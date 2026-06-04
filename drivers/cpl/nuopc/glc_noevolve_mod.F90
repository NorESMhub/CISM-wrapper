module glc_noevolve_mod

  !----------------------------------------------------------------------------
  ! Handles the "noevolve" (data-glacier) path for individual ice sheets within
  ! the CISM NUOPC cap.  A noevolve ice sheet reads static topography/thickness
  ! from a file and computes ice runoff (Fgrg_rofi) from the incoming SMB each
  ! coupling step.  All CISM-specific export fields (heat flux, liquid runoff,
  ! flooding, volume) are zero-filled once at initialization.
  !
  ! Adapted from components/cdeps/dglc/dglc_datamode_noevolve_mod.F90.
  !----------------------------------------------------------------------------

  use ESMF             , only : ESMF_State, ESMF_Mesh, ESMF_DistGrid, ESMF_Field
  use ESMF             , only : ESMF_StateGet, ESMF_FieldGet
  use ESMF             , only : ESMF_FieldBundle, ESMF_FieldBundleCreate, ESMF_FieldCreate
  use ESMF             , only : ESMF_FieldBundleAdd, ESMF_MESHLOC_ELEMENT, ESMF_TYPEKIND_R8
  use ESMF             , only : ESMF_MeshGet, ESMF_DistGridGet
  use ESMF             , only : ESMF_GridComp, ESMF_GridCompGet
  use ESMF             , only : ESMF_VM, ESMF_VMAllreduce, ESMF_REDUCE_SUM
  use ESMF             , only : ESMF_SUCCESS, ESMF_LogWrite, ESMF_LOGMSG_INFO
  use NUOPC            , only : NUOPC_IsConnected
  use shr_kind_mod     , only : r8=>shr_kind_r8, cl=>shr_kind_cl, cs=>shr_kind_cs
  use shr_log_mod      , only : shr_log_error
  use shr_const_mod    , only : SHR_CONST_RHOICE, SHR_CONST_RHOSW, SHR_CONST_REARTH
  use dshr_methods_mod , only : dshr_state_getfldptr, dshr_fldbun_getfldptr, chkerr
  use pio              , only : file_desc_t, io_desc_t, var_desc_t, iosystem_desc_t
  use pio              , only : pio_openfile, pio_inq_varid, pio_inq_varndims, pio_inq_vardimid
  use pio              , only : pio_inq_dimlen, pio_initdecomp, pio_read_darray, pio_double
  use pio              , only : pio_closefile, pio_freedecomp, PIO_BCAST_ERROR, PIO_NOWRITE
  use pio              , only : pio_seterrorhandling
  use shr_pio_mod      , only : shr_pio_getiosys, shr_pio_getiotype
  use glc_constants    , only : stdout
  use glc_communicate  , only : my_task, master_task

  implicit none
  private

  public :: glc_noevolve_init
  public :: glc_noevolve_advance
  public :: glc_noevolve_zero_cism_fields

  !----------------------------------------------------------------------------
  ! Per-ice-sheet pointer type (field data lives in the ESMF field; we hold
  ! Fortran pointers into it for convenience).
  !----------------------------------------------------------------------------
  type icesheet_ptr_t
     real(r8), pointer :: ptr(:) => null()
  end type icesheet_ptr_t

  ! Export field pointers (indexed 1..num_noevolve)
  type(icesheet_ptr_t), allocatable :: Sg_area(:)
  type(icesheet_ptr_t), allocatable :: Sg_topo(:)
  type(icesheet_ptr_t), allocatable :: Sg_ice_covered(:)
  type(icesheet_ptr_t), allocatable :: Sg_icemask(:)
  type(icesheet_ptr_t), allocatable :: Sg_icemask_coupled_fluxes(:)
  type(icesheet_ptr_t), allocatable :: Fgrg_rofi(:)

  ! Import field pointer (SMB, needed every coupling step)
  type(icesheet_ptr_t), allocatable :: Flgl_qice(:)

  ! Field name constants — must match the names used in glc_import_export.F90
  character(len=*), parameter :: fld_out_area       = 'Sg_area'
  character(len=*), parameter :: fld_out_topo       = 'Sg_topo'
  character(len=*), parameter :: fld_out_ice_cov    = 'Sg_ice_covered'
  character(len=*), parameter :: fld_out_icemask    = 'Sg_icemask'
  character(len=*), parameter :: fld_out_icemask_cf = 'Sg_icemask_coupled_fluxes'
  character(len=*), parameter :: fld_out_rofi_ocn   = 'Fgrg_rofi'
  character(len=*), parameter :: fld_in_qice        = 'Flgl_qice'

  ! CISM-specific export fields that are zero for noevolve ice sheets
  character(len=*), parameter :: fld_out_hflx    = 'Flgg_hflx'
  character(len=*), parameter :: fld_out_rofi_si = 'Figg_rofi'
  character(len=*), parameter :: fld_out_rofl    = 'Fgrg_rofl'
  character(len=*), parameter :: fld_out_flood   = 'Flrr_flood'
  character(len=*), parameter :: fld_out_volr    = 'Flrr_volr'
  character(len=*), parameter :: fld_out_volrmch = 'Flrr_volrmch'

  real(r8), parameter :: thk0 = 1._r8  ! thickness scaling (= 1 in modern CISM)

  integer :: num_icesheets_total ! total number of ice sheets (prognostic + noevolve)

  character(len=*), parameter :: u_FILE_u = __FILE__

!===============================================================================
contains
!===============================================================================

  subroutine glc_noevolve_init(NStateExp, NStateImp, meshes, &
       datafiles, nx_global, ny_global, internal_gridsize, rc)

    !---------------------------------------------------------------------------
    ! Read static topography and thickness for each noevolve ice sheet, compute
    ! the time-invariant export fields (area, topo, ice_covered, masks), and
    ! grab a pointer into the import SMB field used on every coupling step.
    !
    ! Cell areas are computed from the user-specified internal grid spacing
    ! (matching dglc datamode_noevolve convention):
    !   Sg_area = (internal_gridsize / SHR_CONST_REARTH)**2   ! radians^2
    !---------------------------------------------------------------------------

    ! input/output variables
    type(ESMF_State)      , intent(inout) :: NStateExp(:)         ! all ice sheets (including prognostic)
    type(ESMF_State)      , intent(inout) :: NStateImp(:)         ! all ice sheets (including prognostic)
    type(ESMF_Mesh)       , intent(in)    :: meshes(:)            ! all ice sheets (including prognostic)
    character(len=*)      , intent(in)    :: datafiles(:)         ! all ice sheets (including prognostic)
    integer               , intent(in)    :: nx_global(:)         ! all ice sheets (including prognostic)
    integer               , intent(in)    :: ny_global(:)         ! all ice sheets (including prognostic)
    real(r8)              , intent(in)    :: internal_gridsize(:) ! [m] all ice sheets (including prognostic)
    integer               , intent(out)   :: rc

    ! local variables
    type(ESMF_DistGrid)    :: distgrid
    type(ESMF_FieldBundle) :: fldbun_noevolve
    type(ESMF_Field)       :: field_tmp
    type(file_desc_t)      :: pioid
    type(io_desc_t)        :: pio_iodesc
    type(var_desc_t)       :: varid
    integer , pointer      :: gindex(:)
    real(r8), pointer      :: topog(:), thck(:)
    integer                :: ns, ng, lsize, ndims, rcode
    integer , allocatable  :: dimid(:)
    real(r8)               :: rhoi, rhoo, eus, lsrf, usrf
    integer                :: io_type
    type(iosystem_desc_t), pointer :: pio_subsystem
    character(len=*), parameter :: subname = '(glc_noevolve_mod:noevolve_init) '
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ! Set module variables
    num_icesheets_total = size(NStateExp)

    allocate(Sg_area(num_icesheets_total))
    allocate(Sg_topo(num_icesheets_total))
    allocate(Sg_ice_covered(num_icesheets_total))
    allocate(Sg_icemask(num_icesheets_total))
    allocate(Sg_icemask_coupled_fluxes(num_icesheets_total))
    allocate(Fgrg_rofi(num_icesheets_total))
    allocate(Flgl_qice(num_icesheets_total))

    rhoi = SHR_CONST_RHOICE
    rhoo = SHR_CONST_RHOSW
    eus  = 0._r8

    ! Get the GLC PIO iosystem from the shared CESM PIO initialization
    pio_subsystem => shr_pio_getiosys('GLC')
    io_type       =  shr_pio_getiotype('GLC')

    ! Loop over ice sheets and initialize only those that are noevolve
    ice_sheet_loop: do ns = 1, num_icesheets_total

       if (trim(get_icesheet_mode(ns)) /= 'noevolve') cycle

       !--- Grab pointers into the ESMF export fields ---
       call dshr_state_getfldptr(NStateExp(ns), fld_out_area, &
            fldptr1=Sg_area(ns)%ptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call dshr_state_getfldptr(NStateExp(ns), fld_out_topo, &
            fldptr1=Sg_topo(ns)%ptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call dshr_state_getfldptr(NStateExp(ns), fld_out_ice_cov, &
            fldptr1=Sg_ice_covered(ns)%ptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call dshr_state_getfldptr(NStateExp(ns), fld_out_icemask, &
            fldptr1=Sg_icemask(ns)%ptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call dshr_state_getfldptr(NStateExp(ns), fld_out_icemask_cf, &
            fldptr1=Sg_icemask_coupled_fluxes(ns)%ptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       call dshr_state_getfldptr(NStateExp(ns), fld_out_rofi_ocn, &
            fldptr1=Fgrg_rofi(ns)%ptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       Fgrg_rofi(ns)%ptr(:) = 0._r8

       !--- Grab pointer into the SMB import field ---
       if (.not. NUOPC_IsConnected(NStateImp(ns), fieldName=fld_in_qice)) then
          call shr_log_error(subname//': '//fld_in_qice// &
               ' must be connected for noevolve ice sheet', rc=rc)
          return
       end if
       call dshr_state_getfldptr(NStateImp(ns), fld_in_qice, &
            fldptr1=Flgl_qice(ns)%ptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       Flgl_qice(ns)%ptr(:) = 0._r8

       !--- Determine local size and global index from mesh ---
       call ESMF_MeshGet(meshes(ns), elementDistGrid=distgrid, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_DistGridGet(distgrid, localDe=0, elementCount=lsize, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       allocate(gindex(lsize))
       call ESMF_DistGridGet(distgrid, localDe=0, seqIndexList=gindex, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       !--- Cell area (radians^2, constant) ---
       !    Computed from the user-specified internal grid spacing (matches the
       !    dglc datamode_noevolve convention).
       !    SHR_CONST_REARTH is the radius of earth in m
       !    model_internal_gridsize is the internal model gridsize in m
       do ng = 1, lsize
          Sg_area(ns)%ptr(ng) = (internal_gridsize(ns) / SHR_CONST_REARTH)**2
       end do

       !--- Build field bundle to hold topg and thk from file ---
       fldbun_noevolve = ESMF_FieldBundleCreate(rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! "ice thickness" ;
       field_noevolve = ESMF_FieldCreate(meshes(ns), ESMF_TYPEKIND_R8, &
            name='thk', meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleAdd(fldbun_noevolve, (/field_noevolve/), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! "bed topography" ;
       field_noevolve = ESMF_FieldCreate(meshes(ns), ESMF_TYPEKIND_R8, &
            name='topg', meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_FieldBundleAdd(fldbun_noevolve, (/field_noevolve/), rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       !--- Open data file, set up PIO decomposition, read topg and thk ---
       inquire(file=trim(datafiles(ns)), exist=exists)
       if (.not.exists) then
          call shr_log_error(' ERROR: model input file '//trim(datafiles(ns))//' does not exist', rc=rc)
          return
       else
          if (my_task == master_task) then
             write(stdout,'(a,a)')' opening file ',trim(datafiles(ns))
          end if
       end if
       rcode = pio_openfile(pio_subsystem, pioid, io_type, trim(datafiles(ns)), PIO_NOWRITE)
       call pio_seterrorhandling(pioid, PIO_BCAST_ERROR)
       rcode = pio_inq_varid(pioid, 'thk', varid)
       rcode = pio_inq_varndims(pioid, varid, ndims)
       allocate(dimid(ndims))
       rcode = pio_inq_vardimid(pioid, varid, dimid(1:ndims))
       deallocate(dimid)
       call pio_initdecomp(pio_subsystem, pio_double, (/nx_global(ns), ny_global(ns)/), gindex, pio_iodesc)

       ! Read in the data into the appropriate field bundle pointers
       ! Note that Sg_ice_covered(ns)%ptr points into the data for
       ! the Sg_ice_covered field in NStateExp(ns)
       ! Note that Sg_topo(ns)%ptr points into the data for
       ! the Sg_topon NStateExp(ns)
       ! Note that topog is bedrock topography

       call dshr_fldbun_getFldPtr(fldbun_noevolve, 'topg', topog, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       rcode = pio_inq_varid(pioid, 'topg', varid)
       call pio_read_darray(pioid, varid, pio_iodesc, topog, rcode)

       call dshr_fldbun_getFldPtr(fldbun_evolve, 'thk', thck, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       rcode = pio_inq_varid(pioid, 'thk', varid)
       call pio_read_darray(pioid, varid, pio_iodesc, thck, rcode)

       call pio_closefile(pioid)
       call pio_freedecomp(pio_subsystem, pio_iodesc)
       deallocate(gindex)

       !--- Compute static mask / topo fields from topg and thk ---
       do ng = 1, lsize
          if (topog(ng) - eus < (-rhoi/rhoo) * thck(ng)) then
             lsrf = (-rhoi/rhoo) * thck(ng)
          else
             lsrf = topog(ng)
          end if
          usrf = max(0._r8, thck(ng) + lsrf)

          if (thk0 * usrf > 0._r8) then
             Sg_icemask(ns)%ptr(ng) = 1._r8
             Sg_icemask_coupled_fluxes(ns)%ptr(ng) = 1._r8
             Sg_topo(ns)%ptr(ng) = thk0 * usrf
             Sg_ice_covered(ns)%ptr(ng) = merge(1._r8, 0._r8, thk0 * thck(ng) > 0._r8)
          else
             Sg_icemask(ns)%ptr(ng) = 0._r8
             Sg_icemask_coupled_fluxes(ns)%ptr(ng) = 0._r8
             Sg_topo(ns)%ptr(ng) = 0._r8
             Sg_ice_covered(ns)%ptr(ng) = 0._r8
          end if
       end do

    end do ice_sheet_loop

    ! Zero-fill the CISM-specific export fields (heat flux, runoff, etc.)
    call glc_noevolve_zero_cism_fields(NStateExp, rc)
    if (ChkErr(rc,__LINE__,u_FILE_u)) return

    call ESMF_LogWrite(subname//' done for '// &
         trim(int_to_str(num_noevolve))//' noevolve ice sheet(s)', ESMF_LOGMSG_INFO)

  end subroutine glc_noevolve_init

  !===============================================================================

  subroutine glc_noevolve_advance(gcomp, rc)

    !---------------------------------------------------------------------------
    ! Compute Fgrg_rofi for each noevolve ice sheet from the imported SMB
    ! (Flgl_qice), conserving total ice mass by redistributing negative fluxes
    ! across positive-SMB cells.  Called each coupling step.
    !---------------------------------------------------------------------------

    ! input/output variables
    type(ESMF_GridComp), intent(in)  :: gcomp
    integer            , intent(out) :: rc

    ! local variables
    type(ESMF_VM) :: vm
    integer       :: ns, ng, lsize
    real(r8)      :: loc_pos(1), Tot_pos(1)
    real(r8)      :: loc_neg(1), Tot_neg(1)
    real(r8)      :: rat
    character(len=*), parameter :: subname = '(glc_noevolve_mod:noevolve_advance) '
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    call ESMF_GridCompGet(gcomp, vm=vm, rc=rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return

    ice_sheet_loop: do ns = 1, num_icesheets_total

       if (trim(get_icesheet_mode(ns)) /= 'noevolve') cycle

       lsize = size(Fgrg_rofi(ns)%ptr)
       Fgrg_rofi(ns)%ptr(:) = 0._r8
       loc_pos(1) = 0._r8
       loc_neg(1) = 0._r8

       ! Accumulate global positive and negative SMB weighted by area
       do ng = 1, lsize
          if (Sg_icemask_coupled_fluxes(ns)%ptr(ng) > 0._r8) then
             if (Flgl_qice(ns)%ptr(ng) > 0._r8) then
                loc_pos(1) = loc_pos(1) + Flgl_qice(ns)%ptr(ng) * Sg_area(ns)%ptr(ng)
             else if (Flgl_qice(ns)%ptr(ng) < 0._r8) then
                loc_neg(1) = loc_neg(1) + Flgl_qice(ns)%ptr(ng) * Sg_area(ns)%ptr(ng)
             end if
          end if
       end do

       call ESMF_VMAllreduce(vm, senddata=loc_pos, recvdata=Tot_pos, count=1, &
            reduceflag=ESMF_REDUCE_SUM, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       call ESMF_VMAllreduce(vm, senddata=loc_neg, recvdata=Tot_neg, count=1, &
            reduceflag=ESMF_REDUCE_SUM, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return

       ! Distribute mass-conserving runoff to ice-mask cells
       do ng = 1, lsize
          if (Sg_icemask_coupled_fluxes(ns)%ptr(ng) > 0._r8) then
             if (abs(Tot_pos(1)) >= abs(Tot_neg(1))) then
                ! More positive than negative: scale positive down
                if (Flgl_qice(ns)%ptr(ng) > 0._r8 .and. Tot_pos(1) /= 0._r8) then
                   rat = Flgl_qice(ns)%ptr(ng) / Tot_pos(1)
                   Fgrg_rofi(ns)%ptr(ng) = Flgl_qice(ns)%ptr(ng) + rat * Tot_neg(1)
                end if
             else
                ! More negative than positive: scale negative down
                if (Flgl_qice(ns)%ptr(ng) < 0._r8 .and. Tot_neg(1) /= 0._r8) then
                   rat = Flgl_qice(ns)%ptr(ng) / Tot_neg(1)
                   Fgrg_rofi(ns)%ptr(ng) = Flgl_qice(ns)%ptr(ng) + rat * Tot_pos(1)
                end if
             end if
          end if
       end do

    end do ice_sheet_loop

  end subroutine glc_noevolve_advance

  !===============================================================================

  subroutine glc_noevolve_zero_cism_fields(NStateExp, rc)

    !---------------------------------------------------------------------------
    ! Zero-fill the CISM-specific export fields (heat flux, ice->seaice runoff,
    ! liquid runoff, flooding, volume) for all noevolve ice sheets.
    ! Called once during InitializeRealize.
    !---------------------------------------------------------------------------

    ! input/output variables
    type(ESMF_State), intent(inout) :: NStateExp(:)  ! (num_noevolve)
    integer         , intent(out)   :: rc

    ! local variables
    real(r8), pointer :: ptr(:)
    integer :: ns
    character(len=*), parameter :: subname = '(glc_noevolve_mod:noevolve_zero_cism_fields) '
    !---------------------------------------------------------------------------

    rc = ESMF_SUCCESS

    ice_sheet_loop: do ns = 1, num_noevolve

       if (trim(get_icesheet_mode(ns)) /= 'noevolve') cycle

       call dshr_state_getfldptr(NStateExp(ns), fld_out_hflx, fldptr1=ptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ptr(:) = 0._r8

       call dshr_state_getfldptr(NStateExp(ns), fld_out_rofi_si, fldptr1=ptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ptr(:) = 0._r8

       call dshr_state_getfldptr(NStateExp(ns), fld_out_rofl, fldptr1=ptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ptr(:) = 0._r8

       call dshr_state_getfldptr(NStateExp(ns), fld_out_flood, fldptr1=ptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ptr(:) = 0._r8

       call dshr_state_getfldptr(NStateExp(ns), fld_out_volr, fldptr1=ptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ptr(:) = 0._r8

       call dshr_state_getfldptr(NStateExp(ns), fld_out_volrmch, fldptr1=ptr, rc=rc)
       if (chkerr(rc,__LINE__,u_FILE_u)) return
       ptr(:) = 0._r8

    end do ice_sheet_loop

  end subroutine glc_noevolve_zero_cism_fields

  !===============================================================================

  function int_to_str(n) result(s)
    integer, intent(in) :: n
    character(len=10)   :: s
    write(s,'(i0)') n
  end function int_to_str

end module glc_noevolve_mod
