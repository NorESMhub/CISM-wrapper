!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module glc_InitMod

!BOP
! !MODULE: glc_InitMod
! !DESCRIPTION:
!  This module contains the glc initialization method and initializes
!  everything needed by a glc simulation.  Primarily it is a driver
!  that calls individual initialization routines for each glc module.
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  SVN:$Id: POP_InitMod.F90 808 2006-04-28 17:06:38Z njn01 $
!  Adapted by William Lipscomb from POP_InitMod.F90

! !USES:

   use glc_kinds_mod
   use glc_communicate, only: my_task, master_task
   use glc_broadcast,   only: broadcast_scalar, broadcast_array
   use glc_time_management, only: iyear0, imonth0, iday0, elapsed_days0,  &
                                  iyear,  imonth,  iday,  elapsed_days,   &
                                  ihour,  iminute, isecond, nsteps_total, &
                                  ymd2eday, eday2ymd, runtype
   use glc_constants, only: nml_in, stdout, zero_gcm_fluxes, test_coupling
   use glc_io,        only: glc_io_read_restart_time
   use glc_files,     only: nml_filename
   use glc_exit_mod
   use shr_kind_mod,  only: CL=>SHR_KIND_CL
   use shr_sys_mod, only: shr_sys_flush

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: glc_initialize

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     module variables
!
!-----------------------------------------------------------------------

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: glc_initialize
! !INTERFACE:

 subroutine glc_initialize(EClock)

! !DESCRIPTION:
!  This routine is the initialization driver that initializes a glc run 
!  by calling individual module initialization routines.
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  same as module

! !USES:
   use glad_main

   use glc_fields, only: glc_allocate_fields, ice_sheet,   &
                         tsfc, qsmb,       &
                         ice_covered, topo,   rofi,   rofl,  hflx,  &
                         ice_sheet_grid_mask

   use glc_override_frac, only: init_glc_frac_overrides
   use glc_constants
   use glc_communicate, only: init_communicate
   use glc_time_management, only: init_time1, init_time2, dtt, ihour
   use glimmer_log
   use glc_route_ice_runoff, only: set_routing
   use glc_history, only : glc_history_init, glc_history_write
   use glc_indexing, only : glc_indexing_init, nx, ny
   use shr_file_mod, only : shr_file_getunit, shr_file_freeunit
   use esmf, only : ESMF_Clock

! !INPUT/OUTPUT PARAMETERS:

   type(ESMF_Clock),     intent(in)    :: EClock

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

  character(fname_length) ::  &
      paramfile        ! Name of the top-level configuration file
 
  character(fname_length) ::  &
      cesm_restart_file  ! Name of the file to be used for a restart
 
  character(CL) :: &
       ice_flux_routing  ! Code for how solid ice should be routed to ocean or sea ice

  ! Scalars which hold information about the global grid --------------
 
  integer (i4) ::  &
      i,j              ! Array index counters

  integer (i4) :: &
      nml_error        ! namelist i/o error flag

  integer (i4) :: &
      nhour_glad      ! number of hours since start of complete glad/CISM run

  integer (i4) :: &
       av_start_time_restart  ! glad averaging start time

  integer (i4) :: &
       days_this_year  ! days since beginning of year

  integer (i4) :: &
       forcing_start_time

  logical :: &
      cesm_restart = .false. ! Logical flag to pass to CISM, telling it to restart
                             ! from a CESM restart

  logical :: &
      cism_debug   = .false. ! Logical flag to pass to CISM, telling it to output extra
                             ! debug diagnostics

  integer :: unit      ! fileunit passed to CISM

  integer :: climate_tstep  ! climate time step (hours)
  
  integer, parameter :: days_in_year = 365
  
  namelist /cism_params/  paramfile, cism_debug, ice_flux_routing, zero_gcm_fluxes, test_coupling
 
! TODO - Write version info?
!-----------------------------------------------------------------------
!  write version information to output log after output redirection
!-----------------------------------------------------------------------
!!   if (my_task == master_task) then
!!      write(stdout,blank_fmt)
!!      write(stdout,ndelim_fmt)
!!      write(stdout,blank_fmt)
!!      write(stdout,'(a)') ' GLC version xxx '
!!      write(stdout,blank_fmt)
!!      call shr_sys_flush(stdout)
!!   endif

!-----------------------------------------------------------------------
!
!  compute time step and initialize time-related quantities
!
!-----------------------------------------------------------------------
 
   call init_time1
 
!-----------------------------------------------------------------------
!
!  output delimiter to log file
!
!-----------------------------------------------------------------------
 
   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      call shr_sys_flush (stdout)
   endif
 
!--------------------------------------------------------------------
! Initialize ice sheet model, grid, and coupling.
! The following code is largely based on CISM.
!-----------------------------------------------------------------------

   paramfile  = 'unknown_paramfile'

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=cism_params,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif
   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call exit_glc(sigAbort,'ERROR reading cism_params nml')
   endif

   call broadcast_scalar(paramfile,         master_task)
   call broadcast_scalar(cism_debug,        master_task)
   call broadcast_scalar(ice_flux_routing,  master_task)
   call broadcast_scalar(zero_gcm_fluxes,   master_task)
   call broadcast_scalar(test_coupling,     master_task)
   call set_routing(ice_flux_routing)

   if (my_task == master_task) then
      write(stdout,*) 'zero_gcm_fluxes: ', zero_gcm_fluxes
      write(stdout,*) 'test_coupling:   ', test_coupling
   end if

   if (verbose .and. my_task==master_task) then
      write (stdout,*) 'paramfile =   ', paramfile
      write (stdout,*) 'dtt =', dtt
      call shr_sys_flush(stdout)
   endif

   ! Set climate time step

   climate_tstep = nint(dtt/3600._r8)   ! convert from sec to integer hours

   if (verbose .and. my_task==master_task) then
      write (stdout,*) 'climate_tstep (hr) =', climate_tstep
      write (stdout,*) 'Set glimmer_unit =', stdout
      write (stdout,*) 'Initialize glad'
   endif

  ! Set glimmer_unit for diagnostic output from CISM. (Log file is already open)
!  call open_log(unit=101)

  call set_glimmer_unit(stdout)
 
  ! Initialize the ice sheet model

  nhour_glad = 0     ! number of hours glad has run since start of complete simulation
                     ! must be set to correct value if reading from a restart file
 
  call init_glc_frac_overrides()

  ! if this is a continuation run, then set up to read restart file and get the restart time
  if (runtype == 'continue') then
    cesm_restart = .true.
    call glc_io_read_restart_time(nhour_glad, av_start_time_restart, cesm_restart_file)
    call ymd2eday (iyear0, imonth0, iday0, elapsed_days0)
    elapsed_days = elapsed_days0 + nhour_glad/24     
    call eday2ymd(elapsed_days, iyear, imonth, iday)
    ihour = 0
    iminute = 0
    isecond = 0
    nsteps_total = nhour_glad / climate_tstep     
    if (verbose .and. my_task==master_task) then
       write(stdout,*) 'Successfully read restart, nhour_glad =', nhour_glad
       write(stdout,*) 'Initial eday/y/m/d:', elapsed_days0, iyear0, imonth0, iday0
       write(stdout,*) 'eday/y/m/d after restart:', elapsed_days, iyear, imonth, iday
       write(stdout,*) 'nsteps_total =', nsteps_total
       write(stdout,*) 'Initialize glad:'
    endif
  endif

  if (verbose .and. my_task==master_task) then
     write(stdout,*) 'Initialize glad, nhour_glad =', nhour_glad
  endif

  unit = shr_file_getUnit()

  call glad_initialize(ice_sheet,                            &
                       climate_tstep,                        &
                       (/paramfile/),                        &
                       daysinyear = days_in_year,            &
                       start_time = nhour_glad,             &
                       gcm_restart = cesm_restart,           &
                       gcm_restart_file = cesm_restart_file, &
                       gcm_debug = cism_debug,               &
                       gcm_fileunit = unit)

  ! TODO(wjs, 2015-03-24) We will need a loop over instances, either here or around the
  ! call to glc_initialize

  if (cesm_restart) then
     forcing_start_time = av_start_time_restart
  else if (test_coupling) then
     ! This assumes that (1) when test_coupling is true, we take a mass balance time step
     ! at the start of every day, and (2) we're starting at the start of a day (i.e.,
     ! ihour0 = iminute0 = isecond0 = 0).
     forcing_start_time = nhour_glad
  else
     ! BUG(wjs, 2017-04-08, https://github.com/NCAR/CISM/issues/1) This assumes that mass
     ! balance time steps always occur on the year boundary - so the current mass balance
     ! time step started at the beginning of the current year. We'd need to generalize
     ! this to allow mid-year mass balance time steps.
     !
     ! This also assumes that ihour0, iminute0 and isecond0 are all 0
     call ymd2eday(year=0, month=imonth0, day=iday0, eday=days_this_year)
     forcing_start_time = nhour_glad - days_this_year * 24
  end if

  call glad_initialize_instance(ice_sheet, instance_index = 1, &
       my_forcing_start_time = forcing_start_time, &
       test_coupling = test_coupling)

  call glc_indexing_init(ice_sheet, instance_index = 1)
  
  call glc_allocate_fields(nx, ny)

  tsfc(:,:) = 0._r8
  qsmb(:,:) = 0._r8
  
  call glad_get_initial_outputs(ice_sheet, instance_index = 1, &
                                ice_covered = ice_covered, &
                                topo = topo, &
                                rofi = rofi, &
                                rofl = rofl, &
                                hflx = hflx, &
                                ice_sheet_grid_mask = ice_sheet_grid_mask)
  
  call glad_initialization_wrapup(ice_sheet)

!TODO - Implement PDD option
 
   call shr_file_freeunit(unit)

! Do the following:
! For each instance, convert ice_sheet%instances(i)%glide_time to hours and compare to nhour_glad.
! If different: Reset params%instances(i)%next_time, params%start_time, params%next_av_start
! Do this here or in initialise_glad?

  ! If restarting (nhour_glad > 0), recompute the year, month, and day
  ! By default, iyear0 = imonth0 = iday0 = 1 (set in namelist file)
  ! Assume that ihour0 = iminute0 = isecond0 = 0
  ! Note that glad does not handle leap years

  ! Set the message level (1 is the default - only fatal errors)
  ! N.B. Must do this after initialization
 
  call glimmer_set_msg_level(6)
 
!-----------------------------------------------------------------------
!
!  finish computing time-related quantities after restart info
!  available (including iyear, imonth, and iday)
!
!-----------------------------------------------------------------------
 
   call init_time2

!-----------------------------------------------------------------------
!
!  initialize history output
!
!-----------------------------------------------------------------------

   call glc_history_init()

!-----------------------------------------------------------------------
!
!  do initial history output
!
!-----------------------------------------------------------------------

   if (.not. cesm_restart) then
      ! TODO loop over instances
      call glc_history_write(ice_sheet%instances(1), EClock, initial_history=.true.)
   end if
   
!-----------------------------------------------------------------------
!
!  output delimiter to log file
!
!-----------------------------------------------------------------------
 
   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,'(" End of GLC initialization")')
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      call shr_sys_flush (stdout)
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine glc_initialize

!***********************************************************************

 end module glc_InitMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
