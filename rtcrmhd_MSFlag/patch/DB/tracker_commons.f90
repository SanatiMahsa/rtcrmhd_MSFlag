module tracker_commons
  
  use amr_parameters, only: ndim, dp, i8b

  ! Tracker generals
  integer,parameter::ntrack_com_max=500
  logical          ::track_just_read=.true.
  logical          ::track_velocity=.false.
  logical          ::track_quantities=.true.     ! Activates the tracking of properties
  logical          ::track2logfile=.true.        ! Activates the dump of tracker properties to logfile

  integer,parameter::ntrack_inout_basic=11        ! Min number of lines written per tracker in I/O files
  integer,parameter::ntrack_inout_measurements=23 ! Properties number of lines written per tracker I/O
  
  ! Maximum number of trackers (could be changed to allocatable)
  integer,parameter::ntracker_max=1
  integer          ::ntracker=0

  ! Tracker print dummy numbers
  integer,parameter::print_alltrack=1
  integer,parameter::print_xtrack=2
  integer,parameter::print_vtrack=3
  integer,parameter::print_measuretrack=4
  integer,parameter::print_idtrack=5
  integer,parameter::print_internaltrack=6
  
  ! ! Dummies to read into tracker structures
  ! real(dp)::max_dist_track_def(ntracker_max)            ! loads to rmax_track
  ! real(dp)::xtrack_center_mass_def(ntracker_max)        ! loads to xpos_def
  ! real(dp)::aexp_track(ntracker_max)                    ! loads to aexp_ini
  ! real(dp)::ntrack_com_use(ntracker_max)                ! loads to ntrack_use
  ! real(dp)::aexp_sink_com(ntracker_max)                 ! loads to aexp_seed_sink
  ! integer ::do_tracker_measure(ntracker_max)=-1         ! loads to quantities
  ! real(dp)::tracker_radphys_measure(ntracker_max)       ! loads to rad_measure(1)
  ! real(dp)::tracker_radcom_measure(ntracker_max)        ! loads to rad_measure(2)  
  
  ! Tracker structure (contains all information for each tracker)
  type part_tracker
     ! Tracked info and configuration
     logical ::active=.false.
     logical ::quantities=.true.                         ! Activates the measurement of tracker quantities
     logical ::merged=.false.
     logical ::has_sink=.false.
     integer ::ntrack_use                                ! Target number of tracked particles
     real(dp)::aexp_ini=0.0d0
     ! ######################
     real(dp)::rmax_track=0.0d0                          ! Maximum radius to link particles to tracker
     real(dp)::rmax_track2=0.0d0                         ! rmax_track squared (for speedup)
     real(dp)::rad_measure(2)=0.0d0                      ! Measurement region radius (physical, comoving)
     real(dp)::rmeasure=0.0d0                            ! Maximum radius to measure tracker properties
     real(dp)::rmeasure2=0.0d0                           ! rmeasure squared (for speedup)
     real(dp)::delta_bound=600.0d0                       ! Boundary thickness in parsecs
     real(dp)::delta_bound_code=0.0d0                    ! Boundary thickness in code units
     real(dp)::rbound2=0.0d0                             ! Boundary radius squared
     real(dp)::rtrack=0.0d0                              ! Distance to farthest particle tracked 
     
     ! Tracker kinematics
     real(dp)::xpos(ndim)=0.5d0
     real(dp)::xpos_def(ndim)=0.5d0
     real(dp)::vel(ndim)=0.0d0
     real(dp)::track_mass=0.0d0

     ! Tracker sinks
     real(dp)::aexp_seed_sink=0.0d0
     
     ! Tracker measurements (only present in myid==1 processor!)
     real(dp)::rvir=0.0d0                               ! Virial radius (currently unused)
     real(dp)::mvir=0.0d0                               ! Virial mass (currently unused)
     ! ######################
     logical::part_boundary_average=.false.             ! Confirms particle flows are averaged
     real(dp)::mdm=0.0d0                                ! Dark matter mass within rmeasure
     real(dp)::mdm_inflow=0.0d0                         ! Inflowing dm mass into rmeasure
     real(dp)::mdm_outflow=0.0d0                        ! Outflowing dm mass out of rmeasure
     real(dp)::mstar=0.0d0                              ! Stars mass within rmeasure
     real(dp)::mstar_inflow=0.0d0                       ! Inflowing stars mass into rmeasure
     real(dp)::mstar_outflow=0.0d0                      ! Outflowing stars mass out of rmeasure
     ! ######################
     real(dp)::mgas=0.0d0                               ! Gas mass within rmeasure
     real(dp)::mgas_inflow=0.0d0                        ! Inflowing gas mass into rmeasure
     real(dp)::mgas_outflow=0.0d0                       ! Outflowing gas mass out of rmeasure
     real(dp)::mmetal=0.0d0                             ! Metal mass within rmeasure
     real(dp)::mmetal_inflow=0.0d0                      ! Inflowing metal mass into rmeasure
     real(dp)::mmetal_outflow=0.0d0                     ! Outflowing metal mass out of rmeasure
     ! ######################
     logical::part_angular_specific=.false.             ! Confirms angular momenta are specific
     real(kind=8),dimension(ndim)::Lgas=0.0d0           ! Gas specific angular momentum
     real(kind=8),dimension(ndim)::Ldm=0.0d0            ! DM specific angular momentum
     real(kind=8),dimension(ndim)::Lstar=0.0d0          ! Stars specific angular momentum 
     real(kind=8),dimension(ndim)::Lgas_in=0.0d0        ! Gas specific angular momentum inflowing
     real(kind=8),dimension(ndim)::Ldm_in=0.0d0         ! DM specific angular momentum inflowing
     real(kind=8),dimension(ndim)::Lstar_in=0.0d0       ! Stars specific angular momentum inflowing
     real(kind=8),dimension(ndim)::Lgas_out=0.0d0       ! Gas specific angular momentum outflowing
     real(kind=8),dimension(ndim)::Ldm_out=0.0d0        ! DM specific angular momentum outflowing
     real(kind=8),dimension(ndim)::Lstar_out=0.0d0      ! Stars specific angular momentum outflowing
     ! ######################
     real(dp)::ethermal=0.0d0                           ! Gas thermal energy within rmeasure
     real(dp)::ekin_bulk=0.0d0                          ! Gas kinetic energy within rmeasure
     real(dp)::ekin=0.0d0                               ! Gas kinetic energy in COM frame within rmeasure
     real(dp)::emag=0.0d0                               ! Gas magnetic energy within rmeasure
     real(dp)::ecrs=0.0d0                               ! Gas cosmic rays energy within rmeasure
     real(dp)::ethermal_in=0.0d0                        ! Inflowing gas thermal energy within rmeasure
     real(dp)::ekin_bulk_in=0.0d0                       ! Inflowing gas kinetic energy within rmeasure
     real(dp)::ekin_in=0.0d0                            ! Inflowing gas kinetic energy in COM frame within rmeasure
     real(dp)::emag_in=0.0d0                            ! Inflowing gas magnetic energy within rmeasure
     real(dp)::ecrs_in=0.0d0                            ! Inflowing gas cosmic rays energy within rmeasure
     real(dp)::ethermal_out=0.0d0                       ! Outflowing gas thermal energy within rmeasure
     real(dp)::ekin_bulk_out=0.0d0                      ! Outflowing gas kinetic energy within rmeasure
     real(dp)::ekin_out=0.0d0                           ! Outflowing gas kinetic energy in COM frame within rmeasure
     real(dp)::emag_out=0.0d0                           ! Outflowing gas magnetic energy within rmeasure
     real(dp)::ecrs_out=0.0d0                           ! Outflowing gas cosmic rays energy within rmeasure
     ! Grid level support####
     real(dp),allocatable,dimension(:)  ::mgas_lev,mgas_inflow_lev,mgas_outflow_lev
     real(dp),allocatable,dimension(:)  ::mmetal_lev,mmetal_inflow_lev,mmetal_outflow_lev
     real(dp),allocatable,dimension(:,:)::Lgas_lev,Lgas_inflow_lev,Lgas_outflow_lev
     real(dp),allocatable,dimension(:)  ::ethermal_lev,ethermal_in_lev,ethermal_out_lev
     real(dp),allocatable,dimension(:)  ::ekin_lev,ekin_in_lev,ekin_out_lev
     real(dp),allocatable,dimension(:)  ::ekin_bulk_lev,ekin_bulk_in_lev,ekin_bulk_out_lev
     real(dp),allocatable,dimension(:)  ::emag_lev,emag_in_lev,emag_out_lev
     real(dp),allocatable,dimension(:)  ::ecrs_lev,ecrs_in_lev,ecrs_out_lev
     ! ######################

     ! Tracker particles
     integer ::ntracked=0                                ! Number of currently tracked particles
     integer ::ntracked_cpu=0                            ! Dummy for number of tracked particles in cpu
     integer(i8b),allocatable,dimension(:)::id_track     ! Particle ID (idp) of tracked particles
     real(dp),allocatable,dimension(:)    ::dist_track   ! Distance to COM for tracked particles
     logical,allocatable,dimension(:)     ::matched_part ! Bool to determine whether particle is processed     
  end type part_tracker
  type(part_tracker)::tracker(ntracker_max)

  ! Tracker sink properties
  logical          ::seed_smbh_com=.false.            ! Sink particles are seeded on trackers
  integer          ::nsink_com=0                      ! Maximum number of sinks to be seeded by trackers
  real(dp)         ::smbh_mass_track_mass=1.0d4       ! Default mass com SMBH seed
  
  
end module tracker_commons
