!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine trackers_revision
   use amr_commons,only: aexp
   use tracker_commons
   implicit none
   real(kind=8),parameter::pc2cm=3.08568d18   
  integer::itrack
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  do itrack=1,ntracker     
     ! Activate trackers when required scale factor is reached
     if (.not. tracker(itrack)%active) then        
        if (aexp.gt.tracker(itrack)%aexp_ini) then
           tracker(itrack)%active=.true.
           tracker(itrack)%xpos=tracker(itrack)%xpos_def
           tracker(itrack)%ntracked=-1
        end if
     end if

     ! Recompute measurement distance for tracker
     if (tracker(itrack)%quantities) then
        tracker(itrack)%rmeasure=tracker(itrack)%rad_measure(1)/aexp+tracker(itrack)%rad_measure(2)
        tracker(itrack)%rmeasure2=tracker(itrack)%rmeasure**2

        ! Convert boundary thickness from physical to code units
        tracker(itrack)%delta_bound_code=tracker(itrack)%delta_bound*pc2cm/scale_l
        tracker(itrack)%rbound2=(tracker(itrack)%rmeasure+tracker(itrack)%delta_bound_code)**2

     else
        ! Set all to zero
        tracker(itrack)%rmeasure2=0.0d0
        tracker(itrack)%rbound2=0.0d0
     end if
     
     ! Recompute maximum tracking distance squared (in case rmax_track is updated)
     tracker(itrack)%rmax_track2=tracker(itrack)%rmax_track**2

     ! Process boundaries and specific quantities final computation
     call process_boundaries_and_specifics(itrack)
  end do
end subroutine trackers_revision
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine process_boundaries_and_specifics(itrack)
  use tracker_commons
  implicit none
  integer,intent(in)::itrack    
  real(dp)::dummy_dp
  ! Average inflow/outflow for particle quantities
  if (.not.tracker(itrack)%part_boundary_average) then
     tracker(itrack)%part_boundary_average=.true.
     dummy_dp=tracker(itrack)%delta_bound_code
     if (dummy_dp.gt.0.0d0) then
        ! Dark matter averaging
        tracker(itrack)%mdm_inflow=tracker(itrack)%mdm_inflow/dummy_dp
        tracker(itrack)%mdm_outflow=tracker(itrack)%mdm_outflow/dummy_dp
        tracker(itrack)%Ldm_in=tracker(itrack)%Ldm_in/dummy_dp
        tracker(itrack)%Ldm_out=tracker(itrack)%Ldm_out/dummy_dp
        ! Stars averaging
        tracker(itrack)%mstar_inflow=tracker(itrack)%mstar_inflow/dummy_dp
        tracker(itrack)%mstar_outflow=tracker(itrack)%mstar_outflow/dummy_dp
        tracker(itrack)%Lstar_in=tracker(itrack)%Lstar_in/dummy_dp
        tracker(itrack)%Lstar_out=tracker(itrack)%Lstar_out/dummy_dp
     end if
  end if
  
  ! Normalise angular momenta of particles
  if (.not.tracker(itrack)%part_angular_specific) then
     tracker(itrack)%part_angular_specific=.true.
     if (tracker(itrack)%mstar.gt.0.0d0) tracker(itrack)%Lstar=tracker(itrack)%Lstar/tracker(itrack)%mstar
     if (tracker(itrack)%mstar_inflow.ne.0.0d0)tracker(itrack)%Lstar_in=tracker(itrack)%Lstar_in/tracker(itrack)%mstar_inflow
     if (tracker(itrack)%mstar_outflow.ne.0.0d0)tracker(itrack)%Lstar_out=tracker(itrack)%Lstar_out/tracker(itrack)%mstar_outflow

     if (tracker(itrack)%mdm.gt.0.0d0) tracker(itrack)%Ldm=tracker(itrack)%Ldm/tracker(itrack)%mdm
     if (tracker(itrack)%mdm_inflow.ne.0.0d0) tracker(itrack)%Ldm_in=tracker(itrack)%Ldm_in/tracker(itrack)%mdm_inflow
     if (tracker(itrack)%mdm_outflow.ne.0.0d0) tracker(itrack)%Ldm_out=tracker(itrack)%Ldm_out/tracker(itrack)%mdm_outflow
  end if
  
end subroutine process_boundaries_and_specifics
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine read_tracker_params(nml_ok)
  use amr_commons
  use tracker_commons
  implicit none
  real(dp),parameter::rmax_def=0.15d0
  logical::nml_ok
  integer::itrack
  ! Dummies to read into tracker structures
  real(dp)::max_dist_track_def(ntracker_max)=0.05         ! loads to rmax_track
  real(dp)::xtrack_center_mass_def(ntracker_max)=0.5d0    ! loads to xpos_def(1)
  real(dp)::ytrack_center_mass_def(ntracker_max)=0.5d0    ! loads to xpos_def(2)
  real(dp)::ztrack_center_mass_def(ntracker_max)=0.5d0    ! loads to xpos_def(3)
  real(dp)::aexp_track(ntracker_max)=1.0d0                ! loads to aexp_ini
  real(dp)::ntrack_com_use(ntracker_max)=500              ! loads to ntrack_use
  real(dp)::aexp_sink_com(ntracker_max)=1.0d0             ! loads to aexp_seed_sink
  integer ::do_tracker_measure(ntracker_max)=-1           ! loads to quantities boolean
  real(dp)::tracker_radphys_measure(ntracker_max)=0.0d0   ! loads to rad_measure(1)
  real(dp)::tracker_radcom_measure(ntracker_max)=0.001d0  ! loads to rad_measure(2) 
  real(dp)::tracker_boundary_delta(ntracker_max)=100.0d0  ! loads to delta_bound (in parsecs)

  
  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  !########################################################
  namelist/tracker_params/track_center_mass,track_center_stars &
       & ,track2logfile,track_quantities,track_velocity &
       !###################################################
       & ,ntracker,aexp_track,ntrack_com_use &
       & ,xtrack_center_mass_def,ytrack_center_mass_def,ztrack_center_mass_def &
       & ,max_dist_track_def &
       !###################################################
       & ,do_tracker_measure &
       & ,tracker_radphys_measure,tracker_radcom_measure,tracker_boundary_delta &
       !###################################################
       & ,aexp_sink_com,nsink_com
  !########################################################

  rewind(1)  
  read(1,NML=tracker_params,END=101)
  goto 102
101 write(*,*)' You need to set up namelist &TRACKER_PARAMS in parameter file' 
  call clean_stop 
102 if (.false.) write(*,*) 
  !--------------------------------------------------
  ! Check number of trackers
  !--------------------------------------------------
  if (ntracker.eq.0) then
     write(*,'("Warning: ntracker has to be initialised. Assuming ntracker=1")')
     ntracker=1
  end if

  !--------------------------------------------------
  ! Check that ntracker < ncpu
  !--------------------------------------------------
  if (ntracker.gt.ncpu) then
     ! This is currently limited by the subroutine: recompute_center_mass_particles
     ! Check the communications for the reduction of the distance list of all cpus
     write(*,'("Error: current tracker implementation does not allow for")')
     write(*,'("ntracker (",I6,") > ncpu (",I6,"). Please reduce ntracker ")') ntracker, ncpu
     stop
  end if

  ! Per tracker checks
  do itrack=1,ntracker
     !--------------------------------------------------
     ! Check number of tracked particles
     !--------------------------------------------------
     if (ntrack_com_use(itrack).lt.1) then
        if (myid.eq.1) write(*,'("Error: more than 0 particles are required for tracking")')
        nml_ok=.false.        
     else if (ntrack_com_use(itrack).gt.ntrack_com_max) then
        if (myid.eq.1) then
           write(*,'("Warning: max num of tracking particles surpassed for ",I3)') itrack
           write(*,'("Reducing ntrack=",I6," to ",I6)') ntrack_com_use(itrack), ntrack_com_max
        end if
        ntrack_com_use(itrack)=ntrack_com_max
     end if

     !--------------------------------------------------
     ! Default tracking distance
     !--------------------------------------------------
     if (max_dist_track_def(itrack).eq.0.0d0) then
        write(*,'("Warning: maximum tracking distance defaulted to ",F7.3)') rmax_def
        max_dist_track_def(itrack)=rmax_def
     end if
  end do
  
  !--------------------------------------------------
  ! Default number of tracker sinks is ntracker
  !--------------------------------------------------
  if (seed_smbh_com.and.(nsink_com.eq.0)) then
     nsink_com=ntracker
  end if

  !--------------------------------------------------
  ! Assign namelist data to corresponding tracker
  !--------------------------------------------------
  call scatter2trackers

  !--------------------------------------------------
  ! Allocate tracker quantities
  !--------------------------------------------------
  do itrack=1,ntracker
     allocate(tracker(itrack)%mgas_lev(1:nlevelmax))
     allocate(tracker(itrack)%mgas_inflow_lev(1:nlevelmax))
     allocate(tracker(itrack)%mgas_outflow_lev(1:nlevelmax))

     allocate(tracker(itrack)%Lgas_lev(1:nlevelmax,1:ndim))
     allocate(tracker(itrack)%Lgas_inflow_lev(1:nlevelmax,1:ndim))
     allocate(tracker(itrack)%Lgas_outflow_lev(1:nlevelmax,1:ndim))

     allocate(tracker(itrack)%mmetal_lev(1:nlevelmax))
     allocate(tracker(itrack)%mmetal_inflow_lev(1:nlevelmax))
     allocate(tracker(itrack)%mmetal_outflow_lev(1:nlevelmax))

     allocate(tracker(itrack)%ethermal_lev(1:nlevelmax))
     allocate(tracker(itrack)%ekin_lev(1:nlevelmax))
     allocate(tracker(itrack)%ekin_bulk_lev(1:nlevelmax))
     allocate(tracker(itrack)%emag_lev(1:nlevelmax))
     allocate(tracker(itrack)%ecrs_lev(1:nlevelmax))

     allocate(tracker(itrack)%ethermal_in_lev(1:nlevelmax))
     allocate(tracker(itrack)%ekin_in_lev(1:nlevelmax))
     allocate(tracker(itrack)%ekin_bulk_in_lev(1:nlevelmax))
     allocate(tracker(itrack)%emag_in_lev(1:nlevelmax))
     allocate(tracker(itrack)%ecrs_in_lev(1:nlevelmax))

     allocate(tracker(itrack)%ethermal_out_lev(1:nlevelmax))
     allocate(tracker(itrack)%ekin_out_lev(1:nlevelmax))
     allocate(tracker(itrack)%ekin_bulk_out_lev(1:nlevelmax))
     allocate(tracker(itrack)%emag_out_lev(1:nlevelmax))
     allocate(tracker(itrack)%ecrs_out_lev(1:nlevelmax))

     tracker(itrack)%mgas_lev=0.0d0; tracker(itrack)%mgas_inflow_lev=0.0d0
     tracker(itrack)%mgas_outflow_lev=0.0d0
     tracker(itrack)%Lgas_lev=0.0d0; tracker(itrack)%Lgas_inflow_lev=0.0d0
     tracker(itrack)%Lgas_outflow_lev=0.0d0
     tracker(itrack)%mmetal_lev=0.0d0; tracker(itrack)%mmetal_inflow_lev=0.0d0
     tracker(itrack)%mmetal_outflow_lev=0.0d0

     tracker(itrack)%ethermal_lev=0.0d0; tracker(itrack)%ethermal_in_lev=0.0d0
     tracker(itrack)%ethermal_out_lev=0.0d0

     tracker(itrack)%ekin_lev=0.0d0; tracker(itrack)%ekin_in_lev=0.0d0
     tracker(itrack)%ekin_out_lev=0.0d0

     tracker(itrack)%ekin_bulk_lev=0.0d0
     tracker(itrack)%ekin_bulk_in_lev=0.0d0
     tracker(itrack)%ekin_bulk_out_lev=0.0d0

     tracker(itrack)%emag_lev=0.0d0; tracker(itrack)%emag_in_lev=0.0d0
     tracker(itrack)%emag_out_lev=0.0d0

     tracker(itrack)%ecrs_lev=0.0d0; tracker(itrack)%ecrs_in_lev=0.0d0
     tracker(itrack)%ecrs_out_lev=0.0d0

  end do

contains
  !###################################################################
  subroutine scatter2trackers
    use tracker_commons
    implicit none
    integer::itrack
    integer::ntrackpart
    do itrack=1,ntracker
       ! Assign initial position and scale factor
       tracker(itrack)%xpos_def(1)=xtrack_center_mass_def(itrack)
       tracker(itrack)%xpos_def(2)=ytrack_center_mass_def(itrack)
       tracker(itrack)%xpos_def(3)=ztrack_center_mass_def(itrack)
       tracker(itrack)%aexp_ini=aexp_track(itrack)
       
       ! Assign desired number of tracked particles
       ntrackpart=ntrack_com_use(itrack)
       tracker(itrack)%ntrack_use=ntrackpart

       ! Assign tracker maximum tracking radius
       tracker(itrack)%rmax_track=max_dist_track_def(itrack)
       tracker(itrack)%rmax_track2=tracker(itrack)%rmax_track**2
       
       ! Assign sink properties
       tracker(itrack)%aexp_seed_sink=aexp_sink_com(itrack)

       ! Assign logical for tracker measurements
       if (do_tracker_measure(itrack).eq.1) then
          tracker(itrack)%quantities=.true.
       else if (do_tracker_measure(itrack).eq.0) then
          tracker(itrack)%quantities=.false.
       end if

       ! Assign measurement region for tracker
       tracker(itrack)%rad_measure(1)=tracker_radphys_measure(itrack)
       tracker(itrack)%rad_measure(2)=tracker_radcom_measure(itrack)
       ! Assign measurement boundary region for tracker
       tracker(itrack)%delta_bound=tracker_boundary_delta(itrack)
       
       ! Per tracker allocations
       allocate(tracker(itrack)%id_track(1:ntrackpart))
       tracker(itrack)%id_track=0
    end do
  end subroutine scatter2trackers
  !###################################################################
End subroutine read_tracker_params
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
! Finds all the required particles for each tracker and computes their
! current center of mass position (and velocity if required)
subroutine compute_track_center_mass_position(seed_call)
  use pm_commons
  use amr_commons
  use tracker_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  ! Routine arguments
  logical,intent(in)::seed_call
  ! Internal dummies
  integer::jlevel
  integer::igrid,jgrid,ipart,jpart,nx_loc
  integer::npart1,ip,info
  logical::bool_dummy
  integer::itrack
  
  ! Tracked particles IDs
  integer,dimension(1:nvector),save::ind_part
  ! Execution type
  logical::do_vels
  ! Center of mass variables  
  real(dp),allocatable,dimension(:,:)::new_center_mass,new_vcom_mass
  real(dp),allocatable,dimension(:)::mass_new_center_mass

  ! Communication vectors
  real(dp)::all_v3(3),send_v3(3),all_v1,send_v1
  
  
  ! Verbose print entering routine
  if(verbose) write(*,111) 
 
  ! Check whether any of the trackers has active particles
  bool_dummy=.false.
  do itrack=1,ntracker
     if (tracker(itrack)%active.and.(tracker(itrack)%ntracked.gt.0)) then
        bool_dummy=.true.
        exit
     end if
  end do
  ! If center of mass is not required, exit
  if (.not.bool_dummy) return
  
  ! Determine whether velocity tracking is required
  require_velocity_tracking: if (track_velocity) then
     do_vels=.true. 
  else
     ! Determine whether velocity tracking is required
     do_vels=.false.

     ! Require center of mass velocity for smbh seeding
     smbh_seeding: if (seed_call) then
        ! Deactivate seeding requirement
        if (nsink.ge.nsink_com) seed_smbh_com=.false.
        
        if (sink.and.seed_smbh_com.and.(nsink.lt.nsink_com)) then
           ! Check whether any of the trackers requires velocity tracking
           do itrack=1,ntracker
              if (aexp.gt.tracker(itrack)%aexp_seed_sink) then
                 do_vels=.true.
                 exit
              end if
           end do
        else if (sink.and.fix_sink2com.and.(nsink.gt.0)) then
           ! Require updated tracker properties to fix sink to tracker
           do_vels=.true.
        else
           ! Seed call does not require updating the trackers
           return
        end if
     end if smbh_seeding
     
  end if require_velocity_tracking


  ! Allocations for center of mass vectors
  allocate(new_center_mass(1:ntracker,1:ndim))
  allocate(mass_new_center_mass(1:ntracker))
  if (do_vels) allocate(new_vcom_mass(1:ntracker,1:ndim))
  
  ! Reset center of mass variables
  new_center_mass=0.0d0
  if (do_vels) new_vcom_mass=0.0d0  
  mass_new_center_mass=0.0d0

  ! Setup matched particles check
  do itrack=1,ntracker
     if (tracker(itrack)%active) then
        allocate(tracker(itrack)%matched_part(tracker(itrack)%ntracked))
        tracker(itrack)%matched_part=.false.
     end if
  end do
  
  Overlevels: do jlevel=levelmin,nlevelmax
     ! Compute maximum time step on active region
     if(numbl(myid,jlevel)>0)then
        ! Loop over grids
        ip=0
        igrid=headl(myid,jlevel)
        do jgrid=1,numbl(myid,jlevel)
           npart1=numbp(igrid)   ! Number of particles in the grid
           if(npart1>0)then
              ! Loop over particles
              ipart=headp(igrid)
              do jpart=1,npart1
                 ! Skip DM particles for stellar tracking
                 !if (track_center_stars.and.(tp(ipart).eq.0.0d0)) then
                 if (track_center_stars.and..not.is_star(idp(ipart),mp(ipart),tp(ipart))) then
                    ipart=nextp(ipart)    
                    cycle
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 if(ip==nvector)then
                    call build_track_center_mass(ind_part,ip)
                    ip=0
                 end if                 
                 ipart=nextp(ipart)    ! Go to next particle
              end do
              ! End loop over particles
           end if
           igrid=next(igrid)   ! Go to next grid
        end do
        ! End loop over grids
        if(ip>0) call build_track_center_mass(ind_part,ip)
     end if
  end do Overlevels


#ifndef WITHOUTMPI   

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Need to figure out which option is faster !!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !------------------------------------------------------
  ! Communications limited to active trackers
  !------------------------------------------------------
  ! Compute local processor center of mass
  do itrack=1,ntracker
     
     ! Skip inactive trackers
     if (.not.tracker(itrack)%active) cycle
     
     ! Communicate center of mass
     all_v3=0.0d0; send_v3=new_center_mass(itrack,1:ndim)
     call MPI_ALLREDUCE(send_v3,all_v3,ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     new_center_mass(itrack,1:ndim)=all_v3

     ! Communicate mass used for center of mass
     all_v1=0.0d0; send_v1=mass_new_center_mass(itrack)
     call MPI_ALLREDUCE(send_v1,all_v1,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     mass_new_center_mass(itrack)=all_v1

     ! Communicate center of mass velocity
     if (do_vels) then
        all_v3=0.0d0; send_v3=new_vcom_mass(itrack,1:ndim)
        call MPI_ALLREDUCE(send_v3,all_v3,ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
        new_vcom_mass(itrack,1:ndim)=all_v3
     end if
     
  end do

  ! !------------------------------------------------------
  ! ! All communications for all trackers used at bulk
  ! !------------------------------------------------------
  !   call MPI_ALLREDUCE(new_center_mass,all_center_mass,ndim*ntracker,&
  !        & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  !   call MPI_ALLREDUCE(mass_new_center_mass,total_center_mass,ntracker,&
  !        & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  !   if (do_vels) then
  !      call MPI_ALLREDUCE(new_vcom_mass,all_vcom_mass,ndim*ntracker,&
  !           & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  !   end if
#endif  
  
  ! Compute center of mass for all processors
  do itrack=1,ntracker
     ! Skip inactive trackers
     if (.not.tracker(itrack)%active) cycle
     
     ! Reduce mass-weighted values
     tracker(itrack)%track_mass=mass_new_center_mass(itrack)
     tracker(itrack)%xpos=new_center_mass(itrack,1:ndim)/(mass_new_center_mass(itrack)+tiny(0.0d0))
     if (do_vels) tracker(itrack)%vel=new_vcom_mass(itrack,1:ndim)/(mass_new_center_mass(itrack)+tiny(0.0d0))
     
     ! Announce new tracker centers to logfile
     if (myid.eq.1) call announce_read_trackers(print_xtrack)  

     ! Tracker deallocations
     deallocate(tracker(itrack)%matched_part)
  end do


  ! Require call to smbh seeding using center of mass
  Seeding_tracker_call: if (seed_call.and.sink.and.seed_smbh_com) then
     do itrack=1,ntracker
        if ((nsink.lt.nsink_com).and.(aexp.gt.tracker(itrack)%aexp_seed_sink)) then
           call seed_smbh_center_mass(itrack)
        end if
     end do

     ! Deactivate seeding requirement
     if (nsink.ge.nsink_com) seed_smbh_com=.false.     
     if (seed_smbh_com) then
        bool_dummy=.true. ! all seeds done
        do itrack=1,ntracker
           if (.not. tracker(itrack)%has_sink) bool_dummy=.false. ! seed missing
        end do
        seed_smbh_com=.not.bool_dummy ! seeds required
     end if
     
  end if Seeding_tracker_call

  ! Routine deallocations
  deallocate(new_center_mass)
  deallocate(mass_new_center_mass)
  if (do_vels) deallocate(new_vcom_mass)
  
  ! Bookkeeping
  track_just_read=.false.

111 format('   Entering compute_track_center_mass_position ')

contains
  !###################################################################
  subroutine build_track_center_mass(ind_part,nn)
    ! use amr_commons,only: myid
    ! use pm_commons
    ! use tracker_commons
    implicit none
    integer::i,ipt,nn
    integer(kind=i8b)::id_check
    integer,dimension(1:nvector)::ind_part
    ! Check for each particle received whether it is one of the tracked particles
    do i=1,nn
       id_check=idp(ind_part(i))
       ! Skip void particles
       if (id_check.eq.0) cycle ! Should this be an 'exit' instruction???
       
       tracker_match: do itrack=1,ntracker
          do ipt=1,tracker(itrack)%ntracked

             ! Skip matched particles
             if (tracker(itrack)%matched_part(ipt)) cycle
             
             matched_part: if (abs(id_check).eq.abs(tracker(itrack)%id_track(ipt))) then                
                ! For matches, compute contribution to center of mass and skip to next tracker
                new_center_mass(itrack,:)=new_center_mass(itrack,:)+mp(ind_part(i))*xp(ind_part(i),:)
                mass_new_center_mass(itrack)=mass_new_center_mass(itrack)+mp(ind_part(i))
                if (do_vels) new_vcom_mass(itrack,:)=new_vcom_mass(itrack,:)+mp(ind_part(i))*vp(ind_part(i),:)
                tracker(itrack)%matched_part(ipt)=.true.
                exit ! To next tracker (itrack)
             end if matched_part
             
          end do
       end do tracker_match
       
    end do
  end subroutine build_track_center_mass
  !###################################################################
end subroutine compute_track_center_mass_position
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine recompute_center_mass_particles
  use pm_commons
  use amr_commons
  use tracker_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  ! Internal dummies
  integer::jlevel,idim
  integer::igrid,jgrid,ipart,jpart,nx_loc
  integer::npart1,info
  integer::i
  logical::required_com
  integer::mpi_err
  integer::extra_min_cpu,min_cpu
  integer::icpu,ipt,jpt
  integer::new_ntrack_center_mass
  integer::itrack
  
  ! Local particle distance minimisation arrays
  integer,allocatable,dimension(:)::ip,qip
  integer,allocatable,dimension(:,:),save::ind_part
  real(dp),allocatable,dimension(:,:)::dist_vector2,qdist_vector2
  real(dp),allocatable,dimension(:,:)::tp_vector,mp_vector
  real(dp),allocatable,dimension(:,:,:)::xp_vector,vp_vector
  
  ! integer::ntrack_mine  
  real(dp)::max_dist_track2
  real(dp)::new_dist,min_max_dist,min_min_dist,rtrack_max
  real(dp)::part_dist2,xpart(3)
  logical ::noinsert
  real(dp),allocatable,dimension(:,:)             :: dist2center

  integer(i8b),dimension(ntrack_com_max)::id_track_center_mass
  
  
  ! Communication variables for tracking
  integer::npart_list
  integer::ntrack_center_mass
  logical,dimension(ncpu) :: required_list_cpus
  real(dp),dimension(ncpu):: all_mindist,all_maxdist
  real(dp),allocatable,dimension(:,:)::sendv,recv
  ! Note to self: for such small arrays, may be better to use the static arrays
  real(dp),allocatable,dimension(:,:)           :: all_dist2center
  integer(i8b),allocatable,dimension(:,:)       :: all_id_track_center_mass
  integer(i8b),allocatable,dimension(:)  :: new_ids_track   
  real(dp),allocatable,dimension(:)      :: new_dist2center
  logical::measure_this,test_part

  if(verbose)write(*,111)

  ! Check whether any of the trackers is active
  required_com=.false.
  do itrack=1,ntracker
     if (tracker(itrack)%active) then
        required_com=.true.
        exit
     end if
  end do
  ! If center of mass is not required, exit
  if (.not.required_com) return
  
  ! Allocate required quantities
  allocate(ip(1:ntracker))
  allocate(ind_part(1:ntracker,1:nvector))
  allocate(dist_vector2(1:ntracker,1:nvector))
  if (track_quantities) then
     ! Reset particle measurements
     tracker(:)%part_angular_specific=.false.
     tracker(:)%part_boundary_average=.false.
     tracker(:)%mstar=0.0d0; tracker(:)%mdm=0.0d0
     tracker(:)%mstar_inflow=0.0d0; tracker(:)%mdm_inflow=0.0d0
     tracker(:)%mstar_outflow=0.0d0; tracker(:)%mdm_outflow=0.0d0
     do itrack=1,ntracker
        tracker(itrack)%Lstar=0.0d0; tracker(itrack)%Ldm=0.0d0
        tracker(itrack)%Lstar_in=0.0d0; tracker(itrack)%Ldm_in=0.0d0
        tracker(itrack)%Lstar_out=0.0d0; tracker(itrack)%Ldm_out=0.0d0
     end do
     allocate(qip(1:ntracker))
     allocate(qdist_vector2(1:ntracker,1:nvector))
     allocate(tp_vector(1:ntracker,1:nvector),mp_vector(1:ntracker,1:nvector))
     allocate(xp_vector(1:ntracker,1:nvector,1:ndim),vp_vector(1:ntracker,1:nvector,1:ndim))
     qip=0; qdist_vector2=0.0d0
  end if

  ip=0; ind_part=0; dist_vector2=0.0d0  

  do itrack=1,ntracker
     ! Skip inactive trackers
     if (.not.tracker(itrack)%active) cycle
     npart_list=tracker(itrack)%ntrack_use
     tracker(itrack)%ntracked=0
     tracker(itrack)%id_track=0
     tracker(itrack)%ntracked_cpu=0

     allocate(tracker(itrack)%dist_track(1:npart_list))
     tracker(itrack)%dist_track=1.0d0
  end do

  !--------------------------------------------------------------------------
  ! Compile arrays of closest particles to each tracker for each processor
  !--------------------------------------------------------------------------  
  Overlevels: do jlevel=levelmin,nlevelmax
     ! Compute maximum time step on active region
     if(numbl(myid,jlevel)>0)then
        ip=0
        qip=0
        igrid=headl(myid,jlevel)
        GridLoop: do jgrid=1,numbl(myid,jlevel)
           npart1=numbp(igrid)   ! Number of particles in the grid
           if(npart1>0)then
              ipart=headp(igrid)
              PartLoop: do jpart=1,npart1
                 ! Skip DM particles for stellar tracking
                 !if ((.not.track_quantities).and.track_center_stars.and.(tp(ipart).eq.0.0d0)) then
                 if ((.not.track_quantities).and.track_center_stars &
                  & .and.(.not.is_star(idp(ipart),mp(ipart),tp(ipart)))) then
                    ipart=nextp(ipart)
                    cycle
                 end if
                                 
                 ! Assignment is done on a per-tracker basis
                 do itrack=1,ntracker
                    ! Skip inactive trackers
                    if (.not.tracker(itrack)%active) cycle
                    measure_this=.false.
                    if (track_quantities.and.tracker(itrack)%quantities) measure_this=.true.
                    
#warning "Warning: track_center_mass does not currently allow for periodic tracking!!"
                    ! part_dist2=sum(( xp(ipart,:)-xtrack_center_mass(:) )**2)
                    xpart=xp(ipart,:)-tracker(itrack)%xpos(:)
                    part_dist2=sum(xpart**2)
                    
                    ! Review only particles within tracking distance
                    if (part_dist2 .lt. tracker(itrack)%rmax_track2) then
                       
                       ! Add particle to tracking list
                       test_part=.true.
                       !if (track_center_stars.and.(tp(ipart).eq.0.0d0)) test_part=.false.
                       if (track_center_stars.and. &
                       & (.not. is_star(idp(ipart),mp(ipart),tp(ipart)))) test_part=.false.
                       if (test_part) then
                          ip(itrack)=ip(itrack)+1
                          ind_part(itrack,ip(itrack))=ipart
                          dist_vector2(itrack,ip(itrack))=part_dist2
                       end if

                       ! Separate measuring particles from tracking particles
                       if (measure_this) then
                          qip(itrack)=qip(itrack)+1
                          qdist_vector2(itrack,qip(itrack))=part_dist2
                          tp_vector(itrack,qip(itrack))=tp(ipart)
                          mp_vector(itrack,qip(itrack))=mp(ipart)
                          vp_vector(itrack,qip(itrack),:)=vp(ipart,:)
                          xp_vector(itrack,qip(itrack),:)=xpart
                       end if
                    end if

                    ! When buffer for tracker is full, review distances list   
                    if(ip(itrack).eq.nvector) then
                       call add2distance_list(itrack,ip(itrack))
                       ip(itrack)=0
                    end if
                    ! If required, measure particle properties for current tracker
                    if(measure_this) then
                       if (qip(itrack).eq.nvector) then
                          call measure_trackpart(itrack,qip(itrack))
                          qip(itrack)=0
                       end if
                    end if
                 end do

                 ! Go to next particle
                 ipart=nextp(ipart)    
              end do PartLoop
              ! End loop over particles
           end if
           igrid=next(igrid)   ! Go to next grid
        end do GridLoop
        
        ! Final buffers processing
        do itrack=1,ntracker
           ! Skip inactive trackers           
           if (.not.tracker(itrack)%active) cycle
           if (ip(itrack).ne.0) call add2distance_list(itrack,ip(itrack))
           ! If required, measure particle properties for current tracker
           if(track_quantities.and.tracker(itrack)%quantities) then
              if (qip(itrack).ne.0) call measure_trackpart(itrack,qip(itrack))
           end if
        end do
        
     end if
   end do Overlevels
  
  ! Routine deallocations
  deallocate(ip,ind_part,dist_vector2)
  if (track_quantities) deallocate(qip,tp_vector,mp_vector,vp_vector,xp_vector,qdist_vector2)

  !--------------------------------------------------------------------------
  ! Gather all arrays of distances onto corresponding processors
  !--------------------------------------------------------------------------
  ! Process iproc receives tracker itrack=iproc. This allows doing all
  ! communications at once and computing the required particles afterwards
  ! The current implementation of this method does not allow for
  ! nproc < ntracker, and is limited on startup by the read_tracker_params
  do itrack=1,ntracker
     ! Skip inactive trackers
     if (.not.tracker(itrack)%active) cycle

     ! Vector length for current tracker
     npart_list=tracker(itrack)%ntrack_use
     ! Check whether current processor is expected to receive this tracker data
     if (myid.eq.itrack) then
        ! Allocate receive arrays
        allocate(all_dist2center(npart_list,ncpu))
        allocate(all_id_track_center_mass(npart_list,ncpu))
     end if
     
#ifndef WITHOUTMPI    
     ! Gather all distances for tracker itrack into processor with myid=itrack
     call MPI_GATHER(tracker(itrack)%dist_track,npart_list,MPI_DOUBLE_PRECISION &
          & ,all_dist2center,npart_list,MPI_DOUBLE_PRECISION &
          & ,itrack-1,MPI_COMM_WORLD,mpi_err)
#ifndef LONGINT
     call MPI_GATHER(tracker(itrack)%id_track,npart_list,MPI_INTEGER &
          & ,all_id_track_center_mass,npart_list,MPI_INTEGER &
          & ,itrack-1,MPI_COMM_WORLD,mpi_err)
#else
     call MPI_GATHER(tracker(itrack)%id_track,npart_list,MPI_INTEGER8 &
          & ,all_id_track_center_mass,npart_list,MPI_INTEGER8 &
          & ,itrack-1,MPI_COMM_WORLD,mpi_err)
#endif
#else
     all_dist2center(itrack,1:npart_list)=tracker(itrack)%dist_track
     all_id_track_center(itrack,1:npart_list)=tracker(itrack)%id_track     
#endif
  end do

#ifndef WITHOUTMPI    
  if (track_quantities) then ! Reduce tracker measurements in lead processor
     allocate(sendv(1:ntracker,1:3),recv(1:ntracker,1:3))
     
     ! Stars communications
     sendv(:,1)=tracker(1:ntracker)%mstar
     !write(*,*)"myid, mstar=",myid,tracker(1:ntracker)%mstar
     sendv(:,2)=tracker(1:ntracker)%mstar_inflow
     sendv(:,3)=tracker(1:ntracker)%mstar_outflow
     call MPI_REDUCE(sendv,recv,ntracker*3,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,info)
     tracker(1:ntracker)%mstar        =recv(:,1)
     tracker(1:ntracker)%mstar_inflow =recv(:,2)
     tracker(1:ntracker)%mstar_outflow=recv(:,3)
     do idim=1,ndim
        sendv(:,1)=tracker(1:ntracker)%Lstar(idim)
        sendv(:,2)=tracker(1:ntracker)%Lstar_in(idim)
        sendv(:,3)=tracker(1:ntracker)%Lstar_out(idim)
        call MPI_REDUCE(sendv,recv,ntracker*3,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,info)
        tracker(1:ntracker)%Lstar(idim)    =recv(:,1)
        tracker(1:ntracker)%Lstar_in(idim) =recv(:,2)
        tracker(1:ntracker)%Lstar_out(idim)=recv(:,3)
     end do
     
     ! DM communications
     sendv(:,1)=tracker(1:ntracker)%mdm
     sendv(:,2)=tracker(1:ntracker)%mdm_inflow
     sendv(:,3)=tracker(1:ntracker)%mdm_outflow
     call MPI_REDUCE(sendv,recv,ntracker*3,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,info)
     tracker(1:ntracker)%mdm        =recv(:,1)
     tracker(1:ntracker)%mdm_inflow =recv(:,2)
     tracker(1:ntracker)%mdm_outflow=recv(:,3)
     do idim=1,ndim
        sendv(:,1)=tracker(1:ntracker)%Ldm(idim)
        sendv(:,2)=tracker(1:ntracker)%Ldm_in(idim)
        sendv(:,3)=tracker(1:ntracker)%Ldm_out(idim)
        call MPI_REDUCE(sendv,recv,ntracker*3,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,info)
        tracker(1:ntracker)%Ldm(idim)    =recv(:,1)
        tracker(1:ntracker)%Ldm_in(idim) =recv(:,2)
        tracker(1:ntracker)%Ldm_out(idim)=recv(:,3)
     end do
     
     deallocate(sendv,recv)
  end if
#endif
  if (myid.eq.1) then
     do itrack=1,ntracker
        call process_boundaries_and_specifics(itrack)
     end do
  end if

  !--------------------------------------------------------------------------
  ! Compute lists intersections to determine required vectors
  !--------------------------------------------------------------------------
  ProcReceivedTracker: if (myid.le.ntracker) then
     ! Assign corresponding tracker depending on processor identity
     itrack=myid
          
     ! Skip inactive trackers
     if (.not.tracker(itrack)%active) go to 17

     ! Vector length for current tracker
     npart_list=tracker(itrack)%ntrack_use
     allocate(new_ids_track(1:npart_list))
     allocate(new_dist2center(1:npart_list))

     ! Tracker properties to seggregate particles
     max_dist_track2=tracker(itrack)%rmax_track2
     
     ! Find the processor with the closest max. distance
     do icpu=1,ncpu
        all_mindist(icpu)=minval(all_dist2center(:,icpu))
        all_maxdist(icpu)=maxval(all_dist2center(:,icpu))
     end do

     min_cpu=-1    
     min_max_dist=1.0
     min_min_dist=1.0
     do icpu=1,ncpu
        if (all_maxdist(icpu).lt.min_max_dist) then
           min_max_dist=all_maxdist(icpu)
           min_cpu=icpu
        end if
        if (all_mindist(icpu).lt.min_min_dist) then
           min_min_dist=all_mindist(icpu)
           extra_min_cpu=icpu
        end if
     end do
     if (min_cpu.eq.-1) then
        if (extra_min_cpu.eq.-1) then
           write(*,'("FATAL! minimised track_center missing minlist: ",L)') track_just_read
           stop
           min_cpu=myid              
        else
           min_cpu=extra_min_cpu
        end if
     end if
     
     ! Mark cpus that intersect with closest list for distance review
     required_list_cpus=.false.
     do icpu=1,ncpu
        if (all_mindist(icpu).lt.min_max_dist) then
           required_list_cpus(icpu)=.true.
       end if
     end do

     ! Build final list of distances
     new_dist2center(:)=all_dist2center(:,min_cpu)
     new_ids_track(:)=all_id_track_center_mass(:,min_cpu)
     do icpu=1,ncpu
        if (icpu.eq.min_cpu) cycle ! Skip initial list
        if (required_list_cpus(icpu).eq..false.) cycle ! Skip non-overlap lists
        
        do ipart=1,npart_list
           
           if (all_dist2center(ipart,icpu).gt.max_dist_track2) exit ! Skip large distances
           ! Check for each particle received whether it should be in the list
           do ipt=1,npart_list
              if (all_dist2center(ipart,icpu).le.new_dist2center(ipt)) then
                 ! Check whether ID of particle already in list and skip if yes
                 noinsert=.false.
                 do jpt=max(ipt-1,1),npart_list
                    if (abs(all_id_track_center_mass(ipart,icpu)).eq.&
                         & abs(new_ids_track(jpt))) then
                       noinsert=.true.
                       exit
                    end if
                 end do
                 if (noinsert) exit
                 ! Insert particle in list, then skip to next
                 do jpt=npart_list,ipt+1,-1
                    new_dist2center(jpt)=new_dist2center(jpt-1)
                 end do
                 new_dist2center(ipt)=all_dist2center(ipart,icpu)
                 new_ids_track(ipt)=all_id_track_center_mass(ipart,icpu)
                 exit
              end if
           end do
        end do
        
     end do
     
     ! Compute number of tracked particles
     new_ntrack_center_mass=npart_list
     do ipart=1,npart_list
        if (new_dist2center(ipart).gt.max_dist_track2) then
           new_ntrack_center_mass=ipart-1
           exit
        end if
     end do

     ! Save the distance to the farthest particle
     rtrack_max=new_dist2center(new_ntrack_center_mass)
     
     ! Deallocate all_cpus lists
     deallocate(all_dist2center)
     deallocate(all_id_track_center_mass)

  end if ProcReceivedTracker
17 continue

  !--------------------------------------------------------------------------  
  ! Scatter final lists of tracked IDs to all processors
  !--------------------------------------------------------------------------
  FinalTrackerScatter: do itrack=1,ntracker
     ! Skip inactive trackers
     if (.not.tracker(itrack)%active) cycle

     ! Vector length for current tracker
     npart_list=tracker(itrack)%ntrack_use
     ! Communicating processor prepares send data in bcast arrays
     if (myid.le.ntracker) then
        if (itrack.eq.myid) then
           id_track_center_mass(1:npart_list)=new_ids_track(1:npart_list)
           ntrack_center_mass=new_ntrack_center_mass
        end if
     end if
     
#ifndef WITHOUTMPI    
     call MPI_BCAST(ntrack_center_mass,1,MPI_INTEGER,itrack-1,MPI_COMM_WORLD,info)

#ifndef LONGINT
     call MPI_BCAST(id_track_center_mass(1:npart_list),npart_list,MPI_INTEGER,itrack-1,MPI_COMM_WORLD,info)
#else
     call MPI_BCAST(id_track_center_mass(1:npart_list),npart_list,MPI_INTEGER8,itrack-1,MPI_COMM_WORLD,info)
#endif
     
     call MPI_BCAST(rtrack_max,1,MPI_DOUBLE_PRECISION,itrack-1,MPI_COMM_WORLD,info)
#else
     id_track_center_mass=new_ids_track
     ntrack_center_mass=new_ntrack_center_mass
#endif

     ! Save received tracker particles onto corresponding tracker
     tracker(itrack)%id_track=id_track_center_mass
     tracker(itrack)%ntracked=ntrack_center_mass
     tracker(itrack)%rtrack=rtrack_max
  end do FinalTrackerScatter
  
  ! Final deallocations
  do itrack=1,ntracker
     ! Skip inactive trackers
     if (.not.tracker(itrack)%active) cycle
     deallocate(tracker(itrack)%dist_track)
  end do
   if (allocated(new_ids_track)) deallocate(new_ids_track)
   if (allocated(new_dist2center)) deallocate(new_dist2center)
   
   
  ! Bookkeeping
  track_just_read=.false.
111 format('   Entering recompute_center_mass_particles ')

contains
  !###################################################################
  subroutine add2distance_list(it,nn)
    implicit none
    integer::i,ipt,jpt,ntr
    integer,intent(in)::it,nn
    logical::noinsert
    integer::nmax_list

#warning "dist_vector2 and ind_part are looped here in row-major. This should be column-major!!!!"    
    ! Check for each particle received whether it should be in the minimum distance list
    do i=1,nn

       ! First matched particle gets inserted regardless, then skip to next
       if (tracker(it)%ntracked_cpu.eq.0) then
          tracker(it)%dist_track(1)=dist_vector2(it,i)
          tracker(it)%id_track(1)=idp(ind_part(it,i))
          tracker(it)%ntracked_cpu=1
          cycle
       end if

       ! For incomplete lists, add particles until full
       nmax_list=min(tracker(it)%ntracked_cpu+1,tracker(it)%ntrack_use)
       ntr=tracker(it)%ntracked_cpu
       ! nmax_list=min(ntrack_mine+1,ntrack_com_max)       
       do ipt=1,nmax_list
          if (dist_vector2(it,i).lt.tracker(it)%dist_track(ipt)) then
             ! Check whether ID of particle already in list and skip if yes
             noinsert=.false.
             do jpt=max(ipt-1,1),ntr
                if (abs(idp(ind_part(it,i))).eq.abs(tracker(it)%id_track(jpt))) then
                   noinsert=.true.                 
                   exit
                end if
             end do
             ! Exit for particles not due to insert
             if (noinsert) exit
             
             ! Insert particle in list, then skip to next
             do jpt=ntr,ipt+1,-1
                ! dist2center(jpt)=dist2center(jpt-1)
                tracker(it)%dist_track(jpt)=tracker(it)%dist_track(jpt-1)
             end do
             tracker(it)%id_track(ipt)=idp(ind_part(it,i)) ! Insert ID
             tracker(it)%dist_track(ipt)=dist_vector2(it,i) ! Insert distance
             tracker(it)%ntracked_cpu=nmax_list ! Increase particle count in list
             exit
          end if
       end do
    end do
  end subroutine add2distance_list
  !###################################################################
  subroutine measure_trackpart(it,nn)
    implicit none
    integer,intent(in)::it,nn
    real(kind=8)::mpart,vpart_rad,mpart_flow
    real(kind=8),dimension(ndim)::xpart,vpart,Lpart
    do i=1,nn       
       sphere_cut: if (qdist_vector2(it,i).lt.tracker(it)%rmeasure2) then ! Contribute to the tracker region
          mpart=mp_vector(it,i)
          if (mpart.eq.0.0d0) cycle
          vpart=vp_vector(it,i,:)-tracker(it)%vel
          xpart=xp_vector(it,i,:) ! Already in tracker CoM frame
          ! Compute specific angular momentum of particle
          Lpart(1) = xpart(2)*vpart(3)-xpart(3)*vpart(2)
          Lpart(2) = xpart(3)*vpart(1)-xpart(1)*vpart(3)
          Lpart(3) = xpart(1)*vpart(2)-xpart(2)*vpart(1)

          if (tp_vector(it,i).eq.0.0d0) then ! Is DM
          !if (tp_vector(it,ip(it)).eq.0.0d0) then ! Is DM
             tracker(it)%mdm=tracker(it)%mdm+mpart
             tracker(it)%Ldm=tracker(it)%Ldm+mpart*Lpart
          else ! Is star
             tracker(it)%mstar=tracker(it)%mstar+mpart
             tracker(it)%Lstar=tracker(it)%Lstar+mpart*Lpart
          end if
       else if (qdist_vector2(it,i).lt.tracker(it)%rbound2) then ! Contribute to the tracker boundary 
          vpart=vp_vector(it,i,:)-tracker(it)%vel
          xpart=xp_vector(it,i,:) ! Already in tracker CoM frame

          ! Compute radial flow of particle
          vpart_rad=sum(vpart*xpart)/sqrt(qdist_vector2(it,i))
          mpart_flow=mpart*vpart_rad

          ! Compute specific angular momentum of particle
          Lpart(1) = xpart(2)*vpart(3)-xpart(3)*vpart(2)
          Lpart(2) = xpart(3)*vpart(1)-xpart(1)*vpart(3)
          Lpart(3) = xpart(1)*vpart(2)-xpart(2)*vpart(1)
          
          ! Contribute to inflow / outflow quantities
          if (vpart_rad.gt.0.0d0) then
             if (tp_vector(it,i).eq.0.0d0) then ! Is DM
             !if (tp_vector(it,ip(it)).eq.0.0d0) then ! Is DM
                tracker(it)%mdm_outflow=tracker(it)%mdm_outflow+mpart_flow
                tracker(it)%Ldm_out=tracker(it)%Ldm_out+mpart_flow*Lpart
             else ! Is star
                tracker(it)%mstar_outflow=tracker(it)%mstar_outflow+mpart_flow
                tracker(it)%Lstar_out=tracker(it)%Lstar_out+mpart_flow*Lpart
             end if
          else
             if (tp_vector(it,i).eq.0.0d0) then ! Is DM
             !if (tp_vector(it,ip(it)).eq.0.0d0) then ! Is DM
                tracker(it)%mdm_inflow=tracker(it)%mdm_inflow+mpart_flow
                tracker(it)%Ldm_in=tracker(it)%Ldm_in+mpart_flow*Lpart
             else ! Is star
                tracker(it)%mstar_inflow=tracker(it)%mstar_inflow+mpart_flow
                tracker(it)%Lstar_in=tracker(it)%Lstar_in+mpart_flow*Lpart
             end if
          end if
       end if sphere_cut
    end do
  end subroutine measure_trackpart
  !###################################################################
end subroutine recompute_center_mass_particles
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine output_tracks(filename)
  use pm_commons
  use tracker_commons
  implicit none
  integer::it
  character(len=80),intent(in)::filename
  character(len=80)::fileloc
  integer::ilun,ntrack_part,dummy
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  if(verbose)write(*,*)'Entering output_track'
  ilun=12
  ! Open trackers file
  fileloc=trim(filename)  
  open(unit=ilun,file=fileloc,form='formatted')
  ! Output trackers header
  write(ilun,'("ntrackers   =",I)') ntracker
  if (track_quantities) then
     write(ilun,'("ntrack_inout=",I)') ntrack_inout_basic+ntrack_inout_measurements
     write(ilun,'("do_quants   =",I)') 1
  else
     write(ilun,'("ntrack_inout=",I)') ntrack_inout_basic
     write(ilun,'("do_quants   =",I)') 0
  end if
  write(ilun,'("unit_l      =",E14.7)') scale_l
  write(ilun,'("unit_d      =",E14.7)') scale_d
  write(ilun,'("unit_t      =",E14.7)') scale_t
  write(ilun,'("")')
  ! Dump all trackers
  do it=1,ntracker
     ! Tracker basics
     write(ilun,'("itrack      =",I)') it
     if (tracker(it)%active) then
        write(ilun,'("active      =",I)') 1
     else
        write(ilun,'("active      =",I)') 0
     end if

     if (track_quantities) then
        if (tracker(it)%quantities) then
           write(ilun,'("quantities  =",I)') 1
        else
           write(ilun,'("quantities  =",I)') 0
        end if
     end if
     ! Tracker position
     write(ilun,'("ntrack      =",I11)') tracker(it)%ntracked
     write(ilun,'("xtrack      =",F11.8)') tracker(it)%xpos(1)
     write(ilun,'("ytrack      =",F11.8)') tracker(it)%xpos(2)
     write(ilun,'("ztrack      =",F11.8)') tracker(it)%xpos(3)
     ! Tracker center of mass velocity
     write(ilun,'("vcom        =",3(E14.7," "))') tracker(it)%vel
     ! Tracker properties     
     write(ilun,'("rmax        =",E14.7)') tracker(it)%rmax_track
     write(ilun,'("rvir        =",E14.7)') tracker(it)%rvir
     write(ilun,'("mvir        =",E14.7)') tracker(it)%mvir
     ! Tracker measurements log
     if (tracker(it)%quantities) then
        write(ilun,'("rphys_log   =",E14.7)') tracker(it)%rad_measure(1)
        write(ilun,'("rcomov_log  =",E14.7)') tracker(it)%rad_measure(2)
        write(ilun,'("delta_bound =",E14.7)') tracker(it)%delta_bound
        write(ilun,'("dbound_code =",E14.7)') tracker(it)%delta_bound_code
        write(ilun,'("rtrack      =",E14.7)') tracker(it)%rtrack
        ! Particle masses budget
        write(ilun,'("mdm         =",3(E14.7," "))') tracker(it)%mdm, tracker(it)%mdm_inflow, tracker(it)%mdm_outflow
        write(ilun,'("Ldm         =",3(E14.7," "))') tracker(it)%Ldm
        write(ilun,'("Ldm_in      =",3(E14.7," "))') tracker(it)%Ldm_in
        write(ilun,'("Ldm_out     =",3(E14.7," "))') tracker(it)%Ldm_out
        write(ilun,'("mstar       =",3(E14.7," "))') tracker(it)%mstar, tracker(it)%mstar_inflow, tracker(it)%mstar_outflow
        write(ilun,'("Lstar       =",3(E14.7," "))') tracker(it)%Lstar
        write(ilun,'("Lstar_in    =",3(E14.7," "))') tracker(it)%Lstar_in
        write(ilun,'("Lstar_out   =",3(E14.7," "))') tracker(it)%Lstar_out
        ! Gas masses budget
        write(ilun,'("mgas        =",3(E14.7," "))') tracker(it)%mgas, tracker(it)%mgas_inflow, tracker(it)%mgas_outflow
        write(ilun,'("Lgas        =",3(E14.7," "))') tracker(it)%Lgas
        write(ilun,'("Lgas_in     =",3(E14.7," "))') tracker(it)%Lgas_in
        write(ilun,'("Lgas_out    =",3(E14.7," "))') tracker(it)%Lgas_out
        write(ilun,'("mmetal      =",3(E14.7," "))') tracker(it)%mmetal, tracker(it)%mmetal_inflow, tracker(it)%mmetal_outflow
        ! Gas energies
        write(ilun,'("ethermal    =",3(E14.7," "))') tracker(it)%ethermal,tracker(it)%ethermal_in,tracker(it)%ethermal_out  
        write(ilun,'("ekin_bulk   =",3(E14.7," "))') tracker(it)%ekin_bulk,tracker(it)%ekin_bulk_in,tracker(it)%ekin_bulk_out
        write(ilun,'("ekin        =",3(E14.7," "))') tracker(it)%ekin,tracker(it)%ekin_in,tracker(it)%ekin_out
        write(ilun,'("emag        =",3(E14.7," "))') tracker(it)%emag,tracker(it)%emag_in,tracker(it)%emag_out
        write(ilun,'("ecrs        =",3(E14.7," "))') tracker(it)%ecrs,tracker(it)%ecrs_in,tracker(it)%ecrs_out
     end if        
  end do
  write(ilun,'("")')
  
  ! Tracker particle IDs  
  write(ilun,'("Tracked particle IDs")') 
  do it=1,ntracker
     ntrack_part=tracker(it)%ntracked
     write(ilun,'("itrack      =",I)') it
     write(ilun,'("ntrack      =",I)') ntrack_part
     if (ntrack_part.gt.0) then
        write(ilun,*) tracker(it)%id_track(1:ntrack_part)
     else
        write(ilun,*) 
     end if
  end do
  close(ilun)  
end subroutine output_tracks
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine input_tracks(filename,myid)
  use pm_commons
  use tracker_commons
  implicit none
  character(len=80)::filename,fileloc
  integer::ilun
  integer,intent(in)::myid
  integer::ntrack_reads,nread_req
  logical::has_quants
  ! Internal dummies
  integer::dummy_int,it,ntrack_part
  real(dp)::dummy_dp,dummy_dp1,dummy_dp2

  ilun=12

  ! Open trackers file
  fileloc=trim(filename)  
  open(unit=ilun,file=fileloc,form='formatted')
  ! Read trackers header
  read(ilun,'("ntrackers   =",I)') ntracker
  read(ilun,'("ntrack_inout=",I)') ntrack_reads     
  read(ilun,'("do_quants   =",I)') dummy_int
  read(ilun,'("unit_l      =",E14.7)') ! scale_l
  read(ilun,'("unit_d      =",E14.7)') ! scale_d
  read(ilun,'("unit_t      =",E14.7)') ! scale_t
  read(ilun,'("")')
  has_quants=.false.
  nread_req=ntrack_inout_basic
  if (track_quantities) nread_req=nread_req+ntrack_inout_measurements
  if (dummy_int.eq.1) has_quants=.true.
  if (ntrack_reads.ne.nread_req) write(*,212) ntrack_reads, nread_req
212 format("Warning: expected nreads differs: ",I3,"/",I3)
  
  ! Read all trackers
  do it=1,ntracker
     ! Tracker basics
     read(ilun,'("itrack      =",I)') ! Skip current tracker number
     read(ilun,'("active      =",I)') dummy_int
     tracker(it)%active=.false.
     if (dummy_int.eq.1) tracker(it)%active=.true.
     read(ilun,'("quantities  =",I)') dummy_int
     if (dummy_int.eq.1) tracker(it)%quantities=.true.
     ! Tracker position
     read(ilun,'("ntrack      =",I11)') tracker(it)%ntracked
     read(ilun,'("xtrack      =",F11.8)') tracker(it)%xpos(1)
     read(ilun,'("ytrack      =",F11.8)') tracker(it)%xpos(2)
     read(ilun,'("ztrack      =",F11.8)') tracker(it)%xpos(3)
     tracker(1)%xpos_def=tracker(1)%xpos
     ! Tracker center of mass velocity
     read(ilun,'("vcom        =",3(E14.7," "))') tracker(it)%vel
     ! Tracker properties
     read(ilun,'("rmax        =",E14.7)') dummy_dp 
     if (tracker(it)%rmax_track.eq.0.0d0) tracker(it)%rmax_track=dummy_dp
     read(ilun,'("rvir        =",E14.7)') tracker(it)%rvir
     read(ilun,'("mvir        =",E14.7)') tracker(it)%mvir
     ! Tracker measurements log
     if (tracker(it)%quantities) then
        read(ilun,'("rphys_log   =",E14.7)') dummy_dp1
        read(ilun,'("rcomov_log  =",E14.7)') dummy_dp2
        if ((tracker(it)%rad_measure(1).eq.0.0d0).and.(tracker(it)%rad_measure(2).eq.0.0d0)) then
           tracker(it)%rad_measure(1)=dummy_dp1
           tracker(it)%rad_measure(2)=dummy_dp2
        end if
        read(ilun,'("delta_bound =",E14.7)') tracker(it)%delta_bound
        read(ilun,'("dbound_code =",E14.7)') tracker(it)%delta_bound_code
        read(ilun,'("rtrack      =",E14.7)') tracker(it)%rtrack
        ! Particle masses budget
        read(ilun,'("mdm         =",3(E14.7," "))') tracker(it)%mdm, tracker(it)%mdm_inflow, tracker(it)%mdm_outflow
        read(ilun,'("Ldm         =",3(E14.7," "))') tracker(it)%Ldm
        read(ilun,'("Ldm_in      =",3(E14.7," "))') tracker(it)%Ldm_in
        read(ilun,'("Ldm_out     =",3(E14.7," "))') tracker(it)%Ldm_out
        read(ilun,'("mstar       =",3(E14.7," "))') tracker(it)%mstar, tracker(it)%mstar_inflow, tracker(it)%mstar_outflow
        read(ilun,'("Lstar       =",3(E14.7," "))') tracker(it)%Lstar
        read(ilun,'("Lstar_in    =",3(E14.7," "))') tracker(it)%Lstar_in
        read(ilun,'("Lstar_out   =",3(E14.7," "))') tracker(it)%Lstar_out
        ! Gas masses budget
        read(ilun,'("mgas        =",3(E14.7," "))') tracker(it)%mgas, tracker(it)%mgas_inflow, tracker(it)%mgas_outflow
        read(ilun,'("Lgas        =",3(E14.7," "))') tracker(it)%Lgas
        read(ilun,'("Lgas_in     =",3(E14.7," "))') tracker(it)%Lgas_in
        read(ilun,'("Lgas_out    =",3(E14.7," "))') tracker(it)%Lgas_out
        read(ilun,'("mmetal      =",3(E14.7," "))') tracker(it)%mmetal, tracker(it)%mmetal_inflow, tracker(it)%mmetal_outflow
        ! Gas energies
        read(ilun,'("ethermal    =",3(E14.7," "))') tracker(it)%ethermal,tracker(it)%ethermal_in,tracker(it)%ethermal_out  
        read(ilun,'("ekin_bulk   =",3(E14.7," "))') tracker(it)%ekin_bulk,tracker(it)%ekin_bulk_in,tracker(it)%ekin_bulk_out
        read(ilun,'("ekin        =",3(E14.7," "))') tracker(it)%ekin,tracker(it)%ekin_in,tracker(it)%ekin_out
        read(ilun,'("emag        =",3(E14.7," "))') tracker(it)%emag,tracker(it)%emag_in,tracker(it)%emag_out
        read(ilun,'("ecrs        =",3(E14.7," "))') tracker(it)%ecrs,tracker(it)%ecrs_in,tracker(it)%ecrs_out
     end if        
  end do
  read(ilun,'("")')


  ! Read tracker particle IDs  
  read(ilun,'("Tracked particle IDs")') 
  do it=1,ntracker
     read(ilun,'("itrack      =",I)') 
     read(ilun,'("ntrack      =",I)') ntrack_part
     
     tracker(it)%ntracked=min(ntrack_com_max,ntrack_part)
     if (ntrack_part.gt.0) then
        read(ilun,*) tracker(it)%id_track(1:ntrack_part)
     else
        read(ilun,*) 
     end if
  end do
  close(ilun)
  track_just_read=.true.
end subroutine input_tracks
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine announce_read_trackers(print_type)
  use amr_commons,only: myid, aexp
  use tracker_commons
  implicit none
  integer,intent(in)::print_type
  ! Internal dummies
  integer::it
  logical::out_xtrack,out_vtrack,out_quanttrack,out_idtrack,out_internal,out_scale_fac
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

  ! Guard clause to allow only prints from main processor
  if (myid.ne.1) return

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  out_xtrack=.false.; out_vtrack=.false.  
  out_quanttrack=.false.; out_idtrack=.false.; out_internal=.false.  
  if ((print_type.eq.print_alltrack).or.(print_type.eq.print_xtrack)) out_xtrack=.true.
  if ((print_type.eq.print_alltrack).or.(print_type.eq.print_vtrack)) out_vtrack=.true.
  if ((print_type.eq.print_alltrack).or.(print_type.eq.print_measuretrack)) out_quanttrack=.true.
  if ((print_type.eq.print_alltrack).or.(print_type.eq.print_measuretrack)) out_scale_fac=.true.
  if (print_type.eq.print_idtrack) out_idtrack=.true.
  if (print_type.eq.print_internaltrack) out_internal=.true.
  
  if (out_scale_fac) then
300  format("a=",F10.7," dfac=",E14.7," lfac=",E14.7," tfac=",E14.7)
     write(*,300) aexp, scale_d, scale_l, scale_t     
  end if
  
  do it=1,ntracker
301 format("Track center of mass (it=",I3," np=",I6,")=",3F10.7)
302 format("Track COM vel    (it=",I3,")=",3(E14.7," "))
     if (out_xtrack) write(*,301) it, tracker(it)%ntracked, tracker(it)%xpos
     if (out_vtrack) write(*,302) it, tracker(it)%vel
     if (track_quantities.and. out_quanttrack) then
303 format("Track region     (it=",I3,")=",3(E14.7," "))
304 format("Track ",A5,"      (it=",I3,")=",3(E14.7," "))
308 format("Track mmetal     (it=",I3,")=",3(E14.7," "))
305 format("Track ",A10," (it=",I3,")=",3(E14.7," "))
307 format("Track part flow  (it=",I3,")=",4(E14.7," "))
401 format("Track L",A8,"  (it=",I3,")=",3(E14.7," "))
402 format("Track sL",A8," (it=",I3,")=",3(E14.7," "))
        ! Particle prints
        write(*,303) it, tracker(it)%rtrack, tracker(it)%rmeasure, tracker(it)%delta_bound_code
        write(*,304) "mdm  ", it,tracker(it)%mdm,tracker(it)%mdm_inflow,tracker(it)%mdm_outflow
        write(*,304) "mstar", it,tracker(it)%mstar,tracker(it)%mstar_inflow,tracker(it)%mstar_outflow
        if (tracker(it)%part_angular_specific) then
           write(*,402) "dm      ",it,tracker(it)%Ldm
           write(*,402) "dm_in   ",it,tracker(it)%Ldm_in 
           write(*,402) "dm_out  ",it,tracker(it)%Ldm_out
           write(*,402) "star    ",it,tracker(it)%Lstar
           write(*,402) "star_in ",it,tracker(it)%Lstar_in
           write(*,402) "star_out",it,tracker(it)%Lstar_out
        else
           write(*,401) "dm      ",it,tracker(it)%Ldm
           write(*,401) "dm_in   ",it,tracker(it)%Ldm_in 
           write(*,401) "dm_out  ",it,tracker(it)%Ldm_out
           write(*,401) "star    ",it,tracker(it)%Lstar
           write(*,401) "star_in ",it,tracker(it)%Lstar_in
           write(*,401) "star_out",it,tracker(it)%Lstar_out
        end if

        ! Gas prints
        write(*,304) "mgas ", it,tracker(it)%mgas,tracker(it)%mgas_inflow,tracker(it)%mgas_outflow
        write(*,402) "gas     ",it,tracker(it)%Lgas
        write(*,402) "gas_in  ",it,tracker(it)%Lgas_in
        write(*,402) "gas_out ",it,tracker(it)%Lgas_out
        write(*,308) it, tracker(it)%mmetal, tracker(it)%mmetal_inflow, tracker(it)%mmetal_outflow
        write(*,305) "ethermal  ",it,tracker(it)%ethermal,tracker(it)%ethermal_in,tracker(it)%ethermal_out
        write(*,305) "ekin_bulk ",it,tracker(it)%ekin_bulk,tracker(it)%ekin_bulk_in,tracker(it)%ekin_bulk_out
        write(*,305) "ekin      ",it,tracker(it)%ekin,tracker(it)%ekin_in,tracker(it)%ekin_out
        write(*,305) "emag      ",it,tracker(it)%emag,tracker(it)%emag_in,tracker(it)%emag_out
        write(*,305) "ecrs      ",it,tracker(it)%ecrs,tracker(it)%ecrs_in,tracker(it)%ecrs_out

     end if
     if (out_idtrack) then
507 format("Tracked particles for tracker ",I3)
        write(*,507) it
        write(*,*) tracker(it)%id_track     
     end if
     if (out_internal) then
508 format("Tracker status active (it=",I3,")=",L2," measure Q=",L2)
509 format("Tracker status merged (it=",I3,")=",L2," has sink=",L2," seed at a=",E14.7)
510 format("Tracker nmax parts    (it=",I3,")=",I4," active from a=",E14.7)
511 format("Tracker track radius  (it=",I3,")=",E14.7," (",E14.7," squared) ")
512 format("Tracker trackmax dist (it=",I3,")=",E14.7," with tracked mass=",E14.7)
513 format("Tracker measure reg   (it=",I3,")=",E14.7," physical / ",E14.7," comoving")
514 format("Tracker int. rad_meas (it=",I3,")=",E14.7," (",E14.7," squared) ")
        write(*,508) it, tracker(it)%active, tracker(it)%quantities
        write(*,509) it, tracker(it)%merged, tracker(it)%has_sink, tracker(it)%aexp_seed_sink
        write(*,510) it, tracker(it)%ntrack_use, tracker(it)%aexp_ini
        write(*,511) it, tracker(it)%rmax_track, tracker(it)%rmax_track2
        write(*,512) it, tracker(it)%rtrack, tracker(it)%track_mass
        write(*,513) it, tracker(it)%rad_measure
        write(*,514) it, tracker(it)%rmeasure, tracker(it)%rmeasure2
     end if
  end do


  
end subroutine announce_read_trackers
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine input_track_old(filename,myid)
  ! This routine serves for backward compatibility and should eventually
  ! be retired!!
  use pm_commons
  use tracker_commons
  implicit none
  character(len=80)::filename,fileloc
  integer::ilun
  integer,intent(in)::myid
  ilun=12
  fileloc=trim(filename)
  open(unit=ilun,file=fileloc,form='formatted')
  read(ilun,'("ntrack      =",I)') tracker(1)%ntracked
  read(ilun,'("xtrack      =",F10.8)') tracker(1)%xpos(1)
  read(ilun,'("ytrack      =",F10.8)') tracker(1)%xpos(2)
  read(ilun,'("ztrack      =",F10.8)') tracker(1)%xpos(3)

  read(ilun,'("")') 
  read(ilun,'("Tracked particle IDs")') 
  read(ilun,*) tracker(1)%id_track(1:tracker(1)%ntracked)
  close(ilun)  
  track_just_read=.true.
  tracker(1)%xpos_def=tracker(1)%xpos
  if (tracker(1)%ntracked.gt.0) tracker(1)%active=.true.
  if (track_quantities.and.(tracker(1)%ntracked.gt.0)) tracker(1)%quantities=.true.
end subroutine input_track_old
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################
subroutine seed_smbh_center_mass(it)
  use amr_commons
  use pm_commons
  use tracker_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,intent(in)::it
  integer::itype,idim,info
  integer::isink,n_newsink,n_newsink_all
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::scale_m
  real(dp)::msun2gram=2d33

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3d0

  ! Set new sink variables to zero 
  msink_new=0d0; mseed_new=0d0; tsink_new=0d0; delta_mass_new=0d0;
  xsink_new=0d0; vsink_new=0d0
  oksink_new=0d0; idsink_new=0; new_born_new=.false.
  n_newsink=0
  
  SeedRequired: if ((nsink.lt.nsink_com).and.(nsink.lt.nsinkmax).and.(myid.eq.1)) then
     
     if (drag_part) then
        v_DFnew=0d0; mass_DFnew=0d0; n_partnew=0d0
        mass_lowspeednew=0d0; fact_fastnew=0d0
     end if

     write(*,*) "New SMBH (m,x,v):",smbh_mass_track_mass, tracker(it)%xpos, tracker(it)%vel
     ! nsink=nsink+1
     nindsink=nindsink+1
     n_newsink=n_newsink+1
     isink=nindsink
     idsink_new(isink)=nindsink
     oksink_new(isink)=1
     !msink_new(isink)=(0.5d0*smbh_mass_track_mass/scale_m)*msun2gram
     msink_new(isink)=(smbh_mass_track_mass/scale_m)*msun2gram
     tracker(it)%has_sink=.true.
     xsink_new(isink,:)=tracker(it)%xpos
     ! xsink_new(isink,1)=xtrack_center_mass(1)
     ! xsink_new(isink,2)=xtrack_center_mass(2)
     ! xsink_new(isink,3)=xtrack_center_mass(3)
     vsink_new(isink,:)=tracker(it)%vel
     !write(*,*) "line1",vsink_new(isink,:) 
     ! vsink_new(isink,1)=vtrack_center_mass(1)
     ! vsink_new(isink,2)=vtrack_center_mass(2)
     ! vsink_new(isink,3)=vtrack_center_mass(3)
     lsink_new(isink,:)=0.0d0
     ! lsink_new(isink,1)=0.0d0
     ! lsink_new(isink,2)=0.0d0
     ! lsink_new(isink,3)=0.0d0
     tsink_new(isink)=t
     new_born_new(isink)=.true.
     mseed_new(isink)=(smbh_mass_track_mass/scale_m)*msun2gram
     if (drag_part) then
        do itype = 1, 2
           v_DFnew(isink, 1, itype) = vsink_new(isink, 1)
           v_DFnew(isink, 2, itype) = vsink_new(isink, 2)
           v_DFnew(isink, 3, itype) = vsink_new(isink, 3)
           !write(*,*) "line2",vsink_new(isink,1), vsink_new(isink, 2), vsink_new(isink, 3) 
           mass_DFnew(isink, itype) = tiny(0d0)
           mass_lowspeednew(isink, itype) = 0d0
           fact_fastnew(isink, itype) = 0d0
           n_partnew(isink, itype) = 0
        end do
     end if
  end if SeedRequired
  
  ! Communicate new number of sinks
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(n_newsink,n_newsink_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
  nsink=nsink+n_newsink_all
  
  ! Communicate new sink to all cores
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(oksink_new,oksink_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(idsink_new,idsink_all,nsinkmax,MPI_INTEGER         ,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(msink_new ,msink_all ,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mseed_new ,mseed_all ,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(tsink_new ,tsink_all ,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  !write(*,'("BMPI Track:",3(1X,1PE14.7))'), xsink_new
  call MPI_ALLREDUCE(xsink_new ,xsink_all ,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  !write(*,'("AMPI Track:",3(1X,1PE14.7))'), xsink_all
  call MPI_ALLREDUCE(vsink_new ,vsink_all ,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(delta_mass_new,delta_mass_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(new_born_new,new_born_all,nsinkmax,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,info)
  !write(*,*) "track MPI",new_born_new,new_born_all,nsinkmax,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,info 
  ! Particles dynamical friction
  if (drag_part) then
     call MPI_ALLREDUCE(v_DFnew,v_DFall,nsinkmax*ndim*2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     !write(*,*) "line30",v_DFnew,v_DFall 
     call MPI_ALLREDUCE(mass_DFnew, mass_DFall,nsinkmax*2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(n_partnew, n_partall,nsinkmax*2,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(mass_lowspeednew, mass_lowspeedall,nsinkmax*2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(fact_fastnew, fact_fastall,nsinkmax*2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  end if
#else
  oksink_all=oksink_new
  idsink_all=idsink_new
  msink_all=msink_new
  mseed_all=mseed_new
  tsink_all=tsink_new
  xsink_all=xsink_new
  vsink_all=vsink_new
  !write(*,*) "line3",vsink_new,vsink_all 
  delta_mass_all=delta_mass_new
  new_born_all=new_born_new
  ! Particles dynamical friction
  if (drag_part) then
     v_DFall(1:nsink, 1:ndim, 1:2) = v_DFnew(1:nsink, 1:ndim, 1:2)
     !write(*,*)"line31",v_DFall(1:nsink, 1:ndim, 1:2)
     mass_DFall(1:nsink, 1:2) = mass_DFnew(1:nsink, 1:2)
     n_partall(1:nsink, 1:2) = n_partnew(1:nsink, 1:2)
     mass_lowspeedall(1:nsink, 1:2) = mass_lowspeednew(1:nsink, 1:2)
     fact_fastall(1:nsink, 1:2) = fact_fastnew(1:nsink, 1:2)
  end if
#endif
  
  do isink=1,nsink
     if(oksink_all(isink)==1)then
        idsink(isink)=idsink_all(isink)
        msink(isink)=msink_all(isink)
        mseed(isink)=mseed_all(isink)
        tsink(isink)=tsink_all(isink)
        xsink(isink,1:ndim)=xsink_all(isink,1:ndim)
        vsink(isink,1:ndim)=vsink_all(isink,1:ndim)
        !write(*,*) "line4",vsink(isink,1:ndim) 
        delta_mass(isink)=delta_mass_all(isink)
        acc_rate(isink)=msink_all(isink)
        new_born(isink)=new_born_all(isink)
        ! Particles dynamical friction 
        if (drag_part) then
           do itype = 1,2
              mass_DF(isink, levelmin, itype) = mass_DFall(isink, itype)
              n_part(isink, levelmin, itype) = n_partall(isink, itype)
              mass_lowspeed(isink, levelmin, itype) = mass_lowspeedall(isink, itype)
              fact_fast(isink, levelmin, itype) = fact_fastall(isink, itype)
              do idim = 1, ndim
                 v_DF(isink,levelmin, idim, itype) = v_DFall(isink, idim, itype)
                 !write(*,*)"line41",v_DF(isink,levelmin, idim, itype)
              end do
           end do
        end if
     end if
  end do
  
end subroutine seed_smbh_center_mass
!#####################################################################
!#####################################################################
!#####################################################################
!#####################################################################


