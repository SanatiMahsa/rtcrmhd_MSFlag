! Patch changes:
! - added mp0 variable for particles
module pm_commons
  use amr_parameters
  use pm_parameters
  use random
  ! Sink particle related arrays
  real(dp),allocatable,dimension(:)::msink,c2sink,oksink_new,oksink_all
  real(dp),allocatable,dimension(:)::r2sink,rinfsink
  real(dp),allocatable,dimension(:)::tsink,tsink_new,tsink_all
  real(dp),allocatable,dimension(:)::msink_new,msink_all
  real(dp),allocatable,dimension(:)::mseed,mseed_new,mseed_all
  real(dp),allocatable,dimension(:)::xmsink
  real(dp),allocatable,dimension(:)::dMsink_overdt,dMBHoverdt
  real(dp),allocatable,dimension(:)::rho_gas,volume_gas,eps_sink
  real(dp),allocatable,dimension(:)::eth_sink,emag_sink,metal_sink
  real(dp),allocatable,dimension(:)::agn_inj_mag,agn_inj_tot
  real(dp),allocatable,dimension(:,:)::vel_gas
  real(dp),allocatable,dimension(:)::delta_mass,delta_mass_new,delta_mass_all
  real(dp),allocatable,dimension(:)::wden,weth,wvol,wdiv,wden_new,weth_new,wvol_new,wdiv_new
  real(dp),allocatable,dimension(:)::wemag,wemag_new
  real(dp),allocatable,dimension(:)::wmetal,wmetal_new
  real(dp),allocatable,dimension(:,:)::lmom_gas,lmom_gas_new,lmom_gas_all
  real(dp),allocatable,dimension(:,:)::wmom,wmom_new
  real(dp),allocatable,dimension(:,:)::vsink,vsink_new,vsink_all
  real(dp),allocatable,dimension(:,:)::fsink,fsink_new,fsink_all
  real(dp),allocatable,dimension(:,:,:)::vsnew,vsold
  real(dp),allocatable,dimension(:,:,:)::fsink_partial,sink_jump
  real(dp),allocatable,dimension(:,:)::lsink,lsink_new,lsink_all!sink angular momentum
  real(dp),allocatable,dimension(:,:)::xsink,xsink_new,xsink_all
  real(dp),allocatable,dimension(:)::acc_rate,acc_lum !sink accretion rate and luminosity
  real(dp),allocatable,dimension(:,:)::weighted_density,weighted_volume,weighted_ethermal,weighted_divergence
  real(dp),allocatable,dimension(:,:)::weighted_emag,weighted_metal
  real(dp),allocatable,dimension(:,:,:)::weighted_momentum
  real(dp),allocatable,dimension(:)::dt_acc                ! maximum timestep allowed by the sink
  real(dp),allocatable,dimension(:)::rho_sink_tff
  integer,allocatable,dimension(:)::idsink,idsink_new,idsink_old,idsink_all
  logical,allocatable,dimension(:,:)::level_sink,level_sink_new
  logical,allocatable,dimension(:)::ok_blast_agn,ok_blast_agn_all,direct_force_sink
  logical,allocatable,dimension(:)::new_born,new_born_all,new_born_new
  integer,allocatable,dimension(:)::idsink_sort
  integer::ncloud_sink,ncloud_sink_massive
  integer::nindsink=0
  integer::sinkint_level=0         ! maximum level currently active is where the global sink variables are updated
  real(dp)::ssoft                  ! sink softening lenght in code units
  real(dp)::rho_min_acc_code


  ! Particles drag on SMBH
  real(dp),allocatable,dimension(:,:,:,:)::v_DFnew_all
  real(dp),allocatable,dimension(:,:,:,:)::v_DF
  real(dp),allocatable,dimension(:,:,:)::v_DFnew,v_DFall
  real(dp),allocatable,dimension(:,:,:)::vrel_sink
  real(dp),allocatable,dimension(:,:,:)::mass_DF,mass_DFnew_all
  real(dp),allocatable,dimension(:,:)::mass_DFnew,mass_DFall
  real(dp),allocatable,dimension(:,:,:)::fact_fast,fact_fastnew_all
  real(dp),allocatable,dimension(:,:)::fact_fastnew,fact_fastall
  real(dp),allocatable,dimension(:,:)::vrel_sink_norm
  real(dp),allocatable,dimension(:,:,:)::mass_lowspeed,mass_lowspeednew_all
  real(dp),allocatable,dimension(:,:)::mass_lowspeednew,mass_lowspeedall
  integer,allocatable,dimension(:,:,:)::n_part,n_partnew_all
  integer,allocatable,dimension(:,:)::n_partnew,n_partall
  integer,allocatable,dimension(:)::sink_cell
  real(dp),allocatable,dimension(:,:,:)::v_background
  real(dp),allocatable,dimension(:,:)  ::fact_fast_background
  real(dp),allocatable,dimension(:,:)  ::n_background
  real(dp),allocatable,dimension(:,:)  ::mass_lowspeed_background
  real(dp),allocatable,dimension(:,:)  ::m_background
  real(dp),allocatable,dimension(:)  ::most_massive_sink
  integer::DF_ncells=5
  
  ! Particles related arrays
  real(dp),allocatable,dimension(:,:)::xp       ! Positions
  real(dp),allocatable,dimension(:,:)::vp       ! Velocities
  real(dp),allocatable,dimension(:)  ::mp       ! Masses
  real(dp),allocatable,dimension(:)  ::mp0      ! Initial masses (for Kimm feedback)
#ifdef OUTPUT_PARTICLE_POTENTIAL
  real(dp),allocatable,dimension(:)  ::ptcl_phi ! Potential of particle added by AP for output purposes 
#endif
#ifdef NTRACEGROUPS
  integer ,allocatable,dimension(:)  ::ptracegroup
#endif
  real(dp),allocatable,dimension(:)  ::tp       ! Birth epoch
  real(dp),allocatable,dimension(:,:)::weightp  ! weight of cloud parts for sink accretion only
  real(dp),allocatable,dimension(:)  ::zp       ! Birth metallicity
  integer ,allocatable,dimension(:)  ::nextp    ! Next particle in list
  integer ,allocatable,dimension(:)  ::prevp    ! Previous particle in list
  integer ,allocatable,dimension(:)  ::levelp   ! Current level of particle
  integer(i8b),allocatable,dimension(:)::idp    ! Identity of particle
  real(dp),allocatable,dimension(:)  ::st_n_tp  ! Gas density at birth epoch         !SD
  real(dp),allocatable,dimension(:)  ::st_n_sn  ! Gas density at SN epoch            !SD
  real(dp),allocatable,dimension(:)  ::st_e_sn  ! SN energy injected                 !SD

  ! Tree related arrays
  integer ,allocatable,dimension(:)  ::headp    ! Head particle in grid
  integer ,allocatable,dimension(:)  ::tailp    ! Tail particle in grid
  integer ,allocatable,dimension(:)  ::numbp    ! Number of particles in grid
  ! Global particle linked lists
  integer::headp_free,tailp_free,numbp_free=0,numbp_free_tot=0
  ! Local and current seed for random number generator
  integer,dimension(IRandNumSize) :: localseed=-1

  ! xtrack related arrays
  ! logical          ::track_just_read=.false.
  ! integer,parameter::ntrack_com_max=500
  ! integer          ::ntrack_com_number=ntrack_com_max
  ! integer          ::ntrack_center_mass=-1
  !real(dp)         ::max_dist_track=0.15d0
  !integer(i8b)     ::id_track_center_mass(ntrack_com_max)
  !integer(i8b)     ::unmatched(ntrack_com_max)
  !real(dp)         ::xtrack_center_mass(3)=0.0d0
  !real(dp)         ::vtrack_center_mass(3)=0.0d0
  !real(dp)         ::mtrack_center_mass
  !real(dp)         ::new_center_mass(3),mass_new_center_mass
  !real(dp)         ::new_vcom_mass(3)

  contains
  function cross(a,b)
    use amr_parameters, only:dp
    real(dp),dimension(1:3)::a,b
    real(dp),dimension(1:3)::cross
    !computes the cross product c= a x b
    cross(1)=a(2)*b(3)-a(3)*b(2)
    cross(2)=a(3)*b(1)-a(1)*b(3)
    cross(3)=a(1)*b(2)-a(2)*b(1)
  end function cross

  function is_cloud_sink(id_part,m_part,tf_part) result(is_cloud_sink_part)
    ! This function receives a particle ID, mass, and formation time
    ! and returns a logical = .true. if the particle 
    ! is a sink cloud particle 
      use amr_parameters, only:dp
      integer,intent(in)::id_part
      real(dp),intent(in)::tf_part, m_part
      logical::is_cloud_sink_part
      if((id_part.lt.0) .and. (tf_part.eq.0.0d0) .and. (abs(id_part).le.nsinkmax_safe))then
        is_cloud_sink_part=.true.
      else
        is_cloud_sink_part=.false.
      end if
    end function is_cloud_sink
  
    function is_star(id_part,m_part,tf_part) result(is_star_part)
      ! This function receives a particle ID, mass and formation time
      ! and returns a logical = .true. if the particle 
      ! is a sink cloud particle 
        use amr_parameters, only:dp
        integer,intent(in)::id_part
        real(dp),intent(in)::tf_part, m_part
        logical::is_star_part
        if((tf_part.ne.0.0d0) .and. (m_part.ne.0.0d0) .and. (abs(id_part).gt.nsinkmax_safe))then
          is_star_part=.true.
        else
          is_star_part=.false.
        end if
      end function is_star
  
      function is_active_star(id_part,m_part,tf_part) result(is_active_star_part)
        ! This function receives a particle ID, mass and formation time
        ! and returns a logical = .true. if the particle 
        ! is a sink cloud particle 
          use amr_parameters, only:dp
          integer,intent(in)::id_part
          real(dp),intent(in)::tf_part, m_part
          logical::is_active_star_part
          is_active_star_part=.false.
          if(is_star(id_part,m_part,tf_part)) then 
            if (id_part.lt.0)then
              is_active_star_part=.true.
            end if
          end if
        end function is_active_star

      function is_dm(id_part,m_part,tf_part) result(is_dm_part)
    ! This function receives a particle ID, mass, and formation time
    ! and returns a logical = .true. if the particle 
    ! is a dark matter particle 
      use amr_parameters, only:dp
      integer,intent(in)::id_part
      real(dp),intent(in)::tf_part, m_part
      logical::is_dm_part
      if((id_part.gt.0) .and. (tf_part.eq.0.0d0))then
        is_dm_part=.true.
      else
        is_dm_part=.false.
      end if
    end function is_dm  
  
end module pm_commons
