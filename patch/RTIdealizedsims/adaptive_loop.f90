subroutine adaptive_loop
  use amr_commons
  use hydro_commons
  use pm_commons
  use poisson_commons
  use cooling_module
#ifdef DICE
  use dice_commons, ONLY:dice_init, dice_mfr, dmf_i, dmf_loop, dmf_ordering
  use amr_parameters, ONLY:ordering      
#endif
#ifdef RT
  use rt_hydro_commons
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer(kind=8)::n_step
  integer::ilevel,idim,ivar,info,tot_pt
  real(kind=8)::tt1,tt2,muspt,muspt_this_step,maxsec,outsec
  real(kind=4)::real_mem,real_mem_tot
  real(kind=8),save::tshar=0.0

#ifndef WITHOUTMPI
  tt1=MPI_WTIME()
  !For output just before runtime end
  if (tshar.eq.0.0) then
     tshar = MPI_WTIME(info)
  end if
#endif

#ifdef DICE
  if (nrestart.eq.0 .and. dice_mfr) then
    call dice_multifile_read
  else
#endif
  call init_amr                      ! Initialize AMR variables
  call init_time                     ! Initialize time variables
  if(hydro)call init_hydro           ! Initialize hydro variables
#ifdef RT
  if(rt.or.neq_chem) &
       & call rt_init_hydro          ! Initialize radiation variables
#endif
  if(poisson)call init_poisson       ! Initialize poisson variables
#ifdef ATON
  if(aton)call init_radiation        ! Initialize radiation variables
#endif
  if(mechanical_feedback) call init_mechanical
  if(nrestart==0)call init_refine    ! Build initial AMR grid

#ifdef grackle
  if(cosmo)then
     ! Compute cooling table at current aexp
  endif
#else  
  if(cooling.and..not.neq_chem) &
       call set_table(dble(aexp))    ! Initialize cooling look up table
#endif
  if(pic)call init_part              ! Initialize particle variables
  if(pic)call init_tree              ! Initialize particle tree
  if(nrestart==0)call init_refine_2  ! Build initial AMR grid again
#ifdef DICE
  end if
  dice_mfr = .FALSE.
#endif



#ifndef WITHOUTMPI
  muspt=0.
  tot_pt=-1
  tt2=MPI_WTIME()
  if(myid==1)write(*,*)'Time elapsed since startup:',tt2-tt1
#endif

  if(myid==1)then
     write(*,*)'Initial mesh structure'
     do ilevel=1,nlevelmax
        if(numbtot(1,ilevel)>0)write(*,999)ilevel,numbtot(1:4,ilevel)
     end do
  end if

  nstep_coarse_old=nstep_coarse

  if(myid==1)write(*,*)'Starting time integration' 
#ifdef DICE
  dice_init = .false.
#endif

  do ! Main time loop
                               call timer('coarse levels','start')

#ifndef WITHOUTMPI
     tt1=MPI_WTIME()
#endif

     if(verbose)write(*,*)'Entering amr_step_coarse'

     epot_tot=0.0D0  ! Reset total potential energy
     ekin_tot=0.0D0  ! Reset total kinetic energy
     mass_tot=0.0D0  ! Reset total mass
     eint_tot=0.0D0  ! Reset total internal energy
#ifdef SOLVERmhd
     emag_tot=0.0D0  ! Reset total magnetic energy     
#endif
#if NCR>0
     ecrs_tot=0.0D0  ! Reset total cosmic rays energy
#endif    

     ! Make new refinements
     if(levelmin.lt.nlevelmax.and.(.not.static.or.(nstep_coarse_old.eq.nstep_coarse.and.restart_remap)))then
        call refine_coarse
        do ilevel=1,levelmin
           call build_comm(ilevel)
           call make_virtual_fine_int(cpu_map(1),ilevel)
           if(hydro)then
#ifdef SOLVERmhd
              do ivar=1,nvar+3
#else
              do ivar=1,nvar
#endif
                 call make_virtual_fine_dp(uold(1,ivar),ilevel)
#ifdef SOLVERmhd
              end do
#else
              end do
#endif
              if(simple_boundary)call make_boundary_hydro(ilevel)
           endif
#ifdef RT
           if(rt)then
              do ivar=1,nrtvartrace
                 call make_virtual_fine_rtdp(rtuold(1,ivar),ilevel)
              end do
              if(simple_boundary)call rt_make_boundary_hydro(ilevel)
           endif
#endif
           if(poisson)then
              call make_virtual_fine_dp(phi(1),ilevel)
              do idim=1,ndim
                 call make_virtual_fine_dp(f(1,idim),ilevel)
              end do
           end if
           if(ilevel<levelmin)call refine_fine(ilevel)
        end do
     endif

     ! Call base level
     call amr_step(levelmin,1)
                               call timer('coarse levels','start')

     if(levelmin.lt.nlevelmax.and.(.not.static.or.(nstep_coarse_old.eq.nstep_coarse.and.restart_remap)))then
        do ilevel=levelmin-1,1,-1
           ! Hydro book-keeping
           if(hydro)then
              call upload_fine(ilevel)
#ifdef SOLVERmhd
              do ivar=1,nvar+3
#else
              do ivar=1,nvar
#endif
                 call make_virtual_fine_dp(uold(1,ivar),ilevel)
#ifdef SOLVERmhd
              end do
#else
              end do
#endif
              if(simple_boundary)call make_boundary_hydro(ilevel)
           end if
#ifdef RT
           ! Radiation book-keeping
           if(rt)then
              call rt_upload_fine(ilevel)
              do ivar=1,nrtvartrace
                 call make_virtual_fine_rtdp(rtuold(1,ivar),ilevel)
              end do
              if(simple_boundary)call rt_make_boundary_hydro(ilevel)
           end if
#endif
           ! Gravity book-keeping
           if(poisson)then
              call make_virtual_fine_dp(phi(1),ilevel)
              do idim=1,ndim
                 call make_virtual_fine_dp(f(1,idim),ilevel)
              end do
           end if
        end do
        
        ! Build refinement map
        do ilevel=levelmin-1,1,-1
           call flag_fine(ilevel,2)
        end do
        call flag_coarse
     endif

     ! New coarse time-step
     nstep_coarse=nstep_coarse+1

#ifndef WITHOUTMPI
     tt2=MPI_WTIME()
     if(mod(nstep_coarse,ncontrol)==0)then
        call getmem(real_mem)
        call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
        if(myid==1)then
           if (tot_pt==0) muspt=0. ! dont count first timestep
           n_step = int(numbtot(1,levelmin),kind=8)*twotondim
           do ilevel=levelmin+1,nlevelmax
             n_step = n_step + int(numbtot(1,ilevel),kind=8)*product(nsubcycle(levelmin:ilevel-1))*(twotondim-1)
           enddo
           muspt_this_step = (tt2-tt1)*1e6/n_step*ncpu
           muspt = muspt + muspt_this_step
           tot_pt = tot_pt + 1
           write(*,'(a,f8.2,a,f12.2,a,f12.2,a)')' Time elapsed since last coarse step:',tt2-tt1 &
          ,' s',muspt_this_step,' mus/pt'  &
          ,muspt / max(tot_pt,1), ' mus/pt (av)'
           call writemem(real_mem_tot)

           write(*,*)'Total running time:', NINT((tt2-tshar)*100.0)*0.01,'s'
        endif
        if(maxruntime.gt.0.d0) then
           maxsec = maxruntime*60.d0*60.d0 !Convert maxruntime from hours to seconds
           outsec = minutes_save*60.d0     !Convert minutes before run time end to seconds
           if(maxsec-outsec.lt.tt2-tshar) then
              output_now=.true.
              if(myid==1) write(*,*) 'Printing snapshot before runtime end!!!'
              !now set maxruntime to a negative number so that we dont keep printing outputs
              maxruntime = -1.0d0
           endif
        endif
     endif
#endif

  end do

999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

end subroutine adaptive_loop


!############################################################################################
!############################################################################################
!############################################################################################
!############################################################################################
!############################################################################################
!############################################################################################
!############################################################################################
!############################################################################################
!############################################################################################
!############################################################################################
!############################################################################################


#ifdef DICE
SUBROUTINE dice_multifile_read

  use amr_commons
  use hydro_commons
  use pm_commons
  use poisson_commons
  use cooling_module
  use dice_commons   
  use amr_parameters, ONLY:ordering  
#ifdef RT
  use rt_hydro_commons
#endif

  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  INTEGER::info,iii
  ! integer ,allocatable,dimension(:)  ::son_dmf
  ! integer ,allocatable,dimension(:)  ::father_dmf

  !<<<<<----- Set up Initial Grids and variables ----->>>>>!
  ! dmf_ordering = ordering
  ! ordering='plana'
  call init_amr                      ! Initialize AMR variables
  call init_time                     ! Initialize time variables
  if(hydro)call init_hydro           ! Initialize hydro variables
#ifdef RT
  if(rt.or.neq_chem) &
    & call rt_init_hydro               ! Initialize radiation variables
#endif
  if(poisson)call init_poisson       ! Initialize poisson variables
#ifdef ATON
  if(aton)call init_radiation        ! Initialize radiation variables
#endif
  if(mechanical_feedback) call init_mechanical
  if(nrestart==0)call init_refine    ! Build initial AMR grid


#ifdef grackle
  if(cosmo)then
    ! Compute cooling table at current aexp
  endif
#else  
  if(cooling.and..not.neq_chem) &
    call set_table(dble(aexp))    ! Initialize cooling look up table
#endif
  !<<<<<----- ---------------------------------- ----->>>>>!

  allocate(readflag(1:dmf_ncell))
  if (dmf_loop.gt.1) then
    allocate(utmp(1:dmf_ncell,1:nvar))
    utmp = 0.0d0
  end if
  readflag=0
! #endif/
  !<<<<<----- LOAD DICE GAS PARTICLES ----->>>>>!

  do dmf_i=1, dmf_loop  
    ! if (dmf_i.eq.3) cycle

    call load_dice_particles
    if (hydro.eqv..false.) call kill_gas_part(1)
    dice_init=.false.
    ! call init_tree
    
    if (dmf_i.gt.1) then
      readflag = 0
      do iii=1, dmf_ncell
        if (uold(iii,1).gt.0.0d0) then
          readflag(iii) = 1
          utmp(iii,:) = uold(iii,:)
        end if
      end do
    end if

    call init_refine_2
    ! deallocate(up)

#ifndef WITHOUTMPI
    call MPI_Barrier(MPI_COMM_WORLD,info)
#endif
    if (dmf_loop.gt.1) then
      if (myid.eq.1) write(*,*) "DICE loop ",dmf_i,' of ',dmf_loop,' completed' 
    else
      if (myid.eq.1) write(*,*) "DICE gas particle read in complete"
    end if
  end do
#ifndef WITHOUTMPI
    call MPI_Barrier(MPI_COMM_WORLD,info)
    ! if (myid.eq.1) print*, 'KEARN BARRIER TEST 2'
#endif
  ! if (myid.eq.1) then
  !   write(*,*) ' Initialising Birth mass of star particles'
  !   write(*,*) ' -----'
  !   write(*,*) ' '
  ! end if

  ! print*, myid,sum(mpb)
  ! allocate(mpb (npartmax))
  ! mpb = 0.0d0
  ! do dmf_i=1, npartmax
  !   if (mp(dmf_i).ne.0.0d0)then
  !     if (tp(dmf_i).lt.0.0d0)then
  !       mpb(dmf_i) = mp(dmf_i)
  !     end if
  !   end if
  ! end do

  ! call load_balance
  dice_init=.false.
#ifndef WITHOUTMPI
    call MPI_Barrier(MPI_COMM_WORLD,info)
    ! if (myid.eq.1) print*, 'KEARN BARRIER TEST 2'
#endif
  if (myid.eq.1) then
    write(*,*) ' Particles Loaded into AMR grid'
    write(*,*) ' -----'
    write(*,*) ' '
  end if

  if (hydro.eqv..false.) then
    T2_star = 0.0d0
    if (myid.eq.1) then
      write(*,*) ''
      write(*,*) ' -----'
      write(*,*) 'Dark Matter Only Run'
      write(*,*) 'No Gas Particles loaded'
      write(*,*) ' -----'
      write(*,*) ' '
    end if
    call MPI_Barrier(MPI_COMM_WORLD,info)
  end if
  ! utmp = 0.0d0
  ! utmp = uold
  !<<<<<----- ---------------------------------- ----->>>>>!



!   !<<<<<----- LOAD DICE STAR AND DM PARTICLES ----->>>>>!
!   call load_dice_particles      
! #ifndef WITHOUTMPI
!   call MPI_Barrier(MPI_COMM_WORLD,info)

! #endif
!   if (myid.eq.1) then
!     write(*,*) ' DICE star and DM Particles Successfully Loaded '
!     write(*,*) ' -----'
!     write(*,*) ' '
!   end if

! #ifndef WITHOUTMPI
!   call MPI_Barrier(MPI_COMM_WORLD,info)
!   ! if (myid.eq.1) print*, 'KEARN BARRIER TEST 2'
! #endif
!   if (myid.eq.1) then
!     write(*,*) ' Adding Stars and DM to AMR grid '
!   end if
!   call particles_to_grid !(levelmin,1)
!   ! call init_refine_3
! #ifndef WITHOUTMPI
!   call MPI_Barrier(MPI_COMM_WORLD,info)
!   ! if (myid.eq.1) print*, 'KEARN BARRIER TEST 2'
! #endif
!   if (myid.eq.1) then
!     write(*,*) ' Particles successfully added to the grid '
!     write(*,*) ' -----'
!     write(*,*) ' '
!   end if
!   !<<<<<----- ------------------------------- ----->>>>>!


  if (dmf_loop.gt.1) deallocate(utmp)
  deallocate(readflag)
  ! ordering =   dmf_ordering
  ! ordering = 'plana'
#ifndef WITHOUTMPI
  call MPI_Barrier(MPI_COMM_WORLD,info)
  ! if (myid.eq.1) print*, 'KEARN BARRIER TEST 2'
#endif
  if (myid.eq.1) then
    write(*,*) ''
    write(*,*) ' DICE IC readin COMPLETE'
    write(*,*) '__________________________________________________'
    write(*,*) ''
  end if
  

END SUBROUTINE dice_multifile_read




















!############################################################################################
!############################################################################################
!############################################################################################
!############################################################################################
!############################################################################################
!############################################################################################
!############################################################################################
!############################################################################################
!############################################################################################
!############################################################################################
!############################################################################################


SUBROUTINE particles_to_grid  !(ilevel,icount)

  use amr_commons
  use hydro_commons
  use pm_commons
  use poisson_commons
  use cooling_module
  use dice_commons, ONLY:dice_init, dice_mfr, dmf_i, dmf_loop, dmf_ordering
  use amr_parameters, ONLY:ordering      
#ifdef RT
  use rt_hydro_commons
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,idim,ivar,info,control_test,i,icount
  real(kind=8)::tt1,tt2
  real(kind=4)::real_mem,real_mem_tot

  icount = 1
  do ilevel=levelmin, nlevelmax
    !--------------------
    ! Poisson source term
    !--------------------
    if(pic) call make_tree_fine(ilevel)
    if(poisson)call rho_fine(ilevel,icount)

    !-------------------------------------------
    ! Sort particles between ilevel and ilevel+1
    !-------------------------------------------
    if(pic)then
      ! Remove particles to finer levels
      call kill_tree_fine(ilevel)
      ! Update boundary conditions for remaining particles
      call virtual_tree_fine(ilevel)
    end if
  end do


  !-------------------------------------------
  ! Make new refinements and update boundaries
  !-------------------------------------------
  if(levelmin.lt.nlevelmax .and..not. static)then
    ilevel = levelmin
    print*, 'check check hello', ilevel, levelmin, nlevelmax
    if(ilevel==levelmin)then
      do i=ilevel,nlevelmax
        if(i>levelmin)then

          !--------------------------
          ! Build communicators
          !--------------------------
          call build_comm(i)

          !--------------------------
          ! Update boundaries
          !--------------------------
          call make_virtual_fine_int(cpu_map(1),i)
          if(hydro)then
#ifdef SOLVERmhd
            do ivar=1,nvar+3
#else
              do ivar=1,nvar
#endif
                call make_virtual_fine_dp(uold(1,ivar),i)
#ifdef SOLVERmhd
              end do
#else
            end do
#endif
            if(simple_boundary)call make_boundary_hydro(i)
          end if
#ifdef RT
          if(rt)then
            do ivar=1,nrtvar
              call make_virtual_fine_dp(rtuold(1,ivar),i)
            end do
            if(simple_boundary)call rt_make_boundary_hydro(i)
          end if
#endif
          if(poisson)then
            call make_virtual_fine_dp(phi(1),i)
            do idim=1,ndim
              call make_virtual_fine_dp(f(1,idim),i)
            end do
            if(simple_boundary)call make_boundary_force(i)
          end if
        end if

        !--------------------------
        ! Refine grids
        !--------------------------
        call refine_fine(i)
      end do
    end if
  end if


  if(levelmin.lt.nlevelmax .and..not. static)then
    do ilevel=levelmin-1,1,-1
      ! Hydro book-keeping
      if(hydro)then
        call upload_fine(ilevel)
#ifdef SOLVERmhd
        do ivar=1,nvar+3
#else
          do ivar=1,nvar
#endif
            call make_virtual_fine_dp(uold(1,ivar),ilevel)
#ifdef SOLVERmhd
          end do
#else
        end do
#endif
        if(simple_boundary)call make_boundary_hydro(ilevel)
      end if

#ifdef RT      
      ! Radiation book-keeping
      if(rt)then
        call rt_upload_fine(ilevel)
        do ivar=1,nrtvar
          call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
        end do
        if(simple_boundary)call rt_make_boundary_hydro(ilevel)
      end if
#endif
    
      ! Gravity book-keeping
      if(poisson)then     
        call make_virtual_fine_dp(phi(1),ilevel)
        do idim=1,ndim
          call make_virtual_fine_dp(f(1,idim),ilevel)
        end do
      end if
    end do
      
    ! Build refinement map
    do ilevel=levelmin-1,1,-1
      call flag_fine(ilevel,2)
    end do
    call flag_coarse
  endif

  if(levelmin.lt.nlevelmax .and..not.static)then
    call refine_coarse

    do ilevel=1,levelmin
      call build_comm(ilevel)
      call make_virtual_fine_int(cpu_map(1),ilevel)
      if(hydro)then
#ifdef SOLVERmhd
        do ivar=1,nvar+3
#else
          do ivar=1,nvar
#endif
            call make_virtual_fine_dp(uold(1,ivar),ilevel)
#ifdef SOLVERmhd
          end do
#else
        end do
#endif
        if(simple_boundary)call make_boundary_hydro(ilevel)
      endif
#ifdef RT
      if(rt)then
        do ivar=1,nrtvar
          call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
        end do
        if(simple_boundary)call rt_make_boundary_hydro(ilevel)
      endif
#endif
      if(poisson)then
        call make_virtual_fine_dp(phi(1),ilevel)
        do idim=1,ndim
          call make_virtual_fine_dp(f(1,idim),ilevel)
        end do
      end if
      if(ilevel<levelmin)call refine_fine(ilevel)
    end do          
  endif















!   use amr_commons
!   use pm_commons
!   use clfind_commons
!   use hydro_commons
!   use poisson_commons
! #ifdef RT
!   use rt_hydro_commons
!   use SED_module
! #endif
!   implicit none
! #ifndef WITHOUTMPI
!   include 'mpif.h'
! #endif
!   integer::ilevel,icount
!   !-------------------------------------------------------------------!
!   ! This routine is the adaptive-mesh/adaptive-time-step main driver. !
!   ! Each routine is called using a specific order, don't change it,   !
!   ! unless you check all consequences first                           !
!   !-------------------------------------------------------------------!
!   integer::i,idim,ivar
!   logical::ok_defrag
!   logical,save::first_step=.true.

!   if(numbtot(1,ilevel)==0)return



!   !-------------------------------------------
!   ! Make new refinements and update boundaries
!   !-------------------------------------------
!   if(levelmin.lt.nlevelmax .and..not. static)then
!     if(ilevel==levelmin.or.icount>1)then
!       do i=ilevel,nlevelmax
!         if(i>levelmin)then

!           !--------------------------
!           ! Build communicators
!           !--------------------------
!           call build_comm(i)

!           !--------------------------
!           ! Update boundaries
!           !--------------------------
!           call make_virtual_fine_int(cpu_map(1),i)
!           if(hydro)then
! #ifdef SOLVERmhd
!             do ivar=1,nvar+3
! #else
!               do ivar=1,nvar
! #endif
!                 call make_virtual_fine_dp(uold(1,ivar),i)
! #ifdef SOLVERmhd
!               end do
! #else
!             end do
! #endif
!             if(simple_boundary)call make_boundary_hydro(i)
!           end if
! #ifdef RT
!           if(rt)then
!             do ivar=1,nrtvar
!               call make_virtual_fine_dp(rtuold(1,ivar),i)
!             end do
!             if(simple_boundary)call rt_make_boundary_hydro(i)
!           end if
! #endif
!           if(poisson)then
!             call make_virtual_fine_dp(phi(1),i)
!             do idim=1,ndim
!               call make_virtual_fine_dp(f(1,idim),i)
!             end do
!             if(simple_boundary)call make_boundary_force(i)
!           end if
!         end if

!         !--------------------------
!         ! Refine grids
!         !--------------------------
!         call refine_fine(i)
!       end do
!     end if
!   end if

!   ! !--------------------------
!   ! ! Load balance
!   ! !--------------------------
!   ! ok_defrag=.false.
!   ! if(levelmin.lt.nlevelmax)then
!   !   if(ilevel==levelmin)then
!   !     if(nremap>0)then
!   !       ! Skip first load balance because it has been performed before file dump
!   !       if(nrestart>0.and.first_step)then
!   !         first_step=.false.
!   !       else
!   !         if(MOD(nstep_coarse,nremap)==0)then
!   !           call load_balance
!   !           call defrag
!   !           ok_defrag=.true.
!   !         endif
!   !       end if
!   !     end if
!   !   endif
!   ! end if

!   !-----------------
!   ! Update sink cloud particle properties
!   !-----------------
!   if(sink)call update_cloud(ilevel)

!   !-----------------
!   ! Particle leakage
!   !-----------------
!   if(pic)call make_tree_fine(ilevel)



!   if(hydro)then
! #ifdef SOLVERmhd
!     do ivar=1,nvar+3
! #else
!       do ivar=1,nvar
! #endif
!         call make_virtual_fine_dp(uold(1,ivar),ilevel)
! #ifdef SOLVERmhd
!       end do
! #else
!     end do
! #endif
!     if(simple_boundary)call make_boundary_hydro(ilevel)
!   endif
! #ifdef RT
!   if(rt)then
!     do ivar=1,nrtvar
!       call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
!     end do
!     if(simple_boundary)call rt_make_boundary_hydro(ilevel)
!   end if
! #endif

! #ifdef SOLVERmhd
!   ! Magnetic diffusion step
!   if(hydro)then
!     if(eta_mag>0d0.and.ilevel==levelmin)then
!       call diffusion
!     endif
!   end if
! #endif

!   !-----------------------
!   ! Compute refinement map
!   !-----------------------
!   if(.not.static) call flag_fine(ilevel,icount)


!   !----------------------------
!   ! Merge finer level particles
!   !----------------------------
!   if(pic)call merge_tree_fine(ilevel)


END SUBROUTINE particles_to_grid

#endif















!############################################################################################
!############################################################################################
!############################################################################################
!############################################################################################
!############################################################################################
!############################################################################################
!############################################################################################
!############################################################################################
!############################################################################################
!############################################################################################
!############################################################################################
