subroutine courant_fine(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
#ifdef TRACKER
  use tracker_commons
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Using the Courant-Friedrich-Levy stability condition,               !
  ! this routine computes the maximum allowed time-step.                !
  !----------------------------------------------------------------------
  integer::i,ivar,idim,ind,ncache,igrid,iskip
  integer::info,nleaf,ngrid,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf

  real(dp)::dt_lev,dx,vol,scale
  real(kind=8)::mass_loc,ekin_loc,eint_loc,emag_loc,dt_loc
#if NCR>0
  integer,parameter::ncomm=5
  real(kind=8)::ecrs_loc,ecrs_all
#else
  integer,parameter::ncomm=4
#endif

  real(kind=8)::mass_all,ekin_all,eint_all,emag_all,dt_all
  real(kind=8),dimension(ncomm)::comm_buffin,comm_buffout
  real(dp),dimension(1:nvector,1:nvar+3),save::uu
  real(dp),dimension(1:nvector,1:ndim),save::gg


#ifdef TRACKER

#ifdef SOLVERmhd
#if NCR>0
  integer,parameter::ncomm_track=30
#else
  integer,parameter::ncomm_track=27
#endif
#else
  integer,parameter::ncomm_track=24
#endif
  real(kind=8),dimension(ncomm_track)::comm_buffin_track,comm_buffout_track
  real(kind=8),dimension(ntracker_max)::mgas_tr,mgas_inflow_tr,mgas_outflow_tr
  real(kind=8),dimension(ntracker_max,1:ndim)::Lgas_tr,Lgas_inflow_tr,Lgas_outflow_tr
  real(kind=8),dimension(ntracker_max)::mmetal_tr,mmetal_inflow_tr,mmetal_outflow_tr
  real(kind=8),dimension(ntracker_max)::ethermal_tr,ekin_bulk_tr,ekin_tr,emag_tr,ecrs_tr
  real(kind=8),dimension(ntracker_max)::ethermal_in_tr,ekin_bulk_in_tr,ekin_in_tr,emag_in_tr,ecrs_in_tr
  real(kind=8),dimension(ntracker_max)::ethermal_out_tr,ekin_bulk_out_tr,ekin_out_tr,emag_out_tr,ecrs_out_tr
  integer::it,icount
  real(dp)::dd_cell2,vrad_cell
  real(dp),dimension(1:ndim)::xpos_cell,vel_cell,vel_cell_COM,Lcell
  real(dp)::mcell,m_metal,dummy_dp
  integer::ix,iy,iz
  real(dp),dimension(1:twotondim,1:ndim)::xc
  real(dp),dimension(1:ndim)::skip_loc
  real(dp),dimension(1:ndim,1:nvector)::xpos
  logical,dimension(1:ntracker,1:nvector)::do_cell_contrib,do_cell_boundary
#endif


#ifdef TRACKER
  if (track_quantities) then

     mgas_tr=0.0d0; mgas_inflow_tr=0.0d0; mgas_outflow_tr=0.0d0
     Lgas_tr=0.0d0; Lgas_inflow_tr=0.0d0; Lgas_outflow_tr=0.0d0
     mmetal_tr=0.0d0; mmetal_inflow_tr=0.0d0; mmetal_outflow_tr=0.0d0

     ethermal_tr=0.0d0; ethermal_in_tr=0.0d0; ethermal_out_tr=0.0d0
     ekin_bulk_tr=0.0d0; ekin_bulk_in_tr=0.0d0; ekin_bulk_out_tr=0.0d0
     ekin_tr=0.0d0; ekin_in_tr=0.0d0; ekin_out_tr=0.0d0
     emag_tr=0.0d0; emag_in_tr=0.0d0; emag_out_tr=0.0d0
     ecrs_tr=0.0d0; ecrs_in_tr=0.0d0; ecrs_out_tr=0.0d0

     if (myid.eq.1) then
        do it=1,ntracker
           tracker(it)%mgas_lev(ilevel)          =0.0d0
           tracker(it)%mgas_inflow_lev(ilevel)   =0.0d0
           tracker(it)%mgas_outflow_lev(ilevel)  =0.0d0
           tracker(it)%Lgas_lev(ilevel,:)        =0.0d0
           tracker(it)%Lgas_inflow_lev(ilevel,:) =0.0d0
           tracker(it)%Lgas_outflow_lev(ilevel,:)=0.0d0

           tracker(it)%mmetal_lev(ilevel)        =0.0d0
           tracker(it)%mmetal_inflow_lev(ilevel) =0.0d0
           tracker(it)%mmetal_outflow_lev(ilevel)=0.0d0

           tracker(it)%ethermal_lev(ilevel)      =0.0d0
           tracker(it)%ethermal_in_lev(ilevel)   =0.0d0
           tracker(it)%ethermal_out_lev(ilevel)  =0.0d0

           tracker(it)%ekin_bulk_lev(ilevel)     =0.0d0
           tracker(it)%ekin_bulk_in_lev(ilevel)  =0.0d0
           tracker(it)%ekin_bulk_out_lev(ilevel) =0.0d0

           tracker(it)%ekin_lev(ilevel)          =0.0d0
           tracker(it)%ekin_in_lev(ilevel)       =0.0d0
           tracker(it)%ekin_out_lev(ilevel)      =0.0d0

           tracker(it)%emag_lev(ilevel)          =0.0d0
           tracker(it)%emag_in_lev(ilevel)       =0.0d0
           tracker(it)%emag_out_lev(ilevel)      =0.0d0

           tracker(it)%ecrs_lev(ilevel)          =0.0d0
           tracker(it)%ecrs_in_lev(ilevel)       =0.0d0
           tracker(it)%ecrs_out_lev(ilevel)      =0.0d0
        end do
     end if
     skip_loc=(/0.0d0,0.0d0,0.0d0/)
     if(ndim>0)skip_loc(1)=dble(icoarse_min)
     if(ndim>1)skip_loc(2)=dble(jcoarse_min)
     if(ndim>2)skip_loc(3)=dble(kcoarse_min)
     ! Cells center position relative to grid center position
     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do
  end if
#endif

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  mass_all=0.0d0; mass_loc=0.0d0
  ekin_all=0.0d0; ekin_loc=0.0d0
  emag_all=0.0d0; emag_loc=0.0d0
  eint_all=0.0d0; eint_loc=0.0d0
#if NCR>0
  ecrs_all=0.0d0; ecrs_loc=0.0d0
#endif
  dt_all=dtnew(ilevel); dt_loc=dt_all

  ! Mesh spacing at that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale
  vol=dx**ndim

  if (ischeme .eq. 1) then
     CALL velocity_fine(ilevel)
     do ivar=1,nvar+3
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     end do
     if(simple_boundary)call make_boundary_hydro(ilevel)
  endif
 
  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Loop over cells
     do ind=1,twotondim        
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=ind_grid(i)+iskip
        end do
        
        ! Gather leaf cells
        nleaf=0
        do i=1,ngrid
           if(son(ind_cell(i))==0)then
              nleaf=nleaf+1
              ind_leaf(nleaf)=ind_cell(i)
#ifdef TRACKER
              ! Gather cells positions
              if (track_quantities) xpos(:,i)=(xg(ind_grid(i),:)+xc(ind,:)-skip_loc(:))*scale
#endif
           end if
        end do

        ! Gather hydro variables
        do ivar=1,nvar+3
           do i=1,nleaf
              uu(i,ivar)=uold(ind_leaf(i),ivar)
           end do
        end do
        
        ! Gather gravitational acceleration
        gg=0.0d0
        if(poisson)then
           do idim=1,ndim
              do i=1,nleaf
                 gg(i,idim)=f(ind_leaf(i),idim)
              end do
           end do
        end if
        
        ! Compute total mass
        do i=1,nleaf
           mass_loc=mass_loc+uu(i,1)*vol
        end do
        
        ! Compute total energy
        do i=1,nleaf
           ekin_loc=ekin_loc+uu(i,5)*vol
           !if(isnan(ekin_loc).or.ekin_loc>1.0d0)then
           !!if(myid==1)then
           ! write(*,'("ekin_loc:",6(1X,1PE14.7))'),i,ekin_loc,uu(i,1),uu(i,2:4)
           ! write(*,'("e_loc:",6(1X,1PE14.7))'),i,uu(i,5),vol,xpos(:,i)
           !endif
        end do
        
        ! Compute total magnetic energy
        do ivar=1,3
           do i=1,nleaf
              emag_loc=emag_loc+0.125d0*(uu(i,5+ivar)+uu(i,nvar+ivar))**2*vol
           end do
        end do
        
        ! Compute total internal energy
        do i=1,nleaf
           eint_loc=eint_loc+uu(i,5)*vol
        end do
        do ivar=1,3
           do i=1,nleaf
              eint_loc=eint_loc-0.5d0*uu(i,1+ivar)**2/uu(i,1)*vol &
                   & -0.125d0*(uu(i,5+ivar)+uu(i,nvar+ivar))**2*vol
           end do
        end do
#if NENER>0
        do ivar=1,nener
           do i=1,nleaf
              eint_loc=eint_loc-uu(i,8+ivar)*vol
           end do
        end do
#endif

        ! Compute cosmic rays energy
#if NCR>0
        ivar=1
        do i=1,nleaf
           ecrs_loc=ecrs_loc+uu(i,8+ivar)*vol
        end do        
#endif

#ifdef TRACKER
        QuantCompute: if (track_quantities) then
           do it=1,ntracker
              ! Skip unrequired trackers
              if (.not.(tracker(it)%active.and.tracker(it)%quantities)) cycle

              cell_contrib: do i=1,nleaf
                 ! Compute cells distance
                 xpos_cell=xpos(:,i)-tracker(it)%xpos
                 dd_cell2=sum( xpos_cell**2 )

                 if (dd_cell2.lt.tracker(it)%rmeasure2) then    ! Contribute to the tracker region
                    ! Compute cell properties
                    mcell=uu(i,1)*vol
                    m_metal=uu(i,imetal)*vol
                    vel_cell=uu(i,2:4)/max(uu(i,1),smallr)
                    vel_cell_COM=vel_cell-tracker(it)%vel
                    Lcell(1) = xpos_cell(2)*vel_cell_COM(3)-xpos_cell(3)*vel_cell_COM(2)
                    Lcell(2) = xpos_cell(3)*vel_cell_COM(1)-xpos_cell(1)*vel_cell_COM(3)
                    Lcell(3) = xpos_cell(1)*vel_cell_COM(2)-xpos_cell(2)*vel_cell_COM(1)
                    Lcell=Lcell*mcell

                    ! Contribute to cell masses
                    mgas_tr(it)=mgas_tr(it)+mcell
                    mmetal_tr(it)=mmetal_tr(it)+m_metal
                    Lgas_tr(it,:)=Lgas_tr(it,:)+Lcell

                    ! Energies contribution
                    ethermal_tr(it)=ethermal_tr(it)+uu(i,5)*vol
ekin_bulk_tr(it)=ekin_bulk_tr(it)+mcell*sum(vel_cell**2)
                    ekin_tr(it)=ekin_tr(it)+mcell*sum(  vel_cell_COM**2  )
#ifdef SOLVERmhd
                    ! Magnetic energy in the cell
                    emag_tr(it)=emag_tr(it)+0.125*vol* ( &
                         & (uu(i,6)+uu(i,nvar+1))**2 + &
                         & (uu(i,7)+uu(i,nvar+2))**2 + &
                         & (uu(i,8)+uu(i,nvar+3))**2 )
#endif
#if NCR>0
                    ! Cosmic ray energy in the cell
                    ecrs_tr(it)=ecrs_tr(it)+uu(i,9)*vol
#endif
                 else if (dd_cell2.lt.tracker(it)%rbound2) then ! Contribute to the tracker boundary
                    ! Compute cell properties
                    mcell=uu(i,1)*vol
                    m_metal=uu(i,imetal)*vol
                    vel_cell=uu(i,2:4)/uu(i,1)/max(uu(i,1),smallr)
                    vel_cell_COM=vel_cell-tracker(it)%vel
                    Lcell(1) = xpos_cell(2)*vel_cell_COM(3)-xpos_cell(3)*vel_cell_COM(2)
                    Lcell(2) = xpos_cell(3)*vel_cell_COM(1)-xpos_cell(1)*vel_cell_COM(3)
                    Lcell(3) = xpos_cell(1)*vel_cell_COM(2)-xpos_cell(2)*vel_cell_COM(1)
                    Lcell=Lcell*mcell
                    vrad_cell=0.0d0
                    if (dd_cell2.ne.0.0d0) vrad_cell=sum(vel_cell_COM*xpos_cell)/sqrt(dd_cell2)

                    if (vrad_cell.gt.0.0d0) then
                       mgas_outflow_tr(it)=mgas_outflow_tr(it)+mcell*vrad_cell
                       mmetal_outflow_tr(it)=mmetal_outflow_tr(it)+m_metal*vrad_cell
Lgas_outflow_tr(it,:)=Lgas_outflow_tr(it,:)+Lcell*vrad_cell

                       ethermal_out_tr(it)=ethermal_out_tr(it)+uu(i,5)*vol*vrad_cell
                       ekin_bulk_out_tr(it)=ekin_bulk_out_tr(it)+mcell*sum(vel_cell**2)*vrad_cell
                       ekin_out_tr(it)=ekin_out_tr(it)+mcell*sum(vel_cell_COM**2)*vrad_cell
#ifdef SOLVERmhd
                       ! Magnetic energy in the cell
                       emag_out_tr(it)=emag_out_tr(it)+0.125*vol*vrad_cell* ( &
                            & (uu(i,6)+uu(i,nvar+1))**2 + &
                            & (uu(i,7)+uu(i,nvar+2))**2 + &
                            & (uu(i,8)+uu(i,nvar+3))**2 )
#endif
#if NCR>0
                       ! Cosmic ray energy in the cell
                       ecrs_out_tr(it)=ecrs_out_tr(it)+uu(i,9)*vol*vrad_cell
#endif
                    else
                       mgas_inflow_tr(it)=mgas_inflow_tr(it)+mcell*vrad_cell
                       mmetal_inflow_tr(it)=mmetal_inflow_tr(it)+m_metal*vrad_cell
                       Lgas_inflow_tr(it,:)=Lgas_inflow_tr(it,:)+Lcell*vrad_cell

                       ethermal_in_tr(it)=ethermal_in_tr(it)+uu(i,5)*vol*vrad_cell
                       ekin_bulk_in_tr(it)=ekin_bulk_in_tr(it)+mcell*sum(vel_cell**2)*vrad_cell
                       ekin_in_tr(it)=ekin_in_tr(it)+mcell*sum(vel_cell_COM**2)*vrad_cell
#ifdef SOLVERmhd
                       ! Magnetic energy in the cell
                       emag_in_tr(it)=emag_in_tr(it)+0.125*vol*vrad_cell* ( &
                            & (uu(i,6)+uu(i,nvar+1))**2 + &
                            & (uu(i,7)+uu(i,nvar+2))**2 + &
                            & (uu(i,8)+uu(i,nvar+3))**2 )
#endif
#if NCR>0
                       ! Cosmic ray energy in the cell
                       ecrs_in_tr(it)=ecrs_in_tr(it)+uu(i,9)*vol*vrad_cell
#endif

                    end if
                 end if
              end do cell_contrib
           end do
        end if QuantCompute
#endif

        
        ! Compute CFL time-step
        if(nleaf>0)then
           call cmpdt(uu,gg,dx,dt_lev,nleaf)
           dt_loc=min(dt_loc,dt_lev)
        end if
        
     end do
     ! End loop over cells
     
  end do
  ! End loop over grids

  ! Compute global quantities
#ifndef WITHOUTMPI
  comm_buffin(1)=mass_loc
  comm_buffin(2)=ekin_loc
  comm_buffin(3)=eint_loc
  comm_buffin(4)=emag_loc
#if NCR>0
  comm_buffin(5)=ecrs_loc
#endif  
  call MPI_ALLREDUCE(comm_buffin,comm_buffout,ncomm,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dt_loc     ,dt_all      ,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
       &MPI_COMM_WORLD,info)
  mass_all=comm_buffout(1)
  ekin_all=comm_buffout(2)
  eint_all=comm_buffout(3)
  emag_all=comm_buffout(4)
#if NCR>0
  ecrs_all=comm_buffout(5)
#endif
#endif
  
#ifdef WITHOUTMPI
  mass_all=mass_loc
  ekin_all=ekin_loc
  eint_all=eint_loc
  emag_all=emag_loc
#if NCR>0
  ecrs_all=ecrs_loc
#endif
  dt_all=dt_loc
#endif

  mass_tot=mass_tot+mass_all
  ekin_tot=ekin_tot+ekin_all
  eint_tot=eint_tot+eint_all
  emag_tot=emag_tot+emag_all
#if NCR>0
  ecrs_tot=ecrs_tot+ecrs_all
#endif
  dtnew(ilevel)=MIN(dtnew(ilevel),dt_all)

  ! Communicate tracker properties to master processor
#ifdef TRACKER
if (track_quantities) then
     do it=1,ntracker
        ! Skip unrequired trackers
        if (.not.(tracker(it)%active.and.tracker(it)%quantities)) cycle
        icount=0
        icount=icount+1; comm_buffin_track(icount)         =mgas_tr(it)
        icount=icount+1; comm_buffin_track(icount)         =mgas_inflow_tr(it)
        icount=icount+1; comm_buffin_track(icount)         =mgas_outflow_tr(it)
        icount=icount+1; comm_buffin_track(icount:icount+2)=Lgas_tr(it,:)
        icount=icount+1; icount=icount+1
        icount=icount+1; comm_buffin_track(icount:icount+2)=Lgas_inflow_tr(it,:)
        icount=icount+1; icount=icount+1
        icount=icount+1; comm_buffin_track(icount:icount+2)=Lgas_outflow_tr(it,:)
        icount=icount+1; icount=icount+1
        icount=icount+1; comm_buffin_track(icount)         =mmetal_tr(it)
        icount=icount+1; comm_buffin_track(icount)         =mmetal_inflow_tr(it)
        icount=icount+1; comm_buffin_track(icount)         =mmetal_outflow_tr(it)
        ! Energies
        icount=icount+1; comm_buffin_track(icount)         =ethermal_tr(it)
        icount=icount+1; comm_buffin_track(icount)         =ethermal_in_tr(it)
        icount=icount+1; comm_buffin_track(icount)         =ethermal_out_tr(it)

        icount=icount+1; comm_buffin_track(icount)         =ekin_bulk_tr(it)
        icount=icount+1; comm_buffin_track(icount)         =ekin_bulk_in_tr(it)
        icount=icount+1; comm_buffin_track(icount)         =ekin_bulk_out_tr(it)

        icount=icount+1; comm_buffin_track(icount)         =ekin_tr(it)
        icount=icount+1; comm_buffin_track(icount)         =ekin_in_tr(it)
        icount=icount+1; comm_buffin_track(icount)         =ekin_out_tr(it)
#ifdef SOLVERmhd
        icount=icount+1; comm_buffin_track(icount)         =emag_tr(it)
        icount=icount+1; comm_buffin_track(icount)         =emag_in_tr(it)
        icount=icount+1; comm_buffin_track(icount)         =emag_out_tr(it)
#if NCR>0
        icount=icount+1; comm_buffin_track(icount)         =ecrs_tr(it)
        icount=icount+1; comm_buffin_track(icount)         =ecrs_in_tr(it)
        icount=icount+1; comm_buffin_track(icount)         =ecrs_out_tr(it)
#endif
#endif
        ! Communicate tracker properties
        call MPI_REDUCE(comm_buffin_track,comm_buffout_track,ncomm_track, &
             & MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,info)

! Return properties to tracker                                                                                                                                                                                           
        master_core: if (myid.eq.1) then
           icount=0
           icount=icount+1; tracker(it)%mgas_lev(ilevel)          =comm_buffout_track(icount)
           icount=icount+1; tracker(it)%mgas_inflow_lev(ilevel)   =comm_buffout_track(icount)
           icount=icount+1; tracker(it)%mgas_outflow_lev(ilevel)  =comm_buffout_track(icount)
           icount=icount+1; tracker(it)%Lgas_lev(ilevel,:)        =comm_buffout_track(icount:icount+2)
           icount=icount+1; icount=icount+1
           icount=icount+1; tracker(it)%Lgas_inflow_lev(ilevel,:) =comm_buffout_track(icount:icount+2)
           icount=icount+1; icount=icount+1
           icount=icount+1; tracker(it)%Lgas_outflow_lev(ilevel,:)=comm_buffout_track(icount:icount+2)
           icount=icount+1; icount=icount+1

           icount=icount+1; tracker(it)%mmetal_lev(ilevel)        =comm_buffout_track(icount)
           icount=icount+1; tracker(it)%mmetal_inflow_lev(ilevel) =comm_buffout_track(icount)
           icount=icount+1; tracker(it)%mmetal_outflow_lev(ilevel)=comm_buffout_track(icount)

           icount=icount+1; tracker(it)%ethermal_lev(ilevel)      =comm_buffout_track(icount)
           icount=icount+1; tracker(it)%ethermal_in_lev(ilevel)   =comm_buffout_track(icount)
           icount=icount+1; tracker(it)%ethermal_out_lev(ilevel)  =comm_buffout_track(icount)

           icount=icount+1; tracker(it)%ekin_bulk_lev(ilevel)     =comm_buffout_track(icount)
           icount=icount+1; tracker(it)%ekin_bulk_in_lev(ilevel)  =comm_buffout_track(icount)
           icount=icount+1; tracker(it)%ekin_bulk_out_lev(ilevel) =comm_buffout_track(icount)

           icount=icount+1; tracker(it)%ekin_lev(ilevel)          =comm_buffout_track(icount)
           icount=icount+1; tracker(it)%ekin_in_lev(ilevel)       =comm_buffout_track(icount)
           icount=icount+1; tracker(it)%ekin_out_lev(ilevel)      =comm_buffout_track(icount)
#ifdef SOLVERmhd
           icount=icount+1; tracker(it)%emag_lev(ilevel)          =comm_buffout_track(icount)
           icount=icount+1; tracker(it)%emag_in_lev(ilevel)       =comm_buffout_track(icount)
           icount=icount+1; tracker(it)%emag_out_lev(ilevel)      =comm_buffout_track(icount)
#if NCR>0
           icount=icount+1; tracker(it)%ecrs_lev(ilevel)          =comm_buffout_track(icount)
           icount=icount+1; tracker(it)%ecrs_in_lev(ilevel)       =comm_buffout_track(icount)
           icount=icount+1; tracker(it)%ecrs_out_lev(ilevel)      =comm_buffout_track(icount)
#endif
#endif

           ! Collapse level quantities
           dummy_dp=tracker(it)%delta_bound_code
           if (dummy_dp.eq.0.0d0) dummy_dp=1.0d99
           tracker(it)%mgas=sum(tracker(it)%mgas_lev(:))
           tracker(it)%mgas_inflow=sum(tracker(it)%mgas_inflow_lev(:))/dummy_dp
           tracker(it)%mgas_outflow=sum(tracker(it)%mgas_outflow_lev(:))/dummy_dp
           do idim=1,ndim
              if (tracker(it)%mgas.gt.0.0d0) then
                 tracker(it)%Lgas(idim)=sum(tracker(it)%Lgas_lev(:,idim))/tracker(it)%mgas
              else
                 tracker(it)%Lgas(idim)=0.0d0
              end if
              if (tracker(it)%mgas_inflow.ne.0.0d0) then
                 tracker(it)%Lgas_in(idim)=sum(tracker(it)%Lgas_inflow_lev(:,idim))/(tracker(it)%mgas_inflow*dummy_dp)
              else
                 tracker(it)%Lgas_in(idim)=0.0d0
end if
              if (tracker(it)%mgas_outflow.ne.0.0d0) then
                 tracker(it)%Lgas_out(idim)=sum(tracker(it)%Lgas_outflow_lev(:,idim))/(tracker(it)%mgas_outflow*dummy_dp)
              else
                 tracker(it)%Lgas_out(idim)=0.0d0
              end if
           end do
           tracker(it)%mmetal=sum(tracker(it)%mmetal_lev(:))
           tracker(it)%mmetal_inflow=sum(tracker(it)%mmetal_inflow_lev(:))/dummy_dp
           tracker(it)%mmetal_outflow=sum(tracker(it)%mmetal_outflow_lev(:))/dummy_dp

           tracker(it)%ethermal=sum(tracker(it)%ethermal_lev(:))
           tracker(it)%ethermal_in=sum(tracker(it)%ethermal_in_lev(:))/dummy_dp
           tracker(it)%ethermal_out=sum(tracker(it)%ethermal_out_lev(:))/dummy_dp

           tracker(it)%ekin_bulk=sum(tracker(it)%ekin_bulk_lev(:))
           tracker(it)%ekin_bulk_in=sum(tracker(it)%ekin_bulk_in_lev(:))/dummy_dp
           tracker(it)%ekin_bulk_out=sum(tracker(it)%ekin_bulk_out_lev(:))/dummy_dp

           tracker(it)%ekin=sum(tracker(it)%ekin_lev(:))
           tracker(it)%ekin_in=sum(tracker(it)%ekin_in_lev(:))/dummy_dp
           tracker(it)%ekin_out=sum(tracker(it)%ekin_out_lev(:))/dummy_dp
#ifdef SOLVERmhd
           tracker(it)%emag=sum(tracker(it)%emag_lev(:))
           tracker(it)%emag_in=sum(tracker(it)%emag_in_lev(:))/dummy_dp
           tracker(it)%emag_out=sum(tracker(it)%emag_out_lev(:))/dummy_dp
#if NCR>0
           tracker(it)%ecrs=sum(tracker(it)%ecrs_lev(:))
           tracker(it)%ecrs_in=sum(tracker(it)%ecrs_in_lev(:))/dummy_dp
tracker(it)%ecrs_out=sum(tracker(it)%ecrs_out_lev(:))/dummy_dp
#endif
#endif
        end if master_core
     end do
  end if
#endif

111 format('   Entering courant_fine for level ',I2)

end subroutine courant_fine
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine velocity_fine(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !----------------------------------------------------------
  ! This routine computes the gravitational acceleration,
  ! the maximum density rho_max, and the potential energy
  !----------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz
  integer::info,ibound,nx_loc,idim,neul=5
  real(dp)::dx,dx_loc,scale,d,u,v,w,A,B,C
  real(kind=8)::rho_max_loc,rho_max_all,epot_loc,epot_all
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:3),save::vv
 
  if(numbtot(1,ilevel)==0)return

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel
  
  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=dble(nx_loc)/boxlen
  dx_loc=dx/scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do
  
  !-------------------------------------
  ! Compute analytical velocity field
  !-------------------------------------
  ncache=active(ilevel)%ngrid
  
  ! Loop over grids by vector sweeps
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Loop over cells
     do ind=1,twotondim
        
        ! Gather cell indices
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        
        ! Gather cell centre positions
        do idim=1,ndim
           do i=1,ngrid
              xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
           end do
        end do
        ! Rescale position from code units to user units
        do idim=1,ndim
           do i=1,ngrid
              xx(i,idim)=(xx(i,idim)-skip_loc(idim))/scale
           end do
        end do
        
        ! Impose analytical velocity field
        call velana(xx,vv,dx_loc,t,ngrid)
        
        ! Impose induction variables
        do i=1,ngrid
           uold(ind_cell(i),1)=1.0
        end do
        do idim=1,3
           do i=1,ngrid
              uold(ind_cell(i),idim+1)=vv(i,idim)
           end do
        end do
        ! Update total energy
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           u=uold(ind_cell(i),2)/d
           v=uold(ind_cell(i),3)/d
           w=uold(ind_cell(i),4)/d
           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))
           uold(ind_cell(i),neul)=1.0+0.5*d*(u**2+v**2+w**2)+0.5*(A**2+B**2+C**2)
        end do
        
     end do
     ! End loop over cells

  end do
  ! End loop over grids

end subroutine velocity_fine
!#########################################################
!#########################################################
!#########################################################
!#########################################################
