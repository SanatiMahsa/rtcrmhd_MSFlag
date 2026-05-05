!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flow  
  use amr_commons
  use hydro_commons, ONLY: nvar, uold
  implicit none

  integer::ilevel,ivar
  
  if(verbose)write(*,*)'Entering init_flow'
  do ilevel=nlevelmax,1,-1
     if(ilevel>=levelmin)call init_flow_fine(ilevel)
     call upload_fine(ilevel)
     do ivar=1,nvar+3
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     end do
     if(simple_boundary)call make_boundary_hydro(ilevel)
  end do
  if(verbose)write(*,*)'Complete init_flow'

end subroutine init_flow
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flow_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  
  integer::i,icell,igrid,ncache,iskip,ngrid,ilun
  integer::ind,idim,ivar,ix,iy,iz,nx_loc
  integer::i1,i2,i3,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
  integer::buf_count,info,nvar_in
  integer ,dimension(1:nvector),save::ind_grid,ind_cell

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_Bco
  real(dp)::dx,rr,vx,vy,vz,bx,by,bz,ek,ei,em,pp,xx1,xx2,xx3,dx_loc,scale,xval
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector)       ,save::vv
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:nvar+3),save::uu

  real(dp),allocatable,dimension(:,:,:)::init_array
  real(kind=4),allocatable,dimension(:,:)  ::init_plane

  logical::error,ok_file1,ok_file2,ok_file3,ok_file
  character(LEN=80)::filename
  character(LEN=5)::nchar,ncharvar

  integer,parameter::tag=1107
   integer::dummy_io,info2

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  ! Initialise received B_ave in comoving units
  scale_Bco = (sqrt(scale_d) * scale_l / scale_t) * (aexp**2)

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel
  
  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do
  
  ! Local constants
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  ncache=active(ilevel)%ngrid

  !--------------------------------------
  ! Compute initial conditions from files
  !--------------------------------------
  filename=TRIM(initfile(ilevel))//'/ic_d'
  INQUIRE(file=filename,exist=ok_file1)
  if(multiple)then
     filename=TRIM(initfile(ilevel))//'/dir_deltab/ic_deltab.00001'
     INQUIRE(file=filename,exist=ok_file2)
  else
     filename=TRIM(initfile(ilevel))//'/ic_deltab'
     INQUIRE(file=filename,exist=ok_file2)
  endif
  ok_file = ok_file1 .or. ok_file2
  if(ok_file)then

     !-------------------------------------------------------------------------
     ! First step: compute level boundaries in terms of initial condition array
     !-------------------------------------------------------------------------
     if(ncache>0)then
     i1_min=n1(ilevel)+1; i1_max=0
     i2_min=n2(ilevel)+1; i2_max=0
     i3_min=n3(ilevel)+1; i3_max=0
     do ind=1,twotondim           
        do i=1,ncache
           igrid=active(ilevel)%igrid(i)
           xx1=xg(igrid,1)+xc(ind,1)-skip_loc(1)
           xx1=(xx1*(dxini(ilevel)/dx)-xoff1(ilevel))/dxini(ilevel)
           xx2=xg(igrid,2)+xc(ind,2)-skip_loc(2)
           xx2=(xx2*(dxini(ilevel)/dx)-xoff2(ilevel))/dxini(ilevel)
           xx3=xg(igrid,3)+xc(ind,3)-skip_loc(3)
           xx3=(xx3*(dxini(ilevel)/dx)-xoff3(ilevel))/dxini(ilevel)
           i1_min=MIN(i1_min,int(xx1)+1)
           i1_max=MAX(i1_max,int(xx1)+1)
           i2_min=MIN(i2_min,int(xx2)+1)
           i2_max=MAX(i2_max,int(xx2)+1)
           i3_min=MIN(i3_min,int(xx3)+1)
           i3_max=MAX(i3_max,int(xx3)+1)
        end do
     end do
     error=.false.
     if(i1_min<1.or.i1_max>n1(ilevel))error=.true.
     if(i2_min<1.or.i2_max>n2(ilevel))error=.true.
     if(i3_min<1.or.i3_max>n3(ilevel))error=.true.
     if(error) then
        write(*,*)'Some grid are outside initial conditions sub-volume'
        write(*,*)'for ilevel=',ilevel
        write(*,*)i1_min,i1_max
        write(*,*)i2_min,i2_max
        write(*,*)i3_min,i3_max
        write(*,*)n1(ilevel),n2(ilevel),n3(ilevel)
        call clean_stop
     end if
     endif

     !-----------------------------------------
     ! Second step: read initial condition file
     !-----------------------------------------
     ! Allocate initial conditions array
     if(ncache>0)allocate(init_array(i1_min:i1_max,i2_min:i2_max,i3_min:i3_max))
     allocate(init_plane(1:n1(ilevel),1:n2(ilevel)))
     ! Loop over input variables
     do ivar=1,nvar+3
        if(cosmo)then
           ! Read baryons initial overdensity and displacement at a=aexp
           if(multiple)then
              call title(myid,nchar)
              if(ivar==1)filename=TRIM(initfile(ilevel))//'/dir_deltab/ic_deltab.'//TRIM(nchar)
              if(ivar==2)filename=TRIM(initfile(ilevel))//'/dir_velcx/ic_velcx.'//TRIM(nchar)
              if(ivar==3)filename=TRIM(initfile(ilevel))//'/dir_velcy/ic_velcy.'//TRIM(nchar)
              if(ivar==4)filename=TRIM(initfile(ilevel))//'/dir_velcz/ic_velcz.'//TRIM(nchar)
              if(ivar==5)filename=TRIM(initfile(ilevel))//'/dir_tempb/ic_tempb.'//TRIM(nchar)
           else
              if(ivar==1)filename=TRIM(initfile(ilevel))//'/ic_deltab'
              if(ivar==2)filename=TRIM(initfile(ilevel))//'/ic_velcx'
              if(ivar==3)filename=TRIM(initfile(ilevel))//'/ic_velcy'
              if(ivar==4)filename=TRIM(initfile(ilevel))//'/ic_velcz'
              if(ivar==5)filename=TRIM(initfile(ilevel))//'/ic_tempb'
           endif
        else
           ! Read primitive variables
           if(ivar==1)filename=TRIM(initfile(ilevel))//'/ic_d'
           if(ivar==2)filename=TRIM(initfile(ilevel))//'/ic_u'
           if(ivar==3)filename=TRIM(initfile(ilevel))//'/ic_v'
           if(ivar==4)filename=TRIM(initfile(ilevel))//'/ic_w'
           if(ivar==5)filename=TRIM(initfile(ilevel))//'/ic_p'
        endif
        if(ivar==6)filename=TRIM(initfile(ilevel))//'/ic_bxleft'
        if(ivar==7)filename=TRIM(initfile(ilevel))//'/ic_byleft'
        if(ivar==8)filename=TRIM(initfile(ilevel))//'/ic_bzleft'
        if(ivar==nvar+1)filename=TRIM(initfile(ilevel))//'/ic_bxright'
        if(ivar==nvar+2)filename=TRIM(initfile(ilevel))//'/ic_byright'
        if(ivar==nvar+3)filename=TRIM(initfile(ilevel))//'/ic_bzright'
        call title(ivar,ncharvar)
#if NCR>0
        if(ivar==9)filename=TRIM(initfile(ilevel))//'/ic_crs'
        if(ivar>9.and.ivar<=nvar)then
           call title(ivar-9,ncharvar)
           filename=TRIM(initfile(ilevel))//'/ic_pvar_'//TRIM(ncharvar)
        end if
#else
        if(ivar>8.and.ivar<=nvar)then
           call title(ivar-8,ncharvar)
           filename=TRIM(initfile(ilevel))//'/ic_pvar_'//TRIM(ncharvar)
        end if
#endif

        INQUIRE(file=filename,exist=ok_file3)
        if(ok_file3)then
           ! Reading the existing file
           if(myid==1)write(*,*)"ivar=",ivar,'Reading file '//TRIM(filename)
           if(multiple)then
              ilun=ncpu+myid+10
              ! Wait for the token
#ifndef WITHOUTMPI
              if(IOGROUPSIZE>0) then
                 if (mod(myid-1,IOGROUPSIZE)/=0) then
                    call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                         & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
                 end if
              endif
#endif
              open(ilun,file=filename,form='unformatted')
              rewind ilun
              read(ilun) ! skip first line
              do i3=1,n3(ilevel)
                 read(ilun) ((init_plane(i1,i2),i1=1,n1(ilevel)),i2=1,n2(ilevel))
                 if(ncache>0)then
                    if(i3.ge.i3_min.and.i3.le.i3_max)then
                       init_array(i1_min:i1_max,i2_min:i2_max,i3) = &
                            & init_plane(i1_min:i1_max,i2_min:i2_max)
                    end if
                 endif
              end do
              close(ilun)
              ! Send the token
#ifndef WITHOUTMPI
              if(IOGROUPSIZE>0) then
                 if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
                    dummy_io=1
                    call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
                         & MPI_COMM_WORLD,info2)
                 end if
              endif
#endif        
           else
              if(myid==1)then
                 open(10,file=filename,form='unformatted')
                 rewind 10
                 read(10) ! skip first line
              endif
              do i3=1,n3(ilevel)
                 if(myid==1)then
                    read(10) ((init_plane(i1,i2),i1=1,n1(ilevel)),i2=1,n2(ilevel))
                 else
                    init_plane=0.0
                 endif
                 buf_count=n1(ilevel)*n2(ilevel)
#ifndef WITHOUTMPI
                 call MPI_BCAST(init_plane,buf_count,MPI_REAL,0,MPI_COMM_WORLD,info)
#endif
                 if(ncache>0)then
                    if(i3.ge.i3_min.and.i3.le.i3_max)then
                       init_array(i1_min:i1_max,i2_min:i2_max,i3) = &
                            & init_plane(i1_min:i1_max,i2_min:i2_max)
                    end if
                 endif
              end do
              if(myid==1)close(10)
           endif
        else
           ! If file doesn't exist, initialize variable to default value
           ! In most cases, this is zero (you can change that if necessary)
           if(myid==1)write(*,*)"ivar=",ivar,'File '//TRIM(filename)//' not found'
           if(myid==1)write(*,*)'Initialize corresponding variable to default value'
           if(ncache>0)then
              init_array=0d0
              ! Default value for metals
              if(cosmo.and.ivar==imetal.and.metal)init_array=z_ave*0.02
              ! Default value for Bz
              if(cosmo.and.ivar==8)init_array=B_ave/scale_Bco
              if(cosmo.and.ivar==nvar+3)init_array=B_ave/scale_Bco
              ! Default value for ionization fraction
              if(cosmo)xval=sqrt(omega_m)/(h0/100.*omega_b) ! From the book of Peebles p. 17
              if(cosmo.and.ivar==ixion.and.aton)init_array=1.2d-5*xval
#if NCR>0
              ! Default value for cosmic rays
              if(ivar==inener)init_array=CR_ave
#endif
           end if
        end if

        if(ncache>0)then
           if(cosmo)then
              ! Rescale initial conditions to code units
              if(.not. cooling)T2_start = 1.356d-2/aexp**2
              if(ivar==1)init_array=(1.0+dfact(ilevel)*init_array)*omega_b/omega_m
              if(ivar==2)init_array=dfact(ilevel)*vfact(1)*dx_loc/dxini(ilevel)*init_array/vfact(ilevel)
              if(ivar==3)init_array=dfact(ilevel)*vfact(1)*dx_loc/dxini(ilevel)*init_array/vfact(ilevel)
              if(ivar==4)init_array=dfact(ilevel)*vfact(1)*dx_loc/dxini(ilevel)*init_array/vfact(ilevel)
#if NCR>0
              ! Initialisation for other NENERs
              ! if(ivar==inener)init_array=(1.0+init_array)*T2_start/scale_T2/(1d0+mu_electron/mu_ion)
              ! if(ivar==5)init_array=(1.0+init_array)*T2_start/scale_T2/(1d0+mu_ion/mu_electron)
              if(ivar==5)init_array=(1.0+init_array)*T2_start/scale_T2
              if(ivar==inener)init_array=(init_array)*T2_start/scale_T2              
#else
              if(ivar==5)init_array=(1.0+init_array)*T2_start/scale_T2
#endif
           end if
           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do i=1,ncache
                 igrid=active(ilevel)%igrid(i)
                 icell=igrid+iskip
                 xx1=xg(igrid,1)+xc(ind,1)-skip_loc(1)
                 xx1=(xx1*(dxini(ilevel)/dx)-xoff1(ilevel))/dxini(ilevel)
                 xx2=xg(igrid,2)+xc(ind,2)-skip_loc(2)
                 xx2=(xx2*(dxini(ilevel)/dx)-xoff2(ilevel))/dxini(ilevel)
                 xx3=xg(igrid,3)+xc(ind,3)-skip_loc(3)
                 xx3=(xx3*(dxini(ilevel)/dx)-xoff3(ilevel))/dxini(ilevel)
                 i1=int(xx1)+1
                 i1=int(xx1)+1
                 i2=int(xx2)+1
                 i2=int(xx2)+1
                 i3=int(xx3)+1
                 i3=int(xx3)+1
                 ! Scatter to corresponding primitive variable
                 uold(icell,ivar)=init_array(i1,i2,i3)
              end do
           end do
           ! End loop over cells
        endif
     end do
     ! End loop over input variables

     ! Deallocate initial conditions array
     if(ncache>0)deallocate(init_array)
     deallocate(init_plane) 

     !-------------------------------------------------------------------
     ! For cosmology runs: compute pressure, prevent negative density 
     !-------------------------------------------------------------------
     if(cosmo)then
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
              ! Prevent negative density
              do i=1,ngrid
                 rr=max(uold(ind_cell(i),1),0.1*omega_b/omega_m)
                 uold(ind_cell(i),1)=rr
              end do
              ! Compute pressure from temperature and density
              do i=1,ngrid
                 uold(ind_cell(i),5)=uold(ind_cell(i),1)*uold(ind_cell(i),5)
#if NENER>0
                  do ivar=0,nener-1
                     uold(ind_cell(i),inener+ivar)=uold(ind_cell(i),1)*uold(ind_cell(i),inener+ivar)
                  end do
#endif
              end do
           end do
           ! End loop over cells        
        end do
        ! End loop over grids
     end if

     !---------------------------------------------------
     ! Third step: compute initial conservative variables
     !---------------------------------------------------
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
           ! Compute total energy
           do i=1,ngrid
              rr=uold(ind_cell(i),1)
              vx=uold(ind_cell(i),2)
              vy=uold(ind_cell(i),3)
              vz=uold(ind_cell(i),4)
              pp=uold(ind_cell(i),5)
              bx=0.5d0*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
              by=0.5d0*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
              bz=0.5d0*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))
              ek=0.5d0*(vx**2+vy**2+vz**2)
              em=0.5d0*(bx**2+by**2+bz**2)
              ei=pp/(gamma-1.0)
#if NCR>0
              ei = ei + uold(ind_cell(i),9)/(gamma_rad(1)-1.0)
#endif
              vv(i)=ei+rr*ek+em
           end do
           ! Scatter to corresponding conservative variable
           do i=1,ngrid
              uold(ind_cell(i),5)=vv(i)
           end do
           ! Compute momentum density
           do ivar=1,ndim
              do i=1,ngrid
                 rr=uold(ind_cell(i),1)
                 vx=uold(ind_cell(i),ivar+1)
                 vv(i)=rr*vx
              end do
              ! Scatter to corresponding conservative variable
              do i=1,ngrid
                 uold(ind_cell(i),ivar+1)=vv(i)
              end do
           end do
           ! Compute passive variable density
#if NVAR > 8

#if NENER > 0
           do ivar=9+nener,nvar
#else
           do ivar=9,nvar
#endif
              do i=1,ngrid
                 rr=uold(ind_cell(i),1)
                 uold(ind_cell(i),ivar)=rr*uold(ind_cell(i),ivar)
              end do
#if NENER > 0
           end do
#else
           end do
#endif
#endif
        end do
        ! End loop over cells
        
     end do
     ! End loop over grids

  !-------------------------------------------------------
  ! Compute initial conditions from subroutine condinit
  !-------------------------------------------------------
  else

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
                 xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
              end do
           end do
           ! Call initial condition routine
           call condinit(xx,uu,dx_loc,ngrid)
           ! Scatter variables
           do ivar=1,nvar+3
              do i=1,ngrid
                 uold(ind_cell(i),ivar)=uu(i,ivar)
              end do
           end do
        end do
        ! End loop over cells
     end do
     ! End loop over grids

  end if
  
111 format('   Entering init_flow_fine for level ',I2)

end subroutine init_flow_fine
!################################################################
!################################################################
!################################################################
!################################################################


#ifdef DICE
subroutine init_flow_fine_dmf(ilevel,rfoverride,gclear)
  use amr_commons
  use hydro_commons
  use cooling_module
  ! use pm_commons, ONLY: up32
  use dice_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel

  logical,intent(in)::rfoverride,gclear

  integer::i,icell,igrid,ncache,iskip,ngrid,ilun
  integer::ind,idim,ivar,ix,iy,iz,nx_loc
  integer::i1,i2,i3,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
  integer::buf_count,info,nvar_in
  integer ,dimension(1:nvector),save::ind_grid,ind_cell

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx,rr,vx,vy,vz,ek,ei,pp,xx1,xx2,xx3,dx_loc,scale,xval
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector)       ,save::vv
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:nvar),save::uu
  real(dp)::axlen

  real(dp),allocatable,dimension(:,:,:)::init_array
  real(kind=4),allocatable,dimension(:,:)  ::init_plane

  logical::error,ok_file1,ok_file2,ok_file3,ok_file
  character(LEN=80)::filename
  character(LEN=5)::nchar,ncharvar

  integer,parameter::tag=1107
  integer::dummy_io,info2

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
    iz=(ind-1)/4
    iy=(ind-1-4*iz)/2
    ix=(ind-1-2*iy-4*iz)
    if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
    if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
    if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Local constants
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  ncache=active(ilevel)%ngrid

#ifndef WITHOUTMPI
  ! write(*,*) 'init_flow_fine check 549', myid
  call MPI_Barrier(MPI_COMM_WORLD,info)
#endif

  do i=1,MAXGAL
    if (ic_mag_scale_B(i) .EQ. 0.0) cycle
    ! renormalise axes
    axlen = SQRT(ic_mag_axis_x(i)**2 + ic_mag_axis_y(i)**2 + ic_mag_axis_z(i)**2)
    ic_mag_axis_x(i) = ic_mag_axis_x(i) / axlen
    ic_mag_axis_y(i) = ic_mag_axis_y(i) / axlen
    ic_mag_axis_z(i) = ic_mag_axis_z(i) / axlen
  enddo



!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!
!---------------------- DICE ADDING GAS TO AMR GRID ----------------------!
  if (gclear) call reset_uold3(ilevel)
    
#ifndef WITHOUTMPI
  call MPI_Barrier(MPI_COMM_WORLD,info)
#endif
  ! call reset_uold(ilevel)
  ! Update the grid using the gas particles read from the Gadget1 file
  ! NGP scheme is used

  if (gclear) call condinit_loc(ilevel,rfoverride)
  ! Reverse update boundaries
  ! call reset_uold3(ilevel)
  do ivar=1,nvar
    call make_virtual_reverse_dp(uold(:,ivar),ilevel)
    ! if (dmf_loop.gt.1) call make_virtual_reverse_dp_dmf(utmp(:,ivar),ilevel) 
  end do
  call make_virtual_reverse_int_dmf(readflag(:),ilevel)

  if (dmf_i.eq.dmf_loop) call init_uold(ilevel)
  do ivar=1,nvar
    call make_virtual_fine_dp(uold(:,ivar),ilevel)
    ! if (dmf_loop.gt.1) call make_virtual_fine_dp(utmp(:,ivar),ilevel)
  end do
  call make_virtual_fine_int(readflag(:),ilevel)
#ifndef WITHOUTMPI
! write(*,*) 'init_flow_fine check 476', myid
call MPI_Barrier(MPI_COMM_WORLD,info)
#endif
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!


111 format('   Entering init_flow_fine for level ',I2)

end subroutine init_flow_fine_dmf
#endif
!################################################################
!################################################################
!################################################################
!################################################################
subroutine region_condinit(x,q,dx,nn)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::nn
  real(dp)::dx
  real(dp),dimension(1:nvector,1:nvar+3)::q
  real(dp),dimension(1:nvector,1:ndim)::x

  integer::i,ivar,k
  real(dp)::vol,r,xn,yn,zn,en

  ! Set some (tiny) default values in case n_region=0
  q(1:nn,1)=smallr
  q(1:nn,2)=0.0d0
  q(1:nn,3)=0.0d0
  q(1:nn,4)=0.0d0
  q(1:nn,5)=smallr*smallc**2/gamma
  q(1:nn,6)=0.0d0
  q(1:nn,7)=0.0d0
  q(1:nn,8)=0.0d0
  q(1:nn,nvar+1)=0.0d0
  q(1:nn,nvar+2)=0.0d0
  q(1:nn,nvar+3)=0.0d0
#if NVAR > 8
  do ivar=9,nvar
     q(1:nn,ivar)=0.0d0
  end do
#endif

  ! Loop over initial conditions regions
  do k=1,nregion
     
     ! For "square" regions only:
     if(region_type(k) .eq. 'square')then
        ! Exponent of choosen norm
        en=exp_region(k)
        do i=1,nn
           ! Compute position in normalized coordinates
           xn=0.0d0; yn=0.0d0; zn=0.0d0
           xn=2.0d0*abs(x(i,1)-x_center(k))/length_x(k)
#if NDIM>1
           yn=2.0d0*abs(x(i,2)-y_center(k))/length_y(k)
#endif
#if NDIM>2
           zn=2.0d0*abs(x(i,3)-z_center(k))/length_z(k)
#endif
           ! Compute cell "radius" relative to region center
           if(exp_region(k)<10)then
              r=(xn**en+yn**en+zn**en)**(1.0/en)
           else
              r=max(xn,yn,zn)
           end if
           ! If cell lies within region,
           ! REPLACE primitive variables by region values
           if(r<1.0)then
              q(i,1)=d_region(k)
              q(i,2)=u_region(k)
              q(i,3)=v_region(k)
              q(i,4)=w_region(k)
              q(i,5)=p_region(k)
              q(i,6)=A_region(k)
              q(i,7)=B_region(k)
              q(i,8)=C_region(k)
              q(i,nvar+1)=A_region(k)
              q(i,nvar+2)=B_region(k)
              q(i,nvar+3)=C_region(k)
#if NENER>0
              do ivar=1,nener
                 q(i,8+ivar)=prad_region(k,ivar)
              enddo
#endif
#if NVAR>8+NENER
              do ivar=9+nener,nvar
                 q(i,ivar)=var_region(k,ivar-8-nener)
              end do
#endif
 
           end if
        end do
     end if
     
     ! For "point" regions only:
     if(region_type(k) .eq. 'point')then
        ! Volume elements
        vol=dx**ndim
        ! Compute CIC weights relative to region center
        do i=1,nn
           xn=1.0; yn=1.0; zn=1.0
           xn=max(1.0-abs(x(i,1)-x_center(k))/dx,0.0_dp)
#if NDIM>1
           yn=max(1.0-abs(x(i,2)-y_center(k))/dx,0.0_dp)
#endif
#if NDIM>2
           zn=max(1.0-abs(x(i,3)-z_center(k))/dx,0.0_dp)
#endif
           r=xn*yn*zn
           ! If cell lies within CIC cloud, 
           ! ADD to primitive variables the region values
           q(i,1)=q(i,1)+d_region(k)*r/vol
           q(i,2)=q(i,2)+u_region(k)*r
           q(i,3)=q(i,3)+v_region(k)*r
           q(i,4)=q(i,4)+w_region(k)*r
           q(i,5)=q(i,5)+p_region(k)*r/vol
#if NENER>0
           do ivar=1,nener
              q(i,8+ivar)=q(i,8+ivar)+prad_region(k,ivar)*r/vol
           enddo
#endif
#if NVAR>8+NENER
           do ivar=9+nener,nvar
              q(i,ivar)=var_region(k,ivar-8-nener)
           end do
#endif
        end do
     end if
  end do

  return
end subroutine region_condinit

#ifdef DICE
subroutine reset_uold3(ilevel)
  use amr_commons
  use hydro_commons
  use dice_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  integer::ncell
  !--------------------------------------------------------------------------
  ! This routine sets array uold to zero before calling
  ! the hydro scheme. uold is set to zero in virtual boundaries as well.
  !--------------------------------------------------------------------------
  integer::i,ivar,irad,ind,icpu,iskip,info, j
  
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
#ifndef WITHOUTMPI
  ! write(*,*) 'ireset_uold2 check 669', myid
  call MPI_Barrier(MPI_COMM_WORLD,info)
#endif

  do ind=1,twotondim
    iskip=ncoarse+(ind-1)*ngridmax
    do ivar=1,nvar
      do i=1,active(ilevel)%ngrid
        j = active(ilevel)%igrid(i)+iskip
        if (readflag(j).eq.0) uold(j,ivar) = 0D0
        if (readflag(j).eq.2) uold(j,ivar) = 0D0
        if (readflag(j).eq.3) uold(j,ivar) = 0D0
        
        if (dmf_i.gt.1) then
          if (readflag(j).eq.3) then
            !if(j > size(utmp,1)) then
            !  write(*,*) 'ERROR1: reset_uold3 - j exceeds utmp bounds, j=', j, 'size(utmp,1)=', size(utmp,1), 'myid=', myid
            !  stop
            !endif
            uold(j,ivar) = utmp(j,ivar)
          end if
        end if
      end do
    end do
  end do

  ! Set uold to 0 for virtual boundary cells
  do icpu=1,ncpu
    do ind=1,twotondim
      iskip=ncoarse+(ind-1)*ngridmax
      do ivar=1,nvar
        do i=1,reception(icpu,ilevel)%ngrid
          j = reception(icpu,ilevel)%igrid(i)+iskip
          if (readflag(j).eq.0) uold(j,ivar) = 0D0
          if (readflag(j).eq.2) uold(j,ivar) = 0D0
          if (readflag(j).eq.3) uold(j,ivar) = 0D0
          
          if (dmf_i.gt.1) then
            if (readflag(j).eq.3) then
              !if(j > size(utmp,1)) then
              !  write(*,*) 'ERROR2: reset_uold3 - Virtual boundary j exceeds utmp bounds, j=', j, 'size(utmp,1)=', size(utmp,1), 'myid=', myid
              !  stop
              !endif
              uold(j,ivar) = utmp(j,ivar)
            end if
          end if
        end do
      end do
    end do
  end do

#ifndef WITHOUTMPI
  call MPI_Barrier(MPI_COMM_WORLD,info)
#endif

111 format('   Entering reset_uold3 for level ',i2)

end subroutine reset_uold3


subroutine init_uold(ilevel)
  use amr_commons
  use hydro_commons
  use dice_commons
  implicit none
  integer::ilevel,info
  !--------------------------------------------------------------------------
  ! This routine sets array unew to its initial value uold before calling
  ! the hydro scheme. unew is set to zero in virtual boundaries.
  !--------------------------------------------------------------------------
  integer::i,ivar,irad,ind,icpu,iskip,idim
  real(dp)::d,u,v,w,e
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,IG_den, IG_TFl

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  ! write(*,'(a,es15.5)') 'IG_rho test 1', IG_rho
  IG_den = IG_rho * (1.989e33) / (dble(3.086e+18)**3)    !----- IG_rho change units from Msol pc^-3 to g cm^-3
  ! write(*,'(a,es15.5)') 'IG_rho test 2', IG_den
  IG_den = IG_den / scale_d                          !----- IG_rho change units from g cm^-3 to internal
  ! write(*,'(a,2es15.5)') 'IG_rho test 3', IG_den, smallr
  IG_den = max(IG_den,smallr)    
  IG_TFl = IG_T2/scale_T2/(gamma-1)
  ! print*, 'IG_den, IG_TFl ',IG_den, IG_TFl
  ! print*, 'IG_den2, IG_TFl2 ',IG_rho/scale_nH, IG_T2/scale_T2/(gamma-1)

  ! Set uold to namelist values for myid cells
  ! do ind=1,twotondim
  !    iskip=ncoarse+(ind-1)*ngridmax

  !    do ivar=nvar,1,-1
  !       do i=1,active(ilevel)%ngrid

  !          if(uold(active(ilevel)%igrid(i)+iskip,1).lt.IG_rho/scale_nH) then
  !             uold(active(ilevel)%igrid(i)+iskip,ivar)    = 0D0
  !             if(ivar.eq.1) then
  !              ! uold(active(ilevel)%igrid(i)+iskip,ivar)   = max(IG_rho/scale_nH,smallr)
  !              uold(active(ilevel)%igrid(i)+iskip,ivar)   = IG_den
  !             end if
  !             if(ivar.eq.ndim+2)then
  !                ! uold(active(ilevel)%igrid(i)+iskip,ivar) = IG_T2/scale_T2/(gamma-1)*max(IG_rho/scale_nH,smallr)
  !                uold(active(ilevel)%igrid(i)+iskip,ivar) = IG_TFl*IG_den
  !             endif
  !             if(metal) then
  !               ! if(ivar.eq.imetal  ) uold(active(ilevel)%igrid(i)+iskip,ivar) = max(IG_rho/scale_nH,smallr)*IG_metal * 0.333d0
  !               ! if(ivar.eq.imetal+1) uold(active(ilevel)%igrid(i)+iskip,ivar) = max(IG_rho/scale_nH,smallr)*IG_metal * 0.667d0
  !               if(ivar.eq.imetal  ) uold(active(ilevel)%igrid(i)+iskip,ivar) = IG_den*IG_metal * 0.333d0
  !               if(ivar.eq.imetal+1) uold(active(ilevel)%igrid(i)+iskip,ivar) = IG_rho*IG_metal * 0.667d0
  !             endif
  !          endif
  !       end do
  !    end do
  ! end do
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        if(uold(active(ilevel)%igrid(i)+iskip,1).lt.IG_den) then

           uold(active(ilevel)%igrid(i)+iskip,:) = 0D0           
           uold(active(ilevel)%igrid(i)+iskip,1) = IG_den
           uold(active(ilevel)%igrid(i)+iskip,5) = IG_TFl*IG_den
           if(metal) then
             if(ivar.eq.imetal  ) uold(active(ilevel)%igrid(i)+iskip,ivar) = IG_den*IG_metal * 0.333d0
             !if(ivar.eq.imetal+1) uold(active(ilevel)%igrid(i)+iskip,ivar) = IG_den*IG_metal * 0.667d0
           endif
        endif
     end do
  end do
  ! Set cell averaged kinetic energy
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        ! Initialisation of the refinement mask
        if((ic_mask_ivar.gt.0).and.(ivar_refine.gt.0).and.(ic_mask_ivar.le.nvar).and.(ivar_refine.le.nvar))then
     ! Switch to K/mu for ic_mask_ivar=ndim+2 case
           if(ic_mask_ivar.eq.ndim+2) then
        u = uold(active(ilevel)%igrid(i)+iskip,ic_mask_ivar)*scale_T2*(gamma-1)
           else
        u = uold(active(ilevel)%igrid(i)+iskip,ic_mask_ivar)
           endif
           if(ic_mask_ivar.gt.1)then
              u = u/uold(active(ilevel)%igrid(i)+iskip,1)
           endif
           if((u.ge.ic_mask_min).and.(u.le.ic_mask_max))then
              uold(active(ilevel)%igrid(i)+iskip,ivar_refine) = 1.0*uold(active(ilevel)%igrid(i)+iskip,1)
           endif
        endif
        e = 0d0
        do idim=1,ndim
           e = e+0.5*uold(active(ilevel)%igrid(i)+iskip,idim+1)**2/uold(active(ilevel)%igrid(i)+iskip,1)
        enddo
        uold(active(ilevel)%igrid(i)+iskip,ndim+2) = uold(active(ilevel)%igrid(i)+iskip,ndim+2)+e
     end do
  end do

!#ifdef SOLVERmhd
!  ! set constant magnetic field
!  CALL mag_constant(ilevel)
!  ! toroidal field
!  CALL mag_compute(ilevel)
!#endif

  ! Set uold to 0 for virtual boundary cells
  do icpu=1,ncpu
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,nvar
        do i=1,reception(icpu,ilevel)%ngrid
           uold(reception(icpu,ilevel)%igrid(i)+iskip,ivar)=0.0
        end do
     end do
  end do
  end do

111 format('   Entering init_uold for level ',i2)

end subroutine init_uold

subroutine condinit_loc(ilevel,rfoverride)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use dice_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------
  ! This routine computes the initial density field at level ilevel using
  ! the CIC scheme from particles that are not entirely in
  ! level ilevel (boundary particles).
  ! Arrays flag1 and flag2 are used as temporary work space.
  !------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,idim,icpu,next_part
  integer::i,ig,ip,npart1,npart2
  real(dp)::dx
  logical,intent(in)::rfoverride

  integer,dimension(1:nvector),save::ind_grid,ind_cell
  integer,dimension(1:nvector),save::ind_part,ind_grid_part
  real(dp),dimension(1:nvector,1:ndim),save::x0

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  ! Loop over cpus
  do icpu=1,ncpu
    ! Loop over grids
    igrid=headl(icpu,ilevel)
    ig=0
    ip=0
    do jgrid=1,numbl(icpu,ilevel)
      npart1=numbp(igrid)  ! Number of particles in the grid
      npart2=0
      ! Count gas particles
      if(npart1>0)then
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
          ! Save next particle   <--- Very important !!!
          next_part=nextp(ipart)
          if(ic_mask_ptype.eq.-1)then
            if(idp(ipart).eq.1)then
              npart2=npart2+1
            end if
          else
            npart2=npart2+1
          endif
          ipart=next_part  ! Go to next particle
        end do
      end if

      ! Gather gas particles∂
      if(npart2>0)then
        ig=ig+1
        ind_grid(ig)=igrid
        ipart=headp(igrid)

        ! Loop over particles
        do jpart=1,npart1
          ! Save next particle   <--- Very important !!!
          next_part=nextp(ipart)
          if(ic_mask_ptype.eq.-1)then
            if(idp(ipart).eq.1)then
              if(ig==0)then
                ig=1
                ind_grid(ig)=igrid
              end if
              ip=ip+1
              ind_part(ip)=ipart
              ind_grid_part(ip)=ig
            endif
          else
            if(ig==0)then
              ig=1
              ind_grid(ig)=igrid
            end if
            ip=ip+1
            ind_part(ip)=ipart
            ind_grid_part(ip)=ig
          endif
          if(ip==nvector)then
            ! Lower left corner of 3x3x3 grid-cube
            do idim=1,ndim
              do i=1,ig
                x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
              end do
            end do

            do i=1,ig
              ind_cell(i)=father(ind_grid(i))
            end do
            
            !if(amr_struct) then
            !  call init_gas_ngp(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
            !else
              call init_gas_cic(ind_cell,ind_part,ind_grid_part,x0,ig,ip,ilevel,rfoverride)
            !endif
            
            ip=0
            ig=0
          end if
          ipart=next_part  ! Go to next particle
        end do
        ! End loop over particles
      end if

      igrid=next(igrid)   ! Go to next grid
    end do
    ! End loop over grids

    if(ip>0)then
      ! Lower left corner of 3x3x3 grid-cube
      do idim=1,ndim
        do i=1,ig
          x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
        end do
      end do
      do i=1,ig
        ind_cell(i)=father(ind_grid(i))
      end do
      !if(amr_struct) then
      !  call init_gas_ngp(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,rfoverride)
      !else
        call init_gas_cic(ind_cell,ind_part,ind_grid_part,x0,ig,ip,ilevel,rfoverride,rfoverride)
      !endif
    end if
  end do

111 format('   Entering condinit_loc for level ',I2)

end subroutine condinit_loc
!==================================================================================
!==================================================================================
!==================================================================================

subroutine init_gas_cic(ind_cell,ind_part,ind_grid_part,x0,ng,np,ilevel,rfoverride)
  use amr_commons
  use pm_commons
  use hydro_commons
  use dice_commons
  use cooling_module
  implicit none
  integer::ng,np,ilevel
  integer ,dimension(1:nvector)::ind_cell,ind_grid_part,ind_part
  real(dp),dimension(1:nvector,1:ndim)::x0
  !------------------------------------------------------------------
  ! This routine computes the initial density field at level ilevel using
  ! the CIC scheme. Only cells that are in level ilevel
  ! are updated by the input particle list.
  !------------------------------------------------------------------
  character(len=128)::testfile
  logical::error,file_exists
  integer::j,ind,idim,nx_loc,kk
  real(dp)::dx,dx_loc,scale
  logical,intent(in)::rfoverride
  ! Grid-based arrays
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle-based arrays
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector,1:ndim),save::xx,dd,dg
  integer ,dimension(1:nvector,1:ndim),save::ig,id,igg,igd,icg,icd
  real(dp),dimension(1:nvector,1:twotondim),save::vol
  integer ,dimension(1:nvector,1:twotondim),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:nvector),save::ethermal
  real(dp),dimension(1:nvector),save::vol_loc
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  ! print*, 'BOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO',MYID
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
   ! write(*,'(a,i,3es15.5,i)') 'check8',myid, sum(up32),minval(up32),maxval(up32),(size(up32))

  ! write(*,'(a,i,3es15.5)') 'check6',myid, sum(up),minval(up),maxval(up)
  ! write(*,'(a,i,3es15.5)') 'check7',myid, sum(up32),minval(up32),maxval(up32)
  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc(1:nvector)=dx_loc**ndim
  
  ! print*, 'Gas_cic', myid, shape(readflag)
  ! Gather neighboring father cells (should be present anytime !)
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Rescale particle position at level ilevel
  do idim=1,ndim
    do j=1,np
      xx(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
    end do
  end do
  do idim=1,ndim
    do j=1,np
      xx(j,idim)=xx(j,idim)-x0(ind_grid_part(j),idim)
    end do
  end do
  do idim=1,ndim
    do j=1,np
      xx(j,idim)=xx(j,idim)/dx
    end do
  end do

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
    do j=1,np
      if(xx(j,idim)<0.5D0.or.xx(j,idim)>5.5D0)error=.true.
    end do
  end do
  if(error)then
    write(*,*)'problem in cic'
    do idim=1,ndim
      do j=1,np
        if(xx(j,idim)<0.5D0.or.xx(j,idim)>5.5D0)then
          write(*,*)xx(j,1:ndim)
        endif
      end do
    end do
    stop
  end if

  ! CIC at level ilevel (dd: right cloud boundary; dg: left cloud boundary)
  do idim=1,ndim
    do j=1,np
      dd(j,idim)=xx(j,idim)+0.5D0
      id(j,idim)=dd(j,idim)
      dd(j,idim)=dd(j,idim)-id(j,idim)
      dg(j,idim)=1.0D0-dd(j,idim)
      ig(j,idim)=id(j,idim)-1
    end do
  end do

  ! Compute cloud volumes
#if NDIM==1
  do j=1,np
    vol(j,1)=dg(j,1)
    vol(j,2)=dd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
    vol(j,1)=dg(j,1)*dg(j,2)
    vol(j,2)=dd(j,1)*dg(j,2)
    vol(j,3)=dg(j,1)*dd(j,2)
    vol(j,4)=dd(j,1)*dd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
    vol(j,1)=dg(j,1)*dg(j,2)*dg(j,3)
    vol(j,2)=dd(j,1)*dg(j,2)*dg(j,3)
    vol(j,3)=dg(j,1)*dd(j,2)*dg(j,3)
    vol(j,4)=dd(j,1)*dd(j,2)*dg(j,3)
    vol(j,5)=dg(j,1)*dg(j,2)*dd(j,3)
    vol(j,6)=dd(j,1)*dg(j,2)*dd(j,3)
    vol(j,7)=dg(j,1)*dd(j,2)*dd(j,3)
    vol(j,8)=dd(j,1)*dd(j,2)*dd(j,3)
  end do
#endif

  ! Compute parent grids
  do idim=1,ndim
    do j=1,np
      igg(j,idim)=ig(j,idim)/2
      igd(j,idim)=id(j,idim)/2
    end do
  end do

#if NDIM==1
  do j=1,np
    kg(j,1)=1+igg(j,1)
    kg(j,2)=1+igd(j,1)
  end do
#endif

#if NDIM==2
  do j=1,np
    kg(j,1)=1+igg(j,1)+3*igg(j,2)
    kg(j,2)=1+igd(j,1)+3*igg(j,2)
    kg(j,3)=1+igg(j,1)+3*igd(j,2)
    kg(j,4)=1+igd(j,1)+3*igd(j,2)
  end do
#endif

#if NDIM==3
  do j=1,np
    kg(j,1)=1+igg(j,1)+3*igg(j,2)+9*igg(j,3)
    kg(j,2)=1+igd(j,1)+3*igg(j,2)+9*igg(j,3)
    kg(j,3)=1+igg(j,1)+3*igd(j,2)+9*igg(j,3)
    kg(j,4)=1+igd(j,1)+3*igd(j,2)+9*igg(j,3)
    kg(j,5)=1+igg(j,1)+3*igg(j,2)+9*igd(j,3)
    kg(j,6)=1+igd(j,1)+3*igg(j,2)+9*igd(j,3)
    kg(j,7)=1+igg(j,1)+3*igd(j,2)+9*igd(j,3)
    kg(j,8)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
#endif

  do ind=1,twotondim
    do j=1,np
      if(nbors_father_cells(ind_grid_part(j),kg(j,ind)).gt.0) then
        igrid(j,ind)=son(nbors_father_cells(ind_grid_part(j),kg(j,ind)))
      else
        igrid(j,ind)=0
      endif
    end do
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
      icg(j,idim)=ig(j,idim)-2*igg(j,idim)
      icd(j,idim)=id(j,idim)-2*igd(j,idim)
     end do
  end do

#if NDIM==1
  do j=1,np
    icell(j,1)=1+icg(j,1)
    icell(j,2)=1+icd(j,1)
  end do
#endif

#if NDIM==2
  do j=1,np
    icell(j,1)=1+icg(j,1)+2*icg(j,2)
    icell(j,2)=1+icd(j,1)+2*icg(j,2)
    icell(j,3)=1+icg(j,1)+2*icd(j,2)
    icell(j,4)=1+icd(j,1)+2*icd(j,2)
  end do
#endif

#if NDIM==3
  do j=1,np
    icell(j,1)=1+icg(j,1)+2*icg(j,2)+4*icg(j,3)
    icell(j,2)=1+icd(j,1)+2*icg(j,2)+4*icg(j,3)
    icell(j,3)=1+icg(j,1)+2*icd(j,2)+4*icg(j,3)
    icell(j,4)=1+icd(j,1)+2*icd(j,2)+4*icg(j,3)
    icell(j,5)=1+icg(j,1)+2*icg(j,2)+4*icd(j,3)
    icell(j,6)=1+icd(j,1)+2*icg(j,2)+4*icd(j,3)
    icell(j,7)=1+icg(j,1)+2*icd(j,2)+4*icd(j,3)
    icell(j,8)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
  end do
#endif

  ! Update mass density and number density fields
  do ind=1,twotondim

    ! Check if particles are entirely in level ilevel
    do j=1,np
      ok(j)=igrid(j,ind)>0
    end do

    ! Compute parent cell adress
    do j=1,np
      if(ok(j))then
        indp(j,ind)=ncoarse+(icell(j,ind)-1)*ngridmax+igrid(j,ind)
      end if
    end do

    do j=1,np
      ok(j)=ok(j).and.(idp(ind_part(j)).eq.1)
    end do
  !=========================================================================!
  !=========================================================================!
  !=========================================================================!
  !=========================================================================!
  !=========================================================================!
  !=========================================================================!
  !=========================================================================!
  !=========================================================================!
  !=========================================================================!
  !=========================================================================!
  !=========================================================================!
  !=========================================================================!
  !---------------------- DICE ADDING GAS TO AMR GRID ----------------------!
    
    ! Update hydro variablesm
    do j=1,np

      if(ok(j)) then
        ! Specific kinetic energy of the gas particle
        if (mp(ind_part(j)).ne.0.0d0 .and. tp(ind_part(j)).eq.0.0d0) print*, 'somwthing is wrong',tp(ind_part(j)), mp(ind_part(j))
        ! if (mp(ind_part(j)).eq.0.0d0) write(*,*) 'gas particle mass is wrong ', mp(ind_part(j))
        ! if (tp(ind_part(j)).le.0.0d0) write(*,*) 'gas particle temp is wrong ', tp(ind_part(j))
        ethermal(j) = tp(ind_part(j))
        kk = indp(j,ind)
        ! if (dmf_i.eq.3) print*, xp(j,:)
        ! print*, xx(j,:)*dx

        ! Update hydro variable in CIC cells

        if (dmf_i.gt.1) then
          if (readflag(kk).eq.1) then
            ! ptmp(j,:)  = uold(kk,:)
            ! utmp(kk,:) = uold(kk,:)
            readflag(kk) = 3
          end if
        end if


        ! ifrfoverride .and. readflag(kk).eq.3) the  = uold(,:) + utmp(kk,:)
        !   if ((j,1).ne.0.0d0) uold(kk,:) = uold(kk,:) + ptmp(j,:)
        ! end if

        uold(kk,1)=uold(kk,1)+mp(ind_part(j))*vol(j,ind)/vol_loc(j)
        if (readflag(kk).eq.0) readflag(kk) = 2 !dmf_i

        !<<<<<----- Set Gas Velocity ----->>>>>!
        do idim=1,ndim
          uold(kk,idim+1)=uold(kk,idim+1) + &
          (mp(ind_part(j))*vol(j,ind) / vol_loc(j))*(vp(ind_part(j),idim))
        end do

        !<<<<<----- Set Gas Temperature ----->>>>>!
        uold(kk,ndim+2)=uold(kk,ndim+2) + &
        ((mp(ind_part(j))*vol(j,ind)) / vol_loc(j)) * (ethermal(j) / scale_T2)


        !<<<<<----- Set Gas Metals ----->>>>>!
        if(metal) then
          ! !<<<<<----- KMG Modification ----->>>>>!
          ! !----- Same metals as OA's model. Assumes Z~ 2O _ 1 Fe (Madau)
          ! !q(i,ndim+3)=0.333*Zgas*0.02d0  !Iron
          uold(kk,imetal) = uold(kk,imetal) + &
          mp(ind_part(j)) * (vol(j,ind) / vol_loc(j)) * (zp(ind_part(j)))! * 0.333) !----- Iron

          !uold(kk,imetal+1) = uold(kk,imetal+1) + &
          !mp(ind_part(j)) * (vol(j,ind) / vol_loc(j)) * (zp2(ind_part(j)))! * 0.667) !----- Oxygen
          ! !<<<<<----- ---------------- ----->>>>>!         
        endif
      endif
          ! Update passive scalar mask
      if(ic_mask_ptype.gt.-1) then
        print*, 'boo'
        uold(indp(j,ind),ivar_refine)=uold(indp(j,ind),ivar_refine) + &
        mp(ind_part(j))*vol(j,ind)/vol_loc(j)*maskp(ind_part(j))
      endif
    end do
    if (dmf_i.eq.3) close(52)
  end do
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!
!=========================================================================!

end subroutine init_gas_cic
#endif

