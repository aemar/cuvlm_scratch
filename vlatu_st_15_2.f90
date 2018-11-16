program vlatunsteady
  use input
  use basictypes
  use maths
  use rotation
  use vtxparticle
  use pg
  use infl
!  use fourdgridtypes
!  use fourdgridsf
  use hyperbolic
  use functions
  use penetration
  implicit none

  integer :: INFO
  real*8,dimension(:,:),allocatable :: wingcirc,COEFFS,panvind

  real*8,dimension(:),allocatable :: RHS,lift,drag,d1ra1
  integer,dimension(:),allocatable :: IPIV
  integer,dimension(:,:),allocatable :: wakesrcpan,wakeneigh,panwakesrc
  real*8,dimension(:,:,:),allocatable :: wakeind,wingind,cmat,panfstream,uwind,duwind
  real*8,dimension(:,:,:),allocatable :: panforce
  real*8,dimension(:,:),allocatable :: pancp
  logical,dimension(:,:),allocatable :: wingprcssd
  logical,dimension(:,:),allocatable :: wakeprcssd

  real*8,dimension(:,:),allocatable :: dummys11,dummys12,dummys21,dummys22,dummys31,dummys32

  real*8,dimension(3) :: d3r,d3r1,d3r2,d3r3,d3r4,d3r5,d3r6
  real*8,dimension(3,3) :: dummy33r
  real*8,dimension(3,4) :: dummy34r
  real*8,dimension(3),parameter :: zerovec = (/0.,0.,0./)
  real*8,parameter :: runity = 1.
  real*8,parameter :: rzero = 0.
  integer,parameter :: izero = 0
  integer,dimension(3) :: dummy3i1
  integer :: nvtx,d1i1,d1i2,dummy1i2
  integer,dimension(4) :: dummy4i1
  real*8,dimension(4) :: d4r1

  integer,dimension(1) :: d1ia1

  logical :: check3

  integer :: nperim,ioin,ioout
  integer :: i,j,k,l,q,count1
  real*8 :: bstabprec,dummy1r1,dummy1r2,dummy1r3,f1,f2,totalA,te,d1r1
  character(len=1),parameter :: cQ = 'Q'
  logical,parameter :: ltrue = .true.
  logical,parameter :: lfalse = .false.

  integer :: dvtx,fvtx,nlevel

  real*8,dimension(:,:,:),allocatable :: acin,inac

  real :: timeon,timeoff

!############# OUTPUT ##############
  character(len=4) :: istring
  character(len=4) :: kstring
!###################################

!###################################
!###################################

  type(point),dimension(:),allocatable,target :: wingpnt
  type(point),dimension(:),allocatable,target :: wakepnt
  type(point),dimension(:),allocatable,target :: vtxpnt
  type(panel),dimension(:),allocatable :: wingpan
  type(panel),dimension(:),allocatable :: wakepan

  type(point), pointer :: pntpntr1, pntpntr2

  logical,dimension(:),allocatable :: dummyl1
  integer,dimension(:),allocatable :: dummyi1
  integer,dimension(:),allocatable :: wkpntsrc

  real*8,parameter :: qtsf = 0.25

  integer,dimension(:),allocatable :: itestarray

  real(kind=8),dimension(12) :: tvec

  real*8 :: lift_seg,area

!###################################
!###################################

  call cpu_time(timeon)

!-----------------------------------------------------------------------------------------------------------------------------------

  open (unit=12, file=trim('input/kosc.dat'), status='old', action='read', iostat=ioin)
    if (ioin == 0) then
      read (12,*), Mosc
      read (12,*), krf
      close(12)
    end if


  !read parameters
  open (unit=10, file=trim('input/vlatinit.dat'), status='old', action='read', iostat=ioin)
    read (10,*), trigger1
    read (10,*), trigger2
    read (10,*), density
    read (10,*), tsteps
    read (10,*), deltat
    read (10,*), uinf
    read (10,*), ntraj
    read (10,*)
    read (10,*), nwpanels
  close(10)

  allocate(lift(tsteps),drag(tsteps))

  !check enough timesteps for two wakepanels.
  if (tsteps < nwpanels+2) then
    print*, 'Cannot run for', tsteps, ',tsteps set to', nwpanels+2 
    tsteps = nwpanels+2
  end if

  !read panels
  np = 0
  npan = 0
  do i=1,ntraj,1
    write(istring,'(i4.4)'), i
    open(unit=11,file=trim('input/'//'vlatin'//trim(istring)//'.dat'),status='old',action='read',iostat=ioin)
    read(11,*) d1i1
    read(11,*) dummy1i2
    np = np + d1i1
    npan = npan + dummy1i2
    if (i==1) then
      npan_wing = dummy1i2
    end if
    close(11)
  end do

  allocate(wingpnt(np),wingpan(npan),d1ra1(npan))
  do i=1,npan,1
    allocate(wingpan(i)%cinfr(3,tsteps),wingpan(i)%dcirc(4,tsteps),wingpan(i)%circ(tsteps))
  end do
  do i=1,np,1
    allocate(wingpnt(i)%set(tsteps))
  end do

  np = 0
  npan = 0

  do k=1,ntraj,1
    write(istring,'(i4.4)'), k
print*, trim('input/'//'vlatin'//trim(istring)//'.dat')
    open(unit=11,file=trim('input/'//'vlatin'//trim(istring)//'.dat'),status='old',action='read',iostat=ioin)
    read(11,*) d1i1
    read(11,*) dummy1i2
    do i=np+1,np+d1i1,1
      read(11,*) wingpnt(i)%acfr
      wingpnt(i)%trj = k
      wingpnt(i)%tag = i
    end do
    do i=npan+1,npan+dummy1i2,1
      read(11,*) wingpan(i)%ptype,dummy4i1
      do j=1,4,1
        wingpan(i)%pnt(j)%p => wingpnt(dummy4i1(j) + np)
      end do
      wingpan(i)%trj = k
    end do
    do i=npan+1,npan+dummy1i2,1
      read(11,*) dummy4i1
      do j=1,4,1
        if (dummy4i1(j)>0) then
          wingpan(i)%neigh(j) = dummy4i1(j) + npan
        else
          wingpan(i)%neigh(j) = dummy4i1(j)
        end if
      end do
    end do
    close(11)
    np = np + d1i1
    npan = npan + dummy1i2
  end do

!  do i=1,np,1
!!    print*, wingpnt(i)%trj, wingpnt(i)%acfr
!    if (wingpnt(i)%trj == 2) then
!      if (wingpnt(i)%acfr(2) == 0. .and. wingpnt(i)%acfr(3) == 0.) then
!        print*, wingpnt(i)%trj, wingpnt(i)%acfr
!      end if
!    end if
!  end do

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

  !find smallest and largest dimension of wingpnts
  dummy1r1 = 0.
  dummy1r2 = 0.
  do i=1,np,1
    do j=1,3,1
      if (wingpnt(i)%acfr(j) < dummy1r1) then
        dummy1r1 = wingpnt(i)%acfr(j)
      end if
      if (wingpnt(i)%acfr(j) > dummy1r2) then
        dummy1r2 = wingpnt(i)%acfr(j)
      end if
    end do
  end do

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

  !determine which panels are wakepanels and which side
  allocate(panwakesrc(npan,4))
  panwakesrc = -1
  nwakepan = 0
  do i=1,npan,1
    do j=1,4,1
      if (wingpan(i)%neigh(j) == 0)then
        nwakepan = nwakepan + 1
        panwakesrc(i,j) = nwakepan
      end if
    end do
  end do

  allocate(wakesrcpan(nwakepan,2),wakeind(3,npan,tsteps),wingind(3,npan,tsteps),wingcirc(npan,tsteps))
  !need IPIV if using LU decomp from lapack
  allocate(COEFFS(npan,npan),RHS(npan),IPIV(npan))

  count1 = 0
  do i=1,npan,1
    do j=1,4,1
      if (wingpan(i)%neigh(j) == 0)then
        count1 = count1+1
        wakesrcpan(count1,1) = i
        wakesrcpan(count1,2) = j
      end if
    end do
  end do

  allocate(wakeneigh(nwakepan,4))
  wakeneigh = -1
  count1 = 0
  do i=1,npan,1
    do j=1,4,1
      if (wingpan(i)%neigh(j) == 0)then
        count1 = count1+1 !out of nwakepan
        if (wingpan(i)%neigh(cycperm(4,1,j)) > 0) then
          if (wingpan(wingpan(i)%neigh(cycperm(4,1,j)))%neigh(j) == 0) then
            wakeneigh(count1,cycperm(4,1,j)) = panwakesrc(wingpan(i)%neigh(cycperm(4,1,j)),j)
          else
            wakeneigh(count1,cycperm(4,1,j)) = 0
          end if
        else if (wingpan(i)%neigh(cycperm(4,1,j)) < 0) then
          wakeneigh(count1,cycperm(4,1,j)) = 0
        else if (wingpan(i)%neigh(cycperm(4,1,j)) == 0) then
          wakeneigh(count1,cycperm(4,1,j)) = panwakesrc(i,cycperm(4,1,j))
        end if
        if (wingpan(i)%neigh(cycperm(4,3,j)) > 0) then
          if (wingpan(wingpan(i)%neigh(cycperm(4,3,j)))%neigh(j) == 0) then
            wakeneigh(count1,cycperm(4,3,j)) = panwakesrc(wingpan(i)%neigh(cycperm(4,3,j)),j)
          else
            wakeneigh(count1,cycperm(4,3,j)) = 0
          end if
        else if (wingpan(i)%neigh(cycperm(4,3,j)) < 0) then
          wakeneigh(count1,cycperm(4,3,j)) = 0
        else if (wingpan(i)%neigh(cycperm(4,3,j)) == 0) then
          wakeneigh(count1,cycperm(4,3,j)) = panwakesrc(i,cycperm(4,3,j))
        end if
        !same wakepanel from next time step (in front)
        wakeneigh(count1,cycperm(4,2,j)) = -5
        !same wakepanel from previous time step (further back)
        wakeneigh(count1,j) = -6
      end if
    end do
  end do
  !neigh = -2 means wake cutoff, vortex horseshoe

!-----------------------------------------------------------------------------------------------------------------------------------
  allocate(wakepan((nwpanels+2)*nwakepan))
  do i=1,(nwpanels+2)*nwakepan,1
    allocate(wakepan(i)%circ(tsteps),wakepan(i)%dcirc(4,tsteps))
    wakepan(i)%circ = 0.
    wakepan(i)%dcirc = 0.
  end do

  do j=1,nwakepan,1
    do k=1,4,1
      do i=2,nwpanels+1,1
        d1i1 = (i-1)*nwakepan + j
        if (wakeneigh(j,k)>0) then
          wakepan(d1i1)%neigh(k) = wakeneigh(j,k) + (i-1)*nwakepan
        else if (wakeneigh(j,k)==-6) then
          wakepan(d1i1)%neigh(k) = d1i1 + nwakepan
        else if (wakeneigh(j,k)==-5) then
          wakepan(d1i1)%neigh(k) = d1i1 - nwakepan
        else if (wakeneigh(j,k)==0) then
          wakepan(d1i1)%neigh(k) = 0
        end if
      end do

      d1i1 = j
      if (wakeneigh(j,k)>0) then
        wakepan(d1i1)%neigh(k) = wakeneigh(j,k)
      else if (wakeneigh(j,k)==-6) then
        wakepan(d1i1)%neigh(k) = d1i1 + nwakepan
      else if (wakeneigh(j,k)==-5) then
        wakepan(d1i1)%neigh(k) = -5
      else if (wakeneigh(j,k)==0) then
        wakepan(d1i1)%neigh(k) = 0
      end if

      d1i1 = (nwpanels+1)*nwakepan + j
      if (wakeneigh(j,k)>0) then
        wakepan(d1i1)%neigh(k) = wakeneigh(j,k) + (nwpanels+1)*nwakepan
      else if (wakeneigh(j,k)==-6) then
        wakepan(d1i1)%neigh(k) = -6
      else if (wakeneigh(j,k)==-5) then
        wakepan(d1i1)%neigh(k) = d1i1 - nwakepan
      else if (wakeneigh(j,k)==0) then
        wakepan(d1i1)%neigh(k) = 0
      end if

    end do
  end do

!do i=1,nwpanels+2,1
!do j=1,nwakepan,1
!d1i1 = (i-1)*nwakepan + j
!print*, wakepan(d1i1)%neigh
!end do
!end do
 
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------   

  allocate(dummyl1(np))
  dummyl1 = .false.

  !determine which wing points become wakepoints
  do i=1,nwakepan,1
    dummyl1(wingpan(wakesrcpan(i,1))%pnt(wakesrcpan(i,2))%p%tag) = .true.
    dummyl1(wingpan(wakesrcpan(i,1))%pnt(cycperm(4,1,wakesrcpan(i,2)))%p%tag) = .true.
  end do

  !count number of wakepoints
  nwkpnt = 0
  do i=1,np,1
    if (dummyl1(i)) then
      nwkpnt = nwkpnt + 1
    end if
  end do

  !allocate wakepnt array and set %tag
  allocate(wakepnt(nwkpnt*(nwpanels+2)))
  do i=1,nwpanels+2,1
    do j=1,nwkpnt,1
      d1i1 = (i-1)*nwkpnt + j
      allocate(wakepnt(d1i1)%set(i:tsteps))
      wakepnt(d1i1)%tag = d1i1
      do k=i,tsteps,1
        wakepnt(d1i1)%set(k)%infr(4) = dble(k)
      end do
    end do
  end do

  !determine connectivity from wingpoints to wakepoints and vice versa
  allocate(dummyi1(np))
  allocate(wkpntsrc(nwkpnt))
  count1 = 0
  do i=1,np,1
    if (dummyl1(i)) then
      count1 = count1 + 1
      wkpntsrc(count1) = i
      dummyi1(i) = count1
    end if
  end do

  !pointer association of first line of wakepnts, which are wingpnts 
  do i=1,nwakepan,1
    !pos2
    wakepan(i)%pnt(cycperm(4,3,wakesrcpan(i,2)))%p => wingpan(wakesrcpan(i,1))%pnt(wakesrcpan(i,2))%p
    !pos1
    wakepan(i)%pnt(cycperm(4,2,wakesrcpan(i,2)))%p => wingpan(wakesrcpan(i,1))%pnt(cycperm(4,1,wakesrcpan(i,2)))%p
  end do
  !do the rear line, which are wakepnts
  do i=1,nwakepan,1
    wakepan(i)%pnt(wakesrcpan(i,2))%p => wakepnt(dummyi1(wingpan(wakesrcpan(i,1))%pnt(wakesrcpan(i,2))%p%tag))
    wakepan(i)%pnt(cycperm(4,1,wakesrcpan(i,2)))%p => &
      & wakepnt(dummyi1(wingpan(wakesrcpan(i,1))%pnt(cycperm(4,1,wakesrcpan(i,2)))%p%tag))
  end do

  !pointer association of all other lines, all wakepnts
  do j=2,nwpanels+2,1
    do i=1,nwakepan,1
      d1i1 = (j-1)*nwakepan + i
      dummy1i2 = (j-2)*nwakepan + i
      wakepan(d1i1)%pnt(cycperm(4,3,wakesrcpan(i,2)))%p => wakepan(dummy1i2)%pnt(wakesrcpan(i,2))%p
      wakepan(d1i1)%pnt(cycperm(4,2,wakesrcpan(i,2)))%p => wakepan(dummy1i2)%pnt(cycperm(4,1,wakesrcpan(i,2)))%p
      wakepan(d1i1)%pnt(cycperm(4,1,wakesrcpan(i,2)))%p => &
        & wakepnt(wakepan(dummy1i2)%pnt(cycperm(4,1,wakesrcpan(i,2)))%p%tag+nwkpnt)
      wakepan(d1i1)%pnt(wakesrcpan(i,2))%p => wakepnt(wakepan(dummy1i2)%pnt(wakesrcpan(i,2))%p%tag+nwkpnt)
    end do
  end do

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

  allocate(cmat(3,npan,npan),panvind(3,npan),uwind(3,npan,tsteps),duwind(3,npan,tsteps),panfstream(3,npan,tsteps))
  allocate(panforce(3,npan,tsteps),pancp(npan,tsteps))

  call panelparam1(npan,wingpan)

  totalA = 0.
  do i=1,npan,1
    totalA = totalA + wingpan(i)%area
  end do

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

  nvtx = tsteps - (nwpanels+2)
  fvtx = npan + (nwpanels+2)*nwakepan
  dvtx = fvtx + nvtx*nwakepan

  allocate(vtxpnt(dvtx))
  !save memory by only allocating non-zero time
  do i=1,npan,1
    allocate(vtxpnt(i)%set(tsteps))
    do j=1,tsteps,1
      vtxpnt(i)%set(j)%act = .false.
      vtxpnt(i)%set(j)%infr(1:3) = 0.
      vtxpnt(i)%set(j)%infr(4) = dble(j)
      vtxpnt(i)%set(j)%minfr = 0.
    end do
  end do
  do i=1,tsteps,1
    do j=1,nwakepan,1
      d1i1 = npan + (i-1)*nwakepan + j
      allocate(vtxpnt(d1i1)%set(i:tsteps))
      do k=i,tsteps,1
        vtxpnt(d1i1)%set(k)%act = .false.
        vtxpnt(d1i1)%set(k)%infr(1:3) = 0.
        vtxpnt(d1i1)%set(k)%infr(4) = dble(k)
        vtxpnt(d1i1)%set(k)%minfr = 0.
      end do
    end do
  end do

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

  wakeind = 0.
  d3r = 0.
  wingind = 0.
  uwind = 0.
  duwind = 0.

  allocate(dummys11(3,nwakepan),dummys12(3,nwakepan),dummys21(3,npan),dummys22(3,npan), &
    & dummys31(3,(nwpanels+2)*nwakepan),dummys32(3,(nwpanels+2)*nwakepan))

  allocate(wingprcssd(4,npan),wakeprcssd(4,(nwpanels+2)*nwakepan))

  wingcirc = 0.
  do i=1,npan,1
    wingpan(i)%circ = 0.
    wingpan(i)%dcirc(:,:) = 0.
  end do

!  allocate(itestarray(npan))

!***********************************************************************************************************************************
!***********************************************************************************************************************************
!***********************************************************************************************************************************
!***********************************************************************************************************************************

  open (unit=39,file='output/rundata.dat',status='replace',action='write',iostat=ioout)

  !start time stepping
  do i=1,tsteps,1
    print*
    print*, 't =', i
    write(39,*)
    write(39,*) 't =', i

    cmat = 0.
    do j=1,npan,1
      wingpan(j)%vind = 0.
    end do

    do j=1,min(i,(nwpanels+2))*nwakepan,1
      do k=1,4,1
        wakepan(j)%pnt(k)%p%vind = 0.
      end do
    end do

    if (i>nwpanels+2) then
      do j=fvtx+1,fvtx+(i-(nwpanels+2))*nwakepan,1
        vtxpnt(j)%gradv = 0.
        vtxpnt(j)%vind = 0.
      end do
    end if

!    do j=1,ntraj,1
!      call rotmatrices(traj(j)%step(i)%t,acin(j,:,:),inac(j,:,:))
!    end do

    do j=1,np,1
!print*, wingpan(j)%trj, i*deltat, wingpnt(j)%acfr
      call traj(wingpnt(j)%trj, i*deltat, wingpnt(j)%acfr, tvec)
      wingpnt(j)%set(i)%infr(1:3) = tvec(1:3)
      wingpnt(j)%set(i)%infr(4) = dble(i)
!print*, wingpnt(j)%set(i)%infr(1:3)
    end do

    do j=1,npan,1
      wingpan(j)%cinfr(:,i) = 0.25*(wingpan(j)%pnt(1)%p%set(i)%infr(1:3)+wingpan(j)%pnt(2)%p%set(i)%infr(1:3)+ &
        & wingpan(j)%pnt(3)%p%set(i)%infr(1:3)+wingpan(j)%pnt(4)%p%set(i)%infr(1:3))
    end do
!    do j=1,npan,1
!      if (wingpan(j)%neigh(3) == 0) then
!        wingpan(j)%cinfr(:,i) = 0.5*(wingpan(j)%pnt(3)%p%set(i)%infr(1:3)+wingpan(j)%pnt(4)%p%set(i)%infr(1:3))
!      end if
!    end do


!get rid of this to increase performance, maybe cycle pointers?
  !copy wake
    do j=1,min(i-1,nwpanels+1),1
      do k=1,nwkpnt,1
        d1i1 = (j-1)*nwkpnt + k
        wakepnt(d1i1+nwkpnt)%set(i)%infr(1:3) = wakepnt(d1i1)%set(i-1)%infr(1:3)
!        wakepnt(d1i1+nwkpnt)%set(i)%infr(4) = wakepnt(d1i1)%set(i-1)%infr(4)
        wakepnt(d1i1+nwkpnt)%set(i)%act = .true.
      end do
      do k=1,nwakepan,1
        d1i1 = (j-1)*nwakepan + k
        wakepan(d1i1+nwakepan)%circ(i) = wakepan(d1i1)%circ(i-1)
        wakepan(d1i1+nwakepan)%dcirc(:,i) = wakepan(d1i1)%dcirc(:,i-1)
      end do
    end do

    !set unknown circ to 1.
    do k=1,nwakepan,1
      wakepan(k)%circ(i) = 0.
      wakepan(k)%dcirc(:,i) = 0.
    end do


    if (i>nwpanels+2) then

      !save old vortex to new vortex
      do j=1,((i-1)-(nwpanels+2))*nwakepan,1
        d1i1 = j + nwakepan
        vtxpnt(fvtx+d1i1)%set(i)%infr(1:3) = vtxpnt(fvtx+j)%set(i-1)%infr(1:3)
        vtxpnt(fvtx+d1i1)%set(i)%minfr = vtxpnt(fvtx+j)%set(i-1)%minfr
        vtxpnt(fvtx+d1i1)%set(i)%act = vtxpnt(fvtx+j)%set(i-1)%act
      end do

      !collapse vortex particle
      do j=1,nwakepan,1
        d1i1 = (nwpanels+1)*nwakepan + j
        call vtxcollapse5(nwpanels,nwakepan,d1i1,i-1,wakepan(d1i1),i,vtxpnt(fvtx+j))
      end do


    end if


    !set dcirc of panel after qtsp to circ
    if (i>1) then
      do k=1,nwakepan,1
        d1i1 = nwakepan+k
        do l=1,4,1
          if (wakepan(d1i1)%neigh(l) <= nwakepan .and. wakepan(d1i1)%neigh(l) > 0) then
            wakepan(d1i1)%dcirc(l,i) = wakepan(d1i1)%circ(i)
          end if
        end do
      end do
      if (i>nwpanels+2) then
        do k=1,nwakepan,1
          d1i1 = nwakepan*(nwpanels+1)+k
          do l=1,4,1
            if (wakepan(d1i1)%neigh(l) == -2) then
              wakepan(d1i1)%dcirc(l,i) = wakepan(d1i1)%circ(i)-wakepan(d1i1)%circ(i-1)
!do I need to check for -6, don't think so...
!            else if (wakepan(d1i1)%neigh(l) == -6) then
!              wakepan(d1i1)%dcirc(l,i) = 0. !wakepan(d1i1)%circ(i)
            end if
          end do
        end do
      else if (i>1) then
!        do k=1,nwakepan,1
!          d1i1 = nwakepan*(i-1)+k
!          do l=1,4,1
!            if (wakepan(d1i1)%neigh(l) == -2) then
!              print*, 'Error 12321'
!!do I need to check for -6, don't think so...
!            else if (wakepan(d1i1)%neigh(l) == -6) then
!              wakepan(d1i1)%dcirc(l,i) = wakepan(d1i1)%circ(i)
!            end if
!          end do
!        end do
      end if
    end if

    !shed wake
    do j=1,nwkpnt,1
      call traj(wingpnt(wkpntsrc(j))%trj, i*deltat, wingpnt(wkpntsrc(j))%acfr, tvec)
!      call ACtoINrot(cross(traj(wingpnt(wkpntsrc(j))%trj)%step(i)%tdot, &
!        & wingpnt(wkpntsrc(j))%acfr-traj(wingpnt(wkpntsrc(j))%trj)%origin),acin(wingpnt(wkpntsrc(j))%trj,:,:),d3r)
!      wakepnt(j)%set(i)%infr(1:3) = wingpnt(wkpntsrc(j))%set(i)%infr(1:3) - &
!        & qtsf*deltat*(traj(wingpnt(wkpntsrc(j))%trj)%step(i)%rdot+d3r)
      wakepnt(j)%set(i)%infr(1:3) = wingpnt(wkpntsrc(j))%set(i)%infr(1:3) - qtsf*deltat*tvec(4:6)
!      wakepnt(j)%set(i)%infr(4) = dble(i)
      wakepnt(j)%trj = wingpnt(wkpntsrc(j))%trj
    end do



  !*********************************************************************************************************************************
  !*********************************************************************************************************************************
    !need panfstream to scale, use this so helicopter works
    do j=1,npan,1
!      call traj(wingpan(j)%trj, i*deltat, wingpan(j)%cinfr(:,1), tvec)
      call traj(wingpan(j)%trj, i*deltat, wingpan(j)%cacfr, tvec)
      panfstream(:,j,i) = -1. * tvec(4:6)
!      call ACtoINrot(cross(traj(wingpan(j)%trj)%step(i)%tdot,wingpan(j)%cacfr-traj(wingpan(j)%trj)%origin), &
!        & acin(wingpan(j)%trj,:,:),d3r)
!      panfstream(:,j,i) = (-1.)*traj(wingpan(j)%trj)%step(i)%rdot - d3r
    end do
  !*********************************************************************************************************************************
  !*********************************************************************************************************************************

    !do for all wingpanel coll wingpnt
    do j=1,npan,1

      do k=1,npan,1
        wingpan(k)%prcssd = .false.
      end do

      d3r2 = wingpan(j)%cinfr(:,i)


      do k=1,np,1
!        call interppnt(d3r2,tsteps,i,1,a,deltat,wingpnt(k),wingpnt(k)%ti(1),wingpnt(k)%ti(2), &
!           & wingpnt(k)%tr(1),wingpnt(k)%tr(2),wingpnt(k)%ffl)
        call traj_interp(d3r2, i*deltat, wingpnt(k)%trj, wingpnt(k)%acfr, wingpnt(k)%te)
!        call traj_interp(d3r2, i*deltat, wingpnt(k)%trj, wingpnt(k)%set(1)%infr(1:3), wingpnt(k)%te)
        call traj(wingpnt(k)%trj, wingpnt(k)%te, wingpnt(k)%acfr, wingpnt(k)%tvec)
!        call traj(wingpnt(k)%trj, wingpnt(k)%te, wingpnt(k)%set(1)%infr(1:3), wingpnt(k)%tvec)

        if (wingpnt(k)%te/deltat >= 1.) then
          wingpnt(k)%ffl = .true.
        else
          wingpnt(k)%ffl = .false.
        end if
          
      end do





      !panel influences, without own panel
      do k=1,npan,1
        if (wingpan(k)%ptype == 'Q') then
          if (k /= j) then
            do q=1,4,1
              if (wingpan(k)%prcssd(q) .eqv. .false.) then
                wingpan(k)%prcssd(q) = .true.
                if (wingpan(k)%neigh(q) /= 0) then
                  if (wingpan(k)%neigh(q) > 0) then
                    wingpan(wingpan(k)%neigh(q))%prcssd(cycperm(4,2,q)) = .true.
                  end if
                  call wingwingcalcseg_3(d3r2,wingpan,cmat,uwind,duwind,i,j,k,q)
                end if
              end if
            end do
          end if
        else
          print*,'T panel!'
        end if
      end do

      !-----------------------------------------------------------------------------------------------------------------------------
      !Wake

      do k=1,(nwpanels+2)*nwakepan,1
        wakepan(k)%prcssd = .false.
      end do

      do k=1,min(i,nwpanels+2),1
        do l=1,nwakepan,1
          d1i1 = (k-1)*nwakepan + l
          do q=1,4,1
            call interppnt(d3r2,tsteps,i,k,a,deltat,wakepan(d1i1)%pnt(q)%p,wakepan(d1i1)%pnt(q)%p%ti(1), &
              & wakepan(d1i1)%pnt(q)%p%ti(2),wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2),wakepan(d1i1)%pnt(q)%p%ffl)
          end do
        end do
      end do


!      do k=1,min(i,nwpanels+2),1
!        do l=1,nwkpnt,1
!          d1i1 = (k-1)*nwkpnt + l
!          call interppnt(d3r2,tsteps,i,k,a,deltat,wakepnt(d1i1),wakepnt(d1i1)%ti(1), &
!            & wakepnt(d1i1)%ti(2),wakepnt(d1i1)%tr(1),wakepnt(d1i1)%tr(2),wakepnt(d1i1)%ffl)
!        end do
!      end do

!      do k=1,nwakepan,1
!        pntpntr1 => wingpan(wakesrcpan(k,1))%pnt(wakesrcpan(k,2))%p
!        pntpntr2 => wingpan(wakesrcpan(k,1))%pnt(cycperm(4,1,wakesrcpan(k,2)))%p
!        call interppnt(d3r2,tsteps,i,1,a,deltat,pntpntr1,pntpntr1%ti(1),pntpntr1%ti(2),pntpntr1%tr(1),pntpntr1%tr(2),pntpntr1%ffl)
!        call interppnt(d3r2,tsteps,i,1,a,deltat,pntpntr2,pntpntr2%ti(1),pntpntr2%ti(2),pntpntr2%tr(1),pntpntr2%tr(2),pntpntr2%ffl)
!      end do


  if ((i == 18) .and. (j == 105)) then
  open (unit=51,file='output/wake_trnsfrm.ply',status='replace',action='write',iostat=ioout)
  write(51,'(A3)') 'ply'
  write(51,'(A16)') 'format ascii 1.0'
  write(51,'(A27)') 'comment Vortex Lattice Mesh'
  write(51,'(A14,1X,I5)') 'element vertex', min(i,nwpanels+2)*nwakepan*4
  write(51,'(A16)') 'property float x'
  write(51,'(A16)') 'property float y'
  write(51,'(A16)') 'property float z'
  write(51,'(A12,1X,I4)') 'element face', min(i,nwpanels+2)*nwakepan
  write(51,'(A38)') 'property list uchar float vertex_index'
  write(51,'(A10)') 'end_header'
  do l=1,min(i,nwpanels+2),1
    do k=1,nwakepan,1
      d1i1 = (l-1)*nwakepan + k
        if ((wakepan(d1i1)%pnt(1)%p%ffl .eqv. .false.) .and. (wakepan(d1i1)%pnt(2)%p%ffl .eqv. .false.) .and. &
          & (wakepan(d1i1)%pnt(3)%p%ffl .eqv. .false.) .and. (wakepan(d1i1)%pnt(4)%p%ffl .eqv. .false.)) then
          do q=1,4,1
!        
!print*, wakepan(d1i1)%pnt(q)%p%ti(1),wakepan(d1i1)%pnt(q)%p%ti(2),wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2)
!print*, wakepan(d1i1)%pnt(q)%p%ffl
          d3r1 = intp(wakepan(d1i1)%pnt(q)%p%set(wakepan(d1i1)%pnt(q)%p%ti(1))%infr(1:3), &
           & wakepan(d1i1)%pnt(q)%p%set(wakepan(d1i1)%pnt(q)%p%ti(2))%infr(1:3), &
           & wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))

!        dummy1r1 = intq(dble(wingpan(k)%pnt(l)%p%ti(1)),dble(wingpan(k)%pnt(l)%p%ti(2)), &
!          & wingpan(k)%pnt(l)%p%tr(1),wingpan(k)%pnt(l)%p%tr(2))
          write (51,'(F10.5,1X,F10.5,1X,F10.5)') d3r1(1),d3r1(2),d3r1(3) !dummy1r1/i
        end do
      end if
    end do
  end do
  do l=1,min(i,nwpanels+2),1
    do k=1,nwakepan,1
      d1i1 = (l-1)*nwakepan + k
       write (51,'(I1,1X,I4,1X,I4,1X,I4,1X,I4)') 4,(d1i1-1)*4,(d1i1-1)*4+1,(d1i1-1)*4+2,(d1i1-1)*4+3
    end do
  end do
  close(51)
  end if



      do l=1,min(i,nwpanels+2),1
        do k=1,nwakepan,1
          d1i1 = (l-1)*nwakepan + k
          do q=1,4,1
            if (wakepan(d1i1)%prcssd(q) .eqv. .false.) then
              wakepan(d1i1)%prcssd(q) = .true.
              if (wakepan(d1i1)%neigh(q) > 0) then
                wakepan(wakepan(d1i1)%neigh(q))%prcssd(cycperm(4,2,q)) = .true.
              end if
              if ((wakepan(d1i1)%neigh(q) /= -5)) then !.and. (wakepan(d1i1)%neigh(q) /= -2) !the -2 leg does not interpolate on the last panel only tstep but that's ok, as don't want starting vortex anyway.


                call wakewingcalcseg_2(wingpan(j),d3r2,cmat, &
                   & wakesrcpan,wakepan,npan,nwpanels,nwakepan,d1i1,deltat,i,j,k,l,q)
              end if
            end if
          end do
        end do
      end do

    end do

    !-----------------------------------------------------------------------------------------------------------------------------
    !Vortex particles

    !do all vtx particles
    if (i>nwpanels+2) then

    do j=1,npan,1

!d3r2 = f1*(wingpan(j)%pnt(3)%p%set(i)%infr(1:3)+wingpan(j)%pnt(4)%p%set(i)%infr(1:3))/2. + &
!     & f2*(wingpan(j)%pnt(1)%p%set(i)%infr(1:3)+wingpan(j)%pnt(2)%p%set(i)%infr(1:3))/2.

      d3r2 = wingpan(j)%cinfr(:,i)

      do k=1,i-(nwpanels+2),1
        dummy1i2 = (nwpanels+2)+k
        do l=1,nwakepan,1
          d1i1 = fvtx + (k-1)*nwakepan + l
          !set stop to earliest possible tstep, anything outside is set to farfield.
          call interppnt(d3r2,tsteps,i,dummy1i2,a,deltat,vtxpnt(d1i1), &
            & vtxpnt(d1i1)%ti(1),vtxpnt(d1i1)%ti(2),vtxpnt(d1i1)%tr(1),vtxpnt(d1i1)%tr(2),vtxpnt(d1i1)%ffl)
        end do
      end do

      do k=1,i-(nwpanels+2),1
        do l=1,nwakepan,1
          d1i2 = fvtx + (k-1)*nwakepan + l

          call vtxwingcalcpart_wake_2(d3r2,wingpan(j),vtxpnt,dvtx,d1i2,deltat,i)

        end do
      end do

    end do
    end if

  !*********************************************************************************************************************************
  !*********************************************************************************************************************************
  !*********************************************************************************************************************************

    do j=1,npan,1
      wakeind(:,j,i) = wingpan(j)%vind
    end do

    print*, 'wakeind(max)', maxval(wakeind), maxloc(wakeind)
    print*, 'wakeind(min)', minval(wakeind), minloc(wakeind)
    write(39,*), 'wakeind(max)', maxval(wakeind), maxloc(wakeind)
    write(39,*), 'wakeind(min)', minval(wakeind), minloc(wakeind)

    !own panel influence
    do j=1,npan,1
!    d3r2 = f1*(wingpan(j)%pnt(3)%p%set(i)%infr(1:3)+wingpan(j)%pnt(4)%p%set(i)%infr(1:3))/2. + &
!      & f2*(wingpan(j)%pnt(1)%p%set(i)%infr(1:3)+wingpan(j)%pnt(2)%p%set(i)%infr(1:3))/2.
      d3r2 = wingpan(j)%cinfr(:,i)
      d3r1 = 0.
      do k=1,4,1
        if (wingpan(j)%neigh(k) /= 0) then
!          d3r1 = d3r1 + velindhypseg_new(i*deltat, d3r2, wingpan(j)%trj, wingpan(j)%pnt(k)%p%set(1)%infr(1:3), &
!            & wingpan(j)%pnt(cycperm(4,1,k))%p%set(1)%infr(1:3), dble(1.), dble(1.), rzero,rzero,zerovec,zerovec,zerovec,zerovec)
          d3r1 = d3r1 + velindhypseg_new(i*deltat, d3r2, wingpan(j)%trj, wingpan(j)%pnt(k)%p%acfr(1:3), &
            & wingpan(j)%pnt(cycperm(4,1,k))%p%acfr(1:3), dble(1.), dble(1.), rzero,rzero,zerovec,zerovec,zerovec,zerovec)
!          if (i>1) then
!            d3r1 = d3r1 + velindhypseg3D(d3r2,wingpan(j)%pnt(k)%p%set(i)%infr(1:3), &
!              & wingpan(j)%pnt(cycperm(4,1,k))%p%set(i)%infr(1:3),dble(1.),dble(1.),rzero, &
!              & traj(wingpan(j)%trj)%step(i)%rdot,panfstream(:,j,i),uwind(:,j,i-1),duwind(:,j,i-1))
!          else
!            d3r1 = d3r1 + velindhypseg3D(d3r2,wingpan(j)%pnt(k)%p%set(i)%infr(1:3), &
!              & wingpan(j)%pnt(cycperm(4,1,k))%p%set(i)%infr(1:3),dble(1.),dble(1.),rzero, &
!              & traj(wingpan(j)%trj)%step(i)%rdot,panfstream(:,j,i),zerovec,zerovec)
!          end if
        end if
      end do
      cmat(:,j,j) = cmat(:,j,j) + d3r1
    end do

    print*, 'cmat(max)', maxval(cmat), maxloc(cmat)
    print*, 'cmat(min)', minval(cmat), minloc(cmat)
    write(39,*), 'cmat(max)', maxval(cmat), maxloc(cmat)
    write(39,*), 'cmat(min)', minval(cmat), minloc(cmat)

    !populate coefficient matrix and RHS vector
    do j=1,npan,1
      d3r2 = cross(wingpan(j)%pnt(3)%p%set(i)%infr(1:3)-wingpan(j)%pnt(1)%p%set(i)%infr(1:3), &
        & wingpan(j)%pnt(2)%p%set(i)%infr(1:3)-wingpan(j)%pnt(4)%p%set(i)%infr(1:3))
!      call ACtoINrot(wingpan(j)%norm,acin(wingpan(j)%trj,:,:),d3r2)
      d3r2 = d3r2/sqrt(dot_product(d3r2,d3r2))
      do k=1,npan,1
        COEFFS(j,k) = dot_product(cmat(:,k,j),d3r2)
      end do
!      call ACtoINrot(cross(traj(wingpan(j)%trj)%step(i)%tdot,wingpan(j)%cacfr-traj(wingpan(j)%trj)%origin), &
!        & acin(wingpan(j)%trj,:,:),d3r)
!      RHS(j) =  (-1.) * dot_product(((-1.)*traj(wingpan(j)%trj)%step(i)%rdot + wakeind(:,j,i) - d3r),d3r2)
      call traj(wingpan(j)%trj, i*deltat, wingpan(j)%cacfr, tvec)
!      call traj(wingpan(j)%trj, i*deltat, wingpan(j)%cinfr(:,i), tvec)
      RHS(j) =  (-1.) * dot_product(((-1.)*tvec(4:6) + wakeind(:,j,i)),d3r2)
    end do

    print*, 'RHS(max)', maxval(RHS), maxloc(RHS)
    print*, 'RHS(min)', minval(RHS), minloc(RHS)
    write(39,*), 'RHS(max)', maxval(RHS), maxloc(RHS)
    write(39,*), 'RHS(min)', minval(RHS), minloc(RHS)

    call DGESV(npan,1,COEFFS,npan,IPIV,RHS,npan,INFO)
    wingcirc(:,i) = RHS(:)

!    if (i>1) then
!      call bicgstabup(npan,bstabprec,COEFFS,RHS,wingcirc(:,i-1),wingcirc(:,i))
!    else if (i==1) then
!      wingcirc(:,i) = 0.
!      call bicgstabup(npan,bstabprec,COEFFS,RHS,wingcirc(:,i),wingcirc(:,i))
!    end if

    print*, 'wingcirc(max)',maxval(wingcirc(:,i)),maxloc(wingcirc(:,i))
    print*, 'wingcirc(min)',minval(wingcirc(:,i)),minloc(wingcirc(:,i))
    write(39,*), 'wingcirc(max)',maxval(wingcirc(:,i)),maxloc(wingcirc(:,i))
    write(39,*), 'wingcirc(min)',minval(wingcirc(:,i)),minloc(wingcirc(:,i))


!    write(istring,'(i4.4)') i
!    open (unit=35,file='output/wing/wing.csv.'//trim(istring),status='replace',action='write',iostat=ioout)
!    write (35,'(A33)') 'x_coord, y_coord, z_coord, scalar'
!    do k=1,npan,1
!      write (35,'(F10.5,1A,1X,F10.5,1A,1X,F10.5,1A,1X,F10.5)') wingpan(k)%cacfr(1),',',wingpan(k)%cacfr(2),',', &
!        & wingpan(k)%cacfr(3),',', COEFFS(315,k) !wingpan(k)%circ(i)
!    end do
!    close(35)


    !save wingcirc to wingpan data structure
    do j=1,npan,1
      wingpan(j)%circ(i) = wingcirc(j,i)
    end do

    !save wakecirc of shed wake panels
    do j=1,nwakepan,1
      wakepan(j)%circ(i) = wingpan(wakesrcpan(j,1))%circ(i)
    end do

    !caculate dcirc for first wakepan and second
  if (i>2) then
    do j=1,nwakepan*2,1
      do k=1,4,1
        if (wakepan(j)%neigh(k) > 0) then
          wakepan(j)%dcirc(k,i) = wakepan(j)%circ(i) - wakepan(wakepan(j)%neigh(k))%circ(i)
        else if (wakepan(j)%neigh(k) == 0) then
          wakepan(j)%dcirc(k,i) = wakepan(j)%circ(i)
        else if (wakepan(j)%neigh(k) == -5) then
          !Kutta condition
          wakepan(j)%dcirc(k,i) = 0.
        else if (wakepan(j)%neigh(k) == -6) then
          wakepan(j)%dcirc(k,i) = wakepan(j)%circ(i)
print*, 'here245'
        else if (wakepan(j)%neigh(k) == -2) then
print*, 'here6354'
          wakepan(j)%dcirc(k,i) = wakepan(j)%circ(i)-wakepan(j)%circ(i-1)
        else
          print*, 'error in main file: 78324'
        end if
      end do
    end do
  else if (i>1) then
     do j=1,nwakepan*2,1
      do k=1,4,1
        if (wakepan(j)%neigh(k) > 0) then
          if (wakepan(j)%neigh(k) > 2.*nwakepan) then
            wakepan(j)%dcirc(k,i) = wakepan(j)%circ(i)
          else
            wakepan(j)%dcirc(k,i) = wakepan(j)%circ(i) - wakepan(wakepan(j)%neigh(k))%circ(i)
          end if
        else if (wakepan(j)%neigh(k) == 0) then
          wakepan(j)%dcirc(k,i) = wakepan(j)%circ(i)
        else if (wakepan(j)%neigh(k) == -5) then
          !Kutta condition
          wakepan(j)%dcirc(k,i) = 0.
        else
          print*, 'error in main file: 78324'
        end if
      end do
    end do
  else
     do j=1,nwakepan,1
      do k=1,4,1
        if (wakepan(j)%neigh(k) > 0) then
          if (wakepan(j)%neigh(k) > nwakepan) then
            wakepan(j)%dcirc(k,i) = wakepan(j)%circ(i)
          else
            wakepan(j)%dcirc(k,i) = wakepan(j)%circ(i) - wakepan(wakepan(j)%neigh(k))%circ(i)
          end if
        else if (wakepan(j)%neigh(k) == 0) then
          wakepan(j)%dcirc(k,i) = wakepan(j)%circ(i)
        else if (wakepan(j)%neigh(k) == -5) then
          !Kutta condition
          wakepan(j)%dcirc(k,i) = 0.
        else
          print*, 'error in main file: 78324'
        end if
      end do
    end do
  end if


    if (i .ge. nwpanels+2) then
      do j=1,nwpanels+2,1
        do k=1,nwakepan,1
          d1i1 = (j-1)*nwakepan+k
          do l=1,4,1
            if (wakepan(d1i1)%neigh(l) == -6) then
              wakepan(d1i1)%dcirc(l,i) = wakepan(d1i1)%circ(i)
            else if (wakepan(d1i1)%neigh(l) == -2) then
              wakepan(d1i1)%dcirc(l,i) = wakepan(d1i1)%circ(i)-wakepan(d1i1)%circ(i-1)
            else if (wakepan(d1i1)%neigh(l) == -5) then
              wakepan(d1i1)%dcirc(l,i) = 0.
            else if (wakepan(d1i1)%neigh(l) == 0) then
              wakepan(d1i1)%dcirc(l,i) = wakepan(d1i1)%circ(i)
            else if (wakepan(d1i1)%neigh(l) > 0) then
              wakepan(d1i1)%dcirc(l,i) = wakepan(d1i1)%circ(i)-wakepan(wakepan(d1i1)%neigh(l))%circ(i)
            end if
          end do
        end do
      end do
    end if


    !calculate induced velocities for panels
    do j=1,npan,1
      do k=1,npan,1
        wingind(:,j,i) = wingind(:,j,i) + cmat(:,k,j) * wingcirc(k,i)
      end do
      !calculate deltawingcirc
      do k=1,4,1
        if (wingpan(j)%neigh(k) > 0) then
          wingpan(j)%dcirc(k,i) = wingcirc(j,i) - wingcirc(wingpan(j)%neigh(k),i)
        else if (wingpan(j)%neigh(k) == 0) then
          wingpan(j)%dcirc(k,i) = 0.
        else if (wingpan(j)%neigh(k) == -1) then
          wingpan(j)%dcirc(k,i) = wingcirc(j,i)
        else
          print*, 'Error 781234'
        end if
      end do
    end do

  !*********************************************************************************************************************************

    !Force calculation
    !add rotational velocity to induced velocity
    uwind(:,:,i) = wingind(:,:,i) + wakeind(:,:,i)

    if (i>1) then
      duwind(:,:,i) = (uwind(:,:,i)-uwind(:,:,i-1))/deltat
    else
      duwind(:,:,i) = uwind(:,:,i)/deltat
    end if

    do j=1,npan,1
      panvind(:,j) = uwind(:,j,i) + panfstream(:,j,i)
    end do


    !call force calculation routine

    call forcecalcM(i,tsteps,deltat,density,np,npan,wingpnt,a,wingpan,wingcirc, &
      & panvind,panfstream,panforce(:,:,i),pancp(:,i))

!    call forcecalc2(i,i,tsteps,deltat,acin(1,:,:),density,np,npan,wingpnt,wingpan,wingcirc, &
!      & panvind,panfstream,panforce(:,:,i),pancp(:,i))
    lift(i) = 0.
    drag(i) = 0.
    do j=1,npan,1
      if (wingpan(j)%trj == 1) then
        lift(i) = lift(i) + panforce(3,j,i)
        drag(i) = drag(i) + panforce(1,j,i)
      end if
    end do
    print*, 'Cl:', lift(i)/(0.5*density*uinf**2*totalA)
    print*, 'Cd:', drag(i)/(0.5*density*uinf**2*totalA)

  !*********************************************************************************************************************************
  !*********************************************************************************************************************************
  !*********************************************************************************************************************************

!do loop below in here


    if (trigger2) then
    !calc influence for wakepoints
    do j=1,min(i,(nwpanels+2))*nwkpnt,1

      !Wingpanels
      do k=1,npan,1
        wingpan(k)%prcssd = .false.
      end do

      do k=1,np,1
        call interppnt(wakepnt(j)%set(i)%infr(1:3),tsteps,i,1,a,deltat,wingpnt(k),wingpnt(k)%ti(1),wingpnt(k)%ti(2), &
          & wingpnt(k)%tr(1),wingpnt(k)%tr(2),wingpnt(k)%ffl)
      end do

      do k=1,npan,1
        if (wingpan(k)%ptype == 'Q') then
          do q=1,4,1
            if ((wingpan(k)%pnt(q)%p%ffl .eqv. .false.) .and. (wingpan(k)%pnt(cycperm(4,1,q))%p%ffl .eqv. .false.)) then
              if (wingpan(k)%prcssd(q) .eqv. .false.) then
                wingpan(k)%prcssd(q) = .true.
                if (wingpan(k)%neigh(q) /= 0) then

!              wakepnt(j)%vind = wakepnt(j)%vind + & 
!                & velindhypseg_new( i*deltat, wakepnt(j)%set(i)%infr(1:3), wingpan(k)%trj, wingpan(k)%pnt(q)%p%acfr(1:3), &
!                & wingpan(k)%pnt(cycperm(4,1,q))%p%acfr(1:3), dble(1.), dble(1.), rzero,rzero,zerovec,zerovec,zerovec,zerovec)

              call wingwakecalcseg_3( wakepnt(j)%set(i)%infr(1:3), wingpan, zerovec, zerovec, i, k, q, wakepnt(j)%vind)

!              wakepnt(j)%vind = wakepnt(j)%vind + &
!                & velindhypseg3D_wake(wakepnt(j)%set(i)%infr(1:3),d3r1,d3r2,dummy1r1,dummy1r2, &
!                & traj(wingpan(k)%trj)%step(i)%rdot,zerovec) 
                !velindlinseg(wakepnt(j)%set(i)%infr(1:3),d3r1,d3r2,dummy1r1,dummy1r2)

              if (wingpan(k)%neigh(q) > 0) then
                wingpan(wingpan(k)%neigh(q))%prcssd(cycperm(4,2,q)) = .true.
              end if

                end if
              end if
!            else
!              if ((wingpan(k)%prcssd(q) .eqv. .false.) .and. (wingpan(k)%neigh(q) /= 0)) then
!                wingpan(k)%prcssd(q) = .true.

!                wakepnt(j)%vind = wakepnt(j)%vind + velindhypseg3D_wake(wakepnt(j)%set(i)%infr(1:3), &
!                  & wingpan(k)%pnt(q)%p%set(1)%infr(1:3),wingpan(k)%pnt(cycperm(4,1,q))%p%set(1)%infr(1:3), &
!                  & wingpan(k)%dcirc(q,1),wingpan(k)%dcirc(q,1),traj(wingpan(k)%trj)%step(i)%rdot)

!                  !velindlinseg(wakepnt(j)%set(i)%infr(1:3), &
!                  !& wingpan(k)%pnt(q)%p%set(1)%infr(1:3),wingpan(k)%pnt(cycperm(4,1,q))%p%set(1)%infr(1:3), &
!                  !& wingpan(k)%dcirc(q,1),wingpan(k)%dcirc(q,1))

!                if (wingpan(k)%neigh(q) > 0) then
!                  wingpan(wingpan(k)%neigh(q))%prcssd(cycperm(4,2,q)) = .true.
!                end if
!              end if
            end if
          end do
        end if
      end do

      !Wakepanels
      do k=1,(nwpanels+2)*nwakepan,1
        wakepan(k)%prcssd = .false.
      end do
      do k=1,min(i,nwpanels+2),1
        do l=1,nwkpnt,1
          d1i1 = (k-1)*nwkpnt + l
          !set stop to earliest possible tstep, anything outside is set to farfield.
          call interppnt(wakepnt(j)%set(i)%infr(1:3),tsteps,i,k,a,deltat,wakepnt(d1i1),wakepnt(d1i1)%ti(1), &
            & wakepnt(d1i1)%ti(2),wakepnt(d1i1)%tr(1),wakepnt(d1i1)%tr(2),wakepnt(d1i1)%ffl)
        end do
      end do

      do l=1,nwkpnt,1
        call interppnt(wakepnt(j)%set(i)%infr(1:3),tsteps,i,1,a,deltat,wingpnt(wkpntsrc(l)),wingpnt(wkpntsrc(l))%ti(1), &
          & wingpnt(wkpntsrc(l))%ti(2),wingpnt(wkpntsrc(l))%tr(1),wingpnt(wkpntsrc(l))%tr(2),wingpnt(wkpntsrc(l))%ffl)
      end do


      do l=1,min(i,nwpanels+2),1
        do k=1,nwakepan,1
          d1i1 = (l-1)*nwakepan + k
          do q=1,4,1
            if ((wakepan(d1i1)%pnt(q)%p%ffl .eqv. .false.) .and. &
              & (wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ffl .eqv. .false.)) then

            if (wakepan(d1i1)%prcssd(q) .eqv. .false.) then
              wakepan(d1i1)%prcssd(q) = .true.
            if (wakepan(d1i1)%neigh(q) /= -5) then

              d3r1 = intp(wakepan(d1i1)%pnt(q)%p%set(wakepan(d1i1)%pnt(q)%p%ti(1))%infr(1:3), &
                & wakepan(d1i1)%pnt(q)%p%set(wakepan(d1i1)%pnt(q)%p%ti(2))%infr(1:3), &
                & wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))

             d3r2 = intp(wakepan(d1i1)%pnt(cycperm(4,1,q))%p%set(wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1))%infr(1:3), &
                & wakepan(d1i1)%pnt(cycperm(4,1,q))%p%set(wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2))%infr(1:3), &
                & wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))

              dummy1r1 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(1)), &
                 & wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2)), &
                 & wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))

              dummy1r2 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)), &
                & wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2)), &
                & wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))

              wakepnt(j)%vind = wakepnt(j)%vind + velindlinseg(wakepnt(j)%set(i)%infr(1:3),d3r1,d3r2,dummy1r1,dummy1r2)
    
              if (wakepan(d1i1)%neigh(q) > nwakepan) then
                wakepan(wakepan(d1i1)%neigh(q))%prcssd(cycperm(4,2,q)) = .true.
              end if

            end if
            end if

            end if
          end do
        end do
      end do

      !Vortices

      do k=1,i-(nwpanels+2),1
        dummy1i2 = (nwpanels+2)+k
        do l=1,nwakepan,1
          d1i1 = fvtx + (k-1)*nwakepan + l
          !set stop to earliest possible tstep, anything outside is set to farfield.
          call interppnt(wakepnt(j)%set(i)%infr(1:3),tsteps,i,dummy1i2,a,deltat,vtxpnt(d1i1),vtxpnt(d1i1)%ti(1), &
            & vtxpnt(d1i1)%ti(2),vtxpnt(d1i1)%tr(1),vtxpnt(d1i1)%tr(2),vtxpnt(d1i1)%ffl)
        end do
      end do

      do k=1,i-(nwpanels+2),1
        dummy1i2 = (nwpanels+2)+k
        do l=1,nwakepan,1
          d1i1 = fvtx + (k-1)*nwakepan + l
          if (vtxpnt(d1i1)%ffl .eqv. .false.) then

            d3r1 = intp(vtxpnt(d1i1)%set(vtxpnt(d1i1)%ti(1))%infr(1:3), &
              & vtxpnt(d1i1)%set(vtxpnt(d1i1)%ti(2))%infr(1:3),vtxpnt(d1i1)%tr(1),vtxpnt(d1i1)%tr(2))

            d3r2 = intp(vtxpnt(d1i1)%set(vtxpnt(d1i1)%ti(1))%minfr, &
              & vtxpnt(d1i1)%set(vtxpnt(d1i1)%ti(2))%minfr,vtxpnt(d1i1)%tr(1),vtxpnt(d1i1)%tr(2))

            call vtxinduced(wakepnt(j)%set(i)%infr(1:3),d3r1,d3r2,d3r)

            wakepnt(j)%vind = wakepnt(j)%vind + d3r

          end if
        end do
      end do


    end do
    end if


  !*********************************************************************************************************************************
  !*********************************************************************************************************************************

    if (trigger2) then

    if (i>nwpanels+2) then
      !calculate vortex influence coefficients
    do j=fvtx+1,fvtx+(i-(nwpanels+2))*nwakepan,1

    if (vtxpnt(j)%set(i)%act) then

      !Wingpanels
      do k=1,npan,1
        wingpan(k)%prcssd = .false.
      end do
      do k=1,np,1
        call interppnt(vtxpnt(j)%set(i)%infr(1:3),tsteps,i,1,a,deltat,wingpnt(k),wingpnt(k)%ti(1),wingpnt(k)%ti(2), &
          & wingpnt(k)%tr(1),wingpnt(k)%tr(2),wingpnt(k)%ffl)
      end do

      do k=1,npan,1
        if (wingpan(k)%ptype == 'Q') then
          do q=1,4,1
            if (wingpan(k)%pnt(q)%p%ffl .eqv. .false.) then
            if (wingpan(k)%pnt(cycperm(4,1,q))%p%ffl .eqv. .false.) then
            if (wingpan(k)%prcssd(q) .eqv. .false.) then
              wingpan(k)%prcssd(q) = .true.
            if (wingpan(k)%neigh(q) /= 0) then

              call wingwakecalcseg_3( vtxpnt(j)%set(i)%infr(1:3), wingpan, zerovec, zerovec, i, k, q, vtxpnt(j)%vind)

!              d3r1 = 0.5*(d3r1+d3r2)
!              dummy1r1 = 0.5*(dummy1r1+dummy1r2)
!              d3r2 = dummy1r1*(d3r2-d3r1)
!              call shearind(vtxpnt(j)%set(i)%infr(1:3),d3r1,d3r2,dummy33r)
!              vtxpnt(j)%gradv = vtxpnt(j)%gradv + dummy33r

              if (wingpan(k)%neigh(q) > 0) then
                wingpan(wingpan(k)%neigh(q))%prcssd(cycperm(4,2,q)) = .true.
              end if

            end if
            end if
            end if
            end if
          end do
        end if
      end do

      !Wakepanels
      do k=1,(nwpanels+2)*nwakepan,1
        wakepan(k)%prcssd = .false.
      end do
      do k=1,min(i,nwpanels+2),1
        do l=1,nwkpnt,1
          d1i1 = (k-1)*nwkpnt + l
          !set stop to earliest possible tstep, anything outside is set to farfield.
          call interppnt(vtxpnt(j)%set(i)%infr(1:3),tsteps,i,k,a,deltat,wakepnt(d1i1),wakepnt(d1i1)%ti(1), &
            & wakepnt(d1i1)%ti(2),wakepnt(d1i1)%tr(1),wakepnt(d1i1)%tr(2),wakepnt(d1i1)%ffl)
        end do
      end do


      do l=1,min(i,nwpanels+2),1
        do k=1,nwakepan,1
          d1i1 = (l-1)*nwakepan + k
          do q=1,4,1
            if (wakepan(d1i1)%pnt(q)%p%ffl .eqv. .false.) then
            if (wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ffl .eqv. .false.) then
            if (wakepan(d1i1)%prcssd(q) .eqv. .false.) then
              wakepan(d1i1)%prcssd(q) = .true.
            if (wakepan(d1i1)%neigh(q) /= -5) then

              d3r1 = intp(wakepan(d1i1)%pnt(q)%p%set(wakepan(d1i1)%pnt(q)%p%ti(1))%infr(1:3), &
                & wakepan(d1i1)%pnt(q)%p%set(wakepan(d1i1)%pnt(q)%p%ti(2))%infr(1:3), &
                & wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))

              d3r2 = intp(wakepan(d1i1)%pnt(cycperm(4,1,q))%p%set(wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1))%infr(1:3), &
                & wakepan(d1i1)%pnt(cycperm(4,1,q))%p%set(wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2))%infr(1:3), &
                & wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))

              dummy1r1 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(1)), &
                 & wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2)), &
                 & wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))

              dummy1r2 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)), &
                & wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2)), &
                & wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))

              vtxpnt(j)%vind = vtxpnt(j)%vind + velindlinseg(vtxpnt(j)%set(i)%infr(1:3),d3r1,d3r2,dummy1r1,dummy1r2)


              d3r1 = 0.5*(d3r1+d3r2)
              dummy1r1 = 0.5*(dummy1r1+dummy1r2)
              d3r2 = dummy1r1*(d3r2-d3r1)
              call shearind(vtxpnt(j)%set(i)%infr(1:3),d3r1,d3r2,dummy33r)
              vtxpnt(j)%gradv = vtxpnt(j)%gradv + dummy33r

              if (wakepan(d1i1)%neigh(q) > nwakepan) then
                wakepan(wakepan(d1i1)%neigh(q))%prcssd(cycperm(4,2,q)) = .true.
              end if

            end if
            end if
            end if
            end if
          end do
        end do
      end do

      !Vortices


      do k=1,i-(nwpanels+2),1
        dummy1i2 = (nwpanels+2)+k
        do l=1,nwakepan,1
          d1i1 = fvtx + (k-1)*nwakepan + l
          !set stop to earliest possible tstep, anything outside is set to farfield.
          call interppnt(vtxpnt(j)%set(i)%infr(1:3),tsteps,i,dummy1i2,a,deltat,vtxpnt(d1i1),vtxpnt(d1i1)%ti(1), &
            & vtxpnt(d1i1)%ti(2),vtxpnt(d1i1)%tr(1),vtxpnt(d1i1)%tr(2),vtxpnt(d1i1)%ffl)
        end do
      end do

      do k=1,i-(nwpanels+2),1
        dummy1i2 = (nwpanels+2)+k
        do l=1,nwakepan,1
          d1i1 = fvtx + (k-1)*nwakepan + l
          if (vtxpnt(d1i1)%ffl .eqv. .false.) then

            d3r1 = intp(vtxpnt(d1i1)%set(vtxpnt(d1i1)%ti(1))%infr(1:3), &
              & vtxpnt(d1i1)%set(vtxpnt(d1i1)%ti(2))%infr(1:3),vtxpnt(d1i1)%tr(1),vtxpnt(d1i1)%tr(2))

            d3r2 = intp(vtxpnt(d1i1)%set(vtxpnt(d1i1)%ti(1))%minfr, &
              & vtxpnt(d1i1)%set(vtxpnt(d1i1)%ti(2))%minfr,vtxpnt(d1i1)%tr(1),vtxpnt(d1i1)%tr(2))

            call vtxinduced(vtxpnt(j)%set(i)%infr(1:3),d3r1,d3r2,d3r)
            vtxpnt(j)%vind = vtxpnt(j)%vind + d3r

            call shearind(vtxpnt(j)%set(i)%infr(1:3),d3r1,d3r2,dummy33r)
            vtxpnt(j)%gradv = vtxpnt(j)%gradv + dummy33r

          end if
        end do
      end do


    end if
    end do
    end if
    end if

!***********************************************************************************************************************************
!***********************************************************************************************************************************
!***********************************************************************************************************************************
!***********************************************************************************************************************************



    if (trigger2) then

      !advance wake panels
      do j=1,min(i,nwpanels+2),1
        do k=1,nwkpnt,1
          d1i1 = (j-1)*nwkpnt + k
          wakepnt(d1i1)%set(i)%infr(1:3) = wakepnt(d1i1)%set(i)%infr(1:3) + deltat*wakepnt(d1i1)%vind
        end do
      end do

      !check vortex penetration
      !call pen_check(fvtx+1, fvtx+(i-(nwpanels+2))*nwakepan, i, vtxpnt, wingpan)

      !advance and rotate vtxs
      do j=1,(i-(nwpanels+2))*nwakepan,1
          if (vtxpnt(fvtx+j)%set(i)%act) then
            call advvtx(deltat,vtxpnt(fvtx+j)%vind,vtxpnt(fvtx+j)%gradv,vtxpnt(fvtx+j)%set(i)%infr(1:3),vtxpnt(fvtx+j)%set(i)%minfr)
          end if
      end do

    end if

  !############# OUTPUT ##############

    write(istring,'(i4.4)') i

    open (unit=34,file=trim('output/wing/wing.ply.'//trim(istring)),status='replace',action='write',iostat=ioout)
    write(34,'(A3)') 'ply'
    write(34,'(A16)') 'format ascii 1.0'
    write(34,'(A27)') 'comment Vortex Lattice Mesh'
    write(34,'(A14,1X,I4)') 'element vertex', npan*4
    write(34,'(A16)') 'property float x'
    write(34,'(A16)') 'property float y'
    write(34,'(A16)') 'property float z'
    write(34,'(A12,1X,I4)') 'element face', npan
    write(34,'(A38)') 'property list uchar float vertex_index'
    write(34,'(A17)') 'property float cp'
    write(34,'(A10)') 'end_header'
    do k=1,npan,1
      do j=1,4,1
        write (34,'(F10.5,1X,F10.5,1X,F10.5)') wingpan(k)%pnt(j)%p%set(i)%infr(1:3)
      end do
    end do
    do k=1,npan,1
      write (34,'(I1,1X,I4,1X,I4,1X,I4,1X,I4,1X,F10.5)') 4,(k-1)*4,(k-1)*4+1,(k-1)*4+2,(k-1)*4+3, pancp(k,i)
    end do
    close(34)

    if (i > nwpanels+1) then

      open (unit=31,file=trim('output/wake/wake.ply.'//trim(istring)),status='replace',action='write',iostat=ioout)
      write(31,'(A3)') 'ply'
      write(31,'(A16)') 'format ascii 1.0'
      write(31,'(A27)') 'comment Vortex Lattice Mesh'
      write(31,'(A14,1X,I4)') 'element vertex', (nwpanels+2)*nwakepan*4
      write(31,'(A16)') 'property float x'
      write(31,'(A16)') 'property float y'
      write(31,'(A16)') 'property float z'
      write(31,'(A12,1X,I4)') 'element face', (nwpanels+2)*nwakepan
      write(31,'(A38)') 'property list uchar float vertex_index'
      write(31,'(A10)') 'end_header'
      do k=1,(nwpanels+2)*nwakepan,1
        do j=1,4,1
          write (31,'(F10.5,1X,F10.5,1X,F10.5)') wakepan(k)%pnt(j)%p%set(i)%infr(1:3)
        end do
      end do
      do k=1,(nwpanels+2)*nwakepan,1
        write (31,'(I1,1X,I4,1X,I4,1X,I4,1X,I4)') 4,(k-1)*4,(k-1)*4+1,(k-1)*4+2,(k-1)*4+3
      end do
      close(31)

    end if

    if (i > nwpanels+2) then
      open (unit=32,file='output/wake/wake.csv.'//trim(istring),status='replace',action='write',iostat=ioout)
      write (32,'(A45)') 'x_coord, y_coord, z_coord, scalar, mx, my, mz'
      do j=fvtx+1,fvtx + (i-(nwpanels+2))*nwakepan,1
          if (vtxpnt(j)%set(i)%act) then
            write (32,'(F10.5,1A,1X,F10.5,1A,1X,F10.5,1A,1X,F10.5,1A,1X,F10.5,1A,1X,F10.5,1A,1X,F10.5)')  &
                & vtxpnt(j)%set(i)%infr(1),',',vtxpnt(j)%set(i)%infr(2),',',vtxpnt(j)%set(i)%infr(3),',',&
                & sqrt(dot_product(vtxpnt(j)%set(i)%minfr,vtxpnt(j)%set(i)%minfr)),',', &
                & vtxpnt(j)%set(i)%minfr(1),',',vtxpnt(j)%set(i)%minfr(2),',',vtxpnt(j)%set(i)%minfr(3)
          end if
      end do
      close(32)
    end if

    write(istring,'(i4.4)') i

    open (unit=35,file='output/wing/wing.csv.'//trim(istring),status='replace',action='write',iostat=ioout)
    write (35,'(A33)') 'x_coord, y_coord, z_coord, scalar'
    do k=1,npan,1
!      write (35,'(F10.5,1A,1X,F10.5,1A,1X,F10.5,1A,1X,F10.5)') wingpan(k)%cinfr(1,i),',',wingpan(k)%cinfr(2,i),',', &
!        & wingpan(k)%cinfr(3,i),',', wingpan(k)%circ(i)

      d3r1 = traj_pos(wingpan(k)%trj, (i+1) * deltat, wingpan(k)%cacfr)

      write (35,'(F10.5,1A,1X,F10.5,1A,1X,F10.5,1A,1X,F10.5)') d3r1(1), ',', d3r1(2), ',', d3r1(3), ',', pancp(k,i) !wingpan(k)%circ(i)

    end do
    close(35)

 !     do k=1,npan,1
 !       if (wingpan(k)%neigh(1) == -1) then
 !         write(kstring,'(i4.4)') k
 !         open (unit=33,file='output/cp/cp.csv.'//trim(kstring)//'.'//trim(istring),status='replace',action='write',iostat=ioout)
 !         write (33,'(A33)') 'x_coord, y_coord, z_coord, scalar'
 !         write (33,'(F10.5,1A,1X,F10.5,1A,1X,F10.5,1A,1X,F10.5)') wingpan(k)%cinfr(1,i),',',wingpan(k)%cinfr(2,i),',', &
 !           & wingpan(k)%cinfr(3,i),',', (-1.)*pancp(k,i)
 !         d1i1 = k
 !         jloop: do j=1,npan,1
 !           if (wingpan(d1i1)%neigh(3) > 0) then
 !             d1i1 = wingpan(d1i1)%neigh(3)
 !             write (33,'(F10.5,1A,1X,F10.5,1A,1X,F10.5,1A,1X,F10.5)') wingpan(d1i1)%cinfr(1,i),',', &
 !               & wingpan(d1i1)%cinfr(2,i),',',wingpan(d1i1)%cinfr(3,i),',', (-1.)*pancp(d1i1,i)
 !           else
 !             exit jloop
 !           end if
 !         end do jloop
 !         close(33)
 !       end if
 !     end do

      open (unit=33,file='output/cp/cp.csv.'//trim(istring),status='replace',action='write',iostat=ioout)
      do k=1,npan,1
        if (wingpan(k)%neigh(1) == -1) then
          write (33,*) k
          write (33,'(A33)') 'x_coord, y_coord, z_coord, scalar'
          write (33,'(F10.5,1A,1X,F10.5,1A,1X,F10.5,1A,1X,F10.5)') wingpan(k)%cinfr(1,i),',',wingpan(k)%cinfr(2,i),',', &
            & wingpan(k)%cinfr(3,i),',', (-1.)*pancp(k,i)
          d1i1 = k
          jloop: do j=1,npan,1
            if (wingpan(d1i1)%neigh(3) > 0) then
              d1i1 = wingpan(d1i1)%neigh(3)
              write (33,'(F10.5,1A,1X,F10.5,1A,1X,F10.5,1A,1X,F10.5)') wingpan(d1i1)%cinfr(1,i),',', &
                & wingpan(d1i1)%cinfr(2,i),',', wingpan(d1i1)%cinfr(3,i),',', pancp(d1i1,i)
            else
              exit jloop
            end if
          end do jloop
          close(33)
        end if
      end do

  !###################################

    write(istring,'(i4.4)') i
    open (unit=34,file='output/cl_sec/clsec'//trim(istring)//'.csv',status='replace',action='write',iostat=ioout)
    do j=1,npan_wing,1
      if (wingpan(j)%neigh(1) == -1) then
        d1i1 = j
        area = wingpan(d1i1)%area
        lift_seg = panforce(3,d1i1,i)
        write(34,'(F12.8)',advance="no") wingpan(d1i1)%cacfr(2)
        do
          if (wingpan(d1i1)%neigh(3) > 0) then
            d1i1 = wingpan(d1i1)%neigh(3)
            area = area + wingpan(d1i1)%area
            lift_seg = lift_seg + panforce(3,d1i1,i)
          else
            exit
          end if
        end do
        write(34,*) ',', lift_seg/(area*0.5*density*dot_product(panfstream(:,j,i),panfstream(:,j,i)))
      end if
    end do
    close(34)

  end do

  call dataout3(panfstream,wingpan,panforce,pancp)

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

  call cpu_time(timeoff)
  print*, 'Runtime:', (timeoff-timeon)/60., 'minutes'
  write(39,*) 'Runtime:', (timeoff-timeon)/60., 'minutes'
  close(39)
end program vlatunsteady
