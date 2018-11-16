module basictypes
  implicit none

  type rmset
    real*8,dimension(4) :: infr
    real*8,dimension(3) :: minfr
    logical :: act
  end type rmset

  type point
    type(rmset),dimension(:),allocatable :: set
    real*8,dimension(3) :: acfr,vind
    real*8,dimension(3,3) :: gradv
    integer :: trj,tag
    integer,dimension(2) :: ti
    real*8,dimension(2) :: tr
    real(kind=8) :: te
    real(kind=8),dimension(12) :: tvec
    logical :: ffl
    character(len=1) :: typ
    integer,dimension(2) :: pseg
    integer :: n
  end type point

  type pntptr
    type(point),pointer :: p => null()
  end type pntptr

  type panel
    type(pntptr),dimension(4) :: pnt
    integer,dimension(4) :: neigh
    logical,dimension(4) :: prcssd
    real*8,dimension(:),allocatable :: circ
    real*8,dimension(:,:),allocatable :: dcirc,cinfr
    real*8,dimension(3) :: cacfr,norm,vind
    real*8 :: area
    character(len=1) :: ptype
    integer :: trj
  end type panel

  type parttraj
    !referenced to origin in inertial frame
    real*8,dimension(3) :: r,rdot,t,tdot
  end type parttraj

  type trajectory
    real*8,dimension(3) :: origin
    type(parttraj),dimension(:),allocatable :: step
  end type trajectory

end module basictypes

!###################################################################################################################################
!###################################################################################################################################

module fourdgridtypes
  use basictypes
  implicit none

  type chld
    type(node),pointer :: child => null()
!    integer,dimension(4) :: i4
  end type chld

  type prnt
    type(node),pointer :: parent => null()
  end type prnt

  type zeronode
    type(rmset),pointer:: vpnt => null()
    type(point),pointer :: vset => null()
    real*8 :: w
    
  end type zeronode

  type node
    type(prnt) :: p
    type(chld),dimension(16) :: c
    logical :: rclcflg,usdflg
    integer :: level
    integer :: nvtx = 0
    type(zeronode),dimension(20) :: bvtx
    real*8,dimension(4) :: centroid
    real*8,dimension(3) :: mentroid
    real*8 :: w
  end type node

  type nodepoint
    type(node),pointer :: npnt => null()
  end type

  type level
    integer,dimension(3) :: nlevelem
    type (nodepoint),dimension(:,:,:,:),allocatable :: elem
  end type level

!##################################################################################################################################

  contains

  !returns true if inside buffer region
  function inoutNF(nlevel,lev,i4calc,buffer)
    implicit none
    logical :: inoutNF
    integer,intent(in) :: nlevel,lev
    integer,dimension(4),intent(in) :: i4calc
    integer,dimension(2,3,0:nlevel),intent(in) :: buffer
    integer :: i
    inoutNF = .false.
    if ((i4calc(1) .gt. buffer(1,1,lev)) .and. (i4calc(1) .le. buffer(2,1,lev))) then
      if ((i4calc(2) .gt. buffer(1,2,lev)) .and. (i4calc(2) .le. buffer(2,2,lev))) then
        if ((i4calc(3) .gt. buffer(1,3,lev)) .and. (i4calc(3) .le. buffer(2,3,lev))) then
          inoutNF = .true.
        end if
      end if
    end if

  end function inoutNF


  function boxcentroid(deltat,lev,i4in)
    use input, only: a
    implicit none
    real*8,dimension(4) :: boxcentroid
    real*8,intent(in) :: deltat
    integer,intent(in) :: lev
    integer,dimension(4),intent(in) :: i4in
    integer :: i
    do i=1,3,1
      boxcentroid(i) = (dble(i4in(i))-sign(dble(0.5),dble(i4in(i)))) * (2**lev) * (a*deltat)
    end do
    boxcentroid(4) = (dble(i4in(4))-sign(dble(0.5),dble(i4in(4)))) * (2**lev)
  end function boxcentroid

  function tcentroid(deltat,lev,i1in)
    implicit none
    real*8 :: tcentroid
    real*8,intent(in) :: deltat
    integer,intent(in) :: lev
    integer,intent(in) :: i1in

    tcentroid = (dble(i1in)-sign(dble(0.5),dble(i1in))) * (2**lev)
  end function tcentroid


  function up(ilower)
    implicit none
    integer :: up
    integer,intent(in) :: ilower
    up = int((dble(ilower)+ 1.5*sign(1,ilower))/2.)
  end function up


  function down(iupper)
    implicit none
    integer :: down
    integer,intent(in) :: iupper
    down = iupper*2
  end function down


  function i4child(pi4in,nchild)
    implicit none
    integer,dimension(4) :: i4child
    integer,dimension(4),intent(in) :: pi4in
    integer,intent(in) :: nchild
    integer,dimension(4) :: d4i1
    integer :: i,d1i1
    d4i1(4) = int(real(nchild-1)/8.)
    d1i1 = nchild - d4i1(4)*8
    d4i1(3) = int(real(d1i1-1)/4.)
    d1i1 = d1i1 - d4i1(3)*4
    d4i1(2) = int(real(d1i1-1)/2.)
    d1i1 = d1i1 - d4i1(2)*2
    d4i1(1) = d1i1-1
    do i=1,4,1
      i4child(i) = down(pi4in(i)) - sign(1,pi4in(i)) + sign(d4i1(i),pi4in(i))
    end do
  end function i4child


  function i4childTLB(nchild)
    implicit none
    integer,dimension(4) :: i4childTLB
    integer,intent(in) :: nchild
    integer,dimension(4) :: d4i1
    integer :: i,d1i1
    d4i1(4) = int(real(nchild-1)/8.)
    d1i1 = nchild - d4i1(4)*8
    d4i1(3) = int(real(d1i1-1)/4.)
    d1i1 = d1i1 - d4i1(3)*4
    d4i1(2) = int(real(d1i1-1)/2.)
    d1i1 = d1i1 - d4i1(2)*2
    d4i1(1) = d1i1-1
    do i=1,4,1
      if (d4i1(i) == 1) then
        i4childTLB(i) = 1
      else if (d4i1(i) == 0) then
        i4childTLB(i) = -1
      else
        print*, 'Error in i4childTLB.'
      end if
    end do
  end function i4childTLB


  function i4parent(i4in)
    implicit none
    integer,dimension(4) :: i4parent
    integer,dimension(4),intent(in) :: i4in
    integer :: i
    do i=1,4,1
      i4parent(i) = up(i4in(i))
    end do
  end function i4parent


  function nchild(pi4in,ci4in)
    implicit none
    integer :: nchild
    integer,dimension(4),intent(in) :: pi4in,ci4in
    integer,dimension(4) :: d4i1
    integer :: i
    do i=1,4,1
      d4i1(i) = idchild(pi4in(i),ci4in(i))
    end do
    nchild = 1 + d4i1(1) + 2*d4i1(2) + 4*d4i1(3) + 8*d4i1(4)
  end function nchild


  function idchild(pi1in,ci1in)
    implicit none
    integer :: idchild
    integer,intent(in) :: pi1in,ci1in
    idchild = abs(ci1in) - (abs(down(pi1in))-1)
  end function idchild


  function TLBchild(ci4in)
    implicit none
    integer :: TLBchild
    integer,dimension(4),intent(in) :: ci4in
    integer,dimension(4) :: d4i1
    integer :: i
    do i=1,4,1
      d4i1(i) = int((ci4in(i)+1)/2)
    end do
    TLBchild = 1 + d4i1(1) + 2*d4i1(2) + 4*d4i1(3) + 8*d4i1(4)
  end function TLBchild


!  function box(deltat,r3in,itin)
!    use input, only: a
!    implicit none
!    integer,dimension(4) :: box
!    real*8,intent(in) :: deltat
!    real*8,dimension(3),intent(in) :: r3in
!    integer,intent(in) :: itin
!    integer :: i
!    do i=1,3,1
!      box(i) = int((r3in(i)/(a*deltat))+sign(dble(1.0),r3in(i)))
!    end do
!    box(4) = int(itin)
!  end function box

  function box4d(deltat,r4in)
    use input, only: a
    implicit none
    integer,dimension(4) :: box4d
    real*8,intent(in) :: deltat
    real*8,dimension(4),intent(in) :: r4in
    integer :: i

    do i=1,3,1
      if (r4in(i) > 0.) then
        if (r4in(i)/(a*deltat) > dble(int(r4in(i)/(a*deltat)))) then
          box4d(i) = int(r4in(i)/(a*deltat)) + 1
        else if (r4in(i)/(a*deltat) == dble(int(r4in(i)/(a*deltat)))) then
          box4d(i) = int(r4in(i)/(a*deltat))
        end if
      else if (r4in(i) == 0.) then
        box4d(i) = -1
      else if (r4in(i) < 0.) then
         if (r4in(i)/(a*deltat) < dble(int(r4in(i)/(a*deltat)))) then
          box4d(i) = int(r4in(i)/(a*deltat)) - 1
        else if (r4in(i)/(a*deltat) == dble(int(r4in(i)/(a*deltat)))) then
          box4d(i) = int(r4in(i)/(a*deltat))
        end if
      end if
    end do

    if (r4in(4) > 0) then
      if (r4in(4) > dble(int(r4in(4)))) then
        box4d(4) = int(r4in(4)) + 1
      else if (r4in(4) == dble(int(r4in(4)))) then
        box4d(4) = int(r4in(4))
      end if
    else if (r4in(4) == 0) then
      box4d(4) = -1
    else if (r4in(4) < 0) then
      if (r4in(4) < dble(int(r4in(4)))) then
        box4d(4) = int(r4in(4)) - 1
      else if (r4in(4) == dble(int(r4in(4)))) then
        box4d(4) = int(r4in(4))
      end if
    end if

!    do i=1,3,1
!      if (abs(r4in(i)/(a*deltat)) > abs(dble(int(r4in(i)/(a*deltat)))) .or. abs(r4in(i)/(a*deltat)) == 0.) then
!        box4d(i) = int(r4in(i)/(a*deltat))+int(sign(dble(1.0),r4in(i)))
!      else
!        box4d(i) = int((r4in(i)/(a*deltat)))
!      end if
!    end do
!    if (abs(r4in(4)) > abs(dble(int(r4in(4))))) then
!      box4d(4) = int(r4in(4)) + int(sign(dble(1.0),r4in(i)))
!    else
!      box4d(4) = int(r4in(4))
!    end if

  end function box4d

  function box1dt(deltat,r1in)
    use input, only: a
    implicit none
    integer :: box1dt
    real*8,intent(in) :: deltat
    real*8,intent(in) :: r1in
    if (r1in > 0) then
      if (r1in > dble(int(r1in))) then
        box1dt = int(r1in) + 1
      else if (r1in == dble(int(r1in))) then
        box1dt = int(r1in)
      end if
    else if (r1in == 0) then
      box1dt = -1
    else if (r1in < 0) then
      if (r1in < dble(int(r1in))) then
        box1dt = int(r1in) - 1
      else if (r1in == dble(int(r1in))) then
        box1dt = int(r1in)
      end if
    end if
  end function box1dt

end module fourdgridtypes

!###################################################################################################################################
!###################################################################################################################################


!module wingpanels
!  implicit none

!  type spoint
!    integer :: trj
!    real*8,dimension(3) :: ac
!    real*8,dimension(:,:),allocatable :: inrt
!    integer,dimension(2) :: ti
!    real*8,dimension(2) :: tr
!    logical :: ffl
!  end type spoint

!  type collpoint
!    real*8,dimension(3) :: pAC,velIN,velAC
!    real*8,dimension(:,:),allocatable :: pIN
!    integer,dimension(3) :: grid
!  end type

!  type wingpanel
!    type(collpoint) :: cpt
!    integer :: trj
!    integer,dimension(4) :: pan,neigh
!    character(len=1) :: ptype
!    real*8 :: area
!    real*8,dimension(3) :: norm
!  end type wingpanel

!end module wingpanels

!module wakepanels
!  implicit none

!  type point
!    real*8,dimension(:,:),allocatable :: pIN
!    real*8,dimension(3) :: pAC,vel
!    integer,dimension(3) :: grid
!    real*8 :: circ
!    integer,dimension(2) :: ti
!    real*8,dimension(2) :: tr
!    logical :: ffl
!  end type point

!  type wakepanel
!    type(point),dimension(4) :: p
!    integer :: trj
!    integer,dimension(4) :: neigh
!    real*8,dimension(3) :: pinf
!    real*8,dimension(:),allocatable :: circ
!    real*8,dimension(:,:),allocatable :: cIN
!  end type wakepanel

!end module wakepanels

!***********************************************************************************************************************************
!***********************************************************************************************************************************

!module gtmodinit
!  implicit none

!  type telem
!    logical :: active,ww
!    integer,dimension(3) :: parent
!    real*8,dimension(8,3) :: gridvtx
!    integer,dimension(8,3) :: grid
!  end type telem

!  type level
!    integer,dimension(3) :: nlevelem
!    type (telem),dimension(:,:,:),allocatable :: levelem
!  end type level

!  type gridstruct
!    integer :: nvtx
!    integer,dimension(:),allocatable :: vtx
!    logical :: ww !,active
!    integer,dimension(3) :: parent
!  end type gridstruct

!  contains
!    subroutine zerogrid(nl,ngrid,grids,grids_active,trees)
!      implicit none

!      integer,intent(in) :: nl
!      integer,dimension(3),intent(in) :: ngrid
!      type (level),dimension(:),intent(inout) :: trees
!      type (gridstruct),dimension(:,:,:),intent(inout) :: grids
!      logical,dimension(:,:,:),intent(inout) :: grids_active
!      integer :: j,k,l,m
!      integer,dimension(3) :: dummy3i


!      do k=1,trees(1)%nlevelem(1),1
!        do l=1,trees(1)%nlevelem(2),1
!          do m=1,trees(1)%nlevelem(3),1
!            if (trees(1)%levelem(k,l,m)%active) then
!              trees(1)%levelem(k,l,m)%active = .false.
!              trees(1)%levelem(k,l,m)%ww = .false.
!              trees(1)%levelem(k,l,m)%gridvtx = 0.
!              dummy3i = trees(1)%levelem(k,l,m)%parent
!              jloop: do j=1,nl-1,1
!                if (trees(j+1)%levelem(dummy3i(1),dummy3i(2),dummy3i(3))%active) then
!                  trees(j+1)%levelem(dummy3i(1),dummy3i(2),dummy3i(3))%active = .false.
!                  trees(j+1)%levelem(dummy3i(1),dummy3i(2),dummy3i(3))%ww = .false.
!                  trees(j+1)%levelem(dummy3i(1),dummy3i(2),dummy3i(3))%gridvtx = 0.
!                  dummy3i = trees(j+1)%levelem(dummy3i(1),dummy3i(2),dummy3i(3))%parent
!                else
!                  exit jloop
!                end if
!              end do jloop
!            end if
!          end do
!        end do
!      end do

!      do k=1,ngrid(1),1
!        do l=1,ngrid(2),1
!          do m=1,ngrid(3),1
!            if (grids_active(k,l,m)) then
!              grids_active(k,l,m) = .false.
!              grids(k,l,m)%ww = .false.
!              grids(k,l,m)%nvtx = 0
!            end if
!          end do
!        end do
!      end do

!    end subroutine zerogrid

!end module gtmodinit

!***********************************************************************************************************************************
!***********************************************************************************************************************************

module vtxparticle
  implicit none

  type vtxpart
    logical,dimension(:),allocatable :: active
    real*8 :: Ring
    real*8,dimension(3) :: posAC,mAC,vel
    real*8,dimension(:,:),allocatable :: posIN,mIN
    integer,dimension(3) :: grid
    real*8,dimension(3,3) :: vtensor
  end type vtxpart

  contains

  subroutine advvtx(deltat,vind,gradvind,vpos,vtx)
    use input, only: avps, pi
    use maths, only: modulus
    implicit none 

    real*8 :: deltat
    real*8,dimension(3),intent(in) :: vind
    real*8,dimension(3,3),intent(in) :: gradvind
    real*8,dimension(3),intent(inout) :: vpos,vtx
    real*8,dimension(3) :: dvtx

    if (modulus(vind) > 340.) then

      vpos = vpos + vind/modulus(vind) * 340. * deltat

    else

      vpos = vpos + vind * deltat

    end if

    !transpose
!    vtx(1) = vtx(1) + deltat * (gradvind(1,1) * vtx(1) + gradvind(2,1) * vtx(2) + gradvind(3,1) * vtx(3))
!    vtx(2) = vtx(2) + deltat * (gradvind(1,2) * vtx(1) + gradvind(2,2) * vtx(2) + gradvind(3,2) * vtx(3))
!    vtx(3) = vtx(3) + deltat * (gradvind(1,3) * vtx(1) + gradvind(2,3) * vtx(2) + gradvind(3,3) * vtx(3))

    !normal
!    vtx(1) = vtx(1) + deltat * (gradvind(1,1) * vtx(1) + gradvind(1,2) * vtx(2) + gradvind(1,3) * vtx(3))
!    vtx(2) = vtx(2) + deltat * (gradvind(2,1) * vtx(1) + gradvind(2,2) * vtx(2) + gradvind(2,3) * vtx(3))
!    vtx(3) = vtx(3) + deltat * (gradvind(3,1) * vtx(1) + gradvind(3,2) * vtx(2) + gradvind(3,3) * vtx(3))

    !mixed

    dvtx(1) = deltat * &
           & (gradvind(1,1) * vtx(1) + (gradvind(1,2)+gradvind(2,1))*0.5 * vtx(2) + (gradvind(1,3)+gradvind(3,1))*0.5 * vtx(3))
    dvtx(2) = deltat * &
           & ((gradvind(2,1)+gradvind(1,2))*0.5 * vtx(1) + gradvind(2,2) * vtx(2) + (gradvind(2,3)+gradvind(3,2))*0.5 * vtx(3))
    dvtx(3) = deltat * &
           & ((gradvind(3,1)+gradvind(1,3))*0.5 * vtx(1) + (gradvind(3,2)+gradvind(2,3))*0.5 * vtx(2) + gradvind(3,3) * vtx(3))

    if (modulus(dvtx+vtx)/(4.*pi*avps**2) > 340.) then
      vtx = (dvtx+vtx)/modulus(dvtx+vtx) * 340. * (4.*pi*avps**2)
    else
      vtx = dvtx + vtx
    end if


!    vtx(1) = vtx(1) + deltat * &
!           & (gradvind(1,1) * vtx(1) + (gradvind(1,2)+gradvind(2,1))*0.5 * vtx(2) + (gradvind(1,3)+gradvind(3,1))*0.5 * vtx(3))
!    vtx(2) = vtx(2) + deltat * &
!           & ((gradvind(2,1)+gradvind(1,2))*0.5 * vtx(1) + gradvind(2,2) * vtx(2) + (gradvind(2,3)+gradvind(3,2))*0.5 * vtx(3))
!    vtx(3) = vtx(3) + deltat * &
!           & ((gradvind(3,1)+gradvind(1,3))*0.5 * vtx(1) + (gradvind(3,2)+gradvind(2,3))*0.5 * vtx(2) + gradvind(3,3) * vtx(3))

  end subroutine advvtx

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

!subroutine findingrid(j,act,pos,gridpos,nl,dbox,ngrid,gridext,grids,grids_active,trees)
!  use maths, only: inc
!  use gtmodinit
!  implicit none

!  type (level),dimension(:),intent(inout) :: trees
!  type (gridstruct),dimension(:,:,:),intent(inout) :: grids
!  logical,dimension(:,:,:),intent(inout) :: grids_active
!  logical,intent(inout) :: act
!  integer,dimension(3),intent(inout) :: gridpos
!  integer,intent(in) :: j,nl
!  integer,dimension(3),intent(in) :: ngrid
!  real*8,dimension(3,2),intent(in) :: gridext
!  real*8,dimension(3),intent(in) :: dbox,pos

!  integer,dimension(3) :: dummy3i1,dummy3i2
!  logical :: check3
!  integer :: l,dummy1i1

!  if (act) then
!    !check whether in same box
!    dummy3i1 = 0
!    if (incheck6(pos,gridpos,dbox,gridext(:,1)) .eqv. .false.) then
!      !if not, find 3D
!      call findpos3D(pos,ngrid,gridext(:,1),dbox,check3,dummy3i1)
!      gridpos = dummy3i1
!      if (check3) then
!        grids_active(dummy3i1(1),dummy3i1(2),dummy3i1(3)) = .true.
!      else
!        act = .false.
!        print*, 'Existing vtx drifted out of bounding domain.'
!      end if
!    else
!      grids_active(gridpos(1),gridpos(2),gridpos(3)) = .true.
!      dummy3i1 = gridpos
!    end if

!    if (act) then
!      !setactive
!      do l=1,3,1
!        dummy3i2(l) = int((gridpos(l)+1)*0.5)
!      end do

!      do l=1,nl,1
!        if (trees(l)%levelem(dummy3i2(1),dummy3i2(2),dummy3i2(3))%active .eqv. .false.) then
!          trees(l)%levelem(dummy3i2(1),dummy3i2(2),dummy3i2(3))%active = .true.
!          dummy3i2 = trees(l)%levelem(dummy3i2(1),dummy3i2(2),dummy3i2(3))%parent
!        else
!          exit
!        end if
!      end do

!      call inc(grids(dummy3i1(1),dummy3i1(2),dummy3i1(3))%nvtx)
!      dummy1i1 = grids(dummy3i1(1),dummy3i1(2),dummy3i1(3))%nvtx
!      grids(dummy3i1(1),dummy3i1(2),dummy3i1(3))%vtx(dummy1i1) = j
!    end if
!  end if

!end subroutine findingrid
!!-----------------------------------------------------------------------------------------------------------------------------------

!subroutine findingridww(j,act,pos,gridpos,nl,dbox,ngrid,gridext,grids,grids_active,trees)
!  use maths, only: inc
!  use gtmodinit
!  implicit none

!  type (level),dimension(:),intent(inout) :: trees
!  type (gridstruct),dimension(:,:,:),intent(inout) :: grids
!  logical,dimension(:,:,:),intent(inout) :: grids_active
!  logical,intent(inout) :: act
!  integer,dimension(3),intent(inout) :: gridpos
!  integer,intent(in) :: j,nl
!  integer,dimension(3),intent(in) :: ngrid
!  real*8,dimension(3,2),intent(in) :: gridext
!  real*8,dimension(3),intent(in) :: dbox,pos

!  integer,dimension(3) :: dummy3i1,dummy3i2
!  logical :: check3
!  integer :: l,dummy1i1

!  if (act) then
!    !check whether in same box
!    dummy3i1 = 0
!    if (incheck6(pos,gridpos,dbox,gridext(:,1)) .eqv. .false.) then
!      !if not, find 3D
!      call findpos3D(pos,ngrid,gridext(:,1),dbox,check3,dummy3i1)
!      gridpos = dummy3i1
!      if (check3) then
!        grids(dummy3i1(1),dummy3i1(2),dummy3i1(3))%ww = .true.
!        grids_active(dummy3i1(1),dummy3i1(2),dummy3i1(3)) = .true.
!      else
!        act = .false.
!        print*, 'Existing vtx drifted out of bounding domain.'
!      end if
!    else
!      grids(gridpos(1),gridpos(2),gridpos(3))%ww = .true.
!      grids_active(gridpos(1),gridpos(2),gridpos(3)) = .true.
!      dummy3i1 = gridpos
!    end if

!    if (act) then
!      !setactive
!      do l=1,3,1
!        dummy3i2(l) = int((gridpos(l)+1)*0.5)
!      end do

!      do l=1,nl,1
!        if (trees(l)%levelem(dummy3i2(1),dummy3i2(2),dummy3i2(3))%ww .eqv. .false.) then
!          trees(l)%levelem(dummy3i2(1),dummy3i2(2),dummy3i2(3))%ww = .true.
!          trees(l)%levelem(dummy3i2(1),dummy3i2(2),dummy3i2(3))%active = .true.
!          dummy3i2 = trees(l)%levelem(dummy3i2(1),dummy3i2(2),dummy3i2(3))%parent
!        else
!          exit
!        end if
!      end do

!      call inc(grids(dummy3i1(1),dummy3i1(2),dummy3i1(3))%nvtx)
!      dummy1i1 = grids(dummy3i1(1),dummy3i1(2),dummy3i1(3))%nvtx
!      grids(dummy3i1(1),dummy3i1(2),dummy3i1(3))%vtx(dummy1i1) = j
!    end if
!  end if

!end subroutine findingridww


!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

subroutine newvtxpart(p1,p2,p3,p4,c1,c2,c3,c4,pout,mout)
  implicit none

  real*8,dimension(4),intent(in) :: p1,p2,p3,p4
  real*8,intent(in) :: c1,c2,c3,c4
  real*8,dimension(4) :: v1,v2,v3,v4
  real*8 :: v1mod,v2mod,v3mod,v4mod,circsum
  real*8,dimension(3),intent(out) :: mout
  real*8,dimension(4),intent(out) :: pout

  v1 = c1*(p2-p1)
  v2 = c2*(p3-p2)
  v3 = c3*(p4-p3)
  v4 = c4*(p1-p4)

  v1mod = sqrt(dot_product(v1(1:3),v1(1:3)))
  v2mod = sqrt(dot_product(v2(1:3),v2(1:3)))
  v3mod = sqrt(dot_product(v3(1:3),v3(1:3)))
  v4mod = sqrt(dot_product(v4(1:3),v4(1:3)))

  mout = v1(1:3) + v2(1:3) + v3(1:3) + v4(1:3)
  circsum = v1mod+v2mod+v3mod+v4mod

!  pout = (p1+p2+p3+p4)*0.25

  pout = 0.5*((p1+p2)*v1mod + (p2+p3)*v2mod + (p3+p4)*v3mod + (p4+p1)*v4mod)/circsum

end subroutine newvtxpart

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!!Collapses a line of wake panels to a vortex particle.
!subroutine vtxsetcollapse(np,npan,points,panels,neigh,ptype,circ,vpos,vtx)
!  use wingpanels
!  implicit none

!  integer,intent(in) :: np,npan
!  type(spoint),dimension(np),intent(in) :: points
!  real*8,dimension(np),intent(in) :: circ
!  integer,dimension(npan,4),intent(in) :: panels,neigh
!  character(len=1),dimension(npan) :: ptype
!  real*8,dimension(npan,4) :: c
!  real*8,dimension(npan,3),intent(out) :: vpos,vtx
!  real*8,dimension(3),parameter :: zerovec = (/0.,0.,0./)
!  integer :: i,j

!  c = 0.

!  do i=1,npan,1
!    do j=1,4,1
!      if (neigh(i,j) > 0) then
!        c(i,j) = 0.5*(circ(i) - circ(neigh(i,j)))
!      else if (neigh(i,j) == 0) then
!        c(i,j) = 0.
!      else if (neigh(i,j) == -1) then
!        c(i,j) = circ(i)
!      else if (neigh(i,j) == -2) then
!        c(i,j) = 0.
!      else if (neigh(i,j) == -3) then
!        !Not sure about this as body junction
!        c(i,j) = circ(i)
!      else if (neigh(i,j) == -4) then
!        c(i,j) = circ(i)
!      end if
!    end do
!  end do

!  do i=1,npan,1
!    if (ptype(i) == 'Q') then
!      call newvtxpart(points(panels(i,1))%inrt(:,t),points(panels(i,2))%inrt(:,t),points(panels(i,3))%inrt(:,t),points(panels(i,4))%inrt(:,t), &
!         & c(i,1),c(i,2),c(i,3),c(i,4),vpos(i,:),vtx(i,:))
!    else if (ptype(i) == 'T') then
!      call newvtxpart(points(panels(i,1))%inrt(:,t),points(panels(i,2))%inrt(:,t),points(panels(i,3))%inrt(:,t),zerovec, &
!         & c(i,1),c(i,2),c(i,3),c(i,4),vpos(i,:),vtx(i,:))
!    end if
!  end do

!end subroutine vtxsetcollapse

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!!Collapses a line of wake panels to a vortex particle.
!subroutine vtxsetcollapse2(t,np,npan,points,wingpans,circ,vpos,vtx)
!  use wingpanels
!  implicit none

!  integer,intent(in) :: np,npan,t
!  type(spoint),dimension(np),intent(in) :: points
!  type(wingpanel),dimension(npan),intent(in) :: wingpans
!  real*8,dimension(npan),intent(in) :: circ

!  real*8,dimension(npan,4) :: c
!  real*8,dimension(3,npan),intent(out) :: vpos,vtx
!  real*8,dimension(3),parameter :: zerovec = (/0.,0.,0./)
!  integer :: i,j

!  c = 0.

!  do i=1,npan,1
!    do j=1,4,1
!      if (wingpans(i)%neigh(j) > 0) then
!        c(i,j) = 0.5*(circ(i) - circ(wingpans(i)%neigh(j)))
!      else if (wingpans(i)%neigh(j) == 0) then
!        c(i,j) = 0.
!      else if (wingpans(i)%neigh(j) == -1) then
!        c(i,j) = circ(i)
!      else if (wingpans(i)%neigh(j) == -2) then
!        c(i,j) = 0.
!      else if (wingpans(i)%neigh(j) == -3) then
!        !Not sure about this as body junction
!        c(i,j) = circ(i)
!      else if (wingpans(i)%neigh(j) == -4) then
!        c(i,j) = circ(i)
!      end if
!    end do
!  end do

!  do i=1,npan,1
!    if (wingpans(i)%ptype == 'Q') then
!      call newvtxpart(points(wingpans(i)%pan(1))%inrt(:,t),points(wingpans(i)%pan(2))%inrt(:,t), &
!        & points(wingpans(i)%pan(3))%inrt(:,t), &
!        & points(wingpans(i)%pan(4))%inrt(:,t),c(i,1),c(i,2),c(i,3),c(i,4),vpos(:,i),vtx(:,i))
!    else if (wingpans(i)%ptype == 'T') then
!      call newvtxpart(points(wingpans(i)%pan(1))%inrt(:,t),points(wingpans(i)%pan(2))%inrt(:,t), &
!        & points(wingpans(i)%pan(3))%inrt(:,t),zerovec, &
!        & c(i,1),c(i,2),c(i,3),c(i,4),vpos(:,i),vtx(:,i))
!    end if
!  end do

!end subroutine vtxsetcollapse2

!!-----------------------------------------------------------------------------------------------------------------------------------
!!-----------------------------------------------------------------------------------------------------------------------------------
!!Collapses a set of panels to vortex particles.
!subroutine vtxwakecollapse6(nwakepan,nwpanels,t,wakepans,vpos,vtx)
!  use maths, only: cycperm
!  use wakepanels
!  implicit none

!  integer,intent(in) :: nwakepan,nwpanels,t
!  type(wakepanel),dimension(:),intent(inout) :: wakepans
!  real*8,dimension(3,(nwpanels+2)*nwakepan),intent(out) :: vpos,vtx
!  integer :: i,j,k,dummy1i1

!  do i=1,min(nwpanels+2,t),1
!    do j=1,nwakepan,1
!      dummy1i1 = (i-1)*nwakepan + j
!      do k=1,4,1
!        !fail safe
!        wakepans(dummy1i1)%p(k)%circ = 0.
!        if (wakepans(dummy1i1)%neigh(k) > 0) then
!          wakepans(dummy1i1)%p(k)%circ = 0.5*(wakepans(dummy1i1)%circ(t) - wakepans(wakepans(dummy1i1)%neigh(k))%circ(t))
!        else if (wakepans(dummy1i1)%neigh(k) == 0 .or. wakepans(dummy1i1)%neigh(k) == -6) then
!          wakepans(dummy1i1)%p(k)%circ = wakepans(dummy1i1)%circ(t)
!        else if (wakepans(dummy1i1)%neigh(k) == -2) then
!          wakepans(dummy1i1)%p(k)%circ = 0.
!        else if (wakepans(dummy1i1)%neigh(k) == -1) then
!          print*, 'check vtxwakecollapse3 in gtmod8 (1), check wake neighbour -1'
!        else if (wakepans(dummy1i1)%neigh(k) == -5) then
!        !see Kutta condition in vtxwingsegring
!          wakepans(dummy1i1)%p(k)%circ = 0.
!        end if
!      end do
!    end do
!  end do

!  do i=1,min(nwpanels+2,t),1
!    do j=1,nwakepan,1
!      dummy1i1 = (i-1)*nwakepan + j
!      call newvtxpart(wakepans(dummy1i1)%p(1)%pIN,wakepans(dummy1i1)%p(2)%pIN,wakepans(dummy1i1)%p(3)%pIN, &
!        & wakepans(dummy1i1)%p(4)%pIN,wakepans(dummy1i1)%p(1)%circ,wakepans(dummy1i1)%p(2)%circ,wakepans(dummy1i1)%p(3)%circ, &
!        & wakepans(dummy1i1)%p(4)%circ,vpos(:,dummy1i1),vtx(:,dummy1i1))
!    end do
!  end do

!end subroutine vtxwakecollapse6

!!Collapses a set of panels to vortex particles.
!subroutine vtxwakecollapse3(nwakepan,nwpanels,tstep,wakepans,vpos,vtx)
!!subroutine vtxwakecollapse3(nwakepan,nwpanels,tstep,ifp,ilp,wakepans,wakeneigh,vpos,vtx)
!  use maths, only: cycperm
!  use wakepanels
!  implicit none

!  integer,intent(in) :: nwakepan,nwpanels,tstep !,ifp,ilp
!  type(wakepanel),dimension(:),intent(inout) :: wakepans
!!  integer,dimension(nwakepan,4),intent(in) :: wakeneigh
!  real*8,dimension(3,(nwpanels+2)*nwakepan),intent(out) :: vpos,vtx
!  integer :: i,j,k,dummy1i1  !,dummy1i2,nlim

!  do i=1,min(nwpanels+2,tstep),1
!    do j=1,nwakepan,1
!      dummy1i1 = (i-1)*nwakepan + j
!      do k=1,4,1
!        !fail safe
!        wakepans(dummy1i1)%p(k)%circ = 0.
!        if (wakepans(dummy1i1)%neigh(k) > 0) then
!          wakepans(dummy1i1)%p(k)%circ = 0.5*(wakepans(dummy1i1)%circ - wakepans(wakepans(dummy1i1)%neigh(k))%circ)
!        else if (wakepans(dummy1i1)%neigh(k) == 0 .or. wakepans(dummy1i1)%neigh(k) == -6) then
!          wakepans(dummy1i1)%p(k)%circ = wakepans(dummy1i1)%circ
!        else if (wakepans(dummy1i1)%neigh(k) == -2) then
!          wakepans(dummy1i1)%p(k)%circ = 0.
!        else if (wakepans(dummy1i1)%neigh(k) == -1) then
!          print*, 'check vtxwakecollapse3 in gtmod8 (1), check wake neighbour -1'
!        else if (wakepans(dummy1i1)%neigh(k) == -5) then
!          wakepans(dummy1i1)%p(k)%circ = 0.
!        end if
!      end do
!    end do
!  end do

!  do i=1,min(nwpanels+2,tstep),1
!    do j=1,nwakepan,1
!      dummy1i1 = (i-1)*nwakepan + j
!      call newvtxpart(wakepans(dummy1i1)%p(1)%pIN,wakepans(dummy1i1)%p(2)%pIN,wakepans(dummy1i1)%p(3)%pIN, &
!        & wakepans(dummy1i1)%p(4)%pIN,wakepans(dummy1i1)%p(1)%circ,wakepans(dummy1i1)%p(2)%circ,wakepans(dummy1i1)%p(3)%circ, &
!        & wakepans(dummy1i1)%p(4)%circ,vpos(:,dummy1i1),vtx(:,dummy1i1))
!    end do
!  end do

!end subroutine vtxwakecollapse3

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!!Collapses a set of panels to vortex particles. except ifp -> dipole
!subroutine vtxwakecollapse5(nwakepan,nwpanels,tstep,ifp,ifpm1,wakepans,vpos,vtx)
!  use maths, only: cycperm
!  use wakepanels
!  implicit none

!  integer,intent(in) :: nwakepan,nwpanels,tstep,ifp,ifpm1
!  type(wakepanel),dimension((nwpanels+2)*nwakepan),intent(inout) :: wakepans
!  real*8,dimension(3,(nwpanels+2)*nwakepan),intent(out) :: vpos,vtx
!  integer :: i,j,k,dummy1i1,ni

!  do i=1,min(nwpanels+2,tstep)-1,1
!    ni = cycperm(nwpanels+2,(nwpanels+2)-i,ifp)
!    do j=1,nwakepan,1
!      dummy1i1 = (ni-1)*nwakepan + j
!      do k=1,4,1
!        !fail safe
!        wakepans(dummy1i1)%p(k)%circ = 0.
!        if (ni /= ifpm1) then
!          if (wakepans(dummy1i1)%neigh(k) > 0) then
!            wakepans(dummy1i1)%p(k)%circ = 0.5*(wakepans(dummy1i1)%circ - wakepans(wakepans(dummy1i1)%neigh(k))%circ)
!          else if (wakepans(dummy1i1)%neigh(k) == 0 .or. wakepans(dummy1i1)%neigh(k) == -6) then
!            wakepans(dummy1i1)%p(k)%circ = wakepans(dummy1i1)%circ
!          else if (wakepans(dummy1i1)%neigh(k) == -2) then
!            wakepans(dummy1i1)%p(k)%circ = 0.
!          else if (wakepans(dummy1i1)%neigh(k) == -1) then
!            print*, 'check vtxwakecollapse5 in gtmod8 (1), check wake neighbour -1'
!          else if (wakepans(dummy1i1)%neigh(k) == -5) then
!            print*, 'check vtxwakecollapse5 in gtmod8 (1), ifp shouldn not be cycled over.'
!          end if
!        else
!          if (wakepans(dummy1i1)%neigh(k) > (ifp-1)*nwakepan .and. wakepans(dummy1i1)%neigh(k) <= ifp*nwakepan) then
!            wakepans(dummy1i1)%p(k)%circ = 0.
!          else if (wakepans(dummy1i1)%neigh(k) > 0) then
!            wakepans(dummy1i1)%p(k)%circ = 0.5*(wakepans(dummy1i1)%circ - wakepans(wakepans(dummy1i1)%neigh(k))%circ)
!          else if (wakepans(dummy1i1)%neigh(k) == 0 .or. wakepans(dummy1i1)%neigh(k) == -6) then
!            wakepans(dummy1i1)%p(k)%circ = wakepans(dummy1i1)%circ
!          else if (wakepans(dummy1i1)%neigh(k) == -2) then
!            wakepans(dummy1i1)%p(k)%circ = 0.
!          else if (wakepans(dummy1i1)%neigh(k) == -1) then
!            print*, 'check vtxwakecollapse5 in gtmod8 (1), check wake neighbour -1'
!          else if (wakepans(dummy1i1)%neigh(k) == -5) then
!            print*, 'check vtxwakecollapse5 in gtmod8 (1), ifp shouldn not be cycled over.'
!          end if
!        end if
!      end do
!    end do
!  end do

!  do i=1,min(nwpanels+2,tstep)-1,1
!    ni = cycperm(nwpanels+2,(nwpanels+2)-i,ifp)
!    do j=1,nwakepan,1
!      dummy1i1 = (ni-1)*nwakepan + j
!      call newvtxpart(wakepans(dummy1i1)%p(1)%pIN,wakepans(dummy1i1)%p(2)%pIN,wakepans(dummy1i1)%p(3)%pIN, &
!        & wakepans(dummy1i1)%p(4)%pIN,wakepans(dummy1i1)%p(1)%circ,wakepans(dummy1i1)%p(2)%circ,wakepans(dummy1i1)%p(3)%circ, &
!        & wakepans(dummy1i1)%p(4)%circ,vpos(:,dummy1i1),vtx(:,dummy1i1))
!    end do
!  end do

!!  do i=1,min(nwpanels+2,tstep)-1,1
!!    ni = cycperm(nwpanels+2,(nwpanels+2)-i,ifp)
!!    do j=1,nwakepan,1
!!      dummy1i1 = (ni-1)*nwakepan + j
!!      do k=1,4,1
!!        wakepans(dummy1i1)%p(k)%circ = 0.
!!        if (wakeneigh(j,k) > 0) then
!!          dummy1i2 = (ni-1)*nwakepan + wakeneigh(j,k)
!!          wakepans(dummy1i1)%p(k)%circ = 0.5*(wakepans(dummy1i1)%circ - wakepans(dummy1i2)%circ)
!!        else if (wakeneigh(j,k) == 0) then
!!          wakepans(dummy1i1)%p(k)%circ = wakepans(dummy1i1)%circ
!!        else if (wakeneigh(j,k) == -1) then
!!          wakepans(dummy1i1)%p(k)%circ = wakepans(dummy1i1)%circ
!!        else if (wakeneigh(j,k) == -2) then
!!          if (ni==ilp) then
!!            wakepans(dummy1i1)%p(k)%circ = 0.
!!          else
!!            dummy1i2 = (cycperm(nwpanels+2,nwpanels+1,ni)-1)*nwakepan + j
!!            wakepans(dummy1i1)%p(k)%circ = 0.5*(wakepans(dummy1i1)%circ - wakepans(dummy1i2)%circ)
!!          end if
!!        else if (wakeneigh(j,k) == -5) then

!!          if (ni==cycperm(nwpanels+2,nwpanels+1,ifp)) then
!!            wakepans(dummy1i1)%p(k)%circ = wakepans(dummy1i1)%circ
!!          else
!!            dummy1i2 = (cycperm(nwpanels+2,1,ni)-1)*nwakepan + j
!!            wakepans(dummy1i1)%p(k)%circ = 0.5*(wakepans(dummy1i1)%circ - wakepans(dummy1i2)%circ)
!!          end if
!!        else if (wakeneigh(j,k) == -6) then
!!          if (ni==ilp) then
!!            wakepans(dummy1i1)%p(k)%circ = wakepans(dummy1i1)%circ
!!          else
!!            dummy1i2 = (cycperm(nwpanels+2,nwpanels+1,ni)-1)*nwakepan + j
!!            wakepans(dummy1i1)%p(k)%circ = 0.5*(wakepans(dummy1i1)%circ - wakepans(dummy1i2)%circ)
!!          end if
!!        end if
!!      end do
!!    end do
!!  end do

!!  do i=1,min(nwpanels+2,tstep)-1,1
!!    ni = cycperm(nwpanels+2,(nwpanels+2)-i,ifp)
!!    do j=1,nwakepan,1
!!      dummy1i1 = (ni-1)*nwakepan + j
!!      call newvtxpart(wakepans(dummy1i1)%p(1)%pIN,wakepans(dummy1i1)%p(2)%pIN,wakepans(dummy1i1)%p(3)%pIN, &
!!        & wakepans(dummy1i1)%p(4)%pIN,wakepans(dummy1i1)%p(1)%circ,wakepans(dummy1i1)%p(2)%circ,wakepans(dummy1i1)%p(3)%circ, &
!!        & wakepans(dummy1i1)%p(4)%circ,vpos(dummy1i1,:),vtx(dummy1i1,:))
!!    end do
!!  end do

!end subroutine vtxwakecollapse5

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!!Collapses a set of panels to vortex particles.
!subroutine vtxwakecollapse2(nwakepan,nwpanels,ifp,ilp,wakepans,wakeneigh,vpos,vtx)
!  use maths, only: cycperm
!  use wakepanels
!  implicit none

!  integer,intent(in) :: nwakepan,nwpanels,ifp,ilp
!  type(wakepanel),dimension(:),intent(in) :: wakepans
!  integer,dimension(nwakepan,4),intent(in) :: wakeneigh
!  real*8,dimension((nwpanels+2)*nwakepan,4) :: c
!  real*8,dimension((nwpanels+2)*nwakepan,3),intent(out) :: vpos,vtx
!  integer :: i,j,k,dummy1i1,dummy1i2

!  c = 0.

!  do i=1,nwpanels+2,1
!    do j=1,nwakepan,1
!      dummy1i1 = (i-1)*nwakepan + j
!      do k=1,4,1
!        if (wakeneigh(j,k) > 0) then
!          dummy1i2 = (i-1)*nwakepan + wakeneigh(j,k)
!          c(dummy1i1,k) = 0.5*(wakepans(dummy1i1)%circ - wakepans(dummy1i2)%circ)
!        else if (wakeneigh(j,k) == 0) then
!          c(dummy1i1,k) = wakepans(dummy1i1)%circ
!        else if (wakeneigh(j,k) == -1) then
!          c(dummy1i1,k) = wakepans(dummy1i1)%circ
!        else if (wakeneigh(j,k) == -2) then
!          if (i==ilp) then
!            c(dummy1i1,k) = 0.
!          else
!            dummy1i2 = (cycperm(nwpanels+2,nwpanels+1,i)-1)*nwakepan + j
!            c(dummy1i1,k) = 0.5*(wakepans(dummy1i1)%circ - wakepans(dummy1i2)%circ)
!          end if
!        else if (wakeneigh(j,k) == -5) then
!          if (i==ifp) then
!            c(dummy1i1,k) = 0.
!          else
!            dummy1i2 = (cycperm(nwpanels+2,1,i)-1)*nwakepan + j
!            c(dummy1i1,k) = 0.5*(wakepans(dummy1i1)%circ - wakepans(dummy1i2)%circ)
!          end if
!        else if (wakeneigh(j,k) == -6) then
!          if (i==ilp) then
!            c(dummy1i1,k) = wakepans(dummy1i1)%circ
!          else
!            dummy1i2 = (cycperm(nwpanels+2,nwpanels+1,i)-1)*nwakepan + j
!            c(dummy1i1,k) = 0.5*(wakepans(dummy1i1)%circ - wakepans(dummy1i2)%circ)
!          end if
!        end if
!      end do
!    end do
!  end do

!  do i=1,nwpanels+2,1
!    do j=1,nwakepan,1
!      dummy1i1 = (i-1)*nwakepan + j
!      call newvtxpart(wakepans(dummy1i1)%p(1)%pIN,wakepans(dummy1i1)%p(2)%pIN,wakepans(dummy1i1)%p(3)%pIN, &
!        & wakepans(dummy1i1)%p(4)%pIN,c(dummy1i1,1),c(dummy1i1,2),c(dummy1i1,3),c(dummy1i1,4),vpos(dummy1i1,:),vtx(dummy1i1,:))
!    end do
!  end do

!end subroutine vtxwakecollapse2

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!!Collapses a set of panels to vortex particles.
!subroutine vtxwakecollapse(nwakepan,nwpanels,count1,firstpan,wakepan,lastpan,wakeneigh,fcirc,wcirc,lcirc,vpos,vtx)
!  use maths, only: cycperm
!  implicit none

!  integer,intent(in) :: nwakepan,nwpanels,count1
!  real*8,dimension(nwakepan,4,3),intent(in) :: firstpan,lastpan
!  real*8,dimension(nwpanels,nwakepan,4,3),intent(in) :: wakepan
!  real*8,dimension(nwakepan),intent(in) :: fcirc,lcirc
!  real*8,dimension(nwpanels,nwakepan),intent(in) :: wcirc
!  integer,dimension(nwakepan,4),intent(inout) :: wakeneigh
!  real*8,dimension((nwpanels+2)*nwakepan,4) :: c
!  real*8,dimension((nwpanels+2)*nwakepan,3),intent(out) :: vpos,vtx
!  integer :: i,j,k

!  c = 0.

!  do i=1,nwakepan,1
!    do j=1,4,1
!      if (wakeneigh(i,j) > 0) then
!        c(i,j) = 0.5*(fcirc(i) - fcirc(wakeneigh(i,j)))
!        c(nwakepan*(nwpanels+1)+i,j) = 0.5*(lcirc(i) - lcirc(wakeneigh(i,j)))
!        c(nwakepan+i,j) = 0.5*(wcirc(count1,i) - wcirc(count1,wakeneigh(i,j)))
!        c(nwakepan*nwpanels+i,j) = 0.5*(wcirc(cycperm(nwpanels,nwpanels-1,count1),i) - &
!                                     & wcirc(cycperm(nwpanels,nwpanels-1,count1),wakeneigh(i,j)))
!      else if (wakeneigh(i,j) == 0) then
!        c(i,j) = fcirc(i)
!        c(nwakepan*(nwpanels+1)+i,j) = lcirc(i)
!        c(nwakepan+i,j) = wcirc(count1,i)
!        c(nwakepan*nwpanels+i,j) =wcirc(cycperm(nwpanels,nwpanels-1,count1),i)
!      else if (wakeneigh(i,j) == -1) then
!        c(i,j) = fcirc(i)
!        c(nwakepan*(nwpanels+1)+i,j) = lcirc(i)
!        c(nwakepan+i,j) = wcirc(count1,i)
!        c(nwakepan*nwpanels+i,j) =wcirc(cycperm(nwpanels,nwpanels-1,count1),i)
!      else if (wakeneigh(i,j) == -2) then
!        c(i,j) = 0.5*(fcirc(i) - wcirc(count1,i))
!        c(nwakepan*(nwpanels+1)+i,j) = 0.
!        c(nwakepan+i,j) = 0.5*(wcirc(count1,i) - wcirc(cycperm(nwpanels,1,count1),i))
!        c(nwakepan*nwpanels+i,j) = 0.5*(wcirc(cycperm(nwpanels,nwpanels-1,count1),i) - &
!                                     & lcirc(i))
!      else if (wakeneigh(i,j) == -5) then
!        c(i,j) = 0.
!        c(nwakepan*(nwpanels+1)+i,j) = 0.5*(lcirc(i) - wcirc(cycperm(nwpanels,nwpanels-1,count1),i))
!        c(nwakepan+i,j) = 0.5*(wcirc(count1,i) - fcirc(i))
!        c(nwakepan*nwpanels+i,j) = 0.5*(wcirc(cycperm(nwpanels,nwpanels-1,count1),i) - &
!                                     & wcirc(cycperm(nwpanels,nwpanels-2,count1),i))
!      else if (wakeneigh(i,j) == -6) then
!        c(i,j) = 0.5*(fcirc(i) - wcirc(count1,i))
!        c(nwakepan*(nwpanels+1)+i,j) = lcirc(i)
!        c(nwakepan+i,j) = 0.5*(wcirc(count1,i) - wcirc(cycperm(nwpanels,1,count1),i))
!        c(nwakepan*nwpanels+i,j) = 0.5*(wcirc(cycperm(nwpanels,nwpanels-1,count1),i) - &
!                                     & lcirc(wakeneigh(i,j)))
!      end if
!    end do
!  end do

!  do k=2,nwpanels-1,1
!  do i=1,nwakepan,1
!    do j=1,4,1
!      if (wakeneigh(i,j) > 0) then
!        c(k*nwakepan+i,j) = 0.5*(wcirc(cycperm(nwpanels,k-1,count1),i) - wcirc(cycperm(nwpanels,k-1,count1),wakeneigh(i,j)))
!      else if (wakeneigh(i,j) == 0) then
!        c(k*nwakepan+i,j) = wcirc(cycperm(nwpanels,k-1,count1),i)
!      else if (wakeneigh(i,j) == -1) then
!        c(k*nwakepan+i,j) = wcirc(cycperm(nwpanels,k-1,count1),i)
!      else if (wakeneigh(i,j) == -2) then
!        c(k*nwakepan+i,j) = 0.5*(wcirc(cycperm(nwpanels,k-1,count1),i) - wcirc(cycperm(nwpanels,k,count1),i))
!      else if (wakeneigh(i,j) == -5) then
!        c(k*nwakepan+i,j) = 0.5*(wcirc(cycperm(nwpanels,k-1,count1),i) - wcirc(cycperm(nwpanels,k-2,count1),i))
!      else if (wakeneigh(i,j) == -6) then
!        c(k*nwakepan+i,j) = 0.5*(wcirc(cycperm(nwpanels,k-1,count1),i) - wcirc(cycperm(nwpanels,k,count1),i))
!      end if
!    end do
!  end do
!  end do

!  do i=1,nwakepan,1
!    call newvtxpart(firstpan(i,1,:),firstpan(i,2,:),firstpan(i,3,:),firstpan(i,4,:), &
!         & c(i,1),c(i,2),c(i,3),c(i,4),vpos(i,:),vtx(i,:))
!    call newvtxpart(lastpan(i,1,:),lastpan(i,2,:),lastpan(i,3,:),lastpan(i,4,:), &
!         & c(nwakepan*(nwpanels+1)+i,1),c(nwakepan*(nwpanels+1)+i,2),c(nwakepan*(nwpanels+1)+i,3),c(nwakepan*(nwpanels+1)+i,4), &
!         & vpos(nwakepan*(nwpanels+1)+i,:),vtx(nwakepan*(nwpanels+1)+i,:))
!    do j=1,nwpanels,1
!      call newvtxpart(wakepan(j,i,1,:),wakepan(j,i,2,:),wakepan(j,i,3,:),wakepan(j,i,4,:), &
!         & c(nwakepan*j+i,1),c(nwakepan*j+i,2),c(nwakepan*j+i,3),c(nwakepan*j+i,4),vpos(nwakepan*j+i,:),vtx(nwakepan*j+i,:))

!    end do
!  end do

!end subroutine vtxwakecollapse

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

!!Collapses a line of wake panels to a vortex particle.
!subroutine vtxcollapse2(nwakepan,wakepans1,wakepans2,wakeneigh,vpos,vtx)
!  use wakepanels
!!  use maths, only: cycperm
!  implicit none

!  integer,intent(in) :: nwakepan
!  type(wakepanel),dimension(:),intent(in) :: wakepans1,wakepans2
!  integer,dimension(nwakepan,4),intent(inout) :: wakeneigh
!  real*8,dimension(nwakepan,4) :: c
!  real*8,dimension(3,nwakepan),intent(out) :: vpos,vtx
!  integer :: i,j

!  c = 0.

!  do i=1,nwakepan,1
!    do j=1,4,1
!      if (wakeneigh(i,j) > 0) then
!        c(i,j) = 0.5*(wakepans1(i)%circ - wakepans1(wakeneigh(i,j))%circ)
!      else if (wakeneigh(i,j) == 0) then
!        c(i,j) = wakepans1(i)%circ
!      else if (wakeneigh(i,j) == -5) then
!        c(i,j) = wakepans1(i)%circ - wakepans2(i)%circ
!      else if (wakeneigh(i,j) == -6) then
!        c(i,j) = wakepans1(i)%circ
!        wakeneigh(i,j) = -2
!      end if
!    end do
!  end do

!  do i=1,nwakepan,1
!    call newvtxpart(wakepans1(i)%p(1)%pIN,wakepans1(i)%p(2)%pIN,wakepans1(i)%p(3)%pIN,wakepans1(i)%p(4)%pIN, &
!      & c(i,1),c(i,2),c(i,3),c(i,4),vpos(:,i),vtx(:,i))
!  end do

!end subroutine vtxcollapse2

!!Collapses a line of wake panels to a vortex particle.
!subroutine vtxcollapse3(t,nwakepan,wakepans1,wakepans2,wakeneigh,vpos,vtx)
!  use wakepanels
!!  use maths, only: cycperm
!  implicit none

!  integer,intent(in) :: t,nwakepan
!  type(wakepanel),dimension(:),intent(in) :: wakepans1,wakepans2
!  integer,dimension(nwakepan,4),intent(inout) :: wakeneigh
!  real*8,dimension(nwakepan,4) :: c
!  real*8,dimension(3,nwakepan),intent(out) :: vpos,vtx
!  integer :: i,j

!  c = 0.

!  do i=1,nwakepan,1
!    do j=1,4,1
!      if (wakeneigh(i,j) > 0) then
!        c(i,j) = 0.5*(wakepans1(i)%circ - wakepans1(wakeneigh(i,j))%circ)
!      else if (wakeneigh(i,j) == 0) then
!        c(i,j) = wakepans1(i)%circ
!      else if (wakeneigh(i,j) == -5) then
!        c(i,j) = wakepans1(i)%circ - wakepans2(i)%circ
!      else if (wakeneigh(i,j) == -6) then
!        c(i,j) = wakepans1(i)%circ
!        wakeneigh(i,j) = -2
!      end if
!    end do
!  end do

!  do i=1,nwakepan,1
!    call newvtxpart(wakepans1(i)%p(1)%pIN(:,t),wakepans1(i)%p(2)%pIN(:,t),wakepans1(i)%p(3)%pIN(:,t),wakepans1(i)%p(4)%pIN(:,t), &
!      & c(i,1),c(i,2),c(i,3),c(i,4),vpos(:,i),vtx(:,i))
!  end do

!end subroutine vtxcollapse3

!Collapses a line of wake panels to a vortex particle.
!subroutine vtxcollapse4(t,nwakepan,nwpanels,wakepans,vpos,vtx)
!  use wakepanels
!!  use maths, only: cycperm
!  implicit none

!  integer,intent(in) :: t,nwakepan,nwpanels
!  type(wakepanel),dimension(:),intent(inout) :: wakepans
!  real*8,dimension(nwakepan,4) :: c
!  real*8,dimension(3,nwakepan),intent(out) :: vpos,vtx
!  integer :: i,j,dummy1i1

!  c = 0.

!  do i=1,nwakepan,1
!    dummy1i1 = (nwpanels+1)*nwakepan + i
!    do j=1,4,1
!      if (wakepans(dummy1i1)%neigh(j) > 0) then
!        c(i,j) = 0.5*(wakepans(dummy1i1)%circ(t) - wakepans(wakepans(dummy1i1)%neigh(j))%circ(t))
!      else if (wakepans(dummy1i1)%neigh(j) == 0) then
!        c(i,j) = wakepans(dummy1i1)%circ(t)
!      else if (wakepans(dummy1i1)%neigh(j) == -5) then
!        c(i,j) = wakepans(dummy1i1)%circ(t)
!      else if (wakepans(dummy1i1)%neigh(j) == -6) then
!        c(i,j) = wakepans(dummy1i1)%circ(t)
!        wakepans(dummy1i1)%neigh(j) = -2
!      end if
!    end do
!  end do

!  do i=1,nwakepan,1
!    dummy1i1 = (nwpanels+1)*nwakepan + i
!    call newvtxpart(wakepans(dummy1i1)%p(1)%pIN(:,t),wakepans(dummy1i1)%p(2)%pIN(:,t),wakepans(dummy1i1)%p(3)%pIN(:,t), &
!      & wakepans(dummy1i1)%p(4)%pIN(:,t),c(i,1),c(i,2),c(i,3),c(i,4),vpos(:,i),vtx(:,i))
!  end do

!end subroutine vtxcollapse4

subroutine vtxsegcollapse(p1,p2,dcirc,vpos,vm)
  implicit none
  real*8,dimension(3),intent(in) :: p1,p2
  real*8,intent(in) :: dcirc
  real*8,dimension(3),intent(out) :: vpos,vm

  vpos = 0.5*(p1+p2)
  vm = dcirc * (p2-p1)

end subroutine vtxsegcollapse

subroutine vtxcollapse5(nwpanels,nwakepan,d1i1,tw,wakepans,tv,vtxpnts)
  use basictypes
!  use maths, only: cycperm
  implicit none

  integer,intent(in) :: nwpanels,nwakepan,d1i1,tw,tv
  type(panel),intent(inout) :: wakepans
  type(point),intent(inout) :: vtxpnts

  integer :: k
  real*8,dimension(4) :: c

  c = 0.
  do k=1,4,1
    if (wakepans%neigh(k) > 0.) then
      if (wakepans%neigh(k) > (nwpanels+1)*nwakepan) then
        c(k) = 0.5*wakepans%dcirc(k,tw)
      else
        c(k) = 0. !wakepans%dcirc(k,tw)
      end if
    else if (wakepans%neigh(k) == 0) then
      c(k) = wakepans%circ(tw)
    else if (wakepans%neigh(k) == -5) then
      c(k) = 0. !wakepans%dcirc(k,tw)
      print*, 'Error in vtxcollapse5'
    else if (wakepans%neigh(k) == -6) then
      c(k) = wakepans%circ(tw)
!print*, d1i1 - (nwpanels+1)*nwakepan
      wakepans%neigh(k) = -2
    else if (wakepans%neigh(k) == -2) then
      c(k) = wakepans%dcirc(k,tw)
    end if
  end do

  call newvtxpart(wakepans%pnt(1)%p%set(tw)%infr,wakepans%pnt(2)%p%set(tw)%infr,wakepans%pnt(3)%p%set(tw)%infr, &
    & wakepans%pnt(4)%p%set(tw)%infr,c(1),c(2),c(3),c(4),vtxpnts%set(tv)%infr,vtxpnts%set(tv)%minfr)

  vtxpnts%set(tv)%infr(4) = dble(tv)
  vtxpnts%set(tv)%act = .true.  

end subroutine vtxcollapse5

subroutine vtxcollapsewing(tw,wakepans,tv,vtxpnts)
  use basictypes
!  use maths, only: cycperm
  implicit none

  integer,intent(in) :: tw,tv
  type(panel),intent(inout) :: wakepans
  type(point),intent(inout) :: vtxpnts

  integer :: k
  real*8,dimension(4) :: c

  c = 0.
  do k=1,4,1
    if (wakepans%neigh(k) > 0.) then
      c(k) = 0.5*wakepans%dcirc(k,tw)
!      c(k) = wakepans%dcirc(k,tw)
    else if (wakepans%neigh(k) == 0) then
      c(k) = wakepans%dcirc(k,tw)
    else if (wakepans%neigh(k) == -5) then
      c(k) = wakepans%dcirc(k,tw)
    else if (wakepans%neigh(k) == -6) then
      c(k) = wakepans%dcirc(k,tw)
      wakepans%neigh(k) = -2
    end if
  end do

  call newvtxpart(wakepans%pnt(1)%p%set(tw)%infr,wakepans%pnt(2)%p%set(tw)%infr,wakepans%pnt(3)%p%set(tw)%infr, &
    & wakepans%pnt(4)%p%set(tw)%infr,c(1),c(2),c(3),c(4),vtxpnts%set(tv)%infr,vtxpnts%set(tv)%minfr)

  vtxpnts%set(tv)%infr(4) = dble(tv)
  vtxpnts%set(tv)%act = .true.  

end subroutine vtxcollapsewing
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

!!Collapses a line of wake panels to a vortex particle.
!subroutine vtxcollapse(nwakepan,lastpanIN,wakecirc,wakecirctm1,wakeneigh,vpos,vtx)
!!  use maths, only: cycperm
!  implicit none

!  integer,intent(in) :: nwakepan
!  real*8,dimension(nwakepan,4,3),intent(in) :: lastpanIN
!  real*8,dimension(nwakepan),intent(in) :: wakecirc,wakecirctm1
!  integer,dimension(nwakepan,4),intent(inout) :: wakeneigh
!  real*8,dimension(nwakepan,4,3) :: wpancirc
!  real*8,dimension(nwakepan,4) :: c
!  real*8,dimension(nwakepan,3),intent(out) :: vpos,vtx
!  integer :: i,j

!  wpancirc = 0.
!  c = 0.

!  do i=1,nwakepan,1
!    do j=1,4,1
!      if (wakeneigh(i,j) > 0) then
!        c(i,j) = 0.5*(wakecirc(i) - wakecirc(wakeneigh(i,j)))
!      else if (wakeneigh(i,j) == 0) then
!        c(i,j) = wakecirc(i)
!      else if (wakeneigh(i,j) == -5) then
!        c(i,j) = wakecirc(i) - wakecirctm1(i)
!      else if (wakeneigh(i,j) == -6) then
!        c(i,j) = wakecirc(i)
!        wakeneigh(i,j) = -2
!      end if
!    end do
!  end do

!  do i=1,nwakepan,1

!call newvtxpart(lastpanIN(i,1,:),lastpanIN(i,2,:),lastpanIN(i,3,:),lastpanIN(i,4,:),c(i,1),c(i,2),c(i,3),c(i,4),vpos(i,:),vtx(i,:))

!  end do

!end subroutine vtxcollapse

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

!logical function incheck6(p,grid,dbox,mingrid)
!  implicit none

!  real*8,dimension(3),intent(in) :: p,dbox,mingrid
!  integer,dimension(3),intent(in) :: grid

!  if (mingrid(1)+(grid(1)-1)*dbox(1) < p(1)) then
!    if (mingrid(1)+grid(1)*dbox(1) > p(1)) then
!      if (mingrid(2)+(grid(2)-1)*dbox(2) < p(2)) then
!        if (mingrid(2)+grid(2)*dbox(2) > p(2)) then
!          if (mingrid(3)+(grid(3)-1)*dbox(3) < p(3)) then
!            if (mingrid(3)+grid(3)*dbox(3) > p(3)) then
!              incheck6 = .true.
!            else
!              incheck6 = .false.
!            end if
!          else
!            incheck6 = .false.
!          end if
!        else
!          incheck6 = .false.
!        end if
!      else
!        incheck6 = .false.
!      end if
!    else
!      incheck6 = .false.
!    end if
!  else
!    incheck6 = .false.
!  end if

!end function incheck6

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

subroutine shearind(pos,vtxpos,vtxm,dummy33r)
  use maths, only: cross,pi
  implicit none

  real*8,dimension(3,3),intent(out) :: dummy33r
  real*8,dimension(3),intent(in) :: pos,vtxpos,vtxm
  real*8,dimension(3) :: partialx,partialy,partialz,r
  real*8 :: r1sq,r2sq,r3sq !,r5,rmod
  real*8 :: sigma,sigmat3,sigmat5,rho,rhosq,r2p1,r2p1t12,r2p1t32,r2p1t52,c1,c2

  r = pos-vtxpos

  if (dot_product(r,r) > 0.1) then

  r1sq = r(1)*r(1)
  r2sq = r(2)*r(2)
  r3sq = r(3)*r(3)

  sigma = 0.4
  sigmat3 = sigma*sigma*sigma
  sigmat5 = sigmat3*sigma*sigma
  rho = sqrt(r1sq+r2sq+r3sq)/sigma
  rhosq = rho*rho
  r2p1 = rhosq+1.
  r2p1t12 = sqrt(r2p1)
  r2p1t32 = r2p1*r2p1t12
  r2p1t52 = r2p1t32*r2p1

  c1 = (rhosq+2.5)/(4.*pi*sigmat3*r2p1t52)
  c2 = (5.*(rhosq+1.5)*r2p1t52 - (rhosq+2.5)*(3.*r2p1t52 + 5.*rhosq*r2p1t32))/(4.*pi*sigmat5*rhosq*r2p1t52*r2p1t52)

  partialx(1) = r1sq*c2 + c1
  partialx(2) = r(1)*r(2)*c2
  partialx(3) = r(1)*r(3)*c2

  partialy(1) = partialx(2)
  partialy(2) = r2sq*c2 + c1
  partialy(3) = r(2)*r(3)*c2

  partialz(1) = partialx(3)
  partialz(2) = partialy(3)
  partialz(3) = r3sq*c2 + c1

  dummy33r(:,1) = cross(vtxm,partialx)
  dummy33r(:,2) = cross(vtxm,partialy)
  dummy33r(:,3) = cross(vtxm,partialz)

!  r1sq = r(1)*r(1)
!  r2sq = r(2)*r(2)
!  r3sq = r(3)*r(3)
!  rmod = sqrt(r1sq+r2sq+r3sq)

!  r5 = 4.*pi*rmod*rmod*rmod*rmod*rmod

!  partialx(1) = (r2sq + r3sq - 2*r1sq) / r5
!  partialx(2) = (-3.*r1sq * r2sq) / r5
!  partialx(3) = (-3.*r1sq * r3sq) / r5

!  partialy(1) = partialx(2)
!  partialy(2) = (r1sq + r3sq - 2*r2sq) / r5
!  partialy(3) = (-3.*r2sq * r3sq) / r5

!  partialz(1) = partialx(3)
!  partialz(2) = partialy(3)
!  partialz(3) = (r1sq + r2sq - 2*r3sq) / r5

!  dummy33r(:,1) = cross(vtxm,partialx)
!  dummy33r(:,2) = cross(vtxm,partialy)
!  dummy33r(:,3) = cross(vtxm,partialz)

  else
    dummy33r = 0.
  end if

end subroutine shearind

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

!function gvtxpos(gridmin,grid,dbox,lev,n)
!  implicit none

!  real*8,dimension(3) :: gvtxpos
!  real*8,dimension(3),intent(in) :: gridmin,dbox
!  integer,dimension(3),intent(in) :: grid
!  integer,intent(in) :: lev,n
!  integer :: i

!  real*8,dimension(3) :: dboxl,rpos

!  dboxl = dbox*2**(lev)

!  do i=1,3,1
!    rpos(i) = gridmin(i) + (grid(i)-1)*dboxl(i)
!  end do
!  
!  gvtxpos = rpos8(rpos,dboxl,n)

!end function gvtxpos

!!-----------------------------------------------------------------------------------------------------------------------------------
!!-----------------------------------------------------------------------------------------------------------------------------------

!function rpos8(rpos,dbox,n)
!  implicit none

!  real*8,dimension(3) :: rpos8
!  real*8,dimension(3),intent(in) :: rpos,dbox
!  integer,intent(in) :: n

!  if (n==1) then
!  rpos8 = 0.
!  else if (n==2) then
!    rpos8(1) = dbox(1)
!    rpos8(2) = 0.
!    rpos8(3) = 0.
!  else if (n==3) then
!    rpos8(1) = 0.
!    rpos8(2) = dbox(2)
!    rpos8(3) = 0.
!  else if (n==4) then
!    rpos8(1) = dbox(1)
!    rpos8(2) = dbox(2)
!    rpos8(3) = 0.
!  else if (n==5) then
!    rpos8(1) = 0.
!    rpos8(2) = 0.
!    rpos8(3) = dbox(3)
!  else if (n==6) then
!    rpos8(1) = dbox(1)
!    rpos8(2) = 0.
!    rpos8(3) = dbox(3)
!  else if (n==7) then
!    rpos8(1) = 0. 
!    rpos8(2) = dbox(2)
!    rpos8(3) = dbox(3)
!  else if (n==8) then
!    rpos8 = dbox
!  end if

!  rpos8 = rpos8 + rpos

!end function rpos8

!!-----------------------------------------------------------------------------------------------------------------------------------
!!-----------------------------------------------------------------------------------------------------------------------------------

end module vtxparticle

!***********************************************************************************************************************************
!***********************************************************************************************************************************

!module gridtree
!  implicit none

!  contains

!!-----------------------------------------------------------------------------------------------------------------------------------
!!-----------------------------------------------------------------------------------------------------------------------------------

!  subroutine treestruct(ngrid,nl,trees)
!    use gtmodinit
!    implicit none

!    integer,dimension(3),intent(in) :: ngrid
!    integer,intent(out) :: nl
!    type (level),dimension(:),allocatable :: trees

!    integer,dimension(1) :: mdim
!    integer,dimension(3) :: dummy3i1

!    integer :: i,j,k,l,m,n,count1,count2,count3

!  mdim = maxloc(ngrid)
!  i = int((ngrid(mdim(1))+1)*0.5)
!  nl = 1
!  k = i
!  
!  if (k>1) then
!    do
!      nl = nl + 1
!      k = int((k+1)*0.5)
!      if (k<=1) exit
!    end do
!  end if

!  allocate(trees(nl))

!  do i=1,3,1
!    trees(1)%nlevelem(i) = int((ngrid(i)+1)*0.5)
!  end do

!  do i=2,nl,1
!    do j=1,3,1
!      trees(i)%nlevelem(j) = int((trees(i-1)%nlevelem(j)+1)*0.5)
!    end do
!  end do

!  do i=1,nl,1
!    allocate(trees(i)%levelem(trees(i)%nlevelem(1),trees(i)%nlevelem(2),trees(i)%nlevelem(3)))
!  end do

!  count1 = 0
!  count2 = 0
!  count3 = 0

!  dummy3i1 = ngrid

!  do l=1,nl,1
!    count1 = 0
!    count2 = 0
!    count3 = 0
!    do i=1,dummy3i1(1),2
!      count1 = int((i+1)*0.5)
!      do j=1,dummy3i1(2),2
!        count2 = int((j+1)*0.5)
!        do k=1,dummy3i1(3),2
!          count3 = int((k+1)*0.5)

!        trees(l)%levelem(count1,count2,count3)%parent(1) = int((count1+1)*0.5)
!        trees(l)%levelem(count1,count2,count3)%parent(2) = int((count2+1)*0.5)
!        trees(l)%levelem(count1,count2,count3)%parent(3) = int((count3+1)*0.5)

!        trees(l)%levelem(count1,count2,count3)%grid(1,1) = i
!        trees(l)%levelem(count1,count2,count3)%grid(1,2) = j
!        trees(l)%levelem(count1,count2,count3)%grid(1,3) = k
!        trees(l)%levelem(count1,count2,count3)%grid(2,1) = i+1
!        trees(l)%levelem(count1,count2,count3)%grid(2,2) = j
!        trees(l)%levelem(count1,count2,count3)%grid(2,3) = k
!        trees(l)%levelem(count1,count2,count3)%grid(3,1) = i
!        trees(l)%levelem(count1,count2,count3)%grid(3,2) = j+1
!        trees(l)%levelem(count1,count2,count3)%grid(3,3) = k
!        trees(l)%levelem(count1,count2,count3)%grid(4,1) = i+1
!        trees(l)%levelem(count1,count2,count3)%grid(4,2) = j+1
!        trees(l)%levelem(count1,count2,count3)%grid(4,3) = k
!        trees(l)%levelem(count1,count2,count3)%grid(5,1) = i
!        trees(l)%levelem(count1,count2,count3)%grid(5,2) = j
!        trees(l)%levelem(count1,count2,count3)%grid(5,3) = k+1
!        trees(l)%levelem(count1,count2,count3)%grid(6,1) = i+1
!        trees(l)%levelem(count1,count2,count3)%grid(6,2) = j
!        trees(l)%levelem(count1,count2,count3)%grid(6,3) = k+1
!        trees(l)%levelem(count1,count2,count3)%grid(7,1) = i
!        trees(l)%levelem(count1,count2,count3)%grid(7,2) = j+1
!        trees(l)%levelem(count1,count2,count3)%grid(7,3) = k+1
!        trees(l)%levelem(count1,count2,count3)%grid(8,1) = i+1
!        trees(l)%levelem(count1,count2,count3)%grid(8,2) = j+1
!        trees(l)%levelem(count1,count2,count3)%grid(8,3) = k+1
!        end do
!      end do
!    end do
!    do i=1,trees(l)%nlevelem(1),1
!      do j=1,trees(l)%nlevelem(2),1
!        do k=1,trees(l)%nlevelem(3),1
!          do m=1,8,1
!            nloop: do n=1,3,1
!              if (trees(l)%levelem(i,j,k)%grid(m,n) > dummy3i1(n)) then
!                trees(l)%levelem(i,j,k)%grid(m,:) = 0
!                exit nloop
!              end if
!            end do nloop
!          end do
!        end do
!      end do
!    end do
!    dummy3i1 = trees(l)%nlevelem
!  end do

!  end subroutine treestruct

!!-----------------------------------------------------------------------------------------------------------------------------------
!!-----------------------------------------------------------------------------------------------------------------------------------

!  !Calculate influence coeffs on vtx box corner at all levels
!  subroutine vtxcorner(dbox,gridext,grids,grids_active,nl,trees,vtx)
!    use gtmodinit
!    use vtxparticle
!    implicit none

!    type (vtxpart),dimension(:),intent(in) :: vtx
!    type (level),dimension(:) :: trees
!    type (gridstruct),dimension(:,:,:),intent(in) :: grids
!    logical,dimension(:,:,:),intent(in) :: grids_active

!    integer,intent(in) :: nl
!    real*8,dimension(3,2),intent(in) :: gridext
!    real*8,dimension(3),intent(in) :: dbox

!    integer :: j,k,l,m,n,o
!    integer :: dummy1i1
!    integer,dimension(3) :: dummy3i1
!    real*8,dimension(3) :: dummy3r
!    real*8,dimension(3),parameter :: zerovec = (/0.,0.,0./)

!    do m=1,trees(1)%nlevelem(3),1
!    do l=1,trees(1)%nlevelem(2),1
!    do k=1,trees(1)%nlevelem(1),1
!      if (trees(1)%levelem(k,l,m)%active) then
!        do n=1,8,1
!          dummy3i1 = trees(1)%levelem(k,l,m)%grid(n,:)
!          if (dummy3i1(1) > 0) then
!          if (grids_active(dummy3i1(1),dummy3i1(2),dummy3i1(3))) then
!            do o=1,grids(dummy3i1(1),dummy3i1(2),dummy3i1(3))%nvtx,1
!              dummy1i1 = grids(dummy3i1(1),dummy3i1(2),dummy3i1(3))%vtx(o)

!              dummy3r(1) = vtx(dummy1i1)%posIN(1) - (gridext(1,1) + (k-1)*(dbox(1)*2))
!              dummy3r(2) = vtx(dummy1i1)%posIN(2) - (gridext(2,1) + (l-1)*(dbox(2)*2))
!              dummy3r(3) = vtx(dummy1i1)%posIN(3) - (gridext(3,1) + (m-1)*(dbox(3)*2))

!              call vtxinterp(dummy3r,vtx(dummy1i1)%mIN,dbox*2,trees(1)%levelem(k,l,m)%gridvtx(:,:))

!            end do
!          end if
!          end if
!        end do
!      end if
!    end do
!    end do
!    end do

!    do j=2,nl,1
!    do m=1,trees(j)%nlevelem(3),1
!    do l=1,trees(j)%nlevelem(2),1
!    do k=1,trees(j)%nlevelem(1),1
!      if (trees(j)%levelem(k,l,m)%active) then
!        do n=1,8,1
!          dummy3i1 = trees(j)%levelem(k,l,m)%grid(n,:)
!          if (dummy3i1(1) > 0) then
!          if (trees(j-1)%levelem(dummy3i1(1),dummy3i1(2),dummy3i1(3))%active) then
!            do o=1,8,1
!              dummy3r = rpos8(rpos8(zerovec,dbox*2**(j-1),n),dbox*2**(j-1),o)

!              call vtxinterp(dummy3r,trees(j-1)%levelem(dummy3i1(1),dummy3i1(2),dummy3i1(3))%gridvtx(o,:),dbox*2**j, &
!                 & trees(j)%levelem(k,l,m)%gridvtx(:,:))

!            end do
!          end if
!          end if
!        end do
!      end if
!    end do
!    end do
!    end do
!    end do

!  end subroutine vtxcorner

!!-----------------------------------------------------------------------------------------------------------------------------------
!!-----------------------------------------------------------------------------------------------------------------------------------

!  !Calculate influence coeffs on vtx box corner at all levels
!  subroutine vtxcorner2(dbox,fvtx,gridext,grids,nl,trees,vtx)
!    use gtmodinit
!    use vtxparticle
!    implicit none

!    type (vtxpart),dimension(:),intent(in) :: vtx
!    type (level),dimension(:) :: trees
!    type (gridstruct),dimension(:,:,:),intent(in) :: grids

!    integer,intent(in) :: nl,fvtx
!    real*8,dimension(3,2),intent(in) :: gridext
!    real*8,dimension(3),intent(in) :: dbox

!    integer :: j,k,l,m,n,o
!    integer :: dummy1i1
!    integer,dimension(3) :: dummy3i1
!    real*8,dimension(3) :: dummy3r
!    real*8,dimension(3),parameter :: zerovec = (/0.,0.,0./)

!    do k=1,trees(1)%nlevelem(1),1
!    do l=1,trees(1)%nlevelem(2),1
!    do m=1,trees(1)%nlevelem(3),1
!      if (trees(1)%levelem(k,l,m)%ww) then
!!      if (trees(1)%levelem(k,l,m)%active) then
!        do n=1,8,1
!          dummy3i1 = trees(1)%levelem(k,l,m)%grid(n,:)
!          if (dummy3i1(1) > 0) then
!            if (grids(dummy3i1(1),dummy3i1(2),dummy3i1(3))%ww) then
!!            if (grids(dummy3i1(1),dummy3i1(2),dummy3i1(3))%active) then
!              do o=1,grids(dummy3i1(1),dummy3i1(2),dummy3i1(3))%nvtx,1
!                if (grids(dummy3i1(1),dummy3i1(2),dummy3i1(3))%vtx(o) <= fvtx) then
!                dummy1i1 = grids(dummy3i1(1),dummy3i1(2),dummy3i1(3))%vtx(o)
!  
!                dummy3r(1) = vtx(dummy1i1)%posIN(1) - (gridext(1,1) + (k-1)*(dbox(1)*2))
!                dummy3r(2) = vtx(dummy1i1)%posIN(2) - (gridext(2,1) + (l-1)*(dbox(2)*2))
!                dummy3r(3) = vtx(dummy1i1)%posIN(3) - (gridext(3,1) + (m-1)*(dbox(3)*2))
!  
!                call vtxinterp(dummy3r,vtx(dummy1i1)%mIN,dbox*2,trees(1)%levelem(k,l,m)%gridvtx(:,:))
!                end if
!              end do
!            end if
!          end if
!        end do
!      end if
!    end do
!    end do
!    end do

!    do j=2,nl,1
!    do k=1,trees(j)%nlevelem(1),1
!    do l=1,trees(j)%nlevelem(2),1
!    do m=1,trees(j)%nlevelem(3),1
!      if (trees(j)%levelem(k,l,m)%ww) then
!!      if (trees(j)%levelem(k,l,m)%active) then
!        trees(j)%levelem(k,l,m)%gridvtx = 0.
!        do n=1,8,1
!          dummy3i1 = trees(j)%levelem(k,l,m)%grid(n,:)
!          if (dummy3i1(1) > 0) then
!          if (trees(j-1)%levelem(dummy3i1(1),dummy3i1(2),dummy3i1(3))%active) then
!            do o=1,8,1
!              dummy3r = rpos8(rpos8(zerovec,dbox*2**(j-1),n),dbox*2**(j-1),o)

!              call vtxinterp(dummy3r,trees(j-1)%levelem(dummy3i1(1),dummy3i1(2),dummy3i1(3))%gridvtx(o,:),dbox*2**j, &
!                 & trees(j)%levelem(k,l,m)%gridvtx(:,:))

!            end do
!          end if
!          end if
!        end do
!      end if
!    end do
!    end do
!    end do
!    end do

!  end subroutine vtxcorner2

!!-----------------------------------------------------------------------------------------------------------------------------------
!!-----------------------------------------------------------------------------------------------------------------------------------

!  !Calculate influence coeffs on vtx box corner at all levels
!  subroutine vtxcorner3(dbox,gridext,grids,grids_active,nl,trees,vtx)
!    use gtmodinit
!    use vtxparticle
!    implicit none

!    type (vtxpart),dimension(:),intent(in) :: vtx
!    type (level),dimension(:) :: trees
!    type (gridstruct),dimension(:,:,:),intent(in) :: grids
!    logical,dimension(:,:,:),intent(in) :: grids_active

!    integer,intent(in) :: nl
!    real*8,dimension(3,2),intent(in) :: gridext
!    real*8,dimension(3),intent(in) :: dbox

!    integer :: j,k,l,m,n,o
!    integer :: dummy1i1
!    integer,dimension(3) :: dummy3i1
!    real*8,dimension(3) :: dummy3r
!    real*8,dimension(3),parameter :: zerovec = (/0.,0.,0./)

!    do m=1,trees(1)%nlevelem(3),1
!    do l=1,trees(1)%nlevelem(2),1
!    do k=1,trees(1)%nlevelem(1),1
!      if (trees(1)%levelem(k,l,m)%ww) then
!!      if (trees(1)%levelem(k,l,m)%active) then
!        trees(1)%levelem(k,l,m)%gridvtx = 0.
!        do n=1,8,1
!          dummy3i1 = trees(1)%levelem(k,l,m)%grid(n,:)
!          if (dummy3i1(1) > 0) then
!!            if (grids(dummy3i1(1),dummy3i1(2),dummy3i1(3))%ww) then
!            if (grids_active(dummy3i1(1),dummy3i1(2),dummy3i1(3))) then
!              do o=1,grids(dummy3i1(1),dummy3i1(2),dummy3i1(3))%nvtx,1
!                dummy1i1 = grids(dummy3i1(1),dummy3i1(2),dummy3i1(3))%vtx(o)
!  
!                dummy3r(1) = vtx(dummy1i1)%posIN(1) - (gridext(1,1) + (k-1)*(dbox(1)*2))
!                dummy3r(2) = vtx(dummy1i1)%posIN(2) - (gridext(2,1) + (l-1)*(dbox(2)*2))
!                dummy3r(3) = vtx(dummy1i1)%posIN(3) - (gridext(3,1) + (m-1)*(dbox(3)*2))
!  
!                call vtxinterp(dummy3r,vtx(dummy1i1)%mIN,dbox*2,trees(1)%levelem(k,l,m)%gridvtx(:,:))
!              end do
!            end if
!          end if
!        end do
!      end if
!    end do
!    end do
!    end do

!    do j=2,nl,1
!    do m=1,trees(j)%nlevelem(3),1
!    do l=1,trees(j)%nlevelem(2),1
!    do k=1,trees(j)%nlevelem(1),1
!      if (trees(j)%levelem(k,l,m)%ww) then
!!      if (trees(j)%levelem(k,l,m)%active) then
!        trees(j)%levelem(k,l,m)%gridvtx = 0.
!        do n=1,8,1
!          dummy3i1 = trees(j)%levelem(k,l,m)%grid(n,:)
!          if (dummy3i1(1) > 0) then
!          if (trees(j-1)%levelem(dummy3i1(1),dummy3i1(2),dummy3i1(3))%active) then
!            do o=1,8,1
!              dummy3r = rpos8(rpos8(zerovec,dbox*2**(j-1),n),dbox*2**(j-1),o)

!              call vtxinterp(dummy3r,trees(j-1)%levelem(dummy3i1(1),dummy3i1(2),dummy3i1(3))%gridvtx(o,:),dbox*2**j, &
!                 & trees(j)%levelem(k,l,m)%gridvtx(:,:))

!            end do
!          end if
!          end if
!        end do
!      end if
!    end do
!    end do
!    end do
!    end do

!  end subroutine vtxcorner3

!!-----------------------------------------------------------------------------------------------------------------------------------
!!-----------------------------------------------------------------------------------------------------------------------------------
!!-----------------------------------------------------------------------------------------------------------------------------------
!!-----------------------------------------------------------------------------------------------------------------------------------

!  subroutine gridboxcalc(j,k,l,m,n,gridext,dbox,vtx,trees,vtxperim)
!    use gtmodinit
!    use vtxparticle
!    implicit none

!    type (vtxpart),dimension(:),intent(inout) :: vtx
!    type (level),dimension(:),intent(in) :: trees

!    integer,intent(in) :: j,k,l,m,n
!    real*8,dimension(3,2),intent(in) :: gridext
!    real*8,dimension(3),intent(in) :: dbox
!    integer,dimension(:,:,:),intent(in) :: vtxperim

!    real*8,dimension(3) :: dummy3r
!    real*8,dimension(3,3) :: dummy33r
!    integer,dimension(3) :: dummy3i2
!    integer :: o,p

!    if (trees(k)%levelem(l,m,n)%active) then
!      do o=1,8,1
!        dummy3i2 = trees(k)%levelem(l,m,n)%grid(o,:)
!        if (dummy3i2(1) > 0) then
!        if (trees(k-1)%levelem(dummy3i2(1),dummy3i2(2),dummy3i2(3))%active) then
!          if (dummy3i2(1) < vtxperim(k-1,1,1) .or. dummy3i2(2) < vtxperim(k-1,1,2) .or. dummy3i2(3) < vtxperim(k-1,1,3) .or. &
!            & dummy3i2(1) > vtxperim(k-1,2,1) .or. dummy3i2(2) > vtxperim(k-1,2,2) .or. dummy3i2(3) > vtxperim(k-1,2,3)) then
!            do p=1,8,1
!              call vtxinduced(vtx(j)%posIN,gvtxpos(gridext(:,1),dummy3i2,dbox,k-1,p), &
!                 & trees(k-1)%levelem(dummy3i2(1),dummy3i2(2),dummy3i2(3))%gridvtx(p,:),dummy3r)
!              vtx(j)%vel = vtx(j)%vel + dummy3r
!              call shearind(vtx(j)%posIN,gvtxpos(gridext(:,1),dummy3i2,dbox,k-1,p), &
!                             & trees(k-1)%levelem(dummy3i2(1),dummy3i2(2),dummy3i2(3))%gridvtx(p,:),dummy33r)
!              vtx(j)%vtensor = vtx(j)%vtensor + dummy33r
!            end do
!          end if
!        end if
!        end if
!      end do
!    end if

!  end subroutine gridboxcalc

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

!  subroutine gridboxcalc2(pos,k,l,m,n,gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!    use gtmodinit
!    use vtxparticle
!    use pg
!    implicit none

!    type (level),dimension(:),intent(in) :: trees
!    integer,intent(in) :: k,l,m,n
!    logical,intent(in) :: mscale
!    real*8,intent(in) :: a
!    real*8,dimension(3),intent(in) :: uinf
!    real*8,dimension(3,2),intent(in) :: gridext
!    real*8,dimension(3),intent(in) :: pos,dbox
!    real*8,dimension(3,3),intent(inout) :: vtens
!    real*8,dimension(3),intent(inout) :: vel
!    integer,dimension(:,:,:),intent(in) :: vtxperim

!    real*8,dimension(3) :: dummy3r
!    real*8,dimension(3,3) :: dummy33r
!    integer,dimension(3) :: dummy3i2
!    integer :: o,p

!    if (trees(k)%levelem(l,m,n)%active) then
!      do o=1,8,1
!        dummy3i2 = trees(k)%levelem(l,m,n)%grid(o,:)
!        if (dummy3i2(1) > 0) then
!        if (trees(k-1)%levelem(dummy3i2(1),dummy3i2(2),dummy3i2(3))%active) then
!          if (dummy3i2(1) < vtxperim(k-1,1,1) .or. dummy3i2(2) < vtxperim(k-1,1,2) .or. dummy3i2(3) < vtxperim(k-1,1,3) .or. &
!            & dummy3i2(1) > vtxperim(k-1,2,1) .or. dummy3i2(2) > vtxperim(k-1,2,2) .or. dummy3i2(3) > vtxperim(k-1,2,3)) then
!            do p=1,8,1
!             if (mscale) then

!                call vtxinduced(pos,rmscale(a,uinf,pos,gvtxpos(gridext(:,1),dummy3i2,dbox,k-1,p)), &
!                  & trees(k-1)%levelem(dummy3i2(1),dummy3i2(2),dummy3i2(3))%gridvtx(p,:),dummy3r)
!                call shearind(pos,rmscale(a,uinf,pos,gvtxpos(gridext(:,1),dummy3i2,dbox,k-1,p)), &
!                  & trees(k-1)%levelem(dummy3i2(1),dummy3i2(2),dummy3i2(3))%gridvtx(p,:),dummy33r)

!              else

!                call vtxinduced(pos,gvtxpos(gridext(:,1),dummy3i2,dbox,k-1,p), &
!                  & trees(k-1)%levelem(dummy3i2(1),dummy3i2(2),dummy3i2(3))%gridvtx(p,:),dummy3r)
!                call shearind(pos,gvtxpos(gridext(:,1),dummy3i2,dbox,k-1,p), &
!                  & trees(k-1)%levelem(dummy3i2(1),dummy3i2(2),dummy3i2(3))%gridvtx(p,:),dummy33r)

!              end if
!              vel = vel + dummy3r
!              vtens = vtens + dummy33r
!            end do
!          end if
!        end if
!        end if
!      end do
!    end if

!  end subroutine gridboxcalc2

!!-----------------------------------------------------------------------------------------------------------------------------------
!!-----------------------------------------------------------------------------------------------------------------------------------

!  subroutine gridboxcalcdipole(pos,k,l,m,n,gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!    use gtmodinit
!    use vtxparticle
!    use pg
!    implicit none

!    type (level),dimension(:),intent(in) :: trees
!    integer,intent(in) :: k,l,m,n
!    logical,intent(in) :: mscale
!    real*8,intent(in) :: a
!    real*8,dimension(3),intent(in) :: uinf
!    real*8,dimension(3,2),intent(in) :: gridext
!    real*8,dimension(3),intent(in) :: pos,dbox
!    real*8,dimension(3,3),intent(inout) :: vtens
!    real*8,dimension(3),intent(inout) :: vel
!    integer,dimension(:,:,:),intent(in) :: vtxperim

!    real*8,dimension(3) :: dummy3r
!    real*8,dimension(3,3) :: dummy33r
!    integer,dimension(3) :: dummy3i2
!    integer :: o,p

!    if (trees(k)%levelem(l,m,n)%active) then
!      do o=1,8,1
!        dummy3i2 = trees(k)%levelem(l,m,n)%grid(o,:)
!        if (dummy3i2(1) > 0) then
!        if (trees(k-1)%levelem(dummy3i2(1),dummy3i2(2),dummy3i2(3))%active) then
!          if (dummy3i2(1) < vtxperim(k-1,1,1) .or. dummy3i2(2) < vtxperim(k-1,1,2) .or. dummy3i2(3) < vtxperim(k-1,1,3) .or. &
!            & dummy3i2(1) > vtxperim(k-1,2,1) .or. dummy3i2(2) > vtxperim(k-1,2,2) .or. dummy3i2(3) > vtxperim(k-1,2,3)) then
!            do p=1,8,1
!              if (mscale) then
!                call dipoleind(pos,rmscale(a,uinf,pos,gvtxpos(gridext(:,1),dummy3i2,dbox,k-1,p)), &
!                  & trees(k-1)%levelem(dummy3i2(1),dummy3i2(2),dummy3i2(3))%gridvtx(p,:),dummy3r)
!                call dipgradvind(pos,rmscale(a,uinf,pos,gvtxpos(gridext(:,1),dummy3i2,dbox,k-1,p)), &
!                  & trees(k-1)%levelem(dummy3i2(1),dummy3i2(2),dummy3i2(3))%gridvtx(p,:),dummy33r)
!              else
!                call dipoleind(pos,gvtxpos(gridext(:,1),dummy3i2,dbox,k-1,p), &
!                  & trees(k-1)%levelem(dummy3i2(1),dummy3i2(2),dummy3i2(3))%gridvtx(p,:),dummy3r)
!                call dipgradvind(pos,gvtxpos(gridext(:,1),dummy3i2,dbox,k-1,p), &
!                  & trees(k-1)%levelem(dummy3i2(1),dummy3i2(2),dummy3i2(3))%gridvtx(p,:),dummy33r)
!              end if
!              vel = vel + dummy3r
!              vtens = vtens + dummy33r
!            end do
!          end if
!        end if
!        end if
!      end do
!    end if

!  end subroutine gridboxcalcdipole


!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

!  subroutine gridboxcalc3(pos,k,l,m,n,gridext,dbox,trees,vtxperim,vel,vtens)
!    use gtmodinit
!    use vtxparticle
!    implicit none

!    type (level),dimension(:),intent(in) :: trees

!    integer,intent(in) :: k,l,m,n
!    real*8,dimension(2,3),intent(in) :: gridext
!    real*8,dimension(3),intent(in) :: pos,dbox
!!    real*8,dimension(3,3),intent(inout) :: vtens
!    real*8,dimension(3),intent(inout) :: vel
!    integer,dimension(:,:,:),intent(in) :: vtxperim

!    real*8,dimension(3) :: dummy3r
!    real*8,dimension(3,3) :: dummy33r
!    integer,dimension(3) :: dummy3i2
!    integer :: o,p

!    if (trees(k)%levelem(l,m,n)%active) then
!      do o=1,8,1
!        dummy3i2 = trees(k)%levelem(l,m,n)%grid(o,:)
!        if (dummy3i2(1) > 0) then
!        if (trees(k-1)%levelem(dummy3i2(1),dummy3i2(2),dummy3i2(3))%active) then
!          if (dummy3i2(1) < vtxperim(k-1,1,1) .or. dummy3i2(2) < vtxperim(k-1,1,2) .or. dummy3i2(3) < vtxperim(k-1,1,3) .or. &
!            & dummy3i2(1) > vtxperim(k-1,2,1) .or. dummy3i2(2) > vtxperim(k-1,2,2) .or. dummy3i2(3) > vtxperim(k-1,2,3)) then
!            do p=1,8,1
!              call vtxinduced(pos,gvtxpos(gridext(1,:),dummy3i2,dbox,k-1,p), &
!                 & trees(k-1)%levelem(dummy3i2(1),dummy3i2(2),dummy3i2(3))%gridvtx(p,:),dummy3r)
!              vel = vel + dummy3r
!!              call shearind(pos,gvtxpos(gridext(1,:),dummy3i2,dbox,k-1,p), &
!!                             & trees(k-1)%levelem(dummy3i2(1),dummy3i2(2),dummy3i2(3))%gridvtx(p,:),dummy33r)
!!              vtens = vtens + dummy33r
!            end do
!          end if
!        end if
!        end if
!      end do
!    end if

!  end subroutine gridboxcalc3

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

!  subroutine ffieldlevel(j,k,gridext,dbox,vtx,trees,vtxperim,tbuffer)
!    use gtmodinit
!    use vtxparticle
!    implicit none

!    type (vtxpart),dimension(:),intent(inout) :: vtx
!    type (level),dimension(:),intent(in) :: trees

!    integer,intent(in) :: j,k
!    real*8,dimension(3,2),intent(in) :: gridext
!    real*8,dimension(3),intent(in) :: dbox
!    integer,dimension(:,:,:),intent(in) :: vtxperim
!    integer,dimension(:),intent(in) :: tbuffer

!    integer :: l,m,n

!    if (vtxperim(k,2,1)-vtxperim(k,1,1) > 2*(tbuffer(k)+1)) then
!      do l=0,tbuffer(k),1
!        do m=vtxperim(k,1,2),vtxperim(k,2,2),1
!          do n=vtxperim(k,1,3),vtxperim(k,2,3),1
!            call gridboxcalc(j,k,vtxperim(k,1,1)+l,m,n,gridext,dbox,vtx,trees,vtxperim)
!            call gridboxcalc(j,k,vtxperim(k,2,1)-l,m,n,gridext,dbox,vtx,trees,vtxperim)
!          end do
!        end do
!      end do
!    else if (vtxperim(k,2,1)-vtxperim(k,1,1) > 0) then
!      do l=vtxperim(k,1,1),vtxperim(k,2,1),1
!        do m=vtxperim(k,1,2),vtxperim(k,2,2),1
!          do n=vtxperim(k,1,3),vtxperim(k,2,3),1
!            call gridboxcalc(j,k,l,m,n,gridext,dbox,vtx,trees,vtxperim)
!          end do
!        end do
!      end do
!    else
!      l = 1
!        do m=vtxperim(k,1,2),vtxperim(k,2,2),1
!          do n=vtxperim(k,1,3),vtxperim(k,2,3),1
!            call gridboxcalc(j,k,vtxperim(k,l,1),m,n,gridext,dbox,vtx,trees,vtxperim)
!          end do
!        end do
!    end if

!    if (vtxperim(k,2,2)-vtxperim(k,1,2) > 2*(tbuffer(k)+1)) then
!      do l=vtxperim(k,1,1)+(tbuffer(k)+1),vtxperim(k,2,1)-(tbuffer(k)+1),1
!        do m=0,tbuffer(k),1
!          do n=vtxperim(k,1,3),vtxperim(k,2,3),1
!            call gridboxcalc(j,k,l,vtxperim(k,1,2)+m,n,gridext,dbox,vtx,trees,vtxperim)
!            call gridboxcalc(j,k,l,vtxperim(k,2,2)-m,n,gridext,dbox,vtx,trees,vtxperim)
!          end do
!        end do
!      end do
!    else if (vtxperim(k,2,2)-vtxperim(k,1,2) > 0) then
!      do l=vtxperim(k,1,1)+(tbuffer(k)+1),vtxperim(k,2,1)-(tbuffer(k)+1),1
!        do m=vtxperim(k,1,2),vtxperim(k,2,2),1
!          do n=vtxperim(k,1,3),vtxperim(k,2,3),1
!            call gridboxcalc(j,k,l,m,n,gridext,dbox,vtx,trees,vtxperim)
!          end do
!        end do
!      end do
!    else
!      do l=vtxperim(k,1,1)+(tbuffer(k)+1),vtxperim(k,2,1)-(tbuffer(k)+1),1
!        m = 1
!          do n=vtxperim(k,1,3),vtxperim(k,2,3),1
!            call gridboxcalc(j,k,l,vtxperim(k,m,2),n,gridext,dbox,vtx,trees,vtxperim)
!          end do
!      end do
!    end if

!    if (vtxperim(k,2,3)-vtxperim(k,1,3) > 2*(tbuffer(k)+1)) then
!      do l=vtxperim(k,1,1)+(tbuffer(k)+1),vtxperim(k,2,1)-(tbuffer(k)+1),1
!        do m=vtxperim(k,1,2)+(tbuffer(k)+1),vtxperim(k,2,2)-(tbuffer(k)+1),1
!          do n=0,tbuffer(k),1
!            call gridboxcalc(j,k,l,m,vtxperim(k,1,3)+n,gridext,dbox,vtx,trees,vtxperim)
!            call gridboxcalc(j,k,l,m,vtxperim(k,2,3)-n,gridext,dbox,vtx,trees,vtxperim)
!          end do
!        end do
!      end do
!    else if (vtxperim(k,2,3)-vtxperim(k,1,3) > 0) then
!      do l=vtxperim(k,1,1)+(tbuffer(k)+1),vtxperim(k,2,1)-(tbuffer(k)+1),1
!        do m=vtxperim(k,1,2)+(tbuffer(k)+1),vtxperim(k,2,2)-(tbuffer(k)+1),1
!          do n=vtxperim(k,1,3),vtxperim(k,2,3),1
!            call gridboxcalc(j,k,l,m,n,gridext,dbox,vtx,trees,vtxperim)
!          end do
!        end do
!      end do
!    else
!      do l=vtxperim(k,1,1)+(tbuffer(k)+1),vtxperim(k,2,1)-(tbuffer(k)+1),1
!        do m=vtxperim(k,1,2)+(tbuffer(k)+1),vtxperim(k,2,2)-(tbuffer(k)+1),1
!          n = 1
!            call gridboxcalc(j,k,l,m,vtxperim(k,n,3),gridext,dbox,vtx,trees,vtxperim)
!        end do
!      end do
!    end if 

!  end subroutine ffieldlevel

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

!  subroutine ffieldlevel2(pos,k,gridext,dbox,trees,vtxperim,tbuffer,mscale,a,uinf,vel,vtens)
!    use gtmodinit
!    use vtxparticle
!    implicit none

!    type (level),dimension(:),intent(in) :: trees
!    logical,intent(in) :: mscale
!    real*8,intent(in) :: a
!    real*8,dimension(3),intent(in) :: uinf
!    integer,intent(in) :: k
!    real*8,dimension(3,2),intent(in) :: gridext
!    real*8,dimension(3),intent(inout) :: vel
!    real*8,dimension(3,3),intent(inout) :: vtens
!    real*8,dimension(3),intent(in) :: pos,dbox
!    integer,dimension(:,:,:),intent(in) :: vtxperim
!    integer,dimension(:),intent(in) :: tbuffer

!    integer :: l,m,n

!    if (vtxperim(k,2,1)-vtxperim(k,1,1) > 2*(tbuffer(k)+1)) then
!      do l=0,tbuffer(k),1
!        do m=vtxperim(k,1,2),vtxperim(k,2,2),1
!          do n=vtxperim(k,1,3),vtxperim(k,2,3),1
!            call gridboxcalc2(pos,k,vtxperim(k,1,1)+l,m,n,gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!            call gridboxcalc2(pos,k,vtxperim(k,2,1)-l,m,n,gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!          end do
!        end do
!      end do
!    else if (vtxperim(k,2,1)-vtxperim(k,1,1) > 0) then
!      do l=vtxperim(k,1,1),vtxperim(k,2,1),1
!        do m=vtxperim(k,1,2),vtxperim(k,2,2),1
!          do n=vtxperim(k,1,3),vtxperim(k,2,3),1
!            call gridboxcalc2(pos,k,l,m,n,gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!          end do
!        end do
!      end do
!    else
!      l = 1
!        do m=vtxperim(k,1,2),vtxperim(k,2,2),1
!          do n=vtxperim(k,1,3),vtxperim(k,2,3),1
!            call gridboxcalc2(pos,k,vtxperim(k,l,1),m,n,gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!          end do
!        end do
!    end if

!    if (vtxperim(k,2,2)-vtxperim(k,1,2) > 2*(tbuffer(k)+1)) then
!      do l=vtxperim(k,1,1)+(tbuffer(k)+1),vtxperim(k,2,1)-(tbuffer(k)+1),1
!        do m=0,tbuffer(k),1
!          do n=vtxperim(k,1,3),vtxperim(k,2,3),1
!            call gridboxcalc2(pos,k,l,vtxperim(k,1,2)+m,n,gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!            call gridboxcalc2(pos,k,l,vtxperim(k,2,2)-m,n,gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!          end do
!        end do
!      end do
!    else if (vtxperim(k,2,2)-vtxperim(k,1,2) > 0) then
!      do l=vtxperim(k,1,1)+(tbuffer(k)+1),vtxperim(k,2,1)-(tbuffer(k)+1),1
!        do m=vtxperim(k,1,2),vtxperim(k,2,2),1
!          do n=vtxperim(k,1,3),vtxperim(k,2,3),1
!            call gridboxcalc2(pos,k,l,m,n,gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!          end do
!        end do
!      end do
!    else
!      do l=vtxperim(k,1,1)+(tbuffer(k)+1),vtxperim(k,2,1)-(tbuffer(k)+1),1
!        m = 1
!          do n=vtxperim(k,1,3),vtxperim(k,2,3),1
!            call gridboxcalc2(pos,k,l,vtxperim(k,m,2),n,gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!          end do
!      end do
!    end if

!    if (vtxperim(k,2,3)-vtxperim(k,1,3) > 2*(tbuffer(k)+1)) then
!      do l=vtxperim(k,1,1)+(tbuffer(k)+1),vtxperim(k,2,1)-(tbuffer(k)+1),1
!        do m=vtxperim(k,1,2)+(tbuffer(k)+1),vtxperim(k,2,2)-(tbuffer(k)+1),1
!          do n=0,tbuffer(k),1
!            call gridboxcalc2(pos,k,l,m,vtxperim(k,1,3)+n,gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!            call gridboxcalc2(pos,k,l,m,vtxperim(k,2,3)-n,gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!          end do
!        end do
!      end do
!    else if (vtxperim(k,2,3)-vtxperim(k,1,3) > 0) then
!      do l=vtxperim(k,1,1)+(tbuffer(k)+1),vtxperim(k,2,1)-(tbuffer(k)+1),1
!        do m=vtxperim(k,1,2)+(tbuffer(k)+1),vtxperim(k,2,2)-(tbuffer(k)+1),1
!          do n=vtxperim(k,1,3),vtxperim(k,2,3),1
!            call gridboxcalc2(pos,k,l,m,n,gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!          end do
!        end do
!      end do
!    else
!      do l=vtxperim(k,1,1)+(tbuffer(k)+1),vtxperim(k,2,1)-(tbuffer(k)+1),1
!        do m=vtxperim(k,1,2)+(tbuffer(k)+1),vtxperim(k,2,2)-(tbuffer(k)+1),1
!          n = 1
!            call gridboxcalc2(pos,k,l,m,vtxperim(k,n,3),gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!        end do
!      end do
!    end if 

!  end subroutine ffieldlevel2

!!-----------------------------------------------------------------------------------------------------------------------------------
!!-----------------------------------------------------------------------------------------------------------------------------------

!  subroutine ffieldleveldipole(pos,k,gridext,dbox,trees,vtxperim,tbuffer,mscale,a,uinf,vel,vtens)
!    use gtmodinit
!    use vtxparticle
!    implicit none

!    type (level),dimension(:),intent(in) :: trees
!    logical,intent(in) :: mscale
!    real*8,intent(in) :: a
!    real*8,dimension(3),intent(in) :: uinf
!    integer,intent(in) :: k
!    real*8,dimension(3,2),intent(in) :: gridext
!    real*8,dimension(3),intent(inout) :: vel
!    real*8,dimension(3,3),intent(inout) :: vtens
!    real*8,dimension(3),intent(in) :: pos,dbox
!    integer,dimension(:,:,:),intent(in) :: vtxperim
!    integer,dimension(:),intent(in) :: tbuffer

!    integer :: l,m,n

!    if (vtxperim(k,2,1)-vtxperim(k,1,1) > 2*(tbuffer(k)+1)) then
!      do l=0,tbuffer(k),1
!        do m=vtxperim(k,1,2),vtxperim(k,2,2),1
!          do n=vtxperim(k,1,3),vtxperim(k,2,3),1
!            call gridboxcalcdipole(pos,k,vtxperim(k,1,1)+l,m,n,gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!            call gridboxcalcdipole(pos,k,vtxperim(k,2,1)-l,m,n,gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!          end do
!        end do
!      end do
!    else if (vtxperim(k,2,1)-vtxperim(k,1,1) > 0) then
!      do l=vtxperim(k,1,1),vtxperim(k,2,1),1
!        do m=vtxperim(k,1,2),vtxperim(k,2,2),1
!          do n=vtxperim(k,1,3),vtxperim(k,2,3),1
!            call gridboxcalcdipole(pos,k,l,m,n,gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!          end do
!        end do
!      end do
!    else
!      l = 1
!        do m=vtxperim(k,1,2),vtxperim(k,2,2),1
!          do n=vtxperim(k,1,3),vtxperim(k,2,3),1
!            call gridboxcalcdipole(pos,k,vtxperim(k,l,1),m,n,gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!          end do
!        end do
!    end if

!    if (vtxperim(k,2,2)-vtxperim(k,1,2) > 2*(tbuffer(k)+1)) then
!      do l=vtxperim(k,1,1)+(tbuffer(k)+1),vtxperim(k,2,1)-(tbuffer(k)+1),1
!        do m=0,tbuffer(k),1
!          do n=vtxperim(k,1,3),vtxperim(k,2,3),1
!            call gridboxcalcdipole(pos,k,l,vtxperim(k,1,2)+m,n,gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!            call gridboxcalcdipole(pos,k,l,vtxperim(k,2,2)-m,n,gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!          end do
!        end do
!      end do
!    else if (vtxperim(k,2,2)-vtxperim(k,1,2) > 0) then
!      do l=vtxperim(k,1,1)+(tbuffer(k)+1),vtxperim(k,2,1)-(tbuffer(k)+1),1
!        do m=vtxperim(k,1,2),vtxperim(k,2,2),1
!          do n=vtxperim(k,1,3),vtxperim(k,2,3),1
!            call gridboxcalcdipole(pos,k,l,m,n,gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!          end do
!        end do
!      end do
!    else
!      do l=vtxperim(k,1,1)+(tbuffer(k)+1),vtxperim(k,2,1)-(tbuffer(k)+1),1
!        m = 1
!          do n=vtxperim(k,1,3),vtxperim(k,2,3),1
!            call gridboxcalcdipole(pos,k,l,vtxperim(k,m,2),n,gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!          end do
!      end do
!    end if

!    if (vtxperim(k,2,3)-vtxperim(k,1,3) > 2*(tbuffer(k)+1)) then
!      do l=vtxperim(k,1,1)+(tbuffer(k)+1),vtxperim(k,2,1)-(tbuffer(k)+1),1
!        do m=vtxperim(k,1,2)+(tbuffer(k)+1),vtxperim(k,2,2)-(tbuffer(k)+1),1
!          do n=0,tbuffer(k),1
!            call gridboxcalcdipole(pos,k,l,m,vtxperim(k,1,3)+n,gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!            call gridboxcalcdipole(pos,k,l,m,vtxperim(k,2,3)-n,gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!          end do
!        end do
!      end do
!    else if (vtxperim(k,2,3)-vtxperim(k,1,3) > 0) then
!      do l=vtxperim(k,1,1)+(tbuffer(k)+1),vtxperim(k,2,1)-(tbuffer(k)+1),1
!        do m=vtxperim(k,1,2)+(tbuffer(k)+1),vtxperim(k,2,2)-(tbuffer(k)+1),1
!          do n=vtxperim(k,1,3),vtxperim(k,2,3),1
!            call gridboxcalcdipole(pos,k,l,m,n,gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!          end do
!        end do
!      end do
!    else
!      do l=vtxperim(k,1,1)+(tbuffer(k)+1),vtxperim(k,2,1)-(tbuffer(k)+1),1
!        do m=vtxperim(k,1,2)+(tbuffer(k)+1),vtxperim(k,2,2)-(tbuffer(k)+1),1
!          n = 1
!            call gridboxcalcdipole(pos,k,l,m,vtxperim(k,n,3),gridext,dbox,trees,vtxperim,mscale,a,uinf,vel,vtens)
!        end do
!      end do
!    end if 

!  end subroutine ffieldleveldipole

!!-----------------------------------------------------------------------------------------------------------------------------------
!!-----------------------------------------------------------------------------------------------------------------------------------

!end module gridtree

!!***********************************************************************************************************************************
!!***********************************************************************************************************************************

