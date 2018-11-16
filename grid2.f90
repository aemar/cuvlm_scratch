!***********************************************************************************************************************************
!***********************************************************************************************************************************

subroutine vtxinterp(rpos,vtxm,dbox,gridvtx)
  implicit none

  !rpos in relation to corner 1 (origin)
  real*8,dimension(3),intent(in) :: rpos,vtxm,dbox
  real*8,dimension(8,3),intent(inout) :: gridvtx

  real*8,dimension(8) :: frac8
  integer :: i

  call vtxinterp8(rpos,dbox,frac8)

  do i=1,8,1
    gridvtx(i,:) = gridvtx(i,:) + frac8(i) * vtxm
  end do

end subroutine vtxinterp

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

subroutine vtxinterp8(rpos,dbox,frac8)
  implicit none

  !rpos in relation to corner 1 (origin)
  real*8,dimension(3),intent(in) :: rpos
  !dbox, may vary for x,y,z
  real*8,dimension(3),intent(in) :: dbox
  !fraction allocated to each corner point
  real*8,dimension(8),intent(out) :: frac8

  frac8(1) = (rpos(1)/dbox(1))*(rpos(2)/dbox(2))*(rpos(3)/dbox(3))
  frac8(2) = (1.-rpos(1)/dbox(1))*(rpos(2)/dbox(2))*(rpos(3)/dbox(3))
  frac8(3) = (rpos(1)/dbox(1))*(1.-rpos(2)/dbox(2))*(rpos(3)/dbox(3))
  frac8(4) = (1.-rpos(1)/dbox(1))*(1.-rpos(2)/dbox(2))*(rpos(3)/dbox(3))
  frac8(5) = (rpos(1)/dbox(1))*(rpos(2)/dbox(2))*(1.-rpos(3)/dbox(3))
  frac8(6) = (1.-rpos(1)/dbox(1))*(rpos(2)/dbox(2))*(1.-rpos(3)/dbox(3))
  frac8(7) = (rpos(1)/dbox(1))*(1.-rpos(2)/dbox(2))*(1.-rpos(3)/dbox(3))
  frac8(8) = (1.-rpos(1)/dbox(1))*(1.-rpos(2)/dbox(2))*(1.-rpos(3)/dbox(3))

end subroutine vtxinterp8

!***********************************************************************************************************************************
!***********************************************************************************************************************************

subroutine findgridext(tsteps,ntraj,traj,np,points,gridext)
  use wingpanels
  implicit none

  integer,intent(in) :: tsteps,ntraj,np
  type(trajectory),dimension(ntraj),intent(in) :: traj
  type(spoint),dimension(np),intent(in) :: points

  real*8,dimension(3,2),intent(out) :: gridext

  real*8,dimension(3) :: dpoints
  real*8 :: maximum,minimum
  integer :: i,j,k

  !Find extensions of model
  do i=1,3,1
    maximum = points(1)%ac(i)+traj(points(1)%trj)%origin(i)
    minimum = points(1)%ac(i)+traj(points(1)%trj)%origin(i)
    do j=2,np,1
      if (points(j)%ac(i)+traj(points(j)%trj)%origin(i) > maximum) then
        maximum = points(j)%ac(i)+traj(points(j)%trj)%origin(i)
      end if
      if (points(j)%ac(i)+traj(points(j)%trj)%origin(i) < minimum) then
        minimum = points(j)%ac(i)+traj(points(j)%trj)%origin(i)
      end if
    end do
    dpoints(i) = maximum - minimum
  end do

  !Find maximum extensions of trajectory and add max model extensions plus safety factor 1.5 so it can't fly outside.
  do i=1,3,1
    maximum = traj(1)%step(1)%r(i)
    minimum = traj(1)%step(1)%r(i)
    do j=1,ntraj,1
      do k=1,tsteps,1
        if (traj(j)%step(k)%r(i) > maximum) then
          maximum = traj(j)%step(k)%r(i)
        end if
        if (traj(j)%step(k)%r(i) < minimum) then
          minimum = traj(j)%step(k)%r(i)
        end if
      end do
    end do
    gridext(i,1) = minimum - maxval(dpoints)*1.5
    gridext(i,2) = maximum + maxval(dpoints)*1.5
  end do

end subroutine findgridext

!***********************************************************************************************************************************
!***********************************************************************************************************************************

subroutine maingrid(gridext,dbox,n)
  implicit none

  real*8,dimension(3,2),intent(inout) :: gridext
  real*8,dimension(3),intent(in) :: dbox

  integer,dimension(3),intent(out) :: n

  integer :: i

  !calculate number of boxes in each direction
  do i=1,3,1
    n(i) = int((gridext(i,2)-gridext(i,1))/dbox(i))
  end do

!do this to check if grid has any influence on different cl outputs
do i=1,3,1
  gridext(i,1) = gridext(i,2) - n(i)*dbox(i)
end do

end subroutine maingrid

!***********************************************************************************************************************************
!***********************************************************************************************************************************

subroutine findpos3D(p,n,mingrid,dbox,check3,pbox)
  implicit none

  real*8,dimension(3),intent(in) :: p,mingrid
  integer,dimension(3),intent(in) :: n
  real*8,dimension(3),intent(in) :: dbox

  logical,intent(out) :: check3
  integer,dimension(3),intent(out) :: pbox

  logical :: check1
  integer :: i

  do i=1,3,1
    call findpos1D(p(i),n(i),mingrid(i),dbox(i),check1,pbox(i))
    if (check1 .eqv. .false.) then
      check3 = .false.
      exit
    end if
    check3 = .true.
  end do

end subroutine findpos3D

!***********************************************************************************************************************************
!***********************************************************************************************************************************

subroutine findpos1D(pi,ni,mini,dbox,check1,pboxi)
  implicit none

  real*8,intent(in) :: pi,mini
  integer,intent(in) :: ni
  real*8,intent(in) :: dbox

  logical,intent(out) :: check1
  integer,intent(out) :: pboxi

  integer :: lower,upper,mid

  !set box limits, midpoint box associated with lower set due to truncation of integer division
  lower = 1
  upper = ni
  mid = int((lower+upper)/2)

  if (pi >= mini .and. pi <= mini + dbox*upper) then
    do
      if (pi >= (lower-1)*dbox+mini .and. pi < mid*dbox+mini) then
        upper = mid
        mid = int((lower+upper)/2)
      else if (pi >= mid*dbox+mini .and. pi <= upper*dbox+mini) then
        lower = mid+1
        mid = int((lower+upper)/2)
      end if
      if (upper-lower == 0) then
        pboxi = lower
        check1 = .true.
        exit
      end if
    end do
  else
    pboxi = 0
    check1 = .false.
  end if

end subroutine findpos1D

!***********************************************************************************************************************************
!***********************************************************************************************************************************
