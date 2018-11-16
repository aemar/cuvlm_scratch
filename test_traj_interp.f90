program test_traj_interp
  use functions
  use maths
  implicit none

  real*8,dimension(3,4) :: vert
  real*8,dimension(3) :: cpt,pout,p0
  real*8,dimension(3) :: vel = 0.
  integer :: npnt = 20
  integer :: i,j,iost

  real*8 :: te

  vert(1,1) = -1.
  vert(2,1) = 0.
  vert(3,1) = 0.
  vert(1,2) = -1.
  vert(2,2) = 1.
  vert(3,2) = 0.
  vert(1,3) = 0.
  vert(2,3) = 1.
  vert(3,3) = 0. 
  vert(1,4) = 0.
  vert(2,4) = 0.
  vert(3,4) = 0.

  cpt = (vert(:,1) + vert(:,2) + vert(:,3) + vert(:,4))/4.

print*, cpt

!  vel(1) = -170.

  open (unit=11,file='test_out.csv',status='replace',action='write',iostat=iost)

  do i=1,4,1
    do j = 1,npnt,1
      p0 = (vert(:,cycperm(4,1,i)) - vert(:,i))*(j-1)/(npnt-1) + vert(:,i)
      call traj_interp(cpt, dble(0.001), 1, p0, te)

      pout = traj_pos( 1, te, p0)

      write (11,'(F10.5,1A,1X,F10.5,1A,1X,F10.5,1A,1X,F10.5)') te,',',pout(1),',',pout(2),',',pout(3)

    end do
  end do


end program test_traj_interp
