module functions
  implicit none

  contains

  subroutine traj(n_trj, t, r_0, tvec)
    use input, only: pi,krf,Mosc !,a
    implicit none
    integer,intent(in) :: n_trj
    real(kind=8),intent(in) :: t
    real(kind=8),dimension(3),intent(in) :: r_0
    real(kind=8),dimension(12),intent(out) :: tvec
    real(kind=8) :: om, alpha_wing
    real(kind=8) :: modr, M_trnsl
    real(kind=8),dimension(3) :: c_rot

real(kind=8) :: a = 340.


    !Create a separate file containing M and k

!Mosc = 0.6
!krf = 1.5


    M_trnsl = -1. * Mosc
    om = krf*Mosc*a/0.5
    alpha_wing =  1. * cos(om*t) * pi/180.

    !X,Y,Z VX,VY,VZ roll,pitch,yaw rollrate,pitchrate,yawrate
    c_rot(1) = 0.
    c_rot(2) = 0.
    c_rot(3) = 0.
    modr = sqrt( (r_0(1)-c_rot(1))**2 + (r_0(3)-c_rot(3))**2 )

    tvec(1:3) = traj_pos(n_trj, t, r_0)

    if (modr > 0.) then
    tvec(4) = -1.* pi/180. *modr * om * sin(om*t) * &
            & cos(alpha_wing + acos( (r_0(3)-c_rot(3))/modr ) * sign( dble(1.), (r_0(1)-c_rot(1)) ) ) + &
            & M_trnsl * a
    tvec(5) = 0.
    tvec(6) = modr * pi/180. * om * sin(om*t) * &
            & sin(alpha_wing + acos( (r_0(3)-c_rot(3))/modr ) * sign( dble(1.), (r_0(1)-c_rot(1)) ) )
    else
    tvec(4) = M_trnsl * a
    tvec(5) = 0.
    tvec(6) = 0.
    end if
    tvec(7) = 0.
    tvec(8) = alpha_wing
    tvec(9:10) = 0.
    tvec(11) = -1.* pi/180. * om * sin(om*t)
    tvec(12) = 0.
    
  end subroutine traj

  subroutine traj_accel(n_trj, t, r_0, taccel)
    use input, only: pi !,a
    implicit none
    integer,intent(in) :: n_trj
    real(kind=8),intent(in) :: t
    real(kind=8),dimension(3),intent(in) :: r_0
    real(kind=8),dimension(3),intent(out) :: taccel
    real(kind=8) :: om
    real(kind=8) :: modr, M_trnsl
    real(kind=8),dimension(3) :: c_rot

real(kind=8) :: a = 340.

    M_trnsl = -0.5
    om = -2.*pi/60. * 1000. !RPM

    !X,Y,Z VX,VY,VZ roll,pitch,yaw rollrate,pitchrate,yawrate
    if (n_trj == 1) then
      taccel(1:3) = 0.
    else !if (n_trj == 2) then
      if (n_trj == 2) then ![3.29, 3.6, 0.21]
        c_rot(1) = -3.29
        c_rot(2) = -3.6
        c_rot(3) = -0.21
      else if (n_trj == 3) then
        c_rot(1) = -3.29
        c_rot(2) = 3.6
        c_rot(3) = -0.21
      else
        print*, 'No trajectory specified - this is a bug'
      end if
      modr = sqrt( (r_0(2)-c_rot(2))**2 + (r_0(3)-c_rot(3))**2 )
      taccel(1) = 0.
      taccel(2) = -1.* modr * om**2 * cos(om*t + acos( (r_0(2)-c_rot(2))/modr ) * sign( dble(1.), (r_0(3)-c_rot(3)) ))
      taccel(3) = -1.* modr * om**2 * sin(om*t + acos( (r_0(2)-c_rot(2))/modr ) * sign( dble(1.), (r_0(3)-c_rot(3)) ))
    end if
    
  end subroutine traj_accel

  function traj_pos(n_trj, t, r_0)
    use input, only: pi,krf,Mosc !,a
    implicit none
    real(kind=8),dimension(3) :: traj_pos
    integer,intent(in) :: n_trj
    real(kind=8),intent(in) :: t
    real(kind=8),dimension(3),intent(in) :: r_0

    real(kind=8) :: om, alpha_wing
    real(kind=8) :: modr, M_trnsl
    real(kind=8),dimension(3) :: c_rot

real(kind=8) :: a = 340.

    M_trnsl = -1. * Mosc
    om = krf*Mosc*a/0.5
    alpha_wing = 1. * cos(om*t) * pi/180.


    !X,Y,Z VX,VY,VZ roll,pitch,yaw rollrate,pitchrate,yawrate
    c_rot(1) = 0.
    c_rot(2) = 0.
    c_rot(3) = 0.
    modr = sqrt( (r_0(1)-c_rot(1))**2 + (r_0(3)-c_rot(3))**2 )
    if (modr > 0.) then
      traj_pos(1) = modr * sin(alpha_wing + acos( (r_0(3)-c_rot(3))/modr ) * sign( dble(1.), (r_0(1)-c_rot(1)) ) ) + &
                  & M_trnsl * a * t + c_rot(1)
      traj_pos(2) = r_0(2)
      traj_pos(3) = modr * cos(alpha_wing + acos( (r_0(3)-c_rot(3))/modr ) * sign( dble(1.), (r_0(1)-c_rot(1)) ) ) + c_rot(3)
    else
      traj_pos(1) = c_rot(1) + M_trnsl * a * t
      traj_pos(2) = r_0(2)
      traj_pos(3) = c_rot(3)
    end if
    
  end function traj_pos

  subroutine traj_interp(pcpt, t, n_trj, r_0, te)
    use input, only: eoc, deltat !, a
    use maths
    implicit none
    real(kind=8),intent(out) :: te
    integer,intent(in) :: n_trj
    real(kind=8),intent(in) :: t
    real(kind=8),dimension(3),intent(in) :: pcpt, r_0

    real(kind=8) :: tn, tnm1, tnp1
    integer :: i

    real(kind=8),dimension(12) :: tvec1, tvec2
    real(kind=8),dimension(3) :: taccel

    real(kind=8) :: f1, f2, fp, fpp

real(kind=8) :: a = 340.


!    tn = t - sqrt(dot_product(pcpt-traj_pos(n_trj, dble(0.), r_0),pcpt-traj_pos(n_trj, dble(0.), r_0)))/a

    tn = 0.
    tnm1 = -1. * deltat

    call traj(n_trj, tn, r_0, tvec1)
    call traj(n_trj, tnm1, r_0, tvec2)

    f1 = (t-tn) - sqrt(dot_product(pcpt-tvec1(1:3),pcpt-tvec1(1:3)))/a
    f2 = (t-tnm1) - sqrt(dot_product(pcpt-tvec2(1:3),pcpt-tvec2(1:3)))/a
    fp = dot_product(pcpt-tvec1(1:3),tvec1(4:6))/(a*sqrt(dot_product(pcpt-tvec1(1:3),pcpt-tvec1(1:3)))) - 1.

!    print*

    do i=1,200,1

      tnp1 = tnm1 - ( (tnm1-tn) / (1. - (f1/f2) * ( ((f1-f2)/(tn-tnm1)) / fp)) )

      if (abs(tnp1-tn) < 0.0000001) exit

      tnm1 = tn
      tn = tnp1

      tvec2 = tvec1
      call traj(n_trj, tn, r_0, tvec1)

      f2 = f1
      f1 = (t-tn) - sqrt(dot_product(pcpt-tvec1(1:3),pcpt-tvec1(1:3)))/a     
      fp = dot_product(pcpt-tvec1(1:3),tvec1(4:6))/(a*sqrt(dot_product(pcpt-tvec1(1:3),pcpt-tvec1(1:3)))) - 1.

!      call traj(n_trj, tn, r_0, tvec2)

!      call traj_accel(n_trj, tn, r_0, taccel)

      !f(t_e) = 0

!      f = (t-tn) - sqrt(dot_product(pcpt-tvec(1:3),pcpt-tvec(1:3)))/a

!      fp = dot_product(pcpt-tvec(1:3),tvec(4:6))/(a*sqrt(dot_product(pcpt-tvec(1:3),pcpt-tvec(1:3)))) - 1.



!      fpp = (dot_product(pcpt-tvec(1:3),taccel)-dot_product(tvec(4:6),tvec(4:6)))/ &
!          & (a*sqrt(dot_product(pcpt-tvec(1:3),pcpt-tvec(1:3)))) + &
!          & dot_product(pcpt-tvec(1:3),tvec(4:6))**2 / (a*sqrt(dot_product(pcpt-tvec(1:3),pcpt-tvec(1:3)))**3)

!      !Halley's method - second order Newton-Raphson
!      tnp1 = tn - 2*f*fp/(2*fp**2 - f*fpp)

!      print*, (t-tn), sqrt(dot_product(pcpt-tvec1(1:3),pcpt-tvec1(1:3)))/a

      !Newton Raphson
!      tnp1 = tn - f/fp


      !Zeroth-order interpolation - does not converge!
!      tnp1 = t - sqrt(dot_product(pcpt-traj_pos(n_trj, tn, r_0),pcpt-traj_pos(n_trj, tn, r_0)))/a

!      if (abs(tnp1-tn) < 0.0000001) exit
      if (i == 400) then
        print*, 'Not converged', tnp1,tn
        print*
      end if
!      tn = tnp1

    end do

    te = tnp1

    

  end subroutine traj_interp



  subroutine traj_interp_o(pcpt, t, n_trj, r_0, te)
    use input, only: eoc !, a
    use maths
    implicit none
    real(kind=8),intent(out) :: te
    integer,intent(in) :: n_trj
    real(kind=8),intent(in) :: t
    real(kind=8),dimension(3),intent(in) :: pcpt, r_0

    real(kind=8) :: tn, tnp1
    integer :: i

    real(kind=8),dimension(12) :: tvec
    real(kind=8),dimension(3) :: taccel

    real(kind=8) :: f, fp, fpp

real(kind=8) :: a = 340.


!    tn = t - sqrt(dot_product(pcpt-traj_pos(n_trj, dble(0.), r_0),pcpt-traj_pos(n_trj, dble(0.), r_0)))/a

    tn = 0.

    print*

    do i=1,200,1

      call traj(n_trj, tn, r_0, tvec)

      call traj_accel(n_trj, tn, r_0, taccel)

      !f(t_e) = 0

      f = (t-tn) - sqrt(dot_product(pcpt-tvec(1:3),pcpt-tvec(1:3)))/a

      fp = dot_product(pcpt-tvec(1:3),tvec(4:6))/(a*sqrt(dot_product(pcpt-tvec(1:3),pcpt-tvec(1:3)))) - 1.



      fpp = (dot_product(pcpt-tvec(1:3),taccel)-dot_product(tvec(4:6),tvec(4:6)))/ &
          & (a*sqrt(dot_product(pcpt-tvec(1:3),pcpt-tvec(1:3)))) + &
          & dot_product(pcpt-tvec(1:3),tvec(4:6))**2 / (a*sqrt(dot_product(pcpt-tvec(1:3),pcpt-tvec(1:3)))**3)

      !Halley's method - second order Newton-Raphson
      tnp1 = tn - 2*f*fp/(2*fp**2 - f*fpp)

      print*, (t-tn), sqrt(dot_product(pcpt-tvec(1:3),pcpt-tvec(1:3)))/a

      !Newton Raphson
!      tnp1 = tn - f/fp


      !Zeroth-order interpolation - does not converge!
!      tnp1 = t - sqrt(dot_product(pcpt-traj_pos(n_trj, tn, r_0),pcpt-traj_pos(n_trj, tn, r_0)))/a

      if (abs(tnp1-tn) < 0.0000001) exit
      if (i == 400) then
        print*, 'Not converged', tnp1,tn
        print*, pcpt
        call traj(n_trj, t, r_0, tvec)
        print*, r_0
        print*, tvec(1:3), modulus(tvec(4:6))/340.
        print*
      end if
      tn = tnp1

    end do

    te = tnp1

    

  end subroutine traj_interp_o



  subroutine hyper_infl(x,rt_1,rt_2,circ,circdot,ntraj,uind)
    implicit none

    real(kind=8),dimension(3),intent(in) :: x,rt_1,rt_2
    real(kind=8),intent(in) :: circ,circdot
    integer,intent(in) :: ntraj
    real(kind=8),dimension(3),intent(out) :: uind

    integer :: nseg,i

    real(kind=8),dimension(3) :: rt

    nseg = 50

    !split into segments
    do i=1,nseg,1
      rt = rt_1 + (rt_2-rt_1)/(nseg-1) * (i-1)
    end do


  end subroutine 



!  subroutine findr_e(x_0, r_start, traj_type, r_e, v_e)
!    use module functions, only: traj1,traj2
!    use module input
!    implicit none
!    real(kind=8),dimension(3),intent(in) :: x_0, r_start
!    real(kind=8),dimension(3),intent(out) :: r_e, v_e
!    real(kind=8),dimension(9) :: tvec

!    real(kind=8) :: t_n,t_np1
!    integer :: i

!!    put some condition here to ensure convergence
!    t_n = -1. * sqrt(dot_product(x_0-r_start,x_0-r_start))/a

!    do i = 1,100,1

!      call traject(t_n, traj_type, r_0(traj_type), c_rot(traj_type), tvec)

!      f_t = a**2 * t_n**2 - dot_product(x_0 - tvec(1:3),x_0 - tvec(1:3))
!      ptf_t = 2. * a**2 * t_n + 2. * dot_product(x_0 - tvec(1:3),tvec(4:6))
!      pt2f_t = 2. * a**2 - 2. * dot_product(tvec(4:6),tvec(4:6)) + 2. * dot_product(x_0 - tvec(1:3),tvec(7:9))

!      t_np1 = t_n - (2. * f_t * ptf_t)/(2. * ptf_t**2 - f_t * pt2f_t)

!      if (abs(t_np1 - t_n) < eoc) then
!        exit
!      else
!        t_n = t_np1
!      end if

!    end do
    

!    call traject(t_np1, traj_type, r_0(traj_type), c_rot(traj_type), tvec)

!    r_e = tvec(1:3)
!    v_e = tvec(4:6)



!  end subroutine findr_e

!  function cross(v1,v2)
!    implicit none
!    real(kind=8),dimension(3) :: cross
!    real(kind=8),dimension(3),intent(in) :: v1,v2

!    cross(1)  = v1(2)*v2(3) - v1(3)*v2(2)
!    cross(2)  = v1(3)*v2(1) - v1(1)*v2(3)
!    cross(3)  = v1(1)*v2(2) - v1(2)*v2(1)

!  end function cross


end module functions
