module hyperbolic
use infl
  implicit none
  contains


  function velindhypseg_new(t,pc,n_trj,p10,p20,c1,c2,dc1,dc2,u1,u2,du1,du2)
    use input, only: a
    use maths, only: pi,cross
    use functions, only: traj_interp, traj_pos, traj
    implicit none
    real(kind=8),dimension(3) :: velindhypseg_new
    integer,intent(in) :: n_trj
    real(kind=8),dimension(3),intent(in) :: pc, p10, p20, u1, u2, du1, du2
    real(kind=8),intent(in) :: t, c1, c2, dc1, dc2
    real(kind=8),dimension(3) :: u, du, om, omhat, dom
    real(kind=8) :: te, c, dc
    real(kind=8),dimension(12) :: tvec

    integer :: i,nseg
    real(kind=8),dimension(3) :: d3r1,d3r2,d3r3,p0

    nseg = 20

    d3r1 = 0.

    do i=1,nseg,1

      p0 = (p20-p10) * (i-1)/(nseg-1) + p10

!  if (dot_product(pc-p0,pc-p0) < 1e-3) print*, dot_product(pc-p0,pc-p0), i

      call traj_interp(pc, t, n_trj, p0, te)

      !omhat: based on config at time of vortex particle!
      d3r2 = traj_pos(n_trj, te, p10)
      d3r3 = traj_pos(n_trj, te, p20)

!      omhat = (d3r3-d3r2)/sqrt(dot_product(d3r3-d3r2,d3r3-d3r2)) /(nseg-1)
      omhat = (d3r3-d3r2)/nseg

      c = (c2-c1)/(nseg-1) + c1
      om = omhat * c

      dc = (dc2-dc1)/(nseg-1) + dc1
      dom = omhat * dc

      u = (u2-u1)/(nseg-1) + u1
      du = (du2-du1)/(nseg-1) + du1

      call traj(n_trj, te, p0, tvec)
      
      d3r1 = d3r1 + uind(om,pc-tvec(1:3),tvec(4:6),dom,u,du)

    end do

    velindhypseg_new = d3r1

  end function velindhypseg_new


  function uind(om,r,vpart,omd,up,upd)
    use input, only: a,pi,eoc
    use maths, only: cross
    implicit none
    real(kind=8),dimension(3) :: uind
    real(kind=8),dimension(3),intent(in) :: om, r, vpart, omd, up, upd
    real(kind=8) :: modr, vdr, fvm

    modr = sqrt(dot_product(r,r))
    vdr = dot_product(vpart,r)
    fvm = 1.-vdr/(a*modr)

    if (modr > eoc) then
      uind = 1./(4.*pi) * (cross(om,r)/(fvm*modr**3) + vdr*cross(om,up)/(a**2*fvm*modr**3) + cross(omd,r)/(a*fvm*modr**2) + &
                         & (cross(om,upd)+cross(omd,up))/(a*fvm*modr))
    else
      uind = 0.
    end if

  end function uind

    function velindhypseg3D(pc,p10,p20,c1,c2,dcirc,vacfr,vfstream,uind,duind)
      use constants, only: a
      use maths, only: pi,cross
      implicit none
      real(kind=8),dimension(3) :: velindhypseg3D
      real(kind=8),dimension(3),intent(in) :: pc,p10,p20,vacfr,vfstream,uind,duind
      real(kind=8),intent(in) :: c1,c2,dcirc
      real(kind=8),dimension(3) :: p1,p2,s,rothat,sm,p1m,vseg,useg,d3r1,d3r2,rk,vtxprt
      real(kind=8) :: M,b,c,d,e,f,g,h,theta,fc,c1m,c2m,beta,f1,f2,f3,x0k,k,dgam,gam
      real(kind=8) :: rkmod,kstep,dxhdk,fvp,fvm,fum,fuvp,furvr
      real(kind=8) :: rm,pm,tm,qm,vdR,udR,udv,udu,vdv
      integer :: i,nseg

      nseg = 50

      M = sqrt(dot_product(vacfr,vacfr))/a
      theta = acos(vacfr(1)/sqrt(dot_product(vacfr,vacfr)))

      vseg = 0.
      vseg(1) = (-1.)*sqrt(dot_product(vfstream,vfstream))

!vseg = -1.*vseg

!      if (theta > 0.) then
!        rothat(1) = 0.
!        rothat(2) = (-1.)*vacfr(3)
!        rothat(3) = vacfr(2)
!        rothat = rothat/sqrt(dot_product(rothat,rothat))
!        p1 = rotate(p10-pc,rothat,theta)
!        p2 = rotate(p20-pc,rothat,theta)
!      else
        p1 = p10-pc
        p2 = p20-pc
!      end if

      beta = sqrt(1.-M**2)
      s = p2-p1
      sm = s

      kstep = 1./dble(nseg)
      dgam = c2-c1
      d3r1 = 0.

      do i=1,nseg,1
        k = (dble(i)-0.5)*kstep
        gam = dgam*k+c1
        g = gam*kstep/(4.*pi)
        rk(2) = s(2)*k+p1(2)
        rk(3) = s(3)*k+p1(3)
        x0k = s(1)*k+p1(1)
        rk(1) = (x0k+M*sqrt(x0k**2+beta**2*(rk(2)**2+rk(3)**2)))/beta**2
        rk = -1.*rk
        rm = sqrt(dot_product(rk,rk))
        vdr = dot_product(vseg,rk)
        vdv = dot_product(vseg,vseg)

        fvm = 1. - (vdR)/(a*rm)
        fvp = 1. + (vdR)/(a*rm)

        f1 = 1./(fvm*rm**3) !+ vdr/(a*fvm*rm**4)
        f2 = 0. !vdr/(a**2*fvm*rm**3) !1./(a*fvm*rm**2) - 1./(a*rm**2)

        d3r1(1) = d3r1(1) + g*(sm(2)*rk(3)-sm(3)*rk(2))*f1
        d3r1(2) = d3r1(2) + g*((sm(3)*rk(1)-sm(1)*rk(3))*f1 ) !+ sm(3)*vseg(1)*f2)
        d3r1(3) = d3r1(3) + g*((sm(1)*rk(2)-sm(2)*rk(1))*f1 ) !- sm(2)*vseg(1)*f2)

!        d3r1 = d3r1 + g*cross(sm,uind) * vdr/(a**2*fvm*rm**3) + dcirc*kstep*cross(sm,rk)/(4.*pi*a*fvm*rm**2) !&
             !& + (dcirc*kstep*cross(sm,uind) + gam*kstep*cross(sm,duind))/(4.*pi*a**2*fvm*rm)



!        d3r1 = d3r1 + 1./a * g*cross(-1. *vseg+utm1,sm) * ((-1.) * vdr/(a*fvm*rm**3))
!        d3r1 = d3r1 - g*cross(sm,vseg-uind) * vdr/(a**2*fvm*rm**3)

!        d3r1 = d3r1 + 2.*g*cross(sm,uind) * vdr/(a**2*fvm*rm**3) !&
!        & + 1./(a*fvm*rm) * (dcirc*kstep*cross(sm,rk))/(4.*pi*rm) + (dcirc*kstep* cross(sm,uind))/(4.*pi*a)/a !((dcirc*kstep*cross(sm,rk))/(4.*pi*rm) + g*cross(sm,duind)/a + (dcirc*kstep* cross(sm,uind))/(4.*pi*a))



      end do

      velindhypseg3D = d3r1




    end function velindhypseg3D

!###################################################################################################################################
!###################################################################################################################################

  function velindhypseg3D_vtx(pc,p10,v1,vvtx) !use vacfr for vvtx, or fraction thereof
    use constants, only: a
    use maths, only: pi,cross
    implicit none
    real(kind=8),dimension(3) :: velindhypseg3D_vtx
      real(kind=8),dimension(3),intent(in) :: pc,p10,v1,vvtx
      real(kind=8) :: modr,vdr,fvm

    if (dot_product(p10-pc,p10-pc) > 0.000001) then
      modr = sqrt(dot_product(p10-pc,p10-pc))
      vdr = dot_product(vvtx,p10-pc)
      fvm = 1. - vdr/(a*modr)
      velindhypseg3D_vtx = cross(v1,p10-pc) / (4.*pi * fvm * modr**3)
    else
      velindhypseg3D_vtx = 0.
    end if

  end function


    function velindhypseg3D_wake(pc,p10,p20,c1,c2,vacfr,utm1)
      use constants, only: a
      use maths, only: pi,cross
      implicit none
      real(kind=8),dimension(3) :: velindhypseg3D_wake
      real(kind=8),dimension(3),intent(in) :: pc,p10,p20,vacfr,utm1
      real(kind=8),intent(in) :: c1,c2
      real(kind=8),dimension(3) :: p1,p2,s,rothat,sm,p1m,useg,d3r1,d3r2,rk,vtxprt
      real(kind=8) :: M,b,c,d,e,f,g,h,theta,fc,c1m,c2m,beta,f1,f2,f3,x0k,k,dgam,gam
      real(kind=8) :: rkmod,kstep,dxhdk,fvp,fvm,fum,fuvp,furvr
      real(kind=8) :: rm,pm,tm,qm,vdR,udR,udv,udu,vdv
      integer :: i,nseg

      nseg = 20

      M = sqrt(dot_product(vacfr,vacfr))/a
      theta = acos((-1.)*vacfr(1)/sqrt(dot_product(vacfr,vacfr)))


!vseg = -1.*vseg

!      if (theta > 0.) then
!        rothat(1) = 0.
!        rothat(2) = (-1.)*vacfr(3)
!        rothat(3) = vacfr(2)
!        rothat = rothat/sqrt(dot_product(rothat,rothat))
!        p1 = rotate(p10-pc,rothat,theta)
!        p2 = rotate(p20-pc,rothat,theta)
!      else
        p1 = p10-pc
        p2 = p20-pc
!      end if

      beta = sqrt(1.-M**2)
      s = p2-p1
      sm = s

      kstep = 1./dble(nseg)
      dgam = c2-c1
      d3r1 = 0.

      do i=1,nseg,1
        k = (dble(i)-0.5)*kstep
        gam = dgam*k+c1
        g = gam*kstep/(4.*pi)
        rk(2) = s(2)*k+p1(2)
        rk(3) = s(3)*k+p1(3)
        x0k = s(1)*k+p1(1)
        rk(1) = (x0k+M*sqrt(x0k**2+beta**2*(rk(2)**2+rk(3)**2)))/beta**2
        rk = -1.*rk
        rm = sqrt(dot_product(rk,rk))

!vdr = dot_product(vacfr,rk)
!fvm = 1. - vdr/(a*rm)

        f1 = 1./(rm**3) !+ vdr/(a*rm**4)

        d3r1(1) = d3r1(1) + g*(sm(2)*rk(3)-sm(3)*rk(2))*f1
        d3r1(2) = d3r1(2) + g*(sm(3)*rk(1)-sm(1)*rk(3))*f1
        d3r1(3) = d3r1(3) + g*(sm(1)*rk(2)-sm(2)*rk(1))*f1

      end do

      velindhypseg3D_wake = d3r1

    contains
      function rotate(rin,ax,ang)
        implicit none
        real(kind=8),dimension(3) :: rotate
        real(kind=8),dimension(3),intent(in) :: rin,ax
        real(kind=8),intent(in) :: ang
        real(kind=8),dimension(3,3) :: rmat

        rmat(1,1) = cos(ang)+ax(1)**2*(1-cos(ang))
        rmat(1,2) = ax(1)*ax(2)*(1-cos(ang))-ax(3)*sin(ang)
        rmat(1,3) = ax(1)*ax(3)*(1-cos(ang))+ax(2)*sin(ang)
        rmat(2,1) = ax(2)*ax(1)*(1-cos(ang))+ax(3)*sin(ang)
        rmat(2,2) = cos(ang)+ax(2)**2*(1-cos(ang))
        rmat(2,3) = ax(2)*ax(3)*(1-cos(ang))-ax(1)*sin(ang)
        rmat(3,1) = ax(3)*ax(1)*(1-cos(ang))-ax(2)*sin(ang)
        rmat(3,2) = ax(3)*ax(2)*(1-cos(ang))+ax(1)*sin(ang)
        rmat(3,3) = cos(ang)+ax(3)**2*(1-cos(ang))
      end function rotate

    end  function velindhypseg3D_wake



end module hyperbolic
