!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!Merges two vortices and finds a new weighted position and vortex moment.

subroutine mergevtx(p1,m1,p2,m2,pout,mout)
  implicit none

  real*8,dimension(3),intent(in) :: p1,m1,p2,m2
  real*8,dimension(3),intent(out) :: pout,mout
  real*8 :: w1,w2,frac1,frac2

  frac1 = sqrt(dot_product(m2,m2))/sqrt(dot_product(m1,m1))
  frac2 = 1./frac1

  w1 = 1./(1.+frac1)
  w2 = 1./(1.+frac2)

  pout = w1*p1 + w2*p2
  mout = m1 + m2

end subroutine

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

logical function incheck(p,box)
  implicit none

  real*8,dimension(8,3),intent(in) :: box
  real*8,dimension(3),intent(in) :: p

  if (dot_product(p-box(1,:),p-box(8,:)) < 0.) then
    if (dot_product(p-box(3,:),p-box(6,:)) < 0.) then
      if (dot_product(p-box(7,:),p-box(2,:)) < 0.) then
        if (dot_product(p-box(5,:),p-box(4,:)) < 0.) then
          incheck = .true.
        else
          incheck = .false.
        end if
      else
        incheck = .false.
      end if
    else
      incheck = .false.
    end if
  else
    incheck = .false.
  end if

end function incheck

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

!Calculates the induced velocity of a vortex particle on a point at position p.
subroutine vtxinduced(p,vpos,vtx,vtxind)
  use maths, only: cross,cycperm,pi,modulus
  implicit none

  real*8,dimension(3),intent(in) :: p,vpos,vtx
  real*8,dimension(3),intent(out) :: vtxind

  real*8 :: rho,rho2,den,modr,sqr2p1

!  if (dot_product(vpos-p,vpos-p) > 10**(-3)) then
!    vtxind = cross(vtx,p-vpos)/(4*pi*sqrt(dot_product(p-vpos,p-vpos))**3)

  !Higher order algebraic smoothing, Winckelmanns
  if (dot_product(vpos-p,vpos-p) > 0.1) then
    modr = sqrt(dot_product(p-vpos,p-vpos))
    rho = modr * 2.5
    rho2 = rho*rho
    sqr2p1 = sqrt((rho2 + 1.))
    den = 4.*pi*sqr2p1*sqr2p1*sqr2p1*sqr2p1*sqr2p1*modr*modr*modr
    vtxind = (cross(vtx,p-vpos) * (rho2*rho) * (rho2 + 2.5))/den
  else
    vtxind = 0.
  end if

  if (modulus(vtxind) > 340.) then
    vtxind = vtxind/modulus(vtxind) * 340.
  end if

end subroutine vtxinduced

!Calculates the induced velocity of a vortex particle on a point at position p.
subroutine vtxinduced_wing(p,vpos,vtx,vtxind)
  use maths, only: cross,cycperm,pi
  implicit none

  real*8,dimension(3),intent(in) :: p,vpos,vtx
  real*8,dimension(3),intent(out) :: vtxind

  real*8 :: rho,rho2,den,modr,sqr2p1

!  if (dot_product(vpos-p,vpos-p) > 10**(-3)) then
!    vtxind = cross(vtx,p-vpos)/(4*pi*sqrt(dot_product(p-vpos,p-vpos))**3)

  !Higher order algebraic smoothing, Winckelmanns
  if (dot_product(vpos-p,vpos-p) > 0.00001) then
    vtxind = cross(vtx,p-vpos)/(4*pi*sqrt(dot_product(p-vpos,p-vpos))**3)
  else
    vtxind = 0.
  end if

end subroutine vtxinduced_wing
