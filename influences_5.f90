subroutine vtxwingcalcpart(pin,vind,gradind,vfstream,vtxpnt,dvtx,d1i1,deltat,i)
  use basictypes
  use vtxparticle
  use pg
  use constants
  implicit none

  real(kind=8),dimension(3),intent(in) :: pin
  real(kind=8),dimension(3),intent(inout) :: vind
  real(kind=8),dimension(3,3),intent(inout) :: gradind
  real(kind=8),dimension(3),intent(in) :: vfstream
  real(kind=8),intent(in) :: deltat
  type(point),dimension(dvtx),intent(in) :: vtxpnt
  integer,intent(in) :: dvtx,d1i1,i

  real(kind=8),dimension(3) :: d3r1,d3r2,d3r4


  if (vtxpnt(d1i1)%ffl) then  
    call ittodummy(pin,vtxpnt(d1i1)%set(i)%infr(1:3),vfstream,deltat,d3r1)
  else
    d3r1 = intp(vtxpnt(d1i1)%set(vtxpnt(d1i1)%ti(1))%infr(1:3), &
         & vtxpnt(d1i1)%set(vtxpnt(d1i1)%ti(2))%infr(1:3),vtxpnt(d1i1)%tr(1),vtxpnt(d1i1)%tr(2))
  end if

  d3r2 = intp(vtxpnt(d1i1)%set(vtxpnt(d1i1)%ti(1))%minfr, &
         & vtxpnt(d1i1)%set(vtxpnt(d1i1)%ti(2))%minfr,vtxpnt(d1i1)%tr(1),vtxpnt(d1i1)%tr(2))

  call vtxinduced(pin,d3r1,d3r2,d3r4)
  vind = vind + d3r4

!  call shearind(pin,d3r1,d3r2,d33r1)
!  gradind = gradind + d33r1

end subroutine vtxwingcalcpart

subroutine vtxwingcalcpart_wake_2(pin,wingpan,vtxpnt,dvtx,d1i1,deltat,i)
  use basictypes
  use vtxparticle
  use pg
  use constants
  use hyperbolic, only: velindhypseg3D_vtx
  use maths, only: cross, modulus
  implicit none

  real(kind=8),dimension(3),intent(in) :: pin
!  real(kind=8),dimension(3),intent(inout) :: vind
  type(panel),intent(inout) :: wingpan
  real(kind=8),intent(in) :: deltat
  integer,intent(in) :: dvtx,d1i1,i
  type(point),dimension(dvtx),intent(in) :: vtxpnt
  real(kind=8),dimension(3) :: norm,d3r1,d3r2,d3r3,d3r4,ahat
  real(kind=8) :: d1r1

  if (vtxpnt(d1i1)%ti(2) >= 1) then


    norm = cross(wingpan%pnt(3)%p%set(i)%infr(1:3)-wingpan%pnt(1)%p%set(i)%infr(1:3), &
               & wingpan%pnt(2)%p%set(i)%infr(1:3)-wingpan%pnt(4)%p%set(i)%infr(1:3))

    !normalise to unit vector
    norm = norm/modulus(norm)  

    d3r1 = intp(vtxpnt(d1i1)%set(vtxpnt(d1i1)%ti(1))%infr(1:3), &
         & vtxpnt(d1i1)%set(vtxpnt(d1i1)%ti(2))%infr(1:3),vtxpnt(d1i1)%tr(1),vtxpnt(d1i1)%tr(2))

    d3r2 = intp(vtxpnt(d1i1)%set(vtxpnt(d1i1)%ti(1))%minfr, &
         & vtxpnt(d1i1)%set(vtxpnt(d1i1)%ti(2))%minfr,vtxpnt(d1i1)%tr(1),vtxpnt(d1i1)%tr(2))

    d1r1 = dot_product(norm, d3r1-pin)

!    if (modulus(d3r1-pin) < 0.5  * sqrt(wingpan%area)) then

!      d3r4 = 0.

    if (abs(d1r1) < 0.3 * sqrt(wingpan%area) .and. modulus(d3r1-pin) < 1.5  * sqrt(wingpan%area)) then

      ahat = (d3r1-pin)/modulus(d3r1-pin)
      d3r3 = pin + 0.5 * sqrt(wingpan%area)/abs(dot_product(ahat,norm)) * ahat
!      d3r3 = pin + sign(dble(1.),d1r1) * 0.3 * sqrt(wingpan%area)/dot_product(ahat,norm) * ahat
      call vtxinduced(pin,d3r3,d3r2,d3r4)
!      d3r4 = 0.
    else

      call vtxinduced(pin,d3r1,d3r2,d3r4)

    end if

    wingpan%vind = wingpan%vind + d3r4

!  call shearind(pin,d3r1,d3r2,d33r1)
!  gradind = gradind + d33r1

  end if

end subroutine vtxwingcalcpart_wake_2

subroutine vtxwingcalcpart_wake(pin,vind,vacfr,vtxpnt,dvtx,d1i1,deltat,i)
  use basictypes
  use vtxparticle
  use pg
  use constants
  use hyperbolic, only: velindhypseg3D_vtx
  implicit none

  real(kind=8),dimension(3),intent(in) :: pin
  real(kind=8),dimension(3),intent(inout) :: vind
  real(kind=8),dimension(3),intent(in) :: vacfr
  real(kind=8),intent(in) :: deltat
  integer,intent(in) :: dvtx,d1i1,i
  type(point),dimension(dvtx),intent(in) :: vtxpnt
  real(kind=8),dimension(3) :: d3r1,d3r2,d3r4

  if (vtxpnt(d1i1)%ffl) then  
    call ittodummy(pin,vtxpnt(d1i1)%set(i)%infr(1:3),(-1.)*vacfr,deltat,d3r1)
  else
    d3r1 = intp(vtxpnt(d1i1)%set(vtxpnt(d1i1)%ti(1))%infr(1:3), &
         & vtxpnt(d1i1)%set(vtxpnt(d1i1)%ti(2))%infr(1:3),vtxpnt(d1i1)%tr(1),vtxpnt(d1i1)%tr(2))
  end if

  d3r2 = intp(vtxpnt(d1i1)%set(vtxpnt(d1i1)%ti(1))%minfr, &
         & vtxpnt(d1i1)%set(vtxpnt(d1i1)%ti(2))%minfr,vtxpnt(d1i1)%tr(1),vtxpnt(d1i1)%tr(2))

  call vtxinduced_wing(pin,d3r1,d3r2,d3r4)
  vind = vind + d3r4

!  vind = vind + velindhypseg3D_vtx(pin,d3r1,d3r2,vacfr)
!  vind = vind + velindhypseg3D_vtx(pin,vtxpnt(d1i1)%set(i)%infr(1:3),vtxpnt(d1i1)%set(i)%minfr,vacfr)

!  call shearind(pin,d3r1,d3r2,d33r1)
!  gradind = gradind + d33r1

end subroutine vtxwingcalcpart_wake

!###################################################################################################################################
!###################################################################################################################################

subroutine wingwingcalcseg_3(pc,pan,cmat,uvind,duvind,i,j,k,q)
  use input
  use basictypes
  use constants
  use maths
  use hyperbolic
  use functions
  implicit none

  type(panel),dimension(npan),intent(inout) :: pan
  real(kind=8),dimension(3,npan,npan),intent(inout) :: cmat
  real(kind=8),dimension(3),intent(in) :: pc
  real(kind=8),dimension(3,npan,tsteps),intent(in) :: uvind,duvind
  integer,intent(in) :: i,j,k,q
  real(kind=8),dimension(3) :: d3r1,d3r2,d3r3,d3r4,u1,u2,du1,du2
  integer :: d1i1
  real(kind=8) :: d1r1,d1r2,d1r3,c1,c2,dc1,dc2

  integer :: nseg, iseg, tmin, n_trj
  real(kind=8) :: t, te, c, dc
  real(kind=8),dimension(3) :: p0, p10, p20
  real(kind=8),dimension(3) :: u, du, om, omhat, dom
  real(kind=8),dimension(12) :: tvec

!  dc1 = 0.
!  dc2 = 0.
!  u1 = 0.
!  u2 = 0.
!  du1 = 0.
!  du2 = 0.


  nseg = 20

  n_trj = pan(k)%trj

  p10 = pan(k)%pnt(q)%p%acfr
  p20 = pan(k)%pnt(cycperm(4,1,q))%p%acfr

  t = i*deltat
  
  do iseg=1,nseg,1

    p0 = (p20-p10)*(iseg-1)/(nseg-1) + p10

    call traj_interp(pc, t, n_trj, p0, te)

    !omhat: based on config at time of vortex particle!
    d3r2 = traj_pos(n_trj, te, p10)
    d3r3 = traj_pos(n_trj, te, p20)

!Comment out the (nseg-1), to increase the Cl, this isn't correct, why does it not show otherwise???

!    omhat = ((d3r3-d3r2)/sqrt(dot_product(d3r3-d3r2,d3r3-d3r2))) / (nseg-1)
   omhat = (d3r3-d3r2)/ (nseg)

    

!      dc = (dc2-dc1)/(nseg-1) + dc1
!      dom = omhat * dc
!      u = (u2-u1)/(nseg-1) + u1
!      du = (du2-du1)/(nseg-1) + du1

      dom = 0.
      u = 0.
      du = 0.

    !check if present time step -> LHS, known -> RHS, split point as necessary and add to cmat or induced velocity

    if (i==1) then

      call traj(n_trj, te, p0, tvec)

      d3r1 = uind(omhat,pc-tvec(1:3),tvec(4:6),dom,u,du) 

      cmat(:,k,j) = cmat(:,k,j) + d3r1
      if (pan(k)%neigh(q) > 0) then
        if (pan(k)%neigh(q) /= j) then
          cmat(:,pan(k)%neigh(q),j) = cmat(:,pan(k)%neigh(q),j) - d3r1
        end if
      end if

    else

      if (te >= deltat .and. te <= (i-1)*deltat) then

        tmin = floor(te/deltat)

        if (pan(k)%neigh(q) /= j) then
          c = (te/deltat-dble(tmin)) * pan(k)%dcirc(q,tmin + 1) + (1. - (te/deltat-dble(tmin))) * pan(k)%dcirc(q,tmin)
        else
          c = (te/deltat-dble(tmin)) * pan(k)%circ(tmin + 1) + (1. - (te/deltat-dble(tmin))) * pan(k)%circ(tmin)
        end if

        om = omhat * c

        call traj(n_trj, te, p0, tvec)
      
        pan(j)%vind = pan(j)%vind + uind(om,pc-tvec(1:3),tvec(4:6),dom,u,du)

      else if (te > (i-1)*deltat) then

        tmin = i-1

        call traj(n_trj, te, p0, tvec)

        d3r1 = uind(omhat,pc-tvec(1:3),tvec(4:6),dom,u,du)

        if (pan(k)%neigh(q) /= j) then
          c = (1. - (te/deltat-tmin)) * pan(k)%dcirc(q,tmin)
        else
          c = (1. - (te/deltat-tmin)) * pan(k)%circ(tmin)
        end if

        pan(j)%vind = pan(j)%vind + c * d3r1

        cmat(:,k,j) = cmat(:,k,j) + (te/deltat-tmin) * d3r1
        if (pan(k)%neigh(q) > 0) then
          if (pan(k)%neigh(q) /= j) then
            cmat(:,pan(k)%neigh(q),j) = cmat(:,pan(k)%neigh(q),j) - (te/deltat-tmin) * d3r1
          end if
        end if

      else if (te < deltat) then

        tmin = min(i-1,5)

        if (pan(k)%neigh(q) /= j) then
          c = pan(k)%dcirc(q,tmin)
        else
          c = pan(k)%circ(tmin)
        end if

        om = omhat * c

        call traj(n_trj, te, p0, tvec)
      
        pan(j)%vind = pan(j)%vind + uind(om,pc-tvec(1:3),tvec(4:6),dom,u,du)

      end if

    end if

  end do

end subroutine wingwingcalcseg_3

!###################################################################################################################################
!###################################################################################################################################

subroutine wingwakecalcseg_3(pc,pan,uvind,duvind,i,k,q,vind)
  use input
  use basictypes
  use constants
  use maths
  use hyperbolic
  use functions
  implicit none

  real(kind=8),dimension(3),intent(inout) :: vind
  type(panel),dimension(npan),intent(inout) :: pan
  real(kind=8),dimension(3),intent(in) :: pc
  real(kind=8),dimension(3,npan,tsteps),intent(in) :: uvind,duvind
  integer,intent(in) :: i,k,q
  real(kind=8),dimension(3) :: d3r1,d3r2,d3r3,d3r4,u1,u2,du1,du2
  integer :: d1i1
  real(kind=8) :: d1r1,d1r2,d1r3,c1,c2,dc1,dc2

  integer :: nseg, iseg, tmin, n_trj
  real(kind=8) :: t, te, c, dc
  real(kind=8),dimension(3) :: p0, p10, p20
  real(kind=8),dimension(3) :: u, du, om, omhat, dom
  real(kind=8),dimension(12) :: tvec

!  dc1 = 0.
!  dc2 = 0.
!  u1 = 0.
!  u2 = 0.
!  du1 = 0.
!  du2 = 0.

  nseg = 20

  n_trj = pan(k)%trj

  p10 = pan(k)%pnt(q)%p%acfr
  p20 = pan(k)%pnt(cycperm(4,1,q))%p%acfr

  t = i*deltat
  
  do iseg=1,nseg,1

    p0 = (p20-p10)*(iseg-1)/(nseg-1) + p10

    call traj_interp(pc, t, n_trj, p0, te)

    !omhat: based on config at time of vortex particle!
    d3r2 = traj_pos(n_trj, te, p10)
    d3r3 = traj_pos(n_trj, te, p20)

!Comment out the (nseg-1), to increase the Cl, this isn't correct, why does it not show otherwise???

!    omhat = ((d3r3-d3r2)/sqrt(dot_product(d3r3-d3r2,d3r3-d3r2))) / (nseg-1)
    omhat = (d3r3-d3r2) / nseg

    

!      dc = (dc2-dc1)/(nseg-1) + dc1
!      dom = omhat * dc
!      u = (u2-u1)/(nseg-1) + u1
!      du = (du2-du1)/(nseg-1) + du1

      dom = 0.
      u = 0.
      du = 0.

    !check if present time step -> LHS, known -> RHS, split point as necessary and add to cmat or induced velocity


    if (te >= deltat) then

      tmin = floor(te/deltat)

      c = (te/deltat-dble(tmin)) * pan(k)%dcirc(q,tmin + 1) + (1. - (te/deltat-dble(tmin))) * pan(k)%dcirc(q,tmin)

      om = omhat * c

      call traj(n_trj, te, p0, tvec)
      
      vind = vind + uind(om,pc-tvec(1:3),tvec(4:6),dom,u,du)

!    else if (te < deltat) then

!      tmin = min(i-1,5)

!      if (pan(k)%neigh(q) /= j) then
!        c = pan(k)%dcirc(q,tmin)
!      else
!        c = pan(k)%circ(tmin)
!      end if

!      om = omhat * c

!      call traj(n_trj, te, p0, tvec)
     
!      wakepnt(j)%vind = wakepnt(j)%vind + uind(om,pc-tvec(1:3),tvec(4:6),dom,u,du)

    end if

  end do

end subroutine wingwakecalcseg_3

!###################################################################################################################################
!###################################################################################################################################

subroutine wingwingcalcseg_2(pcpt,n_trj,pan,cmat,uvind,duvind,i,j,k,q)
  use input
  use basictypes
  use pg
  use infl
  use constants
  use maths
  use hyperbolic
  implicit none

  type(panel),dimension(npan),intent(inout) :: pan
  real(kind=8),dimension(3,npan,npan),intent(inout) :: cmat
  real(kind=8),dimension(3),intent(in) :: pcpt
  real(kind=8),dimension(3,npan,tsteps),intent(in) :: uvind,duvind
  integer,intent(in) :: i,j,k,q,n_trj
  real(kind=8),dimension(3) :: d3r1,d3r2,d3r3,d3r4,u1,u2,du1,du2
  integer :: d1i1
  real(kind=8) :: d1r1,d1r2,d1r3,c1,c2,dc1,dc2
  real(kind=8),parameter :: rzero = 0.
  real(kind=8),parameter :: runity = 1. 
  real(kind=8),dimension(3),parameter :: zerovec = 0.

  integer :: idtime1,idtime2

  dc1 = 0.
  dc2 = 0.
  u1 = 0.
  u2 = 0.
  du1 = 0.
  du2 = 0.

  d1r1 = pan(k)%pnt(q)%p%te/deltat
  d1r2 = pan(k)%pnt(cycperm(4,1,q))%p%te/deltat

  idtime1 = max(1,int(floor(d1r1)))
  idtime2 = max(1,int(floor(d1r2)))

  if (i==1) then

      d3r4 = velindhypseg_new( i*deltat, pcpt, n_trj, pan(k)%pnt(q)%p%acfr, &
        & pan(k)%pnt(cycperm(4,1,q))%p%acfr, dble(1.), dble(1.), dc1, dc2, u1, u2, du1, du2)

      cmat(:,k,j) = cmat(:,k,j) + d3r4
      if (pan(k)%neigh(q) > 0) then
        if (pan(k)%neigh(q) /= j) then
          cmat(:,pan(k)%neigh(q),j) = cmat(:,pan(k)%neigh(q),j) - d3r4
        end if
      end if

  else


  if (dble(i-1)*deltat > pan(k)%pnt(q)%p%te) then

    if (pan(k)%neigh(q) /= j) then
      c1 = (d1r1 - floor(d1r1)) * pan(k)%dcirc(q,idtime1 + 1) + (1. - (d1r1 - floor(d1r1))) * pan(k)%dcirc(q,idtime1)
    else
      c1 = (d1r1 - floor(d1r1)) * pan(k)%circ(idtime1 + 1) + (1. - (d1r1 - floor(d1r1))) * pan(k)%circ(idtime1)
    end if

    if (dble(i-1)*deltat > pan(k)%pnt(cycperm(4,1,q))%p%te) then

      if (pan(k)%neigh(q) /= j) then
        c2 = (d1r2 - floor(d1r2)) * pan(k)%dcirc(q,idtime2 + 1) + (1. - (d1r2 - floor(d1r2))) * pan(k)%dcirc(q,idtime2)
      else
        c2 = (d1r2 - floor(d1r2)) * pan(k)%circ(idtime2 + 1) + (1. - (d1r2 - floor(d1r2))) * pan(k)%circ(idtime2)
      end if

      pan(j)%vind = pan(j)%vind + velindhypseg_new( i*deltat, pcpt, n_trj, pan(k)%pnt(q)%p%acfr, &
        & pan(k)%pnt(cycperm(4,1,q))%p%acfr, c1, c2, dc1, dc2, u1, u2, du1, du2)

    else

      d1i1 = (i-1) - (idtime1+1) !if difference is greater than one timestep

      d1r3 = ints(runity,(1.-(d1r1-floor(d1r1))+dble(d1i1)),(d1r2-floor(d1r2)))

      d3r3 = pan(k)%pnt(cycperm(4,1,q))%p%acfr - d1r3*(pan(k)%pnt(cycperm(4,1,q))%p%acfr-pan(k)%pnt(q)%p%acfr)

      d3r4 = velindhypseg_new( i*deltat, pcpt,n_trj, d3r3, &
        & pan(k)%pnt(cycperm(4,1,q))%p%acfr, rzero, (d1r2-floor(d1r2)), dc1, dc2, u1, u2, du1, du2)

      cmat(:,k,j) = cmat(:,k,j) + d3r4
      if (pan(k)%neigh(q) > 0) then
        if (pan(k)%neigh(q) /= j) then
          cmat(:,pan(k)%neigh(q),j) = cmat(:,pan(k)%neigh(q),j) - d3r4
        end if
      end if

      if (pan(k)%neigh(q) /= j) then
        c2 = pan(k)%dcirc(q,idtime2)
      else
        c2 = pan(k)%circ(idtime2)
      end if

      !add first seg
      pan(j)%vind = pan(j)%vind + velindhypseg_new( i*deltat, pcpt,n_trj, d3r3, pan(k)%pnt(cycperm(4,1,q))%p%acfr, &
        & c2, (1.-(d1r2-floor(d1r2))) * c2, dc1, dc2, u1, u2, du1, du2)

      !add free seg
      pan(j)%vind = pan(j)%vind + velindhypseg_new( i*deltat, pcpt,n_trj, pan(k)%pnt(q)%p%acfr, d3r3, &
        & c1, c2, dc1, dc2, u1, u2, du1, du2)

    end if      

  else

    if (dble(i-1)*deltat > pan(k)%pnt(cycperm(4,1,q))%p%te) then

      d1i1 = (i-1) - (idtime2+1) !if difference is greater than one timestep

      d1r3 = ints(runity,(1.-(d1r2-floor(d1r2))+dble(d1i1)),(d1r1-floor(d1r1)))

      d3r3 = pan(k)%pnt(q)%p%acfr + d1r3*(pan(k)%pnt(cycperm(4,1,q))%p%acfr-pan(k)%pnt(q)%p%acfr)

      d3r4 = velindhypseg_new( i*deltat, pcpt,n_trj, pan(k)%pnt(q)%p%acfr, d3r3, &
        & (d1r1-floor(d1r1)), rzero, dc1, dc2, u1, u2, du1, du2)

      cmat(:,k,j) = cmat(:,k,j) + d3r4
      if (pan(k)%neigh(q) > 0) then
        if (pan(k)%neigh(q) /= j) then
          cmat(:,pan(k)%neigh(q),j) = cmat(:,pan(k)%neigh(q),j) - d3r4
        end if
      end if

      if (pan(k)%neigh(q) /= j) then
        c2 = (d1r2 - floor(d1r2)) * pan(k)%dcirc(q,idtime2 + 1) + (1. - (d1r2 - floor(d1r2))) * pan(k)%dcirc(q,idtime2)
        c1 = pan(k)%dcirc(q,idtime1)
      else
        c2 = (d1r2 - floor(d1r2)) * pan(k)%circ(idtime2 + 1) + (1. - (d1r2 - floor(d1r2))) * pan(k)%circ(idtime2)
        c1 = pan(k)%circ(idtime1)
      end if

      !add first seg
      pan(j)%vind = pan(j)%vind + velindhypseg_new( i*deltat, pcpt, n_trj, pan(k)%pnt(q)%p%acfr, d3r3, &
        & (1.-(d1r1-floor(d1r1))) * c1, c1, dc1, dc2, u1, u2, du1, du2)

      !add free seg
      pan(j)%vind = pan(j)%vind + velindhypseg_new( i*deltat, pcpt,n_trj, d3r3, pan(k)%pnt(cycperm(4,1,q))%p%acfr, &
        & c1, c2, dc1, dc2, u1, u2, du1, du2)

    else

      if (pan(k)%neigh(q) /= j) then
        c1 = (1. - (d1r1 - floor(d1r1))) * pan(k)%dcirc(q,idtime1)
        c2 = (1. - (d1r2 - floor(d1r2))) * pan(k)%dcirc(q,idtime2)
      else
        c1 = (1. - (d1r1 - floor(d1r1))) * pan(k)%circ(idtime1)
        c2 = (1. - (d1r2 - floor(d1r2))) * pan(k)%circ(idtime2)
      end if

      d3r4 = velindhypseg_new( i*deltat, pcpt,n_trj, pan(k)%pnt(q)%p%acfr, &
        & pan(k)%pnt(cycperm(4,1,q))%p%acfr, (d1r1 - floor(d1r1)), (d1r2 - floor(d1r2)), dc1, dc2, u1, u2, du1, du2)

      cmat(:,k,j) = cmat(:,k,j) + d3r4
      if (pan(k)%neigh(q) > 0) then
        if (pan(k)%neigh(q) /= j) then
          cmat(:,pan(k)%neigh(q),j) = cmat(:,pan(k)%neigh(q),j) - d3r4
        end if
      end if

      pan(j)%vind = pan(j)%vind + velindhypseg_new( i*deltat, pcpt,n_trj, pan(k)%pnt(q)%p%acfr, &
        & pan(k)%pnt(cycperm(4,1,q))%p%acfr, c1, c2, dc1, dc2, u1, u2, du1, du2)

    end if

  end if

  end if

end subroutine wingwingcalcseg_2



subroutine wingwingcalcseg(pcpt,pan,cmat,vacfr,vfstream,uvind,duvind,npan,tsteps,deltat,i,j,k,q)
  use basictypes
  use pg
  use infl
  use constants
  use maths
  use hyperbolic
  implicit none

  type(panel),dimension(npan),intent(inout) :: pan
  real(kind=8),dimension(3,npan,npan),intent(inout) :: cmat
  real(kind=8),dimension(3),intent(in) :: pcpt,vacfr
  real(kind=8),dimension(3,npan,tsteps),intent(in) :: vfstream,uvind,duvind
  integer,intent(in) :: npan,tsteps,i,j,k,q
  real(kind=8),intent(in) :: deltat
  real(kind=8),dimension(3) :: d3r1,d3r2,d3r3,d3r4,d3r5,vfstrp,untrp,duntrp
  integer :: d1i1
  real(kind=8) :: d1r1,d1r2,d1r3,f,dcirc
  real(kind=8),parameter :: rzero = 0.
  real(kind=8),parameter :: runity = 1. 
  real(kind=8),dimension(3),parameter :: zerovec = 0.

  if (pan(k)%pnt(q)%p%ffl) then  
!    call ittodummy(pan(j)%cinfr(:,i),pan(k)%pnt(q)%p%set(i)%infr(1:3),vfstream,deltat,d3r1)
    d3r1 = pan(k)%pnt(q)%p%set(i)%infr(1:3)
  else
    d3r1 = intp(pan(k)%pnt(q)%p%set(pan(k)%pnt(q)%p%ti(1))%infr(1:3),pan(k)%pnt(q)%p%set(pan(k)%pnt(q)%p%ti(2))%infr(1:3), &
         & pan(k)%pnt(q)%p%tr(1),pan(k)%pnt(q)%p%tr(2))
    call ittozero(pcpt,d3r1,(-1.)*vacfr,d3r1)
!    d3r1 = pan(k)%pnt(q)%p%set(i)%infr(1:3)
  end if

  if (pan(k)%pnt(cycperm(4,1,q))%p%ffl) then
!    call ittodummy(pan(j)%cinfr(:,i),pan(k)%pnt(cycperm(4,1,q))%p%set(i)%infr(1:3),vfstream,deltat,d3r2)
    d3r2 = pan(k)%pnt(cycperm(4,1,q))%p%set(i)%infr(1:3)
  else
    d3r2 = intp(pan(k)%pnt(cycperm(4,1,q))%p%set(pan(k)%pnt(cycperm(4,1,q))%p%ti(1))%infr(1:3), &
         & pan(k)%pnt(cycperm(4,1,q))%p%set(pan(k)%pnt(cycperm(4,1,q))%p%ti(2))%infr(1:3), &
         & pan(k)%pnt(cycperm(4,1,q))%p%tr(1),pan(k)%pnt(cycperm(4,1,q))%p%tr(2))
    call ittozero(pcpt,d3r2,(-1.)*vacfr,d3r2)
!    d3r2 = pan(k)%pnt(cycperm(4,1,q))%p%set(i)%infr(1:3)
  end if

!need to interpolate panfstream, unind, duvind 

!d3r1 = pan(k)%pnt(q)%p%set(i)%infr(1:3)
!d3r2 = pan(k)%pnt(cycperm(4,1,q))%p%set(i)%infr(1:3)

  if (i > pan(k)%pnt(q)%p%ti(1)) then
    if (i > pan(k)%pnt(cycperm(4,1,q))%p%ti(1))then

      if (pan(k)%neigh(q) /= j) then
        d1r1 = intq(pan(k)%dcirc(q,pan(k)%pnt(q)%p%ti(1)),pan(k)%dcirc(q,pan(k)%pnt(q)%p%ti(2)), &
             & pan(k)%pnt(q)%p%tr(1),pan(k)%pnt(q)%p%tr(2))

        d1r2 = intq(pan(k)%dcirc(q,pan(k)%pnt(cycperm(4,1,q))%p%ti(1)),pan(k)%dcirc(q,pan(k)%pnt(cycperm(4,1,q))%p%ti(2)), &
             & pan(k)%pnt(cycperm(4,1,q))%p%tr(1),pan(k)%pnt(cycperm(4,1,q))%p%tr(2))
      else
        d1r1 = intq(pan(k)%circ(pan(k)%pnt(q)%p%ti(1)),pan(k)%circ(pan(k)%pnt(q)%p%ti(2)), &
             & pan(k)%pnt(q)%p%tr(1),pan(k)%pnt(q)%p%tr(2))

        d1r2 = intq(pan(k)%circ(pan(k)%pnt(cycperm(4,1,q))%p%ti(1)),pan(k)%circ(pan(k)%pnt(cycperm(4,1,q))%p%ti(2)), &
            & pan(k)%pnt(cycperm(4,1,q))%p%tr(1),pan(k)%pnt(cycperm(4,1,q))%p%tr(2))
      end if

      dcirc = 0.5*(pan(k)%circ(pan(k)%pnt(q)%p%ti(1))-pan(k)%circ(pan(k)%pnt(q)%p%ti(2)) + &
            & pan(k)%circ(pan(k)%pnt(cycperm(4,1,q))%p%ti(1))-pan(k)%circ(pan(k)%pnt(cycperm(4,1,q))%p%ti(2)))/deltat
      vfstrp = 0.5*(intp(vfstream(:,k,pan(k)%pnt(q)%p%ti(1)),vfstream(:,k,pan(k)%pnt(q)%p%ti(2)), &
             & pan(k)%pnt(q)%p%tr(1),pan(k)%pnt(q)%p%tr(2)) + &
             & intp(vfstream(:,k,pan(k)%pnt(cycperm(4,1,q))%p%ti(1)),vfstream(:,k,pan(k)%pnt(cycperm(4,1,q))%p%ti(2)), &
             & pan(k)%pnt(cycperm(4,1,q))%p%tr(1),pan(k)%pnt(cycperm(4,1,q))%p%tr(2)))
      untrp = 0.5*(intp(uvind(:,k,pan(k)%pnt(q)%p%ti(1)),uvind(:,k,pan(k)%pnt(q)%p%ti(2)), &
             & pan(k)%pnt(q)%p%tr(1),pan(k)%pnt(q)%p%tr(2)) + &
             & intp(uvind(:,k,pan(k)%pnt(cycperm(4,1,q))%p%ti(1)),uvind(:,k,pan(k)%pnt(cycperm(4,1,q))%p%ti(2)), &
             & pan(k)%pnt(cycperm(4,1,q))%p%tr(1),pan(k)%pnt(cycperm(4,1,q))%p%tr(2)))
      duntrp = 0.5*(intp(duvind(:,k,pan(k)%pnt(q)%p%ti(1)),duvind(:,k,pan(k)%pnt(q)%p%ti(2)), &
             & pan(k)%pnt(q)%p%tr(1),pan(k)%pnt(q)%p%tr(2)) + &
             & intp(duvind(:,k,pan(k)%pnt(cycperm(4,1,q))%p%ti(1)),duvind(:,k,pan(k)%pnt(cycperm(4,1,q))%p%ti(2)), &
             & pan(k)%pnt(cycperm(4,1,q))%p%tr(1),pan(k)%pnt(cycperm(4,1,q))%p%tr(2)))


      pan(j)%vind = pan(j)%vind + velindhypseg3D(pcpt,d3r1,d3r2,d1r1,d1r2,dcirc,vacfr,vfstrp,untrp,duntrp) !velindlinseg(pcpt,d3r1,d3r2,d1r1,d1r2)


    else if (i == pan(k)%pnt(cycperm(4,1,q))%p%ti(1)) then

      d1r1 = ints(runity,pan(k)%pnt(q)%p%tr(1),pan(k)%pnt(q)%p%tr(2))
      d1r2 = ints(runity,pan(k)%pnt(cycperm(4,1,q))%p%tr(1),pan(k)%pnt(cycperm(4,1,q))%p%tr(2))
      d1i1 = (i-1) - pan(k)%pnt(q)%p%ti(1) !if difference is greater than one timestep
      d1r3 = ints(runity,(1.-d1r1)+dble(d1i1),d1r2)
      d3r3 = d3r2 - d1r3*(d3r2-d3r1)

      dcirc = (pan(k)%circ(pan(k)%pnt(q)%p%ti(1))-pan(k)%circ(pan(k)%pnt(q)%p%ti(2)))/deltat
      vfstrp = intp(vfstream(:,k,pan(k)%pnt(q)%p%ti(1)),vfstream(:,k,pan(k)%pnt(q)%p%ti(2)), &
             & pan(k)%pnt(q)%p%tr(1),pan(k)%pnt(q)%p%tr(2))
      untrp = intp(uvind(:,k,pan(k)%pnt(q)%p%ti(1)),uvind(:,k,pan(k)%pnt(q)%p%ti(2)), &
             & pan(k)%pnt(q)%p%tr(1),pan(k)%pnt(q)%p%tr(2))
      duntrp = intp(duvind(:,k,pan(k)%pnt(q)%p%ti(1)),duvind(:,k,pan(k)%pnt(q)%p%ti(2)), &
             & pan(k)%pnt(q)%p%tr(1),pan(k)%pnt(q)%p%tr(2))


      d3r4 = velindhypseg3D(pcpt,d3r3,d3r2,rzero,d1r2,dcirc,vacfr,vfstrp,untrp,duntrp) !velindlinseg(pcpt,d3r3,d3r2,rzero,d1r2)
      cmat(:,k,j) = cmat(:,k,j) + d3r4
      if (pan(k)%neigh(q) > 0) then
        if (pan(k)%neigh(q) /= j) then
          cmat(:,pan(k)%neigh(q),j) = cmat(:,pan(k)%neigh(q),j) - d3r4
        end if
      end if

      !add first seg
      if (pan(k)%neigh(q) /= j) then
        d1r1 = pan(k)%dcirc(q,pan(k)%pnt(cycperm(4,1,q))%p%ti(2))
        d1r2 = (1.-d1r2)*pan(k)%dcirc(q,pan(k)%pnt(cycperm(4,1,q))%p%ti(2))
        d3r4 = velindhypseg3D(pcpt,d3r3,d3r2,d1r1,d1r2,dcirc,vacfr,vfstrp,untrp,duntrp) !velindlinseg(pcpt,d3r3,d3r2,d1r1,d1r2)

        pan(j)%vind = pan(j)%vind + d3r4

        !add free seg
        d1r2 = intq(pan(k)%dcirc(q,pan(k)%pnt(q)%p%ti(1)),pan(k)%dcirc(q,pan(k)%pnt(q)%p%ti(2)), &
             & pan(k)%pnt(q)%p%tr(1),pan(k)%pnt(q)%p%tr(2))

        d3r4 = velindhypseg3D(pcpt,d3r1,d3r3,d1r2,d1r1,dcirc,vacfr,vfstrp,untrp,duntrp) !velindlinseg(pcpt,d3r1,d3r3,d1r2,d1r1)

        pan(j)%vind = pan(j)%vind + d3r4
      else
        d1r1 = pan(k)%circ(pan(k)%pnt(cycperm(4,1,q))%p%ti(2))
        d1r2 = (1.-d1r2)*pan(k)%circ(pan(k)%pnt(cycperm(4,1,q))%p%ti(2))
        d3r4 = velindhypseg3D(pcpt,d3r3,d3r2,d1r1,d1r2,dcirc,vacfr,vfstrp,untrp,duntrp) !velindlinseg(pcpt,d3r3,d3r2,d1r1,d1r2)

        pan(j)%vind = pan(j)%vind + d3r4

        !add free seg
        d1r2 = intq(pan(k)%circ(pan(k)%pnt(q)%p%ti(1)),pan(k)%circ(pan(k)%pnt(q)%p%ti(2)), &
             & pan(k)%pnt(q)%p%tr(1),pan(k)%pnt(q)%p%tr(2))

        d3r4 = velindhypseg3D(pcpt,d3r1,d3r3,d1r2,d1r1,dcirc,vacfr,vfstrp,untrp,duntrp) !velindlinseg(pcpt,d3r1,d3r3,d1r2,d1r1)

        pan(j)%vind = pan(j)%vind + d3r4
      end if

    end if
  else if (i == pan(k)%pnt(q)%p%ti(1)) then

    if (i > pan(k)%pnt(cycperm(4,1,q))%p%ti(1))then

      d1r1 = ints(runity,pan(k)%pnt(q)%p%tr(1),pan(k)%pnt(q)%p%tr(2))
      d1r2 = ints(runity,pan(k)%pnt(cycperm(4,1,q))%p%tr(1),pan(k)%pnt(cycperm(4,1,q))%p%tr(2))
      d1i1 = (i-1) - pan(k)%pnt(cycperm(4,1,q))%p%ti(1)
      d1r3 = ints(runity,(1.-d1r2)+dble(d1i1),d1r1)

      d3r3 = d3r1 + d1r3*(d3r2-d3r1)

      dcirc = (pan(k)%circ(pan(k)%pnt(cycperm(4,1,q))%p%ti(1))-pan(k)%circ(pan(k)%pnt(cycperm(4,1,q))%p%ti(2)))/deltat
      vfstrp = intp(vfstream(:,k,pan(k)%pnt(cycperm(4,1,q))%p%ti(1)),vfstream(:,k,pan(k)%pnt(cycperm(4,1,q))%p%ti(2)), &
             & pan(k)%pnt(cycperm(4,1,q))%p%tr(1),pan(k)%pnt(cycperm(4,1,q))%p%tr(2))
      untrp = intp(uvind(:,k,pan(k)%pnt(cycperm(4,1,q))%p%ti(1)),uvind(:,k,pan(k)%pnt(cycperm(4,1,q))%p%ti(2)), &
             & pan(k)%pnt(cycperm(4,1,q))%p%tr(1),pan(k)%pnt(cycperm(4,1,q))%p%tr(2))
      duntrp = intp(duvind(:,k,pan(k)%pnt(cycperm(4,1,q))%p%ti(1)),duvind(:,k,pan(k)%pnt(cycperm(4,1,q))%p%ti(2)), &
             & pan(k)%pnt(cycperm(4,1,q))%p%tr(1),pan(k)%pnt(cycperm(4,1,q))%p%tr(2))

      d3r4 = velindhypseg3D(pcpt,d3r1,d3r3,d1r1,rzero,dcirc,vacfr,vfstrp,untrp,duntrp) !velindlinseg(pcpt,d3r1,d3r3,d1r1,rzero)
      cmat(:,k,j) = cmat(:,k,j) + d3r4
      if (pan(k)%neigh(q) > 0) then
        if (pan(k)%neigh(q) /= j) then
          cmat(:,pan(k)%neigh(q),j) = cmat(:,pan(k)%neigh(q),j) - d3r4
        end if
      end if

      if (pan(k)%neigh(q) /= j) then
        !add first seg
        d1r1 = (1.-d1r1)*pan(k)%dcirc(q,pan(k)%pnt(q)%p%ti(2))
        d1r2 = pan(k)%dcirc(q,pan(k)%pnt(q)%p%ti(2))
        d3r4 = velindhypseg3D(pcpt,d3r1,d3r3,d1r1,d1r2,dcirc,vacfr,vfstrp,untrp,duntrp) !velindlinseg(pcpt,d3r1,d3r3,d1r1,d1r2)

        pan(j)%vind = pan(j)%vind + d3r4

        !add free seg
        d1r1 = intq(pan(k)%dcirc(q,pan(k)%pnt(cycperm(4,1,q))%p%ti(1)),pan(k)%dcirc(q,pan(k)%pnt(cycperm(4,1,q))%p%ti(2)), &
             & pan(k)%pnt(cycperm(4,1,q))%p%tr(1),pan(k)%pnt(cycperm(4,1,q))%p%tr(2))
        d3r4 = velindhypseg3D(pcpt,d3r3,d3r2,d1r2,d1r1,dcirc,vacfr,vfstrp,untrp,duntrp) !velindlinseg(pcpt,d3r3,d3r2,d1r2,d1r1)

        pan(j)%vind = pan(j)%vind + d3r4
      else
        !add first seg
        d1r1 = (1.-d1r1)*pan(k)%circ(pan(k)%pnt(q)%p%ti(2))
        d1r2 = pan(k)%circ(pan(k)%pnt(q)%p%ti(2))
        d3r4 = velindhypseg3D(pcpt,d3r1,d3r3,d1r1,d1r2,dcirc,vacfr,vfstrp,untrp,duntrp) !velindlinseg(pcpt,d3r1,d3r3,d1r1,d1r2)

        pan(j)%vind = pan(j)%vind + d3r4

        !add free seg
        d1r1 = intq(pan(k)%circ(pan(k)%pnt(cycperm(4,1,q))%p%ti(1)),pan(k)%circ(pan(k)%pnt(cycperm(4,1,q))%p%ti(2)), &
             & pan(k)%pnt(cycperm(4,1,q))%p%tr(1),pan(k)%pnt(cycperm(4,1,q))%p%tr(2))
        d3r4 = velindhypseg3D(pcpt,d3r3,d3r2,d1r2,d1r1,dcirc,vacfr,vfstrp,untrp,duntrp) !velindlinseg(pcpt,d3r3,d3r2,d1r2,d1r1)

        pan(j)%vind = pan(j)%vind + d3r4
      end if

    else if (i == pan(k)%pnt(cycperm(4,1,q))%p%ti(1)) then

      d1r1 = ints(runity,pan(k)%pnt(q)%p%tr(1),pan(k)%pnt(q)%p%tr(2))
      d1r2 = ints(runity,pan(k)%pnt(cycperm(4,1,q))%p%tr(1),pan(k)%pnt(cycperm(4,1,q))%p%tr(2))

      vfstrp = vfstream(:,k,i)
      if (i>1) then
        untrp = uvind(:,k,i-1)
        duntrp = duvind(:,k,i-1)
      else 
        untrp = 0.
        duntrp = 0.
      end if

      dcirc = 0.
      d3r4 = velindhypseg3D(pcpt,d3r1,d3r2,d1r1,d1r2,dcirc,vacfr,vfstrp,untrp,duntrp) !velindlinseg(pcpt,d3r1,d3r2,d1r1,d1r2)

      cmat(:,k,j) = cmat(:,k,j) + d3r4
      if (pan(k)%neigh(q) > 0) then
        if (pan(k)%neigh(q) /= j) then
          cmat(:,pan(k)%neigh(q),j) = cmat(:,pan(k)%neigh(q),j) - d3r4
        end if
      end if

      if (pan(k)%neigh(q) /= j) then
        d1r1 = (1-d1r1)*pan(k)%dcirc(q,pan(k)%pnt(q)%p%ti(2))
        d1r2 = (1-d1r2)*pan(k)%dcirc(q,pan(k)%pnt(cycperm(4,1,q))%p%ti(2))
      else
        d1r1 = (1-d1r1)*pan(k)%circ(pan(k)%pnt(q)%p%ti(2))
        d1r2 = (1-d1r2)*pan(k)%circ(pan(k)%pnt(cycperm(4,1,q))%p%ti(2))
      end if

      if (i>2) then
        dcirc = 0.5*(pan(k)%circ(pan(k)%pnt(q)%p%ti(2))-pan(k)%circ(pan(k)%pnt(q)%p%ti(2)-1) + &
              & pan(k)%circ(pan(k)%pnt(cycperm(4,1,q))%p%ti(2))-pan(k)%circ(pan(k)%pnt(cycperm(4,1,q))%p%ti(2)-1))/deltat
      else 
        dcirc = 0.
      end if

      if (d1r1 /= 0. .or. d1r2 /= 0.) then
        d3r4 = velindhypseg3D(pcpt,d3r1,d3r2,d1r1,d1r2,dcirc,vacfr,vfstrp,untrp,duntrp) !velindlinseg(pcpt,d3r1,d3r2,d1r1,d1r2)
        pan(j)%vind = pan(j)%vind + d3r4
      end if

    end if
  end if

end subroutine wingwingcalcseg


subroutine wakewingcalcseg_2(wingpan,pin,cmat,wakesrcpan,wakepan,npan,nwpanels,nwakepan,d1i1,deltat,i,j,k,l,q)
  use basictypes
  use pg
  use infl
  use constants
  use maths
  use hyperbolic
  use penetration
  use input, only: inner, avps
  implicit none

  real(kind=8),dimension(3),intent(in) :: pin
  type(panel),intent(inout),target :: wingpan
  real(kind=8),dimension(:),pointer :: vind
  real(kind=8),dimension(3,npan,npan),intent(inout) :: cmat
  integer,dimension(nwakepan,2),intent(in) :: wakesrcpan
  type(panel),dimension((nwpanels+2)*nwakepan) :: wakepan
  integer,intent(in) :: npan,nwpanels,nwakepan,d1i1,i,j,k,l,q
  real(kind=8),intent(in) :: deltat
  integer :: d1i2
  real(kind=8),dimension(3) :: d3r1,d3r2,d3r3,d3r4
  real(kind=8) :: d1r1,d1r2,d1r3
  real(kind=8),dimension(3),parameter :: zerovec = 0.
  real(kind=8),parameter :: rzero = 0.
  real(kind=8),parameter :: runity = 1.
  real(kind=8),dimension(3) :: norm

  vind => wingpan%vind

  if ((wakepan(d1i1)%pnt(q)%p%ffl .eqv. .false.) .and. (wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ffl .eqv. .false.)) then  

  if (wakepan(d1i1)%pnt(q)%p%ti(2) >= 1 .and. wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2) >= 1) then
    d3r1 = intp(wakepan(d1i1)%pnt(q)%p%set(wakepan(d1i1)%pnt(q)%p%ti(1))%infr(1:3), &
         & wakepan(d1i1)%pnt(q)%p%set(wakepan(d1i1)%pnt(q)%p%ti(2))%infr(1:3), &
         & wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))

    d3r2 = intp(wakepan(d1i1)%pnt(cycperm(4,1,q))%p%set(wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1))%infr(1:3), &
         & wakepan(d1i1)%pnt(cycperm(4,1,q))%p%set(wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2))%infr(1:3), &
         & wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))
          
    norm = cross(wingpan%pnt(3)%p%set(i)%infr(1:3)-wingpan%pnt(1)%p%set(i)%infr(1:3), &
         & wingpan%pnt(2)%p%set(i)%infr(1:3)-wingpan%pnt(4)%p%set(i)%infr(1:3))

    !normalise to unit vector
    norm = norm/modulus(norm)   

  if (chkpppen(wingpan,i,d3r1,d3r2) .eqv. .false.) then
  if (modulus(d3r1-pin) >  inner * avps .or. (abs(dot_product(norm,d3r1-pin)) > inner * avps &
    & .and. abs(dot_product(norm,d3r2-pin)) > 0.3 * avps)) then


  if (l>2) then

    d1r1 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(1)),wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2)), &
         & wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))
    d1r2 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)), &
         & wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2)), &
         & wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))
    vind = vind + velindlinseg(pin,d3r1,d3r2,d1r1,d1r2)

  else if (l==2) then
    if ((wakepan(d1i1)%neigh(q) > 0) .and. (wakepan(d1i1)%neigh(q) <= nwakepan)) then
      if (i > wakepan(d1i1)%pnt(q)%p%ti(1)) then
        if (i > wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)) then


          d1r1 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(1)),wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2)), &
               & wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))
          d1r2 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)), &
               & wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2)), &
               & wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))
          vind = vind + velindlinseg(pin,d3r1,d3r2,d1r1,d1r2)

        else if (i == wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)) then

          d1r1 = ints(runity,wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))
          d1r2 = ints(runity,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))
          d1i2 = (i-1) - wakepan(d1i1)%pnt(q)%p%ti(1) !if difference is greater than one timestep
          d1r3 = ints(runity,(1.-d1r1)+dble(d1i2),d1r2)
          d3r3 = d3r2 - d1r3*(d3r2-d3r1)

          d3r4 = velindlinseg(pin,d3r3,d3r2,rzero,d1r2)
          cmat(:,wakesrcpan(k,1),j) = cmat(:,wakesrcpan(k,1),j) - d3r4 !-ve because in second line, first line defined other way around
          vind = vind + d3r4*wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1))

          !add first seg
          d1r1 = wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2))
          d1r2 = (1.-d1r2)*wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2))
          d3r4 = velindlinseg(pin,d3r3,d3r2,d1r1,d1r2)
          vind = vind + d3r4

          !add free seg
          d1r2 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(1)),wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2)), &
               & wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))
          d3r4 = velindlinseg(pin,d3r1,d3r3,d1r2,d1r1)
          vind = vind + d3r4

        end if
      else if (i == wakepan(d1i1)%pnt(q)%p%ti(1)) then
        if (i > wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)) then


          d1r1 = ints(runity,wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))
          d1r2 = ints(runity,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))
          d1i2 = (i-1) - wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)
          d1r3 = ints(runity,(1.-d1r2)+dble(d1i2),d1r1)
  
          d3r3 = d3r1 + d1r3*(d3r2-d3r1)

          d3r4 = velindlinseg(pin,d3r1,d3r3,d1r1,rzero)
          cmat(:,wakesrcpan(k,1),j) = cmat(:,wakesrcpan(k,1),j) - d3r4
          vind = vind + d3r4*wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(1))

          !add first seg
          d1r1 = (1.-d1r1)*wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2))
          d1r2 = wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2))

          d3r4 = velindlinseg(pin,d3r1,d3r3,d1r1,d1r2)
  
          vind = vind + d3r4

          !add free seg
          d1r1 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)), &
               & wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2)), &
               & wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))
          d3r4 = velindlinseg(pin,d3r3,d3r2,d1r2,d1r1)
          vind = vind + d3r4

        else if (i == wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)) then

          d1r1 = ints(runity,wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))
          d1r2 = ints(runity,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))

          d3r4 = velindlinseg(pin,d3r1,d3r2,d1r1,d1r2)
          cmat(:,wakesrcpan(k,1),j) = cmat(:,wakesrcpan(k,1),j) - d3r4
          vind = vind + d3r4*wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(1))

          d1r1 = (1-d1r1)*wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2))
          d1r2 = (1-d1r2)*wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2))
          d3r4 = velindlinseg(pin,d3r1,d3r2,d1r1,d1r2)
          vind = vind + d3r4


        end if
      end if

    else

      d1r1 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(1)),wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2)), &
           & wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))
      d1r2 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)), &
           & wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2)), &
           & wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))
      vind = vind + velindlinseg(pin,d3r1,d3r2,d1r1,d1r2)

    end if

  else if (l==1) then

    if (i > wakepan(d1i1)%pnt(q)%p%ti(1)) then
      if (i > wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)) then

        d1r1 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(1)),wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2)), &
             & wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))
        d1r2 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)), &
             & wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2)), &
             & wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))
        vind = vind + velindlinseg(pin,d3r1,d3r2,d1r1,d1r2)

      else if (i == wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)) then

        d1r1 = ints(runity,wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))
        d1r2 = ints(runity,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))
        d1i2 = (i-1) - wakepan(d1i1)%pnt(q)%p%ti(1) !if difference is greater than one timestep
        d1r3 = ints(runity,(1.-d1r1)+dble(d1i2),d1r2)
        d3r3 = d3r2 - d1r3*(d3r2-d3r1)

        d3r4 = velindlinseg(pin,d3r3,d3r2,rzero,d1r2)

        cmat(:,wakesrcpan(k,1),j) = cmat(:,wakesrcpan(k,1),j) + d3r4
        if (wakepan(d1i1)%neigh(q) > 0 .and. wakepan(d1i1)%neigh(q) <= nwakepan) then
          cmat(:,wakesrcpan(wakepan(d1i1)%neigh(q),1),j) = cmat(:,wakesrcpan(wakepan(d1i1)%neigh(q),1),j) - d3r4
        else if (wakepan(d1i1)%neigh(q) > nwakepan) then
          vind = vind - d3r4*wakepan(wakepan(d1i1)%neigh(q))%dcirc(cycperm(4,2,q),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1))
        end if

        !add first seg
        d1r1 = wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2))
        d1r2 = (1.-d1r2)*wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2))
        d3r4 = velindlinseg(pin,d3r3,d3r2,d1r1,d1r2)
        vind = vind + d3r4

        !add free seg
        d1r2 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(1)),wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2)), &
             & wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))
        d3r4 = velindlinseg(pin,d3r1,d3r3,d1r2,d1r1)
        vind = vind + d3r4

      end if
    else if (i == wakepan(d1i1)%pnt(q)%p%ti(1)) then
      if (i > wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)) then


        d1r1 = ints(runity,wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))
        d1r2 = ints(runity,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))
        d1i2 = (i-1) - wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)
        d1r3 = ints(runity,(1.-d1r2)+dble(d1i2),d1r1)
  
        d3r3 = d3r1 + d1r3*(d3r2-d3r1)

!        print*
!        print*, wakepan(d1i1)%pnt(q)%p%ti(1), wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)
!        print*, wakepan(d1i1)%pnt(q)%p%ti(2), wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2)
!        print*, d1r1,d1r2,d1r3
!        print*, sqrt(dot_product(d3r1-d3r3,d3r1-d3r3))

        d3r4 = velindlinseg(pin,d3r1,d3r3,d1r1,rzero)
        cmat(:,wakesrcpan(k,1),j) = cmat(:,wakesrcpan(k,1),j) + d3r4
        if (wakepan(d1i1)%neigh(q) > 0 .and. wakepan(d1i1)%neigh(q) <= nwakepan) then
          cmat(:,wakesrcpan(wakepan(d1i1)%neigh(q),1),j) = cmat(:,wakesrcpan(wakepan(d1i1)%neigh(q),1),j) - d3r4
        else if (wakepan(d1i1)%neigh(q) > nwakepan) then
          vind = vind - d3r4*wakepan(wakepan(d1i1)%neigh(q))%dcirc(cycperm(4,2,q),wakepan(d1i1)%pnt(q)%p%ti(1))
        end if

        !add first seg
        d1r1 = (1.-d1r1)*wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2))
        d1r2 = wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2))
        d3r4 = velindlinseg(pin,d3r1,d3r3,d1r1,d1r2)
        vind = vind + d3r4

        !add free seg
        d1r1 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)), &
             & wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2)), &
             & wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))
        d3r4 = velindlinseg(pin,d3r3,d3r2,d1r2,d1r1)
        vind = vind + d3r4

      else if (i == wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)) then


        d1r1 = ints(runity,wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))
        d1r2 = ints(runity,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))

        d3r4 = velindlinseg(pin,d3r1,d3r2,d1r1,d1r2)
        cmat(:,wakesrcpan(k,1),j) = cmat(:,wakesrcpan(k,1),j) + d3r4
        if (wakepan(d1i1)%neigh(q) > 0 .and. wakepan(d1i1)%neigh(q) <= nwakepan) then
          cmat(:,wakesrcpan(wakepan(d1i1)%neigh(q),1),j) = cmat(:,wakesrcpan(wakepan(d1i1)%neigh(q),1),j) - d3r4
        else if (wakepan(d1i1)%neigh(q) > nwakepan) then
          vind = vind - d3r4*wakepan(wakepan(d1i1)%neigh(q))%dcirc(cycperm(4,2,q),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1))
        end if

        d1r1 = (1-d1r1)*wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2))
        d1r2 = (1-d1r2)*wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2))
        d3r4 = velindlinseg(pin,d3r1,d3r2,d1r1,d1r2)
        vind = vind + d3r4


      end if
    end if
  end if
  end if
  end if

  end if

  end if

end subroutine wakewingcalcseg_2




subroutine wakewingcalcseg(pin,vind,cmat,vacfr,wakesrcpan,wakepan,npan,nwpanels,nwakepan,d1i1,deltat,i,j,k,l,q)
  use basictypes
  use pg
  use infl
  use constants
  use maths
  use hyperbolic
  implicit none

  real(kind=8),dimension(3),intent(in) :: pin
  real(kind=8),dimension(3),intent(inout) :: vind
  real(kind=8),dimension(3,npan,npan),intent(inout) :: cmat
  real(kind=8),dimension(3),intent(in) :: vacfr
  integer,dimension(nwakepan,2),intent(in) :: wakesrcpan
  type(panel),dimension((nwpanels+2)*nwakepan) :: wakepan
  integer,intent(in) :: npan,nwpanels,nwakepan,d1i1,i,j,k,l,q
  real(kind=8),intent(in) :: deltat
  integer :: d1i2
  real(kind=8),dimension(3) :: d3r1,d3r2,d3r3,d3r4
  real(kind=8) :: d1r1,d1r2,d1r3
  real(kind=8),dimension(3),parameter :: zerovec = 0.
  real(kind=8),parameter :: rzero = 0.
  real(kind=8),parameter :: runity = 1.


  if (wakepan(d1i1)%pnt(q)%p%ffl) then  
!    call ittodummy(pan(j)%cinfr(:,i),pan(k)%pnt(q)%p%set(i)%infr(1:3),vacfr,deltat,d3r1)
    d3r1 = wakepan(d1i1)%pnt(q)%p%set(i)%infr(1:3)
  else
    d3r1 = intp(wakepan(d1i1)%pnt(q)%p%set(wakepan(d1i1)%pnt(q)%p%ti(1))%infr(1:3), &
         & wakepan(d1i1)%pnt(q)%p%set(wakepan(d1i1)%pnt(q)%p%ti(2))%infr(1:3), &
         & wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))
    call ittozero(pin,d3r1,(-1.)*vacfr,d3r1)
!    d3r1 = wakepan(d1i1)%pnt(q)%p%set(i)%infr(1:3)
  end if

  if (wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ffl) then
!    call ittodummy(pan(j)%cinfr(:,i),pan(k)%pnt(cycperm(4,1,q))%p%set(i)%infr(1:3),vacfr,deltat,d3r2)
    d3r2 = wakepan(d1i1)%pnt(cycperm(4,1,q))%p%set(i)%infr(1:3)
  else
    d3r2 = intp(wakepan(d1i1)%pnt(cycperm(4,1,q))%p%set(wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1))%infr(1:3), &
         & wakepan(d1i1)%pnt(cycperm(4,1,q))%p%set(wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2))%infr(1:3), &
         & wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))
    call ittozero(pin,d3r2,(-1.)*vacfr,d3r2)
!    d3r2 = wakepan(d1i1)%pnt(cycperm(4,1,q))%p%set(i)%infr(1:3)
  end if


  if (l>2) then

    d1r1 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(1)),wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2)), &
         & wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))
    d1r2 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)), &
         & wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2)), &
         & wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))
    vind = vind + velindhypseg3D_wake(pin,d3r1,d3r2,d1r1,d1r2,vacfr,zerovec)

  else if (l==2) then
    if ((wakepan(d1i1)%neigh(q) > 0) .and. (wakepan(d1i1)%neigh(q) <= nwakepan)) then
      if (i > wakepan(d1i1)%pnt(q)%p%ti(1)) then
        if (i > wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)) then


          d1r1 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(1)),wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2)), &
               & wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))
          d1r2 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)), &
               & wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2)), &
               & wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))
          vind = vind + velindhypseg3D_wake(pin,d3r1,d3r2,d1r1,d1r2,vacfr,zerovec)

        else if (i == wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)) then

          d1r1 = ints(runity,wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))
          d1r2 = ints(runity,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))
          d1i2 = (i-1) - wakepan(d1i1)%pnt(q)%p%ti(1) !if difference is greater than one timestep
          d1r3 = ints(runity,(1.-d1r1)+dble(d1i2),d1r2)
          d3r3 = d3r2 - d1r3*(d3r2-d3r1)

          d3r4 = velindhypseg3D_wake(pin,d3r3,d3r2,rzero,d1r2,vacfr,zerovec)
          cmat(:,wakesrcpan(k,1),j) = cmat(:,wakesrcpan(k,1),j) - d3r4 !-ve because in second line, first line defined other way around
          vind = vind + d3r4*wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1))

          !add first seg
          d1r1 = wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2))
          d1r2 = (1.-d1r2)*wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2))
          d3r4 = velindhypseg3D_wake(pin,d3r3,d3r2,d1r1,d1r2,vacfr,zerovec)
          vind = vind + d3r4

          !add free seg
          d1r2 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(1)),wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2)), &
               & wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))
          d3r4 = velindhypseg3D_wake(pin,d3r1,d3r3,d1r2,d1r1,vacfr,zerovec)
          vind = vind + d3r4

        end if
      else if (i == wakepan(d1i1)%pnt(q)%p%ti(1)) then
        if (i > wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)) then


          d1r1 = ints(runity,wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))
          d1r2 = ints(runity,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))
          d1i2 = (i-1) - wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)
          d1r3 = ints(runity,(1.-d1r2)+dble(d1i2),d1r1)
  
          d3r3 = d3r1 + d1r3*(d3r2-d3r1)

          d3r4 = velindhypseg3D_wake(pin,d3r1,d3r3,d1r1,rzero,vacfr,zerovec)
          cmat(:,wakesrcpan(k,1),j) = cmat(:,wakesrcpan(k,1),j) - d3r4
          vind = vind + d3r4*wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(1))

          !add first seg
          d1r1 = (1.-d1r1)*wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2))
          d1r2 = wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2))

          d3r4 = velindhypseg3D_wake(pin,d3r1,d3r3,d1r1,d1r2,vacfr,zerovec)
  
          vind = vind + d3r4

          !add free seg
          d1r1 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)), &
               & wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2)), &
               & wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))
          d3r4 = velindhypseg3D_wake(pin,d3r3,d3r2,d1r2,d1r1,vacfr,zerovec)
          vind = vind + d3r4

        else if (i == wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)) then

          d1r1 = ints(runity,wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))
          d1r2 = ints(runity,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))

          d3r4 = velindhypseg3D_wake(pin,d3r1,d3r2,d1r1,d1r2,vacfr,zerovec)
          cmat(:,wakesrcpan(k,1),j) = cmat(:,wakesrcpan(k,1),j) - d3r4
          vind = vind + d3r4*wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(1))

          d1r1 = (1-d1r1)*wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2))
          d1r2 = (1-d1r2)*wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2))
          d3r4 = velindhypseg3D_wake(pin,d3r1,d3r2,d1r1,d1r2,vacfr,zerovec)
          vind = vind + d3r4


        end if
      end if

    else

      d1r1 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(1)),wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2)), &
           & wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))
      d1r2 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)), &
           & wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2)), &
           & wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))
      vind = vind + velindhypseg3D_wake(pin,d3r1,d3r2,d1r1,d1r2,vacfr,zerovec)

    end if

  else if (l==1) then

    if (i > wakepan(d1i1)%pnt(q)%p%ti(1)) then
      if (i > wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)) then

        d1r1 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(1)),wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2)), &
             & wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))
        d1r2 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)), &
             & wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2)), &
             & wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))
        vind = vind + velindhypseg3D_wake(pin,d3r1,d3r2,d1r1,d1r2,vacfr,zerovec)

      else if (i == wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)) then

        d1r1 = ints(runity,wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))
        d1r2 = ints(runity,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))
        d1i2 = (i-1) - wakepan(d1i1)%pnt(q)%p%ti(1) !if difference is greater than one timestep
        d1r3 = ints(runity,(1.-d1r1)+dble(d1i2),d1r2)
        d3r3 = d3r2 - d1r3*(d3r2-d3r1)

        d3r4 = velindhypseg3D_wake(pin,d3r3,d3r2,rzero,d1r2,vacfr,zerovec)

        cmat(:,wakesrcpan(k,1),j) = cmat(:,wakesrcpan(k,1),j) + d3r4
        if (wakepan(d1i1)%neigh(q) > 0 .and. wakepan(d1i1)%neigh(q) <= nwakepan) then
          cmat(:,wakesrcpan(wakepan(d1i1)%neigh(q),1),j) = cmat(:,wakesrcpan(wakepan(d1i1)%neigh(q),1),j) - d3r4
        else if (wakepan(d1i1)%neigh(q) > nwakepan) then
          vind = vind - d3r4*wakepan(wakepan(d1i1)%neigh(q))%dcirc(cycperm(4,2,q),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1))
        end if

        !add first seg
        d1r1 = wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2))
        d1r2 = (1.-d1r2)*wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2))
        d3r4 = velindhypseg3D_wake(pin,d3r3,d3r2,d1r1,d1r2,vacfr,zerovec)
        vind = vind + d3r4

        !add free seg
        d1r2 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(1)),wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2)), &
             & wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))
        d3r4 = velindhypseg3D_wake(pin,d3r1,d3r3,d1r2,d1r1,vacfr,zerovec)
        vind = vind + d3r4

      end if
    else if (i == wakepan(d1i1)%pnt(q)%p%ti(1)) then
      if (i > wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)) then


        d1r1 = ints(runity,wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))
        d1r2 = ints(runity,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))
        d1i2 = (i-1) - wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)
        d1r3 = ints(runity,(1.-d1r2)+dble(d1i2),d1r1)
  
        d3r3 = d3r1 + d1r3*(d3r2-d3r1)

        d3r4 = velindhypseg3D_wake(pin,d3r1,d3r3,d1r1,rzero,vacfr,zerovec)
        cmat(:,wakesrcpan(k,1),j) = cmat(:,wakesrcpan(k,1),j) + d3r4
        if (wakepan(d1i1)%neigh(q) > 0 .and. wakepan(d1i1)%neigh(q) <= nwakepan) then
          cmat(:,wakesrcpan(wakepan(d1i1)%neigh(q),1),j) = cmat(:,wakesrcpan(wakepan(d1i1)%neigh(q),1),j) - d3r4
        else if (wakepan(d1i1)%neigh(q) > nwakepan) then
          vind = vind - d3r4*wakepan(wakepan(d1i1)%neigh(q))%dcirc(cycperm(4,2,q),wakepan(d1i1)%pnt(q)%p%ti(1))
        end if

        !add first seg
        d1r1 = (1.-d1r1)*wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2))
        d1r2 = wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2))
        d3r4 = velindhypseg3D_wake(pin,d3r1,d3r3,d1r1,d1r2,vacfr,zerovec)
        vind = vind + d3r4

        !add free seg
        d1r1 = intq(wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)), &
             & wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2)), &
             & wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))
        d3r4 = velindhypseg3D_wake(pin,d3r3,d3r2,d1r2,d1r1,vacfr,zerovec)
        vind = vind + d3r4

      else if (i == wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1)) then


        d1r1 = ints(runity,wakepan(d1i1)%pnt(q)%p%tr(1),wakepan(d1i1)%pnt(q)%p%tr(2))
        d1r2 = ints(runity,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(1),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%tr(2))

        d3r4 = velindhypseg3D_wake(pin,d3r1,d3r2,d1r1,d1r2,vacfr,zerovec)
        cmat(:,wakesrcpan(k,1),j) = cmat(:,wakesrcpan(k,1),j) + d3r4
        if (wakepan(d1i1)%neigh(q) > 0 .and. wakepan(d1i1)%neigh(q) <= nwakepan) then
          cmat(:,wakesrcpan(wakepan(d1i1)%neigh(q),1),j) = cmat(:,wakesrcpan(wakepan(d1i1)%neigh(q),1),j) - d3r4
        else if (wakepan(d1i1)%neigh(q) > nwakepan) then
          vind = vind - d3r4*wakepan(wakepan(d1i1)%neigh(q))%dcirc(cycperm(4,2,q),wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(1))
        end if

        d1r1 = (1-d1r1)*wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(q)%p%ti(2))
        d1r2 = (1-d1r2)*wakepan(d1i1)%dcirc(q,wakepan(d1i1)%pnt(cycperm(4,1,q))%p%ti(2))
        d3r4 = velindhypseg3D_wake(pin,d3r1,d3r2,d1r1,d1r2,vacfr,zerovec)
        vind = vind + d3r4


      end if
    end if
  end if

end subroutine wakewingcalcseg

!###################################################################################################################################
!###################################################################################################################################

!subroutine ittodummy(pcpt,pemit,uinf,deltat,pout)
!  use pg, only:k2,intp
!  use constants
!  implicit none

!  real(kind=8),intent(in) :: deltat
!  real(kind=8),dimension(3),intent(in) :: pcpt,pemit,uinf
!  real(kind=8),dimension(3),intent(out) :: pout
!  integer :: l,tmin,tmax,tmid
!  real(kind=8) :: d1r1,d1r2

!  tmin = 0
!  tmax = 10000
!  tmid = tmin + int((tmax-tmin+1)/2)

!  d1r1 = k2(a,deltat,tmin,pemit-(pcpt-tmin*deltat*uinf))
!  d1r2 = k2(a,deltat,tmax,pemit-(pcpt-tmax*deltat*uinf))

!  if (d1r1 > 0. .and. d1r2 < 0.) then
!    do
!      d1r1 = k2(a,deltat,tmin,pemit-(pcpt-tmin*deltat*uinf))
!      d1r2 = k2(a,deltat,tmid,pemit-(pcpt-tmid*deltat*uinf))
!      if (tmin+1 == tmax) then
!        pout = intp((pemit+tmin*deltat*uinf),(pemit+tmax*deltat*uinf),d1r1,d1r2)
!        exit
!      end if
!      if (d1r1 >= 0. .and. d1r2 < 0.) then
!        tmax = tmid
!        tmid = tmin + int((tmax-tmin+1)/2)
!      else if (d1r1 > 0. .and. d1r2 >= 0.) then
!        tmin = tmid
!        tmid = tmin + int((tmax-tmin+1)/2)
!      end if
!    end do
!  else if (d1r1 <= 0.) then
!    pout = pcpt
!  else if (d1r2 >= 0.) then
!    print*, 'choose larger tmax in ittodummy'
!    pout = (pemit+tmax*deltat*uinf)
!  end if

!!pout = pemit

!!  d1r1 = k2(a,deltat,0,pemit-pcpt)
!!  if (d1r1 > 0.) then
!!    do l=1,10000,1
!!      d1r2 = k2(a,deltat,l,pemit-(pcpt-l*deltat*uinf))
!!!print*, dummy1r2,l
!!      if (d1r2 > 0.) then
!!        d1r1 = d1r2
!!      else
!!!        pout = intp((pcpt-(l-1)*deltat*uinf),(pcpt-l*deltat*uinf),dummy1r1,dummy1r2)
!!        pout = intp((pemit+(l-1)*deltat*uinf),(pemit+l*deltat*uinf),d1r1,d1r2)
!!        exit
!!!print*, 'here'
!!      end if
!!      if (l == 10000) then
!!        print*, 'ittodummy not exiting'
!!!        pout = (pcpt-l*deltat*uinf)
!!        pout = (pemit+l*deltat*uinf)
!!      end if
!!    end do

!!  else
!!    pout = pcpt
!!  end if
!!!    pout = pemit
!end subroutine ittodummy


subroutine ittodummy(pcpt,pemit,uinf,deltat,pout)
  use constants
  implicit none

  real(kind=8),intent(in) :: deltat
  real(kind=8),dimension(3),intent(in) :: pcpt,pemit,uinf
  real(kind=8),dimension(3),intent(out) :: pout
  integer :: l,tmin,tmax,tmid
  real(kind=8) :: d1r1,d1r2,M,beta
!  real(kind=8),dimension(3) :: d3r1,d3r2

!!rotation in here
      M = sqrt(dot_product(uinf,uinf))/a
      beta = sqrt(1.-M**2)

      pout(1) = pcpt(1)+((pemit(1)-pcpt(1)) + &
        & M*sqrt((pemit(1)-pcpt(1))**2+beta**2*((pemit(2)-pcpt(2))**2+(pemit(3)-pcpt(3))**2)))/beta**2
!      pout(1) = pemit(1) +(pemit(1)+M*sqrt(pemit(1)**2+beta**2*(pemit(2)**2+pemit(3)**2)))/beta**2
      pout(2) = pemit(2)
      pout(3) = pemit(3)

!  d1r1 = (dot_product(pemit-pcpt,uinf) + sqrt(dot_product(pemit-pcpt,uinf)**2 + &
!       & (a**2 - dot_product(uinf,uinf))*dot_product(pemit-pcpt,pemit-pcpt)))/(a**2 - dot_product(uinf,uinf))
!  d1r2 = (dot_product(pemit-pcpt,uinf) - sqrt(dot_product(pemit-pcpt,uinf)**2 + &
!       & (a**2 - dot_product(uinf,uinf))*dot_product(pemit-pcpt,pemit-pcpt)))/(a**2 - dot_product(uinf,uinf))

!  d1r1 = ((pemit(1)-pcpt(1))*uinf(1) + sqrt(((pemit(1)-pcpt(1))*uinf(1))**2 + &
!       & (a**2 - dot_product(uinf,uinf))*(pemit(1)-pcpt(1))**2))/(a**2 - dot_product(uinf,uinf))
!  d1r2 = ((pemit(1)-pcpt(1))*uinf(1) - sqrt(((pemit(1)-pcpt(1))*uinf(1))**2 + &
!       & (a**2 - dot_product(uinf,uinf))*(pemit(1)-pcpt(1))**2))/(a**2 - dot_product(uinf,uinf))

!  d3r1 = pemit + d1r1*uinf
!  d3r2 = pemit + d1r2*uinf

!  pout = pemit + max(d1r1,d1r2)*uinf

end subroutine ittodummy

subroutine ittozero(pcpt,phyper,uinf,pzero)
  use constants
  implicit none

  real(kind=8),dimension(3),intent(in) :: pcpt,phyper,uinf
  real(kind=8),dimension(3),intent(out) :: pzero

  pzero = phyper - uinf*(sqrt(dot_product(pcpt-phyper,pcpt-phyper))/a)
!  pzero = phyper - uinf*(abs(pcpt(1)-phyper(1))/a)

end subroutine ittozero

