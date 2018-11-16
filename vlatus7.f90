subroutine panelparam1(npan,wingpans)
  use maths
  use basictypes
  implicit none

  integer,intent(in) :: npan
  type(panel),dimension(npan),intent(inout) :: wingpans
  integer :: i

  !Calculating coordinates for vortex vertices & controlpoints, and the gradient at cpt

  do i=1,npan,1
    if (wingpans(i)%ptype == 'Q') then

      !midpoint for quadrilaterals
      wingpans(i)%cacfr = 0.25 * (wingpans(i)%pnt(1)%p%acfr+wingpans(i)%pnt(2)%p%acfr+ &
        & wingpans(i)%pnt(3)%p%acfr+wingpans(i)%pnt(4)%p%acfr)

      wingpans(i)%norm = cross(wingpans(i)%pnt(4)%p%acfr-wingpans(i)%pnt(2)%p%acfr, &
        & wingpans(i)%pnt(3)%p%acfr-wingpans(i)%pnt(1)%p%acfr)/ &
        & sqrt(dot_product(cross(wingpans(i)%pnt(4)%p%acfr-wingpans(i)%pnt(2)%p%acfr, &
        & wingpans(i)%pnt(3)%p%acfr-wingpans(i)%pnt(1)%p%acfr), &
        & cross(wingpans(i)%pnt(4)%p%acfr-wingpans(i)%pnt(2)%p%acfr,wingpans(i)%pnt(3)%p%acfr-wingpans(i)%pnt(1)%p%acfr)))

      wingpans(i)%area = 0.5*sqrt(dot_product(cross(wingpans(i)%pnt(3)%p%acfr-wingpans(i)%pnt(1)%p%acfr, &
        & wingpans(i)%pnt(4)%p%acfr-wingpans(i)%pnt(2)%p%acfr),cross(wingpans(i)%pnt(3)%p%acfr-wingpans(i)%pnt(1)%p%acfr, &
        & wingpans(i)%pnt(4)%p%acfr-wingpans(i)%pnt(2)%p%acfr)))

    else if (wingpans(i)%ptype == 'T') then

      !centroid for triangles
      wingpans(i)%cacfr = 1./3. * (wingpans(i)%pnt(1)%p%acfr+wingpans(i)%pnt(2)%p%acfr+wingpans(i)%pnt(3)%p%acfr)

      wingpans(i)%norm = cross(wingpans(i)%pnt(3)%p%acfr-wingpans(i)%pnt(1)%p%acfr, &
        & wingpans(i)%pnt(2)%p%acfr-wingpans(i)%pnt(1)%p%acfr)/ &
        & sqrt(dot_product(cross(wingpans(i)%pnt(3)%p%acfr-wingpans(i)%pnt(1)%p%acfr, &
        & wingpans(i)%pnt(2)%p%acfr-wingpans(i)%pnt(1)%p%acfr), &
        & cross(wingpans(i)%pnt(3)%p%acfr-wingpans(i)%pnt(1)%p%acfr,wingpans(i)%pnt(2)%p%acfr-wingpans(i)%pnt(1)%p%acfr)))

      wingpans(i)%area = 0.5*sqrt(dot_product(cross(wingpans(i)%pnt(2)%p%acfr-wingpans(i)%pnt(1)%p%acfr, &
        & wingpans(i)%pnt(3)%p%acfr-wingpans(i)%pnt(1)%p%acfr),cross( wingpans(i)%pnt(2)%p%acfr-wingpans(i)%pnt(1)%p%acfr, &
        & wingpans(i)%pnt(3)%p%acfr-wingpans(i)%pnt(1)%p%acfr)))

    end if
  end do

end subroutine panelparam1





subroutine interppnt(pcpt,tsteps,i,tstop,a,deltat,vpoint,t1,t2,d1r1,d1r2,ff)
  use pg, only:k2,intp,intq
  use fourdgridtypes
  implicit none

  integer,intent(in) :: tsteps,i,tstop
  real*8,intent(in) :: a,deltat
  real*8,dimension(3),intent(in) :: pcpt
  type(point),intent(in) :: vpoint
  integer,intent(out) :: t1,t2
  logical,intent(out) :: ff

  real*8,intent(out) :: d1r1,d1r2
  integer :: l
  real*8 :: d1r3

  real*8,dimension(3) :: d3r1
  real*8 :: ks1,ks2,f1,f2,ls,as,bs,cs,check,dt1,dt2

  l = 1

  d1r1 = k2(a,deltat,0,vpoint%set(i)%infr(1:3)-pcpt)

  ff = .true.

  if (d1r1 >= 0.) then

    do l=i-1,tstop,-1
      d1r2 = k2(a,deltat,i-l,vpoint%set(l)%infr(1:3)-pcpt)

      if (d1r2 >= 0.) then
        d1r1 = d1r2
      else
        t1 = l+1
        t2 = l
        dt1 = dble(t1)
        dt2 = dble(t2)
        as = dot_product(vpoint%set(t1)%infr(1:3)-vpoint%set(t2)%infr(1:3),vpoint%set(t1)%infr(1:3)-vpoint%set(t2)%infr(1:3))/ &
           & (a*deltat)**2 - (dt1-dt2)**2
        bs = 2.*dot_product(vpoint%set(t1)%infr(1:3)-vpoint%set(t2)%infr(1:3),vpoint%set(t2)%infr(1:3)-pcpt)/(a*deltat)**2 + &
           & 2.*(dt1-dt2)*(dble(i)-dt2)
        cs = dot_product(vpoint%set(t2)%infr(1:3)-pcpt,vpoint%set(t2)%infr(1:3)-pcpt)/(a*deltat)**2 - (dble(i)-dt2)**2

        ls = min((-1.*bs+sqrt(bs**2-4.*as*cs))/(2.*as),(-1.*bs-sqrt(bs**2-4.*as*cs))/(2.*as))

        d1r2 = ls
        d1r1 = 1.-ls

        ff = .false.
        exit
      end if

      if (l == tstop) then
        t1 = i !tstop
        t2 = i !tstop
        d1r1 = 0.
        d1r2 = 1.
        ff = .true.
      end if
    end do

    if (i-1 < tstop) then
      t1 = i
      t2 = i !max(max(1,i-1),tstop)
      d1r1 = 0.
      d1r2 = 1.
      ff = .true.
    end if

  else
    t1 = i
    t2 = i !max(max(1,i-1),tstop)
    d1r1 = 0.
    d1r2 = 1.
    ff = .true.
  end if

!if (ff .eqv. .false.) then
!  d3r1 = intp(vpoint%set(t1)%infr(1:3),vpoint%set(t2)%infr(1:3),d1r1,d1r2)
!  print*, d3r1
!  print*, vpoint%set(t1)%infr(1:3)
!  print*, vpoint%set(t2)%infr(1:3)
!  print*, sqrt(dot_product(d3r1-pcpt,d3r1-pcpt))/(a*deltat) - (dble(i)-intq(dble(t1),dble(t2),d1r1,d1r2))
!  print*
!end if

!!**************
!if (i<10) then
!    t1 = i
!    t2 = i
!    d1r1 = 0.
!    d1r2 = 1.
!    ff = .true.
!end if
!!***************


end subroutine interppnt



module infl
  implicit none

  contains

function velindlinseg(pcpt,p1,p2,c1,c2)
  use maths, only:cross
  use input
  implicit none
  real*8,dimension(3) :: velindlinseg
  real*8,dimension(3),intent(in) :: pcpt,p1,p2
  real*8,intent(in) :: c1,c2
  real*8,dimension(3) :: r0,r1,r2,x

  r0 = p2 - p1
  r1 = p1 - pcpt
  r2 = p2 - pcpt

  !calculate the influence of the vortex element confined between points ra & rb (arbitrary orientation) Biot-Savart
  if (dot_product(r1,r1) > eoc .and. dot_product(r2,r2) > eoc .and. dot_product(r0,r0) > eoc) then
    x = p1 - r0 * dot_product(r1,r0)/dot_product(r0,r0)
    if (sqrt(dot_product(x-pcpt,x-pcpt)) > 1.e-6) then
      velindlinseg = (1./(4.*pi)) * (cross(r1,r2)/dot_product(r0,r0)) * &
        & (((c2-c1)*(1./sqrt(dot_product(r2,r2)) - 1./sqrt(dot_product(r1,r1)))) - &
        & (((c2*dot_product(r0,r1) - c1*dot_product(r0,r2))/dot_product(cross(r1,r2),cross(r1,r2))) * &
        & (dot_product(r0,r2)/sqrt(dot_product(r2,r2)) - dot_product(r0,r1)/sqrt(dot_product(r1,r1)))))
    else
      velindlinseg = 0.
    end if
  else
!    print*, 'Error in velindlinseg', r0,r1,r2
    velindlinseg = 0.
  end if

  !prevent any degeneracies
!  x = p1 + r0 * dot_product(r1,r0)/dot_product(r0,r0)
!  if (sqrt(dot_product(x-pcpt,x-pcpt)) <= 1.e-6) then
!    velindlinseg = 0.
!  end if

end function velindlinseg

end module infl




subroutine testvelind()
  use infl
  implicit none

  real*8,dimension(3) :: r1,r2,c,vind
  real*8 :: circ,circ1,circ2

  circ = 0.6
  circ1 = 1.
  circ2 = 0.5

  r1 = (/1.,1.,0./)
  r2 = (/1.,-1.,0./)
  c = (/0.,0.,0./)

  print*, velindlinseg(c,r1,r2,circ1,circ2)

  call inflseg(r1,r2,c,vind)
  vind = circ*vind
  print*, vind
  
end subroutine testvelind
  

!***********************************************************************************************************************************
!***********************************************************************************************************************************

subroutine bicgstabup(n,prec,A,b,x0,xout)
  implicit none

  integer,intent(in) :: n
  real*8,intent(in) :: prec
  real*8,dimension(n,n),intent(in) :: A
  real*8,dimension(n),intent(in) :: b,x0
  real*8,dimension(n),intent(out) :: xout

  real*8 :: alpha,beta,rho,om,rhom1,omm1
  real*8,dimension(n) :: r,rh,v,p,rm1,vm1,pm1,s,t,x,xm1,check
  integer :: i,j

  xm1 = x0
  rm1 = b - matmul(A,x0)
  rh = rm1
  rhom1 = 1.
  alpha = 1.
  omm1 = 1.

  vm1 = 0.
  pm1 = 0.

  do i=1,100,1
    rho = dot_product(rh,rm1)
    beta = rho/rhom1 * alpha/omm1
    p = rm1 + beta*(pm1-omm1*vm1)
    v = matmul(A,p)
    alpha = rho/dot_product(rh,v)
    s = rm1 - alpha*v
    t = matmul(A,s)
    om = dot_product(t,s)/dot_product(t,t)
    x = xm1 + alpha*p + om*s

    do j=1,n,1
      check(j) = abs((alpha*p(j) + om*s(j))/x(j))
    end do
    if (maxval(check) < prec) then
      print*, 'BICGSTABup iterations:', i
      exit
    end if

    r = s - om*t

    rhom1 = rho
    omm1 = om
    rm1 = r
    vm1 = v
    pm1 = p
    xm1 = x
    if (i==100) then
      print*, 'BICGSTABup: precision not achieved', maxval(check)
    end if
  end do

  xout = x

end subroutine bicgstabup

!***********************************************************************************************************************************
!***********************************************************************************************************************************


!!***********************************************************************************************************************************
!!***********************************************************************************************************************************


!***********************************************************************************************************************************
!***********************************************************************************************************************************

subroutine dataout3(panfstream,wingpans,force,pancp)
  use input
  use maths, only: cross
  use rotation, only: ACINrot
  use basictypes
  implicit none

  type(panel),dimension(npan),intent(in) :: wingpans
  real*8,dimension(3,npan,tsteps),intent(in) :: panfstream
  real*8,dimension(3,npan,tsteps),intent(in) :: force
  real*8,dimension(npan,tsteps),intent(in) :: pancp
  real*8,dimension(npan) :: cps

  integer :: i,j,k,ioout,d1i1
  real*8,dimension(3,npan,tsteps) :: forceIN
  real*8,dimension(3,npan) :: forceaux
  real*8,dimension(tsteps) :: lift,drag
  real*8 :: totalA,liftaux,lpg,lkt,lat,M
  real*8, parameter :: g=1.4


  totalA = 0.
  do i=1,npan,1
    totalA = totalA + wingpans(i)%area
  end do

do k=1,tsteps,1
do i=1,npan,1
do j=1,3,1
if (isnan(force(j,i,k))) then
print*, 'nan in force',i,j
end if
end do
end do
end do


  forceIN = force
  drag = 0.
  lift = 0.
  do i=1,tsteps,1
    do j=1,npan,1
      lift(i) = lift(i) + force(3,j,i) !/(0.5*density*dot_product(panfstream(:,j),panfstream(:,j)))
      drag(i) = drag(i) + force(1,j,i) !/(0.5*density*dot_product(panfstream(:,j),panfstream(:,j)))
    end do
if (isnan(drag(i)))print*, 'nanan'
  end do

    print*, 'Cl'
  do i=1,tsteps,1
    print*, lift(i)/(0.5*density*uinf**2*totalA)
!    print*, lift(i)/totalA
  end do

!    print*, 'Cd'
!  do i=1,tsteps,1
!    print*, drag(i)/(0.5*density*uinf**2*totalA)
!  end do

  !---------------------------------------------------------------------------------------------------------------------------------



            





  open (unit=34,file='output/cl.dat',status='replace',action='write',iostat=ioout)

  do i=1,tsteps,1
    write (34,*) lift(i)/(0.5*density*uinf**2*totalA)
  end do

  close(34)


  open (unit=34,file='output/cd.dat',status='replace',action='write',iostat=ioout)

  do i=1,tsteps,1
    write (34,*) drag(i)/(0.5*density*uinf**2*totalA)
  end do

  close(34)




!  open (unit=94,file='cl_compressible.dat',status='replace',action='write',iostat=ioout)

!  write (94,*) 'PG, Karman-Tsien, Laitone'

!  do k=1,19,1

!  M = 1./20. *k 



!!  forceaux = force(:,:,tsteps)

  !PG-scaling
!  do i=1,npan,1
!!    cps(i) = pancp(i,tsteps)/sqrt(1.-M**2)
!    forceaux(:,i) = forceaux(:,i)/sqrt(1.-M**2)
!  end do

!  liftaux = 0.
!  do j=1,npan,1
!    liftaux = liftaux + force(3,j,tsteps)/sqrt(1.-M**2) !forceaux(3,j)
!  end do

!  lpg = liftaux/(0.5*density*uinf**2 *totalA)

!!  forceaux = force(:,:,tsteps)

!!  !Karman-Tsien
!!  do i=1,npan,1
!!    forceaux(:,i) = forceaux(:,i)/(sqrt(1.-M**2)+ (M**2/(1+sqrt(1.-M**2)))*pancp(i,tsteps)/2)
!!!    cps(i) = pancp(i,tsteps)/(sqrt(1.-M**2)+ (M**2/(1+sqrt(1.-M**2)))*pancp(i,tsteps)/2)
!!  end do

!  liftaux = 0.
!  do j=1,npan,1
!    liftaux = liftaux + force(3,j,tsteps)/(sqrt(1.-M**2)+ (M**2/(1+sqrt(1.-M**2)))*pancp(j,tsteps)/2) !forceaux(3,j)
!  end do
!  lkt = liftaux/(0.5*density*uinf**2 *totalA)


!!  forceaux = force(:,:,tsteps)

!!  !Laitone
!!  do i=1,npan,1
!!    forceaux(:,i) = forceaux(:,i)/(sqrt(1.-M**2)+M**2*(1.+ (g-1)*M**2/2.)*sqrt(1.-M**2)*pancp(i,tsteps)/2)
!!!    cps(i) = (-1.)*pancp(i,tsteps)/(sqrt(1.-M**2)+M**2*(1.+ (g-1)*M**2/2.)*sqrt(1.-M**2)*pancp(i,tsteps)/2)
!!  end do

!  liftaux = 0.
!  do j=1,npan,1
!    liftaux = liftaux + force(3,j,tsteps)/(sqrt(1.-M**2)+M**2*(1.+ (g-1)*M**2/2.)*sqrt(1.-M**2)*pancp(j,tsteps)/2) ! forceaux(3,j)
!  end do

!  lat = liftaux/(0.5*density*uinf**2 *totalA)

!    write (94,*) M, ', ', lpg, ', ', lkt, ', ', lat

!  end do

!  close(94)

end subroutine dataout3

!***********************************************************************************************************************************
!***********************************************************************************************************************************
!***********************************************************************************************************************************
!***********************************************************************************************************************************

subroutine forcecalcM(t,tsteps,deltat,density,np,npan,points,d,wingpans,pancirc,panvindin,panfstream,force,cps)
  use pg
  use maths, only: cross,cycperm
  use basictypes
  use rotation
  use constants, only: pi
  use functions
  implicit none

  real*8,intent(in) :: deltat,density,d
  integer,intent(in) :: tsteps,np,npan,t
  type(panel),dimension(npan),intent(in) :: wingpans
  type(point),dimension(np),intent(in) :: points
!  real*8,dimension(3,3),intent(in) :: acin
  real*8,dimension(3,npan),intent(in) :: panvindin
  real*8,dimension(3,npan,tsteps),intent(in) :: panfstream
  real*8,dimension(npan,tsteps),intent(in) :: pancirc
  real*8,dimension(3,npan) :: panvind
  real*8,dimension(3,npan),intent(out) :: force
  real*8,dimension(npan),intent(out) :: cps
  real*8,dimension(npan) :: vindincomp
  real*8 :: dttilde,M,beta
  real*8,dimension(3) :: dummy3r,d3r1,d3r2
  integer :: i,j,dummy1i1
  logical,dimension(npan,4) :: processed

  real*8, parameter :: a = 340.

  !calculate force, then cp for the steady incompressible problem
  processed = .false.
  force = 0.

  do i=1,npan,1
    panvind(:,i) = panvindin(:,i)
  end do

  do i=1,npan,1
    if (wingpans(i)%ptype == 'Q') then
      do j=1,4,1
        if (processed(i,j) .eqv. .false.) then
          d3r1 = wingpans(i)%pnt(j)%p%set(t)%infr(1:3)
          d3r2 = wingpans(i)%pnt(cycperm(4,1,j))%p%set(t)%infr(1:3)
          if (wingpans(i)%neigh(j) > 0) then
            force(:,i) = force(:,i) + 0.5*(pancirc(i,t)-pancirc(wingpans(i)%neigh(j),t))*density*cross(panvind(:,i), &
                       & (d3r2-d3r1))
          else if (wingpans(i)%neigh(j) == 0) then
            !do nothing, see Kutta condition
          else if (wingpans(i)%neigh(j) < 0) then
            force(:,i) = force(:,i) + density*pancirc(i,t)*cross(panvind(:,i), &
                       & (d3r2-d3r1))
          end if
          processed(i,j) = .true.
        end if
      end do
    else if (wingpans(i)%ptype == 'T') then
      do j=1,3,1
        if (processed(i,j) .eqv. .false.) then
          if (wingpans(i)%neigh(j) > 0) then
            force(:,i) = force(:,i) + 0.5*(pancirc(i,t)-pancirc(wingpans(i)%neigh(j),t))*density*cross(panvind(:,i), &
                       & (wingpans(i)%pnt(cycperm(3,1,j))%p%set(t)%infr(1:3)- &
                       & wingpans(i)%pnt(j)%p%set(t)%infr(1:3)))
          else if (wingpans(i)%neigh(j) == 0) then
            !do nothing, see Kutta condition
          else if (wingpans(i)%neigh(j) < 0) then
            force(:,i) = force(:,i) + density*pancirc(i,t)*cross(panvind(:,i), &
                       & (wingpans(i)%pnt(cycperm(3,1,j))%p%set(t)%infr(1:3)- &
                       & wingpans(i)%pnt(j)%p%set(t)%infr(1:3)))
          end if
          processed(i,j) = .true.
        end if
      end do
    end if

    !no unsteady correction?!
!    call ACtoINrot(wingpans(i)%norm,acin,dummy3r)
    d3r2 = cross(wingpans(i)%pnt(3)%p%set(t)%infr(1:3)-wingpans(i)%pnt(1)%p%set(t)%infr(1:3), &
      & wingpans(i)%pnt(2)%p%set(t)%infr(1:3)-wingpans(i)%pnt(4)%p%set(t)%infr(1:3))
    d3r2 = d3r2/sqrt(dot_product(d3r2,d3r2))
    if (t > 1) then
      force(:,i) = force(:,i) + density*wingpans(i)%area*d3r2*(pancirc(i,t)-pancirc(i,t-1))/deltat
    else
      force(:,i) = force(:,i) + density*wingpans(i)%area*d3r2*pancirc(i,t)/deltat
    end if

  end do


do i=1,npan,1
if (isnan(pancirc(i,t))) print*, 'NaN in pancirc'
if (t>1) then
if (isnan(pancirc(i,t-1))) print*, 'NaN in pancirc t-1'
end if
if (isnan(wingpans(i)%area)) print*, 'NaN wingpan area'
do j=1,3,1
if (isnan(force(j,i))) print*, 'NaN in calc force', j,i,t
end do
end do


  do i=1,npan,1
    cps(i) = sqrt(dot_product(force(:,i),force(:,i)))/(wingpans(i)%area* &
      & dot_product(panfstream(:,i,t),panfstream(:,i,t))*0.5*density)
!    call ACtoINrot(wingpans(i)%norm,acin,dummy3r)
    d3r2 = cross(wingpans(i)%pnt(3)%p%set(t)%infr(1:3)-wingpans(i)%pnt(1)%p%set(t)%infr(1:3), &
      & wingpans(i)%pnt(2)%p%set(t)%infr(1:3)-wingpans(i)%pnt(4)%p%set(t)%infr(1:3))
    d3r2 = d3r2/sqrt(dot_product(d3r2,d3r2))
    if (dot_product(d3r2,force(:,i)) > 0.) then
      cps(i) = (-1.) * cps(i)
    end if
  end do

!  !PG Correction
!  do i=1,npan,1
!    M = sqrt(dot_product(panfstream(:,i,t),panfstream(:,i,t)))/a
!    beta = sqrt(1.-M**2)
!    cps(i) = cps(i)/beta
!    force(:,i) = force(:,i)/beta
!  end do

do i=1,npan,1
do j=1,3,1
if (isnan(force(j,i))) then
print*, 'nan in force',i,j
end if
end do
end do


end subroutine forcecalcM

!***********************************************************************************************************************************

!subroutine forcecalc2(t,step,tsteps,deltat,acin,density,np,npan,points,wingpans,pancirc,panvind,panfstream,force,cp)
!  use maths, only: cross,cycperm
!  use wingpanels
!  use rotation
!  implicit none

!  real*8,intent(in) :: deltat,density
!  integer,intent(in) :: np,npan,step,tsteps,t
!  type(spoint),dimension(np),intent(in) :: points
!  type(wingpanel),dimension(npan),intent(in) :: wingpans
!  real*8,dimension(3,3),intent(in) :: acin
!  real*8,dimension(3,npan),intent(in) :: panvind,panfstream
!  real*8,dimension(npan,tsteps),intent(in) :: pancirc

!  real*8,dimension(3,npan),intent(out) :: force
!  real*8,dimension(npan),intent(out) :: cp

!  real*8,dimension(3) :: dummy3r
!  logical,dimension(npan,4) :: processed
!  integer :: i,j

!  processed = .false.
!  force = 0.

!  do i=1,npan,1
!    if (wingpans(i)%ptype == 'Q') then
!      do j=1,4,1
!        if (processed(i,j) .eqv. .false.) then
!          if (wingpans(i)%neigh(j) > 0) then
!            force(:,i) = force(:,i) + 0.5*(pancirc(i,step)-pancirc(wingpans(i)%neigh(j),step))*density*cross(panvind(:,i), &
!                       & (points(wingpans(i)%pan(cycperm(4,1,j)))%inrt(:,t)-points(wingpans(i)%pan(j))%inrt(:,t)))
!          else if (wingpans(i)%neigh(j) == 0) then
!            !do nothing, see Kutta condition
!          else if (wingpans(i)%neigh(j) < 0) then
!            force(:,i) = force(:,i) + density*pancirc(i,step)*cross(panvind(:,i), &
!                       & (points(wingpans(i)%pan(cycperm(4,1,j)))%inrt(:,t)-points(wingpans(i)%pan(j))%inrt(:,t)))
!          end if
!          processed(i,j) = .true.
!        end if
!      end do
!    else if (wingpans(i)%ptype == 'T') then
!      do j=1,3,1
!        if (processed(i,j) .eqv. .false.) then
!          if (wingpans(i)%neigh(j) > 0) then
!            force(:,i) = force(:,i) + 0.5*(pancirc(i,step)-pancirc(wingpans(i)%neigh(j),step))*density*cross(panvind(:,i), &
!                       & (points(wingpans(i)%pan(cycperm(3,1,j)))%inrt(:,t)-points(wingpans(i)%pan(j))%inrt(:,t)))
!          else if (wingpans(i)%neigh(j) == 0) then
!            !do nothing, see Kutta condition
!          else if (wingpans(i)%neigh(j) < 0) then
!            force(:,i) = force(:,i) + density*pancirc(i,step)*cross(panvind(:,i), &
!                       & (points(wingpans(i)%pan(cycperm(3,1,j)))%inrt(:,t)-points(wingpans(i)%pan(j))%inrt(:,t)))
!          end if
!          processed(i,j) = .true.
!        end if
!      end do
!    end if
!    !unsteady correction
!    call ACtoINrot(wingpans(i)%norm,acin,dummy3r)
!    if (step > 1) then
!      force(:,i) = force(:,i) + density*wingpans(i)%area*dummy3r*(pancirc(i,step)-pancirc(i,step-1))/deltat
!    else
!      force(:,i) = force(:,i) + density*wingpans(i)%area*dummy3r*pancirc(i,step)/deltat
!    end if
!  end do

!  do i=1,npan,1
!    cp(i) = (-1.)*sqrt(dot_product(force(:,i),force(:,i)))/(wingpans(i)%area* &
!      & dot_product(panfstream(:,i),panfstream(:,i))*0.5*density)
!    call ACtoINrot(wingpans(i)%norm,acin,dummy3r)
!    if (dot_product(dummy3r,force(:,i)) < 0.) then
!      cp(i) = (-1.) * cp(i)
!    end if
!  end do

!end subroutine forcecalc2

!***********************************************************************************************************************************
!***********************************************************************************************************************************
!***********************************************************************************************************************************
!***********************************************************************************************************************************
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

!subroutine gridparam1(np,npan,points,wingpans)
!  use maths
!  use wingpanels
!  implicit none

!  type(wingpanel),dimension(npan),intent(inout) :: wingpans
!  integer,intent(in) :: np,npan
!  type(spoint),dimension(np),intent(in) :: points
!  integer :: i

!  !Calculating coordinates for vortex vertices & controlpoints, and the gradient at cpt

!  do i=1,npan,1

!    if (wingpans(i)%ptype == 'Q') then

!      !midpoint for quadrilaterals
!      wingpans(i)%cpt%pAC = 0.25 * (points(wingpans(i)%pan(1))%ac+points(wingpans(i)%pan(2))%ac+ &
!        & points(wingpans(i)%pan(3))%ac+points(wingpans(i)%pan(4))%ac)

!      wingpans(i)%norm = cross(points(wingpans(i)%pan(4))%ac-points(wingpans(i)%pan(2))%ac, &
!        & points(wingpans(i)%pan(3))%ac-points(wingpans(i)%pan(1))%ac)/ &
!        & sqrt(dot_product(cross(points(wingpans(i)%pan(4))%ac-points(wingpans(i)%pan(2))%ac, &
!        & points(wingpans(i)%pan(3))%ac-points(wingpans(i)%pan(1))%ac), &
!        & cross(points(wingpans(i)%pan(4))%ac-points(wingpans(i)%pan(2))%ac, &
!        & points(wingpans(i)%pan(3))%ac-points(wingpans(i)%pan(1))%ac)))

!      wingpans(i)%area = 0.5*sqrt(dot_product(cross(points(wingpans(i)%pan(3))%ac-points(wingpans(i)%pan(1))%ac, &
!        & points(wingpans(i)%pan(4))%ac-points(wingpans(i)%pan(2))%ac), &
!        & cross(points(wingpans(i)%pan(3))%ac-points(wingpans(i)%pan(1))%ac, &
!        & points(wingpans(i)%pan(4))%ac-points(wingpans(i)%pan(2))%ac)))

!    else if (wingpans(i)%ptype == 'T') then

!      !centroid for triangles
!      wingpans(i)%cpt%pAC = 1./3. * (points(wingpans(i)%pan(1))%ac+points(wingpans(i)%pan(2))%ac+points(wingpans(i)%pan(3))%ac)

!      wingpans(i)%norm = cross(points(wingpans(i)%pan(3))%ac-points(wingpans(i)%pan(1))%ac, &
!        & points(wingpans(i)%pan(2))%ac-points(wingpans(i)%pan(1))%ac)/ &
!        & sqrt(dot_product(cross(points(wingpans(i)%pan(3))%ac-points(wingpans(i)%pan(1))%ac, &
!        & points(wingpans(i)%pan(2))%ac-points(wingpans(i)%pan(1))%ac), &
!        & cross(points(wingpans(i)%pan(3))%ac-points(wingpans(i)%pan(1))%ac, &
!        & points(wingpans(i)%pan(2))%ac-points(wingpans(i)%pan(1))%ac)))

!      wingpans(i)%area = 0.5*sqrt(dot_product(cross(points(wingpans(i)%pan(2))%ac-points(wingpans(i)%pan(1))%ac, &
!        & points(wingpans(i)%pan(3))%ac-points(wingpans(i)%pan(1))%ac), &
!        & cross(points(wingpans(i)%pan(2))%ac-points(wingpans(i)%pan(1))%ac, &
!        & points(wingpans(i)%pan(3))%ac-points(wingpans(i)%pan(1))%ac)))

!    end if
!  end do

!end subroutine gridparam1

!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------

subroutine ringgrid(np,npan,points,panels,neigh,ptype,pvring,cpt,vrnorm,parea)
  use maths
  implicit none

  integer,intent(in) :: np,npan
  real*8,dimension(np,3),intent(in) :: points
  integer,dimension(npan,4),intent(in) :: panels,neigh
  character(len=1),dimension(npan),intent(in) :: ptype
  real*8,dimension(np,3),intent(out) :: pvring
  real*8,dimension(npan,3),intent(out) :: cpt,vrnorm
  real*8,dimension(npan),intent(out) :: parea
  integer :: i,j

  pvring = points
 
  !Calculating coordinates for vortex vertices & controlpoints, and the gradient at cpt
  do i=1,npan,1
    if (ptype(i) == 'Q') then
      do j=1,4,1
        if (neigh(i,j) == -1) then
          if ((neigh(i,cycperm(4,3,j)) > i .or. neigh(i,cycperm(4,3,j)) <= 0) .and. &
            & (neigh(i,cycperm(4,1,j)) > i .or. neigh(i,cycperm(4,3,j)) <= 0)) then
            pvring(panels(i,j),:) = 0.75 * pvring(panels(i,j),:) + 0.25 * pvring(panels(i,cycperm(4,3,j)),:)
            pvring(panels(i,cycperm(4,1,j)),:) = 0.75*pvring(panels(i,cycperm(4,1,j)),:) + 0.25*pvring(panels(i,cycperm(4,2,j)),:)
          end if
        end if
      end do
    end if
  end do

  do i=1,npan,1
    if (ptype(i) == 'Q') then
      !midpoint for quadrilaterals
      cpt(i,:) = 1./4. * (pvring(panels(i,1),:)+pvring(panels(i,2),:)+pvring(panels(i,3),:)+pvring(panels(i,4),:))
      vrnorm(i,:) = cross(pvring(panels(i,4),:)-pvring(panels(i,2),:),pvring(panels(i,3),:)-pvring(panels(i,1),:))/ &
             & sqrt(dot_product(cross(pvring(panels(i,4),:)-pvring(panels(i,2),:),pvring(panels(i,3),:)-pvring(panels(i,1),:)), &
             & cross(pvring(panels(i,4),:)-pvring(panels(i,2),:),pvring(panels(i,3),:)-pvring(panels(i,1),:))))
      parea(i) = 0.5*sqrt(dot_product(cross(pvring(panels(i,3),:)-pvring(panels(i,1),:),pvring(panels(i,4),:)- &
                 & pvring(panels(i,2),:)),cross(pvring(panels(i,3),:)-pvring(panels(i,1),:),pvring(panels(i,4),:)- &
                 & pvring(panels(i,2),:))))
    else if (ptype(i) == 'T') then
      !centroid for triangles
      cpt(i,:) = 1./3. * (pvring(panels(i,1),:)+pvring(panels(i,2),:)+pvring(panels(i,3),:))
      vrnorm(i,:) = cross(pvring(panels(i,3),:)-pvring(panels(i,1),:),pvring(panels(i,2),:)-pvring(panels(i,1),:))/ &
            & sqrt(dot_product(cross(pvring(panels(i,3),:)-pvring(panels(i,1),:),pvring(panels(i,2),:)-pvring(panels(i,1),:)), &
            & cross(pvring(panels(i,3),:)-pvring(panels(i,1),:),pvring(panels(i,2),:)-pvring(panels(i,1),:))))
      parea(i) = 0.5*sqrt(dot_product(cross(pvring(panels(i,2),:)-pvring(panels(i,1),:),pvring(panels(i,3),:)- &
                 & pvring(panels(i,1),:)),cross(pvring(panels(i,2),:)-pvring(panels(i,1),:),pvring(panels(i,3),:)- &
                 & pvring(panels(i,1),:))))
    end if
  end do

end subroutine ringgrid

!-----------------------------------------------------------------------------------------------------------------------------------

!***********************************************************************************************************************************
!This subroutine calculates the influence coefficient of a vortexringpanel for some arbitrary collocation point.
!Just multiply by its circulation to get the velocity
!***********************************************************************************************************************************

subroutine vortexring(cpt,verta,vertb,vertc,vertd,vtype,circ,c)
  implicit none

  real*8,dimension(3),intent(in) :: cpt,verta,vertb,vertc,vertd
  character(len=1),intent(in) :: vtype
  real*8, intent(in) :: circ
  real*8, dimension(3),intent(out) :: c

  real*8, dimension(3) :: cseg

  interface
    subroutine inflseg(ra,rb,rcpt,cseg)
    implicit none
    real*8,dimension(3),intent(in) :: ra,rb,rcpt
    real*8,dimension(3),intent(out) :: cseg
    real*8,parameter :: pi = 3.14159265359
    real*8, dimension(3) :: r0, r1, r2
    real*8,dimension(3) :: x
    end subroutine inflseg
  end interface

  cseg = 0.
  c = 0.

  if (vtype == 'Q') then
    call inflseg(verta,vertb,cpt,cseg)
    c = cseg
    call inflseg(vertb,vertc,cpt,cseg)
    c = c + cseg
    call inflseg(vertc,vertd,cpt,cseg)
    c = c + cseg
    call inflseg(vertd,verta,cpt,cseg)
    c = c + cseg
  else if (vtype == 'T') then
    call inflseg(verta,vertb,cpt,cseg)
    c = cseg
    call inflseg(vertb,vertc,cpt,cseg)
    c = c + cseg
    call inflseg(vertc,verta,cpt,cseg)
    c = c + cseg
  end if

  c = c * circ

end subroutine vortexring

!-----------------------------------------------------------------------------------------------------------------------------------

subroutine vortexshoe(cpt,verta,vertb,vertc,vertd,vtype,vneigh,circ,c)
  implicit none

  real*8,dimension(3),intent(in) :: cpt,verta,vertb,vertc,vertd
  character(len=1),intent(in) :: vtype
  integer,dimension(4),intent(in) :: vneigh
  real*8, intent(in) :: circ
  real*8, dimension(3),intent(out) :: c

  real*8, dimension(3) :: cseg

  interface
    subroutine inflseg(ra,rb,rcpt,cseg)
    implicit none
    real*8,dimension(3),intent(in) :: ra,rb,rcpt
    real*8,dimension(3),intent(out) :: cseg
    real*8,parameter :: pi = 3.14159265359
    real*8, dimension(3) :: r0, r1, r2
    real*8,dimension(3) :: x
    end subroutine inflseg
  end interface

  cseg = 0.
  c = 0.

  if (vtype == 'Q') then
    if (vneigh(1) /= -2) then
      call inflseg(verta,vertb,cpt,cseg)
      c = cseg
    end if
    if (vneigh(2) /= -2) then
      call inflseg(vertb,vertc,cpt,cseg)
      c = c + cseg
    end if
    if (vneigh(3) /= -2) then
      call inflseg(vertc,vertd,cpt,cseg)
      c = c + cseg
    end if
    if (vneigh(4) /= -2) then
      call inflseg(vertd,verta,cpt,cseg)
      c = c + cseg
    end if
  else if (vtype == 'T') then
    if (vneigh(1) /= -2) then
      call inflseg(verta,vertb,cpt,cseg)
      c = cseg
    end if
    if (vneigh(2) /= -2) then
      call inflseg(vertb,vertc,cpt,cseg)
      c = c + cseg
    end if
    if (vneigh(3) /= -2) then
     call inflseg(vertc,verta,cpt,cseg)
      c = c + cseg
    end if
  end if

  c = c * circ

end subroutine vortexshoe

subroutine vortexsegment(cpt,verta,vertb,circ,vind)
  implicit none

  real*8,dimension(3),intent(in) :: cpt,verta,vertb
  real*8, intent(in) :: circ
  real*8, dimension(3),intent(out) :: vind

  real*8, dimension(3) :: cseg

  interface
    subroutine inflseg(ra,rb,rcpt,cseg)
    implicit none
    real*8,dimension(3),intent(in) :: ra,rb,rcpt
    real*8,dimension(3),intent(out) :: cseg
    real*8,parameter :: pi = 3.14159265359
    real*8, dimension(3) :: r0, r1, r2
    real*8,dimension(3) :: x
    end subroutine inflseg
  end interface

  call inflseg(verta,vertb,cpt,cseg)

  vind = cseg*circ

end subroutine vortexsegment


!-----------------------------------------------------------------------------------------------------------------------------------

subroutine inflinf(ra,rcpt,cseg)
  implicit none

  real*8,dimension(3),intent(in) :: ra,rcpt
  real*8,dimension(3),intent(out) :: cseg

  real*8,parameter :: pi = 3.14159265359
  real*8, dimension(3) :: r1
  real*8,dimension(3) :: x,nx

  r1 = rcpt-ra

  cseg(1) = 0d0
  cseg(2) = (1./(4.*pi)) * (r1(3)/(r1(3)**2 + r1(2)**2)) * (1. + (r1(1) / sqrt(dot_product(r1,r1))))
  cseg(3) = (-1.)/(4.*pi) * (r1(2)/(r1(3)**2 + r1(2)**2)) * (1. + (r1(1) / sqrt(dot_product(r1,r1))))

  nx(1) = 1.
  nx(2) = 0.
  nx(3) = 0.
  
  x = ra+(dot_product(rcpt,nx)-dot_product(ra,nx))*nx

  if (sqrt(dot_product(x-rcpt,x-rcpt)) <= 1.e-6) then
    cseg = 0.
  end if

end subroutine inflinf
!-----------------------------------------------------------------------------------------------------------------------------------

subroutine inflseg(ra,rb,rcpt,cseg)
  use maths
  implicit none

  real*8,dimension(3),intent(in) :: ra,rb,rcpt
  real*8,dimension(3),intent(out) :: cseg

  real*8, dimension(3) :: r0, r1, r2
  real*8,dimension(3) :: x
  
  r0 = rb - ra
!  r1 = rcpt - ra
!  r2 = rcpt -rb
  r1 = ra - rcpt
  r2 = rb - rcpt

  !calculate the influence of the vortex element confined between points ra & rb (arbitrary orientation) Biot-Savart
  if (dot_product(r1,r1) /= 0. .and. dot_product(r2,r2) /= 0.) then
    cseg = (-1./(4.*pi)) * (cross(r1,r2)/dot_product(cross(r1,r2),cross(r1,r2))) * &
         & dot_product(r0,(r1/sqrt(dot_product(r1,r1))) - (r2/sqrt(dot_product(r2,r2))))
  else
    cseg = 0.
  end if

  !prevent any degeneracies
  x = ra + r0 * dot_product(r1,r0)/dot_product(r0,r0)
  if (sqrt(dot_product(x-rcpt,x-rcpt)) <= 1.e-6) then
    cseg = 0.
  end if

end subroutine inflseg
