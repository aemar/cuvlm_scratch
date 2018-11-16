

module fourdgridsf
  implicit none

  contains

  subroutine hextreeinit(deltat,ntsteps,minm,maxm,nlevel,TLB)
    use constants
    use fourdgridtypes
    implicit none
    real*8,intent(in) :: deltat
    integer,intent(in) :: ntsteps
    real*8,intent(in) :: minm,maxm
    type(node) :: TLB
    integer :: nlevel

    integer,dimension(:,:),allocatable :: treelevel
    real*8 :: f
    real*8,dimension(4,2) :: rgsize
    integer,dimension(4,2) :: igsize
    integer :: i,j,n

    f = 1.2
    do i=1,3,1
      !normalised gridsize
      rgsize(i,2) = dble(ntsteps)*f + (maxm-minm)/(a*deltat)
      rgsize(i,1) = (-1) * rgsize(i,2)
      igsize(i,2) = int(rgsize(i,2) + 1.)
      igsize(i,1) = int(rgsize(i,1) - 1.)
    end do

    print*, 'Hextree init:'
    nlevel = 0
    n = igsize(1,2)
    print*, 'Level',nlevel,':',n
    do
      nlevel = nlevel + 1
      n = int(n/2. + 0.5)
      print*, 'Level',nlevel,':',n
      if (n==1) exit
    end do

    allocate(treelevel(0:nlevel,2))

    do j=1,2,1
      treelevel(0,j) = igsize(1,j)
    end do
    
    do i=1,nlevel,1
      do j=1,2,1
        treelevel(i,j) = up(treelevel(i-1,j))
      end do
    end do

    nullify(TLB%p%parent)
    do i=1,16,1
      nullify(TLB%c(i)%child)
    end do
    TLB%centroid = 0.
    TLB%mentroid = 0.
    TLB%rclcflg = .true.
    TLB%usdflg = .false.

  end subroutine hextreeinit

!###################################################################################################################################

  subroutine insertpnt(nlevel,deltat,TLB,vtxpnt,vtxset)
    use basictypes
    use fourdgridtypes
    implicit none
    integer,intent(in) :: nlevel
    real*8,intent(in) :: deltat
    type(node) :: TLB
    type(rmset),target,intent(in) :: vtxpnt
    type(point),target,intent(in) :: vtxset
    integer :: d1i1,i
    integer,dimension(4) :: d4i1
    integer,dimension(4,0:nlevel) :: branch

    branch(:,0) = box4d(deltat,vtxpnt%infr)

    !find all boxes
    do i=1,nlevel,1
      branch(:,i) = i4parent(branch(:,i-1))
    end do

    !check
    if (dot_product(branch(:,nlevel),branch(:,nlevel)) /= 4) then
      print*, 'Error in subroutine "insertpnt" 4dgrid*.f90',dot_product(branch(:,nlevel),branch(:,nlevel)),branch(:,nlevel)
      print*, vtxpnt%infr
    end if

    TLB%rclcflg = .true.

    d1i1 = TLBchild(branch(:,nlevel))

    if (associated(TLB%c(d1i1)%child)) then
      call bstep(nlevel,nlevel,branch,d1i1,TLB,vtxpnt,vtxset)
    else
      call bcnstrct(nlevel,nlevel,branch,d1i1,TLB,vtxpnt,vtxset)
    end if

  end subroutine insertpnt

!###################################################################################################################################


  recursive subroutine bstep(nlevel,lev,branch,ipcin,pnode,vtxpnt,vtxset)
    use fourdgridtypes
    implicit none
    integer,intent(in) :: nlevel,lev,ipcin
    integer,dimension(4,0:nlevel),intent(in) :: branch
    type(node),target :: pnode
    type(rmset),target,intent(in) :: vtxpnt
    type(point),target,intent(in) :: vtxset
    type(node),pointer,save :: cnode
    integer :: d1i1,i

    cnode => pnode%c(ipcin)%child
    cnode%centroid = 0.
    cnode%mentroid = 0.
    cnode%rclcflg = .true.
    if (lev>0) then
      d1i1 = nchild(branch(:,lev),branch(:,lev-1))
      if (associated(cnode%c(d1i1)%child)) then
        call bstep(nlevel,lev-1,branch,d1i1,cnode,vtxpnt,vtxset)
      else
        call bcnstrct(nlevel,lev-1,branch,d1i1,cnode,vtxpnt,vtxset)
      end if
    else if (lev==0) then
      cnode%nvtx = cnode%nvtx+1
      if (cnode%nvtx > 20) then
        do i=1,20,1
          print*, cnode%bvtx(i)%vpnt%infr
        end do
      end if
      cnode%bvtx(cnode%nvtx)%vpnt => vtxpnt
      cnode%bvtx(cnode%nvtx)%vset => vtxset
    end if
  end subroutine bstep


  recursive subroutine bcnstrct(nlevel,lev,branch,ipcin,pnode,vtxpnt,vtxset)
    use fourdgridtypes
    implicit none
    integer,intent(in) :: nlevel,lev,ipcin
    integer,dimension(4,0:nlevel),intent(in) :: branch
    type(node),target :: pnode
    type(rmset),target,intent(in) :: vtxpnt
    type(point),target,intent(in) :: vtxset
    type(node),pointer,save :: cnode
    integer :: d1i1,i

    allocate(pnode%c(ipcin)%child)
    cnode => pnode%c(ipcin)%child
    cnode%p%parent => pnode
    cnode%centroid = 0.
    cnode%mentroid = 0.
    cnode%rclcflg = .true.
    do i=1,16,1
      nullify(cnode%c(i)%child)
    end do
    if (lev>0) then
      d1i1 = nchild(branch(:,lev),branch(:,lev-1))
      call bcnstrct(nlevel,lev-1,branch,d1i1,cnode,vtxpnt,vtxset)
    else if (lev==0) then
      cnode%nvtx = 1
      cnode%bvtx(1)%vpnt => vtxpnt
      cnode%bvtx(1)%vset => vtxset
    end if
  end subroutine bcnstrct


!###################################################################################################################################

  subroutine zerotree(nlevel,TLB)
    use fourdgridtypes
    implicit none
    integer,intent(in) :: nlevel
    type(node) :: TLB

    call zerolevel(nlevel+1,TLB)

  end subroutine zerotree



  recursive subroutine zerolevel(lev,npnt)
    use fourdgridtypes
    implicit none
    integer,intent(in) :: lev
    type(node) :: npnt
    integer :: i

    npnt%mentroid = 0.
    npnt%centroid = 0.
    npnt%p%parent => null()
    npnt%rclcflg = .false.
    if (lev > 0) then
      do i=1,16,1
        if (associated(npnt%c(i)%child)) then
          call zerolevel(lev-1,npnt%c(i)%child)
          deallocate(npnt%c(i)%child)
          npnt%c(i)%child => null()
        end if
      end do
    else if (lev == 0) then
      do i=1,npnt%nvtx,1
        npnt%bvtx(i)%vset => null()
        npnt%bvtx(i)%vpnt => null()
      end do
      npnt%nvtx = 0
    end if

  end subroutine zerolevel

!###################################################################################################################################

  subroutine prune(nlevel,deltat,TLB)
    use fourdgridtypes
    implicit none
    integer,intent(in) :: nlevel
    real*8,intent(in) :: deltat
    type(node) :: TLB

    if (TLB%usdflg) then
      call prunelevel(nlevel,deltat,TLB,nlevel+1,TLB)
    else
      print*, 'No branches to prune. Tree not used.'
    end if

  end subroutine prune



  recursive subroutine prunelevel(nlevel,deltat,TLB,lev,npnt)
    use fourdgridtypes
    implicit none
    integer,intent(in) :: nlevel,lev
    real*8,intent(in) :: deltat 
    type(node) :: TLB,npnt
    type(node),pointer :: cpnt
    type(rmset),pointer,save :: pvtx
    type(point),pointer,save :: pset
    integer :: i,j

    if (lev > 0) then
      do i=16,1,-1
        if (associated(npnt%c(i)%child)) then
          cpnt => npnt%c(i)%child
          if (cpnt%usdflg) then
            call prunelevel(nlevel,deltat,TLB,lev-1,npnt%c(i)%child)
          else
            pvtx => null()
            pset => null()
            call prunebelow(lev-1,npnt%c(i)%child,pvtx,pset)
!            call calclevel(lev-1,npnt%c(i)%child)
            pvtx%infr = npnt%c(i)%child%centroid
            pvtx%minfr = npnt%c(i)%child%mentroid
            pvtx%act = .true.
            do j=1,16,1
              if (associated(cpnt%c(j)%child)) then
                deallocate(cpnt%c(j)%child)
                cpnt%c(j)%child => null()
              end if
            end do
            call insertpnt(nlevel,deltat,TLB,pvtx,pset)
          end if
        end if
      end do
    end if

  end subroutine prunelevel



  recursive subroutine prunebelow(lev,npnt,pvtx,pset)
    use fourdgridtypes
    implicit none
    integer,intent(in) :: lev
    type(node) :: npnt
    type(rmset),pointer :: pvtx
    type(point),pointer :: pset
    
    integer :: i

    if (lev == 0) then
      if (.not. associated(pvtx)) then
        pvtx => npnt%bvtx(1)%vpnt
        pset => npnt%bvtx(1)%vset
      end if
      do i=1,npnt%nvtx,1
        npnt%bvtx(i)%vpnt%act = .false.
      end do
    else
      do i=16,1,-1
        if (associated(npnt%c(i)%child)) then
          call prunebelow(lev-1,npnt%c(i)%child,pvtx,pset)
          deallocate(npnt%c(i)%child)
          npnt%c(i)%child => null()
        end if
      end do
    end if

  end subroutine prunebelow








!###################################################################################################################################

  subroutine zerousdflg(nlevel,TLB)
    use fourdgridtypes
    implicit none
    integer,intent(in) :: nlevel
    type(node) :: TLB

    call zusdflglev(nlevel+1,TLB)

  end subroutine zerousdflg



  recursive subroutine zusdflglev(lev,npnt)
    use fourdgridtypes
    implicit none
    integer,intent(in) :: lev
    type(node) :: npnt
    integer :: i

    npnt%usdflg = .false.
    if (lev .ge. 0) then
      do i=1,16,1
        if (associated(npnt%c(i)%child)) then
          call zusdflglev(lev-1,npnt%c(i)%child)
        end if
      end do
    end if

  end subroutine zusdflglev

!###################################################################################################################################

  subroutine zerotime(t,nlevel,TLB)
    use fourdgridtypes
    implicit none
    integer,intent(in) :: t,nlevel
    type(node) :: TLB
    integer,dimension(0:nlevel) :: tbranch
    integer :: i,i1

    tbranch(0) = t
    do i=1,nlevel,1
      tbranch(i) = up(tbranch(i-1))
    end do

    i1 = int((tbranch(nlevel)+1)/2)
    call zerotimelev(tbranch,i1,nlevel,nlevel+1,TLB)

  end subroutine zerotime



  recursive subroutine zerotimelev(tbranch,i1,nlevel,lev,npnt)
    use fourdgridtypes
    implicit none
    integer,dimension(0:nlevel),intent(in) :: tbranch
    integer,intent(in) :: i1,nlevel,lev
    type(node) :: npnt
    integer :: i,j,d1i1

    npnt%mentroid = 0.
    npnt%centroid = 0.
    npnt%rclcflg = .true.
    if (lev > 1) then
      do i=1+(8*i1),8+(8*i1),1
        if (associated(npnt%c(i)%child)) then
          d1i1 = idchild(tbranch(lev-1),tbranch(lev-2))
          call zerotimelev(tbranch,d1i1,nlevel,lev-1,npnt%c(i)%child)
          jloop: do j=1,16,1
            if (associated(npnt%c(i)%child%c(j)%child)) then
              exit jloop
            else if (j==16) then
              deallocate(npnt%c(i)%child)
              npnt%c(i)%child => null()
            end if
          end do jloop
        end if
      end do
    else if (lev == 1) then
      do i=1+(8*i1),8+(8*i1),1
        if (associated(npnt%c(i)%child)) then
          deallocate(npnt%c(i)%child)
          npnt%c(i)%child => null()
        end if
      end do
    end if

  end subroutine zerotimelev

!###################################################################################################################################

  subroutine setrclcflg(t,nlevel,TLB)
    use fourdgridtypes
    implicit none
    integer,intent(in) :: t,nlevel
    type(node) :: TLB
    integer,dimension(0:nlevel) :: tbranch
    integer :: i,i1

    tbranch(0) = t
    do i=1,nlevel,1
      tbranch(i) = up(tbranch(i-1))
    end do

    i1 = int((tbranch(nlevel)+1)/2)
    call setrclcflglev(tbranch,i1,nlevel,nlevel+1,TLB)

  end subroutine setrclcflg



  recursive subroutine setrclcflglev(tbranch,i1,nlevel,lev,npnt)
    use fourdgridtypes
    implicit none
    integer,dimension(0:nlevel),intent(in) :: tbranch
    integer,intent(in) :: i1,nlevel,lev
    type(node) :: npnt
    integer :: i,j,d1i1

    npnt%mentroid = 0.
    npnt%centroid = 0.
    npnt%rclcflg = .true.
    if (lev > 0) then
      do i=1+(8*i1),8+(8*i1),1
        if (associated(npnt%c(i)%child)) then
          d1i1 = idchild(tbranch(lev-1),tbranch(lev-2))
          call setrclcflglev(tbranch,d1i1,nlevel,lev-1,npnt%c(i)%child)
        end if
      end do
    end if

  end subroutine setrclcflglev

!###################################################################################################################################

  subroutine calctree(nlevel,TLB)
    use fourdgridtypes
    implicit none
    integer,intent(in) :: nlevel
    type(node) :: TLB

    if (TLB%rclcflg) then
      call calclevel(nlevel+1,TLB)
      TLB%rclcflg = .false.
    else
      print*, 'Attempted to recalculate tree - but does not seem necessary'
    end if

!    print*, 'TLB:'
!    print*, 'centroid:', TLB%centroid
!    print*, 'moment:', TLB%mentroid

  end subroutine calctree


!!not weighted at all
!  recursive subroutine calclevel(lev,npnt)
!    use fourdgridtypes
!    implicit none
!    integer,intent(in) :: lev
!    type(node) :: npnt
!    real*8,dimension(3) :: d3r1,d3r2
!    real*8,dimension(4) :: d4r1,d4r2
!    real*8 :: w1,w2,wi,w1t,w2t
!    integer :: i

!    w1 = 0.
!    w2 = 0.
!    w1t = 0.
!    w2t = 0.
!    d3r1 = 0.
!    d3r2 = 0.
!    d4r1 = 0.
!    d4r2 = 0.

!    if (lev > 0) then
!      do i=1,8,1
!        if (associated(npnt%c(i)%child)) then
!          if (npnt%c(i)%child%rclcflg) then
!            call calclevel(lev-1,npnt%c(i)%child)
!            npnt%c(i)%child%rclcflg = .false.
!          end if
!          wi = sqrt(dot_product(npnt%c(i)%child%mentroid,npnt%c(i)%child%mentroid))
!          d4r1(1:3) = d4r1(1:3) + npnt%c(i)%child%centroid(1:3)
!          d4r1(4) = d4r1(4) + npnt%c(i)%child%centroid(4)
!          d3r1 = d3r1 + npnt%c(i)%child%mentroid
!          w1 = w1 + 1.
!          w1t = w1t + 1.
!        end if
!      end do
!      do i=9,16,1
!        if (associated(npnt%c(i)%child)) then
!          if (npnt%c(i)%child%rclcflg) then
!            call calclevel(lev-1,npnt%c(i)%child)
!            npnt%c(i)%child%rclcflg = .false.
!          end if
!          wi = sqrt(dot_product(npnt%c(i)%child%mentroid,npnt%c(i)%child%mentroid))
!          d4r2(1:3) = d4r2(1:3) + npnt%c(i)%child%centroid(1:3)
!          d4r2(4) = d4r2(4) + npnt%c(i)%child%centroid(4)
!          d3r2 = d3r2 + npnt%c(i)%child%mentroid
!          w2 = w2 + 1.
!          w2t = w2t + 1.
!        end if
!      end do
!      npnt%mentroid = 0.
!      npnt%centroid = 0.
!      if (w1 > 0. .and. w2 > 0.) then
!        npnt%mentroid = (d3r1 + d3r2)/2.
!        npnt%centroid(1:3) = (d4r1(1:3)/w1 + d4r2(1:3)/w2)/2.
!        npnt%centroid(4) = (d4r1(4)/w1t + d4r2(4)/w2t)/2.
!      else if (w1 > 0.) then
!        npnt%mentroid = d3r1
!        npnt%centroid(1:3) = d4r1(1:3)/w1
!        npnt%centroid(4) = d4r1(4)/w1t
!      else if (w2 > 0.) then
!        npnt%mentroid = d3r2
!        npnt%centroid(1:3) = d4r2(1:3)/w2
!        npnt%centroid(4) = d4r2(4)/w2t
!      else
!        print*, 'Error in calclevel, level:', lev
!      end if

!    else if (lev == 0) then
!      do i=1,npnt%nvtx,1
!        wi = sqrt(dot_product(npnt%bvtx(i)%vpnt%minfr,npnt%bvtx(i)%vpnt%minfr))
!        d4r1(1:3) = d4r1(1:3) + npnt%bvtx(i)%vpnt%infr(1:3)
!        d4r1(4) = d4r1(4) + npnt%bvtx(i)%vpnt%infr(4) / npnt%nvtx
!        d3r1 = d3r1 + npnt%bvtx(i)%vpnt%minfr
!        w1 = w1 + 1.
!      end do
!      npnt%mentroid = 0.
!      npnt%centroid = 0.
!      if (w1 > 0.) then
!        npnt%mentroid = d3r1
!        npnt%centroid(1:3) = d4r1(1:3)/w1
!        npnt%centroid(4) = d4r1(4)
!      else
!        print*, 'Error in calclevel, level: 0'
!      end if
!    end if

!  end subroutine calclevel


!not weighted in time
  recursive subroutine calclevel(lev,npnt)
    use fourdgridtypes
    implicit none
    integer,intent(in) :: lev
    type(node) :: npnt
    real*8,dimension(3) :: d3r1,d3r2
    real*8,dimension(4) :: d4r1,d4r2
    real*8 :: w1,w2,wi,w1t,w2t
    integer :: i

    w1 = 0.
    w2 = 0.
    w1t = 0.
    w2t = 0.
    d3r1 = 0.
    d3r2 = 0.
    d4r1 = 0.
    d4r2 = 0.

    if (lev > 0) then
      do i=1,8,1
        if (associated(npnt%c(i)%child)) then
          if (npnt%c(i)%child%rclcflg) then
            call calclevel(lev-1,npnt%c(i)%child)
            npnt%c(i)%child%rclcflg = .false.
          end if
          wi = sqrt(dot_product(npnt%c(i)%child%mentroid,npnt%c(i)%child%mentroid))
          d4r1(1:3) = d4r1(1:3) + wi * npnt%c(i)%child%centroid(1:3)
          d4r1(4) = d4r1(4) + npnt%c(i)%child%centroid(4)
          d3r1 = d3r1 + npnt%c(i)%child%mentroid
          w1 = w1 + wi
          w1t = w1t + 1.
        end if
      end do
      do i=9,16,1
        if (associated(npnt%c(i)%child)) then
          if (npnt%c(i)%child%rclcflg) then
            call calclevel(lev-1,npnt%c(i)%child)
            npnt%c(i)%child%rclcflg = .false.
          end if
          wi = sqrt(dot_product(npnt%c(i)%child%mentroid,npnt%c(i)%child%mentroid))
          d4r2(1:3) = d4r2(1:3) + wi * npnt%c(i)%child%centroid(1:3)
          d4r2(4) = d4r2(4) + npnt%c(i)%child%centroid(4)
          d3r2 = d3r2 + npnt%c(i)%child%mentroid
          w2 = w2 + wi
          w2t = w2t + 1.
        end if
      end do
      npnt%mentroid = 0.
      npnt%centroid = 0.
      if (w1 > 0. .and. w2 > 0.) then
        npnt%mentroid = (d3r1 + d3r2)/2.
        npnt%centroid(1:3) = (d4r1(1:3)/w1 + d4r2(1:3)/w2)/2.
        npnt%centroid(4) = (d4r1(4)/w1t + d4r2(4)/w2t)/2.
      else if (w1 > 0.) then
        npnt%mentroid = d3r1
        npnt%centroid(1:3) = d4r1(1:3)/w1
        npnt%centroid(4) = d4r1(4)/w1t
      else if (w2 > 0.) then
        npnt%mentroid = d3r2
        npnt%centroid(1:3) = d4r2(1:3)/w2
        npnt%centroid(4) = d4r2(4)/w2t
      else
        print*, 'Error in calclevel, level:', lev
      end if

    else if (lev == 0) then
      do i=1,npnt%nvtx,1
        wi = sqrt(dot_product(npnt%bvtx(i)%vpnt%minfr,npnt%bvtx(i)%vpnt%minfr))
        d4r1(1:3) = d4r1(1:3) + npnt%bvtx(i)%vpnt%infr(1:3) * wi
        d4r1(4) = d4r1(4) + npnt%bvtx(i)%vpnt%infr(4)
        d3r1 = d3r1 + npnt%bvtx(i)%vpnt%minfr
        w1 = w1 + wi
      end do
      npnt%mentroid = 0.
      npnt%centroid = 0.
      if (w1 > 0.) then
        npnt%mentroid = d3r1
        npnt%centroid(1:3) = d4r1(1:3)/w1
        npnt%centroid(4) = d4r1(4) / npnt%nvtx
      else
        print*, 'Error in calclevel, level: 0'
print*, npnt%nvtx, npnt%mentroid, npnt%centroid
      end if
    end if

  end subroutine calclevel

!weighted in time
!  recursive subroutine calclevel(lev,npnt)
!    use fourdgridtypes
!    implicit none
!    integer,intent(in) :: lev
!    type(node) :: npnt
!    real*8,dimension(3) :: d3r1,d3r2
!    real*8,dimension(4) :: d4r1,d4r2
!    real*8 :: w1,w2,wi
!    integer :: i

!    w1 = 0.
!    w2 = 0.
!    d3r1 = 0.
!    d3r2 = 0.
!    d4r1 = 0.
!    d4r2 = 0.

!    if (lev > 0) then
!      do i=1,8,1
!        if (associated(npnt%c(i)%child)) then
!          if (npnt%c(i)%child%rclcflg) then
!            call calclevel(lev-1,npnt%c(i)%child)
!            npnt%c(i)%child%rclcflg = .false.
!          end if
!          wi = sqrt(dot_product(npnt%c(i)%child%mentroid,npnt%c(i)%child%mentroid))
!          d4r1 = d4r1 + wi * npnt%c(i)%child%centroid
!          d3r1 = d3r1 + npnt%c(i)%child%mentroid
!          w1 = w1 + wi
!        end if
!      end do
!      do i=9,16,1
!        if (associated(npnt%c(i)%child)) then
!          if (npnt%c(i)%child%rclcflg) then
!            call calclevel(lev-1,npnt%c(i)%child)
!            npnt%c(i)%child%rclcflg = .false.
!          end if
!          wi = sqrt(dot_product(npnt%c(i)%child%mentroid,npnt%c(i)%child%mentroid))
!          d4r2 = d4r2 + wi * npnt%c(i)%child%centroid
!          d3r2 = d3r2 + npnt%c(i)%child%mentroid
!          w2 = w2 + wi
!        end if
!      end do
!      npnt%mentroid = 0.
!      npnt%centroid = 0.
!      if (w1 > 0. .and. w2 > 0.) then
!        npnt%mentroid = (d3r1 + d3r2)/2.
!        npnt%centroid = (d4r1/w1 + d4r2/w2)/2.
!      else if (w1 > 0.) then
!        npnt%mentroid = d3r1
!        npnt%centroid = d4r1/w1
!      else if (w2 > 0.) then
!        npnt%mentroid = d3r2
!        npnt%centroid = d4r2/w2
!      else
!        print*, 'Error in calclevel, level:', lev
!      end if

!    else if (lev == 0) then
!      do i=1,npnt%nvtx,1
!        wi = sqrt(dot_product(npnt%bvtx(i)%vpnt%minfr,npnt%bvtx(i)%vpnt%minfr))
!        d4r1 = d4r1 + npnt%bvtx(i)%vpnt%infr * wi
!        d3r1 = d3r1 + npnt%bvtx(i)%vpnt%minfr
!        w1 = w1 + wi
!      end do
!      npnt%mentroid = 0.
!      npnt%centroid = 0.
!      if (w1 > 0.) then
!        npnt%mentroid = d3r1
!        npnt%centroid = d4r1/w1
!      else
!        print*, 'Error in calclevel, level: 0'
!      end if
!    end if

!  end subroutine calclevel

!###################################################################################################################################
!###################################################################################################################################

!  subroutine checktree(tsteps,nlevel,TLB,deltat,branch,buffer,slev)
!    use constants
!    use fourdgridtypes
!    implicit none
!    integer,intent(in) :: tsteps,nlevel,slev
!    type(node) :: TLB
!    real*8,intent(in) :: deltat
!    integer,dimension(2,3,0:nlevel),intent(in) :: buffer
!    integer,dimension(4,0:nlevel),intent(in) :: branch
!    integer :: d1i1,i,j,k,ibit
!    integer,dimension(4) :: d4i1

!    d1i1 = TLBchild(branch(:,nlevel))



!    !check right set of time octree set
!    ibit = int(real(d1i1-0.5)/8.)
!    !nearfield/farfield
!    do i=1+8*ibit,8+8*ibit,1
!      d4i1 = i4childTLB(i)
!      if (inoutNF(nlevel,nlevel,d4i1,buffer)) then
!        if (
!        call checklevel(tsteps,nlevel,npnt,lev,branch,buffer,slev)
!!      else
!!        call calcbox(nlevel,TLB,nlevel,d4i1,branch,deltat,r4in,r3inci,v,vind,gradind,lflag,sflag)
!      end if
!    end do



!  contains

!    recursive subroutine checklevel




!      if (lev > slev) then


!      else


!    end subroutine checklevel

!  end subroutine checktree

!###################################################################################################################################
!###################################################################################################################################

  subroutine calcinfl(it,jp,tsteps,nlevel,TLB,deltat,r4in,r3inci,v,vind,gradind,lflag,sflag,fvtx,cfvtx)
    use constants
    use fourdgridtypes
    implicit none
    integer,intent(in) :: tsteps,nlevel
    type(node) :: TLB
    real*8,intent(in) :: deltat
    real*8,dimension(4),intent(in) :: r4in
    real*8,dimension(3),intent(in) :: r3inci,v
    logical,intent(in) :: lflag,sflag
    real*8,dimension(3) :: vind,d3r
    real*8,dimension(3,3) :: gradind,d33r
    integer :: d1i1,i,j,k,ibit
    integer,dimension(4) :: d4i1
    integer,dimension(2,3,0:nlevel) :: buffer
    integer,dimension(4,0:nlevel) :: branch,vbranch
    integer :: fvtx
    logical,dimension(fvtx) :: cfvtx

    integer,dimension(0:nlevel) :: zerobuff


integer :: it,jp

    vind = 0.
    gradind = 0.

    zerobuff(0) = 2
    zerobuff(1:3) = 2
    zerobuff(3:nlevel) = 0

    branch(:,0) = box4d(deltat,r4in)
    do i=1,3,1
      buffer(1,i,0) = branch(i,0) - zerobuff(0)
      if (buffer(1,i,0)*branch(i,0) <= 0) then
        buffer(1,i,0) = buffer(1,i,0) - 1
      end if
      buffer(2,i,0) = branch(i,0) + zerobuff(0)
      if (buffer(2,i,0)*branch(i,0) <= 0) then
        buffer(2,i,0) = buffer(2,i,0) + 1
      end if
    end do
!    if (zerobuff == 0) then
!      buffer(1,:,0) = branch(1:3,0)
!      buffer(2,:,0) = branch(1:3,0)
!    end if

!    do i=1,nlevel,1
!      branch(:,i) = i4parent(branch(:,i-1))
!      do j=1,3,1
!        buffer(1,j,i) = up(buffer(1,j,i-1)) - zerobuff(i)
!        if (buffer(1,j,i)*up(buffer(1,j,i-1)) <= 0) then
!          buffer(1,j,i) = buffer(1,j,i) - 1
!        end if
!        buffer(2,j,i) = up(buffer(2,j,i-1)) + zerobuff(i)
!        if (buffer(2,j,i)*up(buffer(2,j,i-1)) <= 0) then
!          buffer(2,j,i) = buffer(2,j,i) + 1
!        end if
!      end do

    do i=1,nlevel,1
      branch(:,i) = i4parent(branch(:,i-1))
      do j=1,3,1
        buffer(1,j,i) = branch(j,i) - zerobuff(i)
        if (buffer(1,j,i)*branch(j,i) <= 0) then
          buffer(1,j,i) = buffer(1,j,i) - 1
        end if
        if (buffer(1,j,i) > up(buffer(1,j,i-1))) then
          buffer(1,j,i) = up(buffer(1,j,i-1))
        end if

        buffer(2,j,i) = branch(j,i) + zerobuff(i)
        if (buffer(2,j,i)*branch(j,i) <= 0) then
          buffer(2,j,i) = buffer(2,j,i) + 1
        end if
        if (buffer(2,j,i) < up(buffer(2,j,i-1))) then
          buffer(2,j,i) = up(buffer(2,j,i-1))
        end if
      end do
    end do


    !determine which child particle is in at TLB level
    d1i1 = TLBchild(branch(:,nlevel))


    !check right set of time octree set
    ibit = int(real(d1i1-0.5)/8.)
    !nearfield/farfield
    do i=1+8*ibit,8+8*ibit,1
      if (i /= d1i1) then
        d4i1 = i4childTLB(i)
        if (inoutNF(nlevel,nlevel,d4i1,buffer)) then
          call calcboxF(it,jp,tsteps,nlevel,TLB,nlevel,d4i1,branch,buffer,deltat,r4in,r3inci,v,vind,gradind,lflag,sflag,fvtx,cfvtx)
        else
!          call calcTLB(nlevel,TLB,i,deltat,r4in,r3inci,v,vind,gradind,lflag,sflag)
          call calcbox(nlevel,TLB,nlevel,d4i1,branch,deltat,r4in,r3inci,v,vind,gradind,lflag,sflag)
        end if
      end if
    end do

    if (associated(TLB%c(d1i1)%child)) then
      call calcboxown(it,jp,tsteps,nlevel,TLB,nlevel,TLB%c(d1i1)%child,branch,buffer,deltat,r4in,r3inci,v,vind,gradind, &
         & lflag,sflag,fvtx,cfvtx)
    end if


!    contains


!    subroutine calcTLB(nlevel,TLB,ichld,deltat,r4in,r3inci,v,vind,gradind,lflag,sflag)
!      use vtxparticle
!      use constants
!      use fourdgridtypes
!      use pg
!      implicit none
!      integer,intent(in) :: nlevel,ichld
!      type(node) :: TLB
!      real*8,intent(in) :: deltat
!      real*8,dimension(4),intent(in) :: r4in
!      real*8,dimension(3),intent(in) :: r3inci,v
!      logical,intent(in) :: lflag,sflag
!      real*8,dimension(3) :: vind,d3r
!      real*8,dimension(3,3) :: gradind,d33r
!      integer :: d1i2,j
!      real*8,dimension(4) :: d4r1
!      integer,dimension(4,0:nlevel) :: vbranch
!      type(node),pointer,save :: vnpnt
!      real*8 :: rdt
!print*, 'here'

!      d4r1 = boxcentroid(deltat,nlevel,i4childTLB(ichld))
!      rdt = sqrt(dot_product(d4r1(1:3)-r4in(1:3),d4r1(1:3)-r4in(1:3)))/(a*deltat)
!  !    rdt = abs(d4r1(1)-r4in(1))/(a*deltat)
!      d4r1(4) = r4in(4) - rdt
!      if (d4r1(4) .ge. 1.) then
!        vbranch(:,0) = box4d(deltat,d4r1)
!        !find all boxes
!        do j=1,nlevel,1
!          vbranch(:,j) = i4parent(vbranch(:,j-1))
!        end do
!        d1i2 = TLBchild(vbranch(:,nlevel))
!        if (associated(TLB%c(d1i2)%child)) then
!          if (d1i2 /= d1i1) then
!            vnpnt => TLB%c(d1i2)%child

!!should interpolate here between TLB child boxes


!!            if (sflag) then
!!              call vtxinduced(scaler(a,v,r4in(1:3)),scaler(a,v,vnpnt%centroid(1:3)),vnpnt%mentroid,d3r)
!!              vind = vind + d3r
!!              call shearind(scaler(a,v,r4in(1:3)),scaler(a,v,vnpnt%centroid(1:3)),vnpnt%mentroid,d33r)
!!              gradind = gradind + d33r
!!            else
!              call vtxinduced(r3inci,vnpnt%centroid(1:3),vnpnt%mentroid,d3r)
!              vind = vind + d3r
!              call shearind(r3inci,vnpnt%centroid(1:3),vnpnt%mentroid,d33r)
!              gradind = gradind + d33r
!!            end if
!            if (lflag) then
!              vnpnt%usdflg = .true.
!            end if
!          end if
!        end if
!      end if
!    end subroutine calcTLB  

  end subroutine calcinfl

!###################################################################################################################################

!npnt and plev are at same level
 recursive subroutine calcboxown(it,jp,tsteps,nlevel,TLB,plev,npnt,branch,buffer,deltat,r4in,r3inci,v,vind,gradind, &
    & lflag,sflag,fvtx,cfvtx)
    use constants
    use fourdgridtypes
    implicit none
    integer,intent(in) :: tsteps,nlevel,plev
    type(node) :: TLB,npnt
    integer,dimension(4,0:nlevel),intent(in) :: branch
    integer,dimension(2,3,0:nlevel),intent(in) :: buffer
    real*8,intent(in) :: deltat
    real*8,dimension(4),intent(in) :: r4in
    real*8,dimension(3),intent(in) :: r3inci,v
    logical,intent(in) :: lflag,sflag
    real*8,dimension(3) :: vind
    real*8,dimension(3,3) :: gradind  
    integer :: d1i1,i,ibit
    integer,dimension(4) :: d4i1
    integer :: fvtx
    logical,dimension(fvtx) :: cfvtx

integer :: it,jp

    d1i1 = nchild(branch(:,plev),branch(:,plev-1))

    !check right set of time octree set%n 
    ibit = int(real(d1i1-0.5)/8.)

    !nearest time
    do i=1+8*ibit,8+8*ibit,1
      if (i /= d1i1) then
        d4i1 = i4child(branch(:,plev),i)
        if (plev > 1) then
          if (inoutNF(nlevel,plev-1,d4i1,buffer)) then
            call calcboxF(it,jp,tsteps,nlevel,TLB,plev-1,d4i1,branch,buffer,deltat,r4in,r3inci,v,vind,gradind, &
              & lflag,sflag,fvtx,cfvtx)
          else
            call calcbox(nlevel,TLB,plev-1,d4i1,branch,deltat,r4in,r3inci,v,vind,gradind,lflag,sflag)
          end if
        else
          if (inoutNF(nlevel,plev-1,d4i1,buffer)) then
            call calc1on1(it,jp,tsteps,nlevel,TLB,plev-1,d4i1,branch,deltat,r4in,r3inci,v,vind,gradind,lflag,sflag,fvtx,cfvtx)
          else
            call calcbox(nlevel,TLB,plev-1,d4i1,branch,deltat,r4in,r3inci,v,vind,gradind,lflag,sflag)
          end if
        end if
      end if
    end do


    if (plev>1) then
      if (associated(npnt%c(d1i1)%child)) then
        call calcboxown(it,jp,tsteps,nlevel,TLB,plev-1,npnt%c(d1i1)%child,branch,buffer,deltat,r4in,r3inci,v,vind,gradind, &
           & lflag,sflag,fvtx,cfvtx)
      end if
    else
      if (associated(npnt%c(d1i1)%child)) then
        call calcbvtx(it,jp,tsteps,deltat,npnt%c(d1i1)%child,r4in,r3inci,v,vind,gradind,lflag,sflag,fvtx,cfvtx)
      end if
    end if

  end subroutine calcboxown

!###################################################################################################################################

  recursive subroutine calcboxF(it,jp,tsteps,nlevel,TLB,plev,p4i,branch,buffer,deltat,r4in,r3inci,v,vind,gradind, &
    & lflag,sflag,fvtx,cfvtx)
    use fourdgridtypes
    implicit none
    integer,intent(in) :: tsteps,nlevel,plev
    type(node) :: TLB
    integer,dimension(4),intent(in) :: p4i
    integer,dimension(4,0:nlevel),intent(in) :: branch
    integer,dimension(2,3,0:nlevel),intent(in) :: buffer
    real*8,intent(in) :: deltat
    real*8,dimension(4),intent(in) :: r4in
    real*8,dimension(3),intent(in) :: r3inci,v
    logical,intent(in) :: lflag,sflag
    real*8,dimension(3) :: vind
    real*8,dimension(3,3) :: gradind  
    integer :: d1i1,i,ibit
    integer,dimension(4) :: d4i1
    integer :: fvtx
    logical,dimension(fvtx) :: cfvtx

integer :: it,jp

    !only checking for time here
    d1i1 = nchild(branch(:,plev),branch(:,plev-1))

    !check right set of time octree set
    ibit = int(real(d1i1-0.5)/8.)

    !nearest time


    if (plev > 1) then
      do i=1+8*ibit,8+8*ibit,1
        d4i1 = i4child(p4i,i)
        if (inoutNF(nlevel,plev-1,d4i1,buffer)) then
          call calcboxF(it,jp,tsteps,nlevel,TLB,plev-1,d4i1,branch,buffer,deltat,r4in,r3inci,v,vind,gradind,lflag,sflag,fvtx,cfvtx)
        else
          call calcbox(nlevel,TLB,plev-1,d4i1,branch,deltat,r4in,r3inci,v,vind,gradind,lflag,sflag)
        end if
      end do
    else
      do i=1+8*ibit,8+8*ibit,1
        d4i1 = i4child(p4i,i)
        if (inoutNF(nlevel,plev-1,d4i1,buffer)) then
          !calc 1:1
          call calc1on1(it,jp,tsteps,nlevel,TLB,plev-1,d4i1,branch,deltat,r4in,r3inci,v,vind,gradind,lflag,sflag,fvtx,cfvtx)
        else
          !calc centroid
          call calcbox(nlevel,TLB,plev-1,d4i1,branch,deltat,r4in,r3inci,v,vind,gradind,lflag,sflag)
        end if
      end do
    end if

  end subroutine calcboxF

!###################################################################################################################################

  subroutine calc1on1(it,jp,tsteps,nlevel,TLB,lev,d4i1,branch,deltat,r4in,r3inci,v,vind,gradind,lflag,sflag,fvtx,cfvtx)
    use fourdgridtypes
    implicit none
    integer,intent(in) :: tsteps,nlevel,lev
    type(node) :: TLB
    integer,dimension(4),intent(in) :: d4i1
    real*8,intent(in) :: deltat
    real*8,dimension(4),intent(in) :: r4in
    real*8,dimension(3),intent(in) :: r3inci,v
    integer,dimension(4,0:nlevel),intent(in) :: branch
    real*8,dimension(3) :: vind,d3r1,d3r2
    real*8,dimension(3,3) :: gradind,d33r1,d33r2
    logical,intent(in) :: lflag,sflag
    type(node),pointer,save :: vnpnt1,vnpnt2
    real*8 :: fr1,fr2
    integer,intent(in) :: fvtx
    logical,dimension(fvtx) :: cfvtx
    integer,dimension(4) :: d4i2
    integer,dimension(4,0:nlevel) :: vbranch1
    integer :: d1i1,i

integer :: it,jp

    vnpnt1 => null()
!    vnpnt2 => null()

!    call calcfind_interp(nlevel,TLB,lev,d4i1,deltat,r4in,branch,vnpnt1,vnpnt2,fr1,fr2)
!  
!    !if found then calc

!    if (associated(vnpnt1)) then
!      if (associated(vnpnt2)) then
!        d3r1 = 0.
!        d33r1 = 0.
!        d3r2 = 0.
!        d33r2 = 0.
!!Need to check this.... duplication might be a possibility here
!        call calcbvtx(tsteps,deltat,vnpnt1,r4in,r3inci,v,d3r1,d33r1,lflag,sflag,fvtx,cfvtx)
!        call calcbvtx(tsteps,deltat,vnpnt2,r4in,r3inci,v,d3r2,d33r2,lflag,sflag,fvtx,cfvtx)
!        vind = vind + fr1*d3r1 + fr2*d3r2
!        gradind = gradind + fr1*d33r1 + fr2*d33r2
!      else
!        !calc 1on1 of vnpnt1
!        call calcbvtx(tsteps,deltat,vnpnt1,r4in,r3inci,v,vind,gradind,lflag,sflag,fvtx,cfvtx)
!      end if
!    end if




    !find npnt

    d4i2(1:3) = d4i1(1:3)
!    d4i2(4) = d4i1(4) - 1
    d4i2(4) = box1dt(deltat,r4in(4)) !int(r4in(4)) ! d4i1(4)

    vbranch1(:,0) = d4i2
    do i=1,nlevel,1
      vbranch1(:,i) = i4parent(vbranch1(:,i-1))
    end do

    d1i1 = TLBchild(vbranch1(:,nlevel))
    if (associated(TLB%c(d1i1)%child)) then
      vnpnt1 => TLB%c(d1i1)%child
      do i=nlevel-1,lev,-1
        d1i1 = nchild(vbranch1(:,i+1),vbranch1(:,i))
        if (associated(vnpnt1%c(d1i1)%child)) then
          vnpnt1 => vnpnt1%c(d1i1)%child
        else
          vnpnt1 => null()
          exit
        end if
      end do
    else
      vnpnt1 => null()
    end if

    if (associated(vnpnt1)) then
!if (jp == 351 .and. it == 10) then
!print*
!print*, lev,vbranch1(:,0)
!print*, d4i2
!end if
      call calcbvtx(it,jp,tsteps,deltat,vnpnt1,r4in,r3inci,v,vind,gradind,lflag,sflag,fvtx,cfvtx)
    end if


  end subroutine calc1on1

!###################################################################################################################################

    subroutine calcbvtx(it,jp,tsteps,deltat,npnt,r4in,r3inci,v,vind,gradind,lflag,sflag,fvtx,cfvtx)
      use vtxparticle
      use basictypes
      use constants
      use fourdgridtypes
      use pg
      implicit none
      integer,intent(in) :: tsteps
      type(node) :: npnt
      real*8,intent(in) :: deltat
      real*8,dimension(4),intent(in) :: r4in
      real*8,dimension(3),intent(in) :: r3inci,v
      logical,intent(in) :: lflag,sflag
      real*8,dimension(3) :: vind,d3r
      real*8,dimension(3,3) :: gradind,d33r
      type(point),pointer,save :: vtxset
      type(rmset),pointer,save :: vtxpnt
      integer :: i,d1i1,d1i2,istop
      real*8 :: d1r1,d1r2,rdt
      real*8,dimension(3) :: d3r1,d3r2
      logical :: ffl
      integer :: fvtx
      logical,dimension(fvtx) :: cfvtx

integer :: it,jp

      do i=1,npnt%nvtx,1
        vtxset => npnt%bvtx(i)%vset
        vtxpnt => npnt%bvtx(i)%vpnt
        if (vtxset%typ == 'V') then
          istop = (tsteps-size(vtxset%set)) + 1
          call interppnt(r4in(1:3),tsteps,int(r4in(4)),istop,a,deltat,vtxset,d1i1,d1i2,d1r1,d1r2,ffl)

          if (sflag) then
            if (ffl) then
              call ittodummy(r4in(1:3),vtxset%set(int(r4in(4)))%infr(1:3),v,deltat,d3r1)
            else
              d3r1 = intp(vtxset%set(d1i1)%infr(1:3),vtxset%set(d1i2)%infr(1:3),d1r1,d1r2)
            end if
            d3r2 = intp(vtxset%set(d1i1)%minfr,vtxset%set(d1i2)%minfr,d1r1,d1r2)
            call vtxinduced(r3inci,d3r1,d3r2,d3r)
            vind = vind + d3r
            call shearind(r3inci,d3r1,d3r2,d33r)
            gradind = gradind + d33r
          else
            if (ffl .eqv. .false.) then
              d3r1 = intp(vtxset%set(d1i1)%infr(1:3),vtxset%set(d1i2)%infr(1:3),d1r1,d1r2)
              d3r2 = intp(vtxset%set(d1i1)%minfr,vtxset%set(d1i2)%minfr,d1r1,d1r2)
              call vtxinduced(r3inci,d3r1,d3r2,d3r)
              vind = vind + d3r
              call shearind(r3inci,d3r1,d3r2,d33r)
              gradind = gradind + d33r
            end if
          end if

        else if (vtxset%typ == 'W') then
          cfvtx(vtxset%n) = .true.
        else if (vtxset%typ == 'K') then
          cfvtx(vtxset%n) = .true.
        else
          print*, 'No Vtx particle type specified'
        end if
      end do


      if (lflag) then
        npnt%usdflg = .true.
      end if

    end subroutine calcbvtx

!###################################################################################################################################

!    subroutine calcbvtx(tsteps,deltat,npnt,r4in,v,vind,gradind,lflag,sflag)
!      use vtxparticle
!      use basictypes
!      use constants
!      use fourdgridtypes
!      use pg
!      implicit none
!      integer,intent(in) :: tsteps
!      type(node) :: npnt
!      real*8,intent(in) :: deltat
!      real*8,dimension(4),intent(in) :: r4in
!      real*8,dimension(3),intent(in) :: v
!      logical,intent(in) :: lflag,sflag
!      real*8,dimension(3) :: vind,d3r
!      real*8,dimension(3,3) :: gradind,d33r
!      type(point),pointer,save :: vtxset
!      type(rmset),pointer,save :: vtxpnt
!      integer :: i,d1i1,d1i2,istop
!      real*8 :: d1r1,d1r2,rdt
!      real*8,dimension(3) :: d3r1,d3r2
!      logical :: ffl

!      do i=1,npnt%nvtx,1
!        vtxset => npnt%bvtx(i)%vset
!        vtxpnt => npnt%bvtx(i)%vpnt
!        istop = (tsteps-size(vtxset%set)) + 1
!        call interppnt(r4in(1:3),tsteps,int(r4in(4)),istop,a,deltat,vtxset,d1i1,d1i2,d1r1,d1r2,ffl)

!        if (sflag) then
!          if (ffl) then
!            call ittodummy(r4in(1:3),vtxset%set(r4in(4))%infr(1:3),v,deltat,d3r1)
!          else
!            d3r1 = intp(vtxset%set(d1i1)%infr(1:3),vtxset%set(d1i2)%infr(1:3),d1r1,d1r2)
!          end if
!          d3r2 = intp(vtxset%set(d1i1)%minfr,vtxset%set(d1i2)%minfr,d1r1,d1r2)
!          call vtxinduced(r4in(1:3),d3r1,d3r2,d3r)
!          vind = vind + d3r
!          call shearind(r4in(1:3),d3r1,d3r2,d33r)
!          gradind = gradind + d33r
!        else
!          if (ffl .eqv. .false.) then
!            d3r1 = intp(vtxset%set(d1i1)%infr(1:3),vtxset%set(d1i2)%infr(1:3),d1r1,d1r2)
!            d3r2 = intp(vtxset%set(d1i1)%minfr,vtxset%set(d1i2)%minfr,d1r1,d1r2)
!            call vtxinduced(r4in(1:3),d3r1,d3r2,d3r)
!            vind = vind + d3r
!            call shearind(r4in(1:3),d3r1,d3r2,d33r)
!            gradind = gradind + d33r
!          end if
!        end if
!      end do
!      if (lflag) then
!        npnt%usdflg = .true.
!      end if
!    end subroutine calcbvtx

!###################################################################################################################################

  subroutine calcbox(nlevel,TLB,lev,d4i1,branch,deltat,r4in,r3inci,v,vind,gradind,lflag,sflag)
    use fourdgridtypes
    use infl
    use vtxparticle
    use pg
    use constants
    implicit none
    integer,intent(in) :: nlevel,lev
    type(node) :: TLB
    integer,dimension(4),intent(in) :: d4i1
    real*8,intent(in) :: deltat
    real*8,dimension(4),intent(in) :: r4in
    real*8,dimension(3),intent(in) :: r3inci,v
    integer,dimension(4,0:nlevel),intent(in) :: branch
    real*8,dimension(3) :: vind,d3r1,d3r2
    real*8,dimension(3,3) :: gradind,d33r1,d33r2
    logical,intent(in) :: lflag,sflag
    type(node),pointer,save :: vnpnt1,vnpnt2
    real*8 :: fr1,fr2,rdt,d1r1,d1r2
    integer,dimension(4,0:nlevel) :: vbranch1
    integer :: d1i1,d1i2,j
    real*8,dimension(4) :: d4r1
    integer,dimension(4) :: d4i2
    real*8,dimension(3) :: d3r,zcout
    real*8,dimension(3,3) :: d33r
    real*8,dimension(3),parameter :: zerovec = 0.
    real*8,dimension(3),parameter :: vinf = (/306,0,0/)

    vnpnt1 => null()
    vnpnt2 => null()

!if (lev > 5) print*, 'here'

    call calcfind_interp(nlevel,TLB,lev,d4i1,deltat,r4in,branch,vnpnt1,vnpnt2,zcout,fr1,fr2)

    !if found then calc
    if (associated(vnpnt1)) then
      if (associated(vnpnt2)) then


        d3r1 = intp(vnpnt1%centroid(1:3),vnpnt2%centroid(1:3),fr1,fr2)
        d3r2 = intp(vnpnt1%mentroid,vnpnt2%mentroid,fr1,fr2)
 !       d3r1 = d3r1 !- 0.5 * 2**lev * deltat * v !* sign(dble(1.),dot_product(vnpnt1%centroid(1:3)-r4in(1:3),v))
        call calcvpoint(r3inci,d3r1,d3r2,d3r,d33r)
        vind = vind + d3r
        gradind = gradind + d33r



!        d1r1 = max(1.-(abs(d1r1)/(2**lev)),0.)
!        d1r2 = max(1.-(abs(d1r2)/(2**lev)),0.)
!        d1r1 = abs(fr2)/(abs(fr1)+abs(fr2))
!        d1r2 = abs(fr1)/(abs(fr1)+abs(fr2))

!        call calcvpoint(r3inci,vnpnt1%centroid(1:3),vnpnt1%mentroid,d3r,d33r)
!        vind = vind + d1r1*d3r
!        gradind = gradind + d1r1*d33r
!        call calcvpoint(r3inci,vnpnt2%centroid(1:3),vnpnt2%mentroid,d3r,d33r)
!        vind = vind + d1r2*d3r
!        gradind = gradind + d1r2*d33r




!        d1r1 = vnpnt1%centroid(4) - (r4in(4) - &
!             & sqrt(dot_product(vnpnt1%centroid(1:3)-r4in(1:3),vnpnt1%centroid(1:3)-r4in(1:3)))/(a*deltat))
!        d3r1 = vnpnt1%centroid(1:3) - vinf*d1r1*deltat
!        d1r2 = vnpnt2%centroid(4) - (r4in(4) - &
!             & sqrt(dot_product(vnpnt2%centroid(1:3)-r4in(1:3),vnpnt2%centroid(1:3)-r4in(1:3)))/(a*deltat))
!        d3r2 = vnpnt2%centroid(1:3) - vinf*d1r2*deltat

!        d1r1 = abs(fr2)/(abs(fr1)+abs(fr2))
!        d1r2 = abs(fr1)/(abs(fr1)+abs(fr2))

!        call calcvpoint(r3inci,d3r1,vnpnt1%mentroid,d3r,d33r)
!        vind = vind + d1r1*d3r
!        gradind = gradind + d1r1*d33r
!        call calcvpoint(r3inci,d3r2,vnpnt2%mentroid,d3r,d33r)
!        vind = vind + d1r2*d3r
!        gradind = gradind + d1r2*d33r


        if (lflag) then
          vnpnt1%usdflg = .true.
          vnpnt2%usdflg = .true.
        end if

      else

!        d4r1 = boxcentroid(deltat,lev,d4i1)


!        call rightime(deltat,r4in,vnpnt1%centroid,vinf,d3r1)
!print*, vnpnt1%centroid
!print*, d3r1
!print*


        d3r1 = vnpnt1%centroid(1:3) !- 0.5 * 2**lev * deltat * v !* sign(dble(1.),dot_product(vnpnt1%centroid(1:3)-r4in(1:3),v))
        call calcvpoint(r3inci,d3r1,vnpnt1%mentroid,d3r,d33r)
        vind = vind + d3r
        gradind = gradind + d33r

!        d1r1 = max(1.-(abs(fr1)/(2**lev)),0.)
!!     d1r1 = abs((vnpnt1%centroid(4) - (r4in(4) - &
!!          & sqrt(dot_product(vnpnt1%centroid(1:3)-r4in(1:3),vnpnt1%centroid(1:3)-r4in(1:3)))/(a*deltat)))/(2**lev))
!        vind = vind + d1r1*d3r
!        gradind = gradind + d1r1*d33r

        if (lflag) then
          vnpnt1%usdflg = .true.
        end if



      end if
    else if (sflag) then
      !if present time box non-zero, calc this interpolated back
      vnpnt1 => null()
      !calc centroid coords in 4D at present time
      d4r1 = boxcentroid(deltat,lev,d4i1)
!      d4r1(4) = r4in(4) ! - 1 !one tstep earlier because present t not in tree yet.
      !find tree for centroid at right time

      call findbox(nlevel,deltat,TLB,lev,d4r1,vnpnt1,d4i2)

      if (associated(vnpnt1)) then
        call calccentroid(deltat,lev,d4i2,vnpnt1,r4in,r3inci,v,vind,gradind,lflag,sflag)
      end if
    end if


  end subroutine calcbox


  subroutine rightime(deltat,r4in,e4r1,vrel,d3rout)
    use constants
    implicit none
    real*8,intent(in) :: deltat
    real*8,dimension(4),intent(in) :: r4in,e4r1
    real*8,dimension(3),intent(in) :: vrel
    real*8,dimension(3),intent(out) :: d3rout
    real*8 :: rdt,d1r1,d1r2,tmin,tmid,tmax
    real*8,dimension(4) :: d4r1,d4r2
    integer :: i

    rdt = sqrt(dot_product(r4in(1:3)-e4r1(1:3),r4in(1:3)-e4r1(1:3)))/(a*deltat)
    
    d1r1 = e4r1(4)-(r4in(4)-rdt)

    tmin = 0.
    tmax = d1r1
    tmid = (tmax-tmin)/2.
    d4r1 = e4r1

    do i=1,100,1
      d4r2(1:3) = e4r1(1:3) + deltat*vrel*tmid
      d4r2(4) = e4r1(4) + tmid
      rdt = sqrt(dot_product(r4in(1:3)-d4r2(1:3),r4in(1:3)-d4r2(1:3)))/(a*deltat)
      d1r2 = d4r2(4) - (r4in(4)-rdt)
      if (d1r1*d1r2 > 0) then
        tmin = tmid
!        tmax = tmax
        tmid = (tmax-tmin)/2.
        d4r1 = d4r2
        d1r1 = d1r2
      else
!        tmin = tmin
        tmax = tmid
        tmid = (tmax-tmin)/2.
        d4r1 = d4r2
        d1r1 = d1r2
      end if
    end do
    d3rout = (d4r1(1:3)+d4r2(1:3))/2.

  end subroutine rightime

!###################################################################################################################################

  subroutine calccentroid(deltat,lev,d4i1,npnt,r4in,r3inci,v,vind,gradind,lflag,sflag)
    use vtxparticle
    use constants
    use fourdgridtypes
    use pg
    implicit none
    type(node) :: npnt
    real*8,dimension(4),intent(in) :: r4in
    integer,dimension(4),intent(in) :: d4i1
    real*8,dimension(3),intent(in) :: r3inci,v
    logical,intent(in) :: lflag,sflag
    real*8,dimension(3) :: vind,d3r,d3r1,d3r2
    real*8,dimension(3,3) :: gradind,d33r
    real*8,intent(in) :: deltat
    integer,intent(in) :: lev
    real*8,dimension(4) :: d4r1

    if (sflag) then
      d4r1 = boxcentroid(deltat,lev,d4i1)
      d3r1 = npnt%centroid(1:3) !- (r4in(4) - npnt%centroid(4)) * v * deltat !* sign(dble(1.),dot_product(npnt%centroid(1:3)-r4in(1:3),v)) ! + 0.5 * 2**lev * deltat * v !* sign(dble(1.),dot_product(npnt%centroid(1:3)-r4in(1:3),v))
      call ittodummy(r4in(1:3),d3r1,v,deltat,d3r2)
      !d3r2 = d3r2 !+ 0.5 * 2**lev * deltat * v !* sign(dble(1.),dot_product(npnt%centroid(1:3)-r4in(1:3),v))

      call calcvpoint(r3inci,d3r2,npnt%mentroid,d3r,d33r)
      vind = vind + d3r
      gradind = gradind + d33r
    else
      call calcvpoint(r3inci,npnt%centroid(1:3),npnt%mentroid,d3r,d33r)
      vind = vind + d3r
      gradind = gradind + d33r
    end if
    if (lflag) then
      npnt%usdflg = .true.
    end if
  end subroutine calccentroid

  subroutine calcvpoint(pcpt,centroid,mentroid,vout,gradout)
    use vtxparticle
    implicit none
    real*8,dimension(3),intent(in) :: pcpt,centroid,mentroid
    real*8,dimension(3),intent(out) :: vout
    real*8,dimension(3,3),intent(out) :: gradout

    call vtxinduced_wing(pcpt,centroid,mentroid,vout)
    call shearind(pcpt,centroid,mentroid,gradout)

  end subroutine calcvpoint

!###################################################################################################################################

  subroutine calcfind_interp(nlevel,TLB,lev,d4i1,deltat,r4in,branch,vnpnt1,vnpnt2,zcout,fr1,fr2)
    use constants
    use fourdgridtypes
    implicit none
    integer,intent(in) :: nlevel,lev
    type(node) :: TLB
    integer,dimension(4),intent(in) :: d4i1
    real*8,intent(in) :: deltat
    real*8,dimension(4),intent(in) :: r4in
    integer,dimension(4,0:nlevel),intent(in) :: branch
    type(node),pointer :: vnpnt1,vnpnt2
    integer :: d1i1,d1i2,j,k
    real*8,dimension(4) :: d4r1,d4r2,d4r3
    integer,dimension(4) :: d4i2,d4i3,n4i1,n4i2
    integer,dimension(4,0:nlevel) :: vbranch1,vbranch2
    real*8 :: rdt
    real*8 :: d1r1,d1r2,fr1,fr2
    real*8,dimension(3),intent(out) :: zcout

    real(kind=8) :: dt1,dt2,as,bs,cs,ls

    vnpnt1 => null()
    vnpnt2 => null()
    !box integer parameters fixed in space dimensions, time is free
    !calc centroid coords in 4D at right time
    d4r1 = boxcentroid(deltat,lev,d4i1)
    rdt = sqrt(dot_product(d4r1(1:3)-r4in(1:3),d4r1(1:3)-r4in(1:3)))/(a*deltat)
    d4r1(4) = r4in(4) - rdt
    !d4r1 is the boxcentroid shifted to the right time so emitter has influence on receiver.
    if (d4r1(4) .gt. 0.) then
      call findbox(nlevel,deltat,TLB,lev,d4r1,vnpnt1,n4i1)
    end if

    do k=1,10,1 !some number to limit intterations
      if (associated(vnpnt1)) then
        d4r1 = vnpnt1%centroid
        rdt = sqrt(dot_product(d4r1(1:3)-r4in(1:3),d4r1(1:3)-r4in(1:3)))/(a*deltat)
        d1r1 = d4r1(4) - (r4in(4) - rdt)
        d4i2(1:3) = n4i1(1:3)
        if (d1r1 > 0.) then
          if (n4i1(4) == 1) then
            d4i2(4) = n4i1(4) - 2
          else
            d4i2(4) = n4i1(4) - 1
          end if
        else if (d1r1 < 0.) then
          if (vbranch1(4,lev) == -1) then
            d4i2(4) = n4i1(4) + 2
          else
            d4i2(4) = n4i1(4) + 1
          end if
        else
           d4i2(4) = n4i1(4)
        end if

        if (d4i2(4) /= n4i1(4)) then
          d4r1 = boxcentroid(deltat,lev,d4i2)
          if (d4r1(4) .gt. 0.) then
            call findbox(nlevel,deltat,TLB,lev,d4r1,vnpnt2,n4i2)
          end if

!        dt1 = dble(t1)
!        dt2 = dble(t2)
!        as = dot_product(vpoint%set(t1)%infr(1:3)-vpoint%set(t2)%infr(1:3),vpoint%set(t1)%infr(1:3)-vpoint%set(t2)%infr(1:3))/ &
!           & (a*deltat)**2 - (dt1-dt2)**2
!        bs = 2.*dot_product(vpoint%set(t1)%infr(1:3)-vpoint%set(t2)%infr(1:3),vpoint%set(t2)%infr(1:3)-pcpt)/(a*deltat)**2 + &
!           & 2.*(dt1-dt2)*(dble(i)-dt2)
!        cs = dot_product(vpoint%set(t2)%infr(1:3)-pcpt,vpoint%set(t2)%infr(1:3)-pcpt)/(a*deltat)**2 - (dble(i)-dt2)**2

!        ls = min((-1.*bs+sqrt(bs**2-4.*as*cs))/(2.*as),(-1.*bs-sqrt(bs**2-4.*as*cs))/(2.*as))

!        d1r2 = ls
!        d1r1 = 1.-ls

          if (associated(vnpnt2)) then
            d4r1 = vnpnt2%centroid
            rdt = sqrt(dot_product(d4r1(1:3)-r4in(1:3),d4r1(1:3)-r4in(1:3)))/(a*deltat)
            d1r2 = d4r1(4) - (r4in(4) - rdt)
            if (d1r1*d1r2 <= 0.) then

        if (d1r1 >= 0.) then

        dt1 = vnpnt1%centroid(4)
        dt2 = vnpnt2%centroid(4)
        as = dot_product(vnpnt1%centroid(1:3)-vnpnt2%centroid(1:3),vnpnt1%centroid(1:3)-vnpnt2%centroid(1:3))/ &
           & (a*deltat)**2 - (dt1-dt2)**2
        bs = 2.*dot_product(vnpnt1%centroid(1:3)-vnpnt2%centroid(1:3),vnpnt2%centroid(1:3)-r4in(1:3))/(a*deltat)**2 + &
           & 2.*(dt1-dt2)*(r4in(4)-dt2)
        cs = dot_product(vnpnt2%centroid(1:3)-r4in(1:3),vnpnt2%centroid(1:3)-r4in(1:3))/(a*deltat)**2 - (r4in(4)-dt2)**2

        ls = min((-1.*bs+sqrt(bs**2-4.*as*cs))/(2.*as),(-1.*bs-sqrt(bs**2-4.*as*cs))/(2.*as))

        fr2 = ls
        fr1 = 1.-ls

        else
        dt1 = vnpnt2%centroid(4)
        dt2 = vnpnt1%centroid(4)
        as = dot_product(vnpnt2%centroid(1:3)-vnpnt1%centroid(1:3),vnpnt2%centroid(1:3)-vnpnt1%centroid(1:3))/ &
           & (a*deltat)**2 - (dt1-dt2)**2
        bs = 2.*dot_product(vnpnt2%centroid(1:3)-vnpnt1%centroid(1:3),vnpnt1%centroid(1:3)-r4in(1:3))/(a*deltat)**2 + &
           & 2.*(dt1-dt2)*(r4in(4)-dt2)
        cs = dot_product(vnpnt1%centroid(1:3)-r4in(1:3),vnpnt1%centroid(1:3)-r4in(1:3))/(a*deltat)**2 - (r4in(4)-dt2)**2

        ls = min((-1.*bs+sqrt(bs**2-4.*as*cs))/(2.*as),(-1.*bs-sqrt(bs**2-4.*as*cs))/(2.*as))

        fr1 = ls
        fr2 = 1.-ls
        end if


!              fr1 = d1r1
!              fr2 = d1r2
              exit
            else
              vnpnt1 => vnpnt2
              d1r1 = d1r2
              n4i1 = n4i2
            end if
          else
            exit
          end if
        else
          exit
        end if
      else
        exit
      end if
    end do    
  
  end subroutine calcfind_interp

!###################################################################################################################################

  subroutine findbox(nlevel,deltat,TLB,lev,r4in,npnt,n4iout)
    use constants
    use fourdgridtypes
    implicit none

    integer,intent(in) :: nlevel,lev
    type(node) :: TLB
    real*8,intent(in) :: deltat
    real*8,dimension(4),intent(in) :: r4in
    integer,dimension(4,0:nlevel) :: branch
    type(node),pointer :: npnt
    integer,dimension(4),intent(out) :: n4iout
    integer :: j,d1i2

    npnt => null()
    branch(:,0) = box4d(deltat,r4in)
    do j=1,nlevel,1
      branch(:,j) = i4parent(branch(:,j-1))
    end do
    !find octree box in hextree
    d1i2 = TLBchild(branch(:,nlevel))
    if (associated(TLB%c(d1i2)%child)) then
      npnt => TLB%c(d1i2)%child
      do j=nlevel-1,lev,-1
        d1i2 = nchild(branch(:,j+1),branch(:,j))
        if (associated(npnt%c(d1i2)%child)) then
          npnt => npnt%c(d1i2)%child
        else
          npnt => null()
          exit
        end if
      end do
    else
      npnt => null()
    end if
    n4iout = branch(:,lev)

  end subroutine findbox

!  subroutine calcfind_zero(nlevel,TLB,lev,d4i1,deltat,r4in,branch,npnt,fr1,fr2)
!    use constants
!    use fourdgridtypes
!    implicit none
!    integer,intent(in) :: nlevel,lev
!    type(node) :: TLB
!    integer,dimension(4),intent(in) :: d4i1
!    real*8,intent(in) :: deltat
!    real*8,dimension(4),intent(in) :: r4in
!    integer,dimension(4,0:nlevel),intent(in) :: branch
!    type(node),pointer :: npnt
!    integer :: d1i1,d1i2,j
!    real*8,dimension(4) :: d4r1,d4r2,d4r3
!    integer,dimension(4) :: d4i2
!    integer,dimension(4,0:nlevel) :: vbranch1,vbranch2
!    real*8 :: rdt
!    real*8 :: d1r1,d1r2,fr1,fr2



!    



!    !calc centroid coords in 4D at right time
!    d4r1 = boxcentroid(deltat,lev,d4i1)
!    rdt = sqrt(dot_product(d4r1(1:3)-r4in(1:3),d4r1(1:3)-r4in(1:3)))/(a*deltat)
!!    rdt = abs(d4r1(1)-r4in(1))/(a*deltat)
!    d4r1(4) = r4in(4) - rdt

!    if (d4r1(4) .ge. 1.) then
!      !find tree for centroid at right time
!      vbranch1(:,0) = box4d(deltat,d4r1)
!      do j=1,nlevel,1
!        vbranch1(:,j) = i4parent(vbranch1(:,j-1))
!      end do
!      !find octree box in hextree
!      d1i2 = TLBchild(vbranch1(:,nlevel))
!      if (associated(TLB%c(d1i2)%child)) then
!        vnpnt1 => TLB%c(d1i2)%child
!        do j=nlevel-1,lev,-1
!          d1i2 = nchild(vbranch1(:,j+1),vbranch1(:,j))
!          if (associated(vnpnt1%c(d1i2)%child)) then
!            vnpnt1 => vnpnt1%c(d1i2)%child
!          else
!            vnpnt1 => null()
!            exit
!          end if
!        end do
!      else
!        vnpnt1 => null()
!      end if

!      if (associated(vnpnt1)) then
!        d4r2 = boxcentroid(deltat,lev,vbranch1(:,lev))
!        d1r1 = d4r2(4) - (r4in(4) - rdt)
!        d4i2(1:3) = vbranch1(1:3,lev)
!        if (d1r1 > 0.) then
!          if (vbranch1(4,lev) == 1) then
!            d4i2(4) = vbranch1(4,lev) - 2
!          else
!            d4i2(4) = vbranch1(4,lev) - 1
!          end if
!        else if (d1r1 < 0.) then
!          if (vbranch1(4,lev) == -1) then
!            d4i2(4) = vbranch1(4,lev) + 2
!          else
!            d4i2(4) = vbranch1(4,lev) + 1
!          end if
!        else
!           d4i2(4) = vbranch1(4,lev)
!        end if

!        if (d4i2(4) /= vbranch1(4,lev)) then
!          d4r3 = boxcentroid(deltat,lev,d4i2)
!          d1r2 = d4r3(4) - (r4in(4) - rdt)

!          vbranch2(:,0) = box4d(deltat,d4r3)
!          do j=1,nlevel,1
!            vbranch2(:,j) = i4parent(vbranch2(:,j-1))
!          end do
!          !check it's still inside tree
!          if (abs(vbranch2(4,nlevel)) == 1) then
!            !find octree box in hextree
!            d1i2 = TLBchild(vbranch2(:,nlevel))
!            if (associated(TLB%c(d1i2)%child)) then
!              vnpnt2 => TLB%c(d1i2)%child
!              do j=nlevel-1,lev,-1
!                d1i2 = nchild(vbranch2(:,j+1),vbranch2(:,j))
!                if (associated(vnpnt2%c(d1i2)%child)) then
!                  vnpnt2 => vnpnt2%c(d1i2)%child
!                else
!                  vnpnt2 => null()
!                  exit
!                end if
!              end do
!            end if
!          else
!            vnpnt2 => null()
!          end if
!       else
!          vnpnt2 => null()
!        end if
!      end if

!      if (associated(vnpnt1) .and. associated(vnpnt2)) then


!        fr1 = d1r1 ! abs(d1r1)/(abs(d1r1)+abs(d1r2))
!        fr2 = d1r2 !abs(d1r2)/(abs(d1r1)+abs(d1r2))

!      end if

!    end if

!  end subroutine calcfind_zero

!###################################################################################################################################








!  subroutine vbscanTLB(TLB,vbranch,nlevel,blev,vnpntchild,ltag,d1i1)
!    use basictypes
!    use fourdgridtypes
!    implicit none
!    type(node),intent(in) :: TLB
!    integer,dimension(4,0:nlevel) :: vbranch
!    integer,intent(in) :: nlevel,blev
!    type(node),pointer :: vnpntchild
!    type(node),pointer,save :: npnt
!    logical :: ltag
!    integer :: d1i1

!    d1i1 = TLBchild(vbranch(:,nlevel))
!    if (associated(TLB%c(d1i1)%child)) then
!      npnt => TLB%c(d1i1)%child
!      call vbscan(npnt,vbranch,nlevel,nlevel,blev,vnpntchild,ltag,d1i1)
!    else
!      vnpntchild => null()
!      ltag = .false.
!    end if
!  end subroutine vbscanTLB

!  recursive subroutine vbscan(npnt,vbranch,nlevel,lev,blev,vnpntchild,ltag,d1i1)
!    use basictypes
!    use fourdgridtypes
!    implicit none
!    type(node),pointer :: npnt
!    integer,dimension(4,0:nlevel) :: vbranch
!    integer,intent(in) :: nlevel,lev,blev
!    type(node),pointer :: vnpntchild
!    logical :: ltag
!    integer :: d1i1

!    d1i1 = nchild(vbranch(:,lev),vbranch(:,lev-1))
!    if (associated(npnt%c(d1i1)%child)) then
!      if (lev == blev) then
!        vnpntchild => npnt%c(d1i1)%child
!        ltag = .true.
!      else
!        npnt => npnt%c(d1i1)%child
!        call vbscan(npnt,vbranch,nlevel,lev,blev,vnpntchild,ltag,d1i1)
!      end if
!    else
!      vnpntchild => null()
!      ltag = .false.
!    end if
!  end subroutine vbscan

!###################################################################################################################################
!###################################################################################################################################
!!calc for point
!subroutine calcinfl(nlevel,a,tsteps,deltat,tree,r,t,v,vind,gradind)
!  use vtxparticle
!  use basictypes
!  use fourdgridtypes
!  use pg
!  implicit none

!  integer,intent(in) :: nlevel,tsteps
!  real*8,intent(in) :: a,deltat
!  type(level),dimension(0:nlevel) :: tree
!  real*8,dimension(3),intent(in) :: v
!  real*8,dimension(4),intent(in) :: r
!  integer,intent(in) :: t

!  real*8,dimension(3),intent(inout) :: vind
!  real*8,dimension(3,3),intent(inout) :: gradind

!  integer,dimension(4) :: ipos,lup
!  integer,dimension(4) :: z4i = 0.
!  type(rmset),dimension(:),pointer :: vtxset
!  type(rmset),pointer :: vtxpnt

!  integer :: d1i1,d1i2,istop,i,j
!  real*8 :: d1r1,d1r2
!  real*8,dimension(3) :: d3r1,d3r2,d3r
!  real*8,dimension(3,3) :: d33r
!  logical :: ffl

!  vind = 0.
!  gradind = 0.

!  ipos = box(deltat,r,t)

!  if (associated(tree(0)%elem(ipos(1),ipos(2),ipos(3),ipos(4))%npnt)) then
!    do i=1,tree(0)%elem(ipos(1),ipos(2),ipos(3),ipos(4))%npnt%nvtx,1
!      vtxset => tree(0)%elem(ipos(1),ipos(2),ipos(3),ipos(4))%npnt%bvtx(i)%vset
!      vtxpnt => tree(0)%elem(ipos(1),ipos(2),ipos(3),ipos(4))%npnt%bvtx(i)%vpnt
!      !calculate influence of particles in same box interpolated to time
!      istop = (tsteps-size(vtxset)) + 1
!      call interppnt3(r,v,tsteps,t,istop,a,deltat,vtxset,d1i1,d1i2,d1r1,d1r2,ffl)

!      if (ffl .eqv. .false.) then
!        d3r1 = intp(vtxset(d1i1)%infr(1:3),vtxset(d1i2)%infr(1:3),d1r1,d1r2)
!        d3r2 = intp(vtxset(d1i1)%minfr,vtxset(d1i2)%minfr,d1r1,d1r2)
!        call vtxinduced(scaler(a,v,r(1:3)),scaler(a,v,d3r1),d3r2,d3r)
!        vind = vind + d3r
!        call shearind(scaler(a,v,r(1:3)),scaler(a,v,d3r1),d3r2,d33r)
!        gradind = gradind + d33r
!      end if

!    end do

!    do i=1,4,1
!      lup(i) = up(ipos(i))
!    end do
!    call calcinflupbox(nlevel,tree,0,ipos,lup,r,v,vind,gradind)
!  else
!    do i=1,nlevel,1
!      ipos = lup
!      do j=1,4,1
!        lup(j) = up(ipos(j))
!      end do
!      if (associated(tree(i)%elem(ipos(1),ipos(2),ipos(3),ipos(4))%npnt)) then
!        call calcinflupbox(nlevel,tree,i,ipos,lup,r,v,vind,gradind)
!        exit
!      end if
!    end do
!  end if

!end subroutine calcinfl




!recursive subroutine calcinflupbox(nlevel,tree,lev,lin,lup,r,v,vind,gradind)
!  use vtxparticle
!  use fourdgridtypes
!  use pg
!  use constants
!  implicit none
!  integer,intent(in) :: nlevel
!  type(level),dimension(:) :: tree
!  type(node),pointer :: npnter
!  integer,dimension(4),intent(in) :: lin,lup
!  integer,intent(in) :: lev
!  real*8,dimension(3),intent(in) :: v
!  real*8,dimension(4),intent(in) :: r
!  real*8,dimension(3),intent(inout) :: vind
!  real*8,dimension(3,3),intent(inout) :: gradind

!  integer :: d1i,d1i1,d1i2,thigh,tlow,i,j
!  real*8 :: d1r1,d1r2
!  integer,dimension(4) :: d4i1,d4i2
!  integer,dimension(3) :: lrbox
!  integer,dimension(4) :: lupnxt
!  real*8,dimension(3) :: d3r,d3r1,d3r2
!  real*8,dimension(3,3) :: d33r
!  logical :: ffl



!  npnter => tree(lev+1)%elem(lup(1),lup(2),lup(3),lup(4))%npnt

!  do i=1,4,1
!    d4i1(i) = lin(i) - (down(lup(i))-1)
!  end do
!  d1i = 1 + d4i1(1) + 2*d4i1(2)+ 4*d4i1(3) !+ 8*d4i(4)
!  !only do first aid and interp others
!  do i=1,8,1
!    if (i /= d1i) then
!      if (associated(npnter%c(i)%child)) then
!        !calc vortex influence
!      end if
!      !calc box numbering below for interpolation
!      do j=1,4,1
!        d4i1(j) = down(lup(j))
!      end do
!      d4i2(3) = int(dble(i-1)/4.)
!      d4i2(2) = int(dble(i-1 - d4i1(3)*4)/2.)
!      d4i2(1) = int(dble(i-1 - d4i1(3)*4 - d4i1(2)*2))

!      do j=1,3,1
!        lrbox(j) = d4i2(j) + (down(lup(j))-1)
!      end do

!      thigh = min(down(lup(4)),size(tree(lev)%elem(lrbox(1),lrbox(2),lrbox(3),:)))
!      tlow = max(down(lup(4))-1,1)

!      call interpbox(tree,lev,thigh,tlow,r,lrbox(1:3),d1i1,d1i2,d1r1,d1r2,ffl)

!      if (ffl .eqv. .false.) then
!        d3r1 = intp(tree(lev)%elem(lrbox(1),lrbox(2),lrbox(3),d1i1)%npnt%centroid(1:3), &
!          & tree(lev)%elem(lrbox(1),lrbox(2),lrbox(3),d1i2)%npnt%centroid(1:3),d1r1,d1r2)
!        d3r2 = intp(tree(lev)%elem(lrbox(1),lrbox(2),lrbox(3),d1i1)%npnt%mentroid, &
!          & tree(lev)%elem(lrbox(1),lrbox(2),lrbox(3),d1i2)%npnt%mentroid,d1r1,d1r2)
!        call vtxinduced(scaler(a,v,r(1:3)),scaler(a,v,d3r1),d3r2,d3r)
!        vind = vind + d3r
!        call shearind(scaler(a,v,r(1:3)),scaler(a,v,d3r1),d3r2,d33r)
!        gradind = gradind + d33r
!      end if

!    end if
!  end do

!  if (lev<nlevel) then
!    do i=1,3,1
!      lupnxt(i) = up(lup(i))
!    end do
!    lupnxt(4) = up(d1i1)
!    call calcinflupbox(nlevel,tree,lev+1,lup,lupnxt,r,v,vind,gradind)
!  end if

!end subroutine calcinflupbox



subroutine interpbox(tree,lev,thigh,tlow,r4in,lin3i,d1i1,d1i2,d1r1,d1r2,ff)
  use fourdgridtypes
  implicit none
  type(level),dimension(:) :: tree
  integer,intent(in) :: lev,thigh,tlow
  real*8,dimension(4),intent(in) :: r4in
  integer,dimension(3),intent(in) :: lin3i
  integer :: d1i1,d1i2
  real*8 :: d1r1,d1r2
  logical :: ff
  integer :: i

  d1r1 = k2r(r4in,tree(lev)%elem(lin3i(1),lin3i(2),lin3i(3),thigh)%npnt%centroid)

  if (d1r1 > 0.) then
    do i=thigh-1,tlow,1
      d1r2 = k2r(r4in,tree(lev)%elem(lin3i(1),lin3i(2),lin3i(3),i)%npnt%centroid)
      if (d1r2 > 0.) then
        d1r1 = d1r2
      else
        d1i1 = i+1
        d1i2 = i
        ff = .false.
        exit
      end if

      if (i == tlow) then
        d1i1 = tlow+1
        d1i2 = tlow
        d1r1 = 0.
        d1r2 = -1.
        ff = .false.
      end if
    end do
  else
    d1i1 = thigh
    d1i2 = max(1,thigh-1)
    d1r1 = 0.
    d1r2 = -1.
    ff = .false.
!    ff = .true.
  end if

end subroutine interpbox



function k2r(r4in,r4vtx)
  use constants
  implicit none
  real*8 :: k2r
  real*8,dimension(4),intent(in) :: r4in,r4vtx
  k2r = sqrt(dot_product(r4vtx(1:3)-r4in(1:3),r4vtx(1:3)-r4in(1:3)))/a - dble(r4in(4)-r4vtx(4))
end function k2r





!###################################################################################################################################

end module fourdgridsf
