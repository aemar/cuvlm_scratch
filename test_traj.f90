program test_trajectory
  use functions
  use data_types
  use input
  implicit none

  integer,parameter :: n_part = 1
  integer,parameter :: t_steps = 400
  integer :: nt, i, j, iost
  character(len=4) :: istring
  character(len=6) :: dir1 = 'input/'
  character(len=7) :: dir2 = 'output/'

  type(part),dimension(n_part) :: parts
  real(kind=8),dimension(12) :: tvec

  real(kind=8),dimension(3) :: pcpt,pout
  real(kind=8) :: te

  deltat = 0.001

  !Read parts
  do i=1,n_part,1
    write(istring,'(i4.4)'), i
    open(unit=11,file=trim(dir1//'vlatin'//trim(istring)//'.dat'),status='old',action='read',iostat=iost)
    print*, 'Reading file', trim('vlatin'//trim(istring)//'.dat')
    read(11,*) parts(i)%np
    read(11,*) parts(i)%npan
    allocate(parts(i)%pts(3, parts(i)%np, t_steps),parts(i)%pans(parts(i)%npan))
    do j=1,parts(i)%np,1
      read(11,*) parts(i)%pts(:, j, 1)
    end do
    print*, 'Points read.'
    do j=1,parts(i)%npan,1
      read(11,*) parts(i)%pans(j)%typ, parts(i)%pans(j)%ver
    end do
    print*, 'Panels read.'
    do j=1,parts(i)%npan,1
      read(11,*) parts(i)%pans(j)%ngh
    end do
    print*, 'Neighbours read.'
    close(11)
  end do

  do nt=1,t_steps,1

    write(istring,'(i4.4)') nt
    open (unit=11,file=dir2//'part.csv.'//trim(istring),status='replace',action='write',iostat=iost)
    write (11,'(A33)') 'x coord, y coord, z coord, scalar'

    do i=1,n_part,1

      do j=1,parts(i)%np,1
        call traj(i, (nt-1.)*deltat, parts(i)%pts(:,j,1), tvec)
        write (11,'(F10.5,1A,1X,F10.5,1A,1X,F10.5,1A,1X,F10.5)') &
          & tvec(1),',',tvec(2),',',tvec(3),',',sqrt(dot_product(tvec(4:6),tvec(4:6)))
      end do

    end do

    close(11)

  end do

  pcpt = (/-2.7,-2.,0./)
  open (unit=12,file=dir2//'test_interp_out.csv',status='replace',action='write',iostat=iost)
  write (12,'(A33)') 'x coord, y coord, z coord, scalar'
  do i=1,n_part,1
    do j=1,parts(i)%np,1
      call traj_interp(pcpt, 200.*deltat, i, parts(i)%pts(:,j,1), te)
      pout = traj_pos(i, te, parts(i)%pts(:,j,1))
      write (12,'(F10.5,1A,1X,F10.5,1A,1X,F10.5,1A,1X,F10.5)') pout(1),',',pout(2),',',pout(3),',',te
    end do
  end do
  close(12)


  pcpt = (/-2.7,-2.,0./)
  do nt=1,t_steps,1

    write(istring,'(i4.4)') nt
    open (unit=11,file=dir2//'part_2.csv.'//trim(istring),status='replace',action='write',iostat=iost)
    write (11,'(A33)') 'x coord, y coord, z coord, scalar'

    do i=1,n_part,1

      do j=1,parts(i)%np,1
        call traj_interp(pcpt, nt*deltat, i, parts(i)%pts(:,j,1), te)
        pout = traj_pos(i, te, parts(i)%pts(:,j,1))
        write (11,'(F10.5,1A,1X,F10.5,1A,1X,F10.5,1A,1X,F10.5)') pout(1),',',pout(2),',',pout(3),',',te
      end do

    end do

    close(11)

  end do
  
  


end program test_trajectory
