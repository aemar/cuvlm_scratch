module data_types
  implicit none

  type pan
    character(len=1) :: typ
!    integer,dimension(:),allocatable :: ver,ngh
    integer,dimension(4) :: ver,ngh
    real*8,dimension(3) :: cpt
  end type pan

  type part
    integer :: np,npan
    real(kind=8),dimension(:,:,:),allocatable :: pts
    type(pan),dimension(:),allocatable :: pans
  end type part

end module data_types
