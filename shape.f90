module shape
use precision, only: dp

implicit none

type pyramid
real(kind=dp), dimension(:,:), allocatable :: vol 
real(kind=dp), dimension(:,:), allocatable :: dr 
real(kind=dp), dimension(:,:), allocatable :: r 
real(kind=dp), dimension(:,:), allocatable :: ndens 
real(kind=dp), dimension(:,:), allocatable :: xHI
real(kind=dp), dimension(:,:), allocatable :: temperature
 
end type pyramid

end module shape
