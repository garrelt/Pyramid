module type
	
use precision, only: dp

implicit none
type tv
real(kind=dp), dimension(:,:), allocatable :: transform_volume
end type tv

type position
  real(kind=dp), dimension(:,:), allocatable :: split_position	
end type position

type layer
  type(position), dimension(:), allocatable :: upper_layer
end type 
	
type pyramid
  real(kind=dp) :: volume 
  real(kind=dp), dimension(:,:), allocatable :: cellsize
  real(kind=dp), dimension(:,:), allocatable :: solid_angle
  integer, dimension(:,:), allocatable :: case
  type(tv), dimension(:,:), allocatable :: transform
end type pyramid

type root
  real(kind=dp), dimension(:), allocatable :: with_parent_height
end type root
	
type cartesian
  type(tv), dimension(:,:), allocatable :: transform
end type cartesian

type data
  real(kind=dp), dimension(:,:), allocatable :: HI_density
  real(kind=dp), dimension(:,:), allocatable :: HeI_density
  real(kind=dp), dimension(:,:), allocatable :: HeII_density    
  real(kind=dp), dimension(:,:), allocatable :: HI_column_density_in
  real(kind=dp), dimension(:,:), allocatable :: HeI_column_density_in
  real(kind=dp), dimension(:,:), allocatable :: HeII_column_density_in 
  real(kind=dp), dimension(:,:), allocatable :: HI_column_density_out
  real(kind=dp), dimension(:,:), allocatable :: HeI_column_density_out
  real(kind=dp), dimension(:,:), allocatable :: HeII_column_density_out    
  real(kind=dp), dimension(:,:), allocatable :: HI_photoionization_rate  
  real(kind=dp), dimension(:,:), allocatable :: HeI_photoionization_rate  
  real(kind=dp), dimension(:,:), allocatable :: HeII_photoionization_rate      	
  real(kind=dp), dimension(:,:), allocatable :: photoheating_rate   	  	
end type data

! photrates contains all the photo-ionization rates and heating rates
!type photrates
!   real(kind=dp) :: photo_cell       ! HI photoionization rate of the cell
!   real(kind=dp) :: heat_cell        ! HI heating rate of the cell       
!   real(kind=dp) :: photo_in         ! HI photoionization rate incoming to the cell    
!   real(kind=dp) :: heat_in          ! HI heating rate incoming to the cell
!   real(kind=dp) :: photo_out        ! HI photoionization rate outgoing from the cell
!   real(kind=dp) :: heat_out         ! HI heating rate outgoing from the cell
!end type photrates
type photrates    
   real(kind=dp) :: photo_cell_HI          ! HI photoionization rate of the cell    
   real(kind=dp) :: photo_cell_HeI         ! HeI photoionization rate of the cell    
   real(kind=dp) :: photo_cell_HeII        ! HeII photoionization rate of the cell    
   real(kind=dp) :: heat_cell_HI           ! HI heating rate of the cell       
   real(kind=dp) :: heat_cell_HeI          ! HeI heating rate of the cell    
   real(kind=dp) :: heat_cell_HeII         ! HeII heating rate of the cell          
   real(kind=dp) :: photo_in_HI            ! HI photoionization rate incoming to the cell    
   real(kind=dp) :: photo_in_HeI           ! HeI photoionization rate incoming to the cell
   real(kind=dp) :: photo_in_HeII          ! HeII photoionization rate incoming to the cell
   real(kind=dp) :: heat_in_HI             ! HI heating rate incoming to the cell
   real(kind=dp) :: heat_in_HeI            ! HeI heating rate incoming to the cell
   real(kind=dp) :: heat_in_HeII           ! HeII heating rate incoming to the cell 
   real(kind=dp) :: photo_out_HI           ! HI photoionization rate outgoing from the cell
   real(kind=dp) :: photo_out_HeI          ! HeI photoionization rate outgoing from the cell
   real(kind=dp) :: photo_out_HeII         ! HeII photoionization rate outgoing from the cell 
   real(kind=dp) :: heat_out_HI            ! HI heating rate outgoing from the cell
   real(kind=dp) :: heat_out_HeI           ! HeI heating rate outgoing from the cell
   real(kind=dp) :: heat_out_HeII          ! HeII heating rate outgoing from the cell
   real(kind=dp) :: heat                   ! Total heating rate of the cell
   real(kind=dp) :: photo_in               ! Total photoionization rate incoming to the cell   
end type photrates


end module type
