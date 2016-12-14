module array
	
  use precision, only: dp
  use input, only: level_py, level_sl, partition, layer
  use type, only: pyramid, cartesian, data, root
  
  implicit none
	
  type(pyramid), dimension(1:level_py) :: pyramid_grid
  type(cartesian), dimension(1:level_py) :: cartesian_grid
  type(root), dimension(:,:), allocatable :: parent_position
  type(data), dimension(1:level_py) :: px_py_pz_dx
  type(data), dimension(1:level_py) :: px_py_pz_dy	  
  type(data), dimension(1:level_py) :: px_py_pz_dz	  
  type(data), dimension(1:level_py) :: px_py_nz_dx
  type(data), dimension(1:level_py) :: px_py_nz_dy	  
  type(data), dimension(-level_py:-1) :: px_py_nz_dz	
  type(data), dimension(1:level_py) :: px_ny_pz_dx
  type(data), dimension(-level_py:-1) :: px_ny_pz_dy	  
  type(data), dimension(1:level_py) :: px_ny_pz_dz	  
  type(data), dimension(1:level_py) :: px_ny_nz_dx
  type(data), dimension(-level_py:-1) :: px_ny_nz_dy	  
  type(data), dimension(-level_py:-1) :: px_ny_nz_dz 	   
  type(data), dimension(-level_py:-1) :: nx_py_pz_dx
  type(data), dimension(1:level_py) :: nx_py_pz_dy	  
  type(data), dimension(1:level_py) :: nx_py_pz_dz	  
  type(data), dimension(-level_py:-1) :: nx_py_nz_dx
  type(data), dimension(1:level_py) :: nx_py_nz_dy	  
  type(data), dimension(-level_py:-1) :: nx_py_nz_dz	
  type(data), dimension(-level_py:-1) :: nx_ny_pz_dx
  type(data), dimension(-level_py:-1) :: nx_ny_pz_dy	  
  type(data), dimension(1:level_py) :: nx_ny_pz_dz	  
  type(data), dimension(-level_py:-1) :: nx_ny_nz_dx
  type(data), dimension(-level_py:-1) :: nx_ny_nz_dy	  
  type(data), dimension(-level_py:-1) :: nx_ny_nz_dz 

  real(kind=dp), dimension(1:2*level_py,1:2*level_py,1:2*level_py) :: pyramid_global_number_density_array
  real(kind=dp), dimension(1:2*level_py,1:2*level_py,1:2*level_py) :: pyramid_global_xHI_array
  real(kind=dp), dimension(1:2*level_py,1:2*level_py,1:2*level_py) :: pyramid_global_xHeI_array
  real(kind=dp), dimension(1:2*level_py,1:2*level_py,1:2*level_py) :: pyramid_global_xHeII_array
 
  ! for pyramid ray-tracing
  real(kind=dp), dimension(1:2*level_py,1:2*level_py,1:2*level_py) :: pyramid_global_HI_photoionization_rate_array
  real(kind=dp), dimension(1:2*level_py,1:2*level_py,1:2*level_py) :: pyramid_global_HeI_photoionization_rate_array  
  real(kind=dp), dimension(1:2*level_py,1:2*level_py,1:2*level_py) :: pyramid_global_HeII_photoionization_rate_array  
  real(kind=dp), dimension(1:2*level_py,1:2*level_py,1:2*level_py) :: pyramid_global_photoheating_rate_array
  
  ! for trilinear interpolated 
  real(kind=dp), dimension(1:2*level_py,1:2*level_py,1:2*level_py) :: interpolated_HI_photoionization_rate_array
  real(kind=dp), dimension(1:2*level_py,1:2*level_py,1:2*level_py) :: interpolated_HeI_photoionization_rate_array
  real(kind=dp), dimension(1:2*level_py,1:2*level_py,1:2*level_py) :: interpolated_HeII_photoionization_rate_array
  real(kind=dp), dimension(1:2*level_py,1:2*level_py,1:2*level_py) :: interpolated_photoheating_rate_array
       	
  real(kind=dp), dimension(1:level_py,1:level_py,1:level_py) :: pyramid_source_px_py_pz_HI_density_array
  real(kind=dp), dimension(1:level_py,1:level_py,-level_py:-1) :: pyramid_source_px_py_nz_HI_density_array  
  real(kind=dp), dimension(1:level_py,-level_py:-1,1:level_py) :: pyramid_source_px_ny_pz_HI_density_array
  real(kind=dp), dimension(1:level_py,-level_py:-1,-level_py:-1) :: pyramid_source_px_ny_nz_HI_density_array 
  real(kind=dp), dimension(-level_py:-1,1:level_py,1:level_py) :: pyramid_source_nx_py_pz_HI_density_array
  real(kind=dp), dimension(-level_py:-1,1:level_py,-level_py:-1) :: pyramid_source_nx_py_nz_HI_density_array  
  real(kind=dp), dimension(-level_py:-1,-level_py:-1,1:level_py) :: pyramid_source_nx_ny_pz_HI_density_array
  real(kind=dp), dimension(-level_py:-1,-level_py:-1,-level_py:-1) :: pyramid_source_nx_ny_nz_HI_density_array

  real(kind=dp), dimension(1:level_py,1:level_py,1:level_py) :: pyramid_source_px_py_pz_HeI_density_array
  real(kind=dp), dimension(1:level_py,1:level_py,-level_py:-1) :: pyramid_source_px_py_nz_HeI_density_array  
  real(kind=dp), dimension(1:level_py,-level_py:-1,1:level_py) :: pyramid_source_px_ny_pz_HeI_density_array
  real(kind=dp), dimension(1:level_py,-level_py:-1,-level_py:-1) :: pyramid_source_px_ny_nz_HeI_density_array 
  real(kind=dp), dimension(-level_py:-1,1:level_py,1:level_py) :: pyramid_source_nx_py_pz_HeI_density_array
  real(kind=dp), dimension(-level_py:-1,1:level_py,-level_py:-1) :: pyramid_source_nx_py_nz_HeI_density_array  
  real(kind=dp), dimension(-level_py:-1,-level_py:-1,1:level_py) :: pyramid_source_nx_ny_pz_HeI_density_array
  real(kind=dp), dimension(-level_py:-1,-level_py:-1,-level_py:-1) :: pyramid_source_nx_ny_nz_HeI_density_array
  
  real(kind=dp), dimension(1:level_py,1:level_py,1:level_py) :: pyramid_source_px_py_pz_HeII_density_array
  real(kind=dp), dimension(1:level_py,1:level_py,-level_py:-1) :: pyramid_source_px_py_nz_HeII_density_array  
  real(kind=dp), dimension(1:level_py,-level_py:-1,1:level_py) :: pyramid_source_px_ny_pz_HeII_density_array
  real(kind=dp), dimension(1:level_py,-level_py:-1,-level_py:-1) :: pyramid_source_px_ny_nz_HeII_density_array 
  real(kind=dp), dimension(-level_py:-1,1:level_py,1:level_py) :: pyramid_source_nx_py_pz_HeII_density_array
  real(kind=dp), dimension(-level_py:-1,1:level_py,-level_py:-1) :: pyramid_source_nx_py_nz_HeII_density_array  
  real(kind=dp), dimension(-level_py:-1,-level_py:-1,1:level_py) :: pyramid_source_nx_ny_pz_HeII_density_array
  real(kind=dp), dimension(-level_py:-1,-level_py:-1,-level_py:-1) :: pyramid_source_nx_ny_nz_HeII_density_array
      
  real(kind=dp), dimension(1:level_py,1:level_py,1:level_py) :: pyramid_source_px_py_pz_HI_photoionization_rate_array
  real(kind=dp), dimension(1:level_py,1:level_py,-level_py:-1) :: pyramid_source_px_py_nz_HI_photoionization_rate_array  
  real(kind=dp), dimension(1:level_py,-level_py:-1,1:level_py) :: pyramid_source_px_ny_pz_HI_photoionization_rate_array
  real(kind=dp), dimension(1:level_py,-level_py:-1,-level_py:-1) :: pyramid_source_px_ny_nz_HI_photoionization_rate_array 
  real(kind=dp), dimension(-level_py:-1,1:level_py,1:level_py) :: pyramid_source_nx_py_pz_HI_photoionization_rate_array
  real(kind=dp), dimension(-level_py:-1,1:level_py,-level_py:-1) :: pyramid_source_nx_py_nz_HI_photoionization_rate_array  
  real(kind=dp), dimension(-level_py:-1,-level_py:-1,1:level_py) :: pyramid_source_nx_ny_pz_HI_photoionization_rate_array
  real(kind=dp), dimension(-level_py:-1,-level_py:-1,-level_py:-1) :: pyramid_source_nx_ny_nz_HI_photoionization_rate_array

  real(kind=dp), dimension(1:level_py,1:level_py,1:level_py) :: pyramid_source_px_py_pz_HeI_photoionization_rate_array
  real(kind=dp), dimension(1:level_py,1:level_py,-level_py:-1) :: pyramid_source_px_py_nz_HeI_photoionization_rate_array  
  real(kind=dp), dimension(1:level_py,-level_py:-1,1:level_py) :: pyramid_source_px_ny_pz_HeI_photoionization_rate_array
  real(kind=dp), dimension(1:level_py,-level_py:-1,-level_py:-1) :: pyramid_source_px_ny_nz_HeI_photoionization_rate_array 
  real(kind=dp), dimension(-level_py:-1,1:level_py,1:level_py) :: pyramid_source_nx_py_pz_HeI_photoionization_rate_array
  real(kind=dp), dimension(-level_py:-1,1:level_py,-level_py:-1) :: pyramid_source_nx_py_nz_HeI_photoionization_rate_array  
  real(kind=dp), dimension(-level_py:-1,-level_py:-1,1:level_py) :: pyramid_source_nx_ny_pz_HeI_photoionization_rate_array
  real(kind=dp), dimension(-level_py:-1,-level_py:-1,-level_py:-1) :: pyramid_source_nx_ny_nz_HeI_photoionization_rate_array
  
  real(kind=dp), dimension(1:level_py,1:level_py,1:level_py) :: pyramid_source_px_py_pz_HeII_photoionization_rate_array
  real(kind=dp), dimension(1:level_py,1:level_py,-level_py:-1) :: pyramid_source_px_py_nz_HeII_photoionization_rate_array  
  real(kind=dp), dimension(1:level_py,-level_py:-1,1:level_py) :: pyramid_source_px_ny_pz_HeII_photoionization_rate_array
  real(kind=dp), dimension(1:level_py,-level_py:-1,-level_py:-1) :: pyramid_source_px_ny_nz_HeII_photoionization_rate_array 
  real(kind=dp), dimension(-level_py:-1,1:level_py,1:level_py) :: pyramid_source_nx_py_pz_HeII_photoionization_rate_array
  real(kind=dp), dimension(-level_py:-1,1:level_py,-level_py:-1) :: pyramid_source_nx_py_nz_HeII_photoionization_rate_array  
  real(kind=dp), dimension(-level_py:-1,-level_py:-1,1:level_py) :: pyramid_source_nx_ny_pz_HeII_photoionization_rate_array
  real(kind=dp), dimension(-level_py:-1,-level_py:-1,-level_py:-1) :: pyramid_source_nx_ny_nz_HeII_photoionization_rate_array
    
  real(kind=dp), dimension(1:level_py,1:level_py,1:level_py) :: pyramid_source_px_py_pz_photoheating_rate_array
  real(kind=dp), dimension(1:level_py,1:level_py,-level_py:-1) :: pyramid_source_px_py_nz_photoheating_rate_array  
  real(kind=dp), dimension(1:level_py,-level_py:-1,1:level_py) :: pyramid_source_px_ny_pz_photoheating_rate_array
  real(kind=dp), dimension(1:level_py,-level_py:-1,-level_py:-1) :: pyramid_source_px_ny_nz_photoheating_rate_array 
  real(kind=dp), dimension(-level_py:-1,1:level_py,1:level_py) :: pyramid_source_nx_py_pz_photoheating_rate_array
  real(kind=dp), dimension(-level_py:-1,1:level_py,-level_py:-1) :: pyramid_source_nx_py_nz_photoheating_rate_array  
  real(kind=dp), dimension(-level_py:-1,-level_py:-1,1:level_py) :: pyramid_source_nx_ny_pz_photoheating_rate_array
  real(kind=dp), dimension(-level_py:-1,-level_py:-1,-level_py:-1) :: pyramid_source_nx_ny_nz_photoheating_rate_array
       
  real(kind=dp), dimension(1:level_sl,1:level_sl,1:level_sl) :: short_global_number_density_array
  real(kind=dp), dimension(1:level_sl,1:level_sl,1:level_sl) :: short_global_xHI_array 
  real(kind=dp), dimension(1:level_sl,1:level_sl,1:level_sl) :: short_global_xHeI_array   
  real(kind=dp), dimension(1:level_sl,1:level_sl,1:level_sl) :: short_global_xHeII_array      
  real(kind=dp), dimension(1:level_sl,1:level_sl,1:level_sl) :: short_global_HI_photoionization_rate_array
  real(kind=dp), dimension(1:level_sl,1:level_sl,1:level_sl) :: short_global_HeI_photoionization_rate_array  
  real(kind=dp), dimension(1:level_sl,1:level_sl,1:level_sl) :: short_global_HeII_photoionization_rate_array  
  real(kind=dp), dimension(1:level_sl,1:level_sl,1:level_sl) :: short_global_photoheating_rate_array 
       
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: short_source_HI_density_array
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: short_source_HeI_density_array  
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: short_source_HeII_density_array  
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: short_source_HI_column_density_in
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: short_source_HeI_column_density_in  
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: short_source_HeII_column_density_in  
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: short_source_HI_column_density_out
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: short_source_HeI_column_density_out  
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: short_source_HeII_column_density_out  
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: short_source_HI_photoionization_rate_array
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: short_source_HeI_photoionization_rate_array  
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: short_source_HeII_photoionization_rate_array  
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: short_source_photoheating_rate_array 
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: short_source_shell_volume
  
  real(kind=dp), dimension(1:level_sl,1:level_sl,1:level_sl) :: long_global_number_density_array
  real(kind=dp), dimension(1:level_sl,1:level_sl,1:level_sl) :: long_global_xHI_array 
  real(kind=dp), dimension(1:level_sl,1:level_sl,1:level_sl) :: long_global_xHeI_array  
  real(kind=dp), dimension(1:level_sl,1:level_sl,1:level_sl) :: long_global_xHeII_array  
  real(kind=dp), dimension(1:level_sl,1:level_sl,1:level_sl) :: long_global_HI_photoionization_rate_array
  real(kind=dp), dimension(1:level_sl,1:level_sl,1:level_sl) :: long_global_HeI_photoionization_rate_array  
  real(kind=dp), dimension(1:level_sl,1:level_sl,1:level_sl) :: long_global_HeII_photoionization_rate_array  
  real(kind=dp), dimension(1:level_sl,1:level_sl,1:level_sl) :: long_global_photoheating_rate_array 
    
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: long_source_HI_density_array
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: long_source_HeI_density_array  
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: long_source_HeII_density_array  
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: long_source_HI_column_density_in 
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: long_source_HeI_column_density_in  
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: long_source_HeII_column_density_in  
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: long_source_HI_column_density_out
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: long_source_HeI_column_density_out  
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: long_source_HeII_column_density_out
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: long_source_HI_photoionization_rate_array
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: long_source_HeI_photoionization_rate_array  
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: long_source_HeII_photoionization_rate_array  
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: long_source_photoheating_rate_array
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: long_source_shell_volume

  ! for analytical solution of pyramid model
  real(kind=dp) :: pyramid_analytical_number_density
  real(kind=dp) :: pyramid_analytical_xHI
  real(kind=dp) :: pyramid_analytical_xHeI  
  real(kind=dp) :: pyramid_analytical_xHeII  
    
  real(kind=dp), dimension(1:2*level_py,1:2*level_py,1:2*level_py) :: pyramid_analytical_HI_column_density_in 
  real(kind=dp), dimension(1:2*level_py,1:2*level_py,1:2*level_py) :: pyramid_analytical_HeI_column_density_in  
  real(kind=dp), dimension(1:2*level_py,1:2*level_py,1:2*level_py) :: pyramid_analytical_HeII_column_density_in  
  real(kind=dp), dimension(1:2*level_py,1:2*level_py,1:2*level_py) :: pyramid_analytical_HI_column_density_out
  real(kind=dp), dimension(1:2*level_py,1:2*level_py,1:2*level_py) :: pyramid_analytical_HeI_column_density_out  
  real(kind=dp), dimension(1:2*level_py,1:2*level_py,1:2*level_py) :: pyramid_analytical_HeII_column_density_out
  real(kind=dp), dimension(1:2*level_py,1:2*level_py,1:2*level_py) :: pyramid_analytical_HI_photoionization_rate_array
  real(kind=dp), dimension(1:2*level_py,1:2*level_py,1:2*level_py) :: pyramid_analytical_HeI_photoionization_rate_array  
  real(kind=dp), dimension(1:2*level_py,1:2*level_py,1:2*level_py) :: pyramid_analytical_HeII_photoionization_rate_array  
  real(kind=dp), dimension(1:2*level_py,1:2*level_py,1:2*level_py) :: pyramid_analytical_photoheating_rate_array
  real(kind=dp), dimension(1:2*level_py,1:2*level_py,1:2*level_py) :: pyramid_analytical_shell_volume
         
  ! for analytical solution of cartesian model
  real(kind=dp) :: cartesian_analytical_number_density
  real(kind=dp) :: cartesian_analytical_xHI
  real(kind=dp) :: cartesian_analytical_xHeI  
  real(kind=dp) :: cartesian_analytical_xHeII  
    
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: cartesian_analytical_HI_column_density_in 
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: cartesian_analytical_HeI_column_density_in  
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: cartesian_analytical_HeII_column_density_in  
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: cartesian_analytical_HI_column_density_out
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: cartesian_analytical_HeI_column_density_out  
  real(kind=dp), dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: cartesian_analytical_HeII_column_density_out
  real(kind=dp),dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: &
  cartesian_analytical_HI_photoionization_rate_array
  real(kind=dp),dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: &
  cartesian_analytical_HeI_photoionization_rate_array  
  real(kind=dp),dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: &
  cartesian_analytical_HeII_photoionization_rate_array  
  real(kind=dp),dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: &
  cartesian_analytical_photoheating_rate_array
  real(kind=dp),dimension(-level_py:level_py,-level_py:level_py,-level_py:level_py) :: cartesian_analytical_shell_volume
       
  ! Integrands ( frequency, optical depth )
  real(kind=dp),dimension(:,:), allocatable :: bb_photo_thick_integrand
  real(kind=dp),dimension(:,:), allocatable :: bb_photo_thin_integrand
  real(kind=dp),dimension(:,:), allocatable :: pl_photo_thick_integrand
  real(kind=dp),dimension(:,:), allocatable :: pl_photo_thin_integrand
  real(kind=dp),dimension(:,:), allocatable :: bb_heat_thick_integrand_HI
  real(kind=dp),dimension(:,:), allocatable :: bb_heat_thick_integrand_HeI
  real(kind=dp),dimension(:,:), allocatable :: bb_heat_thick_integrand_HeII
  real(kind=dp),dimension(:,:), allocatable :: bb_heat_thin_integrand_HI
  real(kind=dp),dimension(:,:), allocatable :: bb_heat_thin_integrand_HeI
  real(kind=dp),dimension(:,:), allocatable :: bb_heat_thin_integrand_HeII
  real(kind=dp),dimension(:,:), allocatable :: pl_heat_thick_integrand_HI
  real(kind=dp),dimension(:,:), allocatable :: pl_heat_thick_integrand_HeI
  real(kind=dp),dimension(:,:), allocatable :: pl_heat_thick_integrand_HeII
  real(kind=dp),dimension(:,:), allocatable :: pl_heat_thin_integrand_HI
  real(kind=dp),dimension(:,:), allocatable :: pl_heat_thin_integrand_HeI
  real(kind=dp),dimension(:,:), allocatable :: pl_heat_thin_integrand_HeII

  ! Integration table ( optical depth, sub-bin )
  real(kind=dp),dimension(:,:), target, allocatable :: bb_photo_thick_table  
  real(kind=dp),dimension(:,:), target, allocatable :: bb_photo_thin_table
  real(kind=dp),dimension(:,:), target, allocatable :: pl_photo_thick_table 
  real(kind=dp),dimension(:,:), target, allocatable :: pl_photo_thin_table
  real(kind=dp),dimension(:,:), target, allocatable :: bb_heat_thick_table  
  real(kind=dp),dimension(:,:), target, allocatable :: bb_heat_thin_table
  real(kind=dp),dimension(:,:), target, allocatable :: pl_heat_thick_table  
  real(kind=dp),dimension(:,:), target, allocatable :: pl_heat_thin_table   
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pyramid_creation()

    ! allocation arrays of volume, length...
    integer :: count_x, count_y, count_l
    integer :: count
    integer :: delta_l
    
    do count = 1,level_py
	
      allocate(pyramid_grid(count)%cellsize(1:partition(count),1:partition(count)))
      allocate(pyramid_grid(count)%solid_angle(1:partition(count),1:partition(count)))
      allocate(pyramid_grid(count)%case(1:partition(count),1:partition(count)))
      allocate(pyramid_grid(count)%transform(1:partition(count),1:partition(count)))
	  
    end do
  
    allocate(parent_position(1:level_py,1:partition(level_py)))
  
    do count = 1,level_py
      do count_x = 1,partition(count)
        allocate(parent_position(count,count_x)%with_parent_height(1:count))
      enddo			
    enddo
	
 end subroutine pyramid_creation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cartesian_creation()

    ! allocation arrays of volume, length...
    integer :: count

    do count = 1,level_py
	  
      allocate(cartesian_grid(count)%transform(1:count,1:count))
	  
    end do

  end subroutine cartesian_creation
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  	  
  subroutine domain_creation()

    implicit none

    integer :: count

    do count = 1,level_py
		
      allocate(px_py_pz_dx(count)%HI_density(1:partition(count),1:partition(count)))
      allocate(px_py_pz_dy(count)%HI_density(1:partition(count),1:partition(count)))	  
      allocate(px_py_pz_dz(count)%HI_density(1:partition(count),1:partition(count)))	  
      allocate(px_py_nz_dx(count)%HI_density(1:partition(count),-partition(count):-1))
      allocate(px_py_nz_dy(count)%HI_density(-partition(count):-1,1:partition(count)))	  
      allocate(px_py_nz_dz(-count)%HI_density(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dx(count)%HI_density(-partition(count):-1,1:partition(count)))
      allocate(px_ny_pz_dy(-count)%HI_density(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dz(count)%HI_density(1:partition(count),-partition(count):-1))	  
      allocate(px_ny_nz_dx(count)%HI_density(-partition(count):-1,-partition(count):-1))
      allocate(px_ny_nz_dy(-count)%HI_density(-partition(count):-1,1:partition(count)))	  
      allocate(px_ny_nz_dz(-count)%HI_density(1:partition(count),-partition(count):-1))	
      allocate(nx_py_pz_dx(-count)%HI_density(1:partition(count),1:partition(count)))
      allocate(nx_py_pz_dy(count)%HI_density(1:partition(count),-partition(count):-1))	  
      allocate(nx_py_pz_dz(count)%HI_density(-partition(count):-1,1:partition(count)))  
      allocate(nx_py_nz_dx(-count)%HI_density(1:partition(count),-partition(count):-1))
      allocate(nx_py_nz_dy(count)%HI_density(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_py_nz_dz(-count)%HI_density(-partition(count):-1,1:partition(count)))	  
      allocate(nx_ny_pz_dx(-count)%HI_density(-partition(count):-1,1:partition(count)))
      allocate(nx_ny_pz_dy(-count)%HI_density(1:partition(count),-partition(count):-1))	  
      allocate(nx_ny_pz_dz(count)%HI_density(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dx(-count)%HI_density(-partition(count):-1,-partition(count):-1))
      allocate(nx_ny_nz_dy(-count)%HI_density(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dz(-count)%HI_density(-partition(count):-1,-partition(count):-1))	

      allocate(px_py_pz_dx(count)%HI_column_density_in(1:partition(count),1:partition(count)))
      allocate(px_py_pz_dy(count)%HI_column_density_in(1:partition(count),1:partition(count)))	  
      allocate(px_py_pz_dz(count)%HI_column_density_in(1:partition(count),1:partition(count)))	  
      allocate(px_py_nz_dx(count)%HI_column_density_in(1:partition(count),-partition(count):-1))
      allocate(px_py_nz_dy(count)%HI_column_density_in(-partition(count):-1,1:partition(count)))	  
      allocate(px_py_nz_dz(-count)%HI_column_density_in(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dx(count)%HI_column_density_in(-partition(count):-1,1:partition(count)))
      allocate(px_ny_pz_dy(-count)%HI_column_density_in(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dz(count)%HI_column_density_in(1:partition(count),-partition(count):-1))	  
      allocate(px_ny_nz_dx(count)%HI_column_density_in(-partition(count):-1,-partition(count):-1))
      allocate(px_ny_nz_dy(-count)%HI_column_density_in(-partition(count):-1,1:partition(count)))	  
      allocate(px_ny_nz_dz(-count)%HI_column_density_in(1:partition(count),-partition(count):-1))	
      allocate(nx_py_pz_dx(-count)%HI_column_density_in(1:partition(count),1:partition(count)))
      allocate(nx_py_pz_dy(count)%HI_column_density_in(1:partition(count),-partition(count):-1))	  
      allocate(nx_py_pz_dz(count)%HI_column_density_in(-partition(count):-1,1:partition(count)))  
      allocate(nx_py_nz_dx(-count)%HI_column_density_in(1:partition(count),-partition(count):-1))
      allocate(nx_py_nz_dy(count)%HI_column_density_in(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_py_nz_dz(-count)%HI_column_density_in(-partition(count):-1,1:partition(count)))	  
      allocate(nx_ny_pz_dx(-count)%HI_column_density_in(-partition(count):-1,1:partition(count)))
      allocate(nx_ny_pz_dy(-count)%HI_column_density_in(1:partition(count),-partition(count):-1))	  
      allocate(nx_ny_pz_dz(count)%HI_column_density_in(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dx(-count)%HI_column_density_in(-partition(count):-1,-partition(count):-1))
      allocate(nx_ny_nz_dy(-count)%HI_column_density_in(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dz(-count)%HI_column_density_in(-partition(count):-1,-partition(count):-1))
  
      allocate(px_py_pz_dx(count)%HI_column_density_out(1:partition(count),1:partition(count)))
      allocate(px_py_pz_dy(count)%HI_column_density_out(1:partition(count),1:partition(count)))	  
      allocate(px_py_pz_dz(count)%HI_column_density_out(1:partition(count),1:partition(count)))	  
      allocate(px_py_nz_dx(count)%HI_column_density_out(1:partition(count),-partition(count):-1))
      allocate(px_py_nz_dy(count)%HI_column_density_out(-partition(count):-1,1:partition(count)))	  
      allocate(px_py_nz_dz(-count)%HI_column_density_out(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dx(count)%HI_column_density_out(-partition(count):-1,1:partition(count)))
      allocate(px_ny_pz_dy(-count)%HI_column_density_out(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dz(count)%HI_column_density_out(1:partition(count),-partition(count):-1))	  
      allocate(px_ny_nz_dx(count)%HI_column_density_out(-partition(count):-1,-partition(count):-1))
      allocate(px_ny_nz_dy(-count)%HI_column_density_out(-partition(count):-1,1:partition(count)))	  
      allocate(px_ny_nz_dz(-count)%HI_column_density_out(1:partition(count),-partition(count):-1))	
      allocate(nx_py_pz_dx(-count)%HI_column_density_out(1:partition(count),1:partition(count)))
      allocate(nx_py_pz_dy(count)%HI_column_density_out(1:partition(count),-partition(count):-1))	  
      allocate(nx_py_pz_dz(count)%HI_column_density_out(-partition(count):-1,1:partition(count)))  
      allocate(nx_py_nz_dx(-count)%HI_column_density_out(1:partition(count),-partition(count):-1))
      allocate(nx_py_nz_dy(count)%HI_column_density_out(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_py_nz_dz(-count)%HI_column_density_out(-partition(count):-1,1:partition(count)))	  
      allocate(nx_ny_pz_dx(-count)%HI_column_density_out(-partition(count):-1,1:partition(count)))
      allocate(nx_ny_pz_dy(-count)%HI_column_density_out(1:partition(count),-partition(count):-1))	  
      allocate(nx_ny_pz_dz(count)%HI_column_density_out(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dx(-count)%HI_column_density_out(-partition(count):-1,-partition(count):-1))
      allocate(nx_ny_nz_dy(-count)%HI_column_density_out(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dz(-count)%HI_column_density_out(-partition(count):-1,-partition(count):-1))
  	  	  	  
      allocate(px_py_pz_dx(count)%HI_photoionization_rate(1:partition(count),1:partition(count)))
      allocate(px_py_pz_dy(count)%HI_photoionization_rate(1:partition(count),1:partition(count)))	  
      allocate(px_py_pz_dz(count)%HI_photoionization_rate(1:partition(count),1:partition(count)))	  
      allocate(px_py_nz_dx(count)%HI_photoionization_rate(1:partition(count),-partition(count):-1))
      allocate(px_py_nz_dy(count)%HI_photoionization_rate(-partition(count):-1,1:partition(count)))	  
      allocate(px_py_nz_dz(-count)%HI_photoionization_rate(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dx(count)%HI_photoionization_rate(-partition(count):-1,1:partition(count)))
      allocate(px_ny_pz_dy(-count)%HI_photoionization_rate(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dz(count)%HI_photoionization_rate(1:partition(count),-partition(count):-1))	  
      allocate(px_ny_nz_dx(count)%HI_photoionization_rate(-partition(count):-1,-partition(count):-1))
      allocate(px_ny_nz_dy(-count)%HI_photoionization_rate(-partition(count):-1,1:partition(count)))	  
      allocate(px_ny_nz_dz(-count)%HI_photoionization_rate(1:partition(count),-partition(count):-1))	
      allocate(nx_py_pz_dx(-count)%HI_photoionization_rate(1:partition(count),1:partition(count)))
      allocate(nx_py_pz_dy(count)%HI_photoionization_rate(1:partition(count),-partition(count):-1))	  
      allocate(nx_py_pz_dz(count)%HI_photoionization_rate(-partition(count):-1,1:partition(count)))  
      allocate(nx_py_nz_dx(-count)%HI_photoionization_rate(1:partition(count),-partition(count):-1))
      allocate(nx_py_nz_dy(count)%HI_photoionization_rate(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_py_nz_dz(-count)%HI_photoionization_rate(-partition(count):-1,1:partition(count)))	  
      allocate(nx_ny_pz_dx(-count)%HI_photoionization_rate(-partition(count):-1,1:partition(count)))
      allocate(nx_ny_pz_dy(-count)%HI_photoionization_rate(1:partition(count),-partition(count):-1))	  
      allocate(nx_ny_pz_dz(count)%HI_photoionization_rate(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dx(-count)%HI_photoionization_rate(-partition(count):-1,-partition(count):-1))
      allocate(nx_ny_nz_dy(-count)%HI_photoionization_rate(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dz(-count)%HI_photoionization_rate(-partition(count):-1,-partition(count):-1))				  

      allocate(px_py_pz_dx(count)%HeI_density(1:partition(count),1:partition(count)))
      allocate(px_py_pz_dy(count)%HeI_density(1:partition(count),1:partition(count)))	  
      allocate(px_py_pz_dz(count)%HeI_density(1:partition(count),1:partition(count)))	  
      allocate(px_py_nz_dx(count)%HeI_density(1:partition(count),-partition(count):-1))
      allocate(px_py_nz_dy(count)%HeI_density(-partition(count):-1,1:partition(count)))	  
      allocate(px_py_nz_dz(-count)%HeI_density(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dx(count)%HeI_density(-partition(count):-1,1:partition(count)))
      allocate(px_ny_pz_dy(-count)%HeI_density(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dz(count)%HeI_density(1:partition(count),-partition(count):-1))	  
      allocate(px_ny_nz_dx(count)%HeI_density(-partition(count):-1,-partition(count):-1))
      allocate(px_ny_nz_dy(-count)%HeI_density(-partition(count):-1,1:partition(count)))	  
      allocate(px_ny_nz_dz(-count)%HeI_density(1:partition(count),-partition(count):-1))	
      allocate(nx_py_pz_dx(-count)%HeI_density(1:partition(count),1:partition(count)))
      allocate(nx_py_pz_dy(count)%HeI_density(1:partition(count),-partition(count):-1))	  
      allocate(nx_py_pz_dz(count)%HeI_density(-partition(count):-1,1:partition(count)))  
      allocate(nx_py_nz_dx(-count)%HeI_density(1:partition(count),-partition(count):-1))
      allocate(nx_py_nz_dy(count)%HeI_density(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_py_nz_dz(-count)%HeI_density(-partition(count):-1,1:partition(count)))	  
      allocate(nx_ny_pz_dx(-count)%HeI_density(-partition(count):-1,1:partition(count)))
      allocate(nx_ny_pz_dy(-count)%HeI_density(1:partition(count),-partition(count):-1))	  
      allocate(nx_ny_pz_dz(count)%HeI_density(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dx(-count)%HeI_density(-partition(count):-1,-partition(count):-1))
      allocate(nx_ny_nz_dy(-count)%HeI_density(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dz(-count)%HeI_density(-partition(count):-1,-partition(count):-1))	

      allocate(px_py_pz_dx(count)%HeI_column_density_in(1:partition(count),1:partition(count)))
      allocate(px_py_pz_dy(count)%HeI_column_density_in(1:partition(count),1:partition(count)))	  
      allocate(px_py_pz_dz(count)%HeI_column_density_in(1:partition(count),1:partition(count)))	  
      allocate(px_py_nz_dx(count)%HeI_column_density_in(1:partition(count),-partition(count):-1))
      allocate(px_py_nz_dy(count)%HeI_column_density_in(-partition(count):-1,1:partition(count)))	  
      allocate(px_py_nz_dz(-count)%HeI_column_density_in(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dx(count)%HeI_column_density_in(-partition(count):-1,1:partition(count)))
      allocate(px_ny_pz_dy(-count)%HeI_column_density_in(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dz(count)%HeI_column_density_in(1:partition(count),-partition(count):-1))	  
      allocate(px_ny_nz_dx(count)%HeI_column_density_in(-partition(count):-1,-partition(count):-1))
      allocate(px_ny_nz_dy(-count)%HeI_column_density_in(-partition(count):-1,1:partition(count)))	  
      allocate(px_ny_nz_dz(-count)%HeI_column_density_in(1:partition(count),-partition(count):-1))	
      allocate(nx_py_pz_dx(-count)%HeI_column_density_in(1:partition(count),1:partition(count)))
      allocate(nx_py_pz_dy(count)%HeI_column_density_in(1:partition(count),-partition(count):-1))	  
      allocate(nx_py_pz_dz(count)%HeI_column_density_in(-partition(count):-1,1:partition(count)))  
      allocate(nx_py_nz_dx(-count)%HeI_column_density_in(1:partition(count),-partition(count):-1))
      allocate(nx_py_nz_dy(count)%HeI_column_density_in(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_py_nz_dz(-count)%HeI_column_density_in(-partition(count):-1,1:partition(count)))	  
      allocate(nx_ny_pz_dx(-count)%HeI_column_density_in(-partition(count):-1,1:partition(count)))
      allocate(nx_ny_pz_dy(-count)%HeI_column_density_in(1:partition(count),-partition(count):-1))	  
      allocate(nx_ny_pz_dz(count)%HeI_column_density_in(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dx(-count)%HeI_column_density_in(-partition(count):-1,-partition(count):-1))
      allocate(nx_ny_nz_dy(-count)%HeI_column_density_in(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dz(-count)%HeI_column_density_in(-partition(count):-1,-partition(count):-1))
	  
      allocate(px_py_pz_dx(count)%HeI_column_density_out(1:partition(count),1:partition(count)))
      allocate(px_py_pz_dy(count)%HeI_column_density_out(1:partition(count),1:partition(count)))	  
      allocate(px_py_pz_dz(count)%HeI_column_density_out(1:partition(count),1:partition(count)))	  
      allocate(px_py_nz_dx(count)%HeI_column_density_out(1:partition(count),-partition(count):-1))
      allocate(px_py_nz_dy(count)%HeI_column_density_out(-partition(count):-1,1:partition(count)))	  
      allocate(px_py_nz_dz(-count)%HeI_column_density_out(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dx(count)%HeI_column_density_out(-partition(count):-1,1:partition(count)))
      allocate(px_ny_pz_dy(-count)%HeI_column_density_out(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dz(count)%HeI_column_density_out(1:partition(count),-partition(count):-1))	  
      allocate(px_ny_nz_dx(count)%HeI_column_density_out(-partition(count):-1,-partition(count):-1))
      allocate(px_ny_nz_dy(-count)%HeI_column_density_out(-partition(count):-1,1:partition(count)))	  
      allocate(px_ny_nz_dz(-count)%HeI_column_density_out(1:partition(count),-partition(count):-1))	
      allocate(nx_py_pz_dx(-count)%HeI_column_density_out(1:partition(count),1:partition(count)))
      allocate(nx_py_pz_dy(count)%HeI_column_density_out(1:partition(count),-partition(count):-1))	  
      allocate(nx_py_pz_dz(count)%HeI_column_density_out(-partition(count):-1,1:partition(count)))  
      allocate(nx_py_nz_dx(-count)%HeI_column_density_out(1:partition(count),-partition(count):-1))
      allocate(nx_py_nz_dy(count)%HeI_column_density_out(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_py_nz_dz(-count)%HeI_column_density_out(-partition(count):-1,1:partition(count)))	  
      allocate(nx_ny_pz_dx(-count)%HeI_column_density_out(-partition(count):-1,1:partition(count)))
      allocate(nx_ny_pz_dy(-count)%HeI_column_density_out(1:partition(count),-partition(count):-1))	  
      allocate(nx_ny_pz_dz(count)%HeI_column_density_out(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dx(-count)%HeI_column_density_out(-partition(count):-1,-partition(count):-1))
      allocate(nx_ny_nz_dy(-count)%HeI_column_density_out(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dz(-count)%HeI_column_density_out(-partition(count):-1,-partition(count):-1))
	  	  	  
      allocate(px_py_pz_dx(count)%HeI_photoionization_rate(1:partition(count),1:partition(count)))
      allocate(px_py_pz_dy(count)%HeI_photoionization_rate(1:partition(count),1:partition(count)))	  
      allocate(px_py_pz_dz(count)%HeI_photoionization_rate(1:partition(count),1:partition(count)))	  
      allocate(px_py_nz_dx(count)%HeI_photoionization_rate(1:partition(count),-partition(count):-1))
      allocate(px_py_nz_dy(count)%HeI_photoionization_rate(-partition(count):-1,1:partition(count)))	  
      allocate(px_py_nz_dz(-count)%HeI_photoionization_rate(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dx(count)%HeI_photoionization_rate(-partition(count):-1,1:partition(count)))
      allocate(px_ny_pz_dy(-count)%HeI_photoionization_rate(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dz(count)%HeI_photoionization_rate(1:partition(count),-partition(count):-1))	  
      allocate(px_ny_nz_dx(count)%HeI_photoionization_rate(-partition(count):-1,-partition(count):-1))
      allocate(px_ny_nz_dy(-count)%HeI_photoionization_rate(-partition(count):-1,1:partition(count)))	  
      allocate(px_ny_nz_dz(-count)%HeI_photoionization_rate(1:partition(count),-partition(count):-1))	
      allocate(nx_py_pz_dx(-count)%HeI_photoionization_rate(1:partition(count),1:partition(count)))
      allocate(nx_py_pz_dy(count)%HeI_photoionization_rate(1:partition(count),-partition(count):-1))	  
      allocate(nx_py_pz_dz(count)%HeI_photoionization_rate(-partition(count):-1,1:partition(count)))  
      allocate(nx_py_nz_dx(-count)%HeI_photoionization_rate(1:partition(count),-partition(count):-1))
      allocate(nx_py_nz_dy(count)%HeI_photoionization_rate(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_py_nz_dz(-count)%HeI_photoionization_rate(-partition(count):-1,1:partition(count)))	  
      allocate(nx_ny_pz_dx(-count)%HeI_photoionization_rate(-partition(count):-1,1:partition(count)))
      allocate(nx_ny_pz_dy(-count)%HeI_photoionization_rate(1:partition(count),-partition(count):-1))	  
      allocate(nx_ny_pz_dz(count)%HeI_photoionization_rate(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dx(-count)%HeI_photoionization_rate(-partition(count):-1,-partition(count):-1))
      allocate(nx_ny_nz_dy(-count)%HeI_photoionization_rate(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dz(-count)%HeI_photoionization_rate(-partition(count):-1,-partition(count):-1))				  

      allocate(px_py_pz_dx(count)%HeII_density(1:partition(count),1:partition(count)))
      allocate(px_py_pz_dy(count)%HeII_density(1:partition(count),1:partition(count)))	  
      allocate(px_py_pz_dz(count)%HeII_density(1:partition(count),1:partition(count)))	  
      allocate(px_py_nz_dx(count)%HeII_density(1:partition(count),-partition(count):-1))
      allocate(px_py_nz_dy(count)%HeII_density(-partition(count):-1,1:partition(count)))	  
      allocate(px_py_nz_dz(-count)%HeII_density(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dx(count)%HeII_density(-partition(count):-1,1:partition(count)))
      allocate(px_ny_pz_dy(-count)%HeII_density(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dz(count)%HeII_density(1:partition(count),-partition(count):-1))	  
      allocate(px_ny_nz_dx(count)%HeII_density(-partition(count):-1,-partition(count):-1))
      allocate(px_ny_nz_dy(-count)%HeII_density(-partition(count):-1,1:partition(count)))	  
      allocate(px_ny_nz_dz(-count)%HeII_density(1:partition(count),-partition(count):-1))	
      allocate(nx_py_pz_dx(-count)%HeII_density(1:partition(count),1:partition(count)))
      allocate(nx_py_pz_dy(count)%HeII_density(1:partition(count),-partition(count):-1))	  
      allocate(nx_py_pz_dz(count)%HeII_density(-partition(count):-1,1:partition(count)))  
      allocate(nx_py_nz_dx(-count)%HeII_density(1:partition(count),-partition(count):-1))
      allocate(nx_py_nz_dy(count)%HeII_density(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_py_nz_dz(-count)%HeII_density(-partition(count):-1,1:partition(count)))	  
      allocate(nx_ny_pz_dx(-count)%HeII_density(-partition(count):-1,1:partition(count)))
      allocate(nx_ny_pz_dy(-count)%HeII_density(1:partition(count),-partition(count):-1))	  
      allocate(nx_ny_pz_dz(count)%HeII_density(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dx(-count)%HeII_density(-partition(count):-1,-partition(count):-1))
      allocate(nx_ny_nz_dy(-count)%HeII_density(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dz(-count)%HeII_density(-partition(count):-1,-partition(count):-1))	

      allocate(px_py_pz_dx(count)%HeII_column_density_in(1:partition(count),1:partition(count)))
      allocate(px_py_pz_dy(count)%HeII_column_density_in(1:partition(count),1:partition(count)))	  
      allocate(px_py_pz_dz(count)%HeII_column_density_in(1:partition(count),1:partition(count)))	  
      allocate(px_py_nz_dx(count)%HeII_column_density_in(1:partition(count),-partition(count):-1))
      allocate(px_py_nz_dy(count)%HeII_column_density_in(-partition(count):-1,1:partition(count)))	  
      allocate(px_py_nz_dz(-count)%HeII_column_density_in(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dx(count)%HeII_column_density_in(-partition(count):-1,1:partition(count)))
      allocate(px_ny_pz_dy(-count)%HeII_column_density_in(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dz(count)%HeII_column_density_in(1:partition(count),-partition(count):-1))	  
      allocate(px_ny_nz_dx(count)%HeII_column_density_in(-partition(count):-1,-partition(count):-1))
      allocate(px_ny_nz_dy(-count)%HeII_column_density_in(-partition(count):-1,1:partition(count)))	  
      allocate(px_ny_nz_dz(-count)%HeII_column_density_in(1:partition(count),-partition(count):-1))	
      allocate(nx_py_pz_dx(-count)%HeII_column_density_in(1:partition(count),1:partition(count)))
      allocate(nx_py_pz_dy(count)%HeII_column_density_in(1:partition(count),-partition(count):-1))	  
      allocate(nx_py_pz_dz(count)%HeII_column_density_in(-partition(count):-1,1:partition(count)))  
      allocate(nx_py_nz_dx(-count)%HeII_column_density_in(1:partition(count),-partition(count):-1))
      allocate(nx_py_nz_dy(count)%HeII_column_density_in(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_py_nz_dz(-count)%HeII_column_density_in(-partition(count):-1,1:partition(count)))	  
      allocate(nx_ny_pz_dx(-count)%HeII_column_density_in(-partition(count):-1,1:partition(count)))
      allocate(nx_ny_pz_dy(-count)%HeII_column_density_in(1:partition(count),-partition(count):-1))	  
      allocate(nx_ny_pz_dz(count)%HeII_column_density_in(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dx(-count)%HeII_column_density_in(-partition(count):-1,-partition(count):-1))
      allocate(nx_ny_nz_dy(-count)%HeII_column_density_in(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dz(-count)%HeII_column_density_in(-partition(count):-1,-partition(count):-1))
	  
      allocate(px_py_pz_dx(count)%HeII_column_density_out(1:partition(count),1:partition(count)))
      allocate(px_py_pz_dy(count)%HeII_column_density_out(1:partition(count),1:partition(count)))	  
      allocate(px_py_pz_dz(count)%HeII_column_density_out(1:partition(count),1:partition(count)))	  
      allocate(px_py_nz_dx(count)%HeII_column_density_out(1:partition(count),-partition(count):-1))
      allocate(px_py_nz_dy(count)%HeII_column_density_out(-partition(count):-1,1:partition(count)))	  
      allocate(px_py_nz_dz(-count)%HeII_column_density_out(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dx(count)%HeII_column_density_out(-partition(count):-1,1:partition(count)))
      allocate(px_ny_pz_dy(-count)%HeII_column_density_out(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dz(count)%HeII_column_density_out(1:partition(count),-partition(count):-1))	  
      allocate(px_ny_nz_dx(count)%HeII_column_density_out(-partition(count):-1,-partition(count):-1))
      allocate(px_ny_nz_dy(-count)%HeII_column_density_out(-partition(count):-1,1:partition(count)))	  
      allocate(px_ny_nz_dz(-count)%HeII_column_density_out(1:partition(count),-partition(count):-1))	
      allocate(nx_py_pz_dx(-count)%HeII_column_density_out(1:partition(count),1:partition(count)))
      allocate(nx_py_pz_dy(count)%HeII_column_density_out(1:partition(count),-partition(count):-1))	  
      allocate(nx_py_pz_dz(count)%HeII_column_density_out(-partition(count):-1,1:partition(count)))  
      allocate(nx_py_nz_dx(-count)%HeII_column_density_out(1:partition(count),-partition(count):-1))
      allocate(nx_py_nz_dy(count)%HeII_column_density_out(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_py_nz_dz(-count)%HeII_column_density_out(-partition(count):-1,1:partition(count)))	  
      allocate(nx_ny_pz_dx(-count)%HeII_column_density_out(-partition(count):-1,1:partition(count)))
      allocate(nx_ny_pz_dy(-count)%HeII_column_density_out(1:partition(count),-partition(count):-1))	  
      allocate(nx_ny_pz_dz(count)%HeII_column_density_out(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dx(-count)%HeII_column_density_out(-partition(count):-1,-partition(count):-1))
      allocate(nx_ny_nz_dy(-count)%HeII_column_density_out(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dz(-count)%HeII_column_density_out(-partition(count):-1,-partition(count):-1))
	  	  	  
      allocate(px_py_pz_dx(count)%HeII_photoionization_rate(1:partition(count),1:partition(count)))
      allocate(px_py_pz_dy(count)%HeII_photoionization_rate(1:partition(count),1:partition(count)))	  
      allocate(px_py_pz_dz(count)%HeII_photoionization_rate(1:partition(count),1:partition(count)))	  
      allocate(px_py_nz_dx(count)%HeII_photoionization_rate(1:partition(count),-partition(count):-1))
      allocate(px_py_nz_dy(count)%HeII_photoionization_rate(-partition(count):-1,1:partition(count)))	  
      allocate(px_py_nz_dz(-count)%HeII_photoionization_rate(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dx(count)%HeII_photoionization_rate(-partition(count):-1,1:partition(count)))
      allocate(px_ny_pz_dy(-count)%HeII_photoionization_rate(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dz(count)%HeII_photoionization_rate(1:partition(count),-partition(count):-1))	  
      allocate(px_ny_nz_dx(count)%HeII_photoionization_rate(-partition(count):-1,-partition(count):-1))
      allocate(px_ny_nz_dy(-count)%HeII_photoionization_rate(-partition(count):-1,1:partition(count)))	  
      allocate(px_ny_nz_dz(-count)%HeII_photoionization_rate(1:partition(count),-partition(count):-1))	
      allocate(nx_py_pz_dx(-count)%HeII_photoionization_rate(1:partition(count),1:partition(count)))
      allocate(nx_py_pz_dy(count)%HeII_photoionization_rate(1:partition(count),-partition(count):-1))	  
      allocate(nx_py_pz_dz(count)%HeII_photoionization_rate(-partition(count):-1,1:partition(count)))  
      allocate(nx_py_nz_dx(-count)%HeII_photoionization_rate(1:partition(count),-partition(count):-1))
      allocate(nx_py_nz_dy(count)%HeII_photoionization_rate(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_py_nz_dz(-count)%HeII_photoionization_rate(-partition(count):-1,1:partition(count)))	  
      allocate(nx_ny_pz_dx(-count)%HeII_photoionization_rate(-partition(count):-1,1:partition(count)))
      allocate(nx_ny_pz_dy(-count)%HeII_photoionization_rate(1:partition(count),-partition(count):-1))	  
      allocate(nx_ny_pz_dz(count)%HeII_photoionization_rate(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dx(-count)%HeII_photoionization_rate(-partition(count):-1,-partition(count):-1))
      allocate(nx_ny_nz_dy(-count)%HeII_photoionization_rate(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dz(-count)%HeII_photoionization_rate(-partition(count):-1,-partition(count):-1))				  

      allocate(px_py_pz_dx(count)%photoheating_rate(1:partition(count),1:partition(count)))
      allocate(px_py_pz_dy(count)%photoheating_rate(1:partition(count),1:partition(count)))	  
      allocate(px_py_pz_dz(count)%photoheating_rate(1:partition(count),1:partition(count)))	  
      allocate(px_py_nz_dx(count)%photoheating_rate(1:partition(count),-partition(count):-1))
      allocate(px_py_nz_dy(count)%photoheating_rate(-partition(count):-1,1:partition(count)))	  
      allocate(px_py_nz_dz(-count)%photoheating_rate(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dx(count)%photoheating_rate(-partition(count):-1,1:partition(count)))
      allocate(px_ny_pz_dy(-count)%photoheating_rate(1:partition(count),1:partition(count)))	  
      allocate(px_ny_pz_dz(count)%photoheating_rate(1:partition(count),-partition(count):-1))	  
      allocate(px_ny_nz_dx(count)%photoheating_rate(-partition(count):-1,-partition(count):-1))
      allocate(px_ny_nz_dy(-count)%photoheating_rate(-partition(count):-1,1:partition(count)))	  
      allocate(px_ny_nz_dz(-count)%photoheating_rate(1:partition(count),-partition(count):-1))	
      allocate(nx_py_pz_dx(-count)%photoheating_rate(1:partition(count),1:partition(count)))
      allocate(nx_py_pz_dy(count)%photoheating_rate(1:partition(count),-partition(count):-1))	  
      allocate(nx_py_pz_dz(count)%photoheating_rate(-partition(count):-1,1:partition(count)))  
      allocate(nx_py_nz_dx(-count)%photoheating_rate(1:partition(count),-partition(count):-1))
      allocate(nx_py_nz_dy(count)%photoheating_rate(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_py_nz_dz(-count)%photoheating_rate(-partition(count):-1,1:partition(count)))	  
      allocate(nx_ny_pz_dx(-count)%photoheating_rate(-partition(count):-1,1:partition(count)))
      allocate(nx_ny_pz_dy(-count)%photoheating_rate(1:partition(count),-partition(count):-1))	  
      allocate(nx_ny_pz_dz(count)%photoheating_rate(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dx(-count)%photoheating_rate(-partition(count):-1,-partition(count):-1))
      allocate(nx_ny_nz_dy(-count)%photoheating_rate(-partition(count):-1,-partition(count):-1))	  
      allocate(nx_ny_nz_dz(-count)%photoheating_rate(-partition(count):-1,-partition(count):-1))
	  
    end do

end subroutine domain_creation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pyramid_destruction()

    integer :: k

    ! To deallocate of the arrays
    do k = 1,level_py

      deallocate(pyramid_grid(k)%cellsize)
      deallocate(pyramid_grid(k)%solid_angle)
      deallocate(pyramid_grid(k)%case)
      deallocate(pyramid_grid(k)%transform)
	  	  	      
    end do

  end subroutine pyramid_destruction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cartiesian_destruction()

    integer :: k

    ! To deallocate of the arrays
    do k = 1,level_py

      deallocate(cartesian_grid(k)%transform)
	  	  	      
    end do
	
  end subroutine cartiesian_destruction
	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine domain_destruction()

    implicit none
	
    integer :: count

    ! To deallocate of the arrays
    do count = 1,level_py

      deallocate(px_py_pz_dx(count)%HI_density)
      deallocate(px_py_pz_dy(count)%HI_density)
      deallocate(px_py_pz_dz(count)%HI_density)
      deallocate(px_py_nz_dx(count)%HI_density)
      deallocate(px_py_nz_dy(count)%HI_density)  
      deallocate(px_py_nz_dz(-count)%HI_density)
      deallocate(px_ny_pz_dx(count)%HI_density)
      deallocate(px_ny_pz_dy(-count)%HI_density)
      deallocate(px_ny_pz_dz(count)%HI_density)
      deallocate(px_ny_nz_dx(count)%HI_density)
      deallocate(px_ny_nz_dy(-count)%HI_density)  
      deallocate(px_ny_nz_dz(-count)%HI_density)
      deallocate(nx_py_pz_dx(-count)%HI_density)
      deallocate(nx_py_pz_dy(count)%HI_density)	  
      deallocate(nx_py_pz_dz(count)%HI_density)	  
      deallocate(nx_py_nz_dx(-count)%HI_density)
      deallocate(nx_py_nz_dy(count)%HI_density)  
      deallocate(nx_py_nz_dz(-count)%HI_density)	  
      deallocate(nx_ny_pz_dx(-count)%HI_density)
      deallocate(nx_ny_pz_dy(-count)%HI_density) 
      deallocate(nx_ny_pz_dz(count)%HI_density) 
      deallocate(nx_ny_nz_dx(-count)%HI_density)
      deallocate(nx_ny_nz_dy(-count)%HI_density)	  
      deallocate(nx_ny_nz_dz(-count)%HI_density)

	  deallocate(px_py_pz_dx(count)%HeI_density)
	  deallocate(px_py_pz_dy(count)%HeI_density)
	  deallocate(px_py_pz_dz(count)%HeI_density)
	  deallocate(px_py_nz_dx(count)%HeI_density)
	  deallocate(px_py_nz_dy(count)%HeI_density)  
	  deallocate(px_py_nz_dz(-count)%HeI_density)
	  deallocate(px_ny_pz_dx(count)%HeI_density)
	  deallocate(px_ny_pz_dy(-count)%HeI_density)
	  deallocate(px_ny_pz_dz(count)%HeI_density)
	  deallocate(px_ny_nz_dx(count)%HeI_density)
	  deallocate(px_ny_nz_dy(-count)%HeI_density)  
	  deallocate(px_ny_nz_dz(-count)%HeI_density)
	  deallocate(nx_py_pz_dx(-count)%HeI_density)
	  deallocate(nx_py_pz_dy(count)%HeI_density)	  
	  deallocate(nx_py_pz_dz(count)%HeI_density)	  
	  deallocate(nx_py_nz_dx(-count)%HeI_density)
	  deallocate(nx_py_nz_dy(count)%HeI_density)  
	  deallocate(nx_py_nz_dz(-count)%HeI_density)	  
	  deallocate(nx_ny_pz_dx(-count)%HeI_density)
	  deallocate(nx_ny_pz_dy(-count)%HeI_density) 
	  deallocate(nx_ny_pz_dz(count)%HeI_density) 
	  deallocate(nx_ny_nz_dx(-count)%HeI_density)
	  deallocate(nx_ny_nz_dy(-count)%HeI_density)	  
	  deallocate(nx_ny_nz_dz(-count)%HeI_density)

	  deallocate(px_py_pz_dx(count)%HeII_density)
	  deallocate(px_py_pz_dy(count)%HeII_density)
	  deallocate(px_py_pz_dz(count)%HeII_density)
	  deallocate(px_py_nz_dx(count)%HeII_density)
	  deallocate(px_py_nz_dy(count)%HeII_density)  
	  deallocate(px_py_nz_dz(-count)%HeII_density)
	  deallocate(px_ny_pz_dx(count)%HeII_density)
	  deallocate(px_ny_pz_dy(-count)%HeII_density)
	  deallocate(px_ny_pz_dz(count)%HeII_density)
	  deallocate(px_ny_nz_dx(count)%HeII_density)
	  deallocate(px_ny_nz_dy(-count)%HeII_density)  
	  deallocate(px_ny_nz_dz(-count)%HeII_density)
	  deallocate(nx_py_pz_dx(-count)%HeII_density)
	  deallocate(nx_py_pz_dy(count)%HeII_density)	  
	  deallocate(nx_py_pz_dz(count)%HeII_density)	  
	  deallocate(nx_py_nz_dx(-count)%HeII_density)
	  deallocate(nx_py_nz_dy(count)%HeII_density)  
	  deallocate(nx_py_nz_dz(-count)%HeII_density)	  
	  deallocate(nx_ny_pz_dx(-count)%HeII_density)
	  deallocate(nx_ny_pz_dy(-count)%HeII_density) 
	  deallocate(nx_ny_pz_dz(count)%HeII_density) 
	  deallocate(nx_ny_nz_dx(-count)%HeII_density)
	  deallocate(nx_ny_nz_dy(-count)%HeII_density)	  
	  deallocate(nx_ny_nz_dz(-count)%HeII_density)

	  deallocate(px_py_pz_dx(count)%HI_column_density_in)
	  deallocate(px_py_pz_dy(count)%HI_column_density_in)
	  deallocate(px_py_pz_dz(count)%HI_column_density_in)
	  deallocate(px_py_nz_dx(count)%HI_column_density_in)
	  deallocate(px_py_nz_dy(count)%HI_column_density_in)  
	  deallocate(px_py_nz_dz(-count)%HI_column_density_in)
	  deallocate(px_ny_pz_dx(count)%HI_column_density_in)
	  deallocate(px_ny_pz_dy(-count)%HI_column_density_in)
	  deallocate(px_ny_pz_dz(count)%HI_column_density_in)
	  deallocate(px_ny_nz_dx(count)%HI_column_density_in)
	  deallocate(px_ny_nz_dy(-count)%HI_column_density_in)  
	  deallocate(px_ny_nz_dz(-count)%HI_column_density_in)
	  deallocate(nx_py_pz_dx(-count)%HI_column_density_in)
	  deallocate(nx_py_pz_dy(count)%HI_column_density_in)	  
	  deallocate(nx_py_pz_dz(count)%HI_column_density_in)	  
	  deallocate(nx_py_nz_dx(-count)%HI_column_density_in)
	  deallocate(nx_py_nz_dy(count)%HI_column_density_in)  
	  deallocate(nx_py_nz_dz(-count)%HI_column_density_in)	  
	  deallocate(nx_ny_pz_dx(-count)%HI_column_density_in)
	  deallocate(nx_ny_pz_dy(-count)%HI_column_density_in) 
	  deallocate(nx_ny_pz_dz(count)%HI_column_density_in) 
	  deallocate(nx_ny_nz_dx(-count)%HI_column_density_in)
	  deallocate(nx_ny_nz_dy(-count)%HI_column_density_in)	  
	  deallocate(nx_ny_nz_dz(-count)%HI_column_density_in)

	  deallocate(px_py_pz_dx(count)%HeI_column_density_in)
	  deallocate(px_py_pz_dy(count)%HeI_column_density_in)
	  deallocate(px_py_pz_dz(count)%HeI_column_density_in)
	  deallocate(px_py_nz_dx(count)%HeI_column_density_in)
	  deallocate(px_py_nz_dy(count)%HeI_column_density_in)  
	  deallocate(px_py_nz_dz(-count)%HeI_column_density_in)
	  deallocate(px_ny_pz_dx(count)%HeI_column_density_in)
	  deallocate(px_ny_pz_dy(-count)%HeI_column_density_in)
	  deallocate(px_ny_pz_dz(count)%HeI_column_density_in)
	  deallocate(px_ny_nz_dx(count)%HeI_column_density_in)
	  deallocate(px_ny_nz_dy(-count)%HeI_column_density_in)  
	  deallocate(px_ny_nz_dz(-count)%HeI_column_density_in)
	  deallocate(nx_py_pz_dx(-count)%HeI_column_density_in)
	  deallocate(nx_py_pz_dy(count)%HeI_column_density_in)	  
	  deallocate(nx_py_pz_dz(count)%HeI_column_density_in)	  
	  deallocate(nx_py_nz_dx(-count)%HeI_column_density_in)
	  deallocate(nx_py_nz_dy(count)%HeI_column_density_in)  
	  deallocate(nx_py_nz_dz(-count)%HeI_column_density_in)	  
	  deallocate(nx_ny_pz_dx(-count)%HeI_column_density_in)
	  deallocate(nx_ny_pz_dy(-count)%HeI_column_density_in) 
	  deallocate(nx_ny_pz_dz(count)%HeI_column_density_in) 
	  deallocate(nx_ny_nz_dx(-count)%HeI_column_density_in)
	  deallocate(nx_ny_nz_dy(-count)%HeI_column_density_in)	  
	  deallocate(nx_ny_nz_dz(-count)%HeI_column_density_in)

	  deallocate(px_py_pz_dx(count)%HeII_column_density_in)
	  deallocate(px_py_pz_dy(count)%HeII_column_density_in)
	  deallocate(px_py_pz_dz(count)%HeII_column_density_in)
	  deallocate(px_py_nz_dx(count)%HeII_column_density_in)
	  deallocate(px_py_nz_dy(count)%HeII_column_density_in)  
	  deallocate(px_py_nz_dz(-count)%HeII_column_density_in)
	  deallocate(px_ny_pz_dx(count)%HeII_column_density_in)
	  deallocate(px_ny_pz_dy(-count)%HeII_column_density_in)
	  deallocate(px_ny_pz_dz(count)%HeII_column_density_in)
	  deallocate(px_ny_nz_dx(count)%HeII_column_density_in)
	  deallocate(px_ny_nz_dy(-count)%HeII_column_density_in)  
	  deallocate(px_ny_nz_dz(-count)%HeII_column_density_in)
	  deallocate(nx_py_pz_dx(-count)%HeII_column_density_in)
	  deallocate(nx_py_pz_dy(count)%HeII_column_density_in)	  
	  deallocate(nx_py_pz_dz(count)%HeII_column_density_in)	  
	  deallocate(nx_py_nz_dx(-count)%HeII_column_density_in)
	  deallocate(nx_py_nz_dy(count)%HeII_column_density_in)  
	  deallocate(nx_py_nz_dz(-count)%HeII_column_density_in)	  
	  deallocate(nx_ny_pz_dx(-count)%HeII_column_density_in)
	  deallocate(nx_ny_pz_dy(-count)%HeII_column_density_in) 
	  deallocate(nx_ny_pz_dz(count)%HeII_column_density_in) 
	  deallocate(nx_ny_nz_dx(-count)%HeII_column_density_in)
	  deallocate(nx_ny_nz_dy(-count)%HeII_column_density_in)	  
	  deallocate(nx_ny_nz_dz(-count)%HeII_column_density_in)
	  	  
      deallocate(px_py_pz_dx(count)%HI_column_density_out)
      deallocate(px_py_pz_dy(count)%HI_column_density_out)
      deallocate(px_py_pz_dz(count)%HI_column_density_out)
      deallocate(px_py_nz_dx(count)%HI_column_density_out)
      deallocate(px_py_nz_dy(count)%HI_column_density_out)  
      deallocate(px_py_nz_dz(-count)%HI_column_density_out)
      deallocate(px_ny_pz_dx(count)%HI_column_density_out)
      deallocate(px_ny_pz_dy(-count)%HI_column_density_out)
      deallocate(px_ny_pz_dz(count)%HI_column_density_out)
      deallocate(px_ny_nz_dx(count)%HI_column_density_out)
      deallocate(px_ny_nz_dy(-count)%HI_column_density_out)  
      deallocate(px_ny_nz_dz(-count)%HI_column_density_out)
      deallocate(nx_py_pz_dx(-count)%HI_column_density_out)
      deallocate(nx_py_pz_dy(count)%HI_column_density_out)	  
      deallocate(nx_py_pz_dz(count)%HI_column_density_out)	  
      deallocate(nx_py_nz_dx(-count)%HI_column_density_out)
      deallocate(nx_py_nz_dy(count)%HI_column_density_out)  
      deallocate(nx_py_nz_dz(-count)%HI_column_density_out)	  
      deallocate(nx_ny_pz_dx(-count)%HI_column_density_out)
      deallocate(nx_ny_pz_dy(-count)%HI_column_density_out) 
      deallocate(nx_ny_pz_dz(count)%HI_column_density_out) 
      deallocate(nx_ny_nz_dx(-count)%HI_column_density_out)
      deallocate(nx_ny_nz_dy(-count)%HI_column_density_out)	  
      deallocate(nx_ny_nz_dz(-count)%HI_column_density_out)

	  deallocate(px_py_pz_dx(count)%HeI_column_density_out)
	  deallocate(px_py_pz_dy(count)%HeI_column_density_out)
	  deallocate(px_py_pz_dz(count)%HeI_column_density_out)
	  deallocate(px_py_nz_dx(count)%HeI_column_density_out)
	  deallocate(px_py_nz_dy(count)%HeI_column_density_out)  
	  deallocate(px_py_nz_dz(-count)%HeI_column_density_out)
	  deallocate(px_ny_pz_dx(count)%HeI_column_density_out)
	  deallocate(px_ny_pz_dy(-count)%HeI_column_density_out)
	  deallocate(px_ny_pz_dz(count)%HeI_column_density_out)
	  deallocate(px_ny_nz_dx(count)%HeI_column_density_out)
	  deallocate(px_ny_nz_dy(-count)%HeI_column_density_out)  
	  deallocate(px_ny_nz_dz(-count)%HeI_column_density_out)
	  deallocate(nx_py_pz_dx(-count)%HeI_column_density_out)
	  deallocate(nx_py_pz_dy(count)%HeI_column_density_out)	  
	  deallocate(nx_py_pz_dz(count)%HeI_column_density_out)	  
	  deallocate(nx_py_nz_dx(-count)%HeI_column_density_out)
	  deallocate(nx_py_nz_dy(count)%HeI_column_density_out)  
	  deallocate(nx_py_nz_dz(-count)%HeI_column_density_out)	  
	  deallocate(nx_ny_pz_dx(-count)%HeI_column_density_out)
	  deallocate(nx_ny_pz_dy(-count)%HeI_column_density_out) 
	  deallocate(nx_ny_pz_dz(count)%HeI_column_density_out) 
	  deallocate(nx_ny_nz_dx(-count)%HeI_column_density_out)
	  deallocate(nx_ny_nz_dy(-count)%HeI_column_density_out)	  
	  deallocate(nx_ny_nz_dz(-count)%HeI_column_density_out)

	  deallocate(px_py_pz_dx(count)%HeII_column_density_out)
	  deallocate(px_py_pz_dy(count)%HeII_column_density_out)
	  deallocate(px_py_pz_dz(count)%HeII_column_density_out)
	  deallocate(px_py_nz_dx(count)%HeII_column_density_out)
	  deallocate(px_py_nz_dy(count)%HeII_column_density_out)  
	  deallocate(px_py_nz_dz(-count)%HeII_column_density_out)
	  deallocate(px_ny_pz_dx(count)%HeII_column_density_out)
	  deallocate(px_ny_pz_dy(-count)%HeII_column_density_out)
	  deallocate(px_ny_pz_dz(count)%HeII_column_density_out)
	  deallocate(px_ny_nz_dx(count)%HeII_column_density_out)
	  deallocate(px_ny_nz_dy(-count)%HeII_column_density_out)  
	  deallocate(px_ny_nz_dz(-count)%HeII_column_density_out)
	  deallocate(nx_py_pz_dx(-count)%HeII_column_density_out)
	  deallocate(nx_py_pz_dy(count)%HeII_column_density_out)	  
	  deallocate(nx_py_pz_dz(count)%HeII_column_density_out)	  
	  deallocate(nx_py_nz_dx(-count)%HeII_column_density_out)
	  deallocate(nx_py_nz_dy(count)%HeII_column_density_out)  
	  deallocate(nx_py_nz_dz(-count)%HeII_column_density_out)	  
	  deallocate(nx_ny_pz_dx(-count)%HeII_column_density_out)
	  deallocate(nx_ny_pz_dy(-count)%HeII_column_density_out) 
	  deallocate(nx_ny_pz_dz(count)%HeII_column_density_out) 
	  deallocate(nx_ny_nz_dx(-count)%HeII_column_density_out)
	  deallocate(nx_ny_nz_dy(-count)%HeII_column_density_out)	  
	  deallocate(nx_ny_nz_dz(-count)%HeII_column_density_out)
	  	  
      deallocate(px_py_pz_dx(count)%HI_photoionization_rate)
      deallocate(px_py_pz_dy(count)%HI_photoionization_rate)
      deallocate(px_py_pz_dz(count)%HI_photoionization_rate)
      deallocate(px_py_nz_dx(count)%HI_photoionization_rate)
      deallocate(px_py_nz_dy(count)%HI_photoionization_rate)  
      deallocate(px_py_nz_dz(-count)%HI_photoionization_rate)
      deallocate(px_ny_pz_dx(count)%HI_photoionization_rate)
      deallocate(px_ny_pz_dy(-count)%HI_photoionization_rate)
      deallocate(px_ny_pz_dz(count)%HI_photoionization_rate)
      deallocate(px_ny_nz_dx(count)%HI_photoionization_rate)
      deallocate(px_ny_nz_dy(-count)%HI_photoionization_rate)  
      deallocate(px_ny_nz_dz(-count)%HI_photoionization_rate)
      deallocate(nx_py_pz_dx(-count)%HI_photoionization_rate)
      deallocate(nx_py_pz_dy(count)%HI_photoionization_rate)	  
      deallocate(nx_py_pz_dz(count)%HI_photoionization_rate)	  
      deallocate(nx_py_nz_dx(-count)%HI_photoionization_rate)
      deallocate(nx_py_nz_dy(count)%HI_photoionization_rate)  
      deallocate(nx_py_nz_dz(-count)%HI_photoionization_rate)	  
      deallocate(nx_ny_pz_dx(-count)%HI_photoionization_rate)
      deallocate(nx_ny_pz_dy(-count)%HI_photoionization_rate) 
      deallocate(nx_ny_pz_dz(count)%HI_photoionization_rate) 
      deallocate(nx_ny_nz_dx(-count)%HI_photoionization_rate)
      deallocate(nx_ny_nz_dy(-count)%HI_photoionization_rate)	  
      deallocate(nx_ny_nz_dz(-count)%HI_photoionization_rate)	  

	  deallocate(px_py_pz_dx(count)%HeI_photoionization_rate)
	  deallocate(px_py_pz_dy(count)%HeI_photoionization_rate)
	  deallocate(px_py_pz_dz(count)%HeI_photoionization_rate)
	  deallocate(px_py_nz_dx(count)%HeI_photoionization_rate)
	  deallocate(px_py_nz_dy(count)%HeI_photoionization_rate)  
	  deallocate(px_py_nz_dz(-count)%HeI_photoionization_rate)
	  deallocate(px_ny_pz_dx(count)%HeI_photoionization_rate)
	  deallocate(px_ny_pz_dy(-count)%HeI_photoionization_rate)
	  deallocate(px_ny_pz_dz(count)%HeI_photoionization_rate)
	  deallocate(px_ny_nz_dx(count)%HeI_photoionization_rate)
	  deallocate(px_ny_nz_dy(-count)%HeI_photoionization_rate)  
	  deallocate(px_ny_nz_dz(-count)%HeI_photoionization_rate)
	  deallocate(nx_py_pz_dx(-count)%HeI_photoionization_rate)
	  deallocate(nx_py_pz_dy(count)%HeI_photoionization_rate)	  
	  deallocate(nx_py_pz_dz(count)%HeI_photoionization_rate)	  
	  deallocate(nx_py_nz_dx(-count)%HeI_photoionization_rate)
	  deallocate(nx_py_nz_dy(count)%HeI_photoionization_rate)  
	  deallocate(nx_py_nz_dz(-count)%HeI_photoionization_rate)	  
	  deallocate(nx_ny_pz_dx(-count)%HeI_photoionization_rate)
	  deallocate(nx_ny_pz_dy(-count)%HeI_photoionization_rate) 
	  deallocate(nx_ny_pz_dz(count)%HeI_photoionization_rate) 
	  deallocate(nx_ny_nz_dx(-count)%HeI_photoionization_rate)
	  deallocate(nx_ny_nz_dy(-count)%HeI_photoionization_rate)	  
	  deallocate(nx_ny_nz_dz(-count)%HeI_photoionization_rate)

	  deallocate(px_py_pz_dx(count)%HeII_photoionization_rate)
	  deallocate(px_py_pz_dy(count)%HeII_photoionization_rate)
	  deallocate(px_py_pz_dz(count)%HeII_photoionization_rate)
	  deallocate(px_py_nz_dx(count)%HeII_photoionization_rate)
	  deallocate(px_py_nz_dy(count)%HeII_photoionization_rate)  
	  deallocate(px_py_nz_dz(-count)%HeII_photoionization_rate)
	  deallocate(px_ny_pz_dx(count)%HeII_photoionization_rate)
	  deallocate(px_ny_pz_dy(-count)%HeII_photoionization_rate)
	  deallocate(px_ny_pz_dz(count)%HeII_photoionization_rate)
	  deallocate(px_ny_nz_dx(count)%HeII_photoionization_rate)
	  deallocate(px_ny_nz_dy(-count)%HeII_photoionization_rate)  
	  deallocate(px_ny_nz_dz(-count)%HeII_photoionization_rate)
	  deallocate(nx_py_pz_dx(-count)%HeII_photoionization_rate)
	  deallocate(nx_py_pz_dy(count)%HeII_photoionization_rate)	  
	  deallocate(nx_py_pz_dz(count)%HeII_photoionization_rate)	  
	  deallocate(nx_py_nz_dx(-count)%HeII_photoionization_rate)
	  deallocate(nx_py_nz_dy(count)%HeII_photoionization_rate)  
	  deallocate(nx_py_nz_dz(-count)%HeII_photoionization_rate)	  
	  deallocate(nx_ny_pz_dx(-count)%HeII_photoionization_rate)
	  deallocate(nx_ny_pz_dy(-count)%HeII_photoionization_rate) 
	  deallocate(nx_ny_pz_dz(count)%HeII_photoionization_rate) 
	  deallocate(nx_ny_nz_dx(-count)%HeII_photoionization_rate)
	  deallocate(nx_ny_nz_dy(-count)%HeII_photoionization_rate)	  
	  deallocate(nx_ny_nz_dz(-count)%HeII_photoionization_rate)
	  	  
	  deallocate(px_py_pz_dx(count)%photoheating_rate)
	  deallocate(px_py_pz_dy(count)%photoheating_rate)
	  deallocate(px_py_pz_dz(count)%photoheating_rate)
	  deallocate(px_py_nz_dx(count)%photoheating_rate)
	  deallocate(px_py_nz_dy(count)%photoheating_rate)  
	  deallocate(px_py_nz_dz(-count)%photoheating_rate)
	  deallocate(px_ny_pz_dx(count)%photoheating_rate)
	  deallocate(px_ny_pz_dy(-count)%photoheating_rate)
	  deallocate(px_ny_pz_dz(count)%photoheating_rate)
	  deallocate(px_ny_nz_dx(count)%photoheating_rate)
	  deallocate(px_ny_nz_dy(-count)%photoheating_rate)  
	  deallocate(px_ny_nz_dz(-count)%photoheating_rate)
	  deallocate(nx_py_pz_dx(-count)%photoheating_rate)
	  deallocate(nx_py_pz_dy(count)%photoheating_rate)	  
	  deallocate(nx_py_pz_dz(count)%photoheating_rate)	  
	  deallocate(nx_py_nz_dx(-count)%photoheating_rate)
	  deallocate(nx_py_nz_dy(count)%photoheating_rate)  
	  deallocate(nx_py_nz_dz(-count)%photoheating_rate)	  
	  deallocate(nx_ny_pz_dx(-count)%photoheating_rate)
	  deallocate(nx_ny_pz_dy(-count)%photoheating_rate) 
	  deallocate(nx_ny_pz_dz(count)%photoheating_rate) 
	  deallocate(nx_ny_nz_dx(-count)%photoheating_rate)
	  deallocate(nx_ny_nz_dy(-count)%photoheating_rate)	  
	  deallocate(nx_ny_nz_dz(-count)%photoheating_rate)	  
	  	  
    end do

  end subroutine domain_destruction

end module array
