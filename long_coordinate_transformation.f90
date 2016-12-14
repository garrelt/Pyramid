module long_coordinate_transformation

  use precision, only: dp
  use input, only: level_py, level_sl, source_position
  use array, only: long_global_number_density_array, &
                   long_global_xHI_array, long_global_xHeI_array, long_global_xHeII_array,&
                   long_source_HI_density_array, long_source_HeI_density_array, long_source_HeII_density_array, &
                   long_global_HI_photoionization_rate_array, long_global_HeI_photoionization_rate_array, &
                   long_global_HeII_photoionization_rate_array, long_global_photoheating_rate_array, &				   				   
                   long_source_HI_photoionization_rate_array, long_source_HeI_photoionization_rate_array, &
                   long_source_HeII_photoionization_rate_array, long_source_photoheating_rate_array
				     
contains
       
  subroutine long_global_to_source_transformation()	
	    
    implicit none
		
    integer :: i_source,j_source,k_source
    integer :: i_global,j_global,k_global
    integer :: source_position_x
    integer :: source_position_y		
    integer :: source_position_z
  	integer :: boxsize  
		
    source_position_x = source_position(1,1)
    source_position_y = source_position(2,1)	
    source_position_z = source_position(3,1)
	
    do i_source = -level_py,level_py
      do j_source = -level_py,level_py
        do k_source = -level_py,level_py
			 
          i_global = modulo(i_source+source_position_x-1,level_sl)+1
          j_global = modulo(j_source+source_position_y-1,level_sl)+1
          k_global = modulo(k_source+source_position_z-1,level_sl)+1
          
		  long_source_HI_density_array(i_source,j_source,k_source) = &
		  long_global_number_density_array(i_global,j_global,k_global)* &
		  long_global_xHI_array(i_global,j_global,k_global)
		  long_source_HeI_density_array(i_source,j_source,k_source) = &
		  long_global_number_density_array(i_global,j_global,k_global)* &
		  long_global_xHeI_array(i_global,j_global,k_global)
		  long_source_HeII_density_array(i_source,j_source,k_source) = &
		  long_global_number_density_array(i_global,j_global,k_global)* &
		  long_global_xHeII_array(i_global,j_global,k_global)
				  
        enddo
      enddo
    enddo
	
  end subroutine long_global_to_source_transformation

  subroutine long_source_to_global_transformation()	
	    
    implicit none
		
    integer :: i_source,j_source,k_source
    integer :: i_global,j_global,k_global
    integer :: source_position_x
    integer :: source_position_y		
    integer :: source_position_z
  	integer :: boxsize  
		
    source_position_x = source_position(1,1)
    source_position_y = source_position(2,1)	
    source_position_z = source_position(3,1)
	
    do i_source = -level_py,level_py
      do j_source = -level_py,level_py
        do k_source = -level_py,level_py
			 
          i_global = modulo(i_source+source_position_x-1,level_sl)+1
          j_global = modulo(j_source+source_position_y-1,level_sl)+1
          k_global = modulo(k_source+source_position_z-1,level_sl)+1
		  
		  long_global_HI_photoionization_rate_array(i_global,j_global,k_global) = &
		  	long_source_HI_photoionization_rate_array(i_source,j_source,k_source)
		  long_global_HeI_photoionization_rate_array(i_global,j_global,k_global) = &
		  	long_source_HeI_photoionization_rate_array(i_source,j_source,k_source)
		  long_global_HeII_photoionization_rate_array(i_global,j_global,k_global) = &
		  	long_source_HeII_photoionization_rate_array(i_source,j_source,k_source)						
		  long_global_photoheating_rate_array(i_global,j_global,k_global) = &
		  	long_source_photoheating_rate_array(i_source,j_source,k_source)
							  
        enddo
      enddo
    enddo
	
  end subroutine long_source_to_global_transformation
   
end module long_coordinate_transformation