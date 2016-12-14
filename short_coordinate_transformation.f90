module short_coordinate_transformation

  use precision, only: dp
  use input, only: level_py, level_sl, source_position
  use array, only: short_global_number_density_array, &
                   short_global_xHI_array, short_global_xHeI_array, short_global_xHeII_array,&
                   short_source_HI_density_array, short_source_HeI_density_array, short_source_HeII_density_array, &
                   short_global_HI_photoionization_rate_array, short_global_HeI_photoionization_rate_array, &
                   short_global_HeII_photoionization_rate_array, short_global_photoheating_rate_array, &				   				   
				   short_source_HI_photoionization_rate_array, short_source_HeI_photoionization_rate_array, &
				   short_source_HeII_photoionization_rate_array, short_source_photoheating_rate_array
				   				   				     
contains
       
  subroutine short_global_to_source_transformation()	
	    
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
          
          short_source_HI_density_array(i_source,j_source,k_source) = &
          short_global_number_density_array(i_global,j_global,k_global)* &
		  short_global_xHI_array(i_global,j_global,k_global)
          short_source_HeI_density_array(i_source,j_source,k_source) = &
          short_global_number_density_array(i_global,j_global,k_global)* &
		  short_global_xHeI_array(i_global,j_global,k_global)
          short_source_HeII_density_array(i_source,j_source,k_source) = &
          short_global_number_density_array(i_global,j_global,k_global)* &
		  short_global_xHeII_array(i_global,j_global,k_global)		  
		  				  
        enddo
      enddo
    enddo
	
  end subroutine short_global_to_source_transformation

  subroutine short_source_to_global_transformation()	
	    
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
          
		  short_global_HI_photoionization_rate_array(i_global,j_global,k_global) = &
		  	short_source_HI_photoionization_rate_array(i_source,j_source,k_source)
  		  short_global_HeI_photoionization_rate_array(i_global,j_global,k_global) = &
  		  	short_source_HeI_photoionization_rate_array(i_source,j_source,k_source)
  		  short_global_HeII_photoionization_rate_array(i_global,j_global,k_global) = &
  		  	short_source_HeII_photoionization_rate_array(i_source,j_source,k_source)						
          short_global_photoheating_rate_array(i_global,j_global,k_global) = &
  		  	short_source_photoheating_rate_array(i_source,j_source,k_source)			
							  
        enddo
      enddo
    enddo
	
  end subroutine short_source_to_global_transformation
   
end module short_coordinate_transformation