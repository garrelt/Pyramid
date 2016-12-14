module coordinate_transformation_sl

  use precision, only: dp
  use input, only: level_py, level_sl, source_position
  use array, only: global_sl_number_density_array, global_sl_xHI_array, &
                   source_sl_HI_density_array, global_short_photo_rate_array,&
				   source_short_HI_photo_rate_array,  global_long_photo_rate_array,&
				   source_long_HI_photo_rate_array
				     
contains
       
  subroutine global_to_source_transformation_sl()	
	    
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
          
          source_sl_HI_density_array(i_source,j_source,k_source) = &
          global_sl_number_density_array(i_global,j_global,k_global)* &
		  global_sl_xHI_array(i_global,j_global,k_global)
				  
        enddo
      enddo
    enddo
	
  end subroutine global_to_source_transformation_sl

  subroutine source_to_global_transformation_short()	
	    
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
          
		  global_short_photo_rate_array(i_global,j_global,k_global) = &
		  	source_short_HI_photo_rate_array(i_source,j_source,k_source)

				  
        enddo
      enddo
    enddo
	
  end subroutine source_to_global_transformation_short

  subroutine source_to_global_transformation_long()	
	    
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
          
		  global_long_photo_rate_array(i_global,j_global,k_global) = &
		  	source_long_HI_photo_rate_array(i_source,j_source,k_source)

				  
        enddo
      enddo
    enddo
	
  end subroutine source_to_global_transformation_long
   
end module coordinate_transformation_sl