module pyramid_global_source_transformation

  use precision, only: dp
  use input, only: level_py
  use input, only: source_position
  use array, only: pyramid_global_number_density_array, &
                   pyramid_global_xHI_array, &
                   pyramid_global_xHeI_array, &
                   pyramid_global_xHeII_array, &
                   pyramid_source_px_py_pz_HI_density_array, &
                   pyramid_source_px_py_nz_HI_density_array, &				   
                   pyramid_source_px_ny_pz_HI_density_array, &
                   pyramid_source_px_ny_nz_HI_density_array, &			   
                   pyramid_source_nx_py_pz_HI_density_array, &
                   pyramid_source_nx_py_nz_HI_density_array, &				   
                   pyramid_source_nx_ny_pz_HI_density_array, &					   
                   pyramid_source_nx_ny_nz_HI_density_array, &
                   pyramid_source_px_py_pz_HeI_density_array, &
                   pyramid_source_px_py_nz_HeI_density_array, &				   
                   pyramid_source_px_ny_pz_HeI_density_array, &
                   pyramid_source_px_ny_nz_HeI_density_array, &			   
                   pyramid_source_nx_py_pz_HeI_density_array, &
                   pyramid_source_nx_py_nz_HeI_density_array, &				   
                   pyramid_source_nx_ny_pz_HeI_density_array, &					   
                   pyramid_source_nx_ny_nz_HeI_density_array, &				   
                   pyramid_source_px_py_pz_HeII_density_array, &
                   pyramid_source_px_py_nz_HeII_density_array, &				   
                   pyramid_source_px_ny_pz_HeII_density_array, &
                   pyramid_source_px_ny_nz_HeII_density_array, &			   
                   pyramid_source_nx_py_pz_HeII_density_array, &
                   pyramid_source_nx_py_nz_HeII_density_array, &				   
                   pyramid_source_nx_ny_pz_HeII_density_array, &					   
                   pyramid_source_nx_ny_nz_HeII_density_array, &				   				   		   				   				   				   
                   pyramid_source_px_py_pz_HI_photoionization_rate_array, &
                   pyramid_source_px_py_nz_HI_photoionization_rate_array, &				   
                   pyramid_source_px_ny_pz_HI_photoionization_rate_array, &
                   pyramid_source_px_ny_nz_HI_photoionization_rate_array, &			   
                   pyramid_source_nx_py_pz_HI_photoionization_rate_array, &
                   pyramid_source_nx_py_nz_HI_photoionization_rate_array, &				   
                   pyramid_source_nx_ny_pz_HI_photoionization_rate_array, &					   
                   pyramid_source_nx_ny_nz_HI_photoionization_rate_array, &
                   pyramid_source_px_py_pz_HeI_photoionization_rate_array, &
                   pyramid_source_px_py_nz_HeI_photoionization_rate_array, &				   
                   pyramid_source_px_ny_pz_HeI_photoionization_rate_array, &
                   pyramid_source_px_ny_nz_HeI_photoionization_rate_array, &			   
                   pyramid_source_nx_py_pz_HeI_photoionization_rate_array, &
                   pyramid_source_nx_py_nz_HeI_photoionization_rate_array, &				   
                   pyramid_source_nx_ny_pz_HeI_photoionization_rate_array, &					   
                   pyramid_source_nx_ny_nz_HeI_photoionization_rate_array, &				   
                   pyramid_source_px_py_pz_HeII_photoionization_rate_array, &
                   pyramid_source_px_py_nz_HeII_photoionization_rate_array, &				   
                   pyramid_source_px_ny_pz_HeII_photoionization_rate_array, &
                   pyramid_source_px_ny_nz_HeII_photoionization_rate_array, &			   
                   pyramid_source_nx_py_pz_HeII_photoionization_rate_array, &
                   pyramid_source_nx_py_nz_HeII_photoionization_rate_array, &				   
                   pyramid_source_nx_ny_pz_HeII_photoionization_rate_array, &					   
                   pyramid_source_nx_ny_nz_HeII_photoionization_rate_array, &				   				   								   
                   pyramid_source_px_py_pz_photoheating_rate_array, &
                   pyramid_source_px_py_nz_photoheating_rate_array, &
                   pyramid_source_px_ny_pz_photoheating_rate_array, &
                   pyramid_source_px_ny_nz_photoheating_rate_array, &
                   pyramid_source_nx_py_pz_photoheating_rate_array, &
                   pyramid_source_nx_py_nz_photoheating_rate_array, &
                   pyramid_source_nx_ny_pz_photoheating_rate_array, &
                   pyramid_source_nx_ny_nz_photoheating_rate_array, &
				   
                   pyramid_global_HI_photoionization_rate_array, &
                   pyramid_global_HeI_photoionization_rate_array, &	
                   pyramid_global_HeII_photoionization_rate_array, &				   			   
                   pyramid_global_photoheating_rate_array
				   	   	
				     
contains
	
  subroutine pyramid_global_to_source_species_density_transformation()
	  
    call global_to_source_pos_x_pos_y_pos_z_transformation()
    call global_to_source_pos_x_pos_y_neg_z_transformation()
    call global_to_source_pos_x_neg_y_pos_z_transformation()
    call global_to_source_pos_x_neg_y_neg_z_transformation()
    call global_to_source_neg_x_pos_y_pos_z_transformation()
    call global_to_source_neg_x_pos_y_neg_z_transformation()
    call global_to_source_neg_x_neg_y_pos_z_transformation()
    call global_to_source_neg_x_neg_y_neg_z_transformation()
	
  end subroutine pyramid_global_to_source_species_density_transformation	  

  subroutine pyramid_source_to_global_photo_rate_transformation()
	  
    call source_to_global_pos_x_pos_y_pos_z_transformation()
    call source_to_global_pos_x_pos_y_neg_z_transformation()
    call source_to_global_pos_x_neg_y_pos_z_transformation()
    call source_to_global_pos_x_neg_y_neg_z_transformation()
    call source_to_global_neg_x_pos_y_pos_z_transformation()
    call source_to_global_neg_x_pos_y_neg_z_transformation()
    call source_to_global_neg_x_neg_y_pos_z_transformation()
    call source_to_global_neg_x_neg_y_neg_z_transformation()
	  
  end subroutine pyramid_source_to_global_photo_rate_transformation 
    	  
  subroutine global_to_source_pos_x_pos_y_pos_z_transformation()	
	    
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
    boxsize = 2*level_py

    do i_source = 1,level_py
      do j_source = 1,level_py
        do k_source = 1,level_py
			 
          i_global = modulo(i_source+source_position_x-2,boxsize)+1
          j_global = modulo(j_source+source_position_y-2,boxsize)+1
          k_global = modulo(k_source+source_position_z-2,boxsize)+1
          
          pyramid_source_px_py_pz_HI_density_array(i_source,j_source,k_source) = &
                        pyramid_global_number_density_array(i_global,j_global,k_global) * &
                        pyramid_global_xHI_array(i_global,j_global,k_global)
          pyramid_source_px_py_pz_HeI_density_array(i_source,j_source,k_source) = &
                        pyramid_global_number_density_array(i_global,j_global,k_global) * &
                        pyramid_global_xHeI_array(i_global,j_global,k_global)
          pyramid_source_px_py_pz_HeII_density_array(i_source,j_source,k_source) = &
                        pyramid_global_number_density_array(i_global,j_global,k_global) * &
                        pyramid_global_xHeII_array(i_global,j_global,k_global)
	  
        enddo
      enddo
    enddo

  end subroutine global_to_source_pos_x_pos_y_pos_z_transformation

  subroutine global_to_source_pos_x_pos_y_neg_z_transformation()	
	    
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
    boxsize = 2*level_py
	
    do i_source = 1,level_py
      do j_source = 1,level_py
        do k_source = -1,-level_py,-1
			 
          i_global = modulo(i_source+source_position_x-2,boxsize)+1
          j_global = modulo(j_source+source_position_y-2,boxsize)+1
          k_global = modulo(k_source+source_position_z-1,boxsize)+1
          
          pyramid_source_px_py_nz_HI_density_array(i_source,j_source,k_source) = &
                        pyramid_global_number_density_array(i_global,j_global,k_global) * &
                        pyramid_global_xHI_array(i_global,j_global,k_global)
          pyramid_source_px_py_nz_HeI_density_array(i_source,j_source,k_source) = &
                        pyramid_global_number_density_array(i_global,j_global,k_global) * &
                        pyramid_global_xHeI_array(i_global,j_global,k_global)				  
          pyramid_source_px_py_nz_HeII_density_array(i_source,j_source,k_source) = &
                        pyramid_global_number_density_array(i_global,j_global,k_global) * &
                        pyramid_global_xHeII_array(i_global,j_global,k_global)

        enddo
      enddo
    enddo

  end subroutine global_to_source_pos_x_pos_y_neg_z_transformation
   
  subroutine global_to_source_pos_x_neg_y_pos_z_transformation()	
	    
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
    boxsize = 2*level_py
			
    do i_source = 1,level_py
      do j_source = -1,-level_py,-1
        do k_source = 1,level_py
			 
          i_global = modulo(i_source+source_position_x-2,boxsize)+1
          j_global = modulo(j_source+source_position_y-1,boxsize)+1
          k_global = modulo(k_source+source_position_z-2,boxsize)+1
          
          pyramid_source_px_ny_pz_HI_density_array(i_source,j_source,k_source) = &
                        pyramid_global_number_density_array(i_global,j_global,k_global) * &
                        pyramid_global_xHI_array(i_global,j_global,k_global)
          pyramid_source_px_ny_pz_HeI_density_array(i_source,j_source,k_source) = &
                        pyramid_global_number_density_array(i_global,j_global,k_global) * &
                        pyramid_global_xHeI_array(i_global,j_global,k_global)
          pyramid_source_px_ny_pz_HeII_density_array(i_source,j_source,k_source) = &
                        pyramid_global_number_density_array(i_global,j_global,k_global) * &
                        pyramid_global_xHeII_array(i_global,j_global,k_global)
				  
        enddo
      enddo
    enddo

  end subroutine global_to_source_pos_x_neg_y_pos_z_transformation

  subroutine global_to_source_pos_x_neg_y_neg_z_transformation()	
	    
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
    boxsize = 2*level_py
			
    do i_source = 1,level_py
      do j_source = -1,-level_py,-1
        do k_source = -1,-level_py,-1
			 
          i_global = modulo(i_source+source_position_x-2,boxsize)+1
          j_global = modulo(j_source+source_position_y-1,boxsize)+1
          k_global = modulo(k_source+source_position_z-1,boxsize)+1
          
          pyramid_source_px_ny_nz_HI_density_array(i_source,j_source,k_source) = &
                        pyramid_global_number_density_array(i_global,j_global,k_global) * &
                        pyramid_global_xHI_array(i_global,j_global,k_global)
          pyramid_source_px_ny_nz_HeI_density_array(i_source,j_source,k_source) = &
                        pyramid_global_number_density_array(i_global,j_global,k_global) * &
                        pyramid_global_xHeI_array(i_global,j_global,k_global)
          pyramid_source_px_ny_nz_HeII_density_array(i_source,j_source,k_source) = &
                        pyramid_global_number_density_array(i_global,j_global,k_global) * &
                        pyramid_global_xHeII_array(i_global,j_global,k_global)
				  
        enddo
      enddo
    enddo

  end subroutine global_to_source_pos_x_neg_y_neg_z_transformation

  subroutine global_to_source_neg_x_pos_y_pos_z_transformation()	
	    
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
    boxsize = 2*level_py
			
    do i_source = -1,-level_py,-1
      do j_source = 1,level_py
        do k_source = 1,level_py
			 
          i_global = modulo(i_source+source_position_x-1,boxsize)+1
          j_global = modulo(j_source+source_position_y-2,boxsize)+1
          k_global = modulo(k_source+source_position_z-2,boxsize)+1
          
          pyramid_source_nx_py_pz_HI_density_array(i_source,j_source,k_source) = &
                        pyramid_global_number_density_array(i_global,j_global,k_global) * &
                        pyramid_global_xHI_array(i_global,j_global,k_global)
          pyramid_source_nx_py_pz_HeI_density_array(i_source,j_source,k_source) = &
                        pyramid_global_number_density_array(i_global,j_global,k_global) * &
                        pyramid_global_xHeI_array(i_global,j_global,k_global)
          pyramid_source_nx_py_pz_HeII_density_array(i_source,j_source,k_source) = &
                        pyramid_global_number_density_array(i_global,j_global,k_global) * &
                        pyramid_global_xHeII_array(i_global,j_global,k_global)
				  
        enddo
      enddo
    enddo

  end subroutine global_to_source_neg_x_pos_y_pos_z_transformation

  subroutine global_to_source_neg_x_pos_y_neg_z_transformation()	
	    
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
    boxsize = 2*level_py
			
    do i_source = -1,-level_py,-1
      do j_source = 1,level_py
        do k_source = -1,-level_py,-1
			 
          i_global = modulo(i_source+source_position_x-1,boxsize)+1
          j_global = modulo(j_source+source_position_y-2,boxsize)+1
          k_global = modulo(k_source+source_position_z-1,boxsize)+1
          
          pyramid_source_nx_py_nz_HI_density_array(i_source,j_source,k_source) = &
                        pyramid_global_number_density_array(i_global,j_global,k_global) * &
                        pyramid_global_xHI_array(i_global,j_global,k_global)
          pyramid_source_nx_py_nz_HeI_density_array(i_source,j_source,k_source) = &
                        pyramid_global_number_density_array(i_global,j_global,k_global) * &
                        pyramid_global_xHeI_array(i_global,j_global,k_global)
          pyramid_source_nx_py_nz_HeII_density_array(i_source,j_source,k_source) = &
                        pyramid_global_number_density_array(i_global,j_global,k_global) * &
                        pyramid_global_xHeII_array(i_global,j_global,k_global)
				  
        enddo
      enddo
    enddo

  end subroutine global_to_source_neg_x_pos_y_neg_z_transformation
   
  subroutine global_to_source_neg_x_neg_y_pos_z_transformation()	
	    
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
    boxsize = 2*level_py
			
    do i_source = -1,-level_py,-1
      do j_source = -1,-level_py,-1
        do k_source = 1,level_py
			 
          i_global = modulo(i_source+source_position_x-1,boxsize)+1
          j_global = modulo(j_source+source_position_y-1,boxsize)+1
          k_global = modulo(k_source+source_position_z-2,boxsize)+1
          
          pyramid_source_nx_ny_pz_HI_density_array(i_source,j_source,k_source) = &
                        pyramid_global_number_density_array(i_global,j_global,k_global) * &
                        pyramid_global_xHI_array(i_global,j_global,k_global)
          pyramid_source_nx_ny_pz_HeI_density_array(i_source,j_source,k_source) = &
                        pyramid_global_number_density_array(i_global,j_global,k_global) * &
                        pyramid_global_xHeI_array(i_global,j_global,k_global)
          pyramid_source_nx_ny_pz_HeII_density_array(i_source,j_source,k_source) = &
                        pyramid_global_number_density_array(i_global,j_global,k_global) * &
                        pyramid_global_xHeII_array(i_global,j_global,k_global)
				  
        enddo
      enddo
    enddo

  end subroutine global_to_source_neg_x_neg_y_pos_z_transformation
       
  subroutine global_to_source_neg_x_neg_y_neg_z_transformation()	
	    
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
    boxsize = 2*level_py
	
    do i_source = -1,-level_py,-1
      do j_source = -1,-level_py,-1
        do k_source = -1,-level_py,-1
			 
          i_global = modulo(i_source+source_position_x-1,boxsize)+1
          j_global = modulo(j_source+source_position_y-1,boxsize)+1
          k_global = modulo(k_source+source_position_z-1,boxsize)+1
          
          pyramid_source_nx_ny_nz_HI_density_array(i_source,j_source,k_source) = &
                        pyramid_global_number_density_array(i_global,j_global,k_global) * &
                        pyramid_global_xHI_array(i_global,j_global,k_global)
          pyramid_source_nx_ny_nz_HeI_density_array(i_source,j_source,k_source) = &
                        pyramid_global_number_density_array(i_global,j_global,k_global) * &
                        pyramid_global_xHeI_array(i_global,j_global,k_global)
          pyramid_source_nx_ny_nz_HeII_density_array(i_source,j_source,k_source) = &
                        pyramid_global_number_density_array(i_global,j_global,k_global) * &
                        pyramid_global_xHeII_array(i_global,j_global,k_global)
				  
        enddo
      enddo
    enddo
	
  end subroutine global_to_source_neg_x_neg_y_neg_z_transformation
  
  !!!!!!!!!!!!!!!
  
  subroutine source_to_global_pos_x_pos_y_pos_z_transformation()	
	    
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
    boxsize = 2*level_py
			
    do i_source = 1,level_py
      do j_source = 1,level_py
        do k_source = 1,level_py
			 
          i_global = modulo(i_source+source_position_x-2,boxsize)+1
          j_global = modulo(j_source+source_position_y-2,boxsize)+1
          k_global = modulo(k_source+source_position_z-2,boxsize)+1
          
          pyramid_global_HI_photoionization_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_px_py_pz_HI_photoionization_rate_array(i_source,j_source,k_source)
          pyramid_global_HeI_photoionization_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_px_py_pz_HeI_photoionization_rate_array(i_source,j_source,k_source)
          pyramid_global_HeII_photoionization_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_px_py_pz_HeII_photoionization_rate_array(i_source,j_source,k_source)							  					
          pyramid_global_photoheating_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_px_py_pz_photoheating_rate_array(i_source,j_source,k_source)	
									  
        enddo
      enddo
    enddo

  end subroutine source_to_global_pos_x_pos_y_pos_z_transformation

  subroutine source_to_global_pos_x_pos_y_neg_z_transformation()	
	    
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
    boxsize = 2*level_py
			
    do i_source = 1,level_py
      do j_source = 1,level_py
        do k_source = -1,-level_py,-1
			 
          i_global = modulo(i_source+source_position_x-2,boxsize)+1
          j_global = modulo(j_source+source_position_y-2,boxsize)+1
          k_global = modulo(k_source+source_position_z-1,boxsize)+1
          
          pyramid_global_HI_photoionization_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_px_py_nz_HI_photoionization_rate_array(i_source,j_source,k_source)	
          pyramid_global_HeI_photoionization_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_px_py_nz_HeI_photoionization_rate_array(i_source,j_source,k_source)
          pyramid_global_HeII_photoionization_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_px_py_nz_HeII_photoionization_rate_array(i_source,j_source,k_source)						  
          pyramid_global_photoheating_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_px_py_nz_photoheating_rate_array(i_source,j_source,k_source)
																	  				  
        enddo
      enddo
    enddo

  end subroutine source_to_global_pos_x_pos_y_neg_z_transformation
   
  subroutine source_to_global_pos_x_neg_y_pos_z_transformation()	
	    
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
    boxsize = 2*level_py
			
    do i_source = 1,level_py
      do j_source = -1,-level_py,-1
        do k_source = 1,level_py
			 
          i_global = modulo(i_source+source_position_x-2,boxsize)+1
          j_global = modulo(j_source+source_position_y-1,boxsize)+1
          k_global = modulo(k_source+source_position_z-2,boxsize)+1
          
          pyramid_global_HI_photoionization_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_px_ny_pz_HI_photoionization_rate_array(i_source,j_source,k_source)
          pyramid_global_HeI_photoionization_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_px_ny_pz_HeI_photoionization_rate_array(i_source,j_source,k_source)
          pyramid_global_HeII_photoionization_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_px_ny_pz_HeII_photoionization_rate_array(i_source,j_source,k_source)					
          pyramid_global_photoheating_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_px_ny_pz_photoheating_rate_array(i_source,j_source,k_source)
							  		  		  
        enddo
      enddo
    enddo

  end subroutine source_to_global_pos_x_neg_y_pos_z_transformation

  subroutine source_to_global_pos_x_neg_y_neg_z_transformation()	
	    
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
    boxsize = 2*level_py
			
    do i_source = 1,level_py
      do j_source = -1,-level_py,-1
        do k_source = -1,-level_py,-1
			 
          i_global = modulo(i_source+source_position_x-2,boxsize)+1
          j_global = modulo(j_source+source_position_y-1,boxsize)+1
          k_global = modulo(k_source+source_position_z-1,boxsize)+1
                   
          pyramid_global_HI_photoionization_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_px_ny_nz_HI_photoionization_rate_array(i_source,j_source,k_source)
          pyramid_global_HeI_photoionization_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_px_ny_nz_HeI_photoionization_rate_array(i_source,j_source,k_source)
          pyramid_global_HeII_photoionization_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_px_ny_nz_HeII_photoionization_rate_array(i_source,j_source,k_source)						
          pyramid_global_photoheating_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_px_ny_nz_photoheating_rate_array(i_source,j_source,k_source)
												  			  			  
        enddo
      enddo
    enddo

  end subroutine source_to_global_pos_x_neg_y_neg_z_transformation

  subroutine source_to_global_neg_x_pos_y_pos_z_transformation()	
	    
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
    boxsize = 2*level_py
			
    do i_source = -1,-level_py,-1
      do j_source = 1,level_py
        do k_source = 1,level_py
			 
          i_global = modulo(i_source+source_position_x-1,boxsize)+1
          j_global = modulo(j_source+source_position_y-2,boxsize)+1
          k_global = modulo(k_source+source_position_z-2,boxsize)+1
                    
          pyramid_global_HI_photoionization_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_nx_py_pz_HI_photoionization_rate_array(i_source,j_source,k_source)
          pyramid_global_HeI_photoionization_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_nx_py_pz_HeI_photoionization_rate_array(i_source,j_source,k_source)
          pyramid_global_HeII_photoionization_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_nx_py_pz_HeII_photoionization_rate_array(i_source,j_source,k_source)					
          pyramid_global_photoheating_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_nx_py_pz_photoheating_rate_array(i_source,j_source,k_source)	
													  	  		  
        enddo
      enddo
    enddo

  end subroutine source_to_global_neg_x_pos_y_pos_z_transformation

  subroutine source_to_global_neg_x_pos_y_neg_z_transformation()	
	    
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
    boxsize = 2*level_py
			
    do i_source = -1,-level_py,-1
      do j_source = 1,level_py
        do k_source = -1,-level_py,-1
			 
          i_global = modulo(i_source+source_position_x-1,boxsize)+1
          j_global = modulo(j_source+source_position_y-2,boxsize)+1
          k_global = modulo(k_source+source_position_z-1,boxsize)+1
                   
          pyramid_global_HI_photoionization_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_nx_py_nz_HI_photoionization_rate_array(i_source,j_source,k_source)
          pyramid_global_HeI_photoionization_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_nx_py_nz_HeI_photoionization_rate_array(i_source,j_source,k_source)
          pyramid_global_HeII_photoionization_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_nx_py_nz_HeII_photoionization_rate_array(i_source,j_source,k_source)						
          pyramid_global_photoheating_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_nx_py_nz_photoheating_rate_array(i_source,j_source,k_source)	
								  		  		  
        enddo
      enddo
    enddo

  end subroutine source_to_global_neg_x_pos_y_neg_z_transformation
   
  subroutine source_to_global_neg_x_neg_y_pos_z_transformation()	
	    
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
    boxsize = 2*level_py
			
    do i_source = -1,-level_py,-1
      do j_source = -1,-level_py,-1
        do k_source = 1,level_py
			 
          i_global = modulo(i_source+source_position_x-1,boxsize)+1
          j_global = modulo(j_source+source_position_y-1,boxsize)+1
          k_global = modulo(k_source+source_position_z-2,boxsize)+1
                    
          pyramid_global_HI_photoionization_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_nx_ny_pz_HI_photoionization_rate_array(i_source,j_source,k_source)
          pyramid_global_HeI_photoionization_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_nx_ny_pz_HeI_photoionization_rate_array(i_source,j_source,k_source)
          pyramid_global_HeII_photoionization_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_nx_ny_pz_HeII_photoionization_rate_array(i_source,j_source,k_source)							
          pyramid_global_photoheating_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_nx_ny_pz_photoheating_rate_array(i_source,j_source,k_source)
							  		  		  
        enddo
      enddo
    enddo

  end subroutine source_to_global_neg_x_neg_y_pos_z_transformation
       
  subroutine source_to_global_neg_x_neg_y_neg_z_transformation()	
	    
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
    boxsize = 2*level_py
	
    do i_source = -1,-level_py,-1
      do j_source = -1,-level_py,-1
        do k_source = -1,-level_py,-1
			 
          i_global = modulo(i_source+source_position_x-1,boxsize)+1
          j_global = modulo(j_source+source_position_y-1,boxsize)+1
          k_global = modulo(k_source+source_position_z-1,boxsize)+1
                    
          pyramid_global_HI_photoionization_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_nx_ny_nz_HI_photoionization_rate_array(i_source,j_source,k_source)
          pyramid_global_HeI_photoionization_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_nx_ny_nz_HeI_photoionization_rate_array(i_source,j_source,k_source)
          pyramid_global_HeII_photoionization_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_nx_ny_nz_HeII_photoionization_rate_array(i_source,j_source,k_source)							
          pyramid_global_photoheating_rate_array(i_global,j_global,k_global) = &
                    pyramid_source_nx_ny_nz_photoheating_rate_array(i_source,j_source,k_source)
							  			  		  
        enddo
      enddo
    enddo
	
  end subroutine source_to_global_neg_x_neg_y_neg_z_transformation
  
end module pyramid_global_source_transformation
