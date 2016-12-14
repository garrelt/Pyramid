module pyramid_analytical_photo_rate
	
  use precision, only: dp
  use input, only: level_py, cellsize_py, epsilon, use_which_source
  use array, only: pyramid_analytical_xHI, pyramid_analytical_xHeI, &
                   pyramid_analytical_xHeII, pyramid_analytical_number_density, &
				   pyramid_analytical_HI_column_density_in, &
                   pyramid_analytical_HeI_column_density_in, pyramid_analytical_HeII_column_density_in, &
                   pyramid_analytical_HI_column_density_out, pyramid_analytical_HeI_column_density_out, &
                   pyramid_analytical_HeII_column_density_out, pyramid_analytical_shell_volume, &
                   pyramid_analytical_HI_photoionization_rate_array, pyramid_analytical_HeI_photoionization_rate_array, &				   
                   pyramid_analytical_HeII_photoionization_rate_array, pyramid_analytical_photoheating_rate_array  
  use radiation, only: photoion_shell, photrates
    
contains
	
  subroutine pyramid_analytical_photo_rate_calculation()
		
    implicit none
  		
    type(photrates) :: phi	
    integer :: i,j,k
		
    do i = 1,2*level_py
      do j = 1,2*level_py
        do k = 1,2*level_py
				
call photoion_shell (phi, &
pyramid_analytical_HI_column_density_in(i,j,k),pyramid_analytical_HI_column_density_out(i,j,k),&
pyramid_analytical_HeI_column_density_in(i,j,k),pyramid_analytical_HeI_column_density_out(i,j,k),&
pyramid_analytical_HeII_column_density_in(i,j,k),pyramid_analytical_HeII_column_density_out(i,j,k),&
pyramid_analytical_shell_volume(i,j,k),use_which_source)
								
pyramid_analytical_HI_photoionization_rate_array(i,j,k) = phi%photo_cell_HI/ &
(pyramid_analytical_xHI*pyramid_analytical_number_density)
pyramid_analytical_HeI_photoionization_rate_array(i,j,k) = phi%photo_cell_HeI/ &
(pyramid_analytical_xHeI*pyramid_analytical_number_density)		  
pyramid_analytical_HeII_photoionization_rate_array(i,j,k) = phi%photo_cell_HeII/ &
(pyramid_analytical_xHeII*pyramid_analytical_number_density)		  
pyramid_analytical_photoheating_rate_array(i,j,k) = phi%heat	 	
		  							  
        enddo
      enddo
    enddo
		
  end subroutine pyramid_analytical_photo_rate_calculation 	
	
end module pyramid_analytical_photo_rate
