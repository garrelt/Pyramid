module cartesian_analytical_photo_rate
	
  use precision, only: dp
  use input, only: level_py, cellsize_py, epsilon, use_which_source
  use array, only: cartesian_analytical_xHI, cartesian_analytical_xHeI, &
                   cartesian_analytical_xHeII, cartesian_analytical_number_density, &
                   cartesian_analytical_HI_column_density_in, &
                   cartesian_analytical_HeI_column_density_in, cartesian_analytical_HeII_column_density_in, &
                   cartesian_analytical_HI_column_density_out, cartesian_analytical_HeI_column_density_out, &
                   cartesian_analytical_HeII_column_density_out, cartesian_analytical_shell_volume, &
                   cartesian_analytical_HI_photoionization_rate_array, &
                   cartesian_analytical_HeI_photoionization_rate_array, &				   
                   cartesian_analytical_HeII_photoionization_rate_array, &
                   cartesian_analytical_photoheating_rate_array  
  use radiation, only: photoion_shell, photrates
    
contains
	
  subroutine cartesian_analytical_photo_rate_calculation()
		
    implicit none
  		
    type(photrates) :: phi	
    integer :: i,j,k
		
    do i = -level_py,level_py
      do j = -level_py,level_py
        do k = -level_py,level_py
				
call photoion_shell (phi, &
cartesian_analytical_HI_column_density_in(i,j,k),cartesian_analytical_HI_column_density_out(i,j,k),&
cartesian_analytical_HeI_column_density_in(i,j,k),cartesian_analytical_HeI_column_density_out(i,j,k),&
cartesian_analytical_HeII_column_density_in(i,j,k),cartesian_analytical_HeII_column_density_out(i,j,k),&
cartesian_analytical_shell_volume(i,j,k),use_which_source)
								
cartesian_analytical_HI_photoionization_rate_array(i,j,k) = phi%photo_cell_HI/ &
(cartesian_analytical_xHI*cartesian_analytical_number_density)
cartesian_analytical_HeI_photoionization_rate_array(i,j,k) = phi%photo_cell_HeI/ &
(cartesian_analytical_xHeI*cartesian_analytical_number_density)		  
cartesian_analytical_HeII_photoionization_rate_array(i,j,k) = phi%photo_cell_HeII/ &
(cartesian_analytical_xHeII*cartesian_analytical_number_density)		  
cartesian_analytical_photoheating_rate_array(i,j,k) = phi%heat	 	
		  							  
        enddo
      enddo
    enddo
	
  end subroutine cartesian_analytical_photo_rate_calculation 	
	
end module cartesian_analytical_photo_rate
