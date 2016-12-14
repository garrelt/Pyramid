module analytical_photo_rate
	
  use precision, only: dp
  use input, only: level_py, cellsize_py, epsilon
  use array, only: analytical_xHI, analytical_xHeI, &
                   analytical_xHeII, analytical_number_density, &
				   analytical_HI_column_density_in, &
                   analytical_HeI_column_density_in, analytical_HeII_column_density_in, &
                   analytical_HI_column_density_out, analytical_HeI_column_density_out, &
                   analytical_HeII_column_density_out, analytical_shell_volume, &
                   analytical_HI_photoionization_rate_array, analytical_HeI_photoionization_rate_array, &				   
                   analytical_HeII_photoionization_rate_array, analytical_photoheating_rate_array  
  use radiation, only: photoion_shell, photrates
    
contains
	
  subroutine analytical_photo_rate_calculation()
		
    implicit none
  		
    type(photrates) :: phi	
    integer :: i,j,k
		
    do i = 1,2*level_py
      do j = 1,2*level_py
        do k = 1,2*level_py
				
			call photoion_shell (phi, &
			                     analytical_HI_column_density_in(i,j,k),analytical_HI_column_density_out(i,j,k),&
			                     analytical_HeI_column_density_in(i,j,k),analytical_HeI_column_density_out(i,j,k),&
			                     analytical_HeII_column_density_in(i,j,k),analytical_HeII_column_density_out(i,j,k),&
			                     analytical_shell_volume(i,j,k),'B')
								
			analytical_HI_photoionization_rate_array(i,j,k) = phi%photo_cell_HI/(analytical_xHI*analytical_number_density)
			analytical_HeI_photoionization_rate_array(i,j,k) = phi%photo_cell_HeI/(analytical_xHeI*analytical_number_density)		  
			analytical_HeII_photoionization_rate_array(i,j,k) = phi%photo_cell_HeII/(analytical_xHeII*analytical_number_density)		  
			analytical_photoheating_rate_array(i,j,k) = phi%heat	 	
		  							  
        enddo
      enddo
    enddo
		
  end subroutine analytical_photo_rate_calculation 	
	
end module analytical_photo_rate