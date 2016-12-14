module short_photo_rate
	
  use precision, only: dp
  use input, only: level_py, use_which_source
  use array, only: short_source_HI_density_array, short_source_HeI_density_array, short_source_HeII_density_array, &
                   short_source_HI_column_density_in, short_source_HeI_column_density_in, short_source_HeII_column_density_in, &
                   short_source_HI_column_density_out, short_source_HeI_column_density_out, short_source_HeII_column_density_out, &
                   short_source_HI_photoionization_rate_array, short_source_HeI_photoionization_rate_array, &			   
                   short_source_HeII_photoionization_rate_array, short_source_photoheating_rate_array, &				   
				   short_source_shell_volume
  use radiation, only: photoion_shell, photrates

contains
	
  subroutine short_photo_rate_calculation()
		
    type(photrates) :: phi	
    integer :: i,j,k
		
    do i = -level_py,level_py
      do j = -level_py,level_py
        do k = -level_py,level_py
				
          call photoion_shell (phi, &
                               short_source_HI_column_density_in(i,j,k),short_source_HI_column_density_out(i,j,k),&
                               short_source_HeI_column_density_in(i,j,k),short_source_HeI_column_density_out(i,j,k),&
                               short_source_HeII_column_density_in(i,j,k),short_source_HeII_column_density_out(i,j,k),&
                               short_source_shell_volume(i,j,k),use_which_source)
          short_source_HI_photoionization_rate_array(i,j,k) = phi%photo_cell_HI/short_source_HI_density_array(i,j,k)
          short_source_HeI_photoionization_rate_array(i,j,k) = phi%photo_cell_HeI/short_source_HeI_density_array(i,j,k)		  
          short_source_HeII_photoionization_rate_array(i,j,k) = phi%photo_cell_HeII/short_source_HeII_density_array(i,j,k)		  
          short_source_photoheating_rate_array(i,j,k) = phi%heat	 	 	
		  		  		  				   
        enddo
      enddo
    enddo
	
  end subroutine short_photo_rate_calculation
  	
end module short_photo_rate
