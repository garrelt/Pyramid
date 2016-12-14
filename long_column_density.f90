module long_column_density
	
  use precision, only: dp
  use input, only: level_py, cellsize_sl
  use array, only: long_source_HI_density_array, long_source_HeI_density_array, long_source_HeII_density_array, &
                   long_source_HI_column_density_in, long_source_HeI_column_density_in, long_source_HeII_column_density_in, &
                   long_source_HI_column_density_out, long_source_HeI_column_density_out, long_source_HeII_column_density_out, &
                   long_source_shell_volume
  use long, only: long_characteristic
  
contains	
	
  subroutine long_column_density_calculation()
	  
    implicit none
	
    integer :: i,j,k	
	
	long_source_HI_column_density_in(0,0,0) = 0 
	long_source_HeI_column_density_in(0,0,0) = 0 
	long_source_HeII_column_density_in(0,0,0) = 0 
	long_source_shell_volume(0,0,0) = cellsize_sl*cellsize_sl*cellsize_sl
	long_source_HI_column_density_out(0,0,0) = 0.5*cellsize_sl*long_source_HI_density_array(0,0,0)
	long_source_HeI_column_density_out(0,0,0) = 0.5*cellsize_sl*long_source_HeI_density_array(0,0,0)
	long_source_HeII_column_density_out(0,0,0) = 0.5*cellsize_sl*long_source_HeII_density_array(0,0,0)

    do i = -level_py,level_py
      do j = -level_py,level_py
        do k = -level_py,level_py
          if (i.ne.0 .or. j.ne.0 .or. k.ne.0 ) call long_characteristic(i,j,k)
        enddo
      enddo
    enddo

  end subroutine long_column_density_calculation
  
end module long_column_density
