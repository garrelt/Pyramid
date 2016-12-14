module pyramid_analytical_column_density
	
  use precision, only: dp
  use input, only: level_py, cellsize_py, four_over_three_pi
  use array, only: pyramid_analytical_xHI, pyramid_analytical_xHeI, &
                   pyramid_analytical_xHeII, pyramid_analytical_number_density, &
				   pyramid_analytical_HI_column_density_in, pyramid_analytical_HeI_column_density_in, &
				   pyramid_analytical_HeII_column_density_in, pyramid_analytical_HI_column_density_out, &
				   pyramid_analytical_HeI_column_density_out, pyramid_analytical_HeII_column_density_out, &
				   pyramid_analytical_shell_volume 
				   
  contains
	
  subroutine pyramid_analytical_column_density_calculation
		
    implicit none

  	real(kind=dp) :: absolute_x, absolute_y, absolute_z
	real(kind=dp) :: in_x, in_y, in_z
	real(kind=dp) :: out_x, out_y, out_z
  	real(kind=dp) :: in_distance, out_distance
	integer :: i,j,k
		
    do i = 1,2*level_py
      do j = 1,2*level_py
        do k = 1,2*level_py
		  if (i.le.level_py) absolute_x = real(level_py-i+1)
		  if (i.ge.level_py+1) absolute_x = real(i-level_py)
		  if (j.le.level_py) absolute_y = real(level_py-j+1)
		  if (j.ge.level_py+1) absolute_y = real(j-level_py)		  
		  if (k.le.level_py) absolute_z = real(level_py-k+1)
		  if (k.ge.level_py+1) absolute_z = real(k-level_py)		  
		  if (absolute_z.ge.absolute_x .and. absolute_z.ge.absolute_y) then
                    in_x = (absolute_x-0.5)*(absolute_z-1.0)/(absolute_z-0.5)
                    in_y = (absolute_y-0.5)*(absolute_z-1.0)/(absolute_z-0.5)
                    in_z = absolute_z-1.0	
                    out_x = (absolute_x-0.5)*(absolute_z)/(absolute_z-0.5)
                    out_y = (absolute_y-0.5)*(absolute_z)/(absolute_z-0.5)						
                    out_z = absolute_z
                  elseif (absolute_y.ge.absolute_x .and. absolute_y.ge.absolute_z) then
                    in_x = (absolute_x-0.5)*(absolute_y-1.0)/(absolute_y-0.5)
                    in_y = absolute_y-1.0	
                    in_z = (absolute_z-0.5)*(absolute_y-1.0)/(absolute_y-0.5)				
                    out_x = (absolute_x-0.5)*(absolute_y)/(absolute_y-0.5)
                    out_y = absolute_y
                    out_z = (absolute_z-0.5)*(absolute_y)/(absolute_y-0.5)							  
		  elseif (absolute_x.ge.absolute_y .and. absolute_x.ge.absolute_z) then
                    in_x = absolute_x-1.0
                    in_y = (absolute_y-0.5)*(absolute_x-1.0)/(absolute_x-0.5)
                    in_z = (absolute_z-0.5)*(absolute_x-1.0)/(absolute_x-0.5)	
                    out_x = absolute_x
                    out_y = (absolute_y-0.5)*(absolute_x)/(absolute_x-0.5)
                    out_z = (absolute_z-0.5)*(absolute_x)/(absolute_x-0.5)				  			  
		  endif
		  
		  if (absolute_x.eq.1 .and. absolute_y.eq.1 .and. absolute_z.eq.1) then
                    pyramid_analytical_HI_column_density_in(i,j,k) = 0.0
                    pyramid_analytical_HeI_column_density_in(i,j,k) = 0.0
                    pyramid_analytical_HeII_column_density_in(i,j,k) = 0.0
pyramid_analytical_HI_column_density_out(i,j,k) = cellsize_py*pyramid_analytical_number_density*pyramid_analytical_xHI
pyramid_analytical_HeI_column_density_out(i,j,k) = cellsize_py*pyramid_analytical_number_density*pyramid_analytical_xHeI
pyramid_analytical_HeII_column_density_out(i,j,k) = cellsize_py*pyramid_analytical_number_density*pyramid_analytical_xHeII
pyramid_analytical_shell_volume(i,j,k) = 8.0*cellsize_py*cellsize_py*cellsize_py
          else									 
		    in_distance = in_x*in_x + in_y*in_y + in_z*in_z
		    in_distance = in_distance**0.5
		    in_distance = in_distance*cellsize_py
		    out_distance = out_x*out_x + out_y*out_y + out_z*out_z
		    out_distance = out_distance**0.5
		    out_distance = out_distance*cellsize_py
pyramid_analytical_HI_column_density_in(i,j,k) = in_distance*pyramid_analytical_number_density*pyramid_analytical_xHI
pyramid_analytical_HeI_column_density_in(i,j,k) = in_distance*pyramid_analytical_number_density*pyramid_analytical_xHeI
pyramid_analytical_HeII_column_density_in(i,j,k) = in_distance*pyramid_analytical_number_density*pyramid_analytical_xHeII
pyramid_analytical_HI_column_density_out(i,j,k) = out_distance*pyramid_analytical_number_density*pyramid_analytical_xHI
pyramid_analytical_HeI_column_density_out(i,j,k) = out_distance*pyramid_analytical_number_density*pyramid_analytical_xHeI
pyramid_analytical_HeII_column_density_out(i,j,k) = out_distance*pyramid_analytical_number_density*pyramid_analytical_xHeII
pyramid_analytical_shell_volume(i,j,k) = four_over_three_pi*(out_distance*out_distance*out_distance-&
in_distance*in_distance*in_distance)
          endif
        enddo
      enddo
    enddo	
		
  end subroutine pyramid_analytical_column_density_calculation	
	
end module pyramid_analytical_column_density
