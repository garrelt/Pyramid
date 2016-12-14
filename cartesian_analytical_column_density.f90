module cartesian_analytical_column_density
	
  use precision, only: dp
  use input, only: level_py, cellsize_sl, four_over_three_pi
  use array, only: cartesian_analytical_number_density, &
                   cartesian_analytical_xHI, &
                   cartesian_analytical_xHeI  , &
                   cartesian_analytical_xHeII, &
				   cartesian_analytical_HI_column_density_in, cartesian_analytical_HeI_column_density_in, &
				   cartesian_analytical_HeII_column_density_in, cartesian_analytical_HI_column_density_out, &
				   cartesian_analytical_HeI_column_density_out, cartesian_analytical_HeII_column_density_out, &
				   cartesian_analytical_shell_volume 
				   
  contains
	
  subroutine cartesian_analytical_column_density_calculation
		
    implicit none

  	real(kind=dp) :: absolute_x, absolute_y, absolute_z
	real(kind=dp) :: in_x, in_y, in_z
	real(kind=dp) :: out_x, out_y, out_z
  	real(kind=dp) :: in_distance, out_distance
	integer :: i,j,k
		
	cartesian_analytical_HI_column_density_in(0,0,0) = 0 
	cartesian_analytical_HeI_column_density_in(0,0,0) = 0 
	cartesian_analytical_HeII_column_density_in(0,0,0) = 0 
	cartesian_analytical_shell_volume(0,0,0) = cellsize_sl*cellsize_sl*cellsize_sl
cartesian_analytical_HI_column_density_out(0,0,0) = 0.5*cellsize_sl*cartesian_analytical_number_density*cartesian_analytical_xHI
cartesian_analytical_HeI_column_density_out(0,0,0) = 0.5*cellsize_sl*cartesian_analytical_number_density*cartesian_analytical_xHeI
cartesian_analytical_HeII_column_density_out(0,0,0) = 0.5*cellsize_sl*cartesian_analytical_number_density*cartesian_analytical_xHeII

    do i = -level_py,level_py
      do j = -level_py,level_py
        do k = -level_py,level_py
if (i.ne.0 .or. j.ne.0 .or. k.ne.0) then
absolute_x = abs(i)
absolute_y = abs(j)
absolute_z = abs(k)

if (absolute_z.ge.absolute_x .and. absolute_z.ge.absolute_y) then
in_x = (absolute_x)*(absolute_z-0.5)/(absolute_z)
in_y = (absolute_y)*(absolute_z-0.5)/(absolute_z)
in_z = absolute_z-0.5	
out_x = (absolute_x)*(absolute_z+0.5)/(absolute_z)
out_y = (absolute_y)*(absolute_z+0.5)/(absolute_z)				
out_z = absolute_z+0.5
elseif (absolute_y.ge.absolute_x .and. absolute_y.ge.absolute_z) then
in_x = (absolute_x)*(absolute_y-0.5)/(absolute_y)
in_y = absolute_y-0.5	
in_z = (absolute_z)*(absolute_y-0.5)/(absolute_y)				
out_x = (absolute_x)*(absolute_y+0.5)/(absolute_y)
out_y = absolute_y+0.5
out_z = (absolute_z)*(absolute_y+0.5)/(absolute_y)
elseif (absolute_x.ge.absolute_y .and. absolute_x.ge.absolute_z) then
in_x = absolute_x-0.5
in_y = (absolute_y)*(absolute_x-0.5)/(absolute_x)
in_z = (absolute_z)*(absolute_x-0.5)/(absolute_x)	
out_x = absolute_x+0.5
out_y = (absolute_y)*(absolute_x+0.5)/(absolute_x)
out_z = (absolute_z)*(absolute_x+0.5)/(absolute_x)				  			  
endif

in_distance = in_x*in_x + in_y*in_y + in_z*in_z
in_distance = in_distance**0.5
in_distance = in_distance*cellsize_sl
out_distance = out_x*out_x + out_y*out_y + out_z*out_z
out_distance = out_distance**0.5
out_distance = out_distance*cellsize_sl
cartesian_analytical_HI_column_density_in(i,j,k) = in_distance*cartesian_analytical_number_density*cartesian_analytical_xHI
cartesian_analytical_HeI_column_density_in(i,j,k) = in_distance*cartesian_analytical_number_density*cartesian_analytical_xHeI
cartesian_analytical_HeII_column_density_in(i,j,k) = in_distance*cartesian_analytical_number_density*cartesian_analytical_xHeII
cartesian_analytical_HI_column_density_out(i,j,k) = out_distance*cartesian_analytical_number_density*cartesian_analytical_xHI
cartesian_analytical_HeI_column_density_out(i,j,k) = out_distance*cartesian_analytical_number_density*cartesian_analytical_xHeI
cartesian_analytical_HeII_column_density_out(i,j,k) = out_distance*cartesian_analytical_number_density*cartesian_analytical_xHeII
cartesian_analytical_shell_volume(i,j,k) = four_over_three_pi*(out_distance*out_distance*out_distance-&
in_distance*in_distance*in_distance)

endif

        enddo
      enddo
    enddo	
	
!write(*,*) 'c in ',cartesian_analytical_HI_column_density_in(1,0,0)
!write(*,*) 'c out ',cartesian_analytical_HI_column_density_out(1,0,0)
	
  end subroutine cartesian_analytical_column_density_calculation	
	
end module cartesian_analytical_column_density
