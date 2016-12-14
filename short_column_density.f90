module short_column_density
	
  use precision, only: dp
  use input, only: level_py, cellsize_sl
  use array, only: short_source_HI_density_array, short_source_HeI_density_array, short_source_HeII_density_array, &
                   short_source_HI_column_density_in, short_source_HeI_column_density_in, short_source_HeII_column_density_in, &
                   short_source_HI_column_density_out, short_source_HeI_column_density_out, short_source_HeII_column_density_out, &
                   short_source_shell_volume
  use short, only: short_characteristic
  
contains	
	
  subroutine short_column_density_calculation()
	
    !0D - source position
    short_source_HI_column_density_in(0,0,0) = 0 
    short_source_HeI_column_density_in(0,0,0) = 0 
    short_source_HeII_column_density_in(0,0,0) = 0 
    short_source_shell_volume(0,0,0) = cellsize_sl*cellsize_sl*cellsize_sl
    short_source_HI_column_density_out(0,0,0) = 0.5*cellsize_sl*short_source_HI_density_array(0,0,0)
    short_source_HeI_column_density_out(0,0,0) = 0.5*cellsize_sl*short_source_HeI_density_array(0,0,0)
    short_source_HeII_column_density_out(0,0,0) = 0.5*cellsize_sl*short_source_HeII_density_array(0,0,0)
		  
    !1D +i direction	  
    do i = 1,level_py
      call short_characteristic(i,0,0)
    enddo
	
    !1D -i direction	  
    do i = -1,-level_py,-1
      call short_characteristic(i,0,0)
    enddo

    !1D +j direction	  
    do j = 1,level_py
      call short_characteristic(0,j,0)
    enddo
	
    !1D -j direction	  
    do j = -1,-level_py,-1
      call short_characteristic(0,j,0)
    enddo
	
    !1D +k direction	  
    do k = 1,level_py
      call short_characteristic(0,0,k)
    enddo
	
    !1D -k direction	  
    do k = -1,-level_py,-1
      call short_characteristic(0,0,k)
    enddo
				
    !2D +i+j plane
    do i = 1,level_py
      do j = 1,level_py
        call short_characteristic(i,j,0)						
      enddo
    enddo
	
    !2D +i-j plane
    do i = 1,level_py
      do j = -1,-level_py,-1
        call short_characteristic(i,j,0)						
      enddo
    enddo	

    !2D -i+j plane
    do i = -1,-level_py,-1
      do j = 1,level_py
        call short_characteristic(i,j,0)						
      enddo
    enddo

    !2D -i-j plane
    do i = -1,-level_py,-1
      do j = -1,-level_py,-1
        call short_characteristic(i,j,0)						
      enddo
    enddo

    !2D +i+k plane
    do i = 1,level_py
      do k = 1,level_py
        call short_characteristic(i,0,k)						
      enddo
    enddo
	
    !2D +i-k plane
    do i = 1,level_py
      do k = -1,-level_py,-1
        call short_characteristic(i,0,k)						
      enddo
    enddo	

    !2D -i+k plane
    do i = -1,-level_py,-1
      do k = 1,level_py
        call short_characteristic(i,0,k)						
      enddo
    enddo

    !2D -i-k plane
    do i = -1,-level_py,-1
      do k = -1,-level_py,-1
        call short_characteristic(i,0,k)						
      enddo
    enddo
	
    !2D +j+k plane
    do j = 1,level_py
      do k = 1,level_py
        call short_characteristic(0,j,k)						
      enddo
    enddo
	
    !2D +j-k plane
    do j = 1,level_py
      do k = -1,-level_py,-1
        call short_characteristic(0,j,k)						
      enddo
    enddo	

    !2D -j+k plane
    do j = -1,-level_py,-1
      do k = 1,level_py
        call short_characteristic(0,j,k)						
      enddo
    enddo

    !2D -j-k plane
    do j = -1,-level_py,-1
      do k = -1,-level_py,-1
        call short_characteristic(0,j,k)						
      enddo
    enddo
				
    !3D +i+j+k octant
    do i = 1,level_py
      do j = 1,level_py
        do k = 1,level_py
          call short_characteristic(i,j,k)			
        enddo
      enddo			
    enddo			
			    
    !3D +i+j-k octant
    do i = 1,level_py
      do j = 1,level_py
        do k = -1,-level_py,-1
          call short_characteristic(i,j,k)			
        enddo
      enddo			
    enddo

    !3D +i-j+k octant
    do i = 1,level_py
      do j = -1,-level_py,-1
        do k = 1,level_py
          call short_characteristic(i,j,k)			
        enddo
      enddo			
    enddo			
			    
    !3D +i-j-k octant
    do i = 1,level_py
      do j = -1,-level_py,-1
        do k = -1,-level_py,-1
          call short_characteristic(i,j,k)			
        enddo
      enddo			
    enddo

    !3D -i+j+k octant
    do i = -1,-level_py,-1
      do j = 1,level_py
        do k = 1,level_py
          call short_characteristic(i,j,k)			
        enddo
      enddo			
    enddo			
			    
    !3D -i+j-k octant
    do i = -1,-level_py,-1
      do j = 1,level_py
        do k = -1,-level_py,-1
          call short_characteristic(i,j,k)			
        enddo
      enddo			
    enddo

    !3D -i-j+k octant
    do i = -1,-level_py,-1
      do j = -1,-level_py,-1
        do k = 1,level_py
          call short_characteristic(i,j,k)			
        enddo
      enddo			
    enddo			
			    
    !3D -i-j-k octant
    do i = -1,-level_py,-1
      do j = -1,-level_py,-1
        do k = -1,-level_py,-1
          call short_characteristic(i,j,k)			
        enddo
      enddo			
    enddo
				
  end subroutine short_column_density_calculation
  
end module short_column_density
