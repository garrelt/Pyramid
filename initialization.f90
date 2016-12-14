module array
	
  use precision, only: dp
  use input, only: level, partition, number_density, xHI, xHII, temperature, pi
  use type, only: pyramid, cartesian, data
  use solid_angle, only: ss, cc
  
  implicit none
	
  type(pyramid), dimension(1:level) :: pyramid_grid
  type(cartesian), dimension(1:level) :: cartesian_grid
  type(data), dimension(1:level) :: pos_x_pos_y_pos_z_dir_z
	  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pyramid_creation()

    ! allocation arrays of volume, length...
    integer :: count

    do count = 1,level
	
      allocate(pyramid_grid(count)%volume(1:partition(count),1:partition(count)))
      allocate(pyramid_grid(count)%cellsize(1:partition(count),1:partition(count)))
      allocate(pyramid_grid(count)%distance(1:partition(count),1:partition(count)))
      allocate(pyramid_grid(count)%solid_angle(1:partition(count),1:partition(count)))
      allocate(pyramid_grid(count)%case(1:partition(count),1:partition(count)))
      allocate(pyramid_grid(count)%transform(1:partition(count),1:partition(count)))
	  
    end do

 end subroutine pyramid_creation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cartesian_creation()

    ! allocation arrays of volume, length...
    integer :: count

    do count = 1,level
	  
      allocate(cartesian_grid(count)%transform(1:level,1:level))
	  
    end do

  end subroutine cartesian_creation
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  	  
  subroutine domain_creation()

    ! allocation arrays of number_density, xHI...
    integer :: count

    do count = 1,level
	
      allocate(pos_x_pos_y_pos_z_dir_z(count)%number_density(1:partition(count),1:partition(count)))
      allocate(pos_x_pos_y_pos_z_dir_z(count)%xHI(1:partition(count),1:partition(count)))
      allocate(pos_x_pos_y_pos_z_dir_z(count)%xHII(1:partition(count),1:partition(count)))
      allocate(pos_x_pos_y_pos_z_dir_z(count)%temperature(1:partition(count),1:partition(count)))

    end do

end subroutine domain_creation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pyramid_volume()

    integer :: k, n
    real :: volume

    ! Volume on the same level are the same
    do k = 1,level

      n = partition(k)
      volume = real(3*k*k-3*k+1)/real(3*n*n)
      pyramid_grid(k)%volume = volume
  
    end do

  end subroutine pyramid_volume

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pyramid_cellsize()

    integer :: i,j,k,n

    !The radial length of the pyramid_grid
    do k = 1,level
      do i = 1,partition(k)
        do j = 1,partition(k)

          n = partition(k)
          pyramid_grid(k)%cellsize(i,j) = sqrt(real(i*i+j*j+n*n)+0.5)/real(n)

        end do
      end do
    end do

 end subroutine pyramid_cellsize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pyramid_distance()

    integer :: i,j,k,n

    !Distance of the center of pyramid_grid to the source
    do k = 1,level
		
      n = partition(k)
      do i = 1,partition(k)
        do j = 1,partition(k)

          pyramid_grid(k)%distance(i,j) = sqrt(real(i*i+j*j+n*n-i-j+0.5))*(2*k-1)/real(2*n)

        end do
      end do
    end do

  end subroutine pyramid_distance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pyramid_solid_angle()

    integer :: i,j,k,n

    n = partition(level)
	
    ! solid angle spanned by the pyramid_grid of the toppest level
    do i = 1,partition(level)
      do j = 1,partition(level)
        if (i.ge.j) then
			
          pyramid_grid(level)%solid_angle(i,j) = ss(i-1,j,i-1,n) - ss(i-1,j-1,i-1,n) - &
                                                 ss(i,j,i,n) + ss(i,j-1,i,n) + &
                                                 cc(j,j,i,n) - cc(j,j,i-1,n) + &
                                                 cc(j-1,j-1,i-1,n) - cc(j-1,j-1,i,n)
					
        elseif (i.lt.j) then
			
          pyramid_grid(level)%solid_angle(i,j) = ss(j-1,i,j-1,n) - ss(j-1,i-1,j-1,n) - &
                                                 ss(j,i,j,n) + ss(j,i-1,j,n) + &
                                                 cc(i,i,j,n) - cc(i,i,j-1,n) + &
                                                 cc(i-1,i-1,j-1,n) - cc(i-1,i-1,j,n)
				
        endif	
      end do
    end do
	
    ! solid angle spanned by the pyramid_grid of the other level
    do k = level-1,1,-1
      if (partition(k).eq.partition(k+1)) then
		  
         pyramid_grid(k)%solid_angle = pyramid_grid(k+1)%solid_angle
		 
      else
		  
        do i = 1,partition(k)
          do j = 1,partition(k)
			  
            pyramid_grid(k)%solid_angle(i,j) = pyramid_grid(k+1)%solid_angle(2*i-1,2*j-1) + &
                                               pyramid_grid(k+1)%solid_angle(2*i-1,2*j) + &
                                               pyramid_grid(k+1)%solid_angle(2*i,2*j-1) + &
                                               pyramid_grid(k+1)%solid_angle(2*i,2*j)
				
          end do
        end do  
      endif
    enddo
  
  end subroutine pyramid_solid_angle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine domain_number_density()

	integer :: i,j,k,n

	!Distance of the center of pyramid_grid to the source
	do k = 1,level
	  do i = 1,partition(k)
	    do j = 1,partition(k)
		  
		  pos_x_pos_y_pos_z_dir_z(k)%number_density(i,j) = number_density
		  
	    end do
	  end do
	end do

  end subroutine domain_number_density

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine domain_xHI()

    integer :: i,j,k,n

    ! Ionization fraction of the pyramid_grids
    do k = 1,level
      do i = 1,partition(k)
        do j = 1,partition(k)

          pos_x_pos_y_pos_z_dir_z(k)%xHI(i,j) = xHI

      end do
    end do
  end do

end subroutine domain_xHI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine domain_xHII()

    integer :: i,j,k,n

    ! Ionization fraction of the pyramid_grids
    do k = 1,level
      do i = 1,partition(k)
        do j = 1,partition(k)

          pos_x_pos_y_pos_z_dir_z(k)%xHII(i,j) = xHII

        end do
     end do
    end do

  end subroutine domain_xHII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine domain_temperature()

    integer :: i,j,k,n

    ! Temperature of the pyramid_grids
    do k = 1,level
      do i = 1,partition(k)
        do j = 1,partition(k)

          pos_x_pos_y_pos_z_dir_z(k)%temperature(i,j) = temperature

        end do
      end do
    end do

  end subroutine domain_temperature

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pyramid_destruction()

    integer :: k

    ! To deallocate of the arrays
    do k = 1,level

      deallocate(pyramid_grid(k)%volume)
      deallocate(pyramid_grid(k)%cellsize)
      deallocate(pyramid_grid(k)%distance)
      deallocate(pyramid_grid(k)%solid_angle)
      deallocate(pyramid_grid(k)%case)
      deallocate(pyramid_grid(k)%transform)
	  	  	      
    end do

  end subroutine pyramid_destruction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cartiesian_destruction()

    integer :: k

    ! To deallocate of the arrays
    do k = 1,level

      deallocate(cartesian_grid(k)%transform)
	  	  	      
    end do
	
  end subroutine cartiesian_destruction
	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine domain_destruction()

    integer :: k

    ! To deallocate of the arrays
    do k = 1,level

      deallocate(pos_x_pos_y_pos_z_dir_z(k)%number_density)
      deallocate(pos_x_pos_y_pos_z_dir_z(k)%xHI)
      deallocate(pos_x_pos_y_pos_z_dir_z(k)%xHII)
      deallocate(pos_x_pos_y_pos_z_dir_z(k)%temperature)

    end do

  end subroutine domain_destruction
  
end module array