module pyramid_initialization
	
  use precision, only: dp
  use input, only: level_py, partition, layer, pi, cellsize_py, cellsize_py_cube
  use array, only: pyramid_grid, parent_position
  use solid_angle, only: ss, cc
  
  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pyramid_volume()
    implicit none
    integer :: k, n
    real(kind=dp) :: volume

    ! Volume on the same level are the same
    do k = 1,level_py

      n = partition(k)
      volume = real(3*k*k-3*k+1)/real(3*n*n) * cellsize_py_cube
      pyramid_grid(k)%volume = volume

    end do

  end subroutine pyramid_volume

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pyramid_cellsize()
    implicit none
    integer :: i,j,k,n

    !The radial length of the pyramid_grid
    do k = 1,level_py
      do i = 1,partition(k)
        do j = 1,partition(k)

          n = partition(k)
          pyramid_grid(k)%cellsize(i,j) = cellsize_py*sqrt(real(i*i+j*j+n*n-i-j)+0.5)/real(n)

        end do
      end do
    end do

  end subroutine pyramid_cellsize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pyramid_solid_angle()
    implicit none
    integer :: i,j,k,n
    real(kind=dp) :: four_pi
	
    four_pi = 4*pi
    n = partition(level_py)
	
    ! solid angle spanned by the pyramid_grid of the toppest level
    do i = 1,partition(level_py)
      do j = 1,partition(level_py)
        if (i.ge.j) then
			
          pyramid_grid(level_py)%solid_angle(i,j) = (ss(i-1,j,i-1,n) - ss(i-1,j-1,i-1,n) - &
                                                    ss(i,j,i,n) + ss(i,j-1,i,n) + &
                                                    cc(j,j,i,n) - cc(j,j,i-1,n) + &
                                                    cc(j-1,j-1,i-1,n) - cc(j-1,j-1,i,n))/four_pi
					
        elseif (i.lt.j) then
			
          pyramid_grid(level_py)%solid_angle(i,j) = (ss(j-1,i,j-1,n) - ss(j-1,i-1,j-1,n) - &
                                                    ss(j,i,j,n) + ss(j,i-1,j,n) + &
                                                    cc(i,i,j,n) - cc(i,i,j-1,n) + &
                                                    cc(i-1,i-1,j-1,n) - cc(i-1,i-1,j,n))/four_pi
				
        endif	
      end do
    end do
	
    ! solid angle spanned by the pyramid_grid of the other level
    do k = level_py-1,1,-1
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
  
  subroutine root_setup()

    implicit none
  
    integer :: son_position, son_height, parent_height
  
    do son_height = 1,level_py
      do son_position = 1,partition(son_height)
        do parent_height = 1,son_height		   
          parent_position(son_height,son_position)%with_parent_height(parent_height) = &
                 (son_position-1)/(2**(layer(son_height)-layer(parent_height)))+1 
        enddo
      enddo
    enddo
	 	  
  end subroutine root_setup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pyramid_initialization
