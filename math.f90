module math
	
	use input, only: level
	use array, only: partition
	
	implicit none
	
contains
	! to test if a trapezium cell is divided by more than 4 cubic cells
	subroutine math
		
		integer :: count,i,j,k
		real :: n, y1, y4
		integer :: yy1, yy4
		
		do k=1,level
		  do i=1,partition(k)
		    do j=1,partition(k)
				n = partition(k)
				y1=real(j-1)*real(k-1)/n
				y4=real(j)*real(k)/n
				if (abs(y1-int(y1)).ge.0.0001) then
					yy1=int(y1)+1
				else
					yy1=int(y1)
				endif
				yy4=int(y4)
				if (yy4-yy1.ge.2) write(*,*) 'error'
			enddo
		enddo
		enddo
		
	end subroutine math
	
end module math