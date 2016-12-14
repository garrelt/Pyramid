module case_division
	
  use precision, only: dp
  use input, only: level_py, partition	
  use array, only: pyramid_grid, cartesian_grid
	
  implicit none
	
contains
	
  subroutine pyramid_case

    implicit none
	
    integer :: i,j,k,n
    real(kind=dp) :: i1,i2,i3,i4
    real(kind=dp) :: j1,j2,j3,j4
    integer :: partition_i, partition_j
    integer :: a1,a2,a3
    integer :: b1,b2,b3
    real(kind=dp) :: ha1,ha2
    real(kind=dp) :: hb1,hb2
		
    ! Find out how Cartedian coordinate line partition the pyramid cell
    do k = 1,level_py
      do i = 1,partition(k)
        do j = 1,partition(k)
			
          partition_i = 0
          partition_j = 0
          n = real(partition(k))

          ! x = i1 is the line at the bottom left of pyramidal cell
          i1 = real(i-1)*real(k-1)/n
          ! x = i2 is the line at the top left of pyramidal cell
          i2 = real(i-1)*real(k)/n
          ! x = i3 is the line at the bottom right of pyramidal cell
          i3 = real(i)*real(k-1)/n
          ! x = i4 is the line at the top right of pyramidal cell
          i4 = real(i)*real(k)/n
          ! y = j1 is the line at the bottom down of pyramidal cell
          j1 = real(j-1)*real(k-1)/n
          ! y = j2 is the line at the top down of pyramidal cell
          j2 = real(j-1)*real(k)/n
          ! y = j3 is the line at the bottom up of pyramidal cell
          j3 = real(j)*real(k-1)/n
          ! y = j4 is the line at the top up of pyramidal cell
          j4 = real(j)*real(k)/n
			
          if (int(i4)-int(i1).eq.0) then
            partition_i = 1
            a1 = int(i1) + 1
          elseif (int(i4)-int(i1).eq.1 .and. i4-int(i4).le.0.00000001) then
            partition_i = 1	
            a1 = int(i4)	
          elseif (int(i4)-int(i1).eq.1 .and. i4-int(i4).gt.0.00000001) then
            partition_i = 2
            a1 = int(i4)
            a2 = int(i4) + 1
!come back to this later
           !if (i2.le.i3) then
            if (i1.le.a1 .and. a1.le.i2) then
              ha1 = real(a1)/real(i-1)
            elseif (i3.le.a1 .and. a1.le.i4) then
              ha1 = real(a1)/real(i)
            endif
!           else
!if (i3.lt.a1 .and. a1.lt.i2) write(*,*) 'ewsfewgferg'
!if (a1.lt.i3 .and. i3.lt.i2) write(*,*) 'yyyyrtret43654654'
!if (i3.lt.i2 .and. i2.lt.a1) write(*,*) 'ooooooooooooooooo'


!           end


          elseif (int(i4)-int(i1).eq.2 .and. i4-int(i4).le.0.00000001) then
            partition_i = 2
            a1 = int(i1) + 1
            a2 = int(i4)
            if (i1.le.a1 .and. a1.le.i2) then
              ha1 = real(a1)/real(i-1)
            elseif (i3.le.a1 .and. a1.le.i4) then
              ha1 = real(a1)/real(i)
            endif
          elseif (int(i4)-int(i1).eq.2 .and. i4-int(i4).gt.0.00000001) then
            partition_i = 3
            a1 = int(i1) + 1
            a2 = int(i4)
            a3 = int(i4) + 1	
            ha1 = real(a1)/real(i-1)
            ha2 = real(a2)/real(i)
          endif

          if (int(j4)-int(j1).eq.0) then
            partition_j = 1
            b1 = int(j1) + 1
          elseif (int(j4)-int(j1).eq.1 .and. j4-int(j4).le.0.00000001) then
            partition_j = 1	
            b1 = int(j4)	
          elseif (int(j4)-int(j1).eq.1 .and. j4-int(j4).gt.0.00000001) then
            partition_j = 2
            b1 = int(j4)
            b2 = int(j4) + 1
            if (j1.le.b1 .and. b1.le.j2) then
              hb1 = real(b1)/real(j-1)
            elseif (j3.le.b1 .and. b1.le.j4) then
              hb1 = real(b1)/real(j)	
            endif
          elseif (int(j4)-int(j1).eq.2 .and. j4-int(j4).le.0.00000001) then
            partition_j = 2
            b1 = int(j1) + 1
            b2 = int(j4)
            if (j1.le.b1 .and. b1.le.j2) then
              hb1 = real(b1)/real(j-1)
            elseif (j3.le.b1 .and. b1.le.j4) then
              hb1 = real(b1)/real(j)
            endif	
          elseif (int(j4)-int(j1).eq.2 .and. j4-int(j4).gt.0.00000001) then
            partition_j = 3
            b1 = int(j1) + 1
            b2 = int(j4)
            b3 = int(j4) + 1	
            hb1 = real(b1)/real(j-1)
            hb2 = real(b2)/real(j)	
          endif

          ! There are 42 different partition cases
          if (partition_i.eq.1 .and. partition_j.eq.1) then
            pyramid_grid(k)%case(i,j) = 1
            allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a1,b1:b1))
            pyramid_grid(k)%transform(i,j)%transform_volume = 0.0				
          elseif (partition_i.eq.2 .and. partition_j.eq.1) then
            if (i1.le.a1 .and. a1.le.i2) then
              pyramid_grid(k)%case(i,j) = 2
              allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a2,b1:b1))
              pyramid_grid(k)%transform(i,j)%transform_volume = 0.0				
            elseif (i2.le.a1 .and. a1.le.i3) then
              pyramid_grid(k)%case(i,j) = 3			
              allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a2,b1:b1))
              pyramid_grid(k)%transform(i,j)%transform_volume = 0.0					
            elseif (i3.le.a1 .and. a1.le.i4) then	
              pyramid_grid(k)%case(i,j) = 4
              allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a2,b1:b1))
              pyramid_grid(k)%transform(i,j)%transform_volume = 0.0	
            endif
          elseif (partition_i.eq.1 .and. partition_j.eq.2) then
            if (j1.le.b1 .and. b1.le.j2) then
              pyramid_grid(k)%case(i,j) = 5
              allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a1,b1:b2))
              pyramid_grid(k)%transform(i,j)%transform_volume = 0.0	
            elseif (j2.le.b1 .and. b1.le.j3) then
              pyramid_grid(k)%case(i,j) = 6
              allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a1,b1:b2))
              pyramid_grid(k)%transform(i,j)%transform_volume = 0.0	
            elseif (j3.le.b1 .and. b1.le.j4) then	
              pyramid_grid(k)%case(i,j) = 7
              allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a1,b1:b2))
              pyramid_grid(k)%transform(i,j)%transform_volume = 0.0	
            endif
          elseif (partition_i.eq.2 .and. partition_j.eq.2) then
            if (j1.le.b1 .and. b1.le.j2) then					
              if (i1.le.a1 .and. a1.le.i2) then	
                if (ha1.ge.hb1) then
                  pyramid_grid(k)%case(i,j) = 8
                  allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a2,b1:b2))
                  pyramid_grid(k)%transform(i,j)%transform_volume = 0.0	
                elseif (ha1.lt.hb1) then
                  pyramid_grid(k)%case(i,j) = 9
                  allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a2,b1:b2))
                  pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
                endif							
              elseif (i2.le.a1 .and. a1.le.i3) then		
                pyramid_grid(k)%case(i,j) = 10
                allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a2,b1:b2))
                  pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
              elseif (i3.le.a1 .and. a1.le.i4) then	
                if (ha1.ge.hb1) then
                  pyramid_grid(k)%case(i,j) = 11
                  allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a2,b1:b2))
                  pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
                elseif (ha1.lt.hb1) then
                  pyramid_grid(k)%case(i,j) = 12
                  allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a2,b1:b2))
                  pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
                endif										
              endif								
            elseif (j2.le.b1 .and. b1.le.j3) then							
              if (i1.le.a1 .and. a1.le.i2) then
                pyramid_grid(k)%case(i,j) = 13
                allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a2,b1:b2))
                pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
              elseif (i2.le.a1 .and. a1.le.i3) then		
                pyramid_grid(k)%case(i,j) = 14
                allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a2,b1:b2))
                pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
              elseif (i3.le.a1 .and. a1.le.i4) then
                pyramid_grid(k)%case(i,j) = 15	
                allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a2,b1:b2))
                pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
              endif
            elseif (j3.le.b1 .and. b1.le.j4) then						
              if (i1.le.a1 .and. a1.le.i2) then
                if (ha1.ge.hb1) then
                  pyramid_grid(k)%case(i,j) = 16
                  allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a2,b1:b2))
                  pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
                elseif (ha1.lt.hb1) then
                  pyramid_grid(k)%case(i,j) = 17
                  allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a2,b1:b2))
                  pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
                endif
              elseif (i2.le.a1 .and. a1.le.i3) then		
                pyramid_grid(k)%case(i,j) = 18	
                allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a2,b1:b2))
                pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
              elseif (i3.le.a1 .and. a1.le.i4) then
                if (ha1.ge.hb1) then
                  pyramid_grid(k)%case(i,j) = 19
                  allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a2,b1:b2))
                  pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
                elseif (ha1.lt.hb1) then
                  pyramid_grid(k)%case(i,j) = 20
                  allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a2,b1:b2))
                  pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
                endif
              endif											
            endif				
          elseif (partition_i.eq.3 .and. partition_j.eq.1) then
            pyramid_grid(k)%case(i,j) = 21
            allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a3,b1:b1))
            pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
          elseif (partition_i.eq.1 .and. partition_j.eq.3) then
            pyramid_grid(k)%case(i,j) = 22
            allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a1,b1:b3))
            pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
          elseif (partition_i.eq.3 .and. partition_j.eq.2) then				
            if (j1.le.b1 .and. b1.le.j2) then												
              if (ha1.ge.hb1) then
                pyramid_grid(k)%case(i,j) = 23
                allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a3,b1:b2))
                pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
              elseif (ha2.ge.hb1 .and. hb1.gt.ha1) then
                pyramid_grid(k)%case(i,j) = 24
                allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a3,b1:b2))
                pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
              elseif (hb1.gt.ha2) then
                pyramid_grid(k)%case(i,j) = 25	
                allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a3,b1:b2))
                pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
!write(*,*) 'i1',i1,'i2',i2,'i3',i3,'i4',i4,'a1',a1
              endif												
            elseif (j2.le.b1 .and. b1.le.j3) then
              pyramid_grid(k)%case(i,j) = 26
              allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a3,b1:b2))
              pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
            elseif (j3.le.b1 .and. b1.le.j4) then						
              if (ha1.ge.hb1) then
                pyramid_grid(k)%case(i,j) = 27
                allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a3,b1:b2))
                pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
              elseif (ha2.ge.hb1 .and. hb1.gt.ha1) then
                pyramid_grid(k)%case(i,j) = 28	
                allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a3,b1:b2))
                pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
              elseif (hb1.gt.ha2) then
                pyramid_grid(k)%case(i,j) = 29
                allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a3,b1:b2))
                pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
              endif
            endif
          elseif (partition_i.eq.2 .and. partition_j.eq.3) then				
            if (i1.le.a1 .and. a1.le.i2) then					
              if (hb1.ge.ha1) then
                pyramid_grid(k)%case(i,j) = 30
                allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a2,b1:b3))
                pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
              elseif (hb2.ge.ha1 .and. ha1.gt.hb1) then
                pyramid_grid(k)%case(i,j) = 31
                allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a2,b1:b3))
                pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
              elseif (ha1.gt.hb2) then
                pyramid_grid(k)%case(i,j) = 32
                allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a2,b1:b3))
                pyramid_grid(k)%transform(i,j)%transform_volume = 0.0						
              endif
            elseif (i2.le.a1 .and. a1.le.i3) then					
              pyramid_grid(k)%case(i,j) = 33
              allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a2,b1:b3))
              pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
            elseif (i3.le.a1 .and. a1.le.i4) then					
              if (hb1.ge.ha1) then
                pyramid_grid(k)%case(i,j) = 34
                allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a2,b1:b3))
                pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
              elseif (hb2.ge.ha1 .and. ha1.gt.hb1) then
                pyramid_grid(k)%case(i,j) = 35
                allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a2,b1:b3))
                pyramid_grid(k)%transform(i,j)%transform_volume = 0.0							
              elseif (ha1.gt.hb2) then
                pyramid_grid(k)%case(i,j) = 36
                allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a2,b1:b3))
                pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
              endif													
            endif				
          elseif (partition_i.eq.3 .and. partition_j.eq.3) then								
            if (ha1.ge.hb2) then
              pyramid_grid(k)%case(i,j) = 37
              allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a3,b1:b3))
              pyramid_grid(k)%transform(i,j)%transform_volume = 0.0						
            elseif (ha2.ge.hb2 .and. hb2.gt.ha1 .and. ha1.ge.hb1) then			
              pyramid_grid(k)%case(i,j) = 38
              allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a3,b1:b3))
              pyramid_grid(k)%transform(i,j)%transform_volume = 0.0
            elseif (hb2.gt.ha2 .and. ha1.ge.hb1) then			
              pyramid_grid(k)%case(i,j) = 39
              allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a3,b1:b3))
              pyramid_grid(k)%transform(i,j)%transform_volume = 0.0						
            elseif (hb2.gt.ha2 .and. ha2.ge.hb1 .and. hb1.gt.ha1) then			
              pyramid_grid(k)%case(i,j) = 40
              allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a3,b1:b3))
              pyramid_grid(k)%transform(i,j)%transform_volume = 0.0						
            elseif (ha2.ge.hb2 .and. hb1.gt.ha1) then			
              pyramid_grid(k)%case(i,j) = 41
              allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a3,b1:b3))
              pyramid_grid(k)%transform(i,j)%transform_volume = 0.0						
            elseif (hb1.gt.ha2) then			
              pyramid_grid(k)%case(i,j) = 42
              allocate(pyramid_grid(k)%transform(i,j)%transform_volume(a1:a3,b1:b3))
              pyramid_grid(k)%transform(i,j)%transform_volume = 0.0						
            endif										
          endif	
		  							
        end do
      end do
    end do
	
  end subroutine pyramid_case
  
  subroutine cartesian_case
  
    implicit none
	
    integer :: x,y,z,n
    integer :: min_i, max_i	
    integer :: min_j, max_j
	
    !z=1 is treated separatedly
    allocate(cartesian_grid(1)%transform(1,1)%transform_volume(1,1))
    cartesian_grid(1)%transform(1,1)%transform_volume(1,1) = 0.0
				
    do z = 2,level_py
      n = real(partition(z))
      do x = 1,z
        do y = 1,z
			
          min_i = int(n*real(x-1)/real(z))+1
          if (n*real(x)/real(z-1).gt.int(n*real(x)/real(z-1))) then
            max_i = int(n*real(x)/real(z-1)+1.0)
          elseif (n*real(x)/real(z-1).eq.int(n*real(x)/real(z-1))) then
            max_i = int(n*real(x)/real(z-1))						
          endif

          min_j = int(n*real(y-1)/real(z))+1
          if (n*real(y)/real(z-1).gt.int(n*real(y)/real(z-1))) then
            max_j = int(n*real(y)/real(z-1)+1.0)
          elseif (n*real(y)/real(z-1).eq.int(n*real(y)/real(z-1))) then
            max_j = int(n*real(y)/real(z-1))						
          endif
		
          allocate(cartesian_grid(z)%transform(x,y)%transform_volume(min_i:max_i,min_j:max_j))
          cartesian_grid(z)%transform(x,y)%transform_volume = 0.0		

        enddo
      enddo
    enddo

  end subroutine cartesian_case	
  
end module case_division
