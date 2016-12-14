module geometry_transformation
	
  use precision, only: dp
  use input, only: level_py, partition, cellsize_py_cube
  use array, only: pyramid_grid, cartesian_grid
  use integration, only: volume_integration

contains
	
  subroutine volume_partition

    implicit none
  
    integer :: i,j,k,n
    real(kind=dp) :: i1,i2,i3,i4
    real(kind=dp) :: j1,j2,j3,j4
    integer :: a1,a2,a3
    integer :: b1,b2,b3
    real(kind=dp) :: ha1,ha2
    real(kind=dp) :: hb1,hb2
    real(kind=dp) :: temp
    real(kind=dp) :: a1_dp,a2_dp
    real(kind=dp) :: b1_dp,b2_dp
    real(kind=dp) :: i1s_dp, i2s_dp
    real(kind=dp) :: j1s_dp, j2s_dp	
    real(kind=dp) :: k1s_dp, k2s_dp	
    real(kind=dp) :: relative_difference	
    ! these variables only used to check if the code is correct, must erase when finish
    real(kind=dp) ::sumsum, sumsum2
    integer :: i_primary, i_secondary
    integer ::x,y,z
    integer :: xx,yy,zz
    integer :: min_i,max_i,min_j,max_j

    ! initialization
    a1 = 0
    a2 = 0
    a3 = 0
    b1 = 0
    b2 = 0
    b3 = 0
    ha1 = 0
    ha2 = 0
    hb1 = 0
    hb2 = 0
	
    ! Find out how Cartedian coordinate line partition the pyramid cell
    do k = 1,level_py
      do i = 1,partition(k)
        do j = 1,partition(k)

          n = real(partition(k))
          i1 = real(i-1)*real(k-1)/n
          i2 = real(i-1)*real(k)/n
          i3 = real(i)*real(k-1)/n
          i4 = real(i)*real(k)/n
          j1 = real(j-1)*real(k-1)/n
          j2 = real(j-1)*real(k)/n
          j3 = real(j)*real(k-1)/n
          j4 = real(j)*real(k)/n
			
          if (int(i4)-int(i1).eq.0) then
            a1 = int(i1)+1
          elseif (int(i4)-int(i1).eq.1 .and. i4-int(i4).le.0.00000001) then	
            a1 = int(i4)	
          elseif (int(i4)-int(i1).eq.1 .and. i4-int(i4).gt.0.00000001) then
            a1 = int(i4)
            a2 = int(i4)+1
            if (i1.le.a1 .and. a1.le.i2) then ! comparison fp and int
              ha1 = real(n)*real(a1)/real(i-1)
            elseif (i3.le.a1 .and. a1.le.i4) then
              ha1 = real(n)*real(a1)/real(i)	
            endif
          elseif (int(i4)-int(i1).eq.2 .and. i4-int(i4).le.0.00000001) then
            a1 = int(i1)+1
            a2 = int(i4)
            if (i1.le.a1 .and. a1.le.i2) then
              ha1 = real(n)*real(a1)/real(i-1)
            elseif (i3.le.a1 .and. a1.le.i4) then
              ha1 = real(n)*real(a1)/real(i)	
            endif
          elseif (int(i4)-int(i1).eq.2 .and. i4-int(i4).gt.0.00000001) then
            a1 = int(i1)+1
            a2 = int(i4)
            a3 = int(i4)+1	
            ha1 = real(n)*real(a1)/real(i-1)
            ha2 = real(n)*real(a2)/real(i)	
          endif


          if (int(j4)-int(j1).eq.0) then
            b1 = int(j1)+1
          elseif (int(j4)-int(j1).eq.1 .and. j4-int(j4).le.0.00000001) then	
            b1 = int(j4)	
          elseif (int(j4)-int(j1).eq.1 .and. j4-int(j4).gt.0.00000001) then
            b1 = int(j4)
            b2 = int(j4)+1
            if (j1.le.b1 .and. b1.le.j2) then
              hb1 = real(n)*real(b1)/real(j-1)
            elseif (j3.le.b1 .and. b1.le.j4) then
              hb1 = real(n)*real(b1)/real(j)	
            endif
          elseif (int(j4)-int(j1).eq.2 .and. j4-int(j4).le.0.00000001) then
            b1 = int(j1)+1
            b2 = int(j4)
            if (j1.le.b1 .and. b1.le.j2) then
              hb1 = real(n)*real(b1)/real(j-1)
            elseif (j3.le.b1 .and. b1.le.j4) then
              hb1 = real(n)*real(b1)/real(j)
            endif	
          elseif (int(j4)-int(j1).eq.2 .and. j4-int(j4).gt.0.00000001) then
            b1 = int(j1)+1
            b2 = int(j4)
            b3 = int(j4)+1	
            hb1 = real(n)*real(b1)/real(j-1)
            hb2 = real(n)*real(b2)/real(j)	
          endif
		  		  
          a1_dp = real(a1)
          a2_dp = real(a2)
          b1_dp = real(b1)
          b2_dp = real(b2)		  
          i1s_dp = real(i-1)/n		  
          i2s_dp = real(i)/n
          j1s_dp = real(j-1)/n		  
          j2s_dp = real(j)/n		  
          k1s_dp = real(k-1)		  
          k2s_dp = real(k)		  
		  		  		  
          if (pyramid_grid(k)%case(i,j).eq.1) then
		
            temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,k1s_dp,k2s_dp,1,1)

            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube

relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

          elseif (pyramid_grid(k)%case(i,j).eq.2) then

            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,k1s_dp,ha1,2,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube	

            temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha1,k2s_dp,3,1)+&				
            volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,k1s_dp,ha1,3,2)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube	
	
relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif
	
          elseif (pyramid_grid(k)%case(i,j).eq.3) then
		
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,k1s_dp,k2s_dp,4,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube		
			
            temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,k1s_dp,k2s_dp,5,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube		

relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

          elseif (pyramid_grid(k)%case(i,j).eq.4) then
	
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha1,k2s_dp,6,1)+&		
            volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,k1s_dp,ha1,6,2)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube		
				
            temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha1,k2s_dp,7,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube	
	 
relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif
	
          elseif (pyramid_grid(k)%case(i,j).eq.5) then
	
            temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,hb1,8,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube			
		
            temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,k2s_dp,9,1)+&		
            volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,hb1,9,2)	
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube	
	
relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

          elseif (pyramid_grid(k)%case(i,j).eq.6) then

            temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,k2s_dp,10,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube	
			
            temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,k2s_dp,11,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube			
	
relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif
	
          elseif (pyramid_grid(k)%case(i,j).eq.7) then

            temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,hb1,k2s_dp,12,1)+&
            volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,k1s_dp,hb1,12,2)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube		
			
            temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,hb1,k2s_dp,13,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube	
	
relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif
	
          elseif (pyramid_grid(k)%case(i,j).eq.8) then
				
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,hb1,14,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube		
								
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,ha1,15,1)+&				
            volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,hb1,15,2)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube		
			
            temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,hb1,16,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube		
			
            temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha1,k2s_dp,17,1)+&				
            volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,ha1,17,2)+&
            volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,hb1,17,3)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube	
            cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube	
	
relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

          elseif (pyramid_grid(k)%case(i,j).eq.9) then
				
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,ha1,18,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube	
							
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,ha1,19,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube	
				
            temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha1,hb1,20,1)+&
            volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,ha1,20,2)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube	
				
            temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,k2s_dp,21,1)+&				
            volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha1,hb1,21,2)+&
            volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,ha1,21,3)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube	

relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

          elseif (pyramid_grid(k)%case(i,j).eq.10) then
				
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,hb1,22,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube	
												
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,k2s_dp,23,1)+&
            volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,hb1,23,2)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube	
			
            temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,hb1,24,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube	
				
            temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,k2s_dp,25,1)+&
            volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,hb1,25,2)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube	

relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

          elseif (pyramid_grid(k)%case(i,j).eq.11) then
						
            temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,hb1,26,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube	
							
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha1,k2s_dp,27,1)+&
            volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,ha1,27,2)+&
            volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,hb1,27,3)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube	
	
            temp = 0.0
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube
			
            temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha1,k2s_dp,28,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube	

relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

          elseif (pyramid_grid(k)%case(i,j).eq.12) then
					
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha1,hb1,29,1)+&
            volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,ha1,29,2)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube	
								
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,k2s_dp,30,1)+&
            volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha1,hb1,30,2)+&
            volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,ha1,30,3)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube	
				
            temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha1,hb1,31,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube	
				
            temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,k2s_dp,32,1)+&
            volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha1,hb1,32,2)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube	
	
relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

          elseif (pyramid_grid(k)%case(i,j).eq.13) then
					
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,ha1,33,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube		
							
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,ha1,34,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube		
			
            temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha1,k2s_dp,35,1)+&
            volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,ha1,35,2)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube		
			
            temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha1,k2s_dp,36,1)+&
            volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,ha1,36,2)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube		

relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

          elseif (pyramid_grid(k)%case(i,j).eq.14) then
			
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,k2s_dp,37,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube	
            cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube						
							
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,k2s_dp,38,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube		

            temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,k2s_dp,39,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube		
			
            temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,k2s_dp,40,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube	

relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

          elseif (pyramid_grid(k)%case(i,j).eq.15) then
					
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha1,k2s_dp,41,1)+&
            volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,ha1,41,2)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube				
								
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha1,k2s_dp,42,1)+&
            volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,ha1,42,2)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube		

            temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha1,k2s_dp,43,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube	
			
            temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha1,k2s_dp,44,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube
	
relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

          elseif (pyramid_grid(k)%case(i,j).eq.16) then
				
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,hb1,ha1,45,1)+&
            volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,k1s_dp,hb1,45,2)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube	
            cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube				
							
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,hb1,ha1,46,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube	

            temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha1,k2s_dp,47,1)+&
            volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,hb1,ha1,47,2)+&
            volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,k1s_dp,hb1,47,3)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube			

            temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha1,k2s_dp,48,1)+&				
            volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,hb1,ha1,48,2)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube
	
relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

          elseif (pyramid_grid(k)%case(i,j).eq.17) then
			
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,k1s_dp,ha1,49,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube
		
            temp = 0.0
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube		

            temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,hb1,k2s_dp,50,1)+&
            volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha1,hb1,50,2)+&
            volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,k1s_dp,ha1,50,3)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube				

            temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,hb1,k2s_dp,51,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube		

relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

          elseif (pyramid_grid(k)%case(i,j).eq.18) then
				
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,hb1,k2s_dp,52,1)+&
            volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,k1s_dp,hb1,52,2)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube				
							
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,hb1,k2s_dp,53,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube	

            temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,hb1,k2s_dp,54,1)+&
            volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,k1s_dp,hb1,54,2)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube		
				
            temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,hb1,k2s_dp,55,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube
	
relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

          elseif (pyramid_grid(k)%case(i,j).eq.19) then
						
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha1,k2s_dp,56,1)+&
            volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,hb1,ha1,56,2)+&
            volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,k1s_dp,hb1,56,3)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube				
							
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha1,k2s_dp,57,1)+&
            volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,hb1,ha1,57,2)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube			

            temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha1,k2s_dp,58,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube			
				
            temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha1,k2s_dp,59,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube	

relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

          elseif (pyramid_grid(k)%case(i,j).eq.20) then
				
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,hb1,k2s_dp,60,1)+&
            volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha1,hb1,60,2)+&
            volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,k1s_dp,ha1,60,3)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube					
							
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,hb1,k2s_dp,61,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube			

            temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,hb1,k2s_dp,62,1)+&
            volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha1,hb1,62,2)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube			
				
            temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,hb1,k2s_dp,63,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube	
	
relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

          elseif (pyramid_grid(k)%case(i,j).eq.21) then
											
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,k1s_dp,ha1,64,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube			
			
            temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha2,k2s_dp,65,1)+&
            volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha1,ha2,65,2)+&
            volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,k1s_dp,ha1,65,3)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube	
            cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube			
				
            temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha2,k2s_dp,66,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a3,b1)%transform_volume(i,j) = temp * cellsize_py_cube		
	
relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif
	
          elseif (pyramid_grid(k)%case(i,j).eq.22) then
					
            temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,hb1,67,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube	
            cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube			
			
            temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,hb2,k2s_dp,68,1)+&
            volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,hb2,68,2)+&
            volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,hb1,68,3)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube

            temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,hb2,k2s_dp,69,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3) = temp * cellsize_py_cube	
            cartesian_grid(k)%transform(a1,b3)%transform_volume(i,j) = temp * cellsize_py_cube		

relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

          elseif (pyramid_grid(k)%case(i,j).eq.23) then

            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,hb1,70,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube						
							
            temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,ha1,71,1)+&
            volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,hb1,71,2)
            pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube			
			
            temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,hb1,72,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube	
            cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube			
				
            temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha2,k2s_dp,73,1)+&
            volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha1,ha2,73,2)+&
            volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,ha1,73,3)+&
            volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,hb1,73,4)
            pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube	
			
            temp = 0.0
            pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a3,b1)%transform_volume(i,j) = temp * cellsize_py_cube		
					
            temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha2,k2s_dp,74,1)
            pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2) = temp * cellsize_py_cube
            cartesian_grid(k)%transform(a3,b2)%transform_volume(i,j) = temp * cellsize_py_cube		
	
relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif
		
          elseif (pyramid_grid(k)%case(i,j).eq.24) then

temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,ha1,75,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube					
		
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,ha1,76,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube		
				
temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha1,hb1,77,1)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,ha1,77,2)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube		
				
temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha2,k2s_dp,78,1)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,ha2,78,2)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha1,hb1,78,3)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,ha1,78,4)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube
		
temp = 0.0
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b1)%transform_volume(i,j) = temp * cellsize_py_cube	
					
temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha2,k2s_dp,79,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b2)%transform_volume(i,j) = temp * cellsize_py_cube

relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

elseif (pyramid_grid(k)%case(i,j).eq.25) then
	
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,ha1,80,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube					
			
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,ha1,81,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube			
					
temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha2,hb1,82,1)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha1,ha2,82,2)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,ha1,82,3)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube			
				
temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,k2s_dp,83,1)+&
volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha2,hb1,83,2)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha1,ha2,83,3)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,ha1,83,4)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube			
					
temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha2,hb1,84,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b1)%transform_volume(i,j) = temp * cellsize_py_cube	
				
temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,k2s_dp,85,1)+&
volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha2,hb1,85,2)
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b2)%transform_volume(i,j) = temp * cellsize_py_cube		

relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

elseif (pyramid_grid(k)%case(i,j).eq.26) then
	
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,ha1,86,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube				
		
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,ha1,87,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube		
				
temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha2,k2s_dp,88,1)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha1,ha2,88,2)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,ha1,88,3)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube		
			
temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha2,k2s_dp,89,1)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha1,ha2,89,2)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,ha1,89,3)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube	
cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube		
					
temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha2,k2s_dp,90,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b1)%transform_volume(i,j) = temp * cellsize_py_cube	
			
temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha2,k2s_dp,91,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2) = temp * cellsize_py_cube	
cartesian_grid(k)%transform(a3,b2)%transform_volume(i,j) = temp * cellsize_py_cube	

relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

elseif (pyramid_grid(k)%case(i,j).eq.27) then
	
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,hb1,ha1,92,1)+&	
volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,k1s_dp,hb1,92,2)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube				
		
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,hb1,ha1,93,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube		
				
temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha2,k2s_dp,94,1)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha1,ha2,94,2)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,hb1,ha1,94,3)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,k1s_dp,hb1,94,4)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube			
				
temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha2,k2s_dp,95,1)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha1,ha2,95,2)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,hb1,ha1,95,3)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube		
					
temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha2,k2s_dp,96,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b1)%transform_volume(i,j) = temp * cellsize_py_cube		
				
temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha2,k2s_dp,97,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b2)%transform_volume(i,j) = temp * cellsize_py_cube	

relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

elseif (pyramid_grid(k)%case(i,j).eq.28) then
					
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,k1s_dp,ha1,98,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube
			
temp = 0.0
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube			
					
temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha2,k2s_dp,99,1)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,hb1,ha2,99,2)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha1,hb1,99,3)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,k1s_dp,ha1,99,4)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube		
				
temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha2,k2s_dp,100,1)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,hb1,ha2,100,2)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube		
					
temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha2,k2s_dp,101,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b1)%transform_volume(i,j) = temp * cellsize_py_cube	
					
temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha2,k2s_dp,102,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b2)%transform_volume(i,j) = temp * cellsize_py_cube	

relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

elseif (pyramid_grid(k)%case(i,j).eq.29) then
						
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,k1s_dp,ha1,103,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube
				
temp = 0.0
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp	* cellsize_py_cube	
					
temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,hb1,k2s_dp,104,1)+&
volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha2,hb1,104,2)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha1,ha2,104,3)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,k1s_dp,ha1,104,4)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube		
				
temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,hb1,k2s_dp,105,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube		
				
temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,hb1,k2s_dp,106,1)+&
volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha2,hb1,106,2)
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b1)%transform_volume(i,j) = temp * cellsize_py_cube		
				
temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,hb1,k2s_dp,107,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b2)%transform_volume(i,j) = temp * cellsize_py_cube	
	
relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

elseif (pyramid_grid(k)%case(i,j).eq.30) then
						
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,ha1,108,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube	
			
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,ha1,109,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube	
cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube	
	
temp = 0.0	
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b3)%transform_volume(i,j) = temp * cellsize_py_cube		
				
temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha1,hb1,110,1)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,ha1,110,2)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp
cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp	
			
temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,hb2,k2s_dp,111,1)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,hb2,111,2)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha1,hb1,111,3)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,ha1,111,4)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube	
cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube						
		
temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,hb2,k2s_dp,112,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b3)%transform_volume(i,j) = temp * cellsize_py_cube	

relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

elseif (pyramid_grid(k)%case(i,j).eq.31) then
							
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,hb1,113,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube
				
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,ha1,114,1)+&
volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,hb1,114,2)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube	

temp = 0.0
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b3)%transform_volume(i,j) = temp * cellsize_py_cube		
		
temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,hb1,115,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube	
cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube		
				
temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,hb2,k2s_dp,116,1)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha1,hb2,116,2)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,ha1,116,3)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,hb1,116,4)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube				
		
temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,hb2,k2s_dp,117,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b3)%transform_volume(i,j) = temp * cellsize_py_cube	
	
relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

elseif (pyramid_grid(k)%case(i,j).eq.32) then
						
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,hb1,119,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube	
cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube	
				
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,hb2,ha1,120,1)+&
volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,hb2,120,2)+&
volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,hb1,120,3)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube					
			
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,hb2,ha1,121,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b3)%transform_volume(i,j) = temp * cellsize_py_cube	
				
temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,hb1,122,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube	
					
temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,ha1,k2s_dp,123,1)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,hb2,ha1,123,2)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,hb2,123,3)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,hb1,123,4)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube				
		
temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,ha1,k2s_dp,124,1)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,hb2,ha1,124,2)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b3)%transform_volume(i,j) = temp * cellsize_py_cube		
	
relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

elseif (pyramid_grid(k)%case(i,j).eq.33) then
								
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,hb1,125,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube		
				
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,hb2,k2s_dp,126,1)+&
volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,hb2,126,2)+&
volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,hb1,126,3)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube			
			
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,hb2,k2s_dp,127,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b3)%transform_volume(i,j) = temp * cellsize_py_cube	
			
temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,hb1,129,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube	
				
temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,hb2,k2s_dp,130,1)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,hb2,130,2)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,hb1,130,3)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube				
	
temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,hb2,k2s_dp,131,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3) = temp * cellsize_py_cube	
cartesian_grid(k)%transform(a2,b3)%transform_volume(i,j) = temp * cellsize_py_cube		
	
relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

elseif (pyramid_grid(k)%case(i,j).eq.34) then
					
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha1,hb1,132,1)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,ha1,132,2)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube	
cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube			
					
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,hb2,k2s_dp,133,1)+&
volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,hb2,133,2)+&
volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha1,hb1,133,3)+&	
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,ha1,133,4)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube					
	
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,hb2,k2s_dp,134,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3) = temp * cellsize_py_cube	
cartesian_grid(k)%transform(a1,b3)%transform_volume(i,j) = temp * cellsize_py_cube		
				
temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha1,hb1,135,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube	
				
temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,hb2,k2s_dp,136,1)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,hb2,136,2)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha1,hb1,136,3)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube			
		
temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,hb2,k2s_dp,137,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b3)%transform_volume(i,j) = temp	 * cellsize_py_cube
	
relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

elseif (pyramid_grid(k)%case(i,j).eq.35) then
				
temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,hb1,138,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube	
cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube		
					
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,hb2,k2s_dp,139,1)+&
volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha1,hb2,139,2)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,ha1,139,3)+&	
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,hb1,139,4)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube						
		
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,hb2,k2s_dp,140,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b3)%transform_volume(i,j) = temp * cellsize_py_cube		

temp = 0.0	
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube			
					
temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,hb2,k2s_dp,141,1)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha1,hb2,141,2)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube			
		
temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,hb2,k2s_dp,142,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b3)%transform_volume(i,j) = temp * cellsize_py_cube	
	
relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

elseif (pyramid_grid(k)%case(i,j).eq.36) then
				
temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,hb1,143,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube				
				
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,ha1,k2s_dp,144,1)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,hb2,ha1,144,2)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,hb2,144,3)+&	
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,hb1,144,4)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube					
		
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,ha1,k2s_dp,145,1)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,hb2,ha1,145,2)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b3)%transform_volume(i,j) = temp * cellsize_py_cube			

temp = 0.0
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp	* cellsize_py_cube		
					
temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,ha1,k2s_dp,146,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube				
		
temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,ha1,k2s_dp,147,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b3)%transform_volume(i,j) = temp * cellsize_py_cube		
	
relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

elseif (pyramid_grid(k)%case(i,j).eq.37) then
			
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,hb1,148,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube	

temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,hb2,ha1,149,1)+&
volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,hb2,149,2)+&			
volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,hb1,149,3)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp	* cellsize_py_cube					
		
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,hb2,ha1,150,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b3)%transform_volume(i,j) = temp * cellsize_py_cube	
				
temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,hb1,151,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube

temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,ha2,k2s_dp,152,1)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,ha1,ha2,152,2)+&			
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,hb2,ha1,152,3)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,hb2,152,4)+&		
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,hb1,152,5)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube				

temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,ha2,k2s_dp,153,1)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,ha1,ha2,153,2)+&			
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,hb2,ha1,153,3)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b3)%transform_volume(i,j) = temp * cellsize_py_cube	

temp = 0.0	
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b1)%transform_volume(i,j) = temp * cellsize_py_cube				

temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,ha2,k2s_dp,154,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b2)%transform_volume(i,j) = temp * cellsize_py_cube						
			
temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,ha2,k2s_dp,155,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b3)%transform_volume(i,j) = temp * cellsize_py_cube		
	
relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

elseif (pyramid_grid(k)%case(i,j).eq.38) then
			
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,hb1,156,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube		

temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,ha1,157,1)+&			
volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,hb1,157,2)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube							

temp = 0.0
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b3)%transform_volume(i,j) = temp * cellsize_py_cube		
				
temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,hb1,158,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube		

temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,ha2,k2s_dp,159,1)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,hb2,ha2,159,2)+&			
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha1,hb2,159,3)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,ha1,159,4)+&		
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,hb1,159,5)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube			
	
temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,ha2,k2s_dp,160,1)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,hb2,ha2,160,2)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b3)%transform_volume(i,j) = temp * cellsize_py_cube	
	
temp = 0.0
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b1)%transform_volume(i,j) = temp * cellsize_py_cube		

temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,ha2,k2s_dp,161,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b2)%transform_volume(i,j) = temp * cellsize_py_cube					
			
temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,ha2,k2s_dp,162,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b3)%transform_volume(i,j) = temp * cellsize_py_cube		
	
relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

elseif (pyramid_grid(k)%case(i,j).eq.39) then
					
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,hb1,163,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube			

temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,ha1,164,1)+&			
volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,hb1,164,2)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube					
	
temp = 0.0
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b3)%transform_volume(i,j) = temp * cellsize_py_cube			
				
temp = volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,hb1,165,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube			

temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,hb2,k2s_dp,166,1)+&
volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha2,hb2,166,2)+&			
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha1,ha2,166,3)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,ha1,166,4)+&		
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,hb1,166,5)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube		
	
temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,hb2,k2s_dp,167,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b3)%transform_volume(i,j) = temp * cellsize_py_cube			
	
temp = 0.0
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b1)%transform_volume(i,j) = temp * cellsize_py_cube			

temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,hb2,k2s_dp,168,1)+&
volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha2,hb2,168,2)
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b2)%transform_volume(i,j) = temp * cellsize_py_cube				
		
temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,hb2,k2s_dp,169,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b3)%transform_volume(i,j) = temp * cellsize_py_cube		
	
relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

elseif (pyramid_grid(k)%case(i,j).eq.40) then
				
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,ha1,170,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube	
cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube			
				
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,ha1,171,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube				
	
temp = 0.0
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b3)%transform_volume(i,j) = temp * cellsize_py_cube			
	
temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha1,hb1,172,1)+&					
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,ha1,172,2)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube		

temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,hb2,k2s_dp,173,1)+&
volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha2,hb2,173,2)+&			
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,ha2,173,3)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha1,hb1,173,4)+&		
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,ha1,173,5)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube		

temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,hb2,k2s_dp,174,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b3)%transform_volume(i,j) = temp * cellsize_py_cube		
	
temp = 0.0
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b1)%transform_volume(i,j) = temp * cellsize_py_cube		
	
temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,hb2,k2s_dp,175,1)+&
volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,ha2,hb2,175,2)
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b2)%transform_volume(i,j) = temp * cellsize_py_cube				
		
temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,hb2,k2s_dp,176,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b3)%transform_volume(i,j) = temp * cellsize_py_cube		
	
relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

elseif (pyramid_grid(k)%case(i,j).eq.41) then
						
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,ha1,177,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube		
				
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,ha1,178,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube					

temp = 0.0
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b3)%transform_volume(i,j) = temp * cellsize_py_cube		

temp = volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha1,hb1,179,1)+&					
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,ha1,179,2)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube	

temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,ha2,k2s_dp,180,1)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,hb2,ha2,180,2)+&			
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,hb2,180,3)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha1,hb1,180,4)+&		
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,ha1,180,5)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube				

temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,ha2,k2s_dp,181,1)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,hb2,ha2,181,2)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b3)%transform_volume(i,j) = temp * cellsize_py_cube	

temp = 0.0
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b1)%transform_volume(i,j) = temp * cellsize_py_cube		

temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,ha2,k2s_dp,182,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b2)%transform_volume(i,j) = temp * cellsize_py_cube							
		
temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,ha2,k2s_dp,183,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b3)%transform_volume(i,j) = temp * cellsize_py_cube			
	
relative_difference = (sum(sum(pyramid_grid(k)%transform(i,j)%transform_volume,1),1)-pyramid_grid(k)%volume)&
/pyramid_grid(k)%volume
if(relative_difference.ge.0.0001) then
write(*,*) 'case ',pyramid_grid(k)%case(i,j),' relative difference is ',relative_difference
endif

elseif (pyramid_grid(k)%case(i,j).eq.42) then
				
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,ha1,184,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b1)%transform_volume(i,j) = temp * cellsize_py_cube			
				
temp = volume_integration(0.0_dp,i1s_dp,a1_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,ha1,185,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b2)%transform_volume(i,j) = temp * cellsize_py_cube			

temp = 0.0
pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a1,b3)%transform_volume(i,j) = temp * cellsize_py_cube			
	
temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha2,hb1,186,1)+&					
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha1,ha2,186,2)+&
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,k1s_dp,ha1,186,3)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b1)%transform_volume(i,j) = temp * cellsize_py_cube		

temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,hb2,k2s_dp,187,1)+&
volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,hb2,187,2)+&			
volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha2,hb1,187,3)+&
volume_integration(0.0_dp,i1s_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha1,ha2,187,4)+&		
volume_integration(a1_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,k1s_dp,ha1,187,5)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b2)%transform_volume(i,j) = temp * cellsize_py_cube			

temp = volume_integration(0.0_dp,i1s_dp,a2_dp,0.0_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,hb2,k2s_dp,188,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a2,b3)%transform_volume(i,j) = temp * cellsize_py_cube			
			
temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b1_dp,0.0_dp,ha2,hb1,189,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b1)%transform_volume(i,j) = temp * cellsize_py_cube			
	
temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,b2_dp,0.0_dp,hb2,k2s_dp,190,1)+&
volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,0.0_dp,j1s_dp,0.0_dp,j2s_dp,hb1,hb2,190,2)+&
volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,b1_dp,0.0_dp,0.0_dp,j2s_dp,ha2,hb1,190,3)
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b2)%transform_volume(i,j) = temp * cellsize_py_cube			
		
temp = volume_integration(a2_dp,0.0_dp,0.0_dp,i2s_dp,b2_dp,0.0_dp,0.0_dp,j2s_dp,hb2,k2s_dp,191,1)
pyramid_grid(k)%transform(i,j)%transform_volume(a3,b3) = temp * cellsize_py_cube
cartesian_grid(k)%transform(a3,b3)%transform_volume(i,j) = temp * cellsize_py_cube

endif


        end do
      end do
    end do
		
! the following check if the code calculate transformation correctly, must erase when finish		
sumsum=0.0
sumsum=sumsum+cartesian_grid(1)%transform(1,1)%transform_volume(1,1)

do z=2,level_py
n=real(partition(z))
do x=1,z
do y=1,z
		
min_i=int(n*real(x-1)/real(z))+1
if (n*real(x)/real(z-1).gt.int(n*real(x)/real(z-1))) then
max_i=int(n*real(x)/real(z-1)+1.0)
elseif (n*real(x)/real(z-1).eq.int(n*real(x)/real(z-1))) then
max_i=int(n*real(x)/real(z-1))						
endif

min_j=int(n*real(y-1)/real(z))+1
if (n*real(y)/real(z-1).gt.int(n*real(y)/real(z-1))) then
max_j=int(n*real(y)/real(z-1)+1.0)
elseif (n*real(y)/real(z-1).eq.int(n*real(y)/real(z-1))) then
max_j=int(n*real(y)/real(z-1))						
endif
			
do xx=min_i,max_i
do yy=min_j,max_j
!if (cartesian_grid(z)%transform(x,y)%transform_volume(xx,yy).gt.0.0) then

!sumsum=sumsum+sum(sum(cartesian_grid(z)%transform(x,y)%transform_volume,1),1)
sumsum=sumsum+cartesian_grid(z)%transform(x,y)%transform_volume(xx,yy)
!endif
enddo
enddo
enddo
enddo
enddo		

sumsum2=0.0
do k = 1,level_py
do i = 1,partition(k)
do j = 1,partition(k)
do i_primary = lbound(pyramid_grid(k)%transform(i,j)%transform_volume,1),&
               ubound(pyramid_grid(k)%transform(i,j)%transform_volume,1)
do i_secondary = lbound(pyramid_grid(k)%transform(i,j)%transform_volume,2),&
               ubound(pyramid_grid(k)%transform(i,j)%transform_volume,2)
sumsum2=sumsum2+pyramid_grid(k)%transform(i,j)%transform_volume(i_primary,i_secondary)


enddo
enddo
enddo
enddo
enddo		

write(*,*)'cellsize_py_cube',cellsize_py_cube
write(*,*) 'The cartesian transformation records a volume of ',sumsum
write(*,*) 'The pyramidal transformation records a volume of ',sumsum2
write(*,*) 'The volume of a pyramid is ',level_py*level_py*level_py/3.0 * cellsize_py_cube
write(*,*) 'The relative difference for cartesian is',&
(sumsum-level_py*level_py*level_py/3.0 * cellsize_py_cube)/&
(level_py*level_py*level_py/3.0 * cellsize_py_cube)
write(*,*) 'The relative difference for pyramidal is',&
(sumsum2-level_py*level_py*level_py/3.0 * cellsize_py_cube)/&
(level_py*level_py*level_py/3.0 * cellsize_py_cube)


  end subroutine volume_partition
	
end module geometry_transformation
