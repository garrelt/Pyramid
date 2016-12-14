module transformation
	
  use precision, only: dp
  use input, only: level, partition
  use initialization, only: pyramid_grid, cartesian_grid
  use integration, only: vol_integration

contains
	
  subroutine volume_division

    implicit none
  
    integer :: i,j,k,n
		real :: i1,i2,i3,i4
		real :: j1,j2,j3,j4
		integer :: min_line_i, max_line_i
		integer :: min_line_j, max_line_j
		integer :: partition_i, partition_j
		integer :: a1,a2,a3
		integer :: b1,b2,b3
		real :: ha1,ha2
		real :: hb1,hb2
		integer :: count1, count2
		real :: tv
		real ::sumsum, sumsum2
		integer ::x,y,z
		integer :: level2
		integer :: xx,yy,zz
		integer :: min_i,max_i,min_j,max_j
		
		do k=1,level
		  do i=1,partition(k)
		    do j=1,partition(k)
				partition_i=0
				partition_j=0
				
			n=real(partition(k))
			i1=real(i-1)*real(k-1)/n
			i2=real(i-1)*real(k)/n
			i3=real(i)*real(k-1)/n
			i4=real(i)*real(k)/n
			j1=real(j-1)*real(k-1)/n
			j2=real(j-1)*real(k)/n
			j3=real(j)*real(k-1)/n
			j4=real(j)*real(k)/n
			
			
if (int(i4)-int(i1).eq.0) then
	a1=int(i1)+1
elseif (int(i4)-int(i1).eq.1 .and. i4-int(i4).le.0.00000001) then	
	a1=int(i4)	
elseif (int(i4)-int(i1).eq.1 .and. i4-int(i4).gt.0.00000001) then
	a1=int(i4)
	a2=int(i4)+1
	if (i1.le.a1 .and. a1.le.i2) then
		ha1=real(n)*real(a1)/real(i-1)
	elseif (i3.le.a1 .and. a1.le.i4) then
		ha1=real(n)*real(a1)/real(i)	
	endif
elseif (int(i4)-int(i1).eq.2 .and. i4-int(i4).le.0.00000001) then
	a1=int(i1)+1
	a2=int(i4)
	if (i1.le.a1 .and. a1.le.i2) then
		ha1=real(n)*real(a1)/real(i-1)
	elseif (i3.le.a1 .and. a1.le.i4) then
		ha1=real(n)*real(a1)/real(i)	
	endif
elseif (int(i4)-int(i1).eq.2 .and. i4-int(i4).gt.0.00000001) then
	a1=int(i1)+1
	a2=int(i4)
	a3=int(i4)+1	
	ha1=real(n)*real(a1)/real(i-1)
	ha2=real(n)*real(a2)/real(i)	
endif


if (int(j4)-int(j1).eq.0) then
	b1=int(j1)+1
elseif (int(j4)-int(j1).eq.1 .and. j4-int(j4).le.0.00000001) then	
	b1=int(j4)	
elseif (int(j4)-int(j1).eq.1 .and. j4-int(j4).gt.0.00000001) then
	b1=int(j4)
	b2=int(j4)+1
	if (j1.le.b1 .and. b1.le.j2) then
		hb1=real(n)*real(b1)/real(j-1)
	elseif (j3.le.b1 .and. b1.le.j4) then
		hb1=real(n)*real(b1)/real(j)	
	endif
elseif (int(j4)-int(j1).eq.2 .and. j4-int(j4).le.0.00000001) then
	b1=int(j1)+1
	b2=int(j4)
	if (j1.le.b1 .and. b1.le.j2) then
		hb1=real(n)*real(b1)/real(j-1)
	elseif (j3.le.b1 .and. b1.le.j4) then
		hb1=real(n)*real(b1)/real(j)
	endif	
elseif (int(j4)-int(j1).eq.2 .and. j4-int(j4).gt.0.00000001) then
	b1=int(j1)+1
	b2=int(j4)
	b3=int(j4)+1	
	hb1=real(n)*real(b1)/real(j-1)
	hb2=real(n)*real(b2)/real(j)	
endif




if (pyramid_grid(k)%case(i,j).eq.1) then
		
	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,real(k-1),real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv

elseif (pyramid_grid(k)%case(i,j).eq.2) then

	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,0.0,real(j)/n,real(k-1),ha1)	
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv	

	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,ha1,real(k))+&				
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv	
		
elseif (pyramid_grid(k)%case(i,j).eq.3) then
		
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,0.0,real(j)/n,real(k-1),real(k))		
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv		
			
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,real(k-1),real(k))		
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv		

elseif (pyramid_grid(k)%case(i,j).eq.4) then
	
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,0.0,real(j)/n,ha1,real(k))+&		
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv		
				
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,ha1,real(k))	
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv	
	 	
elseif (pyramid_grid(k)%case(i,j).eq.5) then
	
	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv			
		
	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,hb1,real(k))+&		
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),hb1)	
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv		
	
elseif (pyramid_grid(k)%case(i,j).eq.6) then

	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv	
			
	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv			
		
elseif (pyramid_grid(k)%case(i,j).eq.7) then

	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,hb1,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv		
			
	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,hb1,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv	
		
elseif (pyramid_grid(k)%case(i,j).eq.8) then
				
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv		
								
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,0.0,real(j)/n,hb1,ha1)+&				
	vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv		
			
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv		
			
	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,ha1,real(k))+&				
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,hb1,ha1)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),hb1)	
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv	
	cartesian_grid(a2,b2,k)%transform(i,j)=tv		
	
elseif (pyramid_grid(k)%case(i,j).eq.9) then
				
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv	
							
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv	
				
	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,ha1,hb1)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv	
				
	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,hb1,real(k))+&				
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,ha1,hb1)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),ha1)		
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv	

elseif (pyramid_grid(k)%case(i,j).eq.10) then
				
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv	
												
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,0.0,real(j)/n,hb1,real(k))+&
	vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv	
			
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv	
				
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,hb1,real(k))+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),hb1)	
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv	

elseif (pyramid_grid(k)%case(i,j).eq.11) then
						
	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv	
							
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,0.0,real(j)/n,ha1,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,hb1,ha1)+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv	
	
	tv=0.0
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv
			
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,ha1,real(k))		
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv	

elseif (pyramid_grid(k)%case(i,j).eq.12) then
					
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,ha1,hb1)+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv	
								
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,0.0,real(j)/n,hb1,real(k))+&
	vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,ha1,hb1)+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv	
				
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,ha1,hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv	
				
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,hb1,real(k))+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,ha1,hb1)	
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv	
	
elseif (pyramid_grid(k)%case(i,j).eq.13) then
					
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv		
							
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv		
			
	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,ha1,real(k))+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv		
			
	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,ha1,real(k))+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv		

elseif (pyramid_grid(k)%case(i,j).eq.14) then
			
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,real(k-1),real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv	
	cartesian_grid(a1,b1,k)%transform(i,j)=tv						
							
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,real(k-1),real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv		

	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv		
			
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv	

elseif (pyramid_grid(k)%case(i,j).eq.15) then
					
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,ha1,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv				
								
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,ha1,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv		

	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,ha1,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv		
			
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,ha1,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv
	
elseif (pyramid_grid(k)%case(i,j).eq.16) then
				
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,hb1,ha1)+&
	vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv	
	cartesian_grid(a1,b1,k)%transform(i,j)=tv				
							
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,hb1,ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv	

	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,ha1,real(k))+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,hb1,ha1)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv			

	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,ha1,real(k))+&				
	vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,hb1,ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv	
	
elseif (pyramid_grid(k)%case(i,j).eq.17) then
			
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv
		
	tv=0.0
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv		

	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,hb1,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,ha1,hb1)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv				

	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,hb1,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv		

elseif (pyramid_grid(k)%case(i,j).eq.18) then
				
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,hb1,real(k))+&
	vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv				
							
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,hb1,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv	

	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,hb1,real(k))+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv		
				
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,hb1,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv
	
elseif (pyramid_grid(k)%case(i,j).eq.19) then
						
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,ha1,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,hb1,ha1)+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv				
							
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,ha1,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,hb1,ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv			

	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,ha1,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv			
				
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,ha1,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv	

elseif (pyramid_grid(k)%case(i,j).eq.20) then
				
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,hb1,real(k))+&
	vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,0.0,real(j)/n,ha1,hb1)+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv					
							
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,hb1,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv			

	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,hb1,real(k))+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,ha1,hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv			
				
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,hb1,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv	
	
elseif (pyramid_grid(k)%case(i,j).eq.21) then
											
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv			
			
	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,0.0,real(j-1)/n,0.0,real(j)/n,ha1,real(k))+&
	vol_integration(real(a1),0.0,real(a2),0.0,0.0,real(j-1)/n,0.0,real(j)/n,ha2,ha1)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,real(k-1),ha2)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv	
	cartesian_grid(a2,b1,k)%transform(i,j)=tv			
				
	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,ha2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1)=tv
	cartesian_grid(a3,b1,k)%transform(i,j)=tv		
	
	
elseif (pyramid_grid(k)%case(i,j).eq.22) then
					
	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv	
	cartesian_grid(a1,b1,k)%transform(i,j)=tv			
			
	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b2),0.0,hb1,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b1),0.0,real(b2),0.0,hb2,hb1)+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),hb2)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv				
						
	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b2),0.0,0.0,real(j)/n,hb2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3)=tv	
	cartesian_grid(a1,b3,k)%transform(i,j)=tv		

elseif (pyramid_grid(k)%case(i,j).eq.23) then

	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,real(k-1),hb1)	
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv						
							
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,0.0,real(j)/n,hb1,ha1)+&
	vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv			
			
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv	
	cartesian_grid(a2,b1,k)%transform(i,j)=tv			
				
	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,0.0,real(j-1)/n,0.0,real(j)/n,ha2,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,ha1,ha2)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,hb1,ha1)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),hb1)	
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv	
			
	tv=0.0
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1)=tv
	cartesian_grid(a3,b1,k)%transform(i,j)=tv		
					
	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,ha2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2)=tv
	cartesian_grid(a3,b2,k)%transform(i,j)=tv		
			
elseif (pyramid_grid(k)%case(i,j).eq.24) then

	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,real(k-1),ha1)	
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv					
		
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv		
				
	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,ha1,hb1)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv		
				
	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,0.0,real(j-1)/n,0.0,real(j)/n,ha2,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,hb1,ha2)+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,ha1,hb1)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),ha1)	
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv	
		
	tv=0.0
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1)=tv
	cartesian_grid(a3,b1,k)%transform(i,j)=tv	
					
	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,ha2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2)=tv
	cartesian_grid(a3,b2,k)%transform(i,j)=tv	

elseif (pyramid_grid(k)%case(i,j).eq.25) then
	
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,real(k-1),ha1)	
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv					
			
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv			
					
	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,0.0,real(j-1)/n,real(b1),0.0,ha2,hb1)+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,ha1,ha2)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv			
				
	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,0.0,real(j-1)/n,0.0,real(j)/n,hb1,real(k))+&
	vol_integration(0.0,real(i-1)/n,real(a2),0.0,real(b1),0.0,0.0,real(j)/n,ha2,hb1)+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,ha1,ha2)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),ha1)	
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv			
					
	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,ha2,hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1)=tv
	cartesian_grid(a3,b1,k)%transform(i,j)=tv		
				
	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,hb1,real(k))+&
	vol_integration(real(a2),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,ha2,hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2)=tv
	cartesian_grid(a3,b2,k)%transform(i,j)=tv		

elseif (pyramid_grid(k)%case(i,j).eq.26) then
	
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,real(k-1),ha1)	
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv				
		
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv		
				
	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,0.0,real(j-1)/n,real(b1),0.0,ha2,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,ha1,ha2)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv		
			
	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,real(b1),0.0,0.0,real(j)/n,ha2,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,ha1,ha2)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),ha1)	
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv	
	cartesian_grid(a2,b2,k)%transform(i,j)=tv		
					
	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,ha2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1)=tv
	cartesian_grid(a3,b1,k)%transform(i,j)=tv	
			
	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,ha2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2)=tv	
	cartesian_grid(a3,b2,k)%transform(i,j)=tv	

elseif (pyramid_grid(k)%case(i,j).eq.27) then
	
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,hb1,ha1)+&	
	vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv				
		
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,hb1,ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv		
				
	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,0.0,real(j-1)/n,real(b1),0.0,ha2,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,ha1,ha2)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,hb1,ha1)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv			
				
	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,real(b1),0.0,0.0,real(j)/n,ha2,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,ha1,ha2)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,hb1,ha1)	
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv		
					
	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,ha2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1)=tv
	cartesian_grid(a3,b1,k)%transform(i,j)=tv		
				
	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,ha2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2)=tv
	cartesian_grid(a3,b2,k)%transform(i,j)=tv	

elseif (pyramid_grid(k)%case(i,j).eq.28) then
					
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv
			
	tv=0.0
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv			
					
	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,0.0,real(j-1)/n,real(b1),0.0,ha2,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,hb1,ha2)+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,ha1,hb1)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv		
				
	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,real(b1),0.0,0.0,real(j)/n,ha2,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,hb1,ha2)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv		
					
	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,ha2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1)=tv
	cartesian_grid(a3,b1,k)%transform(i,j)=tv		
					
	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,ha2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2)=tv
	cartesian_grid(a3,b2,k)%transform(i,j)=tv	

elseif (pyramid_grid(k)%case(i,j).eq.29) then
						
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv
				
	tv=0.0
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv		
					
	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,0.0,real(j-1)/n,real(b1),0.0,hb1,real(k))+&
	vol_integration(0.0,real(i-1)/n,real(a2),0.0,0.0,real(j-1)/n,0.0,real(j)/n,ha2,hb1)+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,ha1,ha2)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv		
				
	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,real(b1),0.0,0.0,real(j)/n,hb1,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv		
				
	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,hb1,real(k))+&
	vol_integration(real(a2),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,ha2,hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1)=tv
	cartesian_grid(a3,b1,k)%transform(i,j)=tv		
				
	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,hb1,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2)=tv
	cartesian_grid(a3,b2,k)%transform(i,j)=tv	
	
elseif (pyramid_grid(k)%case(i,j).eq.30) then
						
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv	
			
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv	
	cartesian_grid(a1,b2,k)%transform(i,j)=tv	
	
	tv=0.0	
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3)=tv
	cartesian_grid(a1,b3,k)%transform(i,j)=tv		
				
	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,ha1,hb1)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv	
			
	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b2),0.0,hb2,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,hb1,hb2)+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,ha1,hb1)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv	
	cartesian_grid(a2,b2,k)%transform(i,j)=tv						
		
	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b2),0.0,0.0,real(j)/n,hb2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3)=tv
	cartesian_grid(a2,b3,k)%transform(i,j)=tv	
	
elseif (pyramid_grid(k)%case(i,j).eq.31) then
							
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv	
				
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,0.0,real(j)/n,hb1,ha1)+&
	vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv	

	tv=0.0
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3)=tv
	cartesian_grid(a1,b3,k)%transform(i,j)=tv		
		
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv	
	cartesian_grid(a2,b1,k)%transform(i,j)=tv		
				
	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b2),0.0,hb2,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,ha1,hb2)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,hb1,ha1)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv				
		
	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b2),0.0,0.0,real(j)/n,hb2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3)=tv
	cartesian_grid(a2,b3,k)%transform(i,j)=tv	
	
elseif (pyramid_grid(k)%case(i,j).eq.32) then
						
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv	
	cartesian_grid(a1,b1,k)%transform(i,j)=tv	
				
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b2),0.0,hb2,ha1)+&
	vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,0.0,real(j)/n,hb1,hb2)+&
	vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv					
			
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b2),0.0,0.0,real(j)/n,hb2,ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3)=tv
	cartesian_grid(a1,b3,k)%transform(i,j)=tv	
				
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv	
					
	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b2),0.0,ha1,real(k))+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b2),0.0,hb2,ha1)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,hb1,hb2)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv				
		
	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b2),0.0,0.0,real(j)/n,ha1,real(k))+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,real(b2),0.0,0.0,real(j)/n,hb2,ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3)=tv
	cartesian_grid(a2,b3,k)%transform(i,j)=tv		
	
elseif (pyramid_grid(k)%case(i,j).eq.33) then
								
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv		
				
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b2),0.0,hb2,real(k))+&
	vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,0.0,real(j)/n,hb1,hb2)+&
	vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv			
			
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b2),0.0,0.0,real(j)/n,hb2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3)=tv
	cartesian_grid(a1,b3,k)%transform(i,j)=tv	
			
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv	
				
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b2),0.0,hb2,real(k))+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,hb1,hb2)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv				
	
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,real(b2),0.0,0.0,real(j)/n,hb2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3)=tv	
	cartesian_grid(a2,b3,k)%transform(i,j)=tv		
	
elseif (pyramid_grid(k)%case(i,j).eq.34) then
					
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,ha1,hb1)+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv	
	cartesian_grid(a1,b1,k)%transform(i,j)=tv			
					
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b2),0.0,hb2,real(k))+&
	vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,0.0,real(j)/n,hb1,hb2)+&
	vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,ha1,hb1)+&	
	vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv					
	
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b2),0.0,0.0,real(j)/n,hb2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3)=tv	
	cartesian_grid(a1,b3,k)%transform(i,j)=tv		
				
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,ha1,hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv	
				
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b2),0.0,hb2,real(k))+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,hb1,hb2)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,ha1,hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv			
		
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,real(b2),0.0,0.0,real(j)/n,hb2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3)=tv
	cartesian_grid(a2,b3,k)%transform(i,j)=tv	
	
elseif (pyramid_grid(k)%case(i,j).eq.35) then
				
	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv	
	cartesian_grid(a1,b1,k)%transform(i,j)=tv		
					
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b2),0.0,hb2,real(k))+&
	vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,0.0,real(j)/n,ha1,hb2)+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,hb1,ha1)+&	
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv						
		
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b2),0.0,0.0,real(j)/n,hb2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3)=tv
	cartesian_grid(a1,b3,k)%transform(i,j)=tv		

	tv=0.0	
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv			
					
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b2),0.0,hb2,real(k))+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,ha1,hb2)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv			
		
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,real(b2),0.0,0.0,real(j)/n,hb2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3)=tv
	cartesian_grid(a2,b3,k)%transform(i,j)=tv	
	
elseif (pyramid_grid(k)%case(i,j).eq.36) then
				
	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv				
				
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b2),0.0,ha1,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b2),0.0,hb2,ha1)+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,hb1,hb2)+&	
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv					
		
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b2),0.0,0.0,real(j)/n,ha1,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b2),0.0,0.0,real(j)/n,hb2,ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3)=tv
	cartesian_grid(a1,b3,k)%transform(i,j)=tv			

	tv=0.0
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv			
					
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b2),0.0,ha1,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv				
		
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,real(b2),0.0,0.0,real(j)/n,ha1,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3)=tv
	cartesian_grid(a2,b3,k)%transform(i,j)=tv		
	
elseif (pyramid_grid(k)%case(i,j).eq.37) then
			
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv	

	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b2),0.0,hb2,ha1)+&
	vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,0.0,real(j)/n,hb1,hb2)+&			
	vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv						
		
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b2),0.0,0.0,real(j)/n,hb2,ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3)=tv
	cartesian_grid(a1,b3,k)%transform(i,j)=tv	
				
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv	

	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,0.0,real(j-1)/n,real(b2),0.0,ha2,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b2),0.0,ha1,ha2)+&			
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b2),0.0,hb2,ha1)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,hb1,hb2)+&		
	vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv				

	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,real(b2),0.0,0.0,real(j)/n,ha2,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b2),0.0,0.0,real(j)/n,ha1,ha2)+&			
	vol_integration(real(a1),0.0,0.0,real(i)/n,real(b2),0.0,0.0,real(j)/n,hb2,ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3)=tv
	cartesian_grid(a2,b3,k)%transform(i,j)=tv	

	tv=0.0	
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1)=tv
	cartesian_grid(a3,b1,k)%transform(i,j)=tv				

	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b2),0.0,ha2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2)=tv
	cartesian_grid(a3,b2,k)%transform(i,j)=tv						
			
	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,real(b2),0.0,0.0,real(j)/n,ha2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b3)=tv
	cartesian_grid(a3,b3,k)%transform(i,j)=tv		
	
elseif (pyramid_grid(k)%case(i,j).eq.38) then
			
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv		

	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,0.0,real(j)/n,hb1,ha1)+&			
	vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv							

	tv=0.0
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3)=tv
	cartesian_grid(a1,b3,k)%transform(i,j)=tv		
				
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv		

	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,0.0,real(j-1)/n,real(b2),0.0,ha2,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b2),0.0,hb2,ha2)+&			
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,ha1,hb2)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,hb1,ha1)+&		
	vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv			
	
	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,real(b2),0.0,0.0,real(j)/n,ha2,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b2),0.0,0.0,real(j)/n,hb2,ha2)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3)=tv
	cartesian_grid(a2,b3,k)%transform(i,j)=tv	
	
	tv=0.0
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1)=tv
	cartesian_grid(a3,b1,k)%transform(i,j)=tv		

	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b2),0.0,ha2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2)=tv
	cartesian_grid(a3,b2,k)%transform(i,j)=tv					
			
	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,real(b2),0.0,0.0,real(j)/n,ha2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b3)=tv
	cartesian_grid(a3,b3,k)%transform(i,j)=tv		
	
elseif (pyramid_grid(k)%case(i,j).eq.39) then
					
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv			

	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,0.0,real(j)/n,hb1,ha1)+&			
	vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv					
	
	tv=0.0
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3)=tv
	cartesian_grid(a1,b3,k)%transform(i,j)=tv			
				
	tv=vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv			

	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,0.0,real(j-1)/n,real(b2),0.0,hb2,real(k))+&
	vol_integration(0.0,real(i-1)/n,real(a2),0.0,0.0,real(j-1)/n,0.0,real(j)/n,ha2,hb2)+&			
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,ha1,ha2)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,hb1,ha1)+&		
	vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv		
	
	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,real(b2),0.0,0.0,real(j)/n,hb2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3)=tv
	cartesian_grid(a2,b3,k)%transform(i,j)=tv			
	
	tv=0.0
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1)=tv
	cartesian_grid(a3,b1,k)%transform(i,j)=tv			

	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b2),0.0,hb2,real(k))+&
	vol_integration(real(a2),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,ha2,hb2)
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2)=tv
	cartesian_grid(a3,b2,k)%transform(i,j)=tv				
		
	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,real(b2),0.0,0.0,real(j)/n,hb2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b3)=tv
	cartesian_grid(a3,b3,k)%transform(i,j)=tv		
	
elseif (pyramid_grid(k)%case(i,j).eq.40) then
				
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv	
	cartesian_grid(a1,b1,k)%transform(i,j)=tv			
				
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv				
	
	tv=0.0
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3)=tv
	cartesian_grid(a1,b3,k)%transform(i,j)=tv			
	
	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,ha1,hb1)+&					
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv		

	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,0.0,real(j-1)/n,real(b2),0.0,hb2,real(k))+&
	vol_integration(0.0,real(i-1)/n,real(a2),0.0,0.0,real(j-1)/n,0.0,real(j)/n,ha2,hb2)+&			
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,hb1,ha2)+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,ha1,hb1)+&		
	vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv		

	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,real(b2),0.0,0.0,real(j)/n,hb2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3)=tv
	cartesian_grid(a2,b3,k)%transform(i,j)=tv		
	
	tv=0.0
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1)=tv
	cartesian_grid(a3,b1,k)%transform(i,j)=tv		
	
	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b2),0.0,hb2,real(k))+&
	vol_integration(real(a2),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,ha2,hb2)
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2)=tv
	cartesian_grid(a3,b2,k)%transform(i,j)=tv				
		
	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,real(b2),0.0,0.0,real(j)/n,hb2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b3)=tv
	cartesian_grid(a3,b3,k)%transform(i,j)=tv		
	
elseif (pyramid_grid(k)%case(i,j).eq.41) then
						
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv		
				
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv					

	tv=0.0
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3)=tv
	cartesian_grid(a1,b3,k)%transform(i,j)=tv		

	tv=vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,ha1,hb1)+&					
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv	

	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,0.0,real(j-1)/n,real(b2),0.0,ha2,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b2),0.0,hb2,ha2)+&			
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,hb1,hb2)+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,ha1,hb1)+&		
	vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv				

	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,real(b2),0.0,0.0,real(j)/n,ha2,real(k))+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b2),0.0,0.0,real(j)/n,hb2,ha2)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3)=tv
	cartesian_grid(a2,b3,k)%transform(i,j)=tv	

	tv=0.0
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1)=tv
	cartesian_grid(a3,b1,k)%transform(i,j)=tv		

	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b2),0.0,ha2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2)=tv
	cartesian_grid(a3,b2,k)%transform(i,j)=tv							
		
	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,real(b2),0.0,0.0,real(j)/n,ha2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b3)=tv
	cartesian_grid(a3,b3,k)%transform(i,j)=tv			
	
elseif (pyramid_grid(k)%case(i,j).eq.42) then
				
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,0.0,real(j-1)/n,real(b1),0.0,real(k-1),ha1)		
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b1)=tv
	cartesian_grid(a1,b1,k)%transform(i,j)=tv			
				
	tv=vol_integration(0.0,real(i-1)/n,real(a1),0.0,real(b1),0.0,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b2)=tv
	cartesian_grid(a1,b2,k)%transform(i,j)=tv			

	tv=0.0
	pyramid_grid(k)%transform(i,j)%transform_volume(a1,b3)=tv
	cartesian_grid(a1,b3,k)%transform(i,j)=tv			
	
	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,0.0,real(j-1)/n,real(b1),0.0,ha2,hb1)+&					
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,ha1,ha2)+&
	vol_integration(real(a1),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b1)=tv
	cartesian_grid(a2,b1,k)%transform(i,j)=tv		

	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,0.0,real(j-1)/n,real(b2),0.0,hb2,real(k))+&
	vol_integration(0.0,real(i-1)/n,real(a2),0.0,0.0,real(j-1)/n,0.0,real(j)/n,hb1,hb2)+&			
	vol_integration(0.0,real(i-1)/n,real(a2),0.0,real(b1),0.0,0.0,real(j)/n,ha2,hb1)+&
	vol_integration(0.0,real(i-1)/n,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,ha1,ha2)+&		
	vol_integration(real(a1),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,real(k-1),ha1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b2)=tv
	cartesian_grid(a2,b2,k)%transform(i,j)=tv			

	tv=vol_integration(0.0,real(i-1)/n,real(a2),0.0,real(b2),0.0,0.0,real(j)/n,hb2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a2,b3)=tv
	cartesian_grid(a2,b3,k)%transform(i,j)=tv			
			
	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b1),0.0,ha2,hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b1)=tv
	cartesian_grid(a3,b1,k)%transform(i,j)=tv			
	
	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,0.0,real(j-1)/n,real(b2),0.0,hb2,real(k))+&
	vol_integration(real(a2),0.0,0.0,real(i)/n,0.0,real(j-1)/n,0.0,real(j)/n,hb1,hb2)+&
	vol_integration(real(a2),0.0,0.0,real(i)/n,real(b1),0.0,0.0,real(j)/n,ha2,hb1)
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b2)=tv
	cartesian_grid(a3,b2,k)%transform(i,j)=tv			
		
	tv=vol_integration(real(a2),0.0,0.0,real(i)/n,real(b2),0.0,0.0,real(j)/n,hb2,real(k))
	pyramid_grid(k)%transform(i,j)%transform_volume(a3,b3)=tv
	cartesian_grid(a3,b3,k)%transform(i,j)=tv		
																									
endif

		    end do
		  end do
		end do
		
		
sumsum=0.0
sumsum=sumsum+cartesian_grid(1,1,1)%transform(1,1)
do z=2,level
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
sumsum=sumsum+cartesian_grid(x,y,z)%transform(xx,yy)
enddo
enddo
enddo
enddo
enddo		
write(*,*) sumsum, level*level*level/3.0, &
real(sumsum-level*level*level/3.0)/real(level*level*level/3.0)
		
	!sumsum2=0	
!		do k=1,level2
!		  do i=1,partition(k)
!		    do j=1,partition(k)
!		sumsum2=sumsum2+sum(pyramid_grid(k)%transform(i,j)%transform_volume)
				
!				enddo
!				enddo
!				enddo
			!write(*,*) sumsum2, level2*level2*level2/3.0			
	end subroutine volume_division
	
end module transformation