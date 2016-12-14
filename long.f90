module long
	
  use precision, only: dp
  use input, only: cellsize_sl,pi,four_pi,sigma_HI_at_ion_freq,sigma_HeI_at_ion_freq,sigma_HeII_at_ion_freq,four_over_three_pi
  use array, only: long_source_HI_density_array, long_source_HeI_density_array, long_source_HeII_density_array, &
                   long_source_HI_column_density_in, long_source_HeI_column_density_in, long_source_HeII_column_density_in, &
                   long_source_HI_column_density_out, long_source_HeI_column_density_out, long_source_HeII_column_density_out, &
                   long_source_shell_volume 
   
contains	
	
  subroutine long_characteristic(a,b,c)

    implicit none
	  
    integer, intent(in) :: a
    integer, intent(in) :: b
    integer, intent(in) :: c
    integer :: i,j,k,l
    integer :: sgn_i,sgn_j,sgn_k
    real(kind=dp) :: t, D
    real(kind=dp) :: t_x
    real(kind=dp) :: t_y
    real(kind=dp) :: t_z
    integer :: abc
    real(kind=dp), dimension(:), allocatable :: t_array
    integer, dimension(:), allocatable :: x_array
    integer, dimension(:), allocatable :: y_array
    integer, dimension(:), allocatable :: z_array
    real(kind=dp), dimension(:), allocatable :: p_array
    real(kind=dp) :: distance_square
  	real(kind=dp) :: in_distance, out_distance	
		
    D = sqrt(real(a*a+b*b+c*c))
    sgn_i = sign(1,a)
    sgn_j = sign(1,b)
    sgn_k = sign(1,c)
    abc = abs(a)+ abs(b) + abs(c) + 1
		
    allocate(t_array(0:abc))
    allocate(x_array(1:abc))
    allocate(y_array(1:abc))
    allocate(z_array(1:abc))
    allocate(p_array(1:abc))					
    t_array(0) = 0

    if (a.ne.0 .and. b.ne.0 .and. c.ne.0) then
      i = 1
      j = 1
      k = 1
      t_x = real((i-0.5)*sgn_i)*D/real(a)
      t_y = real((j-0.5)*sgn_j)*D/real(b)
      t_z = real((k-0.5)*sgn_k)*D/real(c)
      do l = 1,abc
        if (t_x.le.t_y .and. t_x.le.t_z) then		
          t_array(l) = t_x
		  x_array(l) = (i-1)*sgn_i
		  y_array(l) = (j-1)*sgn_j
		  z_array(l) = (k-1)*sgn_k
		  p_array(l) = (t_array(l)-t_array(l-1))*cellsize_sl
          i = i+1
          t_x = real((i-0.5)*sgn_i)*D/real(a)
        elseif (t_y.le.t_x .and. t_y.le.t_z) then		
          t_array(l) = t_y
		  x_array(l) = (i-1)*sgn_i
		  y_array(l) = (j-1)*sgn_j
		  z_array(l) = (k-1)*sgn_k
		  p_array(l) = (t_array(l)-t_array(l-1))*cellsize_sl	  		  
          j = j+1
          t_y = real((j-0.5)*sgn_j)*D/real(b)
        elseif (t_z.le.t_x .and. t_z.le.t_y) then		
          t_array(l) = t_z
		  x_array(l) = (i-1)*sgn_i
		  y_array(l) = (j-1)*sgn_j
		  z_array(l) = (k-1)*sgn_k
		  p_array(l) = (t_array(l)-t_array(l-1))*cellsize_sl	  		  
          k = k+1
          t_z = real((k-0.5)*sgn_k)*D/real(c)
        endif
      enddo
    elseif (a.ne.0 .and. b.ne.0 .and. c.eq.0) then
      i = 1
      j = 1
      t_x = real((i-0.5)*sgn_i)*D/real(a)
      t_y = real((j-0.5)*sgn_j)*D/real(b)
      do l = 1,abc
        if (t_x.le.t_y) then		
          t_array(l) = t_x
		  x_array(l) = (i-1)*sgn_i
		  y_array(l) = (j-1)*sgn_j
		  z_array(l) = 0
		  p_array(l) = (t_array(l)-t_array(l-1))*cellsize_sl
          i = i+1
          t_x = real((i-0.5)*sgn_i)*D/real(a)
        elseif (t_y.le.t_x) then		
          t_array(l) = t_y
		  x_array(l) = (i-1)*sgn_i
		  y_array(l) = (j-1)*sgn_j
		  z_array(l) = 0
		  p_array(l) = (t_array(l)-t_array(l-1))*cellsize_sl	  		  
          j = j+1
          t_y = real((j-0.5)*sgn_j)*D/real(b)
        endif
      enddo

	elseif (a.ne.0 .and. b.eq.0 .and. c.ne.0) then
      i = 1
      k = 1
      t_x = real((i-0.5)*sgn_i)*D/real(a)
      t_z = real((k-0.5)*sgn_k)*D/real(c)
      do l = 1,abc
        if (t_x.le.t_z) then		
          t_array(l) = t_x
		  x_array(l) = (i-1)*sgn_i
		  y_array(l) = 0
		  z_array(l) = (k-1)*sgn_k
		  p_array(l) = (t_array(l)-t_array(l-1))*cellsize_sl
          i = i+1
          t_x = real((i-0.5)*sgn_i)*D/real(a)
        elseif (t_z.le.t_x) then		
          t_array(l) = t_z
		  x_array(l) = (i-1)*sgn_i
		  y_array(l) = 0
		  z_array(l) = (k-1)*sgn_k
		  p_array(l) = (t_array(l)-t_array(l-1))*cellsize_sl	  		  
          k = k+1
          t_z = real((k-0.5)*sgn_k)*D/real(c)
        endif
      enddo
    elseif (a.eq.0 .and. b.ne.0 .and. c.ne.0) then
      j = 1
      k = 1
      t_y = real((j-0.5)*sgn_j)*D/real(b)
      t_z = real((k-0.5)*sgn_k)*D/real(c)
      do l = 1,abc
        if (t_y.le.t_z) then		
          t_array(l) = t_y
          x_array(l) = 0
          y_array(l) = (j-1)*sgn_j
          z_array(l) = (k-1)*sgn_k
          p_array(l) = (t_array(l)-t_array(l-1))*cellsize_sl	  		  
          j = j+1
          t_y = real((j-0.5)*sgn_j)*D/real(b)
        elseif (t_z.le.t_y) then		
          t_array(l) = t_z
          x_array(l) = 0
          y_array(l) = (j-1)*sgn_j
          z_array(l) = (k-1)*sgn_k
          p_array(l) = (t_array(l)-t_array(l-1))*cellsize_sl	  		  
          k = k+1
          t_z = real((k-0.5)*sgn_k)*D/real(c)
        endif
      enddo	
    elseif (a.ne.0 .and. b.eq.0 .and. c.eq.0) then
      i = 1
      t_x = real((i-0.5)*sgn_i)*D/real(a)
      do l = 1,abc		
          t_array(l) = t_x
		  x_array(l) = (i-1)*sgn_i
		  y_array(l) = 0
		  z_array(l) = 0
		  p_array(l) = (t_array(l)-t_array(l-1))*cellsize_sl
          i = i+1
          t_x = real((i-0.5)*sgn_i)*D/real(a)
      enddo
    elseif (a.eq.0 .and. b.ne.0 .and. c.eq.0) then
      j = 1
      t_y = real((j-0.5)*sgn_j)*D/real(b)
      do l = 1,abc		
          t_array(l) = t_y
          x_array(l) = 0
          y_array(l) = (j-1)*sgn_j
          z_array(l) = 0
          p_array(l) = (t_array(l)-t_array(l-1))*cellsize_sl		  		  
          j = j+1
          t_y = real((j-0.5)*sgn_j)*D/real(b)
      enddo	
    elseif (a.eq.0 .and. b.eq.0 .and. c.ne.0) then
      k = 1
      t_z = real((k-0.5)*sgn_k)*D/real(c)
      do l = 1,abc		
          t_array(l) = t_z
          x_array(l) = 0
          y_array(l) = 0
          z_array(l) = (k-1)*sgn_k
          p_array(l) = (t_array(l)-t_array(l-1))*cellsize_sl		  		  
          k = k+1
          t_z = real((k-0.5)*sgn_k)*D/real(c)
      enddo		
	endif	
		
    long_source_HI_column_density_in(a,b,c) = 0
    long_source_HeI_column_density_in(a,b,c) = 0
    long_source_HeII_column_density_in(a,b,c) = 0
				
    do l = 1,abc-1
      long_source_HI_column_density_in(a,b,c) = long_source_HI_column_density_in(a,b,c) + &
               p_array(l)*long_source_HI_density_array(x_array(l),y_array(l),z_array(l))
      long_source_HeI_column_density_in(a,b,c) = long_source_HeI_column_density_in(a,b,c) + &
               p_array(l)*long_source_HeI_density_array(x_array(l),y_array(l),z_array(l))
      long_source_HeII_column_density_in(a,b,c) = long_source_HeII_column_density_in(a,b,c) + &
               p_array(l)*long_source_HeII_density_array(x_array(l),y_array(l),z_array(l))
    enddo
	
    long_source_HI_column_density_out(a,b,c) = long_source_HI_column_density_in(a,b,c)+&
             p_array(abc)*long_source_HI_density_array(x_array(abc),y_array(abc),z_array(abc))
    long_source_HeI_column_density_out(a,b,c) = long_source_HeI_column_density_in(a,b,c)+&
             p_array(abc)*long_source_HeI_density_array(x_array(abc),y_array(abc),z_array(abc))			 
    long_source_HeII_column_density_out(a,b,c) = long_source_HeII_column_density_in(a,b,c)+&
             p_array(abc)*long_source_HeII_density_array(x_array(abc),y_array(abc),z_array(abc))			 

out_distance = sum(p_array)
in_distance = out_distance - p_array(abc)
long_source_shell_volume(a,b,c) = four_over_three_pi*(out_distance*out_distance*out_distance-&
in_distance*in_distance*in_distance)

    deallocate(t_array)				
    deallocate(x_array)	
    deallocate(y_array)	
    deallocate(z_array)	
    deallocate(p_array)	
	
  end subroutine long_characteristic

end module long
