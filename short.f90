module short
	
  use precision, only: dp
  use input, only: cellsize_sl,pi,four_pi,sigma_HI_at_ion_freq,sigma_HeI_at_ion_freq,sigma_HeII_at_ion_freq,four_over_three_pi
  use array, only: short_source_HI_density_array, short_source_HeI_density_array, short_source_HeII_density_array, &
                   short_source_HI_column_density_in, short_source_HeI_column_density_in, short_source_HeII_column_density_in, &
                   short_source_HI_column_density_out, short_source_HeI_column_density_out, short_source_HeII_column_density_out, &
                   short_source_shell_volume
   
contains	
	
  subroutine short_characteristic(i,j,k)

    implicit none
	
    integer, intent(in) :: i
    integer, intent(in) :: j
    integer, intent(in) :: k
    real(kind=dp),parameter :: sqrt3 = sqrt(3.0)
    real(kind=dp),parameter :: sqrt2 = sqrt(2.0)
    integer :: abs_i,abs_j,abs_k
    integer :: im,jm,km
    integer :: ip,imp,jp,jmp,kp,kmp
    integer :: sgn_i,sgn_j,sgn_k
    real(kind=dp) :: ratio,xc,yc,zc,dx,dy,dz,s1,s2,s3,s4
    real(kind=dp) :: c1_HI,c2_HI,c3_HI,c4_HI
    real(kind=dp) :: c1_HeI,c2_HeI,c3_HeI,c4_HeI
    real(kind=dp) :: c1_HeII,c2_HeII,c3_HeII,c4_HeII		
    real(kind=dp) :: dxp,dyp,dzp
    real(kind=dp) :: w1_HI,w2_HI,w3_HI,w4_HI
    real(kind=dp) :: w1_HeI,w2_HeI,w3_HeI,w4_HeI	
    real(kind=dp) :: w1_HeII,w2_HeII,w3_HeII,w4_HeII	
    real(kind=dp) :: di,dj,dk,path,distance_square
    real(kind=dp) :: column_density_HI,column_density_HeI,column_density_HeII
  	real(kind=dp) :: absolute_x, absolute_y, absolute_z
	real(kind=dp) :: in_x, in_y, in_z
	real(kind=dp) :: out_x, out_y, out_z
  	real(kind=dp) :: in_distance, out_distance

    abs_i = abs(i)
    abs_j = abs(j)
    abs_k = abs(k)
    sgn_i = sign(1,i)
    sgn_j = sign(1,j)
    sgn_k = sign(1,k)
    im = i-sgn_i
    jm = j-sgn_j
    km = k-sgn_k
    di = real(i)
    dj = real(j)
    dk = real(k)

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
short_source_shell_volume(i,j,k) = four_over_three_pi*(out_distance*out_distance*out_distance-&
in_distance*in_distance*in_distance)




    if (abs_k .ge. abs_j.and.abs_k .ge. abs_i) then
       
      ratio = (real(km)+sgn_k*0.5)/dk
      xc = ratio*di 
      yc = ratio*dj 
      dx = 2.0*abs(xc-(real(im)+0.5*sgn_i)) ! distances from c-point to
      dy = 2.0*abs(yc-(real(jm)+0.5*sgn_j)) ! the corners.
      s1 = (1.-dx)*(1.-dy)    ! interpolation weights of
      s2 = (1.-dy)*dx         ! corner points to c-point
      s3 = (1.-dx)*dy
      s4 = dx*dy
      c1_HI = short_source_HI_column_density_out(im,jm,km)   
      c2_HI = short_source_HI_column_density_out(i,jm,km)     
      c3_HI = short_source_HI_column_density_out(im,j,km)
      c4_HI = short_source_HI_column_density_out(i,j,km)
      c1_HeI = short_source_HeI_column_density_out(im,jm,km)   
      c2_HeI = short_source_HeI_column_density_out(i,jm,km)     
      c3_HeI = short_source_HeI_column_density_out(im,j,km)
      c4_HeI = short_source_HeI_column_density_out(i,j,km)
      c1_HeII = short_source_HeII_column_density_out(im,jm,km)   
      c2_HeII = short_source_HeII_column_density_out(i,jm,km)     
      c3_HeII = short_source_HeII_column_density_out(im,j,km)
      c4_HeII = short_source_HeII_column_density_out(i,j,km)	  	  
      w1_HI = s1*weight_function(c1_HI,1)
      w2_HI = s2*weight_function(c2_HI,1)
      w3_HI = s3*weight_function(c3_HI,1)
      w4_HI = s4*weight_function(c4_HI,1)
      w1_HeI = s1*weight_function(c1_HeI,2)
      w2_HeI = s2*weight_function(c2_HeI,2)
      w3_HeI = s3*weight_function(c3_HeI,2)
      w4_HeI = s4*weight_function(c4_HeI,2)
      w1_HeII = s1*weight_function(c1_HeII,3)
      w2_HeII = s2*weight_function(c2_HeII,3)
      w3_HeII = s3*weight_function(c3_HeII,3)
      w4_HeII = s4*weight_function(c4_HeII,3)	  	  
      column_density_HI = (c1_HI*w1_HI+c2_HI*w2_HI+c3_HI*w3_HI+c4_HI*w4_HI)/(w1_HI+w2_HI+w3_HI+w4_HI) 
      column_density_HeI = (c1_HeI*w1_HeI+c2_HeI*w2_HeI+c3_HeI*w3_HeI+c4_HeI*w4_HeI)/(w1_HeI+w2_HeI+w3_HeI+w4_HeI)
      column_density_HeII = (c1_HeII*w1_HeII+c2_HeII*w2_HeII+c3_HeII*w3_HeII+c4_HeII*w4_HeII)/(w1_HeII+w2_HeII+w3_HeII+w4_HeII)
	  	  	  
      if (abs_k .eq. 1 .and. (abs_i .eq. 1 .or. abs_j .eq. 1)) then
        if (abs_i .eq. 1 .and. abs_j .eq. 1) then
          column_density_HI = sqrt3*column_density_HI
          column_density_HeI = sqrt3*column_density_HeI
          column_density_HeII = sqrt3*column_density_HeII		  		  
        else if (abs_i .eq. 1 .and. abs_j .eq. 0) then
          column_density_HI = sqrt2*column_density_HI
          column_density_HeI = sqrt2*column_density_HeI
          column_density_HeII = sqrt2*column_density_HeII		  		  
        else if (abs_i .eq. 0 .and. abs_j .eq. 1) then
          column_density_HI = sqrt2*column_density_HI
          column_density_HeI = sqrt2*column_density_HeI
          column_density_HeII = sqrt2*column_density_HeII
        endif
      endif

      short_source_HI_column_density_in(i,j,k) = column_density_HI
      short_source_HeI_column_density_in(i,j,k) = column_density_HeI
      short_source_HeII_column_density_in(i,j,k) = column_density_HeII	  	  
      path = cellsize_sl*sqrt((di*di+dj*dj)/(dk*dk)+1.0) 
      short_source_HI_column_density_out(i,j,k) = short_source_HI_column_density_in(i,j,k) + &	  
                                                  short_source_HI_density_array(i,j,k)*path	
      short_source_HeI_column_density_out(i,j,k) = short_source_HeI_column_density_in(i,j,k) + &	  
                                                  short_source_HeI_density_array(i,j,k)*path	
      short_source_HeII_column_density_out(i,j,k) = short_source_HeII_column_density_in(i,j,k) + &	  
                                                  short_source_HeII_density_array(i,j,k)*path
      !distance_square = cellsize_sl*cellsize_sl*(di*di+dj*dj+dk*dk)						
      !short_source_shell_volume(i,j,k) = four_pi*distance_square*path
	  					  
    elseif (abs_j .ge. abs_i.and.abs_j .ge. abs_k) then
          
      ratio = (real(jm)+sgn_j*0.5)/dj
      zc = ratio*dk
      xc = ratio*di
      dz = 2.0*abs(zc-(real(km)+0.5*sgn_k))
      dx = 2.0*abs(xc-(real(im)+0.5*sgn_i))
      s1 = (1.-dx)*(1.-dz)
      s2 = (1.-dz)*dx
      s3 = (1.-dx)*dz
      s4 = dx*dz
      c1_HI = short_source_HI_column_density_out(im,jm,km)
      c2_HI = short_source_HI_column_density_out(i,jm,km)
      c3_HI = short_source_HI_column_density_out(im,jm,k)
      c4_HI = short_source_HI_column_density_out(i,jm,k)
      c1_HeI = short_source_HeI_column_density_out(im,jm,km)
      c2_HeI = short_source_HeI_column_density_out(i,jm,km)
      c3_HeI = short_source_HeI_column_density_out(im,jm,k)
      c4_HeI = short_source_HeI_column_density_out(i,jm,k)
      c1_HeII = short_source_HeII_column_density_out(im,jm,km)
      c2_HeII = short_source_HeII_column_density_out(i,jm,km)
      c3_HeII = short_source_HeII_column_density_out(im,jm,k)
      c4_HeII = short_source_HeII_column_density_out(i,jm,k)	  	  
      w1_HI = s1*weight_function(c1_HI,1)
      w2_HI = s2*weight_function(c2_HI,1)
      w3_HI = s3*weight_function(c3_HI,1)
      w4_HI = s4*weight_function(c4_HI,1)
      w1_HeI = s1*weight_function(c1_HeI,2)
      w2_HeI = s2*weight_function(c2_HeI,2)
      w3_HeI = s3*weight_function(c3_HeI,2)
      w4_HeI = s4*weight_function(c4_HeI,2)
      w1_HeII = s1*weight_function(c1_HeII,3)
      w2_HeII = s2*weight_function(c2_HeII,3)
      w3_HeII = s3*weight_function(c3_HeII,3)
      w4_HeII = s4*weight_function(c4_HeII,3)	
      column_density_HI = (c1_HI*w1_HI+c2_HI*w2_HI+c3_HI*w3_HI+c4_HI*w4_HI)/(w1_HI+w2_HI+w3_HI+w4_HI) 
      column_density_HeI = (c1_HeI*w1_HeI+c2_HeI*w2_HeI+c3_HeI*w3_HeI+c4_HeI*w4_HeI)/(w1_HeI+w2_HeI+w3_HeI+w4_HeI)
      column_density_HeII = (c1_HeII*w1_HeII+c2_HeII*w2_HeII+c3_HeII*w3_HeII+c4_HeII*w4_HeII)/(w1_HeII+w2_HeII+w3_HeII+w4_HeII)
       
      if (abs_j .eq. 1 .and. (abs_i .eq. 1 .or. abs_k .eq. 1)) then
        if (abs_i .eq. 1 .and. abs_k .eq. 1) then
            column_density_HI = sqrt3*column_density_HI
            column_density_HeI = sqrt3*column_density_HeI
            column_density_HeII = sqrt3*column_density_HeII
        else if (abs_i .eq. 1 .and. abs_k .eq. 0) then
            column_density_HI = sqrt2*column_density_HI
            column_density_HeI = sqrt2*column_density_HeI
            column_density_HeII = sqrt2*column_density_HeII
        else if (abs_i .eq. 0 .and. abs_k .eq. 1) then
            column_density_HI = sqrt2*column_density_HI
            column_density_HeI = sqrt2*column_density_HeI
            column_density_HeII = sqrt2*column_density_HeII
        endif
      endif

      short_source_HI_column_density_in(i,j,k) = column_density_HI
      short_source_HeI_column_density_in(i,j,k) = column_density_HeI
      short_source_HeII_column_density_in(i,j,k) = column_density_HeII	
      path = cellsize_sl*sqrt((di*di+dk*dk)/(dj*dj)+1.0)	  
      short_source_HI_column_density_out(i,j,k) = short_source_HI_column_density_in(i,j,k) + &	  
                                                  short_source_HI_density_array(i,j,k)*path	
      short_source_HeI_column_density_out(i,j,k) = short_source_HeI_column_density_in(i,j,k) + &	  
                                                  short_source_HeI_density_array(i,j,k)*path	
      short_source_HeII_column_density_out(i,j,k) = short_source_HeII_column_density_in(i,j,k) + &	  
                                                  short_source_HeII_density_array(i,j,k)*path
      !distance_square = cellsize_sl*cellsize_sl*(di*di+dj*dj+dk*dk)							
      !short_source_shell_volume(i,j,k) = four_pi*distance_square*path
												  	   
    elseif(abs_i .ge. abs_j .and. abs_i .ge. abs_k) then

      ratio = (real(im)+sgn_i*0.5)/di
      zc = ratio*dk
      yc = ratio*dj
      dz = 2.0*abs(zc-(real(km)+0.5*sgn_k))
      dy = 2.0*abs(yc-(real(jm)+0.5*sgn_j))
      s1 = (1.-dz)*(1.-dy)
      s2 = (1.-dz)*dy
      s3 = (1.-dy)*dz
      s4 = dy*dz
      c1_HI = short_source_HI_column_density_out(im,jm,km)
      c2_HI = short_source_HI_column_density_out(im,j,km)
      c3_HI = short_source_HI_column_density_out(im,jm,k)
      c4_HI = short_source_HI_column_density_out(im,j,k)
      c1_HeI = short_source_HeI_column_density_out(im,jm,km)
      c2_HeI = short_source_HeI_column_density_out(im,j,km)
      c3_HeI = short_source_HeI_column_density_out(im,jm,k)
      c4_HeI = short_source_HeI_column_density_out(im,j,k)
      c1_HeII = short_source_HeII_column_density_out(im,jm,km)
      c2_HeII = short_source_HeII_column_density_out(im,j,km)
      c3_HeII = short_source_HeII_column_density_out(im,jm,k)
      c4_HeII = short_source_HeII_column_density_out(im,j,k)	  	  
      w1_HI = s1*weight_function(c1_HI,1)
      w2_HI = s2*weight_function(c2_HI,1)
      w3_HI = s3*weight_function(c3_HI,1)
      w4_HI = s4*weight_function(c4_HI,1)
      w1_HeI = s1*weight_function(c1_HeI,2)
      w2_HeI = s2*weight_function(c2_HeI,2)
      w3_HeI = s3*weight_function(c3_HeI,2)
      w4_HeI = s4*weight_function(c4_HeI,2)
      w1_HeII = s1*weight_function(c1_HeII,3)
      w2_HeII = s2*weight_function(c2_HeII,3)
      w3_HeII = s3*weight_function(c3_HeII,3)
      w4_HeII = s4*weight_function(c4_HeII,3)	
      column_density_HI = (c1_HI*w1_HI+c2_HI*w2_HI+c3_HI*w3_HI+c4_HI*w4_HI)/(w1_HI+w2_HI+w3_HI+w4_HI) 
      column_density_HeI = (c1_HeI*w1_HeI+c2_HeI*w2_HeI+c3_HeI*w3_HeI+c4_HeI*w4_HeI)/(w1_HeI+w2_HeI+w3_HeI+w4_HeI)
      column_density_HeII = (c1_HeII*w1_HeII+c2_HeII*w2_HeII+c3_HeII*w3_HeII+c4_HeII*w4_HeII)/(w1_HeII+w2_HeII+w3_HeII+w4_HeII)
       
      if (abs_i .eq. 1 .and. (abs_j .eq. 1 .or. abs_k .eq. 1)) then
        if (abs_j .eq. 1 .and. abs_k .eq. 1) then
            column_density_HI = sqrt3*column_density_HI
            column_density_HeI = sqrt3*column_density_HeI
            column_density_HeII = sqrt3*column_density_HeII
        else if (abs_j .eq. 1 .and. abs_k .eq. 0) then
            column_density_HI = sqrt2*column_density_HI
            column_density_HeI = sqrt2*column_density_HeI
            column_density_HeII = sqrt2*column_density_HeII
        else if (abs_j .eq. 0 .and. abs_k .eq. 1) then
            column_density_HI = sqrt2*column_density_HI
            column_density_HeI = sqrt2*column_density_HeI
            column_density_HeII = sqrt2*column_density_HeII
        endif
      endif

      short_source_HI_column_density_in(i,j,k) = column_density_HI
      short_source_HeI_column_density_in(i,j,k) = column_density_HeI
      short_source_HeII_column_density_in(i,j,k) = column_density_HeII
      path = cellsize_sl*sqrt(1.0+(dj*dj+dk*dk)/(di*di))  
      short_source_HI_column_density_out(i,j,k) = short_source_HI_column_density_in(i,j,k) + &	  
                                                  short_source_HI_density_array(i,j,k)*path	
      short_source_HeI_column_density_out(i,j,k) = short_source_HeI_column_density_in(i,j,k) + &	  
                                                  short_source_HeI_density_array(i,j,k)*path	
      short_source_HeII_column_density_out(i,j,k) = short_source_HeII_column_density_in(i,j,k) + &	  
                                                  short_source_HeII_density_array(i,j,k)*path
      !distance_square = cellsize_sl*cellsize_sl*(di*di+dj*dj+dk*dk)
      !short_source_shell_volume(i,j,k) = four_pi*distance_square*path
												  	  
    end if
    
  end subroutine short_characteristic

  real(kind=dp) function weight_function (cd,id)
  
    implicit none
  
    real(kind=dp),intent(in) :: cd
    integer,intent(in) :: id
	real(kind=dp) :: sigma
    real(kind=dp),parameter :: minweight = 1.0_dp/0.6_dp

    if (id .eq. 1) then
		sigma = sigma_HI_at_ion_freq
		elseif (id .eq. 2) then
		sigma = sigma_HeI_at_ion_freq	
		elseif (id .eq. 3) then
		sigma = sigma_HeII_at_ion_freq	
	endif
    weight_function = 1.0/max(0.6_dp,cd*sigma)
 
  end function weight_function
	
end module short
