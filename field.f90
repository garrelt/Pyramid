module field
	
  use precision, only: dp
  use input, only: level_py, level_sl, number_density, xHI, xHeI, xHeII, cellsize_py_cube, &
                   use_which_field
  use array, only: pyramid_global_number_density_array, pyramid_global_xHI_array, &
                   pyramid_global_xHeI_array, pyramid_global_xHeII_array, &  
                   short_global_number_density_array, short_global_xHI_array, &
                   short_global_xHeI_array, short_global_xHeII_array, &				   
                   long_global_number_density_array, long_global_xHI_array, &
                   long_global_xHeI_array, long_global_xHeII_array, &
                   pyramid_analytical_number_density, pyramid_analytical_xHI, & 
                   pyramid_analytical_xHeI, pyramid_analytical_xHeII, &  
                   cartesian_analytical_number_density, cartesian_analytical_xHI, & 
                   cartesian_analytical_xHeI, cartesian_analytical_xHeII
				   				   				
contains	
	
  subroutine pyramid_field_generation()
	  
    implicit none	  
    integer :: i,j,k
    real(kind=dp) :: r2
    
  if (use_which_field .eq. 'C') then
    do i = 1,2*level_py  
      do j = 1,2*level_py
        do k = 1,2*level_py
r2 = (i-level_py-0.5)**2 + (j-level_py-0.5)**2 + (k-level_py-0.5)**2
            !pyramid_global_number_density_array(i,j,k) = number_density	
	pyramid_global_number_density_array(i,j,k) = number_density! / r2			 
        enddo
      enddo
    enddo
  elseif (use_which_field .eq. 'O') then
    do i = 1,2*level_py  
      do j = 1,2*level_py
        do k = 1,2*level_py

            pyramid_global_number_density_array(i,j,k) = number_density	
					 
        enddo
      enddo
    enddo
! obstacle 
! cp pyramid_3D_HI_photoionization.dat obstacle_pyramid_HI_photoionization.dat
! cp short_3D_HI_photoionization.dat obstacle_short_HI_photoionization.dat
! cp long_3D_HI_photoionization.dat obstacle_long_HI_photoionization.dat
!pyramid_global_number_density_array((level_py*4)/3,level_py+10,level_py) = number_density*1e3
!pyramid_global_number_density_array((level_py*4)/3+1,level_py+10,level_py) = number_density*1e3
!pyramid_global_number_density_array((level_py*4)/3,level_py+11,level_py) = number_density*1e3 
!pyramid_global_number_density_array((level_py*4)/3+1,level_py+11,level_py) = number_density*1e3
!pyramid_global_number_density_array((level_py*2)/3,level_py+10,level_py) = number_density*1e3
!pyramid_global_number_density_array((level_py*2)/3+7,level_py-30,level_py) = number_density*1e3
!pyramid_global_number_density_array((level_py*2)/3+8,level_py-30,level_py) = number_density*1e3
!pyramid_global_number_density_array((level_py*2)/3+9,level_py-30,level_py) = number_density*1e3
!pyramid_global_number_density_array((level_py*2)/3+10,level_py-30,level_py) = number_density*1e3
!pyramid_global_number_density_array((level_py*2)/3+11,level_py-30,level_py) = number_density*1e3
!pyramid_global_number_density_array((level_py*2)/3+12,level_py-30,level_py) = number_density*1e3
!pyramid_global_number_density_array((level_py*2)/3+13,level_py-30,level_py) = number_density*1e3
!pyramid_global_number_density_array((level_py*2)/3+14,level_py-30,level_py) = number_density*1e3
!pyramid_global_number_density_array((level_py*2)/3+15,level_py-30,level_py) = number_density*1e3
!pyramid_global_number_density_array((level_py*2)/3+16,level_py-30,level_py) = number_density*1e3
!pyramid_global_number_density_array((level_py*2)/3+17,level_py-30,level_py) = number_density*1e3
!pyramid_global_number_density_array((level_py*2)/3+18,level_py-30,level_py) = number_density*1e3
!pyramid_global_number_density_array((level_py*2)/3+19,level_py-30,level_py) = number_density*1e3
!pyramid_global_number_density_array((level_py*2)/3+20,level_py-30,level_py) = number_density*1e3
  endif



    pyramid_global_xHI_array = xHI
    pyramid_global_xHeI_array = xHeI
    pyramid_global_xHeII_array = xHeII

  end subroutine pyramid_field_generation
	
  subroutine short_field_generation()
	   
    implicit none
    integer :: i,j,k
real(kind=dp) :: r2
  if (use_which_field .eq. 'C') then	      
      do i = 1,level_sl	  
        do j = 1,level_sl
           do k = 1,level_sl

              if (i == level_py+1 .and. j==level_py+1 .and. k==level_py+1) then
                 r2 = 0.25* 0.25
              else
                 r2 = (i-level_py-1)**2 + (j-level_py-1)**2 + (k-level_py-1)**2
              endif
              short_global_number_density_array(i,j,k) = number_density!/r2			
         enddo
       enddo
     enddo
  elseif (use_which_field .eq. 'O') then 
      do i = 1,level_sl	  
        do j = 1,level_sl
           do k = 1,level_sl
                   short_global_number_density_array(i,j,k) = number_density			
         enddo
       enddo
     enddo
!short_global_number_density_array(level_sl*2/3,level_sl/2+10,level_py+1) = number_density*1e3
!short_global_number_density_array(level_sl*2/3+1,level_sl/2+10,level_py+1) = number_density*1e3
!short_global_number_density_array(level_sl*2/3,level_sl/2+11,level_py+1) = number_density*1e3
!short_global_number_density_array(level_sl*2/3+1,level_sl/2+11,level_py+1) = number_density*1e3
!short_global_number_density_array(level_sl/3,level_sl/2+10,level_py+1) = number_density*1e3
!short_global_number_density_array(level_sl/3+7,level_sl/2-30,level_py+1) = number_density*1e3
!short_global_number_density_array(level_sl/3+8,level_sl/2-30,level_py+1) = number_density*1e3
!short_global_number_density_array(level_sl/3+9,level_sl/2-30,level_py+1) = number_density*1e3
!short_global_number_density_array(level_sl/3+10,level_sl/2-30,level_py+1) = number_density*1e3
!short_global_number_density_array(level_sl/3+11,level_sl/2-30,level_py+1) = number_density*1e3
!short_global_number_density_array(level_sl/3+12,level_sl/2-30,level_py+1) = number_density*1e3
!short_global_number_density_array(level_sl/3+13,level_sl/2-30,level_py+1) = number_density*1e3
!short_global_number_density_array(level_sl/3+14,level_sl/2-30,level_py+1) = number_density*1e3
!short_global_number_density_array(level_sl/3+15,level_sl/2-30,level_py+1) = number_density*1e3
!short_global_number_density_array(level_sl/3+16,level_sl/2-30,level_py+1) = number_density*1e3
!short_global_number_density_array(level_sl/3+17,level_sl/2-30,level_py+1) = number_density*1e3
!short_global_number_density_array(level_sl/3+18,level_sl/2-30,level_py+1) = number_density*1e3
!short_global_number_density_array(level_sl/3+19,level_sl/2-30,level_py+1) = number_density*1e3
!short_global_number_density_array(level_sl/3+20,level_sl/2-30,level_py+1) = number_density*1e3
  endif

    short_global_xHI_array = xHI
    short_global_xHeI_array = xHeI
    short_global_xHeII_array = xHeII
			
  end subroutine short_field_generation

  subroutine long_field_generation()
	  
   implicit none	  
   integer :: i,j,k
real(kind=dp) :: r2
  if (use_which_field .eq. 'C') then	      
     do i = 1,level_sl	  
       do j = 1,level_sl
         do k = 1,level_sl		 
            if (i==level_py+1 .and. j==level_py+1 .and. k==level_py+1) then
               r2 = 0.25*0.25
            else
               r2 = (i-level_py-1)**2 + (j-level_py-1)**2 + (k-level_py-1)**2
               endif
           long_global_number_density_array(i,j,k) = number_density!/r2
		   
        enddo
      enddo
    enddo
  elseif (use_which_field .eq. 'O') then 
     do i = 1,level_sl	  
       do j = 1,level_sl
         do k = 1,level_sl		 
			 
           long_global_number_density_array(i,j,k) = number_density
		   
        enddo
      enddo
    enddo
!long_global_number_density_array(level_sl*2/3,level_sl/2+10,level_py+1) = number_density*1e3
!long_global_number_density_array(level_sl*2/3+1,level_sl/2+10,level_py+1) = number_density*1e3
!long_global_number_density_array(level_sl*2/3,level_sl/2+11,level_py+1) = number_density*1e3
!long_global_number_density_array(level_sl*2/3+1,level_sl/2+11,level_py+1) = number_density*1e3
!long_global_number_density_array(level_sl/3,level_sl/2+10,level_py+1) = number_density*1e3
!long_global_number_density_array(level_sl/3+7,level_sl/2-30,level_py+1) = number_density*1e3
!long_global_number_density_array(level_sl/3+8,level_sl/2-30,level_py+1) = number_density*1e3
!long_global_number_density_array(level_sl/3+9,level_sl/2-30,level_py+1) = number_density*1e3
!long_global_number_density_array(level_sl/3+10,level_sl/2-30,level_py+1) = number_density*1e3
!long_global_number_density_array(level_sl/3+11,level_sl/2-30,level_py+1) = number_density*1e3
!long_global_number_density_array(level_sl/3+12,level_sl/2-30,level_py+1) = number_density*1e3
!long_global_number_density_array(level_sl/3+13,level_sl/2-30,level_py+1) = number_density*1e3
!long_global_number_density_array(level_sl/3+14,level_sl/2-30,level_py+1) = number_density*1e3
!long_global_number_density_array(level_sl/3+15,level_sl/2-30,level_py+1) = number_density*1e3
!long_global_number_density_array(level_sl/3+16,level_sl/2-30,level_py+1) = number_density*1e3
!long_global_number_density_array(level_sl/3+17,level_sl/2-30,level_py+1) = number_density*1e3
!long_global_number_density_array(level_sl/3+18,level_sl/2-30,level_py+1) = number_density*1e3
!long_global_number_density_array(level_sl/3+19,level_sl/2-30,level_py+1) = number_density*1e3
!long_global_number_density_array(level_sl/3+20,level_sl/2-30,level_py+1) = number_density*1e3
  endif 

    long_global_xHI_array = xHI
    long_global_xHeI_array = xHeI
    long_global_xHeII_array = xHeII
				
  end subroutine long_field_generation
 
  subroutine pyramid_analytical_field_generation()
	   
    implicit none
	
    pyramid_analytical_number_density = number_density
    pyramid_analytical_xHI = xHI 
    pyramid_analytical_xHeI = xHeI
    pyramid_analytical_xHeII = xHeII  
	  
  end subroutine pyramid_analytical_field_generation	     	

  subroutine cartesian_analytical_field_generation()
	   
    implicit none
	
    cartesian_analytical_number_density = number_density
    cartesian_analytical_xHI = xHI 
    cartesian_analytical_xHeI = xHeI
    cartesian_analytical_xHeII = xHeII  
	  
  end subroutine cartesian_analytical_field_generation	

end module field
