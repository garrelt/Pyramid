module output
  
  use precision, only: dp
  use array, only: pyramid_global_HI_photoionization_rate_array, &
                   pyramid_global_HeI_photoionization_rate_array, &
                   pyramid_global_HeII_photoionization_rate_array, &
                   pyramid_global_photoheating_rate_array, & 
                   short_global_HI_photoionization_rate_array, short_global_HeI_photoionization_rate_array, &
                   short_global_HeII_photoionization_rate_array, short_global_photoheating_rate_array, &				   	   
                   long_global_HI_photoionization_rate_array, long_global_HeI_photoionization_rate_array, &				   
                   long_global_HeII_photoionization_rate_array, long_global_photoheating_rate_array, &				   
                   pyramid_analytical_HI_photoionization_rate_array, pyramid_analytical_HeI_photoionization_rate_array, &
                   pyramid_analytical_HeII_photoionization_rate_array, pyramid_analytical_photoheating_rate_array, &
                   cartesian_analytical_HI_photoionization_rate_array, &
                   cartesian_analytical_HeI_photoionization_rate_array, &				   
                   cartesian_analytical_HeII_photoionization_rate_array, &
                   cartesian_analytical_photoheating_rate_array
  use array, only : interpolated_HI_photoionization_rate_array
  use array, only : interpolated_HeI_photoionization_rate_array  
  use array, only : interpolated_HeII_photoionization_rate_array  
  use array, only : interpolated_photoheating_rate_array
  use array, only : long_global_HI_photoionization_rate_array
  use array, only : long_global_HeI_photoionization_rate_array  
  use array, only : long_global_HeII_photoionization_rate_array  
  use array, only : long_global_photoheating_rate_array 
				   				   
  use input, only: level_py, box_size,level_sl
  
  implicit none

contains

  subroutine output_global()

    implicit none

    integer :: i,j,k
    real(kind=dp) :: r
    real(kind=dp), dimension (:,:,:), allocatable :: pyramid_global_HI_photoionization_rate_r2_array
    real(kind=dp), dimension (:,:,:), allocatable :: short_global_HI_photoionization_rate_r2_array
    real(kind=dp), dimension (:,:,:), allocatable :: long_global_HI_photoionization_rate_r2_array
    allocate (pyramid_global_HI_photoionization_rate_r2_array(1:2*level_py,1:2*level_py,1:2*level_py))
    allocate (short_global_HI_photoionization_rate_r2_array(1:level_sl,1:level_sl,1:level_sl))
    allocate (long_global_HI_photoionization_rate_r2_array(1:level_sl,1:level_sl,1:level_sl))
    do i = 1, 2*level_py
       do j = 1, 2*level_py
          do k = 1, 2*level_py
             r = ((i-level_py-0.5)**2+(j-level_py-0.5)**2+(k-level_py-0.5)**2)*(box_size/real(2*level_py))**2
             pyramid_global_HI_photoionization_rate_r2_array(i,j,k)=pyramid_global_HI_photoionization_rate_array(i,j,k)*r
          enddo
       enddo
    enddo

    do i = 1, level_sl
       do j = 1, level_sl
          do k = 1, level_sl
             r = ((i-level_py-1)**2+(j-level_py-1)**2+(k-level_py-1)**2)*(box_size/real(level_sl))**2
             if (i.eq.level_py+1 .and. j.eq.level_py+1 .and. k.eq.level_py+1) then
                r = (box_size/real(level_sl) * 0.25)**2
                endif
             short_global_HI_photoionization_rate_r2_array(i,j,k)=short_global_HI_photoionization_rate_array(i,j,k)*r
             long_global_HI_photoionization_rate_r2_array(i,j,k)=long_global_HI_photoionization_rate_array(i,j,k)*r
          enddo
       enddo
    enddo

    open(unit=41,file='pyramid_diagonal_photo_rate.dat',form="formatted",status="unknown")
    do i=1,level_py
      r = abs(sqrt(3.0)*real(i-0.5-level_py)*box_size/real(2*level_py))
      write(41,*) r, &
                  pyramid_global_HI_photoionization_rate_array(i,i,i), &
                  pyramid_global_HeI_photoionization_rate_array(i,i,i), &				  
                  pyramid_global_HeII_photoionization_rate_array(i,i,i), &				  
                  pyramid_global_photoheating_rate_array(i,i,i) 
    enddo
    close(41)

    open(unit=41,file='pyramid_diagonal_photo_rate_from3DData.dat',form="formatted",status="unknown")
    do i=1,level_py
       r = sqrt(3.0)*real(i-0.5-level_py)*box_size/real(2*level_py)
       write(41,*) r, pyramid_global_HI_photoionization_rate_r2_array(i,i,i)
    enddo
    close(41)
    
    open(unit=42,file='short_diagonal_photo_rate.dat',form="formatted",status="unknown")
    do i=1,level_py
      r = abs(sqrt(3.0)*real(i-(level_sl+1)/2)*box_size/real(level_sl))
      write(42,*) r, &
                  short_global_HI_photoionization_rate_array(i,i,i), &
                  short_global_HeI_photoionization_rate_array(i,i,i), &				  
                  short_global_HeII_photoionization_rate_array(i,i,i), &
                  short_global_photoheating_rate_array(i,i,i) 			  				  		  
    enddo
    close(42)

    open(unit=42,file='short_diagonal_photo_rate_from3DData.dat',form="formatted",status="unknown")
    do i=1,level_py
       r = sqrt(3.0)*real(i-(level_sl+1)/2)*box_size/real(level_sl)
       write(42,*) r, short_global_HI_photoionization_rate_r2_array(i,i,i)
    enddo
    close(42)
    
    open(unit=43,file='long_diagonal_photo_rate.dat',form="formatted",status="unknown")
    do i=1,level_py
      r = abs(sqrt(3.0)*real(i-(level_sl+1)/2)*box_size/real(level_sl))
      write(43,*) r, &
                  long_global_HI_photoionization_rate_array(i,i,i), &	
                  long_global_HeI_photoionization_rate_array(i,i,i), &					  
                  long_global_HeII_photoionization_rate_array(i,i,i), &					  
                  long_global_photoheating_rate_array(i,i,i)				  	  
    enddo
    close(43)

    open(unit=43,file='long_diagonal_photo_rate_from3DData.dat',form="formatted",status="unknown")
    do i=1,level_py
       r = sqrt(3.0)*real(i-(level_sl+1)/2)*box_size/real(level_sl)
       write(43,*) r, long_global_HI_photoionization_rate_r2_array(i,i,i)
    enddo
    close(43)
    
    open(unit=44,file='pyramid_analytical_diagonal_photo_rate.dat',form="formatted",status="unknown")
    do i=1,2*level_py
      write(44,*) sqrt(3.0)*real(i-0.5-level_py)*box_size/real(2*level_py), &
                  pyramid_analytical_HI_photoionization_rate_array(i,i,i), &
                  pyramid_analytical_HeI_photoionization_rate_array(i,i,i), &				  
                  pyramid_analytical_HeII_photoionization_rate_array(i,i,i), &				  
                  pyramid_analytical_photoheating_rate_array(i,i,i)
    enddo
    close(44)

    open(unit=45,file='cartesian_analytical_diagonal_photo_rate.dat',form="formatted",status="unknown")
    do i=-level_py,level_py !1,level_sl
      write(45,*) sqrt(3.0)*real(i-(level_sl+1)/2)*box_size/real(level_sl), &
                  cartesian_analytical_HI_photoionization_rate_array(i,i,i), &
                  cartesian_analytical_HeI_photoionization_rate_array(i,i,i), &				  
                  cartesian_analytical_HeII_photoionization_rate_array(i,i,i), &				  
                  cartesian_analytical_photoheating_rate_array(i,i,i)
    enddo
    close(45)
		
    open(unit=46,file='pyramid_3D_HI_photoionization.dat',form="unformatted",status="unknown")
    write(46) 2*level_py,2*level_py,2*level_py
    write(46) (((log10(pyramid_global_HI_photoionization_rate_array(i,j,k)),i=1,2*level_py),&
         j=1,2*level_py),k=1,2*level_py)
    close(46)
    write(*,*) "near source py rate", log10(pyramid_global_HI_photoionization_rate_array(level_py,level_py,level_py))
    write(*,*) "far source py rate", log10(pyramid_global_HI_photoionization_rate_array(1,1,1))
    
    open(unit=47,file='pyramid_3D_HeI_photoionization.dat',form="unformatted",status="unknown")
    write(47) 2*level_py,2*level_py,2*level_py
    write(47) (((log10(pyramid_global_HeI_photoionization_rate_array(i,j,k)),i=1,2*level_py),&
	      j=1,2*level_py),k=1,2*level_py)
    close(47)

    open(unit=48,file='pyramid_3D_HeII_photoionization.dat',form="unformatted",status="unknown")
    write(48) 2*level_py,2*level_py,2*level_py
    write(48) (((log10(pyramid_global_HeII_photoionization_rate_array(i,j,k)),i=1,2*level_py),&
	      j=1,2*level_py),k=1,2*level_py)
    close(48)

    open(unit=49,file='pyramid_3D_photoheating.dat',form="unformatted",status="unknown")
    write(49) 2*level_py,2*level_py,2*level_py
    write(49) (((log10(pyramid_global_photoheating_rate_array(i,j,k)),i=1,2*level_py),&
	      j=1,2*level_py),k=1,2*level_py)
    close(49)

    open(unit=50,file='short_3D_HI_photoionization.dat',form="unformatted",status="unknown")
    write(50) level_sl,level_sl,level_sl
    write(50) (((log10(short_global_HI_photoionization_rate_array(i,j,k)),i=1,level_sl),&
         j=1,level_sl),k=1,level_sl)
    close(50)
    write(*,*) "near source sh rate",log10(short_global_HI_photoionization_rate_r2_array(level_py+1,level_py+1,level_py+1))
    write(*,*) "far source sh rate", log10(short_global_HI_photoionization_rate_r2_array(1,1,1))

    open(unit=51,file='short_3D_HeI_photoionization.dat',form="unformatted",status="unknown")
    write(51) level_sl,level_sl,level_sl
    write(51) (((log10(short_global_HeI_photoionization_rate_array(i,j,k)),i=1,level_sl),&
	      j=1,level_sl),k=1,level_sl)
    close(51)

    open(unit=52,file='short_3D_HeII_photoionization.dat',form="unformatted",status="unknown")
    write(52) level_sl,level_sl,level_sl
    write(52) (((log10(short_global_HeII_photoionization_rate_array(i,j,k)),i=1,level_sl),&
	      j=1,level_sl),k=1,level_sl)
    close(52)

    open(unit=53,file='short_3D_photoheating.dat',form="unformatted",status="unknown")
    write(53) level_sl,level_sl,level_sl
    write(53) (((log10(short_global_photoheating_rate_array(i,j,k)),i=1,level_sl),&
	      j=1,level_sl),k=1,level_sl)	
    close(53)

    open(unit=54,file='long_3D_HI_photoionization.dat',form="unformatted",status="unknown")
    write(54) level_sl,level_sl,level_sl
    write(54) (((log10(long_global_HI_photoionization_rate_array(i,j,k)),i=1,level_sl),&
               j=1,level_sl),k=1,level_sl)
    write(*,*) "near source lo rate",log10(long_global_HI_photoionization_rate_array(level_py+1,level_py+1,level_py+1))
        write(*,*) "far source lo rate", log10(long_global_HI_photoionization_rate_array(1,1,1))

    open(unit=55,file='long_3D_HeI_photoionization.dat',form="unformatted",status="unknown")
    write(55) level_sl,level_sl,level_sl
    write(55) (((log10(long_global_HeI_photoionization_rate_array(i,j,k)),i=1,level_sl),&
	      j=1,level_sl),k=1,level_sl)
    close(55)	

    open(unit=56,file='long_3D_HeII_photoionization.dat',form="unformatted",status="unknown")
    write(56) level_sl,level_sl,level_sl
    write(56) (((log10(long_global_HeII_photoionization_rate_array(i,j,k)),i=1,level_sl),&
	      j=1,level_sl),k=1,level_sl)
    close(56)

    open(unit=57,file='long_3D_photoheating.dat',form="unformatted",status="unknown")
    write(57) level_sl,level_sl,level_sl
    write(57) (((log10(long_global_photoheating_rate_array(i,j,k)),i=1,level_sl),&
	      j=1,level_sl),k=1,level_sl)	
    close(57)

    open(unit=62,file='interpolated_3D_HI_photoionization.dat',form="unformatted",status="unknown")
    write(62) 2*level_py,2*level_py,2*level_py
    write(62) (((log10(interpolated_HI_photoionization_rate_array(i,j,k)),i=1,2*level_py),&
	      j=1,2*level_py),k=1,2*level_py)
    close(62)

    open(unit=63,file='interpolated_3D_HeI_photoionization.dat',form="unformatted",status="unknown")
    write(63) 2*level_py,2*level_py,2*level_py
    write(63) (((log10(interpolated_HeI_photoionization_rate_array(i,j,k)),i=1,2*level_py),&
	      j=1,2*level_py),k=1,2*level_py)
    close(63)

    open(unit=64,file='interpolated_3D_HeII_photoionization.dat',form="unformatted",status="unknown")
    write(64) 2*level_py,2*level_py,2*level_py
    write(64) (((log10(interpolated_HeII_photoionization_rate_array(i,j,k)),i=1,2*level_py),&
	      j=1,2*level_py),k=1,2*level_py)
    close(64)

    open(unit=65,file='interpolated_3D_photoheating.dat',form="unformatted",status="unknown")
    write(65) 2*level_py,2*level_py,2*level_py
    write(65) (((log10(interpolated_photoheating_rate_array(i,j,k)),i=1,2*level_py),&
	      j=1,2*level_py),k=1,2*level_py)
    close(65)

    open(unit=66,file='pyramid_analytical_3D_HI_photoionization.dat',form="unformatted",status="unknown")
    write(66) 2*level_py,2*level_py,2*level_py
    write(66) (((log10(pyramid_analytical_HI_photoionization_rate_array(i,j,k)),i=1,2*level_py),&
	      j=1,2*level_py),k=1,2*level_py)
    close(66)
	
    open(unit=67,file='pyramid_analytical_3D_HeI_photoionization.dat',form="unformatted",status="unknown")
    write(67) 2*level_py,2*level_py,2*level_py
    write(67) (((log10(pyramid_analytical_HeI_photoionization_rate_array(i,j,k)),i=1,2*level_py),&
	      j=1,2*level_py),k=1,2*level_py)
    close(67)
	
    open(unit=68,file='pyramid_analytical_3D_HeII_photoionization.dat',form="unformatted",status="unknown")
    write(68) 2*level_py,2*level_py,2*level_py
    write(68) (((log10(pyramid_analytical_HeII_photoionization_rate_array(i,j,k)),i=1,2*level_py),&
	      j=1,2*level_py),k=1,2*level_py)
    close(68)

    open(unit=69,file='pyramid_analytical_3D_photoheating.dat',form="unformatted",status="unknown")
    write(69) 2*level_py,2*level_py,2*level_py
    write(69) (((log10(pyramid_analytical_photoheating_rate_array(i,j,k)),i=1,2*level_py),&
	      j=1,2*level_py),k=1,2*level_py)
    close(69)

    open(unit=70,file='cartesian_analytical_3D_HI_photoionization.dat',form="unformatted",status="unknown")
    write(70) level_sl,level_sl,level_sl
    write(70) (((log10(cartesian_analytical_HI_photoionization_rate_array(i,j,k)),i=-level_py,level_py),&
	      j=-level_py,level_py),k=-level_py,level_py)
    close(70)

    open(unit=71,file='cartesian_analytical_3D_HeI_photoionization.dat',form="unformatted",status="unknown")
    write(71) level_sl,level_sl,level_sl
    write(71) (((log10(cartesian_analytical_HI_photoionization_rate_array(i,j,k)),i=-level_py,level_py),&
	      j=-level_py,level_py),k=-level_py,level_py)
    close(71)
	
    open(unit=72,file='cartesian_analytical_3D_HeII_photoionization.dat',form="unformatted",status="unknown")
    write(72) level_sl,level_sl,level_sl
    write(72) (((log10(cartesian_analytical_HI_photoionization_rate_array(i,j,k)),i=-level_py,level_py),&
              j=-level_py,level_py),k=-level_py,level_py)
    close(72)
	
    open(unit=73,file='cartesian_analytical_3D_photoheating.dat',form="unformatted",status="unknown")
    write(73) level_sl,level_sl,level_sl
    write(73) (((log10(cartesian_analytical_photoheating_rate_array(i,j,k)),i=-level_py,level_py),&
              j=-level_py,level_py),k=-level_py,level_py)
    close(73)

    open(unit=74,file='pyramid_analytical_vs_pyramid_HI_photoionization.dat',form="unformatted",status="unknown")
    write(74) 2*level_py,2*level_py,2*level_py
    write(74) ((((pyramid_global_HI_photoionization_rate_array(i,j,k)-&
              pyramid_analytical_HI_photoionization_rate_array(i,j,k))/&
              pyramid_analytical_HI_photoionization_rate_array(i,j,k),&
              i=1,2*level_py),j=1,2*level_py),k=1,2*level_py)
    close(74)

    open(unit=75,file='pyramid_analytical_vs_pyramid_HeI_photoionization.dat',form="unformatted",status="unknown")
    write(75) 2*level_py,2*level_py,2*level_py
    write(75) ((((pyramid_global_HeI_photoionization_rate_array(i,j,k)-&
              pyramid_analytical_HeI_photoionization_rate_array(i,j,k))/&
              pyramid_analytical_HeI_photoionization_rate_array(i,j,k),&
              i=1,2*level_py),j=1,2*level_py),k=1,2*level_py)
    close(75)

    open(unit=76,file='pyramid_analytical_vs_pyramid_HeII_photoionization.dat',form="unformatted",status="unknown")
    write(76) 2*level_py,2*level_py,2*level_py
    write(76) ((((pyramid_global_HeII_photoionization_rate_array(i,j,k)-&
              pyramid_analytical_HeII_photoionization_rate_array(i,j,k))/&
              pyramid_analytical_HeII_photoionization_rate_array(i,j,k),&
              i=1,2*level_py),j=1,2*level_py),k=1,2*level_py)
    close(76)

    open(unit=77,file='pyramid_analytical_vs_pyramid_photoheating.dat',form="unformatted",status="unknown")
    write(77) 2*level_py,2*level_py,2*level_py
    write(77) ((((pyramid_global_photoheating_rate_array(i,j,k)-&
              pyramid_analytical_photoheating_rate_array(i,j,k))/&
              pyramid_analytical_photoheating_rate_array(i,j,k),&
              i=1,2*level_py),j=1,2*level_py),k=1,2*level_py)
    close(77)

    open(unit=78,file='interpolated_vs_pyramid_HI_photoionization.dat',form="unformatted",status="unknown")
    write(78) 2*level_py,2*level_py,2*level_py
    write(78) ((((pyramid_global_HI_photoionization_rate_array(i,j,k)-&
              interpolated_HI_photoionization_rate_array(i,j,k))/&
              interpolated_HI_photoionization_rate_array(i,j,k),&
              i=1,2*level_py),j=1,2*level_py),k=1,2*level_py)   
    close(78)

    open(unit=79,file='interpolated_vs_pyramid_HeI_photoionization.dat',form="unformatted",status="unknown")
    write(79) 2*level_py,2*level_py,2*level_py
    write(79) ((((pyramid_global_HeI_photoionization_rate_array(i,j,k)-&
              interpolated_HeI_photoionization_rate_array(i,j,k))/&
              interpolated_HeI_photoionization_rate_array(i,j,k),&
              i=1,2*level_py),j=1,2*level_py),k=1,2*level_py)   
    close(79)

    open(unit=80,file='interpolated_vs_pyramid_HeII_photoionization.dat',form="unformatted",status="unknown")
    write(80) 2*level_py,2*level_py,2*level_py
    write(80) ((((pyramid_global_HeII_photoionization_rate_array(i,j,k)-&
              interpolated_HeII_photoionization_rate_array(i,j,k))/&
              interpolated_HeII_photoionization_rate_array(i,j,k),&
              i=1,2*level_py),j=1,2*level_py),k=1,2*level_py)   
    close(80)

    open(unit=81,file='interpolated_vs_pyramid_photoheating.dat',form="unformatted",status="unknown")
    write(81) 2*level_py,2*level_py,2*level_py
    write(81) ((((pyramid_global_photoheating_rate_array(i,j,k)-&
              interpolated_photoheating_rate_array(i,j,k))/&
              interpolated_photoheating_rate_array(i,j,k),&
              i=1,2*level_py),j=1,2*level_py),k=1,2*level_py)   
    close(81)
	
    open(unit=82,file='pyramid_analytical_vs_interpolated_HI_photoionization.dat',form="unformatted",status="unknown")
    write(82) 2*level_py,2*level_py,2*level_py
    write(82) ((((pyramid_analytical_HI_photoionization_rate_array(i,j,k)-&
              interpolated_HI_photoionization_rate_array(i,j,k))/&
              interpolated_HI_photoionization_rate_array(i,j,k),&
              i=1,2*level_py),j=1,2*level_py),k=1,2*level_py)
    close(82)

    open(unit=83,file='pyramid_analytical_vs_interpolated_HeI_photoionization.dat',form="unformatted",status="unknown")
    write(83) 2*level_py,2*level_py,2*level_py
    write(83) ((((pyramid_analytical_HeI_photoionization_rate_array(i,j,k)-&
              interpolated_HeI_photoionization_rate_array(i,j,k))/&
              interpolated_HeI_photoionization_rate_array(i,j,k),&
              i=1,2*level_py),j=1,2*level_py),k=1,2*level_py)
    close(83)

    open(unit=84,file='pyramid_analytical_vs_interpolated_HeII_photoionization.dat',form="unformatted",status="unknown")
    write(84) 2*level_py,2*level_py,2*level_py
    write(84) ((((pyramid_analytical_HeII_photoionization_rate_array(i,j,k)-&
              interpolated_HeII_photoionization_rate_array(i,j,k))/&
              interpolated_HeII_photoionization_rate_array(i,j,k),&
              i=1,2*level_py),j=1,2*level_py),k=1,2*level_py)
    close(84)

    open(unit=85,file='pyramid_analytical_vs_interpolated_photoheating.dat',form="unformatted",status="unknown")
    write(85) 2*level_py,2*level_py,2*level_py
    write(85) ((((pyramid_analytical_photoheating_rate_array(i,j,k)-&
              interpolated_photoheating_rate_array(i,j,k))/&
              interpolated_photoheating_rate_array(i,j,k),&
              i=1,2*level_py),j=1,2*level_py),k=1,2*level_py)
    close(85)

!write(*,*) 'pyramid ', pyramid_global_HI_photoionization_rate_array(50,50,49)
!write(*,*) 'analytical ', pyramid_analytical_HI_photoionization_rate_array(50,50,49)
!write(*,*) 'interpolated ', interpolated_HI_photoionization_rate_array(50,50,49)
!write(*,*) 'but ', (interpolated_HI_photoionization_rate_array(50,50,49)-&
!              pyramid_analytical_HI_photoionization_rate_array(50,50,49))/&
!              pyramid_analytical_HI_photoionization_rate_array(50,50,49)
  end subroutine output_global

end module output
