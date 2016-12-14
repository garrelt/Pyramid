module pyramid_species_density
	
  use precision, only: dp
  use input, only: level_py, partition, number_density, xHI
  use array, only: pyramid_grid, &
                   px_py_pz_dx, px_py_pz_dy, px_py_pz_dz, &
				   px_py_nz_dx, px_py_nz_dy, px_py_nz_dz, &
				   px_ny_pz_dx, px_ny_pz_dy, px_ny_pz_dz, &	  
				   px_ny_nz_dx, px_ny_nz_dy, px_ny_nz_dz, & 	   
				   nx_py_pz_dx, nx_py_pz_dy, nx_py_pz_dz, &	  
				   nx_py_nz_dx, nx_py_nz_dy, nx_py_nz_dz, &	
				   nx_ny_pz_dx, nx_ny_pz_dy, nx_ny_pz_dz, &	  
				   nx_ny_nz_dx, nx_ny_nz_dy, nx_ny_nz_dz, & 			   	   
                   pyramid_source_px_py_pz_HI_density_array, &
                   pyramid_source_px_py_nz_HI_density_array, &
                   pyramid_source_px_ny_pz_HI_density_array, &
                   pyramid_source_px_ny_nz_HI_density_array, &
                   pyramid_source_nx_py_pz_HI_density_array, &
                   pyramid_source_nx_py_nz_HI_density_array, &
                   pyramid_source_nx_ny_pz_HI_density_array, &
                   pyramid_source_nx_ny_nz_HI_density_array 
				   
  implicit none
	  
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pyramid_source_to_domain_HI_density_transformation ()
	
    call px_py_pz_HI_density_calculation ()	
    call px_py_nz_HI_density_calculation ()	
    call px_ny_pz_HI_density_calculation ()	
    call px_ny_nz_HI_density_calculation ()	
    call nx_py_pz_HI_density_calculation ()	
    call nx_py_nz_HI_density_calculation ()	
    call nx_ny_pz_HI_density_calculation ()	
    call nx_ny_nz_HI_density_calculation ()	
			
  end subroutine pyramid_source_to_domain_HI_density_transformation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine px_py_pz_HI_density_calculation ()
	
  	integer :: i_primary,i_secondary,i_level_py,n
    integer :: a,b
    real :: HI_number_dx,HI_number_dy,HI_number_dz
	
  	!Distance of the center of pyramid_grid to the source
  	do i_level_py = 1,level_py
  	  do i_primary = 1,partition(i_level_py)
  	    do i_secondary = 1,partition(i_level_py)
		  
          HI_number_dx = 0
          HI_number_dy = 0
          HI_number_dz = 0	
		  	  
          do a = lbound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,1),&
	             ubound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,1)	
            do b = lbound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,2),&
                   ubound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,2)

              HI_number_dx = HI_number_dx + &
                             pyramid_source_px_py_pz_HI_density_array(i_level_py,a,b) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HI_number_dy = HI_number_dy + &
                             pyramid_source_px_py_pz_HI_density_array(b,i_level_py,a) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)									 	   
              HI_number_dz = HI_number_dz + &
                             pyramid_source_px_py_pz_HI_density_array(a,b,i_level_py) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
		  
            enddo		  
          enddo		  

          px_py_pz_dx(i_level_py)%HI_density(i_primary,i_secondary) = &
                                  HI_number_dx/pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
          px_py_pz_dy(i_level_py)%HI_density(i_primary,i_secondary) = &
                                  HI_number_dy/pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
          px_py_pz_dz(i_level_py)%HI_density(i_primary,i_secondary) = &
                                  HI_number_dz/pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
		  		  
        end do
  	  end do
  	end do

  end subroutine px_py_pz_HI_density_calculation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine px_py_nz_HI_density_calculation ()

  	integer :: i_primary,i_secondary,i_level_py,n
    integer :: a,b
    real :: HI_number_dx,HI_number_dy,HI_number_dz
	
  	!Distance of the center of pyramid_grid to the source
  	do i_level_py = 1,level_py
  	  do i_primary = 1,partition(i_level_py)
  	    do i_secondary = 1,partition(i_level_py)
		  
          HI_number_dx = 0
          HI_number_dy = 0
          HI_number_dz = 0	
		  	  
          do a = lbound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,1),&
	             ubound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,1)	
            do b = lbound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,2),&
                   ubound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,2)

              HI_number_dx = HI_number_dx + &
                             pyramid_source_px_py_nz_HI_density_array(i_level_py,a,-b) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HI_number_dy = HI_number_dy + &
                             pyramid_source_px_py_nz_HI_density_array(b,i_level_py,-a) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)		
	          HI_number_dz = HI_number_dz + &
                             pyramid_source_px_py_nz_HI_density_array(a,b,-i_level_py) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
		  
            enddo		  
          enddo		  

          px_py_nz_dx(i_level_py)%HI_density(i_primary,-i_secondary) = &
                                  HI_number_dx/pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
          px_py_nz_dy(i_level_py)%HI_density(-i_primary,i_secondary) = &
                                  HI_number_dy/pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
          px_py_nz_dz(-i_level_py)%HI_density(i_primary,i_secondary) = &
                                  HI_number_dz/pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
		  		  
        end do
  	  end do
  	end do

  end subroutine px_py_nz_HI_density_calculation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine px_ny_pz_HI_density_calculation ()
	
    integer :: i_primary,i_secondary,i_level_py,n
    integer :: a,b
    real :: HI_number_dx,HI_number_dy,HI_number_dz
	
    !Distance of the center of pyramid_grid to the source
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)
		  
          HI_number_dx = 0
          HI_number_dy = 0
          HI_number_dz = 0	
		  	  
          do a = lbound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,1),&
  	             ubound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,1)	
            do b = lbound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,2),&
                   ubound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,2)

               HI_number_dx = HI_number_dx + &
                              pyramid_source_px_ny_pz_HI_density_array(i_level_py,-a,b) * &
                              pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
               HI_number_dy = HI_number_dy + &
                              pyramid_source_px_ny_pz_HI_density_array(b,-i_level_py,a) * &
                              pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)									 	   
               HI_number_dz = HI_number_dz + &
                              pyramid_source_px_ny_pz_HI_density_array(a,-b,i_level_py) * &
                              pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
		  
            enddo		  
          enddo		  

          px_ny_pz_dx(i_level_py)%HI_density(-i_primary,i_secondary) = &
                                  HI_number_dx/pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
          px_ny_pz_dy(-i_level_py)%HI_density(i_primary,i_secondary) = &
                                  HI_number_dy/pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
          px_ny_pz_dz(i_level_py)%HI_density(i_primary,-i_secondary) = &
                                  HI_number_dz/pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
		  		  
        end do
      end do
    end do

  end subroutine px_ny_pz_HI_density_calculation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine px_ny_nz_HI_density_calculation ()
	
  	integer :: i_primary,i_secondary,i_level_py,n
    integer :: a,b
    real :: HI_number_dx,HI_number_dy,HI_number_dz
	
    !Distance of the center of pyramid_grid to the source
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)
		  
          HI_number_dx = 0
          HI_number_dy = 0
          HI_number_dz = 0	
		  	  
          do a = lbound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,1),&
                 ubound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,1)	
            do b = lbound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,2),&
                   ubound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,2)

              HI_number_dx = HI_number_dx + &
                             pyramid_source_px_ny_nz_HI_density_array(i_level_py,-a,-b) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HI_number_dy = HI_number_dy + &
                             pyramid_source_px_ny_nz_HI_density_array(b,-i_level_py,-a) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)									 	   
              HI_number_dz = HI_number_dz + &
                             pyramid_source_px_ny_nz_HI_density_array(a,-b,-i_level_py) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
		  
            enddo		  
          enddo		  

          px_ny_nz_dx(i_level_py)%HI_density(-i_primary,-i_secondary) = &
                                  HI_number_dx/pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
          px_ny_nz_dy(-i_level_py)%HI_density(-i_primary,i_secondary) = &
                                  HI_number_dy/pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
          px_ny_nz_dz(-i_level_py)%HI_density(i_primary,-i_secondary) = &
                                  HI_number_dz/pyramid_grid(i_level_py)%volume(i_primary,i_secondary)

        end do
  	  end do
  	end do

  end subroutine px_ny_nz_HI_density_calculation
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nx_py_pz_HI_density_calculation ()
	
   	integer :: i_primary,i_secondary,i_level_py,n
    integer :: a,b
    real :: HI_number_dx,HI_number_dy,HI_number_dz
	
    !Distance of the center of pyramid_grid to the source
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)
		  
          HI_number_dx = 0
          HI_number_dy = 0
          HI_number_dz = 0	
		  	  
          do a = lbound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,1),&
                 ubound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,1)	
            do b = lbound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,2),&
                   ubound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,2)

              HI_number_dx = HI_number_dx + &
                             pyramid_source_nx_py_pz_HI_density_array(-i_level_py,a,b) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HI_number_dy = HI_number_dy + &
                             pyramid_source_nx_py_pz_HI_density_array(-b,i_level_py,a) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)									 	   
              HI_number_dz = HI_number_dz + &
                             pyramid_source_nx_py_pz_HI_density_array(-a,b,i_level_py) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)

            enddo		  
          enddo		  

          nx_py_pz_dx(-i_level_py)%HI_density(i_primary,i_secondary) = &
                                  HI_number_dx/pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
          nx_py_pz_dy(i_level_py)%HI_density(i_primary,-i_secondary) = &
                                  HI_number_dy/pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
          nx_py_pz_dz(i_level_py)%HI_density(-i_primary,i_secondary) = &
                                  HI_number_dz/pyramid_grid(i_level_py)%volume(i_primary,i_secondary)

        end do
      end do
    end do

  end subroutine nx_py_pz_HI_density_calculation
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nx_py_nz_HI_density_calculation ()

    integer :: i_primary,i_secondary,i_level_py,n
    integer :: a,b
    real :: HI_number_dx,HI_number_dy,HI_number_dz
	
    !Distance of the center of pyramid_grid to the source
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)
		  
          HI_number_dx = 0
          HI_number_dy = 0
          HI_number_dz = 0	
		  	  
          do a = lbound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,1),&
  	             ubound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,1)	
            do b = lbound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,2),&
                   ubound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,2)

              HI_number_dx = HI_number_dx + &
                             pyramid_source_nx_py_nz_HI_density_array(-i_level_py,a,-b) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HI_number_dy = HI_number_dy + &
                             pyramid_source_nx_py_nz_HI_density_array(-b,i_level_py,-a) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)		
  	          HI_number_dz = HI_number_dz + &
                             pyramid_source_nx_py_nz_HI_density_array(-a,b,-i_level_py) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
		  
            enddo		  
          enddo		  

          nx_py_nz_dx(-i_level_py)%HI_density(i_primary,-i_secondary) = &
                                  HI_number_dx/pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
          nx_py_nz_dy(i_level_py)%HI_density(-i_primary,-i_secondary) = &
                                  HI_number_dy/pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
          nx_py_nz_dz(-i_level_py)%HI_density(-i_primary,i_secondary) = &
                                  HI_number_dz/pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
		  		  
        end do
      end do
    end do

  end subroutine nx_py_nz_HI_density_calculation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nx_ny_pz_HI_density_calculation ()
	
    integer :: i_primary,i_secondary,i_level_py,n
    integer :: a,b
    real :: HI_number_dx,HI_number_dy,HI_number_dz
	
    !Distance of the center of pyramid_grid to the source
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)
		  
          HI_number_dx = 0
          HI_number_dy = 0
          HI_number_dz = 0	
		  	  
          do a = lbound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,1),&
                 ubound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,1)	
            do b = lbound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,2),&
                   ubound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,2)

               HI_number_dx = HI_number_dx + &
                              pyramid_source_nx_ny_pz_HI_density_array(-i_level_py,-a,b) * &
                              pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
               HI_number_dy = HI_number_dy + &
                              pyramid_source_nx_ny_pz_HI_density_array(-b,-i_level_py,a) * &
                              pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)									 	   
               HI_number_dz = HI_number_dz + &
                              pyramid_source_nx_ny_pz_HI_density_array(-a,-b,i_level_py) * &
                              pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
		  
            enddo		  
          enddo		  

          nx_ny_pz_dx(-i_level_py)%HI_density(-i_primary,i_secondary) = &
                                  HI_number_dx/pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
          nx_ny_pz_dy(-i_level_py)%HI_density(i_primary,-i_secondary) = &
                                  HI_number_dy/pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
          nx_ny_pz_dz(i_level_py)%HI_density(-i_primary,-i_secondary) = &
                                  HI_number_dz/pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
		  		  
        end do
      end do
    end do

  end subroutine nx_ny_pz_HI_density_calculation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nx_ny_nz_HI_density_calculation ()
	
    integer :: i_primary,i_secondary,i_level_py,n
    integer :: a,b
    real :: HI_number_dx,HI_number_dy,HI_number_dz
	
    !Distance of the center of pyramid_grid to the source
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)
		  
          HI_number_dx = 0
          HI_number_dy = 0
          HI_number_dz = 0	
		  	  
          do a = lbound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,1),&
                 ubound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,1)	
            do b = lbound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,2),&
                   ubound(pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume,2)

              HI_number_dx = HI_number_dx + &
                             pyramid_source_nx_ny_nz_HI_density_array(-i_level_py,-a,-b) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HI_number_dy = HI_number_dy + &
                             pyramid_source_nx_ny_nz_HI_density_array(-b,-i_level_py,-a) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)									 	   
              HI_number_dz = HI_number_dz + &
                             pyramid_source_nx_ny_nz_HI_density_array(-a,-b,-i_level_py) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
		  
            enddo		  
          enddo		  

          nx_ny_nz_dx(-i_level_py)%HI_density(-i_primary,-i_secondary) = &
                                  HI_number_dx/pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
          nx_ny_nz_dy(-i_level_py)%HI_density(-i_primary,-i_secondary) = &
                                  HI_number_dy/pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
          nx_ny_nz_dz(-i_level_py)%HI_density(-i_primary,-i_secondary) = &
                                  HI_number_dz/pyramid_grid(i_level_py)%volume(i_primary,i_secondary)

        end do
      end do
    end do

  end subroutine nx_ny_nz_HI_density_calculation			

end module pyramid_species_density