module pyramid_source_domain_transformation
	
  use precision, only: dp
  use input, only: level_py, partition, cellsize_py_cube, &
                   ionization_weighted_by_atom, heating_weighted_by_atom
  use array, only: cartesian_grid, pyramid_grid, &
                   px_py_pz_dx, px_py_pz_dy, px_py_pz_dz, &
                   px_py_nz_dx, px_py_nz_dy, px_py_nz_dz, &
                   px_ny_pz_dx, px_ny_pz_dy, px_ny_pz_dz, &	  
                   px_ny_nz_dx, px_ny_nz_dy, px_ny_nz_dz, & 	   
                   nx_py_pz_dx, nx_py_pz_dy, nx_py_pz_dz, &	  
                   nx_py_nz_dx, nx_py_nz_dy, nx_py_nz_dz, &	
                   nx_ny_pz_dx, nx_ny_pz_dy, nx_ny_pz_dz, &	  
                   nx_ny_nz_dx, nx_ny_nz_dy, nx_ny_nz_dz, &
                   pyramid_source_px_py_pz_HI_photoionization_rate_array, &
                   pyramid_source_px_py_nz_HI_photoionization_rate_array, &  
                   pyramid_source_px_ny_pz_HI_photoionization_rate_array, &
                   pyramid_source_px_ny_nz_HI_photoionization_rate_array, &
                   pyramid_source_nx_py_pz_HI_photoionization_rate_array, &
                   pyramid_source_nx_py_nz_HI_photoionization_rate_array, &
                   pyramid_source_nx_ny_pz_HI_photoionization_rate_array, &
                   pyramid_source_nx_ny_nz_HI_photoionization_rate_array, &				   
                   pyramid_source_px_py_pz_HeI_photoionization_rate_array, &
                   pyramid_source_px_py_nz_HeI_photoionization_rate_array, &  
                   pyramid_source_px_ny_pz_HeI_photoionization_rate_array, &
                   pyramid_source_px_ny_nz_HeI_photoionization_rate_array, &
                   pyramid_source_nx_py_pz_HeI_photoionization_rate_array, &
                   pyramid_source_nx_py_nz_HeI_photoionization_rate_array, &
                   pyramid_source_nx_ny_pz_HeI_photoionization_rate_array, &
                   pyramid_source_nx_ny_nz_HeI_photoionization_rate_array, &				   				   
                   pyramid_source_px_py_pz_HeII_photoionization_rate_array, &
                   pyramid_source_px_py_nz_HeII_photoionization_rate_array, &  
                   pyramid_source_px_ny_pz_HeII_photoionization_rate_array, &
                   pyramid_source_px_ny_nz_HeII_photoionization_rate_array, &
                   pyramid_source_nx_py_pz_HeII_photoionization_rate_array, &
                   pyramid_source_nx_py_nz_HeII_photoionization_rate_array, &
                   pyramid_source_nx_ny_pz_HeII_photoionization_rate_array, &
                   pyramid_source_nx_ny_nz_HeII_photoionization_rate_array, &				   				   
                   pyramid_source_px_py_pz_photoheating_rate_array, &
                   pyramid_source_px_py_nz_photoheating_rate_array, &
                   pyramid_source_px_ny_pz_photoheating_rate_array, &
                   pyramid_source_px_ny_nz_photoheating_rate_array, &
                   pyramid_source_nx_py_pz_photoheating_rate_array, &
                   pyramid_source_nx_py_nz_photoheating_rate_array, &
                   pyramid_source_nx_ny_pz_photoheating_rate_array, &
                   pyramid_source_nx_ny_nz_photoheating_rate_array, &
                   pyramid_source_px_py_pz_HI_density_array, &
                   pyramid_source_px_py_nz_HI_density_array, &
                   pyramid_source_px_ny_pz_HI_density_array, &
                   pyramid_source_px_ny_nz_HI_density_array, &
                   pyramid_source_nx_py_pz_HI_density_array, &
                   pyramid_source_nx_py_nz_HI_density_array, &
                   pyramid_source_nx_ny_pz_HI_density_array, &
                   pyramid_source_nx_ny_nz_HI_density_array, &				   
                   pyramid_source_px_py_pz_HeI_density_array, &
                   pyramid_source_px_py_nz_HeI_density_array, &
                   pyramid_source_px_ny_pz_HeI_density_array, &
                   pyramid_source_px_ny_nz_HeI_density_array, &
                   pyramid_source_nx_py_pz_HeI_density_array, &
                   pyramid_source_nx_py_nz_HeI_density_array, &
                   pyramid_source_nx_ny_pz_HeI_density_array, &
                   pyramid_source_nx_ny_nz_HeI_density_array, &				   
                   pyramid_source_px_py_pz_HeII_density_array, &
                   pyramid_source_px_py_nz_HeII_density_array, &
                   pyramid_source_px_ny_pz_HeII_density_array, &
                   pyramid_source_px_ny_nz_HeII_density_array, &
                   pyramid_source_nx_py_pz_HeII_density_array, &
                   pyramid_source_nx_py_nz_HeII_density_array, &
                   pyramid_source_nx_ny_pz_HeII_density_array, &
                   pyramid_source_nx_ny_nz_HeII_density_array				   
				   				   
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pyramid_source_to_domain_species_density_transformation ()
	
    call px_py_pz_species_density_calculation ()	
    call px_py_nz_species_density_calculation ()	
    call px_ny_pz_species_density_calculation ()	
    call px_ny_nz_species_density_calculation ()	
    call nx_py_pz_species_density_calculation ()	
    call nx_py_nz_species_density_calculation ()	
    call nx_ny_pz_species_density_calculation ()	
    call nx_ny_nz_species_density_calculation ()	

  end subroutine pyramid_source_to_domain_species_density_transformation
	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  	
  subroutine pyramid_domain_to_source_photo_rate_transformation()
 
    call source_photo_rate_initialization()   
    call px_py_pz_photo_rate_calculation() 
    call px_py_nz_photo_rate_calculation() 
    call px_ny_pz_photo_rate_calculation() 
    call px_ny_nz_photo_rate_calculation()
    call nx_py_pz_photo_rate_calculation() 
    call nx_py_nz_photo_rate_calculation()
    call nx_ny_pz_photo_rate_calculation()
    call nx_ny_nz_photo_rate_calculation()
    call photoionization_rate_per_atom() 
    call photoheating_rate_per_volume()
	  	
  end subroutine pyramid_domain_to_source_photo_rate_transformation
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine px_py_pz_species_density_calculation ()
	
    integer :: i_primary,i_secondary,i_level_py,n
    integer :: a,b
    real(kind=dp) :: HI_number_dx,HI_number_dy,HI_number_dz
    real(kind=dp) :: HeI_number_dx,HeI_number_dy,HeI_number_dz
    real(kind=dp) :: HeII_number_dx,HeII_number_dy,HeII_number_dz	

    !Distance of the center of pyramid_grid to the source
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)
		  
          HI_number_dx = 0
          HI_number_dy = 0
          HI_number_dz = 0	
          HeI_number_dx = 0
          HeI_number_dy = 0
          HeI_number_dz = 0
          HeII_number_dx = 0
          HeII_number_dy = 0
          HeII_number_dz = 0
		  		  		  	  							  
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
              HeI_number_dx = HeI_number_dx + &
                             pyramid_source_px_py_pz_HeI_density_array(i_level_py,a,b) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeI_number_dy = HeI_number_dy + &
                             pyramid_source_px_py_pz_HeI_density_array(b,i_level_py,a) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeI_number_dz = HeI_number_dz + &
                             pyramid_source_px_py_pz_HeI_density_array(a,b,i_level_py) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeII_number_dx = HeII_number_dx + &
                             pyramid_source_px_py_pz_HeII_density_array(i_level_py,a,b) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeII_number_dy = HeII_number_dy + &
                             pyramid_source_px_py_pz_HeII_density_array(b,i_level_py,a) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeII_number_dz = HeII_number_dz + &
                             pyramid_source_px_py_pz_HeII_density_array(a,b,i_level_py) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
																						  
            enddo		  
          enddo		  

          px_py_pz_dx(i_level_py)%HI_density(i_primary,i_secondary) = &
                                  HI_number_dx/pyramid_grid(i_level_py)%volume
          px_py_pz_dy(i_level_py)%HI_density(i_primary,i_secondary) = &
                                  HI_number_dy/pyramid_grid(i_level_py)%volume
          px_py_pz_dz(i_level_py)%HI_density(i_primary,i_secondary) = &
                                  HI_number_dz/pyramid_grid(i_level_py)%volume
          px_py_pz_dx(i_level_py)%HeI_density(i_primary,i_secondary) = &
                                  HeI_number_dx/pyramid_grid(i_level_py)%volume
          px_py_pz_dy(i_level_py)%HeI_density(i_primary,i_secondary) = &
                                  HeI_number_dy/pyramid_grid(i_level_py)%volume
          px_py_pz_dz(i_level_py)%HeI_density(i_primary,i_secondary) = &
                                  HeI_number_dz/pyramid_grid(i_level_py)%volume
          px_py_pz_dx(i_level_py)%HeII_density(i_primary,i_secondary) = &
                                  HeII_number_dx/pyramid_grid(i_level_py)%volume
          px_py_pz_dy(i_level_py)%HeII_density(i_primary,i_secondary) = &
                                  HeII_number_dy/pyramid_grid(i_level_py)%volume
          px_py_pz_dz(i_level_py)%HeII_density(i_primary,i_secondary) = &
                                  HeII_number_dz/pyramid_grid(i_level_py)%volume
        end do
      end do
    end do

  end subroutine px_py_pz_species_density_calculation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine px_py_nz_species_density_calculation ()

    integer :: i_primary,i_secondary,i_level_py,n
    integer :: a,b
    real(kind=dp) :: HI_number_dx,HI_number_dy,HI_number_dz
    real(kind=dp) :: HeI_number_dx,HeI_number_dy,HeI_number_dz
    real(kind=dp) :: HeII_number_dx,HeII_number_dy,HeII_number_dz	
	
    !Distance of the center of pyramid_grid to the source
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)
		  
          HI_number_dx = 0
          HI_number_dy = 0
          HI_number_dz = 0	
          HeI_number_dx = 0
          HeI_number_dy = 0
          HeI_number_dz = 0
          HeII_number_dx = 0
          HeII_number_dy = 0
          HeII_number_dz = 0
		  	  
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
              HeI_number_dx = HeI_number_dx + &
                             pyramid_source_px_py_nz_HeI_density_array(i_level_py,a,-b) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeI_number_dy = HeI_number_dy + &
                             pyramid_source_px_py_nz_HeI_density_array(b,i_level_py,-a) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)		
              HeI_number_dz = HeI_number_dz + &
                             pyramid_source_px_py_nz_HeI_density_array(a,b,-i_level_py) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeII_number_dx = HeII_number_dx + &
                             pyramid_source_px_py_nz_HeII_density_array(i_level_py,a,-b) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeII_number_dy = HeII_number_dy + &
                             pyramid_source_px_py_nz_HeII_density_array(b,i_level_py,-a) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)		
              HeII_number_dz = HeII_number_dz + &
                             pyramid_source_px_py_nz_HeII_density_array(a,b,-i_level_py) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
							 																								   													  
            enddo		  
          enddo		  

          px_py_nz_dx(i_level_py)%HI_density(i_primary,-i_secondary) = &
                                  HI_number_dx/pyramid_grid(i_level_py)%volume
          px_py_nz_dy(i_level_py)%HI_density(-i_primary,i_secondary) = &
                                  HI_number_dy/pyramid_grid(i_level_py)%volume
          px_py_nz_dz(-i_level_py)%HI_density(i_primary,i_secondary) = &
                                  HI_number_dz/pyramid_grid(i_level_py)%volume
          px_py_nz_dx(i_level_py)%HeI_density(i_primary,-i_secondary) = &
                                  HeI_number_dx/pyramid_grid(i_level_py)%volume
          px_py_nz_dy(i_level_py)%HeI_density(-i_primary,i_secondary) = &
                                  HeI_number_dy/pyramid_grid(i_level_py)%volume
          px_py_nz_dz(-i_level_py)%HeI_density(i_primary,i_secondary) = &
                                  HeI_number_dz/pyramid_grid(i_level_py)%volume
          px_py_nz_dx(i_level_py)%HeII_density(i_primary,-i_secondary) = &
                                  HeII_number_dx/pyramid_grid(i_level_py)%volume
          px_py_nz_dy(i_level_py)%HeII_density(-i_primary,i_secondary) = &
                                  HeII_number_dy/pyramid_grid(i_level_py)%volume
          px_py_nz_dz(-i_level_py)%HeII_density(i_primary,i_secondary) = &
                                  HeII_number_dz/pyramid_grid(i_level_py)%volume
		  		  
        end do
  	  end do
  	end do

  end subroutine px_py_nz_species_density_calculation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine px_ny_pz_species_density_calculation ()
	
    integer :: i_primary,i_secondary,i_level_py,n
    integer :: a,b
    real(kind=dp) :: HI_number_dx,HI_number_dy,HI_number_dz
    real(kind=dp) :: HeI_number_dx,HeI_number_dy,HeI_number_dz
    real(kind=dp) :: HeII_number_dx,HeII_number_dy,HeII_number_dz	
	
    !Distance of the center of pyramid_grid to the source
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)
		  
          HI_number_dx = 0
          HI_number_dy = 0
          HI_number_dz = 0	
          HeI_number_dx = 0
          HeI_number_dy = 0
          HeI_number_dz = 0
          HeII_number_dx = 0
          HeII_number_dy = 0
          HeII_number_dz = 0
		  	  
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
               HeI_number_dx = HeI_number_dx + &
                              pyramid_source_px_ny_pz_HeI_density_array(i_level_py,-a,b) * &
                              pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
               HeI_number_dy = HeI_number_dy + &
                              pyramid_source_px_ny_pz_HeI_density_array(b,-i_level_py,a) * &
                              pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
               HeI_number_dz = HeI_number_dz + &
                              pyramid_source_px_ny_pz_HeI_density_array(a,-b,i_level_py) * &
                              pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
               HeII_number_dx = HeII_number_dx + &
                              pyramid_source_px_ny_pz_HeII_density_array(i_level_py,-a,b) * &
                              pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
               HeII_number_dy = HeII_number_dy + &
                              pyramid_source_px_ny_pz_HeII_density_array(b,-i_level_py,a) * &
                              pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
               HeII_number_dz = HeII_number_dz + &
                              pyramid_source_px_ny_pz_HeII_density_array(a,-b,i_level_py) * &
                              pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
											 		  
            enddo		  
          enddo		  

          px_ny_pz_dx(i_level_py)%HI_density(-i_primary,i_secondary) = &
                                  HI_number_dx/pyramid_grid(i_level_py)%volume
          px_ny_pz_dy(-i_level_py)%HI_density(i_primary,i_secondary) = &
                                  HI_number_dy/pyramid_grid(i_level_py)%volume
          px_ny_pz_dz(i_level_py)%HI_density(i_primary,-i_secondary) = &
                                  HI_number_dz/pyramid_grid(i_level_py)%volume
          px_ny_pz_dx(i_level_py)%HeI_density(-i_primary,i_secondary) = &
                                  HeI_number_dx/pyramid_grid(i_level_py)%volume
          px_ny_pz_dy(-i_level_py)%HeI_density(i_primary,i_secondary) = &
                                  HeI_number_dy/pyramid_grid(i_level_py)%volume
          px_ny_pz_dz(i_level_py)%HeI_density(i_primary,-i_secondary) = &
                                  HeI_number_dz/pyramid_grid(i_level_py)%volume
          px_ny_pz_dx(i_level_py)%HeII_density(-i_primary,i_secondary) = &
                                  HeII_number_dx/pyramid_grid(i_level_py)%volume
          px_ny_pz_dy(-i_level_py)%HeII_density(i_primary,i_secondary) = &
                                  HeII_number_dy/pyramid_grid(i_level_py)%volume
          px_ny_pz_dz(i_level_py)%HeII_density(i_primary,-i_secondary) = &
                                  HeII_number_dz/pyramid_grid(i_level_py)%volume
		  		  
        end do
      end do
    end do

  end subroutine px_ny_pz_species_density_calculation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine px_ny_nz_species_density_calculation ()
	
    integer :: i_primary,i_secondary,i_level_py,n
    integer :: a,b
    real(kind=dp) :: HI_number_dx,HI_number_dy,HI_number_dz
    real(kind=dp) :: HeI_number_dx,HeI_number_dy,HeI_number_dz
    real(kind=dp) :: HeII_number_dx,HeII_number_dy,HeII_number_dz	
	
    !Distance of the center of pyramid_grid to the source
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)
		  
          HI_number_dx = 0
          HI_number_dy = 0
          HI_number_dz = 0	
          HeI_number_dx = 0
          HeI_number_dy = 0
          HeI_number_dz = 0
          HeII_number_dx = 0
          HeII_number_dy = 0
          HeII_number_dz = 0
		  	  
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
              HeI_number_dx = HeI_number_dx + &
                             pyramid_source_px_ny_nz_HeI_density_array(i_level_py,-a,-b) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeI_number_dy = HeI_number_dy + &
                             pyramid_source_px_ny_nz_HeI_density_array(b,-i_level_py,-a) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeI_number_dz = HeI_number_dz + &
                             pyramid_source_px_ny_nz_HeI_density_array(a,-b,-i_level_py) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeII_number_dx = HeII_number_dx + &
                             pyramid_source_px_ny_nz_HeII_density_array(i_level_py,-a,-b) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeII_number_dy = HeII_number_dy + &
                             pyramid_source_px_ny_nz_HeII_density_array(b,-i_level_py,-a) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeII_number_dz = HeII_number_dz + &
                             pyramid_source_px_ny_nz_HeII_density_array(a,-b,-i_level_py) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)

            enddo		  
          enddo		  

          px_ny_nz_dx(i_level_py)%HI_density(-i_primary,-i_secondary) = &
                                  HI_number_dx/pyramid_grid(i_level_py)%volume
          px_ny_nz_dy(-i_level_py)%HI_density(-i_primary,i_secondary) = &
                                  HI_number_dy/pyramid_grid(i_level_py)%volume
          px_ny_nz_dz(-i_level_py)%HI_density(i_primary,-i_secondary) = &
                                  HI_number_dz/pyramid_grid(i_level_py)%volume
          px_ny_nz_dx(i_level_py)%HeI_density(-i_primary,-i_secondary) = &
                                  HeI_number_dx/pyramid_grid(i_level_py)%volume
          px_ny_nz_dy(-i_level_py)%HeI_density(-i_primary,i_secondary) = &
                                  HeI_number_dy/pyramid_grid(i_level_py)%volume
          px_ny_nz_dz(-i_level_py)%HeI_density(i_primary,-i_secondary) = &
                                  HeI_number_dz/pyramid_grid(i_level_py)%volume
          px_ny_nz_dx(i_level_py)%HeII_density(-i_primary,-i_secondary) = &
                                  HeII_number_dx/pyramid_grid(i_level_py)%volume
          px_ny_nz_dy(-i_level_py)%HeII_density(-i_primary,i_secondary) = &
                                  HeII_number_dy/pyramid_grid(i_level_py)%volume
          px_ny_nz_dz(-i_level_py)%HeII_density(i_primary,-i_secondary) = &
                                  HeII_number_dz/pyramid_grid(i_level_py)%volume

        end do
      end do
    end do

  end subroutine px_ny_nz_species_density_calculation
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nx_py_pz_species_density_calculation ()
	
    integer :: i_primary,i_secondary,i_level_py,n
    integer :: a,b
    real(kind=dp) :: HI_number_dx,HI_number_dy,HI_number_dz
    real(kind=dp) :: HeI_number_dx,HeI_number_dy,HeI_number_dz
    real(kind=dp) :: HeII_number_dx,HeII_number_dy,HeII_number_dz	
	
    !Distance of the center of pyramid_grid to the source
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)
		  
          HI_number_dx = 0
          HI_number_dy = 0
          HI_number_dz = 0	
          HeI_number_dx = 0
          HeI_number_dy = 0
          HeI_number_dz = 0
          HeII_number_dx = 0
          HeII_number_dy = 0
          HeII_number_dz = 0
		  	  
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
              HeI_number_dx = HeI_number_dx + &
                             pyramid_source_nx_py_pz_HeI_density_array(-i_level_py,a,b) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeI_number_dy = HeI_number_dy + &
                             pyramid_source_nx_py_pz_HeI_density_array(-b,i_level_py,a) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeI_number_dz = HeI_number_dz + &
                             pyramid_source_nx_py_pz_HeI_density_array(-a,b,i_level_py) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeII_number_dx = HeII_number_dx + &
                             pyramid_source_nx_py_pz_HeII_density_array(-i_level_py,a,b) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeII_number_dy = HeII_number_dy + &
                             pyramid_source_nx_py_pz_HeII_density_array(-b,i_level_py,a) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeII_number_dz = HeII_number_dz + &
                             pyramid_source_nx_py_pz_HeII_density_array(-a,b,i_level_py) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
														   											
            enddo		  
          enddo		  

          nx_py_pz_dx(-i_level_py)%HI_density(i_primary,i_secondary) = &
                                  HI_number_dx/pyramid_grid(i_level_py)%volume
          nx_py_pz_dy(i_level_py)%HI_density(i_primary,-i_secondary) = &
                                  HI_number_dy/pyramid_grid(i_level_py)%volume
          nx_py_pz_dz(i_level_py)%HI_density(-i_primary,i_secondary) = &
                                  HI_number_dz/pyramid_grid(i_level_py)%volume
          nx_py_pz_dx(-i_level_py)%HeI_density(i_primary,i_secondary) = &
                                  HeI_number_dx/pyramid_grid(i_level_py)%volume
          nx_py_pz_dy(i_level_py)%HeI_density(i_primary,-i_secondary) = &
                                  HeI_number_dy/pyramid_grid(i_level_py)%volume
          nx_py_pz_dz(i_level_py)%HeI_density(-i_primary,i_secondary) = &
                                  HeI_number_dz/pyramid_grid(i_level_py)%volume
          nx_py_pz_dx(-i_level_py)%HeII_density(i_primary,i_secondary) = &
                                  HeII_number_dx/pyramid_grid(i_level_py)%volume
          nx_py_pz_dy(i_level_py)%HeII_density(i_primary,-i_secondary) = &
                                  HeII_number_dy/pyramid_grid(i_level_py)%volume
          nx_py_pz_dz(i_level_py)%HeII_density(-i_primary,i_secondary) = &
                                  HeII_number_dz/pyramid_grid(i_level_py)%volume

        end do
      end do
    end do

  end subroutine nx_py_pz_species_density_calculation
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nx_py_nz_species_density_calculation ()

    integer :: i_primary,i_secondary,i_level_py,n
    integer :: a,b
    real(kind=dp) :: HI_number_dx,HI_number_dy,HI_number_dz
    real(kind=dp) :: HeI_number_dx,HeI_number_dy,HeI_number_dz
    real(kind=dp) :: HeII_number_dx,HeII_number_dy,HeII_number_dz	
	
    !Distance of the center of pyramid_grid to the source
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)
		  
          HI_number_dx = 0
          HI_number_dy = 0
          HI_number_dz = 0	
          HeI_number_dx = 0
          HeI_number_dy = 0
          HeI_number_dz = 0
          HeII_number_dx = 0
          HeII_number_dy = 0
          HeII_number_dz = 0
		  	  
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
              HeI_number_dx = HeI_number_dx + &
                             pyramid_source_nx_py_nz_HeI_density_array(-i_level_py,a,-b) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeI_number_dy = HeI_number_dy + &
                             pyramid_source_nx_py_nz_HeI_density_array(-b,i_level_py,-a) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)		
              HeI_number_dz = HeI_number_dz + &
                             pyramid_source_nx_py_nz_HeI_density_array(-a,b,-i_level_py) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeII_number_dx = HeII_number_dx + &
                             pyramid_source_nx_py_nz_HeII_density_array(-i_level_py,a,-b) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeII_number_dy = HeII_number_dy + &
                             pyramid_source_nx_py_nz_HeII_density_array(-b,i_level_py,-a) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)		
              HeII_number_dz = HeII_number_dz + &
                             pyramid_source_nx_py_nz_HeII_density_array(-a,b,-i_level_py) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
														   													  
            enddo		  
          enddo		  

          nx_py_nz_dx(-i_level_py)%HI_density(i_primary,-i_secondary) = &
                                  HI_number_dx/pyramid_grid(i_level_py)%volume
          nx_py_nz_dy(i_level_py)%HI_density(-i_primary,-i_secondary) = &
                                  HI_number_dy/pyramid_grid(i_level_py)%volume
          nx_py_nz_dz(-i_level_py)%HI_density(-i_primary,i_secondary) = &
                                  HI_number_dz/pyramid_grid(i_level_py)%volume
          nx_py_nz_dx(-i_level_py)%HeI_density(i_primary,-i_secondary) = &
                                  HeI_number_dx/pyramid_grid(i_level_py)%volume
          nx_py_nz_dy(i_level_py)%HeI_density(-i_primary,-i_secondary) = &
                                  HeI_number_dy/pyramid_grid(i_level_py)%volume
          nx_py_nz_dz(-i_level_py)%HeI_density(-i_primary,i_secondary) = &
                                  HeI_number_dz/pyramid_grid(i_level_py)%volume
          nx_py_nz_dx(-i_level_py)%HeII_density(i_primary,-i_secondary) = &
                                  HeII_number_dx/pyramid_grid(i_level_py)%volume
          nx_py_nz_dy(i_level_py)%HeII_density(-i_primary,-i_secondary) = &
                                  HeII_number_dy/pyramid_grid(i_level_py)%volume
          nx_py_nz_dz(-i_level_py)%HeII_density(-i_primary,i_secondary) = &
                                  HeII_number_dz/pyramid_grid(i_level_py)%volume
		  		  
        end do
      end do
    end do

  end subroutine nx_py_nz_species_density_calculation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nx_ny_pz_species_density_calculation ()
	
    integer :: i_primary,i_secondary,i_level_py,n
    integer :: a,b
    real(kind=dp) :: HI_number_dx,HI_number_dy,HI_number_dz
    real(kind=dp) :: HeI_number_dx,HeI_number_dy,HeI_number_dz
    real(kind=dp) :: HeII_number_dx,HeII_number_dy,HeII_number_dz	
	
    !Distance of the center of pyramid_grid to the source
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)
		  
          HI_number_dx = 0
          HI_number_dy = 0
          HI_number_dz = 0	
          HeI_number_dx = 0
          HeI_number_dy = 0
          HeI_number_dz = 0
          HeII_number_dx = 0
          HeII_number_dy = 0
          HeII_number_dz = 0
		  	  
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
               HeI_number_dx = HeI_number_dx + &
                              pyramid_source_nx_ny_pz_HeI_density_array(-i_level_py,-a,b) * &
                              pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
               HeI_number_dy = HeI_number_dy + &
                              pyramid_source_nx_ny_pz_HeI_density_array(-b,-i_level_py,a) * &
                              pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
               HeI_number_dz = HeI_number_dz + &
                              pyramid_source_nx_ny_pz_HeI_density_array(-a,-b,i_level_py) * &
                              pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
               HeII_number_dx = HeII_number_dx + &
                              pyramid_source_nx_ny_pz_HeII_density_array(-i_level_py,-a,b) * &
                              pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
               HeII_number_dy = HeII_number_dy + &
                              pyramid_source_nx_ny_pz_HeII_density_array(-b,-i_level_py,a) * &
                              pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
               HeII_number_dz = HeII_number_dz + &
                              pyramid_source_nx_ny_pz_HeII_density_array(-a,-b,i_level_py) * &
                              pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
																										 		  
            enddo		  
          enddo		  

          nx_ny_pz_dx(-i_level_py)%HI_density(-i_primary,i_secondary) = &
                                  HI_number_dx/pyramid_grid(i_level_py)%volume
          nx_ny_pz_dy(-i_level_py)%HI_density(i_primary,-i_secondary) = &
                                  HI_number_dy/pyramid_grid(i_level_py)%volume
          nx_ny_pz_dz(i_level_py)%HI_density(-i_primary,-i_secondary) = &
                                  HI_number_dz/pyramid_grid(i_level_py)%volume
          nx_ny_pz_dx(-i_level_py)%HeI_density(-i_primary,i_secondary) = &
                                  HeI_number_dx/pyramid_grid(i_level_py)%volume
          nx_ny_pz_dy(-i_level_py)%HeI_density(i_primary,-i_secondary) = &
                                  HeI_number_dy/pyramid_grid(i_level_py)%volume
          nx_ny_pz_dz(i_level_py)%HeI_density(-i_primary,-i_secondary) = &
                                  HeI_number_dz/pyramid_grid(i_level_py)%volume
          nx_ny_pz_dx(-i_level_py)%HeII_density(-i_primary,i_secondary) = &
                                  HeII_number_dx/pyramid_grid(i_level_py)%volume
          nx_ny_pz_dy(-i_level_py)%HeII_density(i_primary,-i_secondary) = &
                                  HeII_number_dy/pyramid_grid(i_level_py)%volume
          nx_ny_pz_dz(i_level_py)%HeII_density(-i_primary,-i_secondary) = &
                                  HeII_number_dz/pyramid_grid(i_level_py)%volume
		  		  
        end do
      end do
    end do

  end subroutine nx_ny_pz_species_density_calculation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine nx_ny_nz_species_density_calculation ()
	
    integer :: i_primary,i_secondary,i_level_py,n
    integer :: a,b
    real(kind=dp) :: HI_number_dx,HI_number_dy,HI_number_dz
    real(kind=dp) :: HeI_number_dx,HeI_number_dy,HeI_number_dz
    real(kind=dp) :: HeII_number_dx,HeII_number_dy,HeII_number_dz	
	
    !Distance of the center of pyramid_grid to the source
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)
		  
          HI_number_dx = 0
          HI_number_dy = 0
          HI_number_dz = 0	
          HeI_number_dx = 0
          HeI_number_dy = 0
          HeI_number_dz = 0
          HeII_number_dx = 0
          HeII_number_dy = 0
          HeII_number_dz = 0
		  	  
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
              HeI_number_dx = HeI_number_dx + &
                             pyramid_source_nx_ny_nz_HeI_density_array(-i_level_py,-a,-b) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeI_number_dy = HeI_number_dy + &
                             pyramid_source_nx_ny_nz_HeI_density_array(-b,-i_level_py,-a) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeI_number_dz = HeI_number_dz + &
                             pyramid_source_nx_ny_nz_HeI_density_array(-a,-b,-i_level_py) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeII_number_dx = HeII_number_dx + &
                             pyramid_source_nx_ny_nz_HeII_density_array(-i_level_py,-a,-b) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeII_number_dy = HeII_number_dy + &
                             pyramid_source_nx_ny_nz_HeII_density_array(-b,-i_level_py,-a) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
              HeII_number_dz = HeII_number_dz + &
                             pyramid_source_nx_ny_nz_HeII_density_array(-a,-b,-i_level_py) * &
                             pyramid_grid(i_level_py)%transform(i_primary,i_secondary)%transform_volume(a,b)
														   													  
            enddo		  
          enddo		  

          nx_ny_nz_dx(-i_level_py)%HI_density(-i_primary,-i_secondary) = &
                                  HI_number_dx/pyramid_grid(i_level_py)%volume
          nx_ny_nz_dy(-i_level_py)%HI_density(-i_primary,-i_secondary) = &
                                  HI_number_dy/pyramid_grid(i_level_py)%volume
          nx_ny_nz_dz(-i_level_py)%HI_density(-i_primary,-i_secondary) = &
                                  HI_number_dz/pyramid_grid(i_level_py)%volume
          nx_ny_nz_dx(-i_level_py)%HeI_density(-i_primary,-i_secondary) = &
                                  HeI_number_dx/pyramid_grid(i_level_py)%volume
          nx_ny_nz_dy(-i_level_py)%HeI_density(-i_primary,-i_secondary) = &
                                  HeI_number_dy/pyramid_grid(i_level_py)%volume
          nx_ny_nz_dz(-i_level_py)%HeI_density(-i_primary,-i_secondary) = &
                                  HeI_number_dz/pyramid_grid(i_level_py)%volume
          nx_ny_nz_dx(-i_level_py)%HeII_density(-i_primary,-i_secondary) = &
                                  HeII_number_dx/pyramid_grid(i_level_py)%volume
          nx_ny_nz_dy(-i_level_py)%HeII_density(-i_primary,-i_secondary) = &
                                  HeII_number_dy/pyramid_grid(i_level_py)%volume
          nx_ny_nz_dz(-i_level_py)%HeII_density(-i_primary,-i_secondary) = &
                                  HeII_number_dz/pyramid_grid(i_level_py)%volume

        end do
      end do
    end do

  end subroutine nx_ny_nz_species_density_calculation	
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  subroutine source_photo_rate_initialization()
	  
    pyramid_source_px_py_pz_HI_photoionization_rate_array = 0
    pyramid_source_px_py_nz_HI_photoionization_rate_array = 0
    pyramid_source_px_ny_pz_HI_photoionization_rate_array = 0
    pyramid_source_px_ny_nz_HI_photoionization_rate_array = 0
    pyramid_source_nx_py_pz_HI_photoionization_rate_array = 0
    pyramid_source_nx_py_nz_HI_photoionization_rate_array = 0
    pyramid_source_nx_ny_pz_HI_photoionization_rate_array = 0
    pyramid_source_nx_ny_nz_HI_photoionization_rate_array = 0  

    pyramid_source_px_py_pz_HeI_photoionization_rate_array = 0
    pyramid_source_px_py_nz_HeI_photoionization_rate_array = 0
    pyramid_source_px_ny_pz_HeI_photoionization_rate_array = 0
    pyramid_source_px_ny_nz_HeI_photoionization_rate_array = 0
    pyramid_source_nx_py_pz_HeI_photoionization_rate_array = 0
    pyramid_source_nx_py_nz_HeI_photoionization_rate_array = 0
    pyramid_source_nx_ny_pz_HeI_photoionization_rate_array = 0
    pyramid_source_nx_ny_nz_HeI_photoionization_rate_array = 0  
	
    pyramid_source_px_py_pz_HeII_photoionization_rate_array = 0
    pyramid_source_px_py_nz_HeII_photoionization_rate_array = 0
    pyramid_source_px_ny_pz_HeII_photoionization_rate_array = 0
    pyramid_source_px_ny_nz_HeII_photoionization_rate_array = 0
    pyramid_source_nx_py_pz_HeII_photoionization_rate_array = 0
    pyramid_source_nx_py_nz_HeII_photoionization_rate_array = 0
    pyramid_source_nx_ny_pz_HeII_photoionization_rate_array = 0
    pyramid_source_nx_ny_nz_HeII_photoionization_rate_array = 0  
		
    pyramid_source_px_py_pz_photoheating_rate_array = 0
    pyramid_source_px_py_nz_photoheating_rate_array = 0
    pyramid_source_px_ny_pz_photoheating_rate_array = 0
    pyramid_source_px_ny_nz_photoheating_rate_array = 0
    pyramid_source_nx_py_pz_photoheating_rate_array = 0
    pyramid_source_nx_py_nz_photoheating_rate_array = 0
    pyramid_source_nx_ny_pz_photoheating_rate_array = 0
    pyramid_source_nx_ny_nz_photoheating_rate_array = 0 
				  	   
  end subroutine source_photo_rate_initialization
  
  subroutine px_py_pz_photo_rate_calculation()
	
    implicit none

    integer :: i_level_py,a,b
    integer :: i_primary,i_secondary
    real(kind=dp) :: HI_photoionization_rate_dx,HI_photoionization_rate_dy,HI_photoionization_rate_dz
    real(kind=dp) :: HeI_photoionization_rate_dx,HeI_photoionization_rate_dy,HeI_photoionization_rate_dz
    real(kind=dp) :: HeII_photoionization_rate_dx,HeII_photoionization_rate_dy,HeII_photoionization_rate_dz	
    real(kind=dp) :: photoheating_rate_dx,photoheating_rate_dy,photoheating_rate_dz
				
    do i_level_py=1,level_py
      do a=1,i_level_py
        do b=1,i_level_py
			  
          HI_photoionization_rate_dx = 0
          HI_photoionization_rate_dy = 0
          HI_photoionization_rate_dz = 0
          HeI_photoionization_rate_dx = 0
          HeI_photoionization_rate_dy = 0
          HeI_photoionization_rate_dz = 0
          HeII_photoionization_rate_dx = 0
          HeII_photoionization_rate_dy = 0
          HeII_photoionization_rate_dz = 0									
          photoheating_rate_dx = 0
          photoheating_rate_dy = 0
          photoheating_rate_dz = 0
																		
          do i_primary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1),&
                         ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1)	
            do i_secondary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2),&
                             ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2)
 
              if (cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) .gt. 0) then

                if (ionization_weighted_by_atom.eqv..true.) then
                  HI_photoionization_rate_dx = HI_photoionization_rate_dx + &
                                           px_py_pz_dx(i_level_py)%HI_photoionization_rate(i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_py_pz_HI_density_array(i_level_py,a,b)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_py_pz_dx(i_level_py)%HI_density(i_primary,i_secondary))
                  HI_photoionization_rate_dy = HI_photoionization_rate_dy + &
                                           px_py_pz_dy(i_level_py)%HI_photoionization_rate(i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_py_pz_HI_density_array(b,i_level_py,a)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_py_pz_dy(i_level_py)%HI_density(i_primary,i_secondary))
                  HI_photoionization_rate_dz = HI_photoionization_rate_dz + &
                                           px_py_pz_dz(i_level_py)%HI_photoionization_rate(i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_py_pz_HI_density_array(a,b,i_level_py)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_py_pz_dz(i_level_py)%HI_density(i_primary,i_secondary))
                  HeI_photoionization_rate_dx = HeI_photoionization_rate_dx + &
                                           px_py_pz_dx(i_level_py)%HeI_photoionization_rate(i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_py_pz_HeI_density_array(i_level_py,a,b)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_py_pz_dx(i_level_py)%HeI_density(i_primary,i_secondary))
                  HeI_photoionization_rate_dy = HeI_photoionization_rate_dy + &
                                           px_py_pz_dy(i_level_py)%HeI_photoionization_rate(i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_py_pz_HeI_density_array(b,i_level_py,a)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_py_pz_dy(i_level_py)%HeI_density(i_primary,i_secondary))
                  HeI_photoionization_rate_dz = HeI_photoionization_rate_dz + &
                                           px_py_pz_dz(i_level_py)%HeI_photoionization_rate(i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_py_pz_HeI_density_array(a,b,i_level_py)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_py_pz_dz(i_level_py)%HeI_density(i_primary,i_secondary))
                  HeII_photoionization_rate_dx = HeII_photoionization_rate_dx + &
                                           px_py_pz_dx(i_level_py)%HeII_photoionization_rate(i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_py_pz_HeII_density_array(i_level_py,a,b)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_py_pz_dx(i_level_py)%HeII_density(i_primary,i_secondary))
                  HeII_photoionization_rate_dy = HeII_photoionization_rate_dy + &
                                           px_py_pz_dy(i_level_py)%HeII_photoionization_rate(i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_py_pz_HeII_density_array(b,i_level_py,a)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_py_pz_dy(i_level_py)%HeII_density(i_primary,i_secondary))
                  HeII_photoionization_rate_dz = HeII_photoionization_rate_dz + &
                                           px_py_pz_dz(i_level_py)%HeII_photoionization_rate(i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_py_pz_HeII_density_array(a,b,i_level_py)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_py_pz_dz(i_level_py)%HeII_density(i_primary,i_secondary))

                else
                  HI_photoionization_rate_dx = HI_photoionization_rate_dx + &
                                           px_py_pz_dx(i_level_py)%HI_photoionization_rate(i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                           pyramid_grid(i_level_py)%volume
                  HI_photoionization_rate_dy = HI_photoionization_rate_dy + &
                                           px_py_pz_dy(i_level_py)%HI_photoionization_rate(i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                           pyramid_grid(i_level_py)%volume
                  HI_photoionization_rate_dz = HI_photoionization_rate_dz + &
                                           px_py_pz_dz(i_level_py)%HI_photoionization_rate(i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                           pyramid_grid(i_level_py)%volume
                  HeI_photoionization_rate_dx = HeI_photoionization_rate_dx + &
                                            px_py_pz_dx(i_level_py)%HeI_photoionization_rate(i_primary,i_secondary) * &
                                            cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                            pyramid_grid(i_level_py)%volume
                  HeI_photoionization_rate_dy = HeI_photoionization_rate_dy + &
                                            px_py_pz_dy(i_level_py)%HeI_photoionization_rate(i_primary,i_secondary) * &
                                            cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                            pyramid_grid(i_level_py)%volume
                  HeI_photoionization_rate_dz = HeI_photoionization_rate_dz + &
                                            px_py_pz_dz(i_level_py)%HeI_photoionization_rate(i_primary,i_secondary) * &
                                            cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                            pyramid_grid(i_level_py)%volume
                  HeII_photoionization_rate_dx = HeII_photoionization_rate_dx + &
                                             px_py_pz_dx(i_level_py)%HeII_photoionization_rate(i_primary,i_secondary) * &
                                             cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                             pyramid_grid(i_level_py)%volume
                  HeII_photoionization_rate_dy = HeII_photoionization_rate_dy + &
                                             px_py_pz_dy(i_level_py)%HeII_photoionization_rate(i_primary,i_secondary) * &
                                             cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                             pyramid_grid(i_level_py)%volume
                  HeII_photoionization_rate_dz = HeII_photoionization_rate_dz + &
                                             px_py_pz_dz(i_level_py)%HeII_photoionization_rate(i_primary,i_secondary) * &
                                             cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                             pyramid_grid(i_level_py)%volume

                endif

                if (heating_weighted_by_atom.eqv..true.) then

                  photoheating_rate_dx = photoheating_rate_dx + &
                                     px_py_pz_dx(i_level_py)%photoheating_rate(i_primary,i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                     pyramid_source_px_py_pz_HI_density_array(i_level_py,a,b)/ &
                                     (pyramid_grid(i_level_py)%volume * &
                                     px_py_pz_dx(i_level_py)%HI_density(i_primary,i_secondary))
                  photoheating_rate_dy = photoheating_rate_dy + &
                                     px_py_pz_dy(i_level_py)%photoheating_rate(i_primary,i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                     pyramid_source_px_py_pz_HI_density_array(b,i_level_py,a)/ &
                                     (pyramid_grid(i_level_py)%volume * &
                                     px_py_pz_dy(i_level_py)%HI_density(i_primary,i_secondary))
                  photoheating_rate_dz = photoheating_rate_dz + &
                                     px_py_pz_dz(i_level_py)%photoheating_rate(i_primary,i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                     pyramid_source_px_py_pz_HI_density_array(a,b,i_level_py)/ &
                                     (pyramid_grid(i_level_py)%volume * &
                                     px_py_pz_dz(i_level_py)%HI_density(i_primary,i_secondary))
                else

                  photoheating_rate_dx = photoheating_rate_dx + &
                                     px_py_pz_dx(i_level_py)%photoheating_rate(i_primary,i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                     pyramid_grid(i_level_py)%volume
                  photoheating_rate_dy = photoheating_rate_dy + &
                                     px_py_pz_dy(i_level_py)%photoheating_rate(i_primary,i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                     pyramid_grid(i_level_py)%volume
                  photoheating_rate_dz = photoheating_rate_dz + &
                                     px_py_pz_dz(i_level_py)%photoheating_rate(i_primary,i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                     pyramid_grid(i_level_py)%volume

                endif
				
              endif

            enddo
          enddo
		  
          pyramid_source_px_py_pz_HI_photoionization_rate_array(i_level_py,a,b) = &
                          pyramid_source_px_py_pz_HI_photoionization_rate_array(i_level_py,a,b)+HI_photoionization_rate_dx
          pyramid_source_px_py_pz_HI_photoionization_rate_array(b,i_level_py,a) = &
                          pyramid_source_px_py_pz_HI_photoionization_rate_array(b,i_level_py,a)+HI_photoionization_rate_dy
          pyramid_source_px_py_pz_HI_photoionization_rate_array(a,b,i_level_py) = &
                          pyramid_source_px_py_pz_HI_photoionization_rate_array(a,b,i_level_py)+HI_photoionization_rate_dz
          pyramid_source_px_py_pz_HeI_photoionization_rate_array(i_level_py,a,b) = &
                          pyramid_source_px_py_pz_HeI_photoionization_rate_array(i_level_py,a,b)+HeI_photoionization_rate_dx
          pyramid_source_px_py_pz_HeI_photoionization_rate_array(b,i_level_py,a) = &
                          pyramid_source_px_py_pz_HeI_photoionization_rate_array(b,i_level_py,a)+HeI_photoionization_rate_dy
          pyramid_source_px_py_pz_HeI_photoionization_rate_array(a,b,i_level_py) = &
                          pyramid_source_px_py_pz_HeI_photoionization_rate_array(a,b,i_level_py)+HeI_photoionization_rate_dz
          pyramid_source_px_py_pz_HeII_photoionization_rate_array(i_level_py,a,b) = &
                          pyramid_source_px_py_pz_HeII_photoionization_rate_array(i_level_py,a,b)+HeII_photoionization_rate_dx
          pyramid_source_px_py_pz_HeII_photoionization_rate_array(b,i_level_py,a) = &
                          pyramid_source_px_py_pz_HeII_photoionization_rate_array(b,i_level_py,a)+HeII_photoionization_rate_dy
          pyramid_source_px_py_pz_HeII_photoionization_rate_array(a,b,i_level_py) = &
                          pyramid_source_px_py_pz_HeII_photoionization_rate_array(a,b,i_level_py)+HeII_photoionization_rate_dz
          pyramid_source_px_py_pz_photoheating_rate_array(i_level_py,a,b) = &
                          pyramid_source_px_py_pz_photoheating_rate_array(i_level_py,a,b)+photoheating_rate_dx
          pyramid_source_px_py_pz_photoheating_rate_array(b,i_level_py,a) = &
                          pyramid_source_px_py_pz_photoheating_rate_array(b,i_level_py,a)+photoheating_rate_dy
          pyramid_source_px_py_pz_photoheating_rate_array(a,b,i_level_py) = &
                          pyramid_source_px_py_pz_photoheating_rate_array(a,b,i_level_py)+photoheating_rate_dz

        enddo
      enddo		
    enddo

  end subroutine px_py_pz_photo_rate_calculation

  subroutine px_py_nz_photo_rate_calculation()
	
    integer :: i_level_py,a,b
    integer :: i_primary,i_secondary
    real(kind=dp) :: HI_photoionization_rate_dx,HI_photoionization_rate_dy,HI_photoionization_rate_dz
    real(kind=dp) :: HeI_photoionization_rate_dx,HeI_photoionization_rate_dy,HeI_photoionization_rate_dz
    real(kind=dp) :: HeII_photoionization_rate_dx,HeII_photoionization_rate_dy,HeII_photoionization_rate_dz	
    real(kind=dp) :: photoheating_rate_dx,photoheating_rate_dy,photoheating_rate_dz
		
    do i_level_py=1,level_py
      do a=1,i_level_py
        do b=1,i_level_py
			  
          HI_photoionization_rate_dx = 0
          HI_photoionization_rate_dy = 0
          HI_photoionization_rate_dz = 0
          HeI_photoionization_rate_dx = 0
          HeI_photoionization_rate_dy = 0
          HeI_photoionization_rate_dz = 0
          HeII_photoionization_rate_dx = 0
          HeII_photoionization_rate_dy = 0
          HeII_photoionization_rate_dz = 0
          photoheating_rate_dx = 0
          photoheating_rate_dy = 0
          photoheating_rate_dz = 0
									
          do i_primary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1),&
  	                     ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1)	
            do i_secondary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2),&
                             ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2)
 
              if (cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) .gt. 0) then

                if (ionization_weighted_by_atom.eqv..true.) then

                  HI_photoionization_rate_dx = HI_photoionization_rate_dx + &
                                           px_py_nz_dx(i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_py_nz_HI_density_array(i_level_py,a,-b)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_py_nz_dx(i_level_py)%HI_density(i_primary,-i_secondary))
                  HI_photoionization_rate_dy = HI_photoionization_rate_dy + &
                                           px_py_nz_dy(i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_py_nz_HI_density_array(b,i_level_py,-a)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_py_nz_dy(i_level_py)%HI_density(-i_primary,i_secondary))
                  HI_photoionization_rate_dz = HI_photoionization_rate_dz + &
                                           px_py_nz_dz(-i_level_py)%HI_photoionization_rate(i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_py_nz_HI_density_array(a,b,-i_level_py)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_py_nz_dz(-i_level_py)%HI_density(i_primary,i_secondary))
                  HeI_photoionization_rate_dx = HeI_photoionization_rate_dx + &
                                           px_py_nz_dx(i_level_py)%HeI_photoionization_rate(i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_py_nz_HeI_density_array(i_level_py,a,-b)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_py_nz_dx(i_level_py)%HeI_density(i_primary,-i_secondary))
                  HeI_photoionization_rate_dy = HeI_photoionization_rate_dy + &
                                           px_py_nz_dy(i_level_py)%HeI_photoionization_rate(-i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_py_nz_HeI_density_array(b,i_level_py,-a)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_py_nz_dy(i_level_py)%HeI_density(-i_primary,i_secondary))
                  HeI_photoionization_rate_dz = HeI_photoionization_rate_dz + &
                                           px_py_nz_dz(-i_level_py)%HeI_photoionization_rate(i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_py_nz_HeI_density_array(a,b,-i_level_py)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_py_nz_dz(-i_level_py)%HeI_density(i_primary,i_secondary))
                  HeII_photoionization_rate_dx = HeII_photoionization_rate_dx + &
                                           px_py_nz_dx(i_level_py)%HeII_photoionization_rate(i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_py_nz_HeII_density_array(i_level_py,a,-b)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_py_nz_dx(i_level_py)%HeII_density(i_primary,-i_secondary))
                  HeII_photoionization_rate_dy = HeII_photoionization_rate_dy + &
                                           px_py_nz_dy(i_level_py)%HeII_photoionization_rate(-i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_py_nz_HeII_density_array(b,i_level_py,-a)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_py_nz_dy(i_level_py)%HeII_density(-i_primary,i_secondary))
                  HeII_photoionization_rate_dz = HeII_photoionization_rate_dz + &
                                           px_py_nz_dz(-i_level_py)%HeII_photoionization_rate(i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_py_nz_HeII_density_array(a,b,-i_level_py)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_py_nz_dz(-i_level_py)%HeII_density(i_primary,i_secondary))


                else

                  HI_photoionization_rate_dx = HI_photoionization_rate_dx + &
			                  px_py_nz_dx(i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
                  HI_photoionization_rate_dy = HI_photoionization_rate_dy + &
			                  px_py_nz_dy(i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume							  
                  HI_photoionization_rate_dz = HI_photoionization_rate_dz + &
			                  px_py_nz_dz(-i_level_py)%HI_photoionization_rate(i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume						  
                  HeI_photoionization_rate_dx = HeI_photoionization_rate_dx + &
			                  px_py_nz_dx(i_level_py)%HeI_photoionization_rate(i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
                  HeI_photoionization_rate_dy = HeI_photoionization_rate_dy + &
			                  px_py_nz_dy(i_level_py)%HeI_photoionization_rate(-i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume							  
                  HeI_photoionization_rate_dz = HeI_photoionization_rate_dz + &
			                  px_py_nz_dz(-i_level_py)%HeI_photoionization_rate(i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
                  HeII_photoionization_rate_dx = HeII_photoionization_rate_dx + &
			                  px_py_nz_dx(i_level_py)%HeII_photoionization_rate(i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
                  HeII_photoionization_rate_dy = HeII_photoionization_rate_dy + &
			                  px_py_nz_dy(i_level_py)%HeII_photoionization_rate(-i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume							  
                  HeII_photoionization_rate_dz = HeII_photoionization_rate_dz + &
			                  px_py_nz_dz(-i_level_py)%HeII_photoionization_rate(i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
					
                endif
						               
                if (heating_weighted_by_atom.eqv..true.) then

                  photoheating_rate_dx = photoheating_rate_dx + &
                                     px_py_nz_dx(i_level_py)%photoheating_rate(i_primary,-i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                     pyramid_source_px_py_nz_HI_density_array(i_level_py,a,-b)/ &
                                     (pyramid_grid(i_level_py)%volume * &
                                     px_py_nz_dx(i_level_py)%HI_density(i_primary,-i_secondary))
                  photoheating_rate_dy = photoheating_rate_dy + &
                                     px_py_nz_dy(i_level_py)%photoheating_rate(-i_primary,i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                     pyramid_source_px_py_nz_HI_density_array(b,i_level_py,-a)/ &
                                     (pyramid_grid(i_level_py)%volume * &
                                     px_py_nz_dy(i_level_py)%HI_density(-i_primary,i_secondary))
                  photoheating_rate_dz = photoheating_rate_dz + &
                                     px_py_nz_dz(-i_level_py)%photoheating_rate(i_primary,i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                     pyramid_source_px_py_nz_HI_density_array(a,b,-i_level_py)/ &
                                     (pyramid_grid(i_level_py)%volume * &
                                     px_py_nz_dz(-i_level_py)%HI_density(i_primary,i_secondary))
                else

                  photoheating_rate_dx = photoheating_rate_dx + &
                                     px_py_nz_dx(i_level_py)%photoheating_rate(i_primary,-i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                     pyramid_grid(i_level_py)%volume
                  photoheating_rate_dy = photoheating_rate_dy + &
                                     px_py_nz_dy(i_level_py)%photoheating_rate(-i_primary,i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                     pyramid_grid(i_level_py)%volume
                  photoheating_rate_dz = photoheating_rate_dz + &
                                     px_py_nz_dz(-i_level_py)%photoheating_rate(i_primary,i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                     pyramid_grid(i_level_py)%volume	

                endif					  							  			  
              endif				  
            enddo
          enddo
		  
          pyramid_source_px_py_nz_HI_photoionization_rate_array(i_level_py,a,-b) = &
		             pyramid_source_px_py_nz_HI_photoionization_rate_array(i_level_py,a,-b)+HI_photoionization_rate_dx
          pyramid_source_px_py_nz_HI_photoionization_rate_array(b,i_level_py,-a) = &
		             pyramid_source_px_py_nz_HI_photoionization_rate_array(b,i_level_py,-a)+HI_photoionization_rate_dy
          pyramid_source_px_py_nz_HI_photoionization_rate_array(a,b,-i_level_py) = &
		             pyramid_source_px_py_nz_HI_photoionization_rate_array(a,b,-i_level_py)+HI_photoionization_rate_dz
          pyramid_source_px_py_nz_HeI_photoionization_rate_array(i_level_py,a,-b) = &
		             pyramid_source_px_py_nz_HeI_photoionization_rate_array(i_level_py,a,-b)+HeI_photoionization_rate_dx
          pyramid_source_px_py_nz_HeI_photoionization_rate_array(b,i_level_py,-a) = &
		             pyramid_source_px_py_nz_HeI_photoionization_rate_array(b,i_level_py,-a)+HeI_photoionization_rate_dy
          pyramid_source_px_py_nz_HeI_photoionization_rate_array(a,b,-i_level_py) = &
		             pyramid_source_px_py_nz_HeI_photoionization_rate_array(a,b,-i_level_py)+HeI_photoionization_rate_dz
          pyramid_source_px_py_nz_HeII_photoionization_rate_array(i_level_py,a,-b) = &
		             pyramid_source_px_py_nz_HeII_photoionization_rate_array(i_level_py,a,-b)+HeII_photoionization_rate_dx
          pyramid_source_px_py_nz_HeII_photoionization_rate_array(b,i_level_py,-a) = &
		             pyramid_source_px_py_nz_HeII_photoionization_rate_array(b,i_level_py,-a)+HeII_photoionization_rate_dy
          pyramid_source_px_py_nz_HeII_photoionization_rate_array(a,b,-i_level_py) = &
		             pyramid_source_px_py_nz_HeII_photoionization_rate_array(a,b,-i_level_py)+HeII_photoionization_rate_dz
          pyramid_source_px_py_nz_photoheating_rate_array(i_level_py,a,-b) = &
		             pyramid_source_px_py_nz_photoheating_rate_array(i_level_py,a,-b)+photoheating_rate_dx
          pyramid_source_px_py_nz_photoheating_rate_array(b,i_level_py,-a) = &
		             pyramid_source_px_py_nz_photoheating_rate_array(b,i_level_py,-a)+photoheating_rate_dy
          pyramid_source_px_py_nz_photoheating_rate_array(a,b,-i_level_py) = &
		             pyramid_source_px_py_nz_photoheating_rate_array(a,b,-i_level_py)+photoheating_rate_dz

        enddo
      enddo		
    enddo
 
  end subroutine px_py_nz_photo_rate_calculation

  subroutine px_ny_pz_photo_rate_calculation()
	
    integer :: i_level_py,a,b
    integer :: i_primary,i_secondary
    real(kind=dp) :: HI_photoionization_rate_dx,HI_photoionization_rate_dy,HI_photoionization_rate_dz
    real(kind=dp) :: HeI_photoionization_rate_dx,HeI_photoionization_rate_dy,HeI_photoionization_rate_dz
    real(kind=dp) :: HeII_photoionization_rate_dx,HeII_photoionization_rate_dy,HeII_photoionization_rate_dz	
    real(kind=dp) :: photoheating_rate_dx,photoheating_rate_dy,photoheating_rate_dz
		
    do i_level_py=1,level_py
      do a=1,i_level_py
        do b=1,i_level_py
			  
          HI_photoionization_rate_dx = 0
          HI_photoionization_rate_dy = 0
          HI_photoionization_rate_dz = 0
          HeI_photoionization_rate_dx = 0
          HeI_photoionization_rate_dy = 0
          HeI_photoionization_rate_dz = 0
          HeII_photoionization_rate_dx = 0
          HeII_photoionization_rate_dy = 0
          HeII_photoionization_rate_dz = 0
          photoheating_rate_dx = 0
          photoheating_rate_dy = 0
          photoheating_rate_dz = 0
									
          do i_primary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1),&
  	                     ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1)	
            do i_secondary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2),&
                             ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2)
 
              if (cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) .gt. 0) then

                if (ionization_weighted_by_atom.eqv..true.) then

                  HI_photoionization_rate_dx = HI_photoionization_rate_dx + &
                                           px_ny_pz_dx(i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_ny_pz_HI_density_array(i_level_py,-a,b)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_ny_pz_dx(i_level_py)%HI_density(-i_primary,i_secondary))
                  HI_photoionization_rate_dy = HI_photoionization_rate_dy + &
                                           px_ny_pz_dy(-i_level_py)%HI_photoionization_rate(i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_ny_pz_HI_density_array(b,-i_level_py,a)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_ny_pz_dy(-i_level_py)%HI_density(i_primary,i_secondary))
                  HI_photoionization_rate_dz = HI_photoionization_rate_dz + &
                                           px_ny_pz_dz(i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_ny_pz_HI_density_array(a,-b,i_level_py)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_ny_pz_dz(i_level_py)%HI_density(i_primary,-i_secondary))
                  HeI_photoionization_rate_dx = HeI_photoionization_rate_dx + &
                                           px_ny_pz_dx(i_level_py)%HeI_photoionization_rate(-i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_ny_pz_HeI_density_array(i_level_py,-a,b)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_ny_pz_dx(i_level_py)%HeI_density(-i_primary,i_secondary))
                  HeI_photoionization_rate_dy = HeI_photoionization_rate_dy + &
                                           px_ny_pz_dy(-i_level_py)%HeI_photoionization_rate(i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_ny_pz_HeI_density_array(b,-i_level_py,a)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_ny_pz_dy(-i_level_py)%HeI_density(i_primary,i_secondary))
                  HeI_photoionization_rate_dz = HeI_photoionization_rate_dz + &
                                           px_ny_pz_dz(i_level_py)%HeI_photoionization_rate(i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_ny_pz_HeI_density_array(a,-b,i_level_py)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_ny_pz_dz(i_level_py)%HeI_density(i_primary,-i_secondary))
                  HeII_photoionization_rate_dx = HeII_photoionization_rate_dx + &
                                           px_ny_pz_dx(i_level_py)%HeII_photoionization_rate(-i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_ny_pz_HeII_density_array(i_level_py,-a,b)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_ny_pz_dx(i_level_py)%HeII_density(-i_primary,i_secondary))
                  HeII_photoionization_rate_dy = HeII_photoionization_rate_dy + &
                                           px_ny_pz_dy(-i_level_py)%HeII_photoionization_rate(i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_ny_pz_HeII_density_array(b,-i_level_py,a)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_ny_pz_dy(-i_level_py)%HeII_density(i_primary,i_secondary))
                  HeII_photoionization_rate_dz = HeII_photoionization_rate_dz + &
                                           px_ny_pz_dz(i_level_py)%HeII_photoionization_rate(i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_ny_pz_HeII_density_array(a,-b,i_level_py)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_ny_pz_dz(i_level_py)%HeII_density(i_primary,-i_secondary))

                else

                  HI_photoionization_rate_dx = HI_photoionization_rate_dx + &
			                  px_ny_pz_dx(i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
                  HI_photoionization_rate_dy = HI_photoionization_rate_dy + &
			                  px_ny_pz_dy(-i_level_py)%HI_photoionization_rate(i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume							  
                  HI_photoionization_rate_dz = HI_photoionization_rate_dz + &
			                  px_ny_pz_dz(i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume						  
                  HeI_photoionization_rate_dx = HeI_photoionization_rate_dx + &
			                  px_ny_pz_dx(i_level_py)%HeI_photoionization_rate(-i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
                  HeI_photoionization_rate_dy = HeI_photoionization_rate_dy + &
			                  px_ny_pz_dy(-i_level_py)%HeI_photoionization_rate(i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume							  
                  HeI_photoionization_rate_dz = HeI_photoionization_rate_dz + &
			                  px_ny_pz_dz(i_level_py)%HeI_photoionization_rate(i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume	
                  HeII_photoionization_rate_dx = HeII_photoionization_rate_dx + &
			                  px_ny_pz_dx(i_level_py)%HeII_photoionization_rate(-i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
                  HeII_photoionization_rate_dy = HeII_photoionization_rate_dy + &
			                  px_ny_pz_dy(-i_level_py)%HeII_photoionization_rate(i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume							  
                  HeII_photoionization_rate_dz = HeII_photoionization_rate_dz + &
			                  px_ny_pz_dz(i_level_py)%HeII_photoionization_rate(i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
					
                endif

                if (heating_weighted_by_atom.eqv..true.) then

                  photoheating_rate_dx = photoheating_rate_dx + &
                                     px_ny_pz_dx(i_level_py)%photoheating_rate(-i_primary,i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                     pyramid_source_px_ny_pz_HI_density_array(i_level_py,-a,b)/ &
                                     (pyramid_grid(i_level_py)%volume * &
                                     px_ny_pz_dx(i_level_py)%HI_density(-i_primary,i_secondary))
                  photoheating_rate_dy = photoheating_rate_dy + &
                                     px_ny_pz_dy(-i_level_py)%photoheating_rate(i_primary,i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                     pyramid_source_px_ny_pz_HI_density_array(b,-i_level_py,a)/ &
                                     (pyramid_grid(i_level_py)%volume * &
                                     px_ny_pz_dy(-i_level_py)%HI_density(i_primary,i_secondary))
                  photoheating_rate_dz = photoheating_rate_dz + &
                                     px_ny_pz_dz(i_level_py)%photoheating_rate(i_primary,-i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                     pyramid_source_px_ny_pz_HI_density_array(a,-b,i_level_py)/ &
                                     (pyramid_grid(i_level_py)%volume * &
                                     px_ny_pz_dz(i_level_py)%HI_density(i_primary,-i_secondary))
                else

                  photoheating_rate_dx = photoheating_rate_dx + &
                                     px_ny_pz_dx(i_level_py)%photoheating_rate(-i_primary,i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                     pyramid_grid(i_level_py)%volume
                  photoheating_rate_dy = photoheating_rate_dy + &
                                     px_ny_pz_dy(-i_level_py)%photoheating_rate(i_primary,i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                     pyramid_grid(i_level_py)%volume
                  photoheating_rate_dz = photoheating_rate_dz + &
                                     px_ny_pz_dz(i_level_py)%photoheating_rate(i_primary,-i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                     pyramid_grid(i_level_py)%volume

                endif					  	
              endif
            enddo
          enddo
		  
          pyramid_source_px_ny_pz_HI_photoionization_rate_array(i_level_py,-a,b) = &
		             pyramid_source_px_ny_pz_HI_photoionization_rate_array(i_level_py,-a,b)+HI_photoionization_rate_dx
          pyramid_source_px_ny_pz_HI_photoionization_rate_array(b,-i_level_py,a) = &
		             pyramid_source_px_ny_pz_HI_photoionization_rate_array(b,-i_level_py,a)+HI_photoionization_rate_dy
          pyramid_source_px_ny_pz_HI_photoionization_rate_array(a,-b,i_level_py) = &
		             pyramid_source_px_ny_pz_HI_photoionization_rate_array(a,-b,i_level_py)+HI_photoionization_rate_dz
          pyramid_source_px_ny_pz_HeI_photoionization_rate_array(i_level_py,-a,b) = &
		             pyramid_source_px_ny_pz_HeI_photoionization_rate_array(i_level_py,-a,b)+HeI_photoionization_rate_dx
          pyramid_source_px_ny_pz_HeI_photoionization_rate_array(b,-i_level_py,a) = &
		             pyramid_source_px_ny_pz_HeI_photoionization_rate_array(b,-i_level_py,a)+HeI_photoionization_rate_dy
          pyramid_source_px_ny_pz_HeI_photoionization_rate_array(a,-b,i_level_py) = &
		             pyramid_source_px_ny_pz_HeI_photoionization_rate_array(a,-b,i_level_py)+HeI_photoionization_rate_dz
          pyramid_source_px_ny_pz_HeII_photoionization_rate_array(i_level_py,-a,b) = &
		             pyramid_source_px_ny_pz_HeII_photoionization_rate_array(i_level_py,-a,b)+HeII_photoionization_rate_dx
          pyramid_source_px_ny_pz_HeII_photoionization_rate_array(b,-i_level_py,a) = &
		             pyramid_source_px_ny_pz_HeII_photoionization_rate_array(b,-i_level_py,a)+HeII_photoionization_rate_dy
          pyramid_source_px_ny_pz_HeII_photoionization_rate_array(a,-b,i_level_py) = &
		             pyramid_source_px_ny_pz_HeII_photoionization_rate_array(a,-b,i_level_py)+HeII_photoionization_rate_dz
          pyramid_source_px_ny_pz_photoheating_rate_array(i_level_py,-a,b) = &
		             pyramid_source_px_ny_pz_photoheating_rate_array(i_level_py,-a,b)+photoheating_rate_dx
          pyramid_source_px_ny_pz_photoheating_rate_array(b,-i_level_py,a) = &
		             pyramid_source_px_ny_pz_photoheating_rate_array(b,-i_level_py,a)+photoheating_rate_dy
          pyramid_source_px_ny_pz_photoheating_rate_array(a,-b,i_level_py) = &
		             pyramid_source_px_ny_pz_photoheating_rate_array(a,-b,i_level_py)+photoheating_rate_dz
   
        enddo
      enddo		
    enddo

  end subroutine px_ny_pz_photo_rate_calculation

  subroutine px_ny_nz_photo_rate_calculation()
	
    integer :: i_level_py,a,b
    integer :: i_primary,i_secondary
    real(kind=dp) :: HI_photoionization_rate_dx,HI_photoionization_rate_dy,HI_photoionization_rate_dz
    real(kind=dp) :: HeI_photoionization_rate_dx,HeI_photoionization_rate_dy,HeI_photoionization_rate_dz
    real(kind=dp) :: HeII_photoionization_rate_dx,HeII_photoionization_rate_dy,HeII_photoionization_rate_dz	
    real(kind=dp) :: photoheating_rate_dx,photoheating_rate_dy,photoheating_rate_dz
		
    do i_level_py=1,level_py
	  do a=1,i_level_py
        do b=1,i_level_py
			  
          HI_photoionization_rate_dx = 0
          HI_photoionization_rate_dy = 0
          HI_photoionization_rate_dz = 0
          HeI_photoionization_rate_dx = 0
          HeI_photoionization_rate_dy = 0
          HeI_photoionization_rate_dz = 0
          HeII_photoionization_rate_dx = 0
          HeII_photoionization_rate_dy = 0
          HeII_photoionization_rate_dz = 0
          photoheating_rate_dx = 0
          photoheating_rate_dy = 0
          photoheating_rate_dz = 0
									
          do i_primary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1),&
  	                     ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1)	
            do i_secondary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2),&
                             ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2)
 
if (cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) .gt. 0) then

if (ionization_weighted_by_atom.eqv..true.) then

              HI_photoionization_rate_dx = HI_photoionization_rate_dx + &
                                           px_ny_nz_dx(i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_ny_nz_HI_density_array(i_level_py,-a,-b)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_ny_nz_dx(i_level_py)%HI_density(-i_primary,-i_secondary))
              HI_photoionization_rate_dy = HI_photoionization_rate_dy + &
                                           px_ny_nz_dy(-i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_ny_nz_HI_density_array(b,-i_level_py,-a)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_ny_nz_dy(-i_level_py)%HI_density(-i_primary,i_secondary))
              HI_photoionization_rate_dz = HI_photoionization_rate_dz + &
                                           px_ny_nz_dz(-i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_ny_nz_HI_density_array(a,-b,-i_level_py)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_ny_nz_dz(-i_level_py)%HI_density(i_primary,-i_secondary))
              HeI_photoionization_rate_dx = HeI_photoionization_rate_dx + &
                                           px_ny_nz_dx(i_level_py)%HeI_photoionization_rate(-i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_ny_nz_HeI_density_array(i_level_py,-a,-b)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_ny_nz_dx(i_level_py)%HeI_density(-i_primary,-i_secondary))
              HeI_photoionization_rate_dy = HeI_photoionization_rate_dy + &
                                           px_ny_nz_dy(-i_level_py)%HeI_photoionization_rate(-i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_ny_nz_HeI_density_array(b,-i_level_py,-a)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_ny_nz_dy(-i_level_py)%HeI_density(-i_primary,i_secondary))
              HeI_photoionization_rate_dz = HeI_photoionization_rate_dz + &
                                           px_ny_nz_dz(-i_level_py)%HeI_photoionization_rate(i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_ny_nz_HeI_density_array(a,-b,-i_level_py)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_ny_nz_dz(-i_level_py)%HeI_density(i_primary,-i_secondary))
              HeII_photoionization_rate_dx = HeII_photoionization_rate_dx + &
                                           px_ny_nz_dx(i_level_py)%HeII_photoionization_rate(-i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_ny_nz_HeII_density_array(i_level_py,-a,-b)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_ny_nz_dx(i_level_py)%HeII_density(-i_primary,-i_secondary))
              HeII_photoionization_rate_dy = HeII_photoionization_rate_dy + &
                                           px_ny_nz_dy(-i_level_py)%HeII_photoionization_rate(-i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_ny_nz_HeII_density_array(b,-i_level_py,-a)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_ny_nz_dy(-i_level_py)%HeII_density(-i_primary,i_secondary))
              HeII_photoionization_rate_dz = HeII_photoionization_rate_dz + &
                                           px_ny_nz_dz(-i_level_py)%HeII_photoionization_rate(i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_px_ny_nz_HeII_density_array(a,-b,-i_level_py)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           px_ny_nz_dz(-i_level_py)%HeII_density(i_primary,-i_secondary))


else

              HI_photoionization_rate_dx = HI_photoionization_rate_dx + &
			                  px_ny_nz_dx(i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) * &
                              cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                              pyramid_grid(i_level_py)%volume
              HI_photoionization_rate_dy = HI_photoionization_rate_dy + &
			                  px_ny_nz_dy(-i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume							  
              HI_photoionization_rate_dz = HI_photoionization_rate_dz + &
			                  px_ny_nz_dz(-i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume						  
              HeI_photoionization_rate_dx = HeI_photoionization_rate_dx + &
			                  px_ny_nz_dx(i_level_py)%HeI_photoionization_rate(-i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
              HeI_photoionization_rate_dy = HeI_photoionization_rate_dy + &
			                  px_ny_nz_dy(-i_level_py)%HeI_photoionization_rate(-i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume							  
              HeI_photoionization_rate_dz = HeI_photoionization_rate_dz + &
			                  px_ny_nz_dz(-i_level_py)%HeI_photoionization_rate(i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume	
              HeII_photoionization_rate_dx = HeII_photoionization_rate_dx + &
			                  px_ny_nz_dx(i_level_py)%HeII_photoionization_rate(-i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
              HeII_photoionization_rate_dy = HeII_photoionization_rate_dy + &
			                  px_ny_nz_dy(-i_level_py)%HeII_photoionization_rate(-i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume							  
              HeII_photoionization_rate_dz = HeII_photoionization_rate_dz + &
			                  px_ny_nz_dz(-i_level_py)%HeII_photoionization_rate(i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
 					
endif

if (heating_weighted_by_atom.eqv..true.) then

              photoheating_rate_dx = photoheating_rate_dx + &
                                     px_ny_nz_dx(i_level_py)%photoheating_rate(-i_primary,-i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                     pyramid_source_px_ny_nz_HI_density_array(i_level_py,-a,-b)/ &
                                     (pyramid_grid(i_level_py)%volume * &
                                     px_ny_nz_dx(i_level_py)%HI_density(-i_primary,-i_secondary))
              photoheating_rate_dy = photoheating_rate_dy + &
                                     px_ny_nz_dy(-i_level_py)%photoheating_rate(-i_primary,i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                     pyramid_source_px_ny_nz_HI_density_array(b,-i_level_py,-a)/ &
                                     (pyramid_grid(i_level_py)%volume * &
                                     px_ny_nz_dy(-i_level_py)%HI_density(-i_primary,i_secondary))
              photoheating_rate_dz = photoheating_rate_dz + &
                                     px_ny_nz_dz(-i_level_py)%photoheating_rate(i_primary,-i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                     pyramid_source_px_ny_nz_HI_density_array(a,-b,-i_level_py)/ &
                                     (pyramid_grid(i_level_py)%volume * &
                                     px_ny_nz_dz(-i_level_py)%HI_density(i_primary,-i_secondary))


else

              photoheating_rate_dx = photoheating_rate_dx + &
                                     px_ny_nz_dx(i_level_py)%photoheating_rate(-i_primary,-i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                     pyramid_grid(i_level_py)%volume
              photoheating_rate_dy = photoheating_rate_dy + &
                                     px_ny_nz_dy(-i_level_py)%photoheating_rate(-i_primary,i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                     pyramid_grid(i_level_py)%volume
              photoheating_rate_dz = photoheating_rate_dz + &
                                     px_ny_nz_dz(-i_level_py)%photoheating_rate(i_primary,-i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                     pyramid_grid(i_level_py)%volume



endif
						
endif							  				  				  						  
            enddo
          enddo
		  
		  pyramid_source_px_ny_nz_HI_photoionization_rate_array(i_level_py,-a,-b) = &
		             pyramid_source_px_ny_nz_HI_photoionization_rate_array(i_level_py,-a,-b)+HI_photoionization_rate_dx
          pyramid_source_px_ny_nz_HI_photoionization_rate_array(b,-i_level_py,-a) = &
		             pyramid_source_px_ny_nz_HI_photoionization_rate_array(b,-i_level_py,-a)+HI_photoionization_rate_dy
          pyramid_source_px_ny_nz_HI_photoionization_rate_array(a,-b,-i_level_py) = &
		             pyramid_source_px_ny_nz_HI_photoionization_rate_array(a,-b,-i_level_py)+HI_photoionization_rate_dz
		  pyramid_source_px_ny_nz_HeI_photoionization_rate_array(i_level_py,-a,-b) = &
		             pyramid_source_px_ny_nz_HeI_photoionization_rate_array(i_level_py,-a,-b)+HeI_photoionization_rate_dx
		  pyramid_source_px_ny_nz_HeI_photoionization_rate_array(b,-i_level_py,-a) = &
		             pyramid_source_px_ny_nz_HeI_photoionization_rate_array(b,-i_level_py,-a)+HeI_photoionization_rate_dy
		  pyramid_source_px_ny_nz_HeI_photoionization_rate_array(a,-b,-i_level_py) = &
		             pyramid_source_px_ny_nz_HeI_photoionization_rate_array(a,-b,-i_level_py)+HeI_photoionization_rate_dz
		  pyramid_source_px_ny_nz_HeII_photoionization_rate_array(i_level_py,-a,-b) = &
		             pyramid_source_px_ny_nz_HeII_photoionization_rate_array(i_level_py,-a,-b)+HeII_photoionization_rate_dx
		  pyramid_source_px_ny_nz_HeII_photoionization_rate_array(b,-i_level_py,-a) = &
		             pyramid_source_px_ny_nz_HeII_photoionization_rate_array(b,-i_level_py,-a)+HeII_photoionization_rate_dy
		  pyramid_source_px_ny_nz_HeII_photoionization_rate_array(a,-b,-i_level_py) = &
		             pyramid_source_px_ny_nz_HeII_photoionization_rate_array(a,-b,-i_level_py)+HeII_photoionization_rate_dz
		  pyramid_source_px_ny_nz_photoheating_rate_array(i_level_py,-a,-b) = &
		             pyramid_source_px_ny_nz_photoheating_rate_array(i_level_py,-a,-b)+photoheating_rate_dx
		  pyramid_source_px_ny_nz_photoheating_rate_array(b,-i_level_py,-a) = &
		             pyramid_source_px_ny_nz_photoheating_rate_array(b,-i_level_py,-a)+photoheating_rate_dy
		  pyramid_source_px_ny_nz_photoheating_rate_array(a,-b,-i_level_py) = &
		             pyramid_source_px_ny_nz_photoheating_rate_array(a,-b,-i_level_py)+photoheating_rate_dz
															 					 
        enddo
      enddo		
    enddo
	  
  end subroutine px_ny_nz_photo_rate_calculation

  subroutine nx_py_pz_photo_rate_calculation()
	
    integer :: i_level_py,a,b
    integer :: i_primary,i_secondary
    real(kind=dp) :: HI_photoionization_rate_dx,HI_photoionization_rate_dy,HI_photoionization_rate_dz
    real(kind=dp) :: HeI_photoionization_rate_dx,HeI_photoionization_rate_dy,HeI_photoionization_rate_dz
    real(kind=dp) :: HeII_photoionization_rate_dx,HeII_photoionization_rate_dy,HeII_photoionization_rate_dz	
    real(kind=dp) :: photoheating_rate_dx,photoheating_rate_dy,photoheating_rate_dz
		
    do i_level_py=1,level_py
	  do a=1,i_level_py
        do b=1,i_level_py
			  
          HI_photoionization_rate_dx = 0
          HI_photoionization_rate_dy = 0
          HI_photoionization_rate_dz = 0
          HeI_photoionization_rate_dx = 0
          HeI_photoionization_rate_dy = 0
          HeI_photoionization_rate_dz = 0
          HeII_photoionization_rate_dx = 0
          HeII_photoionization_rate_dy = 0
          HeII_photoionization_rate_dz = 0
          photoheating_rate_dx = 0
          photoheating_rate_dy = 0
          photoheating_rate_dz = 0
									
          do i_primary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1),&
  	                     ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1)	
            do i_secondary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2),&
                             ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2)

if (cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) .gt. 0) then

if (ionization_weighted_by_atom.eqv..true.) then

              HI_photoionization_rate_dx = HI_photoionization_rate_dx + &
                                           nx_py_pz_dx(-i_level_py)%HI_photoionization_rate(i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_py_pz_HI_density_array(-i_level_py,a,b)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_py_pz_dx(-i_level_py)%HI_density(i_primary,i_secondary))
              HI_photoionization_rate_dy = HI_photoionization_rate_dy + &
                                           nx_py_pz_dy(i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_py_pz_HI_density_array(-b,i_level_py,a)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_py_pz_dy(i_level_py)%HI_density(i_primary,-i_secondary))
              HI_photoionization_rate_dz = HI_photoionization_rate_dz + &
                                           nx_py_pz_dz(i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_py_pz_HI_density_array(-a,b,i_level_py)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_py_pz_dz(i_level_py)%HI_density(-i_primary,i_secondary))
              HeI_photoionization_rate_dx = HeI_photoionization_rate_dx + &
                                           nx_py_pz_dx(-i_level_py)%HeI_photoionization_rate(i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_py_pz_HeI_density_array(-i_level_py,a,b)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_py_pz_dx(-i_level_py)%HeI_density(i_primary,i_secondary))
              HeI_photoionization_rate_dy = HeI_photoionization_rate_dy + &
                                           nx_py_pz_dy(i_level_py)%HeI_photoionization_rate(i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_py_pz_HeI_density_array(-b,i_level_py,a)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_py_pz_dy(i_level_py)%HeI_density(i_primary,-i_secondary))
              HeI_photoionization_rate_dz = HeI_photoionization_rate_dz + &
                                           nx_py_pz_dz(i_level_py)%HeI_photoionization_rate(-i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_py_pz_HeI_density_array(-a,b,i_level_py)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_py_pz_dz(i_level_py)%HeI_density(-i_primary,i_secondary))
              HeII_photoionization_rate_dx = HeII_photoionization_rate_dx + &
                                           nx_py_pz_dx(-i_level_py)%HeII_photoionization_rate(i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_py_pz_HeII_density_array(-i_level_py,a,b)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_py_pz_dx(-i_level_py)%HeII_density(i_primary,i_secondary))
              HeII_photoionization_rate_dy = HeII_photoionization_rate_dy + &
                                           nx_py_pz_dy(i_level_py)%HeII_photoionization_rate(i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_py_pz_HeII_density_array(-b,i_level_py,a)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_py_pz_dy(i_level_py)%HeII_density(i_primary,-i_secondary))
              HeII_photoionization_rate_dz = HeII_photoionization_rate_dz + &
                                           nx_py_pz_dz(i_level_py)%HeII_photoionization_rate(-i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_py_pz_HeII_density_array(-a,b,i_level_py)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_py_pz_dz(i_level_py)%HeII_density(-i_primary,i_secondary))


else

              HI_photoionization_rate_dx = HI_photoionization_rate_dx + &
			                  nx_py_pz_dx(-i_level_py)%HI_photoionization_rate(i_primary,i_secondary) * &
                              cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                              pyramid_grid(i_level_py)%volume
              HI_photoionization_rate_dy = HI_photoionization_rate_dy + &
			                  nx_py_pz_dy(i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume							  
              HI_photoionization_rate_dz = HI_photoionization_rate_dz + &
			                  nx_py_pz_dz(i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume						  
              HeI_photoionization_rate_dx = HeI_photoionization_rate_dx + &
			                  nx_py_pz_dx(-i_level_py)%HeI_photoionization_rate(i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
              HeI_photoionization_rate_dy = HeI_photoionization_rate_dy + &
			                  nx_py_pz_dy(i_level_py)%HeI_photoionization_rate(i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume							  
              HeI_photoionization_rate_dz = HeI_photoionization_rate_dz + &
			                  nx_py_pz_dz(i_level_py)%HeI_photoionization_rate(-i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume	
              HeII_photoionization_rate_dx = HeII_photoionization_rate_dx + &
			                  nx_py_pz_dx(-i_level_py)%HeII_photoionization_rate(i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
              HeII_photoionization_rate_dy = HeII_photoionization_rate_dy + &
			                  nx_py_pz_dy(i_level_py)%HeII_photoionization_rate(i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume							  
              HeII_photoionization_rate_dz = HeII_photoionization_rate_dz + &
			                  nx_py_pz_dz(i_level_py)%HeII_photoionization_rate(-i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
					
endif

if (heating_weighted_by_atom.eqv..true.) then

              photoheating_rate_dx = photoheating_rate_dx + &
                                     nx_py_pz_dx(-i_level_py)%photoheating_rate(i_primary,i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                     pyramid_source_nx_py_pz_HI_density_array(-i_level_py,a,b)/ &
                                     (pyramid_grid(i_level_py)%volume * &
                                     nx_py_pz_dx(-i_level_py)%HI_density(i_primary,i_secondary))
              photoheating_rate_dy = photoheating_rate_dy + &
                                     nx_py_pz_dy(i_level_py)%photoheating_rate(i_primary,-i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                     pyramid_source_nx_py_pz_HI_density_array(-b,i_level_py,a)/ &
                                     (pyramid_grid(i_level_py)%volume * &
                                     nx_py_pz_dy(i_level_py)%HI_density(i_primary,-i_secondary))
              photoheating_rate_dz = photoheating_rate_dz + &
                                     nx_py_pz_dz(i_level_py)%photoheating_rate(-i_primary,i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                     pyramid_source_nx_py_pz_HI_density_array(-a,b,i_level_py)/ &
                                     (pyramid_grid(i_level_py)%volume * &
                                     nx_py_pz_dz(i_level_py)%HI_density(-i_primary,i_secondary))

else
              photoheating_rate_dx = photoheating_rate_dx + &
                                     nx_py_pz_dx(-i_level_py)%photoheating_rate(i_primary,i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                     pyramid_grid(i_level_py)%volume
              photoheating_rate_dy = photoheating_rate_dy + &
                                     nx_py_pz_dy(i_level_py)%photoheating_rate(i_primary,-i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                     pyramid_grid(i_level_py)%volume
              photoheating_rate_dz = photoheating_rate_dz + &
                                     nx_py_pz_dz(i_level_py)%photoheating_rate(-i_primary,i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                     pyramid_grid(i_level_py)%volume	



endif	

endif											  										  				  							  
            enddo
          enddo
		  
		  pyramid_source_nx_py_pz_HI_photoionization_rate_array(-i_level_py,a,b) = &
		             pyramid_source_nx_py_pz_HI_photoionization_rate_array(-i_level_py,a,b)+HI_photoionization_rate_dx
          pyramid_source_nx_py_pz_HI_photoionization_rate_array(-b,i_level_py,a) = &
		             pyramid_source_nx_py_pz_HI_photoionization_rate_array(-b,i_level_py,a)+HI_photoionization_rate_dy
          pyramid_source_nx_py_pz_HI_photoionization_rate_array(-a,b,i_level_py) = &
		             pyramid_source_nx_py_pz_HI_photoionization_rate_array(-a,b,i_level_py)+HI_photoionization_rate_dz
          pyramid_source_nx_py_pz_HeI_photoionization_rate_array(-i_level_py,a,b) = &
		             pyramid_source_nx_py_pz_HeI_photoionization_rate_array(-i_level_py,a,b)+HeI_photoionization_rate_dx
          pyramid_source_nx_py_pz_HeI_photoionization_rate_array(-b,i_level_py,a) = &
		             pyramid_source_nx_py_pz_HeI_photoionization_rate_array(-b,i_level_py,a)+HeI_photoionization_rate_dy
          pyramid_source_nx_py_pz_HeI_photoionization_rate_array(-a,b,i_level_py) = &
		             pyramid_source_nx_py_pz_HeI_photoionization_rate_array(-a,b,i_level_py)+HeI_photoionization_rate_dz
          pyramid_source_nx_py_pz_HeII_photoionization_rate_array(-i_level_py,a,b) = &
		             pyramid_source_nx_py_pz_HeII_photoionization_rate_array(-i_level_py,a,b)+HeII_photoionization_rate_dx
          pyramid_source_nx_py_pz_HeII_photoionization_rate_array(-b,i_level_py,a) = &
		             pyramid_source_nx_py_pz_HeII_photoionization_rate_array(-b,i_level_py,a)+HeII_photoionization_rate_dy
          pyramid_source_nx_py_pz_HeII_photoionization_rate_array(-a,b,i_level_py) = &
		             pyramid_source_nx_py_pz_HeII_photoionization_rate_array(-a,b,i_level_py)+HeII_photoionization_rate_dz
		  pyramid_source_nx_py_pz_photoheating_rate_array(-i_level_py,a,b) = &
		             pyramid_source_nx_py_pz_photoheating_rate_array(-i_level_py,a,b)+photoheating_rate_dx
		  pyramid_source_nx_py_pz_photoheating_rate_array(-b,i_level_py,a) = &
		             pyramid_source_nx_py_pz_photoheating_rate_array(-b,i_level_py,a)+photoheating_rate_dy
		  pyramid_source_nx_py_pz_photoheating_rate_array(-a,b,i_level_py) = &
		             pyramid_source_nx_py_pz_photoheating_rate_array(-a,b,i_level_py)+photoheating_rate_dz
															 					 
        enddo
      enddo		
    enddo
	  
  end subroutine nx_py_pz_photo_rate_calculation

  subroutine nx_py_nz_photo_rate_calculation()
	
    integer :: i_level_py,a,b
    integer :: i_primary,i_secondary
    real(kind=dp) :: HI_photoionization_rate_dx,HI_photoionization_rate_dy,HI_photoionization_rate_dz
    real(kind=dp) :: HeI_photoionization_rate_dx,HeI_photoionization_rate_dy,HeI_photoionization_rate_dz
    real(kind=dp) :: HeII_photoionization_rate_dx,HeII_photoionization_rate_dy,HeII_photoionization_rate_dz	
    real(kind=dp) :: photoheating_rate_dx,photoheating_rate_dy,photoheating_rate_dz
		
    do i_level_py=1,level_py
	  do a=1,i_level_py
        do b=1,i_level_py
			  
          HI_photoionization_rate_dx = 0
          HI_photoionization_rate_dy = 0
          HI_photoionization_rate_dz = 0
          HeI_photoionization_rate_dx = 0
          HeI_photoionization_rate_dy = 0
          HeI_photoionization_rate_dz = 0
          HeII_photoionization_rate_dx = 0
          HeII_photoionization_rate_dy = 0
          HeII_photoionization_rate_dz = 0
          photoheating_rate_dx = 0
          photoheating_rate_dy = 0
          photoheating_rate_dz = 0
									
          do i_primary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1),&
  	                     ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1)	
            do i_secondary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2),&
                             ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2)
 
if (cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) .gt. 0) then

if (ionization_weighted_by_atom.eqv..true.) then

              HI_photoionization_rate_dx = HI_photoionization_rate_dx + &
                                           nx_py_nz_dx(-i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_py_nz_HI_density_array(-i_level_py,a,-b)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_py_nz_dx(-i_level_py)%HI_density(i_primary,-i_secondary))
              HI_photoionization_rate_dy = HI_photoionization_rate_dy + &
                                           nx_py_nz_dy(i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_py_nz_HI_density_array(-b,i_level_py,-a)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_py_nz_dy(i_level_py)%HI_density(-i_primary,-i_secondary))
              HI_photoionization_rate_dz = HI_photoionization_rate_dz + &
                                           nx_py_nz_dz(-i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_py_nz_HI_density_array(-a,b,-i_level_py)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_py_nz_dz(-i_level_py)%HI_density(-i_primary,i_secondary))
              HeI_photoionization_rate_dx = HeI_photoionization_rate_dx + &
                                           nx_py_nz_dx(-i_level_py)%HeI_photoionization_rate(i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_py_nz_HeI_density_array(-i_level_py,a,-b)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_py_nz_dx(-i_level_py)%HeI_density(i_primary,-i_secondary))
              HeI_photoionization_rate_dy = HeI_photoionization_rate_dy + &
                                           nx_py_nz_dy(i_level_py)%HeI_photoionization_rate(-i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_py_nz_HeI_density_array(-b,i_level_py,-a)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_py_nz_dy(i_level_py)%HeI_density(-i_primary,-i_secondary))
              HeI_photoionization_rate_dz = HeI_photoionization_rate_dz + &
                                           nx_py_nz_dz(-i_level_py)%HeI_photoionization_rate(-i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_py_nz_HeI_density_array(-a,b,-i_level_py)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_py_nz_dz(-i_level_py)%HeI_density(-i_primary,i_secondary))
              HeII_photoionization_rate_dx = HeII_photoionization_rate_dx + &
                                           nx_py_nz_dx(-i_level_py)%HeII_photoionization_rate(i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_py_nz_HeII_density_array(-i_level_py,a,-b)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_py_nz_dx(-i_level_py)%HeII_density(i_primary,-i_secondary))
              HeII_photoionization_rate_dy = HeII_photoionization_rate_dy + &
                                           nx_py_nz_dy(i_level_py)%HeII_photoionization_rate(-i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_py_nz_HeII_density_array(-b,i_level_py,-a)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_py_nz_dy(i_level_py)%HeII_density(-i_primary,-i_secondary))
              HeII_photoionization_rate_dz = HeII_photoionization_rate_dz + &
                                           nx_py_nz_dz(-i_level_py)%HeII_photoionization_rate(-i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_py_nz_HeII_density_array(-a,b,-i_level_py)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_py_nz_dz(-i_level_py)%HeII_density(-i_primary,i_secondary))


else

              HI_photoionization_rate_dx = HI_photoionization_rate_dx + &
			                  nx_py_nz_dx(-i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) * &
                              cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                              pyramid_grid(i_level_py)%volume
              HI_photoionization_rate_dy = HI_photoionization_rate_dy + &
			                  nx_py_nz_dy(i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume							  
              HI_photoionization_rate_dz = HI_photoionization_rate_dz + &
			                  nx_py_nz_dz(-i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume	
              HeI_photoionization_rate_dx = HeI_photoionization_rate_dx + &
			                  nx_py_nz_dx(-i_level_py)%HeI_photoionization_rate(i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
              HeI_photoionization_rate_dy = HeI_photoionization_rate_dy + &
			                  nx_py_nz_dy(i_level_py)%HeI_photoionization_rate(-i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume							  
              HeI_photoionization_rate_dz = HeI_photoionization_rate_dz + &
			                  nx_py_nz_dz(-i_level_py)%HeI_photoionization_rate(-i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
              HeII_photoionization_rate_dx = HeII_photoionization_rate_dx + &
			                  nx_py_nz_dx(-i_level_py)%HeII_photoionization_rate(i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
              HeII_photoionization_rate_dy = HeII_photoionization_rate_dy + &
			                  nx_py_nz_dy(i_level_py)%HeII_photoionization_rate(-i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume							  
              HeII_photoionization_rate_dz = HeII_photoionization_rate_dz + &
			                  nx_py_nz_dz(-i_level_py)%HeII_photoionization_rate(-i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
					
endif

if (heating_weighted_by_atom.eqv..true.) then

              photoheating_rate_dx = photoheating_rate_dx + &
                                     nx_py_nz_dx(-i_level_py)%photoheating_rate(i_primary,-i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                     pyramid_source_nx_py_nz_HI_density_array(-i_level_py,a,-b)/ &
                                     (pyramid_grid(i_level_py)%volume * &
                                     nx_py_nz_dx(-i_level_py)%HI_density(i_primary,-i_secondary))
              photoheating_rate_dy = photoheating_rate_dy + &
                                     nx_py_nz_dy(i_level_py)%photoheating_rate(-i_primary,-i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                     pyramid_source_nx_py_nz_HI_density_array(-b,i_level_py,-a)/ &
                                     (pyramid_grid(i_level_py)%volume * &
                                     nx_py_nz_dy(i_level_py)%HI_density(-i_primary,-i_secondary))
              photoheating_rate_dz = photoheating_rate_dz + &
                                     nx_py_nz_dz(-i_level_py)%photoheating_rate(-i_primary,i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                     pyramid_source_nx_py_nz_HI_density_array(-a,b,-i_level_py)/ &
                                     (pyramid_grid(i_level_py)%volume * &
                                     nx_py_nz_dz(-i_level_py)%HI_density(-i_primary,i_secondary))
	else
              photoheating_rate_dx = photoheating_rate_dx + &
                                     nx_py_nz_dx(-i_level_py)%photoheating_rate(i_primary,-i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                     pyramid_grid(i_level_py)%volume
              photoheating_rate_dy = photoheating_rate_dy + &
                                     nx_py_nz_dy(i_level_py)%photoheating_rate(-i_primary,-i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                     pyramid_grid(i_level_py)%volume
              photoheating_rate_dz = photoheating_rate_dz + &
                                     nx_py_nz_dz(-i_level_py)%photoheating_rate(-i_primary,i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                     pyramid_grid(i_level_py)%volume



endif

endif

            enddo
          enddo
		  
		  pyramid_source_nx_py_nz_HI_photoionization_rate_array(-i_level_py,a,-b) = &
		             pyramid_source_nx_py_nz_HI_photoionization_rate_array(-i_level_py,a,-b)+HI_photoionization_rate_dx
          pyramid_source_nx_py_nz_HI_photoionization_rate_array(-b,i_level_py,-a) = &
		             pyramid_source_nx_py_nz_HI_photoionization_rate_array(-b,i_level_py,-a)+HI_photoionization_rate_dy
          pyramid_source_nx_py_nz_HI_photoionization_rate_array(-a,b,-i_level_py) = &
		             pyramid_source_nx_py_nz_HI_photoionization_rate_array(-a,b,-i_level_py)+HI_photoionization_rate_dz
          pyramid_source_nx_py_nz_HeI_photoionization_rate_array(-i_level_py,a,-b) = &
		             pyramid_source_nx_py_nz_HeI_photoionization_rate_array(-i_level_py,a,-b)+HeI_photoionization_rate_dx
          pyramid_source_nx_py_nz_HeI_photoionization_rate_array(-b,i_level_py,-a) = &
		             pyramid_source_nx_py_nz_HeI_photoionization_rate_array(-b,i_level_py,-a)+HeI_photoionization_rate_dy
          pyramid_source_nx_py_nz_HeI_photoionization_rate_array(-a,b,-i_level_py) = &
		             pyramid_source_nx_py_nz_HeI_photoionization_rate_array(-a,b,-i_level_py)+HeI_photoionization_rate_dz
          pyramid_source_nx_py_nz_HeII_photoionization_rate_array(-i_level_py,a,-b) = &
		             pyramid_source_nx_py_nz_HeII_photoionization_rate_array(-i_level_py,a,-b)+HeII_photoionization_rate_dx
          pyramid_source_nx_py_nz_HeII_photoionization_rate_array(-b,i_level_py,-a) = &
		             pyramid_source_nx_py_nz_HeII_photoionization_rate_array(-b,i_level_py,-a)+HeII_photoionization_rate_dy
          pyramid_source_nx_py_nz_HeII_photoionization_rate_array(-a,b,-i_level_py) = &
		             pyramid_source_nx_py_nz_HeII_photoionization_rate_array(-a,b,-i_level_py)+HeII_photoionization_rate_dz
		  pyramid_source_nx_py_nz_photoheating_rate_array(-i_level_py,a,-b) = &
		             pyramid_source_nx_py_nz_photoheating_rate_array(-i_level_py,a,-b)+photoheating_rate_dx
		  pyramid_source_nx_py_nz_photoheating_rate_array(-b,i_level_py,-a) = &
		             pyramid_source_nx_py_nz_photoheating_rate_array(-b,i_level_py,-a)+photoheating_rate_dy
		  pyramid_source_nx_py_nz_photoheating_rate_array(-a,b,-i_level_py) = &
		             pyramid_source_nx_py_nz_photoheating_rate_array(-a,b,-i_level_py)+photoheating_rate_dz
															 					 
        enddo
      enddo		
    enddo

  end subroutine nx_py_nz_photo_rate_calculation

  subroutine nx_ny_pz_photo_rate_calculation()
	
    integer :: i_level_py,a,b
    integer :: i_primary,i_secondary
    real(kind=dp) :: HI_photoionization_rate_dx,HI_photoionization_rate_dy,HI_photoionization_rate_dz
    real(kind=dp) :: HeI_photoionization_rate_dx,HeI_photoionization_rate_dy,HeI_photoionization_rate_dz
    real(kind=dp) :: HeII_photoionization_rate_dx,HeII_photoionization_rate_dy,HeII_photoionization_rate_dz	
    real(kind=dp) :: photoheating_rate_dx,photoheating_rate_dy,photoheating_rate_dz
		
    do i_level_py=1,level_py
	  do a=1,i_level_py
        do b=1,i_level_py
			  
          HI_photoionization_rate_dx = 0
          HI_photoionization_rate_dy = 0
          HI_photoionization_rate_dz = 0
          HeI_photoionization_rate_dx = 0
          HeI_photoionization_rate_dy = 0
          HeI_photoionization_rate_dz = 0
          HeII_photoionization_rate_dx = 0
          HeII_photoionization_rate_dy = 0
          HeII_photoionization_rate_dz = 0
          photoheating_rate_dx = 0
          photoheating_rate_dy = 0
          photoheating_rate_dz = 0
									
          do i_primary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1),&
  	                     ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1)	
            do i_secondary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2),&
                             ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2)
 
if (cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) .gt. 0) then

if (ionization_weighted_by_atom.eqv..true.) then

              HI_photoionization_rate_dx = HI_photoionization_rate_dx + &
                                           nx_ny_pz_dx(-i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_ny_pz_HI_density_array(-i_level_py,-a,b)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_ny_pz_dx(-i_level_py)%HI_density(-i_primary,i_secondary))
              HI_photoionization_rate_dy = HI_photoionization_rate_dy + &
                                           nx_ny_pz_dy(-i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_ny_pz_HI_density_array(-b,-i_level_py,a)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_ny_pz_dy(-i_level_py)%HI_density(i_primary,-i_secondary))
              HI_photoionization_rate_dz = HI_photoionization_rate_dz + &
                                           nx_ny_pz_dz(i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_ny_pz_HI_density_array(-a,-b,i_level_py)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_ny_pz_dz(i_level_py)%HI_density(-i_primary,-i_secondary))
              HeI_photoionization_rate_dx = HeI_photoionization_rate_dx + &
                                           nx_ny_pz_dx(-i_level_py)%HeI_photoionization_rate(-i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_ny_pz_HeI_density_array(-i_level_py,-a,b)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_ny_pz_dx(-i_level_py)%HeI_density(-i_primary,i_secondary))
              HeI_photoionization_rate_dy = HeI_photoionization_rate_dy + &
                                           nx_ny_pz_dy(-i_level_py)%HeI_photoionization_rate(i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_ny_pz_HeI_density_array(-b,-i_level_py,a)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_ny_pz_dy(-i_level_py)%HeI_density(i_primary,-i_secondary))
              HeI_photoionization_rate_dz = HeI_photoionization_rate_dz + &
                                           nx_ny_pz_dz(i_level_py)%HeI_photoionization_rate(-i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_ny_pz_HeI_density_array(-a,-b,i_level_py)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_ny_pz_dz(i_level_py)%HeI_density(-i_primary,-i_secondary))
              HeII_photoionization_rate_dx = HeII_photoionization_rate_dx + &
                                           nx_ny_pz_dx(-i_level_py)%HeII_photoionization_rate(-i_primary,i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_ny_pz_HeII_density_array(-i_level_py,-a,b)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_ny_pz_dx(-i_level_py)%HeII_density(-i_primary,i_secondary))
              HeII_photoionization_rate_dy = HeII_photoionization_rate_dy + &
                                           nx_ny_pz_dy(-i_level_py)%HeII_photoionization_rate(i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_ny_pz_HeII_density_array(-b,-i_level_py,a)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_ny_pz_dy(-i_level_py)%HeII_density(i_primary,-i_secondary))
              HeII_photoionization_rate_dz = HeII_photoionization_rate_dz + &
                                           nx_ny_pz_dz(i_level_py)%HeII_photoionization_rate(-i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_ny_pz_HeII_density_array(-a,-b,i_level_py)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_ny_pz_dz(i_level_py)%HeII_density(-i_primary,-i_secondary))


else

              HI_photoionization_rate_dx = HI_photoionization_rate_dx + &
			                  nx_ny_pz_dx(-i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) * &
                              cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                              pyramid_grid(i_level_py)%volume
              HI_photoionization_rate_dy = HI_photoionization_rate_dy + &
			                  nx_ny_pz_dy(-i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume							  
              HI_photoionization_rate_dz = HI_photoionization_rate_dz + &
			                  nx_ny_pz_dz(i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume						  
              HeI_photoionization_rate_dx = HeI_photoionization_rate_dx + &
			                  nx_ny_pz_dx(-i_level_py)%HeI_photoionization_rate(-i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
              HeI_photoionization_rate_dy = HeI_photoionization_rate_dy + &
			                  nx_ny_pz_dy(-i_level_py)%HeI_photoionization_rate(i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume							  
              HeI_photoionization_rate_dz = HeI_photoionization_rate_dz + &
			                  nx_ny_pz_dz(i_level_py)%HeI_photoionization_rate(-i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
              HeII_photoionization_rate_dx = HeII_photoionization_rate_dx + &
			                  nx_ny_pz_dx(-i_level_py)%HeII_photoionization_rate(-i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
              HeII_photoionization_rate_dy = HeII_photoionization_rate_dy + &
			                  nx_ny_pz_dy(-i_level_py)%HeII_photoionization_rate(i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume							  
              HeII_photoionization_rate_dz = HeII_photoionization_rate_dz + &
			                  nx_ny_pz_dz(i_level_py)%HeII_photoionization_rate(-i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
					
endif

if (heating_weighted_by_atom.eqv..true.) then

              photoheating_rate_dx = photoheating_rate_dx + &
                                     nx_ny_pz_dx(-i_level_py)%photoheating_rate(-i_primary,i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                     pyramid_source_nx_ny_pz_HI_density_array(-i_level_py,-a,b)/ &
                                     (pyramid_grid(i_level_py)%volume * &
                                     nx_ny_pz_dx(-i_level_py)%HI_density(-i_primary,i_secondary))
              photoheating_rate_dy = photoheating_rate_dy + &
                                     nx_ny_pz_dy(-i_level_py)%photoheating_rate(i_primary,-i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                     pyramid_source_nx_ny_pz_HI_density_array(-b,-i_level_py,a)/ &
                                     (pyramid_grid(i_level_py)%volume * &
                                     nx_ny_pz_dy(-i_level_py)%HI_density(i_primary,-i_secondary))
              photoheating_rate_dz = photoheating_rate_dz + &
                                     nx_ny_pz_dz(i_level_py)%photoheating_rate(-i_primary,-i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                     pyramid_source_nx_ny_pz_HI_density_array(-a,-b,i_level_py)/ &
                                     (pyramid_grid(i_level_py)%volume * &
                                     nx_ny_pz_dz(i_level_py)%HI_density(-i_primary,-i_secondary))
else
              photoheating_rate_dx = photoheating_rate_dx + &
                                     nx_ny_pz_dx(-i_level_py)%photoheating_rate(-i_primary,i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                     pyramid_grid(i_level_py)%volume
              photoheating_rate_dy = photoheating_rate_dy + &
                                     nx_ny_pz_dy(-i_level_py)%photoheating_rate(i_primary,-i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                     pyramid_grid(i_level_py)%volume
              photoheating_rate_dz = photoheating_rate_dz + &
                                     nx_ny_pz_dz(i_level_py)%photoheating_rate(-i_primary,-i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                     pyramid_grid(i_level_py)%volume



endif
											  	endif									  				  							  
            enddo
          enddo
		  
		  pyramid_source_nx_ny_pz_HI_photoionization_rate_array(-i_level_py,-a,b) = &
		             pyramid_source_nx_ny_pz_HI_photoionization_rate_array(-i_level_py,-a,b)+HI_photoionization_rate_dx
          pyramid_source_nx_ny_pz_HI_photoionization_rate_array(-b,-i_level_py,a) = &
		             pyramid_source_nx_ny_pz_HI_photoionization_rate_array(-b,-i_level_py,a)+HI_photoionization_rate_dy
          pyramid_source_nx_ny_pz_HI_photoionization_rate_array(-a,-b,i_level_py) = &
		             pyramid_source_nx_ny_pz_HI_photoionization_rate_array(-a,-b,i_level_py)+HI_photoionization_rate_dz
          pyramid_source_nx_ny_pz_HeI_photoionization_rate_array(-i_level_py,-a,b) = &
		             pyramid_source_nx_ny_pz_HeI_photoionization_rate_array(-i_level_py,-a,b)+HeI_photoionization_rate_dx
          pyramid_source_nx_ny_pz_HeI_photoionization_rate_array(-b,-i_level_py,a) = &
		             pyramid_source_nx_ny_pz_HeI_photoionization_rate_array(-b,-i_level_py,a)+HeI_photoionization_rate_dy
          pyramid_source_nx_ny_pz_HeI_photoionization_rate_array(-a,-b,i_level_py) = &
		             pyramid_source_nx_ny_pz_HeI_photoionization_rate_array(-a,-b,i_level_py)+HeI_photoionization_rate_dz
          pyramid_source_nx_ny_pz_HeII_photoionization_rate_array(-i_level_py,-a,b) = &
		             pyramid_source_nx_ny_pz_HeII_photoionization_rate_array(-i_level_py,-a,b)+HeII_photoionization_rate_dx
          pyramid_source_nx_ny_pz_HeII_photoionization_rate_array(-b,-i_level_py,a) = &
		             pyramid_source_nx_ny_pz_HeII_photoionization_rate_array(-b,-i_level_py,a)+HeII_photoionization_rate_dy
          pyramid_source_nx_ny_pz_HeII_photoionization_rate_array(-a,-b,i_level_py) = &
		             pyramid_source_nx_ny_pz_HeII_photoionization_rate_array(-a,-b,i_level_py)+HeII_photoionization_rate_dz
		  pyramid_source_nx_ny_pz_photoheating_rate_array(-i_level_py,-a,b) = &
		             pyramid_source_nx_ny_pz_photoheating_rate_array(-i_level_py,-a,b)+photoheating_rate_dx
		  pyramid_source_nx_ny_pz_photoheating_rate_array(-b,-i_level_py,a) = &
		             pyramid_source_nx_ny_pz_photoheating_rate_array(-b,-i_level_py,a)+photoheating_rate_dy
		  pyramid_source_nx_ny_pz_photoheating_rate_array(-a,-b,i_level_py) = &
		             pyramid_source_nx_ny_pz_photoheating_rate_array(-a,-b,i_level_py)+photoheating_rate_dz
															 					 
        enddo
      enddo		
    enddo
	  
  end subroutine nx_ny_pz_photo_rate_calculation

  subroutine nx_ny_nz_photo_rate_calculation()
	
    integer :: i_level_py,a,b
    integer :: i_primary,i_secondary
    real(kind=dp) :: HI_photoionization_rate_dx,HI_photoionization_rate_dy,HI_photoionization_rate_dz
    real(kind=dp) :: HeI_photoionization_rate_dx,HeI_photoionization_rate_dy,HeI_photoionization_rate_dz
    real(kind=dp) :: HeII_photoionization_rate_dx,HeII_photoionization_rate_dy,HeII_photoionization_rate_dz	
    real(kind=dp) :: photoheating_rate_dx,photoheating_rate_dy,photoheating_rate_dz
		
    do i_level_py=1,level_py
	  do a=1,i_level_py
        do b=1,i_level_py
			  
          HI_photoionization_rate_dx = 0
          HI_photoionization_rate_dy = 0
          HI_photoionization_rate_dz = 0
          HeI_photoionization_rate_dx = 0
          HeI_photoionization_rate_dy = 0
          HeI_photoionization_rate_dz = 0
          HeII_photoionization_rate_dx = 0
          HeII_photoionization_rate_dy = 0
          HeII_photoionization_rate_dz = 0
          photoheating_rate_dx = 0
          photoheating_rate_dy = 0
          photoheating_rate_dz = 0
									
          do i_primary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1),&
  	                     ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1)	
            do i_secondary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2),&
                             ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2)
 
if (cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) .gt. 0) then

if (ionization_weighted_by_atom.eqv..true.) then

              HI_photoionization_rate_dx = HI_photoionization_rate_dx + &
                                           nx_ny_nz_dx(-i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_ny_nz_HI_density_array(-i_level_py,-a,-b)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_ny_nz_dx(-i_level_py)%HI_density(-i_primary,-i_secondary))
              HI_photoionization_rate_dy = HI_photoionization_rate_dy + &
                                           nx_ny_nz_dy(-i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_ny_nz_HI_density_array(-b,-i_level_py,-a)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_ny_nz_dy(-i_level_py)%HI_density(-i_primary,-i_secondary))
              HI_photoionization_rate_dz = HI_photoionization_rate_dz + &
                                           nx_ny_nz_dz(-i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_ny_nz_HI_density_array(-a,-b,-i_level_py)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_ny_nz_dz(-i_level_py)%HI_density(-i_primary,-i_secondary))
              HeI_photoionization_rate_dx = HeI_photoionization_rate_dx + &
                                           nx_ny_nz_dx(-i_level_py)%HeI_photoionization_rate(-i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_ny_nz_HeI_density_array(-i_level_py,-a,-b)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_ny_nz_dx(-i_level_py)%HeI_density(-i_primary,-i_secondary))
              HeI_photoionization_rate_dy = HeI_photoionization_rate_dy + &
                                           nx_ny_nz_dy(-i_level_py)%HeI_photoionization_rate(-i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_ny_nz_HeI_density_array(-b,-i_level_py,-a)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_ny_nz_dy(-i_level_py)%HeI_density(-i_primary,-i_secondary))
              HeI_photoionization_rate_dz = HeI_photoionization_rate_dz + &
                                           nx_ny_nz_dz(-i_level_py)%HeI_photoionization_rate(-i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_ny_nz_HeI_density_array(-a,-b,-i_level_py)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_ny_nz_dz(-i_level_py)%HeI_density(-i_primary,-i_secondary))
              HeII_photoionization_rate_dx = HeII_photoionization_rate_dx + &
                                           nx_ny_nz_dx(-i_level_py)%HeII_photoionization_rate(-i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_ny_nz_HeII_density_array(-i_level_py,-a,-b)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_ny_nz_dx(-i_level_py)%HeII_density(-i_primary,-i_secondary))
              HeII_photoionization_rate_dy = HeII_photoionization_rate_dy + &
                                           nx_ny_nz_dy(-i_level_py)%HeII_photoionization_rate(-i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_ny_nz_HeII_density_array(-b,-i_level_py,-a)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_ny_nz_dy(-i_level_py)%HeII_density(-i_primary,-i_secondary))
              HeII_photoionization_rate_dz = HeII_photoionization_rate_dz + &
                                           nx_ny_nz_dz(-i_level_py)%HeII_photoionization_rate(-i_primary,-i_secondary) * &
                                           cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                           pyramid_source_nx_ny_nz_HeII_density_array(-a,-b,-i_level_py)/ &
                                           (pyramid_grid(i_level_py)%volume * &
                                           nx_ny_nz_dz(-i_level_py)%HeII_density(-i_primary,-i_secondary))

else

              HI_photoionization_rate_dx = HI_photoionization_rate_dx + &
			                  nx_ny_nz_dx(-i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) * &
                              cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                              pyramid_grid(i_level_py)%volume
              HI_photoionization_rate_dy = HI_photoionization_rate_dy + &
			                  nx_ny_nz_dy(-i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume							  
              HI_photoionization_rate_dz = HI_photoionization_rate_dz + &
			                  nx_ny_nz_dz(-i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume						  
              HeI_photoionization_rate_dx = HeI_photoionization_rate_dx + &
			                  nx_ny_nz_dx(-i_level_py)%HeI_photoionization_rate(-i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
              HeI_photoionization_rate_dy = HeI_photoionization_rate_dy + &
			                  nx_ny_nz_dy(-i_level_py)%HeI_photoionization_rate(-i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume							  
              HeI_photoionization_rate_dz = HeI_photoionization_rate_dz + &
			                  nx_ny_nz_dz(-i_level_py)%HeI_photoionization_rate(-i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
              HeII_photoionization_rate_dx = HeII_photoionization_rate_dx + &
			                  nx_ny_nz_dx(-i_level_py)%HeII_photoionization_rate(-i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
              HeII_photoionization_rate_dy = HeII_photoionization_rate_dy + &
			                  nx_ny_nz_dy(-i_level_py)%HeII_photoionization_rate(-i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume							  
              HeII_photoionization_rate_dz = HeII_photoionization_rate_dz + &
			                  nx_ny_nz_dz(-i_level_py)%HeII_photoionization_rate(-i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume
					
endif

if (heating_weighted_by_atom.eqv..true.) then

              photoheating_rate_dx = photoheating_rate_dx + &
                                     nx_ny_nz_dx(-i_level_py)%photoheating_rate(-i_primary,-i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                     pyramid_source_nx_ny_nz_HI_density_array(-i_level_py,-a,-b)/ &
                                     (pyramid_grid(i_level_py)%volume * &
                                     nx_ny_nz_dx(-i_level_py)%HI_density(-i_primary,-i_secondary))
              photoheating_rate_dy = photoheating_rate_dy + &
                                     nx_ny_nz_dy(-i_level_py)%photoheating_rate(-i_primary,-i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                     pyramid_source_nx_ny_nz_HI_density_array(-b,-i_level_py,-a)/ &
                                     (pyramid_grid(i_level_py)%volume * &
                                     nx_ny_nz_dy(-i_level_py)%HI_density(-i_primary,-i_secondary))
              photoheating_rate_dz = photoheating_rate_dz + &
                                     nx_ny_nz_dz(-i_level_py)%photoheating_rate(-i_primary,-i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) * &
                                     pyramid_source_nx_ny_nz_HI_density_array(-a,-b,-i_level_py)/ &
                                     (pyramid_grid(i_level_py)%volume * &
                                     nx_ny_nz_dz(-i_level_py)%HI_density(-i_primary,-i_secondary))
else
              photoheating_rate_dx = photoheating_rate_dx + &
                                     nx_ny_nz_dx(-i_level_py)%photoheating_rate(-i_primary,-i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                     pyramid_grid(i_level_py)%volume
              photoheating_rate_dy = photoheating_rate_dy + &
                                     nx_ny_nz_dy(-i_level_py)%photoheating_rate(-i_primary,-i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                     pyramid_grid(i_level_py)%volume
              photoheating_rate_dz = photoheating_rate_dz + &
                                     nx_ny_nz_dz(-i_level_py)%photoheating_rate(-i_primary,-i_secondary) * &
                                     cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                                     pyramid_grid(i_level_py)%volume	



endif		
											  endif									  				  							  
            enddo
          enddo
		  
		  pyramid_source_nx_ny_nz_HI_photoionization_rate_array(-i_level_py,-a,-b) = &
		             pyramid_source_nx_ny_nz_HI_photoionization_rate_array(-i_level_py,-a,-b)+HI_photoionization_rate_dx
          pyramid_source_nx_ny_nz_HI_photoionization_rate_array(-b,-i_level_py,-a) = &
		             pyramid_source_nx_ny_nz_HI_photoionization_rate_array(-b,-i_level_py,-a)+HI_photoionization_rate_dy
          pyramid_source_nx_ny_nz_HI_photoionization_rate_array(-a,-b,-i_level_py) = &
		             pyramid_source_nx_ny_nz_HI_photoionization_rate_array(-a,-b,-i_level_py)+HI_photoionization_rate_dz
		  pyramid_source_nx_ny_nz_HeI_photoionization_rate_array(-i_level_py,-a,-b) = &
		             pyramid_source_nx_ny_nz_HeI_photoionization_rate_array(-i_level_py,-a,-b)+HeI_photoionization_rate_dx
		  pyramid_source_nx_ny_nz_HeI_photoionization_rate_array(-b,-i_level_py,-a) = &
		             pyramid_source_nx_ny_nz_HeI_photoionization_rate_array(-b,-i_level_py,-a)+HeI_photoionization_rate_dy
		  pyramid_source_nx_ny_nz_HeI_photoionization_rate_array(-a,-b,-i_level_py) = &
		             pyramid_source_nx_ny_nz_HeI_photoionization_rate_array(-a,-b,-i_level_py)+HeI_photoionization_rate_dz
		  pyramid_source_nx_ny_nz_HeII_photoionization_rate_array(-i_level_py,-a,-b) = &
		             pyramid_source_nx_ny_nz_HeII_photoionization_rate_array(-i_level_py,-a,-b)+HeII_photoionization_rate_dx
		  pyramid_source_nx_ny_nz_HeII_photoionization_rate_array(-b,-i_level_py,-a) = &
		             pyramid_source_nx_ny_nz_HeII_photoionization_rate_array(-b,-i_level_py,-a)+HeII_photoionization_rate_dy
		  pyramid_source_nx_ny_nz_HeII_photoionization_rate_array(-a,-b,-i_level_py) = &
		             pyramid_source_nx_ny_nz_HeII_photoionization_rate_array(-a,-b,-i_level_py)+HeII_photoionization_rate_dz
		  pyramid_source_nx_ny_nz_photoheating_rate_array(-i_level_py,-a,-b) = &
		             pyramid_source_nx_ny_nz_photoheating_rate_array(-i_level_py,-a,-b)+photoheating_rate_dx
		  pyramid_source_nx_ny_nz_photoheating_rate_array(-b,-i_level_py,-a) = &
		             pyramid_source_nx_ny_nz_photoheating_rate_array(-b,-i_level_py,-a)+photoheating_rate_dy
		  pyramid_source_nx_ny_nz_photoheating_rate_array(-a,-b,-i_level_py) = &
		             pyramid_source_nx_ny_nz_photoheating_rate_array(-a,-b,-i_level_py)+photoheating_rate_dz
															 					 
        enddo
      enddo		
    enddo
		  
  end subroutine nx_ny_nz_photo_rate_calculation
  
  subroutine photoionization_rate_per_atom()
	  
	pyramid_source_px_py_pz_HI_photoionization_rate_array=&
	pyramid_source_px_py_pz_HI_photoionization_rate_array/(pyramid_source_px_py_pz_HI_density_array*cellsize_py_cube)	
	pyramid_source_px_py_nz_HI_photoionization_rate_array=&
	pyramid_source_px_py_nz_HI_photoionization_rate_array/(pyramid_source_px_py_nz_HI_density_array*cellsize_py_cube)	
	pyramid_source_px_ny_pz_HI_photoionization_rate_array=&
	pyramid_source_px_ny_pz_HI_photoionization_rate_array/(pyramid_source_px_ny_pz_HI_density_array*cellsize_py_cube)	
	pyramid_source_px_ny_nz_HI_photoionization_rate_array=&
	pyramid_source_px_ny_nz_HI_photoionization_rate_array/(pyramid_source_px_ny_nz_HI_density_array*cellsize_py_cube)	
	pyramid_source_nx_py_pz_HI_photoionization_rate_array=&
	pyramid_source_nx_py_pz_HI_photoionization_rate_array/(pyramid_source_nx_py_pz_HI_density_array*cellsize_py_cube)	
	pyramid_source_nx_py_nz_HI_photoionization_rate_array=&
	pyramid_source_nx_py_nz_HI_photoionization_rate_array/(pyramid_source_nx_py_nz_HI_density_array*cellsize_py_cube)	
	pyramid_source_nx_ny_pz_HI_photoionization_rate_array=&
	pyramid_source_nx_ny_pz_HI_photoionization_rate_array/(pyramid_source_nx_ny_pz_HI_density_array*cellsize_py_cube)	
	pyramid_source_nx_ny_nz_HI_photoionization_rate_array=&
	pyramid_source_nx_ny_nz_HI_photoionization_rate_array/(pyramid_source_nx_ny_nz_HI_density_array*cellsize_py_cube)		  	  
	pyramid_source_px_py_pz_HeI_photoionization_rate_array=&
	pyramid_source_px_py_pz_HeI_photoionization_rate_array/(pyramid_source_px_py_pz_HeI_density_array*cellsize_py_cube)	
	pyramid_source_px_py_nz_HeI_photoionization_rate_array=&
	pyramid_source_px_py_nz_HeI_photoionization_rate_array/(pyramid_source_px_py_nz_HeI_density_array*cellsize_py_cube)	
	pyramid_source_px_ny_pz_HeI_photoionization_rate_array=&
	pyramid_source_px_ny_pz_HeI_photoionization_rate_array/(pyramid_source_px_ny_pz_HeI_density_array*cellsize_py_cube)	
	pyramid_source_px_ny_nz_HeI_photoionization_rate_array=&
	pyramid_source_px_ny_nz_HeI_photoionization_rate_array/(pyramid_source_px_ny_nz_HeI_density_array*cellsize_py_cube)	
	pyramid_source_nx_py_pz_HeI_photoionization_rate_array=&
	pyramid_source_nx_py_pz_HeI_photoionization_rate_array/(pyramid_source_nx_py_pz_HeI_density_array*cellsize_py_cube)	
	pyramid_source_nx_py_nz_HeI_photoionization_rate_array=&
	pyramid_source_nx_py_nz_HeI_photoionization_rate_array/(pyramid_source_nx_py_nz_HeI_density_array*cellsize_py_cube)	
	pyramid_source_nx_ny_pz_HeI_photoionization_rate_array=&
	pyramid_source_nx_ny_pz_HeI_photoionization_rate_array/(pyramid_source_nx_ny_pz_HeI_density_array*cellsize_py_cube)	
	pyramid_source_nx_ny_nz_HeI_photoionization_rate_array=&
	pyramid_source_nx_ny_nz_HeI_photoionization_rate_array/(pyramid_source_nx_ny_nz_HeI_density_array*cellsize_py_cube)
	pyramid_source_px_py_pz_HeII_photoionization_rate_array=&
	pyramid_source_px_py_pz_HeII_photoionization_rate_array/(pyramid_source_px_py_pz_HeII_density_array*cellsize_py_cube)	
	pyramid_source_px_py_nz_HeII_photoionization_rate_array=&
	pyramid_source_px_py_nz_HeII_photoionization_rate_array/(pyramid_source_px_py_nz_HeII_density_array*cellsize_py_cube)	
	pyramid_source_px_ny_pz_HeII_photoionization_rate_array=&
	pyramid_source_px_ny_pz_HeII_photoionization_rate_array/(pyramid_source_px_ny_pz_HeII_density_array*cellsize_py_cube)	
	pyramid_source_px_ny_nz_HeII_photoionization_rate_array=&
	pyramid_source_px_ny_nz_HeII_photoionization_rate_array/(pyramid_source_px_ny_nz_HeII_density_array*cellsize_py_cube)	
	pyramid_source_nx_py_pz_HeII_photoionization_rate_array=&
	pyramid_source_nx_py_pz_HeII_photoionization_rate_array/(pyramid_source_nx_py_pz_HeII_density_array*cellsize_py_cube)	
	pyramid_source_nx_py_nz_HeII_photoionization_rate_array=&
	pyramid_source_nx_py_nz_HeII_photoionization_rate_array/(pyramid_source_nx_py_nz_HeII_density_array*cellsize_py_cube)	
	pyramid_source_nx_ny_pz_HeII_photoionization_rate_array=&
	pyramid_source_nx_ny_pz_HeII_photoionization_rate_array/(pyramid_source_nx_ny_pz_HeII_density_array*cellsize_py_cube)	
	pyramid_source_nx_ny_nz_HeII_photoionization_rate_array=&
	pyramid_source_nx_ny_nz_HeII_photoionization_rate_array/(pyramid_source_nx_ny_nz_HeII_density_array*cellsize_py_cube)
			  	  
  end subroutine photoionization_rate_per_atom

  subroutine photoheating_rate_per_volume()
	  
	pyramid_source_px_py_pz_photoheating_rate_array=&
	pyramid_source_px_py_pz_photoheating_rate_array/cellsize_py_cube	
	pyramid_source_px_py_nz_photoheating_rate_array=&
	pyramid_source_px_py_nz_photoheating_rate_array/cellsize_py_cube
	pyramid_source_px_ny_pz_photoheating_rate_array=&
	pyramid_source_px_ny_pz_photoheating_rate_array/cellsize_py_cube	
	pyramid_source_px_ny_nz_photoheating_rate_array=&
	pyramid_source_px_ny_nz_photoheating_rate_array/cellsize_py_cube
	pyramid_source_nx_py_pz_photoheating_rate_array=&
	pyramid_source_nx_py_pz_photoheating_rate_array/cellsize_py_cube	
	pyramid_source_nx_py_nz_photoheating_rate_array=&
	pyramid_source_nx_py_nz_photoheating_rate_array/cellsize_py_cube	
	pyramid_source_nx_ny_pz_photoheating_rate_array=&
	pyramid_source_nx_ny_pz_photoheating_rate_array/cellsize_py_cube	
	pyramid_source_nx_ny_nz_photoheating_rate_array=&
	pyramid_source_nx_ny_nz_photoheating_rate_array/cellsize_py_cube		  	  
			  	  
  end subroutine photoheating_rate_per_volume
    
end module pyramid_source_domain_transformation
