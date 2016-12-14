module global_photo_rate
	
  use precision, only: dp
  use input, only: level_py, partition, number_density, xHI, cellsize_py_cube
  use array, only: cartesian_grid, pyramid_grid, &
                   px_py_pz_dx, px_py_pz_dy, px_py_pz_dz, &
  				   px_py_nz_dx, px_py_nz_dy, px_py_nz_dz, &
  				   px_ny_pz_dx, px_ny_pz_dy, px_ny_pz_dz, &	  
  				   px_ny_nz_dx, px_ny_nz_dy, px_ny_nz_dz, & 	   
  				   nx_py_pz_dx, nx_py_pz_dy, nx_py_pz_dz, &	  
  				   nx_py_nz_dx, nx_py_nz_dy, nx_py_nz_dz, &	
  				   nx_ny_pz_dx, nx_ny_pz_dy, nx_ny_pz_dz, &	  
  				   nx_ny_nz_dx, nx_ny_nz_dy, nx_ny_nz_dz, &
                   pyramid_source_px_py_pz_HI_photo_rate_array, &
                   pyramid_source_px_py_nz_HI_photo_rate_array, &  
                   pyramid_source_px_ny_pz_HI_photo_rate_array, &
                   pyramid_source_px_ny_nz_HI_photo_rate_array, &
                   pyramid_source_nx_py_pz_HI_photo_rate_array, &
                   pyramid_source_nx_py_nz_HI_photo_rate_array, &
                   pyramid_source_nx_ny_pz_HI_photo_rate_array, &
                   pyramid_source_nx_ny_nz_HI_photo_rate_array, &
				   pyramid_source_px_py_pz_HI_density_array, &
				   pyramid_source_px_py_nz_HI_density_array, &
				   pyramid_source_px_ny_pz_HI_density_array, &
				   pyramid_source_px_ny_nz_HI_density_array, &
				   pyramid_source_nx_py_pz_HI_density_array, &
				   pyramid_source_nx_py_nz_HI_density_array, &
				   pyramid_source_nx_ny_pz_HI_density_array, &
				   pyramid_source_nx_ny_nz_HI_density_array
contains
	
  subroutine global_HI_photo_rate_calculation()

    call source_photo_rate_initialization() 
    call px_py_pz_HI_photo_rate_calculation()
    call px_py_nz_HI_photo_rate_calculation()
    call px_ny_pz_HI_photo_rate_calculation()
    call px_ny_nz_HI_photo_rate_calculation()
    call nx_py_pz_HI_photo_rate_calculation()
    call nx_py_nz_HI_photo_rate_calculation()
    call nx_ny_pz_HI_photo_rate_calculation()
    call nx_ny_nz_HI_photo_rate_calculation()
  call photo_rate_per_atom()
  	
  end subroutine global_HI_photo_rate_calculation
  
  subroutine source_photo_rate_initialization()
	  
    pyramid_source_px_py_pz_HI_photo_rate_array = 0
    pyramid_source_px_py_nz_HI_photo_rate_array = 0
    pyramid_source_px_ny_pz_HI_photo_rate_array = 0
    pyramid_source_px_ny_nz_HI_photo_rate_array = 0
    pyramid_source_nx_py_pz_HI_photo_rate_array = 0
    pyramid_source_nx_py_nz_HI_photo_rate_array = 0
    pyramid_source_nx_ny_pz_HI_photo_rate_array = 0
    pyramid_source_nx_ny_nz_HI_photo_rate_array = 0  
	  	   
  end subroutine source_photo_rate_initialization
  
  subroutine px_py_pz_HI_photo_rate_calculation()
	
    integer :: i_level_py,a,b
    integer :: i_primary,i_secondary
    real(kind=dp) :: photo_rate_dx,photo_rate_dy,photo_rate_dz
	
    do i_level_py=1,level_py
	  do a=1,i_level_py
        do b=1,i_level_py
			  
			photo_rate_dx = 0
			photo_rate_dy = 0
			photo_rate_dz = 0
									
          do i_primary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1),&
  	                     ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1)	
            do i_secondary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2),&
                             ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2)
 
              photo_rate_dx = photo_rate_dx + &
			                  px_py_pz_dx(i_level_py)%HI_photoionization_rate(i_primary,i_secondary) * &
                              cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                              pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
              photo_rate_dy = photo_rate_dy + &
			                  px_py_pz_dy(i_level_py)%HI_photoionization_rate(i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume(i_primary,i_secondary)							  
              photo_rate_dz = photo_rate_dz + &
			                  px_py_pz_dz(i_level_py)%HI_photoionization_rate(i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume(i_primary,i_secondary)						  
				  							  
            enddo
          enddo
		  
		  pyramid_source_px_py_pz_HI_photo_rate_array(i_level_py,a,b) = &
		             pyramid_source_px_py_pz_HI_photo_rate_array(i_level_py,a,b)+photo_rate_dx
          pyramid_source_px_py_pz_HI_photo_rate_array(b,i_level_py,a) = &
		             pyramid_source_px_py_pz_HI_photo_rate_array(b,i_level_py,a)+photo_rate_dy
          pyramid_source_px_py_pz_HI_photo_rate_array(a,b,i_level_py) = &
		             pyramid_source_px_py_pz_HI_photo_rate_array(a,b,i_level_py)+photo_rate_dz
					 
        enddo
      enddo		
    enddo

  end subroutine px_py_pz_HI_photo_rate_calculation

  subroutine px_py_nz_HI_photo_rate_calculation()
	
    integer :: i_level_py,a,b
    integer :: i_primary,i_secondary
    real(kind=dp) :: photo_rate_dx,photo_rate_dy,photo_rate_dz
	
    do i_level_py=1,level_py
	  do a=1,i_level_py
        do b=1,i_level_py
			  
			photo_rate_dx = 0
			photo_rate_dy = 0
			photo_rate_dz = 0
									
          do i_primary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1),&
  	                     ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1)	
            do i_secondary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2),&
                             ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2)
 
              photo_rate_dx = photo_rate_dx + &
			                  px_py_nz_dx(i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) * &
                              cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                              pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
              photo_rate_dy = photo_rate_dy + &
			                  px_py_nz_dy(i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume(i_primary,i_secondary)							  
              photo_rate_dz = photo_rate_dz + &
			                  px_py_nz_dz(-i_level_py)%HI_photoionization_rate(i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume(i_primary,i_secondary)						  
				  							  
            enddo
          enddo
		  
		  pyramid_source_px_py_nz_HI_photo_rate_array(i_level_py,a,-b) = &
		             pyramid_source_px_py_nz_HI_photo_rate_array(i_level_py,a,-b)+photo_rate_dx
          pyramid_source_px_py_nz_HI_photo_rate_array(b,i_level_py,-a) = &
		             pyramid_source_px_py_nz_HI_photo_rate_array(b,i_level_py,-a)+photo_rate_dy
          pyramid_source_px_py_nz_HI_photo_rate_array(a,b,-i_level_py) = &
		             pyramid_source_px_py_nz_HI_photo_rate_array(a,b,-i_level_py)+photo_rate_dz
					 
        enddo
      enddo		
    enddo
 
  end subroutine px_py_nz_HI_photo_rate_calculation

  subroutine px_ny_pz_HI_photo_rate_calculation()
	
    integer :: i_level_py,a,b
    integer :: i_primary,i_secondary
    real(kind=dp) :: photo_rate_dx,photo_rate_dy,photo_rate_dz
	
    do i_level_py=1,level_py
	  do a=1,i_level_py
        do b=1,i_level_py
			  
			photo_rate_dx = 0
			photo_rate_dy = 0
			photo_rate_dz = 0
									
          do i_primary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1),&
  	                     ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1)	
            do i_secondary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2),&
                             ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2)
 
              photo_rate_dx = photo_rate_dx + &
			                  px_ny_pz_dx(i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) * &
                              cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                              pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
              photo_rate_dy = photo_rate_dy + &
			                  px_ny_pz_dy(-i_level_py)%HI_photoionization_rate(i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume(i_primary,i_secondary)							  
              photo_rate_dz = photo_rate_dz + &
			                  px_ny_pz_dz(i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume(i_primary,i_secondary)						  
				  							  
            enddo
          enddo
		  
		  pyramid_source_px_ny_pz_HI_photo_rate_array(i_level_py,-a,b) = &
		             pyramid_source_px_ny_pz_HI_photo_rate_array(i_level_py,-a,b)+photo_rate_dx
          pyramid_source_px_ny_pz_HI_photo_rate_array(b,-i_level_py,a) = &
		             pyramid_source_px_ny_pz_HI_photo_rate_array(b,-i_level_py,a)+photo_rate_dy
          pyramid_source_px_ny_pz_HI_photo_rate_array(a,-b,i_level_py) = &
		             pyramid_source_px_ny_pz_HI_photo_rate_array(a,-b,i_level_py)+photo_rate_dz
					 
        enddo
      enddo		
    enddo
	  
  end subroutine px_ny_pz_HI_photo_rate_calculation

  subroutine px_ny_nz_HI_photo_rate_calculation()
	
    integer :: i_level_py,a,b
    integer :: i_primary,i_secondary
    real(kind=dp) :: photo_rate_dx,photo_rate_dy,photo_rate_dz
	
    do i_level_py=1,level_py
	  do a=1,i_level_py
        do b=1,i_level_py
			  
			photo_rate_dx = 0
			photo_rate_dy = 0
			photo_rate_dz = 0
									
          do i_primary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1),&
  	                     ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1)	
            do i_secondary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2),&
                             ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2)
 
              photo_rate_dx = photo_rate_dx + &
			                  px_ny_nz_dx(i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) * &
                              cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                              pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
              photo_rate_dy = photo_rate_dy + &
			                  px_ny_nz_dy(-i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume(i_primary,i_secondary)							  
              photo_rate_dz = photo_rate_dz + &
			                  px_ny_nz_dz(-i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume(i_primary,i_secondary)						  
				  							  
            enddo
          enddo
		  
		  pyramid_source_px_ny_nz_HI_photo_rate_array(i_level_py,-a,-b) = &
		             pyramid_source_px_ny_nz_HI_photo_rate_array(i_level_py,-a,-b)+photo_rate_dx
          pyramid_source_px_ny_nz_HI_photo_rate_array(b,-i_level_py,-a) = &
		             pyramid_source_px_ny_nz_HI_photo_rate_array(b,-i_level_py,-a)+photo_rate_dy
          pyramid_source_px_ny_nz_HI_photo_rate_array(a,-b,-i_level_py) = &
		             pyramid_source_px_ny_nz_HI_photo_rate_array(a,-b,-i_level_py)+photo_rate_dz
					 
        enddo
      enddo		
    enddo
	  
  end subroutine px_ny_nz_HI_photo_rate_calculation

  subroutine nx_py_pz_HI_photo_rate_calculation()
	
    integer :: i_level_py,a,b
    integer :: i_primary,i_secondary
    real(kind=dp) :: photo_rate_dx,photo_rate_dy,photo_rate_dz
	
    do i_level_py=1,level_py
	  do a=1,i_level_py
        do b=1,i_level_py
			  
			photo_rate_dx = 0
			photo_rate_dy = 0
			photo_rate_dz = 0
									
          do i_primary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1),&
  	                     ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1)	
            do i_secondary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2),&
                             ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2)
 
              photo_rate_dx = photo_rate_dx + &
			                  nx_py_pz_dx(-i_level_py)%HI_photoionization_rate(i_primary,i_secondary) * &
                              cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                              pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
              photo_rate_dy = photo_rate_dy + &
			                  nx_py_pz_dy(i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume(i_primary,i_secondary)							  
              photo_rate_dz = photo_rate_dz + &
			                  nx_py_pz_dz(i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume(i_primary,i_secondary)						  
				  							  
            enddo
          enddo
		  
		  pyramid_source_nx_py_pz_HI_photo_rate_array(-i_level_py,a,b) = &
		             pyramid_source_nx_py_pz_HI_photo_rate_array(-i_level_py,a,b)+photo_rate_dx
          pyramid_source_nx_py_pz_HI_photo_rate_array(-b,i_level_py,a) = &
		             pyramid_source_nx_py_pz_HI_photo_rate_array(-b,i_level_py,a)+photo_rate_dy
          pyramid_source_nx_py_pz_HI_photo_rate_array(-a,b,i_level_py) = &
		             pyramid_source_nx_py_pz_HI_photo_rate_array(-a,b,i_level_py)+photo_rate_dz
					 
        enddo
      enddo		
    enddo
	  
  end subroutine nx_py_pz_HI_photo_rate_calculation

  subroutine nx_py_nz_HI_photo_rate_calculation()
	
    integer :: i_level_py,a,b
    integer :: i_primary,i_secondary
    real(kind=dp) :: photo_rate_dx,photo_rate_dy,photo_rate_dz
	
    do i_level_py=1,level_py
	  do a=1,i_level_py
        do b=1,i_level_py
			  
			photo_rate_dx = 0
			photo_rate_dy = 0
			photo_rate_dz = 0
									
          do i_primary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1),&
  	                     ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1)	
            do i_secondary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2),&
                             ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2)
 
              photo_rate_dx = photo_rate_dx + &
			                  nx_py_nz_dx(-i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) * &
                              cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                              pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
              photo_rate_dy = photo_rate_dy + &
			                  nx_py_nz_dy(i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume(i_primary,i_secondary)							  
              photo_rate_dz = photo_rate_dz + &
			                  nx_py_nz_dz(-i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume(i_primary,i_secondary)						  
				  							  
            enddo
          enddo
		  
		  pyramid_source_nx_py_nz_HI_photo_rate_array(-i_level_py,a,-b) = &
		             pyramid_source_nx_py_nz_HI_photo_rate_array(-i_level_py,a,-b)+photo_rate_dx
          pyramid_source_nx_py_nz_HI_photo_rate_array(-b,i_level_py,-a) = &
		             pyramid_source_nx_py_nz_HI_photo_rate_array(-b,i_level_py,-a)+photo_rate_dy
          pyramid_source_nx_py_nz_HI_photo_rate_array(-a,b,-i_level_py) = &
		             pyramid_source_nx_py_nz_HI_photo_rate_array(-a,b,-i_level_py)+photo_rate_dz
					 
        enddo
      enddo		
    enddo

  end subroutine nx_py_nz_HI_photo_rate_calculation

  subroutine nx_ny_pz_HI_photo_rate_calculation()
	
    integer :: i_level_py,a,b
    integer :: i_primary,i_secondary
    real(kind=dp) :: photo_rate_dx,photo_rate_dy,photo_rate_dz
	
    do i_level_py=1,level_py
	  do a=1,i_level_py
        do b=1,i_level_py
			  
			photo_rate_dx = 0
			photo_rate_dy = 0
			photo_rate_dz = 0
									
          do i_primary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1),&
  	                     ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1)	
            do i_secondary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2),&
                             ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2)
 
              photo_rate_dx = photo_rate_dx + &
			                  nx_ny_pz_dx(-i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) * &
                              cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                              pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
              photo_rate_dy = photo_rate_dy + &
			                  nx_ny_pz_dy(-i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume(i_primary,i_secondary)							  
              photo_rate_dz = photo_rate_dz + &
			                  nx_ny_pz_dz(i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume(i_primary,i_secondary)						  
				  							  
            enddo
          enddo
		  
		  pyramid_source_nx_ny_pz_HI_photo_rate_array(-i_level_py,-a,b) = &
		             pyramid_source_nx_ny_pz_HI_photo_rate_array(-i_level_py,-a,b)+photo_rate_dx
          pyramid_source_nx_ny_pz_HI_photo_rate_array(-b,-i_level_py,a) = &
		             pyramid_source_nx_ny_pz_HI_photo_rate_array(-b,-i_level_py,a)+photo_rate_dy
          pyramid_source_nx_ny_pz_HI_photo_rate_array(-a,-b,i_level_py) = &
		             pyramid_source_nx_ny_pz_HI_photo_rate_array(-a,-b,i_level_py)+photo_rate_dz
					 
        enddo
      enddo		
    enddo
	  
  end subroutine nx_ny_pz_HI_photo_rate_calculation

  subroutine nx_ny_nz_HI_photo_rate_calculation()
	
    integer :: i_level_py,a,b
    integer :: i_primary,i_secondary
    real(kind=dp) :: photo_rate_dx,photo_rate_dy,photo_rate_dz
	
    do i_level_py=1,level_py
	  do a=1,i_level_py
        do b=1,i_level_py
			  
			photo_rate_dx = 0
			photo_rate_dy = 0
			photo_rate_dz = 0
									
          do i_primary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1),&
  	                     ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,1)	
            do i_secondary = lbound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2),&
                             ubound(cartesian_grid(i_level_py)%transform(a,b)%transform_volume,2)
 
              photo_rate_dx = photo_rate_dx + &
			                  nx_ny_nz_dx(-i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) * &
                              cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
                              pyramid_grid(i_level_py)%volume(i_primary,i_secondary)
              photo_rate_dy = photo_rate_dy + &
			                  nx_ny_nz_dy(-i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume(i_primary,i_secondary)							  
              photo_rate_dz = photo_rate_dz + &
			                  nx_ny_nz_dz(-i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) * &
			                  cartesian_grid(i_level_py)%transform(a,b)%transform_volume(i_primary,i_secondary) / &
			                  pyramid_grid(i_level_py)%volume(i_primary,i_secondary)						  
				  							  
            enddo
          enddo
		  
		  pyramid_source_nx_ny_nz_HI_photo_rate_array(-i_level_py,-a,-b) = &
		             pyramid_source_nx_ny_nz_HI_photo_rate_array(-i_level_py,-a,-b)+photo_rate_dx
          pyramid_source_nx_ny_nz_HI_photo_rate_array(-b,-i_level_py,-a) = &
		             pyramid_source_nx_ny_nz_HI_photo_rate_array(-b,-i_level_py,-a)+photo_rate_dy
          pyramid_source_nx_ny_nz_HI_photo_rate_array(-a,-b,-i_level_py) = &
		             pyramid_source_nx_ny_nz_HI_photo_rate_array(-a,-b,-i_level_py)+photo_rate_dz
					 
        enddo
      enddo		
    enddo
		  
  end subroutine nx_ny_nz_HI_photo_rate_calculation
  
  subroutine photo_rate_per_atom()
	  
	pyramid_source_px_py_pz_HI_photo_rate_array=&
	pyramid_source_px_py_pz_HI_photo_rate_array/(pyramid_source_px_py_pz_HI_density_array*cellsize_py_cube)	
	pyramid_source_px_py_nz_HI_photo_rate_array=&
	pyramid_source_px_py_nz_HI_photo_rate_array/(pyramid_source_px_py_nz_HI_density_array*cellsize_py_cube)	
	pyramid_source_px_ny_pz_HI_photo_rate_array=&
	pyramid_source_px_ny_pz_HI_photo_rate_array/(pyramid_source_px_ny_pz_HI_density_array*cellsize_py_cube)	
	pyramid_source_px_ny_nz_HI_photo_rate_array=&
	pyramid_source_px_ny_nz_HI_photo_rate_array/(pyramid_source_px_ny_nz_HI_density_array*cellsize_py_cube)	
	pyramid_source_nx_py_pz_HI_photo_rate_array=&
	pyramid_source_nx_py_pz_HI_photo_rate_array/(pyramid_source_nx_py_pz_HI_density_array*cellsize_py_cube)	
	pyramid_source_nx_py_nz_HI_photo_rate_array=&
	pyramid_source_nx_py_nz_HI_photo_rate_array/(pyramid_source_nx_py_nz_HI_density_array*cellsize_py_cube)	
	pyramid_source_nx_ny_pz_HI_photo_rate_array=&
	pyramid_source_nx_ny_pz_HI_photo_rate_array/(pyramid_source_nx_ny_pz_HI_density_array*cellsize_py_cube)	
	pyramid_source_nx_ny_nz_HI_photo_rate_array=&
	pyramid_source_nx_ny_nz_HI_photo_rate_array/(pyramid_source_nx_ny_nz_HI_density_array*cellsize_py_cube)		  	  

		  	  
  end subroutine photo_rate_per_atom
end module global_photo_rate