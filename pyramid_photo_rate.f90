module pyramid_photo_rate
	
    use precision, only: dp
    use input, only: level_py, partition ,volume_factor, use_which_source, pi
    use array, only: pyramid_grid, &
                     px_py_pz_dx, px_py_pz_dy, px_py_pz_dz, &
                     px_py_nz_dx, px_py_nz_dy, px_py_nz_dz, &
                     px_ny_pz_dx, px_ny_pz_dy, px_ny_pz_dz, &	  
                     px_ny_nz_dx, px_ny_nz_dy, px_ny_nz_dz, & 	   
                     nx_py_pz_dx, nx_py_pz_dy, nx_py_pz_dz, &	  
                     nx_py_nz_dx, nx_py_nz_dy, nx_py_nz_dz, &	
                     nx_ny_pz_dx, nx_ny_pz_dy, nx_ny_pz_dz, &	  
                     nx_ny_nz_dx, nx_ny_nz_dy, nx_ny_nz_dz		   	   
    use radiation, only: photoion_solid_angle, photrates
		
contains
	
  subroutine pyramid_domain_photo_rate_calculation()
		
    call px_py_pz_dx_photo_rate_calculation ()
    call px_py_pz_dy_photo_rate_calculation ()
    call px_py_pz_dz_photo_rate_calculation ()	
    call px_py_nz_dx_photo_rate_calculation ()
    call px_py_nz_dy_photo_rate_calculation ()	
    call px_py_nz_dz_photo_rate_calculation ()		
    call px_ny_pz_dx_photo_rate_calculation ()	
    call px_ny_pz_dy_photo_rate_calculation ()	
    call px_ny_pz_dz_photo_rate_calculation ()	
    call px_ny_nz_dx_photo_rate_calculation ()
    call px_ny_nz_dy_photo_rate_calculation ()
    call px_ny_nz_dz_photo_rate_calculation ()
    call nx_py_pz_dx_photo_rate_calculation ()
    call nx_py_pz_dy_photo_rate_calculation ()
    call nx_py_pz_dz_photo_rate_calculation ()
    call nx_py_nz_dx_photo_rate_calculation ()	
    call nx_py_nz_dy_photo_rate_calculation ()	
    call nx_py_nz_dz_photo_rate_calculation ()	
    call nx_ny_pz_dx_photo_rate_calculation ()
    call nx_ny_pz_dy_photo_rate_calculation ()
    call nx_ny_pz_dz_photo_rate_calculation ()	
    call nx_ny_nz_dx_photo_rate_calculation ()
    call nx_ny_nz_dy_photo_rate_calculation ()
    call nx_ny_nz_dz_photo_rate_calculation ()		

  end subroutine pyramid_domain_photo_rate_calculation	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  subroutine px_py_pz_dx_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
    real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_out
    real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_out		
    type(photrates) :: phi		
    real(kind=dp) :: i,j,k,x,y,z,n	

    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          i = real(i_level_py)
          j = real(i_primary)
          k = real(i_secondary)
          n = real(partition(i_level_py))

          x = (2.0*i-1.0)/2.0
          y = (2.0*j-1.0)*(2.0*i-1.0)/(4.0*n)
          z = (2.0*k-1.0)*(2.0*i-1.0)/(4.0*n)

          HI_column_density_in = px_py_pz_dx(i_level_py)%HI_column_density_in(i_primary,i_secondary)
          HeI_column_density_in = px_py_pz_dx(i_level_py)%HeI_column_density_in(i_primary,i_secondary)
          HeII_column_density_in = px_py_pz_dx(i_level_py)%HeII_column_density_in(i_primary,i_secondary)
          HI_column_density_out = px_py_pz_dx(i_level_py)%HI_column_density_out(i_primary,i_secondary)	  
          HeI_column_density_out = px_py_pz_dx(i_level_py)%HeI_column_density_out(i_primary,i_secondary)	  
          HeII_column_density_out = px_py_pz_dx(i_level_py)%HeII_column_density_out(i_primary,i_secondary)  	  
          !call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
	  !	HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,correction (x,y,z))
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
		  HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,1.0_dp)
          px_py_pz_dx(i_level_py)%HI_photoionization_rate(i_primary,i_secondary) = 1.0!phi%photo_cell_HI
          px_py_pz_dx(i_level_py)%HeI_photoionization_rate(i_primary,i_secondary) = phi%photo_cell_HeI
          px_py_pz_dx(i_level_py)%HeII_photoionization_rate(i_primary,i_secondary) = phi%photo_cell_HeII
          px_py_pz_dx(i_level_py)%photoheating_rate(i_primary,i_secondary) = phi%heat

        enddo
      enddo	  
    enddo
write(*,*)'py HI rate near source is ',  px_py_pz_dx(1)%HI_photoionization_rate(1,1)  
  end subroutine px_py_pz_dx_photo_rate_calculation

  subroutine px_py_pz_dy_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
    real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_out
    real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_out
    type(photrates) :: phi		
    real(kind=dp) :: i,j,k,x,y,z,n		
	
		
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          i = real(i_secondary)
          j = real(i_level_py)
          k = real(i_primary)
          n = real(partition(i_level_py))

          x = (2.0*i-1.0)*(2.0*j-1.0)/(4.0*n)
          y = (2.0*j-1.0)/2.0
          z = (2.0*k-1.0)*(2.0*j-1.0)/(4.0*n)	

          HI_column_density_in = px_py_pz_dy(i_level_py)%HI_column_density_in(i_primary,i_secondary)
          HeI_column_density_in = px_py_pz_dy(i_level_py)%HeI_column_density_in(i_primary,i_secondary)
          HeII_column_density_in = px_py_pz_dy(i_level_py)%HeII_column_density_in(i_primary,i_secondary)
          HI_column_density_out = px_py_pz_dy(i_level_py)%HI_column_density_out(i_primary,i_secondary)	  
          HeI_column_density_out = px_py_pz_dy(i_level_py)%HeI_column_density_out(i_primary,i_secondary)	  
          HeII_column_density_out = px_py_pz_dy(i_level_py)%HeII_column_density_out(i_primary,i_secondary)	  
          !call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
	  !	HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,correction (x,y,z))
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
		  HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,1.0_dp)
          px_py_pz_dy(i_level_py)%HI_photoionization_rate(i_primary,i_secondary) = 1.0!phi%photo_cell_HI
          px_py_pz_dy(i_level_py)%HeI_photoionization_rate(i_primary,i_secondary) = phi%photo_cell_HeI
          px_py_pz_dy(i_level_py)%HeII_photoionization_rate(i_primary,i_secondary) = phi%photo_cell_HeII
          px_py_pz_dy(i_level_py)%photoheating_rate(i_primary,i_secondary) = phi%heat

        enddo
      enddo	
    enddo
  	  
  end subroutine px_py_pz_dy_photo_rate_calculation

  subroutine px_py_pz_dz_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
    real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_out
    real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_out
    type(photrates) :: phi			
    real(kind=dp) :: i,j,k,x,y,z,n	
	
		
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          i = real(i_primary)
          j = real(i_secondary)
          k = real(i_level_py)
          n = real(partition(i_level_py))

          x = (2.0*i-1.0)*(2.0*k-1.0)/(4.0*n)
          y = (2.0*j-1.0)*(2.0*k-1.0)/(4.0*n)
          z = (2.0*k-1.0)/2.0	

          HI_column_density_in = px_py_pz_dz(i_level_py)%HI_column_density_in(i_primary,i_secondary)
          HeI_column_density_in = px_py_pz_dz(i_level_py)%HeI_column_density_in(i_primary,i_secondary)
          HeII_column_density_in = px_py_pz_dz(i_level_py)%HeII_column_density_in(i_primary,i_secondary)
          HI_column_density_out = px_py_pz_dz(i_level_py)%HI_column_density_out(i_primary,i_secondary)	  
          HeI_column_density_out = px_py_pz_dz(i_level_py)%HeI_column_density_out(i_primary,i_secondary)	  
          HeII_column_density_out = px_py_pz_dz(i_level_py)%HeII_column_density_out(i_primary,i_secondary)	  
          !call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
	  !	HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,correction (x,y,z))
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
		  HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,1.0_dp)
          px_py_pz_dz(i_level_py)%HI_photoionization_rate(i_primary,i_secondary) = 1.0!phi%photo_cell_HI
          px_py_pz_dz(i_level_py)%HeI_photoionization_rate(i_primary,i_secondary) = phi%photo_cell_HeI
          px_py_pz_dz(i_level_py)%HeII_photoionization_rate(i_primary,i_secondary) = phi%photo_cell_HeII
          px_py_pz_dz(i_level_py)%photoheating_rate(i_primary,i_secondary) = phi%heat

        enddo
      enddo	
    enddo
  	  
  end subroutine px_py_pz_dz_photo_rate_calculation
    
  subroutine px_py_nz_dx_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
    real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_out
    real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_out
    type(photrates) :: phi			
    real(kind=dp) :: i,j,k,x,y,z,n	
	

    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          i = real(i_level_py)
          j = real(i_primary)
          k = real(i_secondary)
          n = real(partition(i_level_py))

          x = (2.0*i-1.0)/2.0
          y = (2.0*j-1.0)*(2.0*i-1.0)/(4.0*n)
          z = -(2.0*k-1.0)*(2.0*i-1.0)/(4.0*n)

          HI_column_density_in = px_py_nz_dx(i_level_py)%HI_column_density_in(i_primary,-i_secondary)
          HeI_column_density_in = px_py_nz_dx(i_level_py)%HeI_column_density_in(i_primary,-i_secondary)
          HeII_column_density_in = px_py_nz_dx(i_level_py)%HeII_column_density_in(i_primary,-i_secondary)
          HI_column_density_out = px_py_nz_dx(i_level_py)%HI_column_density_out(i_primary,-i_secondary)	  
          HeI_column_density_out = px_py_nz_dx(i_level_py)%HeI_column_density_out(i_primary,-i_secondary)	  
          HeII_column_density_out = px_py_nz_dx(i_level_py)%HeII_column_density_out(i_primary,-i_secondary)  	  
          !call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
	  !	HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,correction (x,y,z))
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
		  HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,1.0_dp)
          px_py_nz_dx(i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) = 1.0!phi%photo_cell_HI
          px_py_nz_dx(i_level_py)%HeI_photoionization_rate(i_primary,-i_secondary) = phi%photo_cell_HeI
          px_py_nz_dx(i_level_py)%HeII_photoionization_rate(i_primary,-i_secondary) = phi%photo_cell_HeII
          px_py_nz_dx(i_level_py)%photoheating_rate(i_primary,-i_secondary) = phi%heat

        enddo
      enddo	  
    enddo
  	  
  end subroutine px_py_nz_dx_photo_rate_calculation

  subroutine px_py_nz_dy_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
    real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_out
    real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_out
    type(photrates) :: phi			
    real(kind=dp) :: i,j,k,x,y,z,n	
	
		
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          i = real(i_secondary)
          j = real(i_level_py)
          k = real(i_primary)
          n = real(partition(i_level_py))

          x = (2.0*i-1.0)*(2.0*j-1.0)/(4.0*n)
          y = (2.0*j-1.0)/2.0
          z = -(2.0*k-1.0)*(2.0*j-1.0)/(4.0*n)	

          HI_column_density_in = px_py_nz_dy(i_level_py)%HI_column_density_in(-i_primary,i_secondary)
          HeI_column_density_in = px_py_nz_dy(i_level_py)%HeI_column_density_in(-i_primary,i_secondary)
          HeII_column_density_in = px_py_nz_dy(i_level_py)%HeII_column_density_in(-i_primary,i_secondary)
          HI_column_density_out = px_py_nz_dy(i_level_py)%HI_column_density_out(-i_primary,i_secondary)	  
          HeI_column_density_out = px_py_nz_dy(i_level_py)%HeI_column_density_out(-i_primary,i_secondary)	  
          HeII_column_density_out = px_py_nz_dy(i_level_py)%HeII_column_density_out(-i_primary,i_secondary)  	  
          !call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
	  !	HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,correction (x,y,z))
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
		  HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,1.0_dp)
          px_py_nz_dy(i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) = 1.0!phi%photo_cell_HI
          px_py_nz_dy(i_level_py)%HeI_photoionization_rate(-i_primary,i_secondary) = phi%photo_cell_HeI
          px_py_nz_dy(i_level_py)%HeII_photoionization_rate(-i_primary,i_secondary) = phi%photo_cell_HeII
          px_py_nz_dy(i_level_py)%photoheating_rate(-i_primary,i_secondary) = phi%heat

        enddo
      enddo	  
    enddo
  	  
  end subroutine px_py_nz_dy_photo_rate_calculation

  subroutine px_py_nz_dz_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
    real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_out
    real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_out
    type(photrates) :: phi			
    real(kind=dp) :: i,j,k,x,y,z,n	
	
		
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          i = real(i_primary)
          j = real(i_secondary)
          k = real(i_level_py)
          n = real(partition(i_level_py))

          x = (2.0*i-1.0)*(2.0*k-1.0)/(4.0*n)
          y = (2.0*j-1.0)*(2.0*k-1.0)/(4.0*n)
          z = -(2.0*k-1.0)/2.0	

          HI_column_density_in = px_py_nz_dz(-i_level_py)%HI_column_density_in(i_primary,i_secondary)
          HeI_column_density_in = px_py_nz_dz(-i_level_py)%HeI_column_density_in(i_primary,i_secondary)
          HeII_column_density_in = px_py_nz_dz(-i_level_py)%HeII_column_density_in(i_primary,i_secondary)
          HI_column_density_out = px_py_nz_dz(-i_level_py)%HI_column_density_out(i_primary,i_secondary)	  
          HeI_column_density_out = px_py_nz_dz(-i_level_py)%HeI_column_density_out(i_primary,i_secondary)	  
          HeII_column_density_out = px_py_nz_dz(-i_level_py)%HeII_column_density_out(i_primary,i_secondary)  	  
          !call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
	  !	HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,correction (x,y,z))
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
		  HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,1.0_dp)
          px_py_nz_dz(-i_level_py)%HI_photoionization_rate(i_primary,i_secondary) = 1.0!phi%photo_cell_HI
          px_py_nz_dz(-i_level_py)%HeI_photoionization_rate(i_primary,i_secondary) = phi%photo_cell_HeI
          px_py_nz_dz(-i_level_py)%HeII_photoionization_rate(i_primary,i_secondary) = phi%photo_cell_HeII
          px_py_nz_dz(-i_level_py)%photoheating_rate(i_primary,i_secondary) = phi%heat

        enddo
      enddo	  
    enddo
  	  
  end subroutine px_py_nz_dz_photo_rate_calculation
  
  subroutine px_ny_pz_dx_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
    real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_out
    real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_out
    type(photrates) :: phi			
    real(kind=dp) :: i,j,k,x,y,z,n	
	
		
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          i = real(i_level_py)
          j = real(i_primary)
          k = real(i_secondary)
          n = real(partition(i_level_py))

          x = (2.0*i-1.0)/2.0
          y = -(2.0*j-1.0)*(2.0*i-1.0)/(4.0*n)
          z = (2.0*k-1.0)*(2.0*i-1.0)/(4.0*n)

          HI_column_density_in = px_ny_pz_dx(i_level_py)%HI_column_density_in(-i_primary,i_secondary)
          HeI_column_density_in = px_ny_pz_dx(i_level_py)%HeI_column_density_in(-i_primary,i_secondary)
          HeII_column_density_in = px_ny_pz_dx(i_level_py)%HeII_column_density_in(-i_primary,i_secondary)
          HI_column_density_out = px_ny_pz_dx(i_level_py)%HI_column_density_out(-i_primary,i_secondary)	  
          HeI_column_density_out = px_ny_pz_dx(i_level_py)%HeI_column_density_out(-i_primary,i_secondary)	  
          HeII_column_density_out = px_ny_pz_dx(i_level_py)%HeII_column_density_out(-i_primary,i_secondary)  	  
          !call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
	  !	HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,correction (x,y,z))
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
		  HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,1.0_dp)
          px_ny_pz_dx(i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) = 1.0!phi%photo_cell_HI
          px_ny_pz_dx(i_level_py)%HeI_photoionization_rate(-i_primary,i_secondary) = phi%photo_cell_HeI
          px_ny_pz_dx(i_level_py)%HeII_photoionization_rate(-i_primary,i_secondary) = phi%photo_cell_HeII
          px_ny_pz_dx(i_level_py)%photoheating_rate(-i_primary,i_secondary) = phi%heat

        enddo
      enddo	  
    enddo
  	  
  end subroutine px_ny_pz_dx_photo_rate_calculation

  subroutine px_ny_pz_dy_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
    real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_out
    real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_out
    type(photrates) :: phi			
    real(kind=dp) :: i,j,k,x,y,z,n	
	
		
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          i = real(i_secondary)
          j = real(i_level_py)
          k = real(i_primary)
          n = real(partition(i_level_py))

          x = (2.0*i-1.0)*(2.0*j-1.0)/(4.0*n)
          y = -(2.0*j-1.0)/2.0
          z = (2.0*k-1.0)*(2.0*j-1.0)/(4.0*n)	

          HI_column_density_in = px_ny_pz_dy(-i_level_py)%HI_column_density_in(i_primary,i_secondary)
          HeI_column_density_in = px_ny_pz_dy(-i_level_py)%HeI_column_density_in(i_primary,i_secondary)
          HeII_column_density_in = px_ny_pz_dy(-i_level_py)%HeII_column_density_in(i_primary,i_secondary)
          HI_column_density_out = px_ny_pz_dy(-i_level_py)%HI_column_density_out(i_primary,i_secondary)	  
          HeI_column_density_out = px_ny_pz_dy(-i_level_py)%HeI_column_density_out(i_primary,i_secondary)	  
          HeII_column_density_out = px_ny_pz_dy(-i_level_py)%HeII_column_density_out(i_primary,i_secondary)  	  
          !call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
	  !	HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,correction (x,y,z))
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
		  HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,1.0_dp)
          px_ny_pz_dy(-i_level_py)%HI_photoionization_rate(i_primary,i_secondary) = 1.0!phi%photo_cell_HI
          px_ny_pz_dy(-i_level_py)%HeI_photoionization_rate(i_primary,i_secondary) = phi%photo_cell_HeI
          px_ny_pz_dy(-i_level_py)%HeII_photoionization_rate(i_primary,i_secondary) = phi%photo_cell_HeII
          px_ny_pz_dy(-i_level_py)%photoheating_rate(i_primary,i_secondary) = phi%heat

        enddo
      enddo	  
    enddo
  	  
  end subroutine px_ny_pz_dy_photo_rate_calculation

  subroutine px_ny_pz_dz_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
    real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_out
    real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_out
    type(photrates) :: phi			
    real(kind=dp) :: i,j,k,x,y,z,n	
	
		
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          i = real(i_primary)
          j = real(i_secondary)
          k = real(i_level_py)
          n = real(partition(i_level_py))

          x = (2.0*i-1.0)*(2.0*k-1.0)/(4.0*n)
          y = -(2.0*j-1.0)*(2.0*k-1.0)/(4.0*n)
          z = (2.0*k-1.0)/2.0	

          HI_column_density_in = px_ny_pz_dz(i_level_py)%HI_column_density_in(i_primary,-i_secondary)
          HeI_column_density_in = px_ny_pz_dz(i_level_py)%HeI_column_density_in(i_primary,-i_secondary)
          HeII_column_density_in = px_ny_pz_dz(i_level_py)%HeII_column_density_in(i_primary,-i_secondary)
          HI_column_density_out = px_ny_pz_dz(i_level_py)%HI_column_density_out(i_primary,-i_secondary)	  
          HeI_column_density_out = px_ny_pz_dz(i_level_py)%HeI_column_density_out(i_primary,-i_secondary)	  
          HeII_column_density_out = px_ny_pz_dz(i_level_py)%HeII_column_density_out(i_primary,-i_secondary)  	  
          !call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
	  !	HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,correction (x,y,z))
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
		  HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,1.0_dp)
          px_ny_pz_dz(i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) = 1.0!phi%photo_cell_HI
          px_ny_pz_dz(i_level_py)%HeI_photoionization_rate(i_primary,-i_secondary) = phi%photo_cell_HeI
          px_ny_pz_dz(i_level_py)%HeII_photoionization_rate(i_primary,-i_secondary) = phi%photo_cell_HeII
          px_ny_pz_dz(i_level_py)%photoheating_rate(i_primary,-i_secondary) = phi%heat

        enddo
      enddo	  
    enddo
  	  
  end subroutine px_ny_pz_dz_photo_rate_calculation
    
  subroutine px_ny_nz_dx_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
    real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_out
    real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_out
    type(photrates) :: phi			
    real(kind=dp) :: i,j,k,x,y,z,n	
	
		
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          i = real(i_level_py)
          j = real(i_primary)
          k = real(i_secondary)
          n = real(partition(i_level_py))

          x = (2.0*i-1.0)/2.0
          y = -(2.0*j-1.0)*(2.0*i-1.0)/(4.0*n)
          z = -(2.0*k-1.0)*(2.0*i-1.0)/(4.0*n)

          HI_column_density_in = px_ny_nz_dx(i_level_py)%HI_column_density_in(-i_primary,-i_secondary)
          HeI_column_density_in = px_ny_nz_dx(i_level_py)%HeI_column_density_in(-i_primary,-i_secondary)
          HeII_column_density_in = px_ny_nz_dx(i_level_py)%HeII_column_density_in(-i_primary,-i_secondary)
          HI_column_density_out = px_ny_nz_dx(i_level_py)%HI_column_density_out(-i_primary,-i_secondary)	  
          HeI_column_density_out = px_ny_nz_dx(i_level_py)%HeI_column_density_out(-i_primary,-i_secondary)	  
          HeII_column_density_out = px_ny_nz_dx(i_level_py)%HeII_column_density_out(-i_primary,-i_secondary)  	  
          !call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
	  !	HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,correction (x,y,z))
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
		  HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,1.0_dp)
          px_ny_nz_dx(i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) = 1.0!phi%photo_cell_HI
          px_ny_nz_dx(i_level_py)%HeI_photoionization_rate(-i_primary,-i_secondary) = phi%photo_cell_HeI
          px_ny_nz_dx(i_level_py)%HeII_photoionization_rate(-i_primary,-i_secondary) = phi%photo_cell_HeII
          px_ny_nz_dx(i_level_py)%photoheating_rate(-i_primary,-i_secondary) = phi%heat

        enddo
      enddo	  
    enddo
  	  
  end subroutine px_ny_nz_dx_photo_rate_calculation

  subroutine px_ny_nz_dy_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
    real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_out
    real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_out
    type(photrates) :: phi			
    real(kind=dp) :: i,j,k,x,y,z,n	
	
		
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          i = real(i_secondary)
          j = real(i_level_py)
          k = real(i_primary)
          n = real(partition(i_level_py))

          x = (2.0*i-1.0)*(2.0*j-1.0)/(4.0*n)
          y = -(2.0*j-1.0)/2.0
          z = -(2.0*k-1.0)*(2.0*j-1.0)/(4.0*n)	

          HI_column_density_in = px_ny_nz_dy(-i_level_py)%HI_column_density_in(-i_primary,i_secondary)
          HeI_column_density_in = px_ny_nz_dy(-i_level_py)%HeI_column_density_in(-i_primary,i_secondary)
          HeII_column_density_in = px_ny_nz_dy(-i_level_py)%HeII_column_density_in(-i_primary,i_secondary)
          HI_column_density_out = px_ny_nz_dy(-i_level_py)%HI_column_density_out(-i_primary,i_secondary)	  
          HeI_column_density_out = px_ny_nz_dy(-i_level_py)%HeI_column_density_out(-i_primary,i_secondary)	  
          HeII_column_density_out = px_ny_nz_dy(-i_level_py)%HeII_column_density_out(-i_primary,i_secondary)  	  
          !call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
	  !	HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,correction (x,y,z))
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
		  HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,1.0_dp)
          px_ny_nz_dy(-i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) = 1.0!phi%photo_cell_HI
          px_ny_nz_dy(-i_level_py)%HeI_photoionization_rate(-i_primary,i_secondary) = phi%photo_cell_HeI
          px_ny_nz_dy(-i_level_py)%HeII_photoionization_rate(-i_primary,i_secondary) = phi%photo_cell_HeII
          px_ny_nz_dy(-i_level_py)%photoheating_rate(-i_primary,i_secondary) = phi%heat

        enddo
      enddo	  
    enddo
  	  
  end subroutine px_ny_nz_dy_photo_rate_calculation

  subroutine px_ny_nz_dz_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
    real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_out
    real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_out
    type(photrates) :: phi			
    real(kind=dp) :: i,j,k,x,y,z,n	
	
		
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          i = real(i_primary)
          j = real(i_secondary)
          k = real(i_level_py)
          n = real(partition(i_level_py))

          x = (2.0*i-1.0)*(2.0*k-1.0)/(4.0*n)
          y = -(2.0*j-1.0)*(2.0*k-1.0)/(4.0*n)
          z = -(2.0*k-1.0)/2.0	

          HI_column_density_in = px_ny_nz_dz(-i_level_py)%HI_column_density_in(i_primary,-i_secondary)
          HeI_column_density_in = px_ny_nz_dz(-i_level_py)%HeI_column_density_in(i_primary,-i_secondary)
          HeII_column_density_in = px_ny_nz_dz(-i_level_py)%HeII_column_density_in(i_primary,-i_secondary)
          HI_column_density_out = px_ny_nz_dz(-i_level_py)%HI_column_density_out(i_primary,-i_secondary)	  
          HeI_column_density_out = px_ny_nz_dz(-i_level_py)%HeI_column_density_out(i_primary,-i_secondary)	  
          HeII_column_density_out = px_ny_nz_dz(-i_level_py)%HeII_column_density_out(i_primary,-i_secondary)  	  
          !call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
	  !	HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,correction (x,y,z))
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
		  HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,1.0_dp)
          px_ny_nz_dz(-i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) = 1.0!phi%photo_cell_HI
          px_ny_nz_dz(-i_level_py)%HeI_photoionization_rate(i_primary,-i_secondary) = phi%photo_cell_HeI
          px_ny_nz_dz(-i_level_py)%HeII_photoionization_rate(i_primary,-i_secondary) = phi%photo_cell_HeII
          px_ny_nz_dz(-i_level_py)%photoheating_rate(i_primary,-i_secondary) = phi%heat

        enddo
      enddo	  
    enddo
  	  
  end subroutine px_ny_nz_dz_photo_rate_calculation
  
  subroutine nx_py_pz_dx_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
    real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_out
    real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_out
    type(photrates) :: phi			
    real(kind=dp) :: i,j,k,x,y,z,n	
	
		
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          i = real(i_level_py)
          j = real(i_primary)
          k = real(i_secondary)
          n = real(partition(i_level_py))

          x = -(2.0*i-1.0)/2.0
          y = (2.0*j-1.0)*(2.0*i-1.0)/(4.0*n)
          z = (2.0*k-1.0)*(2.0*i-1.0)/(4.0*n)

          HI_column_density_in = nx_py_pz_dx(-i_level_py)%HI_column_density_in(i_primary,i_secondary)
          HeI_column_density_in = nx_py_pz_dx(-i_level_py)%HeI_column_density_in(i_primary,i_secondary)
          HeII_column_density_in = nx_py_pz_dx(-i_level_py)%HeII_column_density_in(i_primary,i_secondary)
          HI_column_density_out = nx_py_pz_dx(-i_level_py)%HI_column_density_out(i_primary,i_secondary)	  
          HeI_column_density_out = nx_py_pz_dx(-i_level_py)%HeI_column_density_out(i_primary,i_secondary)	  
          HeII_column_density_out = nx_py_pz_dx(-i_level_py)%HeII_column_density_out(i_primary,i_secondary)  	  
          !call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
	  !	HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,correction (x,y,z))
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
		  HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,1.0_dp)
          nx_py_pz_dx(-i_level_py)%HI_photoionization_rate(i_primary,i_secondary) = 1.0!phi%photo_cell_HI
          nx_py_pz_dx(-i_level_py)%HeI_photoionization_rate(i_primary,i_secondary) = phi%photo_cell_HeI
          nx_py_pz_dx(-i_level_py)%HeII_photoionization_rate(i_primary,i_secondary) = phi%photo_cell_HeII
          nx_py_pz_dx(-i_level_py)%photoheating_rate(i_primary,i_secondary) = phi%heat

        enddo
      enddo	  
    enddo
  	  
  end subroutine nx_py_pz_dx_photo_rate_calculation

  subroutine nx_py_pz_dy_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
    real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_out
    real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_out
    type(photrates) :: phi		
    real(kind=dp) :: i,j,k,x,y,z,n	
		
		
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          i = real(i_secondary)
          j = real(i_level_py)
          k = real(i_primary)
          n = real(partition(i_level_py))

          x = -(2.0*i-1.0)*(2.0*j-1.0)/(4.0*n)
          y = (2.0*j-1.0)/2.0
          z = (2.0*k-1.0)*(2.0*j-1.0)/(4.0*n)	

          HI_column_density_in = nx_py_pz_dy(i_level_py)%HI_column_density_in(i_primary,-i_secondary)
          HeI_column_density_in = nx_py_pz_dy(i_level_py)%HeI_column_density_in(i_primary,-i_secondary)
          HeII_column_density_in = nx_py_pz_dy(i_level_py)%HeII_column_density_in(i_primary,-i_secondary)
          HI_column_density_out = nx_py_pz_dy(i_level_py)%HI_column_density_out(i_primary,-i_secondary)	  
          HeI_column_density_out = nx_py_pz_dy(i_level_py)%HeI_column_density_out(i_primary,-i_secondary)	  
          HeII_column_density_out = nx_py_pz_dy(i_level_py)%HeII_column_density_out(i_primary,-i_secondary)  	  
          !call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
	  !	HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,correction (x,y,z))
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
		  HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,1.0_dp)
          nx_py_pz_dy(i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) = 1.0!phi%photo_cell_HI
          nx_py_pz_dy(i_level_py)%HeI_photoionization_rate(i_primary,-i_secondary) = phi%photo_cell_HeI
          nx_py_pz_dy(i_level_py)%HeII_photoionization_rate(i_primary,-i_secondary) = phi%photo_cell_HeII
          nx_py_pz_dy(i_level_py)%photoheating_rate(i_primary,-i_secondary) = phi%heat

        enddo
      enddo	  
    enddo
  	  
  end subroutine nx_py_pz_dy_photo_rate_calculation

  subroutine nx_py_pz_dz_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
    real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_out
    real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_out
    type(photrates) :: phi			
    real(kind=dp) :: i,j,k,x,y,z,n	
	
		
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          i = real(i_primary)
          j = real(i_secondary)
          k = real(i_level_py)
          n = real(partition(i_level_py))

          x = -(2.0*i-1.0)*(2.0*k-1.0)/(4.0*n)
          y = (2.0*j-1.0)*(2.0*k-1.0)/(4.0*n)
          z = (2.0*k-1.0)/2.0	

          HI_column_density_in = nx_py_pz_dz(i_level_py)%HI_column_density_in(-i_primary,i_secondary)
          HeI_column_density_in = nx_py_pz_dz(i_level_py)%HeI_column_density_in(-i_primary,i_secondary)
          HeII_column_density_in = nx_py_pz_dz(i_level_py)%HeII_column_density_in(-i_primary,i_secondary)
          HI_column_density_out = nx_py_pz_dz(i_level_py)%HI_column_density_out(-i_primary,i_secondary)	  
          HeI_column_density_out = nx_py_pz_dz(i_level_py)%HeI_column_density_out(-i_primary,i_secondary)	  
          HeII_column_density_out = nx_py_pz_dz(i_level_py)%HeII_column_density_out(-i_primary,i_secondary)  	  
          !call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
	  !	HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,correction (x,y,z))
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
		  HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,1.0_dp)
          nx_py_pz_dz(i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) = 1.0!phi%photo_cell_HI
          nx_py_pz_dz(i_level_py)%HeI_photoionization_rate(-i_primary,i_secondary) = phi%photo_cell_HeI
          nx_py_pz_dz(i_level_py)%HeII_photoionization_rate(-i_primary,i_secondary) = phi%photo_cell_HeII
          nx_py_pz_dz(i_level_py)%photoheating_rate(-i_primary,i_secondary) = phi%heat

        enddo
      enddo	  
    enddo
  	  
  end subroutine nx_py_pz_dz_photo_rate_calculation
    
  subroutine nx_py_nz_dx_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
    real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_out
    real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_out
    type(photrates) :: phi			
    real(kind=dp) :: i,j,k,x,y,z,n
		
		
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          i = real(i_level_py)
          j = real(i_primary)
          k = real(i_secondary)
          n = real(partition(i_level_py))

          x = -(2.0*i-1.0)/2.0
          y = (2.0*j-1.0)*(2.0*i-1.0)/(4.0*n)
          z = -(2.0*k-1.0)*(2.0*i-1.0)/(4.0*n)

          HI_column_density_in = nx_py_nz_dx(-i_level_py)%HI_column_density_in(i_primary,-i_secondary)
          HeI_column_density_in = nx_py_nz_dx(-i_level_py)%HeI_column_density_in(i_primary,-i_secondary)
          HeII_column_density_in = nx_py_nz_dx(-i_level_py)%HeII_column_density_in(i_primary,-i_secondary)
          HI_column_density_out = nx_py_nz_dx(-i_level_py)%HI_column_density_out(i_primary,-i_secondary)	  
          HeI_column_density_out = nx_py_nz_dx(-i_level_py)%HeI_column_density_out(i_primary,-i_secondary)	  
          HeII_column_density_out = nx_py_nz_dx(-i_level_py)%HeII_column_density_out(i_primary,-i_secondary)  	  
          !call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
	  !	HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,correction (x,y,z))
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
		  HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,1.0_dp)
          nx_py_nz_dx(-i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) = 1.0!phi%photo_cell_HI
          nx_py_nz_dx(-i_level_py)%HeI_photoionization_rate(i_primary,-i_secondary) = phi%photo_cell_HeI
          nx_py_nz_dx(-i_level_py)%HeII_photoionization_rate(i_primary,-i_secondary) = phi%photo_cell_HeII
          nx_py_nz_dx(-i_level_py)%photoheating_rate(i_primary,-i_secondary) = phi%heat

        enddo
      enddo	  
    enddo
  	  
  end subroutine nx_py_nz_dx_photo_rate_calculation

  subroutine nx_py_nz_dy_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
    real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_out
    real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_out
    type(photrates) :: phi			
    real(kind=dp) :: i,j,k,x,y,z,n	
	
		
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          i = real(i_secondary)
          j = real(i_level_py)
          k = real(i_primary)
          n = real(partition(i_level_py))

          x = -(2.0*i-1.0)*(2.0*j-1.0)/(4.0*n)
          y = (2.0*j-1.0)/2.0
          z = -(2.0*k-1.0)*(2.0*j-1.0)/(4.0*n)	

          HI_column_density_in = nx_py_nz_dy(i_level_py)%HI_column_density_in(-i_primary,-i_secondary)
          HeI_column_density_in = nx_py_nz_dy(i_level_py)%HeI_column_density_in(-i_primary,-i_secondary)
          HeII_column_density_in = nx_py_nz_dy(i_level_py)%HeII_column_density_in(-i_primary,-i_secondary)
          HI_column_density_out = nx_py_nz_dy(i_level_py)%HI_column_density_out(-i_primary,-i_secondary)	  
          HeI_column_density_out = nx_py_nz_dy(i_level_py)%HeI_column_density_out(-i_primary,-i_secondary)	  
          HeII_column_density_out = nx_py_nz_dy(i_level_py)%HeII_column_density_out(-i_primary,-i_secondary)  	  
          !call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
	  !	HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,correction (x,y,z))
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
		  HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,1.0_dp)
          nx_py_nz_dy(i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) = 1.0!phi%photo_cell_HI
          nx_py_nz_dy(i_level_py)%HeI_photoionization_rate(-i_primary,-i_secondary) = phi%photo_cell_HeI
          nx_py_nz_dy(i_level_py)%HeII_photoionization_rate(-i_primary,-i_secondary) = phi%photo_cell_HeII
          nx_py_nz_dy(i_level_py)%photoheating_rate(-i_primary,-i_secondary) = phi%heat

        enddo
      enddo	  
    enddo
  	  
  end subroutine nx_py_nz_dy_photo_rate_calculation

  subroutine nx_py_nz_dz_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
    real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_out
    real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_out
    type(photrates) :: phi			
    real(kind=dp) :: i,j,k,x,y,z,n
		
		
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          i = real(i_primary)
          j = real(i_secondary)
          k = real(i_level_py)
          n = real(partition(i_level_py))

          x = -(2.0*i-1.0)*(2.0*k-1.0)/(4.0*n)
          y = (2.0*j-1.0)*(2.0*k-1.0)/(4.0*n)
          z = -(2.0*k-1.0)/2.0	

          HI_column_density_in = nx_py_nz_dz(-i_level_py)%HI_column_density_in(-i_primary,i_secondary)
          HeI_column_density_in = nx_py_nz_dz(-i_level_py)%HeI_column_density_in(-i_primary,i_secondary)
          HeII_column_density_in = nx_py_nz_dz(-i_level_py)%HeII_column_density_in(-i_primary,i_secondary)
          HI_column_density_out = nx_py_nz_dz(-i_level_py)%HI_column_density_out(-i_primary,i_secondary)	  
          HeI_column_density_out = nx_py_nz_dz(-i_level_py)%HeI_column_density_out(-i_primary,i_secondary)	  
          HeII_column_density_out = nx_py_nz_dz(-i_level_py)%HeII_column_density_out(-i_primary,i_secondary)  	  
          !call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
	  !	HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,correction (x,y,z))
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
		  HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,1.0_dp)
          nx_py_nz_dz(-i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) = 1.0!phi%photo_cell_HI
          nx_py_nz_dz(-i_level_py)%HeI_photoionization_rate(-i_primary,i_secondary) = phi%photo_cell_HeI
          nx_py_nz_dz(-i_level_py)%HeII_photoionization_rate(-i_primary,i_secondary) = phi%photo_cell_HeII
          nx_py_nz_dz(-i_level_py)%photoheating_rate(-i_primary,i_secondary) = phi%heat

        enddo
      enddo	  
    enddo
  	  
  end subroutine nx_py_nz_dz_photo_rate_calculation
  
  subroutine nx_ny_pz_dx_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
    real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_out
    real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_out
    type(photrates) :: phi			
    real(kind=dp) :: i,j,k,x,y,z,n
		
		
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          i = real(i_level_py)
          j = real(i_primary)
          k = real(i_secondary)
          n = real(partition(i_level_py))

          x = -(2.0*i-1.0)/2.0
          y = -(2.0*j-1.0)*(2.0*i-1.0)/(4.0*n)
          z = (2.0*k-1.0)*(2.0*i-1.0)/(4.0*n)

          HI_column_density_in = nx_ny_pz_dx(-i_level_py)%HI_column_density_in(-i_primary,i_secondary)
          HeI_column_density_in = nx_ny_pz_dx(-i_level_py)%HeI_column_density_in(-i_primary,i_secondary)
          HeII_column_density_in = nx_ny_pz_dx(-i_level_py)%HeII_column_density_in(-i_primary,i_secondary)
          HI_column_density_out = nx_ny_pz_dx(-i_level_py)%HI_column_density_out(-i_primary,i_secondary)	  
          HeI_column_density_out = nx_ny_pz_dx(-i_level_py)%HeI_column_density_out(-i_primary,i_secondary)	  
          HeII_column_density_out = nx_ny_pz_dx(-i_level_py)%HeII_column_density_out(-i_primary,i_secondary)  	  
          !call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
	  !	HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,correction (x,y,z))
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
		  HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,1.0_dp)
          nx_ny_pz_dx(-i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) = 1.0!phi%photo_cell_HI
          nx_ny_pz_dx(-i_level_py)%HeI_photoionization_rate(-i_primary,i_secondary) = phi%photo_cell_HeI
          nx_ny_pz_dx(-i_level_py)%HeII_photoionization_rate(-i_primary,i_secondary) = phi%photo_cell_HeII
          nx_ny_pz_dx(-i_level_py)%photoheating_rate(-i_primary,i_secondary) = phi%heat

        enddo
      enddo	  
    enddo
  	  
  end subroutine nx_ny_pz_dx_photo_rate_calculation

  subroutine nx_ny_pz_dy_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
    real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_out
    real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_out
    type(photrates) :: phi			
    real(kind=dp) :: i,j,k,x,y,z,n
		
		
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          i = real(i_secondary)
          j = real(i_level_py)
          k = real(i_primary)
          n = real(partition(i_level_py))

          x = -(2.0*i-1.0)*(2.0*j-1.0)/(4.0*n)
          y = -(2.0*j-1.0)/2.0
          z = (2.0*k-1.0)*(2.0*j-1.0)/(4.0*n)	

          HI_column_density_in = nx_ny_pz_dy(-i_level_py)%HI_column_density_in(i_primary,-i_secondary)
          HeI_column_density_in = nx_ny_pz_dy(-i_level_py)%HeI_column_density_in(i_primary,-i_secondary)
          HeII_column_density_in = nx_ny_pz_dy(-i_level_py)%HeII_column_density_in(i_primary,-i_secondary)
          HI_column_density_out = nx_ny_pz_dy(-i_level_py)%HI_column_density_out(i_primary,-i_secondary)	  
          HeI_column_density_out = nx_ny_pz_dy(-i_level_py)%HeI_column_density_out(i_primary,-i_secondary)	  
          HeII_column_density_out = nx_ny_pz_dy(-i_level_py)%HeII_column_density_out(i_primary,-i_secondary)  	  
          !call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
	  !	HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,correction (x,y,z))
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
		  HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,1.0_dp)
          nx_ny_pz_dy(-i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) = 1.0!phi%photo_cell_HI
          nx_ny_pz_dy(-i_level_py)%HeI_photoionization_rate(i_primary,-i_secondary) = phi%photo_cell_HeI
          nx_ny_pz_dy(-i_level_py)%HeII_photoionization_rate(i_primary,-i_secondary) = phi%photo_cell_HeII
          nx_ny_pz_dy(-i_level_py)%photoheating_rate(i_primary,-i_secondary) = phi%heat

        enddo
      enddo	  
    enddo
  	  
  end subroutine nx_ny_pz_dy_photo_rate_calculation

  subroutine nx_ny_pz_dz_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
    real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_out
    real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_out
    type(photrates) :: phi			
    real(kind=dp) :: i,j,k,x,y,z,n	

		
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          i = real(i_primary)
          j = real(i_secondary)
          k = real(i_level_py)
          n = real(partition(i_level_py))

          x = -(2.0*i-1.0)*(2.0*k-1.0)/(4.0*n)
          y = -(2.0*j-1.0)*(2.0*k-1.0)/(4.0*n)
          z = (2.0*k-1.0)/2.0	

          HI_column_density_in = nx_ny_pz_dz(i_level_py)%HI_column_density_in(-i_primary,-i_secondary)
          HeI_column_density_in = nx_ny_pz_dz(i_level_py)%HeI_column_density_in(-i_primary,-i_secondary)
          HeII_column_density_in = nx_ny_pz_dz(i_level_py)%HeII_column_density_in(-i_primary,-i_secondary)
          HI_column_density_out = nx_ny_pz_dz(i_level_py)%HI_column_density_out(-i_primary,-i_secondary)	  
          HeI_column_density_out = nx_ny_pz_dz(i_level_py)%HeI_column_density_out(-i_primary,-i_secondary)	  
          HeII_column_density_out = nx_ny_pz_dz(i_level_py)%HeII_column_density_out(-i_primary,-i_secondary)  	  
          !call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
	  !	HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,correction (x,y,z))
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
		  HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,1.0_dp)
          nx_ny_pz_dz(i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) = 1.0! phi%photo_cell_HI
          nx_ny_pz_dz(i_level_py)%HeI_photoionization_rate(-i_primary,-i_secondary) = phi%photo_cell_HeI
          nx_ny_pz_dz(i_level_py)%HeII_photoionization_rate(-i_primary,-i_secondary) = phi%photo_cell_HeII
          nx_ny_pz_dz(i_level_py)%photoheating_rate(-i_primary,-i_secondary) = phi%heat

        enddo
      enddo	  
    enddo
  	  
  end subroutine nx_ny_pz_dz_photo_rate_calculation
    
  subroutine nx_ny_nz_dx_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
    real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_out
    real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_out
    type(photrates) :: phi			
    real(kind=dp) :: i,j,k,x,y,z,n	

		
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          i = real(i_level_py)
          j = real(i_primary)
          k = real(i_secondary)
          n = real(partition(i_level_py))

          x = -(2.0*i-1.0)/2.0
          y = -(2.0*j-1.0)*(2.0*i-1.0)/(4.0*n)
          z = -(2.0*k-1.0)*(2.0*i-1.0)/(4.0*n)

          HI_column_density_in = nx_ny_nz_dx(-i_level_py)%HI_column_density_in(-i_primary,-i_secondary)
          HeI_column_density_in = nx_ny_nz_dx(-i_level_py)%HeI_column_density_in(-i_primary,-i_secondary)
          HeII_column_density_in = nx_ny_nz_dx(-i_level_py)%HeII_column_density_in(-i_primary,-i_secondary)
          HI_column_density_out = nx_ny_nz_dx(-i_level_py)%HI_column_density_out(-i_primary,-i_secondary)	  
          HeI_column_density_out = nx_ny_nz_dx(-i_level_py)%HeI_column_density_out(-i_primary,-i_secondary)	  
          HeII_column_density_out = nx_ny_nz_dx(-i_level_py)%HeII_column_density_out(-i_primary,-i_secondary)  	  
          !call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
	  !	HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,correction (x,y,z))
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
		  HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,1.0_dp)
          nx_ny_nz_dx(-i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) = 1.0!phi%photo_cell_HI
          nx_ny_nz_dx(-i_level_py)%HeI_photoionization_rate(-i_primary,-i_secondary) = phi%photo_cell_HeI
          nx_ny_nz_dx(-i_level_py)%HeII_photoionization_rate(-i_primary,-i_secondary) = phi%photo_cell_HeII
          nx_ny_nz_dx(-i_level_py)%photoheating_rate(-i_primary,-i_secondary) = phi%heat

        enddo
      enddo	  
    enddo
  	  
  end subroutine nx_ny_nz_dx_photo_rate_calculation

  subroutine nx_ny_nz_dy_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
    real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_out
    real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_out
    type(photrates) :: phi			
    real(kind=dp) :: i,j,k,x,y,z,n		
		
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          i = real(i_secondary)
          j = real(i_level_py)
          k = real(i_primary)
          n = real(partition(i_level_py))

          x = -(2.0*i-1.0)*(2.0*j-1.0)/(4.0*n)
          y = -(2.0*j-1.0)/2.0
          z = -(2.0*k-1.0)*(2.0*j-1.0)/(4.0*n)	

          HI_column_density_in = nx_ny_nz_dy(-i_level_py)%HI_column_density_in(-i_primary,-i_secondary)
          HeI_column_density_in = nx_ny_nz_dy(-i_level_py)%HeI_column_density_in(-i_primary,-i_secondary)
          HeII_column_density_in = nx_ny_nz_dy(-i_level_py)%HeII_column_density_in(-i_primary,-i_secondary)
          HI_column_density_out = nx_ny_nz_dy(-i_level_py)%HI_column_density_out(-i_primary,-i_secondary)	  
          HeI_column_density_out = nx_ny_nz_dy(-i_level_py)%HeI_column_density_out(-i_primary,-i_secondary)	  
          HeII_column_density_out = nx_ny_nz_dy(-i_level_py)%HeII_column_density_out(-i_primary,-i_secondary)  	  
          !call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
	  !	HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,correction (x,y,z))
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
		  HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,1.0_dp)
          nx_ny_nz_dy(-i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) = 1.0! phi%photo_cell_HI
          nx_ny_nz_dy(-i_level_py)%HeI_photoionization_rate(-i_primary,-i_secondary) = phi%photo_cell_HeI
          nx_ny_nz_dy(-i_level_py)%HeII_photoionization_rate(-i_primary,-i_secondary) = phi%photo_cell_HeII
          nx_ny_nz_dy(-i_level_py)%photoheating_rate(-i_primary,-i_secondary) = phi%heat

        enddo
      enddo	  
    enddo
  	  
  end subroutine nx_ny_nz_dy_photo_rate_calculation

  subroutine nx_ny_nz_dz_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
    real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_out
    real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_out
    type(photrates) :: phi			
    real(kind=dp) :: i,j,k,x,y,z,n
		
    do i_level_py = 1,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          i = real(i_primary)
          j = real(i_secondary)
          k = real(i_level_py)
          n = real(partition(i_level_py))

          x = -(2.0*i-1.0)*(2.0*k-1.0)/(4.0*n)
          y = -(2.0*j-1.0)*(2.0*k-1.0)/(4.0*n)
          z = -(2.0*k-1.0)/2.0	

          HI_column_density_in = nx_ny_nz_dz(-i_level_py)%HI_column_density_in(-i_primary,-i_secondary)
          HeI_column_density_in = nx_ny_nz_dz(-i_level_py)%HeI_column_density_in(-i_primary,-i_secondary)
          HeII_column_density_in = nx_ny_nz_dz(-i_level_py)%HeII_column_density_in(-i_primary,-i_secondary)
          HI_column_density_out = nx_ny_nz_dz(-i_level_py)%HI_column_density_out(-i_primary,-i_secondary)	  
          HeI_column_density_out = nx_ny_nz_dz(-i_level_py)%HeI_column_density_out(-i_primary,-i_secondary)	  
          HeII_column_density_out = nx_ny_nz_dz(-i_level_py)%HeII_column_density_out(-i_primary,-i_secondary)  	  
          !call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
	  !	HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,correction (x,y,z))
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,HeI_column_density_in,HeI_column_density_out,&
		  HeII_column_density_in,HeII_column_density_out,solid_angle,use_which_source,1.0_dp)
          nx_ny_nz_dz(-i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) = 1.0!phi%photo_cell_HI
          nx_ny_nz_dz(-i_level_py)%HeI_photoionization_rate(-i_primary,-i_secondary) = phi%photo_cell_HeI
          nx_ny_nz_dz(-i_level_py)%HeII_photoionization_rate(-i_primary,-i_secondary) = phi%photo_cell_HeII
          nx_ny_nz_dz(-i_level_py)%photoheating_rate(-i_primary,-i_secondary) = phi%heat

        enddo
      enddo	  
    enddo
  	  
  end subroutine nx_ny_nz_dz_photo_rate_calculation
  
  real(kind=dp) function correction (x,y,z)
    
    implicit none

    real(kind=dp), intent(in) :: x,y,z
    real(kind=dp) :: r, theta, phi
    integer :: cc = 5

    r = sqrt(x*x + y*y + z*z)
    theta = acos (z/r)
    phi = acos (x/(sqrt(x*x + y*y)))
    if (y .lt. 0.0) then
      phi = 2*pi-phi
    endif

  if (cc.eq.1) then
    ! anisotropic 1 cp pyramid_3D_HI_photoionization.dat anisotropic1_HI_photoionization.dat
    correction = (1.0 - abs(0.999* sin(6*theta)))
  else if (cc.eq.2) then
    ! anisotropic 2 cp pyramid_3D_HI_photoionization.dat anisotropic2_HI_photoionization.dat
    correction = (1.0 - abs(0.999* sin(6*phi)))
  else if (cc.eq.3) then
    ! anisotropic 3 cp pyramid_3D_HI_photoionization.dat anisotropic3_HI_photoionization.dat
    correction = (1.0 - abs(0.999* sin(6*(phi+theta))))
  else if (cc.eq.4) then
    ! anisotropic 4 cp pyramid_3D_HI_photoionization.dat anisotropic4_HI_photoionization.dat
    correction = 0.001
    if ( (theta.le.pi/12.0) .or. (theta.ge.3.0*pi/12.0 .and. theta.le.5.0*pi/12.0) .or. &
         (theta.ge.7.0*pi/12.0 .and. theta.le.9.0*pi/12.0) .or. (theta.ge.11.0*pi/12.0)) then
      if ( (phi.le.pi/12.0) .or. (phi.ge.3.0*pi/12.0 .and. phi.le.5.0*pi/12.0) .or. &
           (phi.ge.7.0*pi/12.0 .and. phi.le.9.0*pi/12.0) .or. (phi.ge.11.0*pi/12.0 .and. phi.le.13.0*pi/12.0) .or. &
           (phi.ge.15.0*pi/12.0 .and. phi.le.17.0*pi/12.0) .or. (phi.ge.19.0*pi/12.0)) then
        correction = 1.0
      endif
    endif
  else
    correction = 1.0
  endif

  end function correction


end module pyramid_photo_rate
