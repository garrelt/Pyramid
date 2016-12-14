module pyramid_photo_rate
	
    use precision, only: dp
    use input, only: level_py, partition, number_density, xHI, pi, volume_factor
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
	
  subroutine pyramid_HI_photo_rate_calculation()
		
    call px_py_pz_dx_HI_photo_rate_calculation ()
    call px_py_pz_dy_HI_photo_rate_calculation ()
    call px_py_pz_dz_HI_photo_rate_calculation ()	
    call px_py_nz_dx_HI_photo_rate_calculation ()
    call px_py_nz_dy_HI_photo_rate_calculation ()	
    call px_py_nz_dz_HI_photo_rate_calculation ()		
    call px_ny_pz_dx_HI_photo_rate_calculation ()	
    call px_ny_pz_dy_HI_photo_rate_calculation ()	
    call px_ny_pz_dz_HI_photo_rate_calculation ()	
    call px_ny_nz_dx_HI_photo_rate_calculation ()
    call px_ny_nz_dy_HI_photo_rate_calculation ()
    call px_ny_nz_dz_HI_photo_rate_calculation ()
    call nx_py_pz_dx_HI_photo_rate_calculation ()
    call nx_py_pz_dy_HI_photo_rate_calculation ()
    call nx_py_pz_dz_HI_photo_rate_calculation ()
    call nx_py_nz_dx_HI_photo_rate_calculation ()	
    call nx_py_nz_dy_HI_photo_rate_calculation ()	
    call nx_py_nz_dz_HI_photo_rate_calculation ()	
    call nx_ny_pz_dx_HI_photo_rate_calculation ()
    call nx_ny_pz_dy_HI_photo_rate_calculation ()
    call nx_ny_pz_dz_HI_photo_rate_calculation ()	
    call nx_ny_nz_dx_HI_photo_rate_calculation ()
    call nx_ny_nz_dy_HI_photo_rate_calculation ()
    call nx_ny_nz_dz_HI_photo_rate_calculation ()		
				
  end subroutine pyramid_HI_photo_rate_calculation		
		
  subroutine px_py_pz_dx_HI_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    type(photrates) :: phi	
		!volume錯撚哂 volume計岩之餘仲要用solid angle去tune番個luminosity 
    ! i_level_py = 1
    solid_angle = pyramid_grid(1)%solid_angle(1,1)
    HI_column_density_in = 0
    HI_column_density_out = px_py_pz_dx(1)%HI_column_density_out(1,1)
    call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
    px_py_pz_dx(1)%HI_photoionization_rate(1,1) = phi%photo_cell
    px_py_pz_dx(1)%HI_heating_rate(1,1) = phi%heat_cell
	
   	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          if (partition(i_level_py).eq.partition(i_level_py-1)) then	
            HI_column_density_in = px_py_pz_dx(i_level_py-1)%HI_column_density_out(i_primary,i_secondary)
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2	 
            HI_column_density_in = px_py_pz_dx(i_level_py-1)%HI_column_density_out(o_primary,o_secondary)	 
          endif	 
          HI_column_density_out = px_py_pz_dx(i_level_py)%HI_column_density_out(i_primary,i_secondary)	  
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
          px_py_pz_dx(i_level_py)%HI_photoionization_rate(i_primary,i_secondary) = phi%photo_cell
          px_py_pz_dx(i_level_py)%HI_heating_rate(i_primary,i_secondary) = phi%heat_cell

        enddo
      enddo	  
    enddo
  
  end subroutine px_py_pz_dx_HI_photo_rate_calculation

  subroutine px_py_pz_dy_HI_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    type(photrates) :: phi	
		
    ! i_level_py = 1
    solid_angle = pyramid_grid(1)%solid_angle(1,1)
    HI_column_density_in = 0
    HI_column_density_out = px_py_pz_dy(1)%HI_column_density_out(1,1)
    call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
    px_py_pz_dy(1)%HI_photoionization_rate(1,1) = phi%photo_cell
    px_py_pz_dy(1)%HI_heating_rate(1,1) = phi%heat_cell
	
   	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          if (partition(i_level_py).eq.partition(i_level_py-1)) then	
            HI_column_density_in = px_py_pz_dy(i_level_py-1)%HI_column_density_out(i_primary,i_secondary)
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2	 
            HI_column_density_in = px_py_pz_dy(i_level_py-1)%HI_column_density_out(o_primary,o_secondary)	 
          endif	 
          HI_column_density_out = px_py_pz_dy(i_level_py)%HI_column_density_out(i_primary,i_secondary)	  
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
          px_py_pz_dy(i_level_py)%HI_photoionization_rate(i_primary,i_secondary) = phi%photo_cell
          px_py_pz_dy(i_level_py)%HI_heating_rate(i_primary,i_secondary) = phi%heat_cell

        enddo
      enddo	
    enddo
  	  
  end subroutine px_py_pz_dy_HI_photo_rate_calculation

  subroutine px_py_pz_dz_HI_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    type(photrates) :: phi	
		
    ! i_level_py = 1
    solid_angle = pyramid_grid(1)%solid_angle(1,1)
    HI_column_density_in = 0
    HI_column_density_out = px_py_pz_dz(1)%HI_column_density_out(1,1)
    call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
    px_py_pz_dz(1)%HI_photoionization_rate(1,1) = phi%photo_cell
    px_py_pz_dz(1)%HI_heating_rate(1,1) = phi%heat_cell
	
   	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          if (partition(i_level_py).eq.partition(i_level_py-1)) then	
            HI_column_density_in = px_py_pz_dz(i_level_py-1)%HI_column_density_out(i_primary,i_secondary)
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2	 
            HI_column_density_in = px_py_pz_dz(i_level_py-1)%HI_column_density_out(o_primary,o_secondary)	 
          endif	 
          HI_column_density_out = px_py_pz_dz(i_level_py)%HI_column_density_out(i_primary,i_secondary)	  
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
          px_py_pz_dz(i_level_py)%HI_photoionization_rate(i_primary,i_secondary) = phi%photo_cell
          px_py_pz_dz(i_level_py)%HI_heating_rate(i_primary,i_secondary) = phi%heat_cell

        enddo
      enddo	
    enddo
  	  
  end subroutine px_py_pz_dz_HI_photo_rate_calculation
    
  subroutine px_py_nz_dx_HI_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    type(photrates) :: phi	
		
    ! i_level_py = 1
    solid_angle = pyramid_grid(1)%solid_angle(1,1)
    HI_column_density_in = 0
    HI_column_density_out = px_py_nz_dx(1)%HI_column_density_out(1,-1)
    call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
    px_py_nz_dx(1)%HI_photoionization_rate(1,-1) = phi%photo_cell
    px_py_nz_dx(1)%HI_heating_rate(1,-1) = phi%heat_cell
	
   	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          if (partition(i_level_py).eq.partition(i_level_py-1)) then	
            HI_column_density_in = px_py_nz_dx(i_level_py-1)%HI_column_density_out(i_primary,-i_secondary)
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2	 
            HI_column_density_in = px_py_nz_dx(i_level_py-1)%HI_column_density_out(o_primary,-o_secondary)	 
          endif	 
          HI_column_density_out = px_py_nz_dx(i_level_py)%HI_column_density_out(i_primary,-i_secondary)	  
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
          px_py_nz_dx(i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) = phi%photo_cell
          px_py_nz_dx(i_level_py)%HI_heating_rate(i_primary,-i_secondary) = phi%heat_cell

        enddo
      enddo	 
    enddo
  	  
  end subroutine px_py_nz_dx_HI_photo_rate_calculation

  subroutine px_py_nz_dy_HI_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    type(photrates) :: phi	
		
    ! i_level_py = 1
    solid_angle = pyramid_grid(1)%solid_angle(1,1)
    HI_column_density_in = 0
    HI_column_density_out = px_py_nz_dy(1)%HI_column_density_out(-1,1)
    call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
    px_py_nz_dy(1)%HI_photoionization_rate(-1,1) = phi%photo_cell
    px_py_nz_dy(1)%HI_heating_rate(-1,1) = phi%heat_cell
	
   	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          if (partition(i_level_py).eq.partition(i_level_py-1)) then	
            HI_column_density_in = px_py_nz_dy(i_level_py-1)%HI_column_density_out(-i_primary,i_secondary)
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2	 
            HI_column_density_in = px_py_nz_dy(i_level_py-1)%HI_column_density_out(-o_primary,o_secondary)	 
          endif	 
          HI_column_density_out = px_py_nz_dy(i_level_py)%HI_column_density_out(-i_primary,i_secondary)	  
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
          px_py_nz_dy(i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) = phi%photo_cell
          px_py_nz_dy(i_level_py)%HI_heating_rate(-i_primary,i_secondary) = phi%heat_cell

        enddo
      enddo
    enddo
  	  
  end subroutine px_py_nz_dy_HI_photo_rate_calculation

  subroutine px_py_nz_dz_HI_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    type(photrates) :: phi	
		
    ! i_level_py = 1
    solid_angle = pyramid_grid(1)%solid_angle(1,1)
    HI_column_density_in = 0
    HI_column_density_out = px_py_nz_dz(-1)%HI_column_density_out(1,1)
    call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
    px_py_nz_dz(-1)%HI_photoionization_rate(1,1) = phi%photo_cell
    px_py_nz_dz(-1)%HI_heating_rate(1,1) = phi%heat_cell
	
   	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          if (partition(i_level_py).eq.partition(i_level_py-1)) then	
            HI_column_density_in = px_py_nz_dz(-i_level_py+1)%HI_column_density_out(i_primary,i_secondary)
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2	 
            HI_column_density_in = px_py_nz_dz(-i_level_py+1)%HI_column_density_out(o_primary,o_secondary)	 
          endif	 
          HI_column_density_out = px_py_nz_dz(-i_level_py)%HI_column_density_out(i_primary,i_secondary)	  
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
          px_py_nz_dz(-i_level_py)%HI_photoionization_rate(i_primary,i_secondary) = phi%photo_cell
          px_py_nz_dz(-i_level_py)%HI_heating_rate(i_primary,i_secondary) = phi%heat_cell

        enddo
      enddo	
    enddo
  	  
  end subroutine px_py_nz_dz_HI_photo_rate_calculation
  
  subroutine px_ny_pz_dx_HI_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    type(photrates) :: phi	
		
    ! i_level_py = 1
    solid_angle = pyramid_grid(1)%solid_angle(1,1)
    HI_column_density_in = 0
    HI_column_density_out = px_ny_pz_dx(1)%HI_column_density_out(-1,1)
    call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
    px_ny_pz_dx(1)%HI_photoionization_rate(-1,1) = phi%photo_cell
    px_ny_pz_dx(1)%HI_heating_rate(-1,1) = phi%heat_cell
	
   	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          if (partition(i_level_py).eq.partition(i_level_py-1)) then	
            HI_column_density_in = px_ny_pz_dx(i_level_py-1)%HI_column_density_out(-i_primary,i_secondary)
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2	 
            HI_column_density_in = px_ny_pz_dx(i_level_py-1)%HI_column_density_out(-o_primary,o_secondary)	 
          endif	 
          HI_column_density_out = px_ny_pz_dx(i_level_py)%HI_column_density_out(-i_primary,i_secondary)	  
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
          px_ny_pz_dx(i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) = phi%photo_cell
          px_ny_pz_dx(i_level_py)%HI_heating_rate(-i_primary,i_secondary) = phi%heat_cell

        enddo
      enddo	  
    enddo
  	  
  end subroutine px_ny_pz_dx_HI_photo_rate_calculation

  subroutine px_ny_pz_dy_HI_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    type(photrates) :: phi	
		
    ! i_level_py = 1
    solid_angle = pyramid_grid(1)%solid_angle(1,1)
    HI_column_density_in = 0
    HI_column_density_out = px_ny_pz_dy(-1)%HI_column_density_out(1,1)
    call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
    px_ny_pz_dy(-1)%HI_photoionization_rate(1,1) = phi%photo_cell
    px_ny_pz_dy(-1)%HI_heating_rate(1,1) = phi%heat_cell
	
   	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          if (partition(i_level_py).eq.partition(i_level_py-1)) then	
            HI_column_density_in = px_ny_pz_dy(-i_level_py+1)%HI_column_density_out(i_primary,i_secondary)
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2	 
            HI_column_density_in = px_ny_pz_dy(-i_level_py+1)%HI_column_density_out(o_primary,o_secondary)	 
          endif	 
          HI_column_density_out = px_ny_pz_dy(-i_level_py)%HI_column_density_out(i_primary,i_secondary)	  
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
          px_ny_pz_dy(-i_level_py)%HI_photoionization_rate(i_primary,i_secondary) = phi%photo_cell
          px_ny_pz_dy(-i_level_py)%HI_heating_rate(i_primary,i_secondary) = phi%heat_cell

        enddo
      enddo	
    enddo
  	  
  end subroutine px_ny_pz_dy_HI_photo_rate_calculation

  subroutine px_ny_pz_dz_HI_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    type(photrates) :: phi	
		
    ! i_level_py = 1
    solid_angle = pyramid_grid(1)%solid_angle(1,1)
    HI_column_density_in = 0
    HI_column_density_out = px_ny_pz_dz(1)%HI_column_density_out(1,-1)
    call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
    px_ny_pz_dz(1)%HI_photoionization_rate(1,-1) = phi%photo_cell
    px_ny_pz_dz(1)%HI_heating_rate(1,-1) = phi%heat_cell
	
   	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          if (partition(i_level_py).eq.partition(i_level_py-1)) then	
            HI_column_density_in = px_ny_pz_dz(i_level_py-1)%HI_column_density_out(i_primary,-i_secondary)
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2	 
            HI_column_density_in = px_ny_pz_dz(i_level_py-1)%HI_column_density_out(o_primary,-o_secondary)	 
          endif	 
          HI_column_density_out = px_ny_pz_dz(i_level_py)%HI_column_density_out(i_primary,-i_secondary)	  
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
          px_ny_pz_dz(i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) = phi%photo_cell
          px_ny_pz_dz(i_level_py)%HI_heating_rate(i_primary,-i_secondary) = phi%heat_cell

        enddo
      enddo	
    enddo
  	  
  end subroutine px_ny_pz_dz_HI_photo_rate_calculation
    
  subroutine px_ny_nz_dx_HI_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    type(photrates) :: phi	
		
    ! i_level_py = 1
    solid_angle = pyramid_grid(1)%solid_angle(1,1)
    HI_column_density_in = 0
    HI_column_density_out = px_ny_nz_dx(1)%HI_column_density_out(-1,-1)
    call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
    px_ny_nz_dx(1)%HI_photoionization_rate(-1,-1) = phi%photo_cell
    px_ny_nz_dx(1)%HI_heating_rate(-1,-1) = phi%heat_cell

   	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          if (partition(i_level_py).eq.partition(i_level_py-1)) then	
            HI_column_density_in = px_ny_nz_dx(i_level_py-1)%HI_column_density_out(-i_primary,-i_secondary)
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2	 
            HI_column_density_in = px_ny_nz_dx(i_level_py-1)%HI_column_density_out(-o_primary,-o_secondary)	 
          endif	 
          HI_column_density_out = px_ny_nz_dx(i_level_py)%HI_column_density_out(-i_primary,-i_secondary)	 
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
          px_ny_nz_dx(i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) = phi%photo_cell
          px_ny_nz_dx(i_level_py)%HI_heating_rate(-i_primary,-i_secondary) = phi%heat_cell

        enddo
      enddo	
    enddo
  	  
  end subroutine px_ny_nz_dx_HI_photo_rate_calculation

  subroutine px_ny_nz_dy_HI_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    type(photrates) :: phi	
		
    ! i_level_py = 1
    solid_angle = pyramid_grid(1)%solid_angle(1,1)
    HI_column_density_in = 0
    HI_column_density_out = px_ny_nz_dy(-1)%HI_column_density_out(-1,1)
    call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
    px_ny_nz_dy(-1)%HI_photoionization_rate(-1,1) = phi%photo_cell
    px_ny_nz_dy(-1)%HI_heating_rate(-1,1) = phi%heat_cell
	
   	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          if (partition(i_level_py).eq.partition(i_level_py-1)) then	
            HI_column_density_in = px_ny_nz_dy(-i_level_py+1)%HI_column_density_out(-i_primary,i_secondary)
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2	 
            HI_column_density_in = px_ny_nz_dy(-i_level_py+1)%HI_column_density_out(-o_primary,o_secondary)	 
          endif	 
          HI_column_density_out = px_ny_nz_dy(-i_level_py)%HI_column_density_out(-i_primary,i_secondary)	  
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
          px_ny_nz_dy(-i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) = phi%photo_cell
          px_ny_nz_dy(-i_level_py)%HI_heating_rate(-i_primary,i_secondary) = phi%heat_cell

        enddo
      enddo	 
    enddo
  	  
  end subroutine px_ny_nz_dy_HI_photo_rate_calculation

  subroutine px_ny_nz_dz_HI_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    type(photrates) :: phi	
		
    ! i_level_py = 1
    solid_angle = pyramid_grid(1)%solid_angle(1,1)
    HI_column_density_in = 0
    HI_column_density_out = px_ny_nz_dz(-1)%HI_column_density_out(1,-1)
    call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
    px_ny_nz_dz(-1)%HI_photoionization_rate(1,-1) = phi%photo_cell
    px_ny_nz_dz(-1)%HI_heating_rate(1,-1) = phi%heat_cell
	
   	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          if (partition(i_level_py).eq.partition(i_level_py-1)) then	
            HI_column_density_in = px_ny_nz_dz(-i_level_py+1)%HI_column_density_out(i_primary,-i_secondary)
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2	 
            HI_column_density_in = px_ny_nz_dz(-i_level_py+1)%HI_column_density_out(o_primary,-o_secondary)	 
          endif	 
          HI_column_density_out = px_ny_nz_dz(-i_level_py)%HI_column_density_out(i_primary,-i_secondary)	  
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
          px_ny_nz_dz(-i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) = phi%photo_cell
          px_ny_nz_dz(-i_level_py)%HI_heating_rate(i_primary,-i_secondary) = phi%heat_cell

        enddo
      enddo	
    enddo
  	  
  end subroutine px_ny_nz_dz_HI_photo_rate_calculation
  
  subroutine nx_py_pz_dx_HI_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    type(photrates) :: phi	
		
    ! i_level_py = 1
    solid_angle = pyramid_grid(1)%solid_angle(1,1)
    HI_column_density_in = 0
    HI_column_density_out = nx_py_pz_dx(-1)%HI_column_density_out(1,1)
	!write(*,*)HI_column_density_out ! it is zero!!!!
    call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
    nx_py_pz_dx(-1)%HI_photoionization_rate(1,1) = phi%photo_cell
    nx_py_pz_dx(-1)%HI_heating_rate(1,1) = phi%heat_cell
	!write(*,*)nx_py_pz_dx(-1)%HI_photoionization_rate(1,1)
   	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          if (partition(i_level_py).eq.partition(i_level_py-1)) then	
            HI_column_density_in = nx_py_pz_dx(-i_level_py+1)%HI_column_density_out(i_primary,i_secondary)
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2	 
            HI_column_density_in = nx_py_pz_dx(-i_level_py+1)%HI_column_density_out(o_primary,o_secondary)	 
          endif	 
          HI_column_density_out = nx_py_pz_dx(-i_level_py)%HI_column_density_out(i_primary,i_secondary)	  
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
          nx_py_pz_dx(-i_level_py)%HI_photoionization_rate(i_primary,i_secondary) = phi%photo_cell
          nx_py_pz_dx(-i_level_py)%HI_heating_rate(i_primary,i_secondary) = phi%heat_cell

        enddo
      enddo 	  	  
    enddo
  	  
  end subroutine nx_py_pz_dx_HI_photo_rate_calculation

  subroutine nx_py_pz_dy_HI_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    type(photrates) :: phi	
		
    ! i_level_py = 1
    solid_angle = pyramid_grid(1)%solid_angle(1,1)
    HI_column_density_in = 0
    HI_column_density_out = nx_py_pz_dy(1)%HI_column_density_out(1,-1)
    call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
    nx_py_pz_dy(1)%HI_photoionization_rate(1,-1) = phi%photo_cell
    nx_py_pz_dy(1)%HI_heating_rate(1,-1) = phi%heat_cell
	
   	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          if (partition(i_level_py).eq.partition(i_level_py-1)) then	
            HI_column_density_in = nx_py_pz_dy(i_level_py-1)%HI_column_density_out(i_primary,-i_secondary)
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2	 
            HI_column_density_in = nx_py_pz_dy(i_level_py-1)%HI_column_density_out(o_primary,-o_secondary)	 
          endif	 
          HI_column_density_out = nx_py_pz_dy(i_level_py)%HI_column_density_out(i_primary,-i_secondary)	  
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
          nx_py_pz_dy(i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) = phi%photo_cell
          nx_py_pz_dy(i_level_py)%HI_heating_rate(i_primary,-i_secondary) = phi%heat_cell

        enddo
      enddo	
    enddo
  	  
  end subroutine nx_py_pz_dy_HI_photo_rate_calculation

  subroutine nx_py_pz_dz_HI_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    type(photrates) :: phi	
		
    ! i_level_py = 1
    solid_angle = pyramid_grid(1)%solid_angle(1,1)
    HI_column_density_in = 0
    HI_column_density_out = nx_py_pz_dz(1)%HI_column_density_out(-1,1)
    call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
    nx_py_pz_dz(1)%HI_photoionization_rate(-1,1) = phi%photo_cell
    nx_py_pz_dz(1)%HI_heating_rate(-1,1) = phi%heat_cell
	
   	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          if (partition(i_level_py).eq.partition(i_level_py-1)) then	
            HI_column_density_in = nx_py_pz_dz(i_level_py-1)%HI_column_density_out(-i_primary,i_secondary)
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2	 
            HI_column_density_in = nx_py_pz_dz(i_level_py-1)%HI_column_density_out(-o_primary,o_secondary)	 
          endif	 
          HI_column_density_out = nx_py_pz_dz(i_level_py)%HI_column_density_out(-i_primary,i_secondary)	  
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
          nx_py_pz_dz(i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) = phi%photo_cell
          nx_py_pz_dz(i_level_py)%HI_heating_rate(-i_primary,i_secondary) = phi%heat_cell

        enddo
      enddo	
    enddo
  	  
  end subroutine nx_py_pz_dz_HI_photo_rate_calculation
    
  subroutine nx_py_nz_dx_HI_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    type(photrates) :: phi	
		
    ! i_level_py = 1
    solid_angle = pyramid_grid(1)%solid_angle(1,1)
    HI_column_density_in = 0
    HI_column_density_out = nx_py_nz_dx(-1)%HI_column_density_out(1,-1)
    call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
    nx_py_nz_dx(-1)%HI_photoionization_rate(1,-1) = phi%photo_cell
    nx_py_nz_dx(-1)%HI_heating_rate(1,-1) = phi%heat_cell
	
   	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          if (partition(i_level_py).eq.partition(i_level_py-1)) then	
            HI_column_density_in = nx_py_nz_dx(-i_level_py+1)%HI_column_density_out(i_primary,-i_secondary)
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2	 
            HI_column_density_in = nx_py_nz_dx(-i_level_py+1)%HI_column_density_out(o_primary,-o_secondary)	 
          endif	 
          HI_column_density_out = nx_py_nz_dx(-i_level_py)%HI_column_density_out(i_primary,-i_secondary)	  
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
          nx_py_nz_dx(-i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) = phi%photo_cell
          nx_py_nz_dx(-i_level_py)%HI_heating_rate(i_primary,-i_secondary) = phi%heat_cell

        enddo
      enddo	 
    enddo
  	  
  end subroutine nx_py_nz_dx_HI_photo_rate_calculation

  subroutine nx_py_nz_dy_HI_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    type(photrates) :: phi	
		
    ! i_level_py = 1
    solid_angle = pyramid_grid(1)%solid_angle(1,1)
    HI_column_density_in = 0
    HI_column_density_out = nx_py_nz_dy(1)%HI_column_density_out(-1,-1)
    call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
    nx_py_nz_dy(1)%HI_photoionization_rate(-1,-1) = phi%photo_cell
    nx_py_nz_dy(1)%HI_heating_rate(-1,-1) = phi%heat_cell
	
   	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          if (partition(i_level_py).eq.partition(i_level_py-1)) then	
            HI_column_density_in = nx_py_nz_dy(i_level_py-1)%HI_column_density_out(-i_primary,-i_secondary)
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2	 
            HI_column_density_in = nx_py_nz_dy(i_level_py-1)%HI_column_density_out(-o_primary,-o_secondary)	 
          endif	 
          HI_column_density_out = nx_py_nz_dy(i_level_py)%HI_column_density_out(-i_primary,-i_secondary)	  
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
          nx_py_nz_dy(i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) = phi%photo_cell
          nx_py_nz_dy(i_level_py)%HI_heating_rate(-i_primary,-i_secondary) = phi%heat_cell

        enddo
      enddo	
    enddo
  	  
  end subroutine nx_py_nz_dy_HI_photo_rate_calculation

  subroutine nx_py_nz_dz_HI_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    type(photrates) :: phi	
		
    ! i_level_py = 1
    solid_angle = pyramid_grid(1)%solid_angle(1,1)
    HI_column_density_in = 0
    HI_column_density_out = nx_py_nz_dz(-1)%HI_column_density_out(-1,1)
    call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
    nx_py_nz_dz(-1)%HI_photoionization_rate(-1,1) = phi%photo_cell
    nx_py_nz_dz(-1)%HI_heating_rate(-1,1) = phi%heat_cell
	
   	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          if (partition(i_level_py).eq.partition(i_level_py-1)) then	
            HI_column_density_in = nx_py_nz_dz(-i_level_py+1)%HI_column_density_out(-i_primary,i_secondary)
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2	 
            HI_column_density_in = nx_py_nz_dz(-i_level_py+1)%HI_column_density_out(-o_primary,o_secondary)	 
          endif	 
          HI_column_density_out = nx_py_nz_dz(-i_level_py)%HI_column_density_out(-i_primary,i_secondary)	  
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
          nx_py_nz_dz(-i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) = phi%photo_cell
          nx_py_nz_dz(-i_level_py)%HI_heating_rate(-i_primary,i_secondary) = phi%heat_cell

        enddo
      enddo	
    enddo
  	  
  end subroutine nx_py_nz_dz_HI_photo_rate_calculation
  
  subroutine nx_ny_pz_dx_HI_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    type(photrates) :: phi	
		
    ! i_level_py = 1
    solid_angle = pyramid_grid(1)%solid_angle(1,1)
    HI_column_density_in = 0
    HI_column_density_out = nx_ny_pz_dx(-1)%HI_column_density_out(-1,1)
    call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
    nx_ny_pz_dx(-1)%HI_photoionization_rate(-1,1) = phi%photo_cell
    nx_ny_pz_dx(-1)%HI_heating_rate(-1,1) = phi%heat_cell
	
   	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          if (partition(i_level_py).eq.partition(i_level_py-1)) then	
            HI_column_density_in = nx_ny_pz_dx(-i_level_py+1)%HI_column_density_out(-i_primary,i_secondary)
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2	 
            HI_column_density_in = nx_ny_pz_dx(-i_level_py+1)%HI_column_density_out(-o_primary,o_secondary)	 
          endif	 
          HI_column_density_out = nx_ny_pz_dx(-i_level_py)%HI_column_density_out(-i_primary,i_secondary)	  
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
          nx_ny_pz_dx(-i_level_py)%HI_photoionization_rate(-i_primary,i_secondary) = phi%photo_cell
          nx_ny_pz_dx(-i_level_py)%HI_heating_rate(-i_primary,i_secondary) = phi%heat_cell

        enddo
      enddo	  
    enddo
  	  
  end subroutine nx_ny_pz_dx_HI_photo_rate_calculation

  subroutine nx_ny_pz_dy_HI_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    type(photrates) :: phi	
		
    ! i_level_py = 1
    solid_angle = pyramid_grid(1)%solid_angle(1,1)
    HI_column_density_in = 0
    HI_column_density_out = nx_ny_pz_dy(-1)%HI_column_density_out(1,-1)
    call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
    nx_ny_pz_dy(-1)%HI_photoionization_rate(1,-1) = phi%photo_cell
    nx_ny_pz_dy(-1)%HI_heating_rate(1,-1) = phi%heat_cell
	
   	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          if (partition(i_level_py).eq.partition(i_level_py-1)) then	
            HI_column_density_in = nx_ny_pz_dy(-i_level_py+1)%HI_column_density_out(i_primary,-i_secondary)
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2	 
            HI_column_density_in = nx_ny_pz_dy(-i_level_py+1)%HI_column_density_out(o_primary,-o_secondary)	 
          endif	 
          HI_column_density_out = nx_ny_pz_dy(-i_level_py)%HI_column_density_out(i_primary,-i_secondary)	  
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
          nx_ny_pz_dy(-i_level_py)%HI_photoionization_rate(i_primary,-i_secondary) = phi%photo_cell
          nx_ny_pz_dy(-i_level_py)%HI_heating_rate(i_primary,-i_secondary) = phi%heat_cell

        enddo
      enddo	
    enddo
  	  
  end subroutine nx_ny_pz_dy_HI_photo_rate_calculation

  subroutine nx_ny_pz_dz_HI_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    type(photrates) :: phi	
		
    ! i_level_py = 1
    solid_angle = pyramid_grid(1)%solid_angle(1,1)
    HI_column_density_in = 0
    HI_column_density_out = nx_ny_pz_dz(1)%HI_column_density_out(-1,-1)
    call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
    nx_ny_pz_dz(1)%HI_photoionization_rate(-1,-1) = phi%photo_cell
    nx_ny_pz_dz(1)%HI_heating_rate(-1,-1) = phi%heat_cell
	
   	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          if (partition(i_level_py).eq.partition(i_level_py-1)) then	
            HI_column_density_in = nx_ny_pz_dz(i_level_py-1)%HI_column_density_out(-i_primary,-i_secondary)
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2	 
            HI_column_density_in = nx_ny_pz_dz(i_level_py-1)%HI_column_density_out(-o_primary,-o_secondary)	 
          endif	 
          HI_column_density_out = nx_ny_pz_dz(i_level_py)%HI_column_density_out(-i_primary,-i_secondary)	  
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
          nx_ny_pz_dz(i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) = phi%photo_cell
          nx_ny_pz_dz(i_level_py)%HI_heating_rate(-i_primary,-i_secondary) = phi%heat_cell

        enddo
      enddo	
    enddo
  	  
  end subroutine nx_ny_pz_dz_HI_photo_rate_calculation
    
  subroutine nx_ny_nz_dx_HI_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    type(photrates) :: phi	
		
    ! i_level_py = 1
    solid_angle = pyramid_grid(1)%solid_angle(1,1)
    HI_column_density_in = 0
    HI_column_density_out = nx_ny_nz_dx(-1)%HI_column_density_out(-1,-1)
    call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
    nx_ny_nz_dx(-1)%HI_photoionization_rate(-1,-1) = phi%photo_cell
    nx_ny_nz_dx(-1)%HI_heating_rate(-1,-1) = phi%heat_cell
	
   	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          if (partition(i_level_py).eq.partition(i_level_py-1)) then	
            HI_column_density_in = nx_ny_nz_dx(-i_level_py+1)%HI_column_density_out(-i_primary,-i_secondary)
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2	 
            HI_column_density_in = nx_ny_nz_dx(-i_level_py+1)%HI_column_density_out(-o_primary,-o_secondary)	 
          endif	 
          HI_column_density_out = nx_ny_nz_dx(-i_level_py)%HI_column_density_out(-i_primary,-i_secondary)	  
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
          nx_ny_nz_dx(-i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) = phi%photo_cell
          nx_ny_nz_dx(-i_level_py)%HI_heating_rate(-i_primary,-i_secondary) = phi%heat_cell

        enddo
      enddo	 
    enddo
  	  
  end subroutine nx_ny_nz_dx_HI_photo_rate_calculation

  subroutine nx_ny_nz_dy_HI_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    type(photrates) :: phi	
		
    ! i_level_py = 1
    solid_angle = pyramid_grid(1)%solid_angle(1,1)
    HI_column_density_in = 0
    HI_column_density_out = nx_ny_nz_dy(-1)%HI_column_density_out(-1,-1)
    call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
    nx_ny_nz_dy(-1)%HI_photoionization_rate(-1,-1) = phi%photo_cell
    nx_ny_nz_dy(-1)%HI_heating_rate(-1,-1) = phi%heat_cell
	
   	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          if (partition(i_level_py).eq.partition(i_level_py-1)) then	
            HI_column_density_in = nx_ny_nz_dy(-i_level_py+1)%HI_column_density_out(-i_primary,-i_secondary)
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2	 
            HI_column_density_in = nx_ny_nz_dy(-i_level_py+1)%HI_column_density_out(-o_primary,-o_secondary)	 
          endif	 
          HI_column_density_out = nx_ny_nz_dy(-i_level_py)%HI_column_density_out(-i_primary,-i_secondary)	  
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
          nx_ny_nz_dy(-i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) = phi%photo_cell
          nx_ny_nz_dy(-i_level_py)%HI_heating_rate(-i_primary,-i_secondary) = phi%heat_cell

        enddo
      enddo	
    enddo
  	  
  end subroutine nx_ny_nz_dy_HI_photo_rate_calculation

  subroutine nx_ny_nz_dz_HI_photo_rate_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    real(kind=dp) :: solid_angle
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_out
    type(photrates) :: phi	
		
    ! i_level_py = 1
    solid_angle = pyramid_grid(1)%solid_angle(1,1)
    HI_column_density_in = 0
    HI_column_density_out = nx_ny_nz_dz(-1)%HI_column_density_out(-1,-1)
    call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
    nx_ny_nz_dz(-1)%HI_photoionization_rate(-1,-1) = phi%photo_cell
    nx_ny_nz_dz(-1)%HI_heating_rate(-1,-1) = phi%heat_cell
	
   	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)

          solid_angle = pyramid_grid(i_level_py)%solid_angle(i_primary,i_secondary)
          if (partition(i_level_py).eq.partition(i_level_py-1)) then	
            HI_column_density_in = nx_ny_nz_dz(-i_level_py+1)%HI_column_density_out(-i_primary,-i_secondary)
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2	 
            HI_column_density_in = nx_ny_nz_dz(-i_level_py+1)%HI_column_density_out(-o_primary,-o_secondary)	 
          endif	 
          HI_column_density_out = nx_ny_nz_dz(-i_level_py)%HI_column_density_out(-i_primary,-i_secondary)	  
          call photoion_solid_angle(phi,HI_column_density_in,HI_column_density_out,solid_angle)
          nx_ny_nz_dz(-i_level_py)%HI_photoionization_rate(-i_primary,-i_secondary) = phi%photo_cell
          nx_ny_nz_dz(-i_level_py)%HI_heating_rate(-i_primary,-i_secondary) = phi%heat_cell

        enddo
      enddo	  
    enddo
  	  
  end subroutine nx_ny_nz_dz_HI_photo_rate_calculation
  
end module pyramid_photo_rate