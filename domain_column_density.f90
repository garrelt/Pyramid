module pyramid_column_density
	
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
  	     		     nx_ny_nz_dx, nx_ny_nz_dy, nx_ny_nz_dz		   	   
		
contains
	
  subroutine pyramid_HI_column_density_calculation()
		
    call px_py_pz_dx_HI_column_density_calculation ()
    call px_py_pz_dy_HI_column_density_calculation ()
    call px_py_pz_dz_HI_column_density_calculation ()	
    call px_py_nz_dx_HI_column_density_calculation ()
    call px_py_nz_dy_HI_column_density_calculation ()	
    call px_py_nz_dz_HI_column_density_calculation ()		
    call px_ny_pz_dx_HI_column_density_calculation ()	
    call px_ny_pz_dy_HI_column_density_calculation ()	
    call px_ny_pz_dz_HI_column_density_calculation ()	
    call px_ny_nz_dx_HI_column_density_calculation ()
    call px_ny_nz_dy_HI_column_density_calculation ()
    call px_ny_nz_dz_HI_column_density_calculation ()
    call nx_py_pz_dx_HI_column_density_calculation ()
    call nx_py_pz_dy_HI_column_density_calculation ()
    call nx_py_pz_dz_HI_column_density_calculation ()
    call nx_py_nz_dx_HI_column_density_calculation ()	
    call nx_py_nz_dy_HI_column_density_calculation ()	
    call nx_py_nz_dz_HI_column_density_calculation ()	
    call nx_ny_pz_dx_HI_column_density_calculation ()
    call nx_ny_pz_dy_HI_column_density_calculation ()
    call nx_ny_pz_dz_HI_column_density_calculation ()	
    call nx_ny_nz_dx_HI_column_density_calculation ()
    call nx_ny_nz_dy_HI_column_density_calculation ()
    call nx_ny_nz_dz_HI_column_density_calculation ()		
				
  end subroutine pyramid_HI_column_density_calculation		
		
  subroutine px_py_pz_dx_HI_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	
	! i_level_py = 1
	HI_column_density_in = 0.0
    HI_column_density_cell = px_py_pz_dx(1)%HI_density(1,1) * pyramid_grid(1)%cellsize(1,1)
	px_py_pz_dx(1)%HI_column_density_out(1,1) = HI_column_density_in + HI_column_density_cell
	
	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)	
   	      if (partition(i_level_py).eq.partition(i_level_py-1)) then
		  	HI_column_density_in = px_py_pz_dx(i_level_py-1)%HI_column_density_out(i_primary,i_secondary)
            HI_column_density_cell = px_py_pz_dx(i_level_py)%HI_density(i_primary,i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	px_py_pz_dx(i_level_py)%HI_column_density_out(i_primary,i_secondary) = HI_column_density_in + &
                                                                                HI_column_density_cell 		
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2
		  	HI_column_density_in = px_py_pz_dx(i_level_py-1)%HI_column_density_out(o_primary,o_secondary)
            HI_column_density_cell = px_py_pz_dx(i_level_py)%HI_density(i_primary,i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	px_py_pz_dx(i_level_py)%HI_column_density_out(i_primary,i_secondary) = HI_column_density_in + &
		  			  			  			  			  			  	        HI_column_density_cell			
          endif				  
		enddo	  
      enddo	  
    enddo	  

  end subroutine px_py_pz_dx_HI_column_density_calculation
  
  subroutine px_py_pz_dy_HI_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	
	! i_level_py = 1
	HI_column_density_in = 0.0
    HI_column_density_cell = px_py_pz_dy(1)%HI_density(1,1) * pyramid_grid(1)%cellsize(1,1)
	px_py_pz_dy(1)%HI_column_density_out(1,1) = HI_column_density_in + HI_column_density_cell
	
	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)	
   	      if (partition(i_level_py).eq.partition(i_level_py-1)) then
		  	HI_column_density_in = px_py_pz_dy(i_level_py-1)%HI_column_density_out(i_primary,i_secondary)
            HI_column_density_cell = px_py_pz_dy(i_level_py)%HI_density(i_primary,i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	px_py_pz_dy(i_level_py)%HI_column_density_out(i_primary,i_secondary) = HI_column_density_in + &
                                                                                HI_column_density_cell 		
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2
		  	HI_column_density_in = px_py_pz_dy(i_level_py-1)%HI_column_density_out(o_primary,o_secondary)
            HI_column_density_cell = px_py_pz_dy(i_level_py)%HI_density(i_primary,i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	px_py_pz_dy(i_level_py)%HI_column_density_out(i_primary,i_secondary) = HI_column_density_in + &
		  			  			  			  			  			  	        HI_column_density_cell			
          endif				  
		enddo	  
      enddo	  
    enddo	  
	  
  end subroutine px_py_pz_dy_HI_column_density_calculation
  
  subroutine px_py_pz_dz_HI_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	
	! i_level_py = 1
	HI_column_density_in = 0.0
    HI_column_density_cell = px_py_pz_dz(1)%HI_density(1,1) * pyramid_grid(1)%cellsize(1,1)
	px_py_pz_dz(1)%HI_column_density_out(1,1) = HI_column_density_in + HI_column_density_cell
	
	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)	
   	      if (partition(i_level_py).eq.partition(i_level_py-1)) then
		  	HI_column_density_in = px_py_pz_dz(i_level_py-1)%HI_column_density_out(i_primary,i_secondary)
            HI_column_density_cell = px_py_pz_dz(i_level_py)%HI_density(i_primary,i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	px_py_pz_dz(i_level_py)%HI_column_density_out(i_primary,i_secondary) = HI_column_density_in + &
                                                                                HI_column_density_cell 		
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2
		  	HI_column_density_in = px_py_pz_dz(i_level_py-1)%HI_column_density_out(o_primary,o_secondary)
            HI_column_density_cell = px_py_pz_dz(i_level_py)%HI_density(i_primary,i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	px_py_pz_dz(i_level_py)%HI_column_density_out(i_primary,i_secondary) = HI_column_density_in + &
		  			  			  			  			  			  	        HI_column_density_cell			
          endif				  
		enddo	  
      enddo	  
    enddo	  
	  
  end subroutine px_py_pz_dz_HI_column_density_calculation
     
  subroutine px_py_nz_dx_HI_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	
	! i_level_py = 1
	HI_column_density_in = 0.0
    HI_column_density_cell = px_py_nz_dx(1)%HI_density(1,-1) * pyramid_grid(1)%cellsize(1,1)
	px_py_nz_dx(1)%HI_column_density_out(1,-1) = HI_column_density_in + HI_column_density_cell
	
	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)	
   	      if (partition(i_level_py).eq.partition(i_level_py-1)) then
		  	HI_column_density_in = px_py_nz_dx(i_level_py-1)%HI_column_density_out(i_primary,-i_secondary)
            HI_column_density_cell = px_py_nz_dx(i_level_py)%HI_density(i_primary,-i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	px_py_nz_dx(i_level_py)%HI_column_density_out(i_primary,-i_secondary) = HI_column_density_in + &
                                                                                HI_column_density_cell 		
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2
		  	HI_column_density_in = px_py_nz_dx(i_level_py-1)%HI_column_density_out(o_primary,-o_secondary)
            HI_column_density_cell = px_py_nz_dx(i_level_py)%HI_density(i_primary,-i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	px_py_nz_dx(i_level_py)%HI_column_density_out(i_primary,-i_secondary) = HI_column_density_in + &
		  			  			  			  			  			  	        HI_column_density_cell			
          endif				  
		enddo	  
      enddo	  
    enddo	  
	  
  end subroutine px_py_nz_dx_HI_column_density_calculation
  
  subroutine px_py_nz_dy_HI_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	
	! i_level_py = 1
	HI_column_density_in = 0.0
    HI_column_density_cell = px_py_nz_dy(1)%HI_density(-1,1) * pyramid_grid(1)%cellsize(1,1)
	px_py_nz_dy(1)%HI_column_density_out(-1,1) = HI_column_density_in + HI_column_density_cell
	
	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)	
   	      if (partition(i_level_py).eq.partition(i_level_py-1)) then
		  	HI_column_density_in = px_py_nz_dy(i_level_py-1)%HI_column_density_out(-i_primary,i_secondary)
            HI_column_density_cell = px_py_nz_dy(i_level_py)%HI_density(-i_primary,i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	px_py_nz_dy(i_level_py)%HI_column_density_out(-i_primary,i_secondary) = HI_column_density_in + &
                                                                                HI_column_density_cell 		
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2
		  	HI_column_density_in = px_py_nz_dy(i_level_py-1)%HI_column_density_out(-o_primary,o_secondary)
            HI_column_density_cell = px_py_nz_dy(i_level_py)%HI_density(-i_primary,i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	px_py_nz_dy(i_level_py)%HI_column_density_out(-i_primary,i_secondary) = HI_column_density_in + &
		  			  			  			  			  			  	        HI_column_density_cell			
          endif				  
		enddo	  
      enddo	  
    enddo	  
	  
  end subroutine px_py_nz_dy_HI_column_density_calculation
  
  subroutine px_py_nz_dz_HI_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	
	! i_level_py = 1
	HI_column_density_in = 0.0
    HI_column_density_cell = px_py_nz_dz(-1)%HI_density(1,1) * pyramid_grid(1)%cellsize(1,1)
	px_py_nz_dz(-1)%HI_column_density_out(1,1) = HI_column_density_in + HI_column_density_cell
	
	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)	
   	      if (partition(i_level_py).eq.partition(i_level_py-1)) then
		  	HI_column_density_in = px_py_nz_dz(-i_level_py+1)%HI_column_density_out(i_primary,i_secondary)
            HI_column_density_cell = px_py_nz_dz(-i_level_py)%HI_density(i_primary,i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	px_py_nz_dz(-i_level_py)%HI_column_density_out(i_primary,i_secondary) = HI_column_density_in + &
                                                                                HI_column_density_cell 		
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2
		  	HI_column_density_in = px_py_nz_dz(-i_level_py+1)%HI_column_density_out(o_primary,o_secondary)
            HI_column_density_cell = px_py_nz_dz(-i_level_py)%HI_density(i_primary,i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	px_py_nz_dz(-i_level_py)%HI_column_density_out(i_primary,i_secondary) = HI_column_density_in + &
		  			  			  			  			  			  	        HI_column_density_cell			
          endif				  
		enddo	  
      enddo	  
    enddo	  
	  
  end subroutine px_py_nz_dz_HI_column_density_calculation
  	 
  subroutine px_ny_pz_dx_HI_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	
	! i_level_py = 1
	HI_column_density_in = 0.0
    HI_column_density_cell = px_ny_pz_dx(1)%HI_density(-1,1) * pyramid_grid(1)%cellsize(1,1)
	px_ny_pz_dx(1)%HI_column_density_out(-1,1) = HI_column_density_in + HI_column_density_cell
	
	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)	
   	      if (partition(i_level_py).eq.partition(i_level_py-1)) then
		  	HI_column_density_in = px_ny_pz_dx(i_level_py-1)%HI_column_density_out(-i_primary,i_secondary)
            HI_column_density_cell = px_ny_pz_dx(i_level_py)%HI_density(-i_primary,i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	px_ny_pz_dx(i_level_py)%HI_column_density_out(-i_primary,i_secondary) = HI_column_density_in + &
                                                                                HI_column_density_cell 		
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2
		  	HI_column_density_in = px_ny_pz_dx(i_level_py-1)%HI_column_density_out(-o_primary,o_secondary)
            HI_column_density_cell = px_ny_pz_dx(i_level_py)%HI_density(-i_primary,i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	px_ny_pz_dx(i_level_py)%HI_column_density_out(-i_primary,i_secondary) = HI_column_density_in + &
		  			  			  			  			  			  	        HI_column_density_cell			
          endif				  
		enddo	  
      enddo	  
    enddo	  
	  
  end subroutine px_ny_pz_dx_HI_column_density_calculation
  
  subroutine px_ny_pz_dy_HI_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	
	! i_level_py = 1
	HI_column_density_in = 0.0
    HI_column_density_cell = px_ny_pz_dy(-1)%HI_density(1,1) * pyramid_grid(1)%cellsize(1,1)
	px_ny_pz_dy(-1)%HI_column_density_out(1,1) = HI_column_density_in + HI_column_density_cell
	
	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)	
   	      if (partition(i_level_py).eq.partition(i_level_py-1)) then
		  	HI_column_density_in = px_ny_pz_dy(-i_level_py+1)%HI_column_density_out(i_primary,i_secondary)
            HI_column_density_cell = px_ny_pz_dy(-i_level_py)%HI_density(i_primary,i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	px_ny_pz_dy(-i_level_py)%HI_column_density_out(i_primary,i_secondary) = HI_column_density_in + &
                                                                                HI_column_density_cell 		
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2
		  	HI_column_density_in = px_ny_pz_dy(-i_level_py+1)%HI_column_density_out(o_primary,o_secondary)
            HI_column_density_cell = px_ny_pz_dy(-i_level_py)%HI_density(i_primary,i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	px_ny_pz_dy(-i_level_py)%HI_column_density_out(i_primary,i_secondary) = HI_column_density_in + &
		  			  			  			  			  			  	        HI_column_density_cell			
          endif				  
		enddo	  
      enddo	  
    enddo	  
	  
  end subroutine px_ny_pz_dy_HI_column_density_calculation
  
  subroutine px_ny_pz_dz_HI_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	
	! i_level_py = 1
	HI_column_density_in = 0.0
    HI_column_density_cell = px_ny_pz_dz(1)%HI_density(1,-1) * pyramid_grid(1)%cellsize(1,1)
	px_ny_pz_dz(1)%HI_column_density_out(1,-1) = HI_column_density_in + HI_column_density_cell
	
	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)	
   	      if (partition(i_level_py).eq.partition(i_level_py-1)) then
		  	HI_column_density_in = px_ny_pz_dz(i_level_py-1)%HI_column_density_out(i_primary,-i_secondary)
            HI_column_density_cell = px_ny_pz_dz(i_level_py)%HI_density(i_primary,-i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	px_ny_pz_dz(i_level_py)%HI_column_density_out(i_primary,-i_secondary) = HI_column_density_in + &
                                                                                HI_column_density_cell 		
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2
		  	HI_column_density_in = px_ny_pz_dz(i_level_py-1)%HI_column_density_out(o_primary,-o_secondary)
            HI_column_density_cell = px_ny_pz_dz(i_level_py)%HI_density(i_primary,-i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	px_ny_pz_dz(i_level_py)%HI_column_density_out(i_primary,-i_secondary) = HI_column_density_in + &
		  			  			  			  			  			  	        HI_column_density_cell			
          endif				  
		enddo	  
      enddo	  
    enddo	  
	  
  end subroutine px_ny_pz_dz_HI_column_density_calculation
     
  subroutine px_ny_nz_dx_HI_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	
	! i_level_py = 1
	HI_column_density_in = 0.0
    HI_column_density_cell = px_ny_nz_dx(1)%HI_density(-1,-1) * pyramid_grid(1)%cellsize(1,1)
	px_ny_nz_dx(1)%HI_column_density_out(-1,-1) = HI_column_density_in + HI_column_density_cell
	
	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)	
   	      if (partition(i_level_py).eq.partition(i_level_py-1)) then
		  	HI_column_density_in = px_ny_nz_dx(i_level_py-1)%HI_column_density_out(-i_primary,-i_secondary)
            HI_column_density_cell = px_ny_nz_dx(i_level_py)%HI_density(-i_primary,-i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	px_ny_nz_dx(i_level_py)%HI_column_density_out(-i_primary,-i_secondary) = HI_column_density_in + &
                                                                                HI_column_density_cell 		
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2
		  	HI_column_density_in = px_ny_nz_dx(i_level_py-1)%HI_column_density_out(-o_primary,-o_secondary)
            HI_column_density_cell = px_ny_nz_dx(i_level_py)%HI_density(-i_primary,-i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	px_ny_nz_dx(i_level_py)%HI_column_density_out(-i_primary,-i_secondary) = HI_column_density_in + &
		  			  			  			  			  			  	        HI_column_density_cell			
          endif				  
		enddo	  
      enddo	  
    enddo	  
	  
  end subroutine px_ny_nz_dx_HI_column_density_calculation
  
  subroutine px_ny_nz_dy_HI_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	
	! i_level_py = 1
	HI_column_density_in = 0.0
    HI_column_density_cell = px_ny_nz_dy(-1)%HI_density(-1,1) * pyramid_grid(1)%cellsize(1,1)
	px_ny_nz_dy(-1)%HI_column_density_out(-1,1) = HI_column_density_in + HI_column_density_cell
	
	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)	
   	      if (partition(i_level_py).eq.partition(i_level_py-1)) then
		  	HI_column_density_in = px_ny_nz_dy(-i_level_py+1)%HI_column_density_out(-i_primary,i_secondary)
            HI_column_density_cell = px_ny_nz_dy(-i_level_py)%HI_density(-i_primary,i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	px_ny_nz_dy(-i_level_py)%HI_column_density_out(-i_primary,i_secondary) = HI_column_density_in + &
                                                                                HI_column_density_cell 		
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2
		  	HI_column_density_in = px_ny_nz_dy(-i_level_py+1)%HI_column_density_out(-o_primary,o_secondary)
            HI_column_density_cell = px_ny_nz_dy(-i_level_py)%HI_density(-i_primary,i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	px_ny_nz_dy(-i_level_py)%HI_column_density_out(-i_primary,i_secondary) = HI_column_density_in + &
		  			  			  			  			  			  	        HI_column_density_cell			
          endif				  
		enddo	  
      enddo	  
    enddo	  
	  
  end subroutine px_ny_nz_dy_HI_column_density_calculation
  
  subroutine px_ny_nz_dz_HI_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	
	! i_level_py = 1
	HI_column_density_in = 0.0
    HI_column_density_cell = px_ny_nz_dz(-1)%HI_density(1,-1) * pyramid_grid(1)%cellsize(1,1)
	px_ny_nz_dz(-1)%HI_column_density_out(1,-1) = HI_column_density_in + HI_column_density_cell
	
	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)	
   	      if (partition(i_level_py).eq.partition(i_level_py-1)) then
		  	HI_column_density_in = px_ny_nz_dz(-i_level_py+1)%HI_column_density_out(i_primary,-i_secondary)
            HI_column_density_cell = px_ny_nz_dz(-i_level_py)%HI_density(i_primary,-i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	px_ny_nz_dz(-i_level_py)%HI_column_density_out(i_primary,-i_secondary) = HI_column_density_in + &
                                                                                HI_column_density_cell 		
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2
		  	HI_column_density_in = px_ny_nz_dz(-i_level_py+1)%HI_column_density_out(o_primary,-o_secondary)
            HI_column_density_cell = px_ny_nz_dz(-i_level_py)%HI_density(i_primary,-i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	px_ny_nz_dz(-i_level_py)%HI_column_density_out(i_primary,-i_secondary) = HI_column_density_in + &
		  			  			  			  			  			  	        HI_column_density_cell			
          endif				  
		enddo	  
      enddo	  
    enddo	  
	  
  end subroutine px_ny_nz_dz_HI_column_density_calculation

  subroutine nx_py_pz_dx_HI_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	
	! i_level_py = 1
	HI_column_density_in = 0.0
    HI_column_density_cell = nx_py_pz_dx(-1)%HI_density(1,1) * pyramid_grid(1)%cellsize(1,1)
	nx_py_pz_dx(-1)%HI_column_density_out(1,1) = HI_column_density_in + HI_column_density_cell

	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)	
   	      if (partition(i_level_py).eq.partition(i_level_py-1)) then
		  	HI_column_density_in = nx_py_pz_dx(-i_level_py+1)%HI_column_density_out(i_primary,i_secondary)
            HI_column_density_cell = nx_py_pz_dx(-i_level_py)%HI_density(i_primary,i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	nx_py_pz_dx(-i_level_py)%HI_column_density_out(i_primary,i_secondary) = HI_column_density_in + &
                                                                                HI_column_density_cell 		
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2
		  	HI_column_density_in = nx_py_pz_dx(-i_level_py+1)%HI_column_density_out(o_primary,o_secondary)
            HI_column_density_cell = nx_py_pz_dx(-i_level_py)%HI_density(i_primary,i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	nx_py_pz_dx(-i_level_py)%HI_column_density_out(i_primary,i_secondary) = HI_column_density_in + &
		  			  			  			  			  			  	        HI_column_density_cell			
          endif				  
		enddo	  
      enddo	  
    enddo	  
	  
  end subroutine nx_py_pz_dx_HI_column_density_calculation
  
  subroutine nx_py_pz_dy_HI_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	
	! i_level_py = 1
	HI_column_density_in = 0.0
    HI_column_density_cell = nx_py_pz_dy(1)%HI_density(1,-1) * pyramid_grid(1)%cellsize(1,1)
	nx_py_pz_dy(1)%HI_column_density_out(1,-1) = HI_column_density_in + HI_column_density_cell
	
	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)	
   	      if (partition(i_level_py).eq.partition(i_level_py-1)) then
		  	HI_column_density_in = nx_py_pz_dy(i_level_py-1)%HI_column_density_out(i_primary,-i_secondary)
            HI_column_density_cell = nx_py_pz_dy(i_level_py)%HI_density(i_primary,-i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	nx_py_pz_dy(i_level_py)%HI_column_density_out(i_primary,-i_secondary) = HI_column_density_in + &
                                                                                HI_column_density_cell 		
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2
		  	HI_column_density_in = nx_py_pz_dy(i_level_py-1)%HI_column_density_out(o_primary,-o_secondary)
            HI_column_density_cell = nx_py_pz_dy(i_level_py)%HI_density(i_primary,-i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	nx_py_pz_dy(i_level_py)%HI_column_density_out(i_primary,-i_secondary) = HI_column_density_in + &
		  			  			  			  			  			  	        HI_column_density_cell			
          endif				  
		enddo	  
      enddo	  
    enddo	  
	  
  end subroutine nx_py_pz_dy_HI_column_density_calculation
  
  subroutine nx_py_pz_dz_HI_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	
	! i_level_py = 1
	HI_column_density_in = 0.0
    HI_column_density_cell = nx_py_pz_dz(1)%HI_density(-1,1) * pyramid_grid(1)%cellsize(1,1)
	nx_py_pz_dz(1)%HI_column_density_out(-1,1) = HI_column_density_in + HI_column_density_cell
	
	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)	
   	      if (partition(i_level_py).eq.partition(i_level_py-1)) then
		  	HI_column_density_in = nx_py_pz_dz(i_level_py-1)%HI_column_density_out(-i_primary,i_secondary)
            HI_column_density_cell = nx_py_pz_dz(i_level_py)%HI_density(-i_primary,i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	nx_py_pz_dz(i_level_py)%HI_column_density_out(-i_primary,i_secondary) = HI_column_density_in + &
                                                                                HI_column_density_cell 		
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2
		  	HI_column_density_in = nx_py_pz_dz(i_level_py-1)%HI_column_density_out(-o_primary,o_secondary)
            HI_column_density_cell = nx_py_pz_dz(i_level_py)%HI_density(-i_primary,i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	nx_py_pz_dz(i_level_py)%HI_column_density_out(-i_primary,i_secondary) = HI_column_density_in + &
		  			  			  			  			  			  	        HI_column_density_cell			
          endif				  
		enddo	  
      enddo	  
    enddo	  
	  
  end subroutine nx_py_pz_dz_HI_column_density_calculation
     
  subroutine nx_py_nz_dx_HI_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	
	! i_level_py = 1
	HI_column_density_in = 0.0
    HI_column_density_cell = nx_py_nz_dx(-1)%HI_density(1,-1) * pyramid_grid(1)%cellsize(1,1)
	nx_py_nz_dx(-1)%HI_column_density_out(1,-1) = HI_column_density_in + HI_column_density_cell
	
	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)	
   	      if (partition(i_level_py).eq.partition(i_level_py-1)) then
		  	HI_column_density_in = nx_py_nz_dx(-i_level_py+1)%HI_column_density_out(i_primary,-i_secondary)
            HI_column_density_cell = nx_py_nz_dx(-i_level_py)%HI_density(i_primary,-i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	nx_py_nz_dx(-i_level_py)%HI_column_density_out(i_primary,-i_secondary) = HI_column_density_in + &
                                                                                HI_column_density_cell 		
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2
		  	HI_column_density_in = nx_py_nz_dx(-i_level_py+1)%HI_column_density_out(o_primary,-o_secondary)
            HI_column_density_cell = nx_py_nz_dx(-i_level_py)%HI_density(i_primary,-i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	nx_py_nz_dx(-i_level_py)%HI_column_density_out(i_primary,-i_secondary) = HI_column_density_in + &
		  			  			  			  			  			  	        HI_column_density_cell			
          endif				  
		enddo	  
      enddo	  
    enddo	  
	  
  end subroutine nx_py_nz_dx_HI_column_density_calculation
  
  subroutine nx_py_nz_dy_HI_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	
	! i_level_py = 1
	HI_column_density_in = 0.0
    HI_column_density_cell = nx_py_nz_dy(1)%HI_density(-1,-1) * pyramid_grid(1)%cellsize(1,1)
	nx_py_nz_dy(1)%HI_column_density_out(-1,-1) = HI_column_density_in + HI_column_density_cell
	
	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)	
   	      if (partition(i_level_py).eq.partition(i_level_py-1)) then
		  	HI_column_density_in = nx_py_nz_dy(i_level_py-1)%HI_column_density_out(-i_primary,-i_secondary)
            HI_column_density_cell = nx_py_nz_dy(i_level_py)%HI_density(-i_primary,-i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	nx_py_nz_dy(i_level_py)%HI_column_density_out(-i_primary,-i_secondary) = HI_column_density_in + &
                                                                                HI_column_density_cell 		
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2
		  	HI_column_density_in = nx_py_nz_dy(i_level_py-1)%HI_column_density_out(-o_primary,-o_secondary)
            HI_column_density_cell = nx_py_nz_dy(i_level_py)%HI_density(-i_primary,-i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	nx_py_nz_dy(i_level_py)%HI_column_density_out(-i_primary,-i_secondary) = HI_column_density_in + &
		  			  			  			  			  			  	        HI_column_density_cell			
          endif				  
		enddo	  
      enddo	  
    enddo	  
	  
  end subroutine nx_py_nz_dy_HI_column_density_calculation
  
  subroutine nx_py_nz_dz_HI_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	
	! i_level_py = 1
	HI_column_density_in = 0.0
    HI_column_density_cell = nx_py_nz_dz(-1)%HI_density(-1,1) * pyramid_grid(1)%cellsize(1,1)
	nx_py_nz_dz(-1)%HI_column_density_out(-1,1) = HI_column_density_in + HI_column_density_cell
	
	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)	
   	      if (partition(i_level_py).eq.partition(i_level_py-1)) then
		  	HI_column_density_in = nx_py_nz_dz(-i_level_py+1)%HI_column_density_out(-i_primary,i_secondary)
            HI_column_density_cell = nx_py_nz_dz(-i_level_py)%HI_density(-i_primary,i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	nx_py_nz_dz(-i_level_py)%HI_column_density_out(-i_primary,i_secondary) = HI_column_density_in + &
                                                                                HI_column_density_cell 		
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2
		  	HI_column_density_in = nx_py_nz_dz(-i_level_py+1)%HI_column_density_out(-o_primary,o_secondary)
            HI_column_density_cell = nx_py_nz_dz(-i_level_py)%HI_density(-i_primary,i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	nx_py_nz_dz(-i_level_py)%HI_column_density_out(-i_primary,i_secondary) = HI_column_density_in + &
		  			  			  			  			  			  	        HI_column_density_cell			
          endif				  
		enddo	  
      enddo	  
    enddo	  
	  
  end subroutine nx_py_nz_dz_HI_column_density_calculation
  	 
  subroutine nx_ny_pz_dx_HI_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	
	! i_level_py = 1
	HI_column_density_in = 0.0
    HI_column_density_cell = nx_ny_pz_dx(-1)%HI_density(-1,1) * pyramid_grid(1)%cellsize(1,1)
	nx_ny_pz_dx(-1)%HI_column_density_out(-1,1) = HI_column_density_in + HI_column_density_cell
	
	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)	
   	      if (partition(i_level_py).eq.partition(i_level_py-1)) then
		  	HI_column_density_in = nx_ny_pz_dx(-i_level_py+1)%HI_column_density_out(-i_primary,i_secondary)
            HI_column_density_cell = nx_ny_pz_dx(-i_level_py)%HI_density(-i_primary,i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	nx_ny_pz_dx(-i_level_py)%HI_column_density_out(-i_primary,i_secondary) = HI_column_density_in + &
                                                                                HI_column_density_cell 		
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2
		  	HI_column_density_in = nx_ny_pz_dx(-i_level_py+1)%HI_column_density_out(-o_primary,o_secondary)
            HI_column_density_cell = nx_ny_pz_dx(-i_level_py)%HI_density(-i_primary,i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	nx_ny_pz_dx(-i_level_py)%HI_column_density_out(-i_primary,i_secondary) = HI_column_density_in + &
		  			  			  			  			  			  	        HI_column_density_cell			
          endif				  
		enddo	  
      enddo	  
    enddo	  
	  
  end subroutine nx_ny_pz_dx_HI_column_density_calculation
  
  subroutine nx_ny_pz_dy_HI_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	
	! i_level_py = 1
	HI_column_density_in = 0.0
    HI_column_density_cell = nx_ny_pz_dy(-1)%HI_density(1,-1) * pyramid_grid(1)%cellsize(1,1)
	nx_ny_pz_dy(-1)%HI_column_density_out(1,-1) = HI_column_density_in + HI_column_density_cell
	
	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)	
   	      if (partition(i_level_py).eq.partition(i_level_py-1)) then
		  	HI_column_density_in = nx_ny_pz_dy(-i_level_py+1)%HI_column_density_out(i_primary,-i_secondary)
            HI_column_density_cell = nx_ny_pz_dy(-i_level_py)%HI_density(i_primary,-i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	nx_ny_pz_dy(-i_level_py)%HI_column_density_out(i_primary,-i_secondary) = HI_column_density_in + &
                                                                                HI_column_density_cell 		
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2
		  	HI_column_density_in = nx_ny_pz_dy(-i_level_py+1)%HI_column_density_out(o_primary,-o_secondary)
            HI_column_density_cell = nx_ny_pz_dy(-i_level_py)%HI_density(i_primary,-i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	nx_ny_pz_dy(-i_level_py)%HI_column_density_out(i_primary,-i_secondary) = HI_column_density_in + &
		  			  			  			  			  			  	        HI_column_density_cell			
          endif				  
		enddo	  
      enddo	  
    enddo	  
	  
  end subroutine nx_ny_pz_dy_HI_column_density_calculation
  
  subroutine nx_ny_pz_dz_HI_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	
	! i_level_py = 1
	HI_column_density_in = 0.0
    HI_column_density_cell = nx_ny_pz_dz(1)%HI_density(-1,-1) * pyramid_grid(1)%cellsize(1,1)
	nx_ny_pz_dz(1)%HI_column_density_out(-1,-1) = HI_column_density_in + HI_column_density_cell
	
	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)	
   	      if (partition(i_level_py).eq.partition(i_level_py-1)) then
		  	HI_column_density_in = nx_ny_pz_dz(i_level_py-1)%HI_column_density_out(-i_primary,-i_secondary)
            HI_column_density_cell = nx_ny_pz_dz(i_level_py)%HI_density(-i_primary,-i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	nx_ny_pz_dz(i_level_py)%HI_column_density_out(-i_primary,-i_secondary) = HI_column_density_in + &
                                                                                HI_column_density_cell 		
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2
		  	HI_column_density_in = nx_ny_pz_dz(i_level_py-1)%HI_column_density_out(-o_primary,-o_secondary)
            HI_column_density_cell = nx_ny_pz_dz(i_level_py)%HI_density(-i_primary,-i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	nx_ny_pz_dz(i_level_py)%HI_column_density_out(-i_primary,-i_secondary) = HI_column_density_in + &
		  			  			  			  			  			  	        HI_column_density_cell			
          endif				  
		enddo	  
      enddo	  
    enddo	  
	  
  end subroutine nx_ny_pz_dz_HI_column_density_calculation
     
  subroutine nx_ny_nz_dx_HI_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	
	! i_level_py = 1
	HI_column_density_in = 0.0
    HI_column_density_cell = nx_ny_nz_dx(-1)%HI_density(-1,-1) * pyramid_grid(1)%cellsize(1,1)
	nx_ny_nz_dx(-1)%HI_column_density_out(-1,-1) = HI_column_density_in + HI_column_density_cell
	
	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)	
   	      if (partition(i_level_py).eq.partition(i_level_py-1)) then
		  	HI_column_density_in = nx_ny_nz_dx(-i_level_py+1)%HI_column_density_out(-i_primary,-i_secondary)
            HI_column_density_cell = nx_ny_nz_dx(-i_level_py)%HI_density(-i_primary,-i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	nx_ny_nz_dx(-i_level_py)%HI_column_density_out(-i_primary,-i_secondary) = HI_column_density_in + &
                                                                                HI_column_density_cell 		
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2
		  	HI_column_density_in = nx_ny_nz_dx(-i_level_py+1)%HI_column_density_out(-o_primary,-o_secondary)
            HI_column_density_cell = nx_ny_nz_dx(-i_level_py)%HI_density(-i_primary,-i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	nx_ny_nz_dx(-i_level_py)%HI_column_density_out(-i_primary,-i_secondary) = HI_column_density_in + &
		  			  			  			  			  			  	        HI_column_density_cell			
          endif				  
		enddo	  
      enddo	  
    enddo	  
	  
  end subroutine nx_ny_nz_dx_HI_column_density_calculation
  
  subroutine nx_ny_nz_dy_HI_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	
	! i_level_py = 1
	HI_column_density_in = 0.0
    HI_column_density_cell = nx_ny_nz_dy(-1)%HI_density(-1,-1) * pyramid_grid(1)%cellsize(1,1)
	nx_ny_nz_dy(-1)%HI_column_density_out(-1,-1) = HI_column_density_in + HI_column_density_cell
	
	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)	
   	      if (partition(i_level_py).eq.partition(i_level_py-1)) then
		  	HI_column_density_in = nx_ny_nz_dy(-i_level_py+1)%HI_column_density_out(-i_primary,-i_secondary)
            HI_column_density_cell = nx_ny_nz_dy(-i_level_py)%HI_density(-i_primary,-i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	nx_ny_nz_dy(-i_level_py)%HI_column_density_out(-i_primary,-i_secondary) = HI_column_density_in + &
                                                                                HI_column_density_cell 		
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2
		  	HI_column_density_in = nx_ny_nz_dy(-i_level_py+1)%HI_column_density_out(-o_primary,-o_secondary)
            HI_column_density_cell = nx_ny_nz_dy(-i_level_py)%HI_density(-i_primary,-i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	nx_ny_nz_dy(-i_level_py)%HI_column_density_out(-i_primary,-i_secondary) = HI_column_density_in + &
		  			  			  			  			  			  	        HI_column_density_cell			
          endif				  
		enddo	  
      enddo	  
    enddo	  
	  
  end subroutine nx_ny_nz_dy_HI_column_density_calculation
  
  subroutine nx_ny_nz_dz_HI_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	
	! i_level_py = 1
	HI_column_density_in = 0.0
    HI_column_density_cell = nx_ny_nz_dz(-1)%HI_density(-1,-1) * pyramid_grid(1)%cellsize(1,1)
	nx_ny_nz_dz(-1)%HI_column_density_out(-1,-1) = HI_column_density_in + HI_column_density_cell
	
	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)	
   	      if (partition(i_level_py).eq.partition(i_level_py-1)) then
		  	HI_column_density_in = nx_ny_nz_dz(-i_level_py+1)%HI_column_density_out(-i_primary,-i_secondary)
            HI_column_density_cell = nx_ny_nz_dz(-i_level_py)%HI_density(-i_primary,-i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	nx_ny_nz_dz(-i_level_py)%HI_column_density_out(-i_primary,-i_secondary) = HI_column_density_in + &
                                                                                HI_column_density_cell 		
          else
            o_primary = (i_primary+1)/2
            o_secondary = (i_secondary+1)/2
		  	HI_column_density_in = nx_ny_nz_dz(-i_level_py+1)%HI_column_density_out(-o_primary,-o_secondary)
            HI_column_density_cell = nx_ny_nz_dz(-i_level_py)%HI_density(-i_primary,-i_secondary) * &
			                         pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  	nx_ny_nz_dz(-i_level_py)%HI_column_density_out(-i_primary,-i_secondary) = HI_column_density_in + &
		  			  			  			  			  			  	        HI_column_density_cell			
          endif				  
		enddo	  
      enddo	  
    enddo	  
	  
  end subroutine nx_ny_nz_dz_HI_column_density_calculation
    	 
end module pyramid_column_density