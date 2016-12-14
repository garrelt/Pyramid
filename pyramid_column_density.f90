module pyramid_column_density
	
    use precision, only: dp
    use input, only: level_py, partition, layer
    use array, only: pyramid_grid, parent_position, &
                     px_py_pz_dx, px_py_pz_dy, px_py_pz_dz, &
                     px_py_nz_dx, px_py_nz_dy, px_py_nz_dz, &
                     px_ny_pz_dx, px_ny_pz_dy, px_ny_pz_dz, &	  
                     px_ny_nz_dx, px_ny_nz_dy, px_ny_nz_dz, & 	   
                     nx_py_pz_dx, nx_py_pz_dy, nx_py_pz_dz, &	  
                     nx_py_nz_dx, nx_py_nz_dy, nx_py_nz_dz, &	
                     nx_ny_pz_dx, nx_ny_pz_dy, nx_ny_pz_dz, &	  
                     nx_ny_nz_dx, nx_ny_nz_dy, nx_ny_nz_dz		   	   
		
contains
	
  subroutine pyramid_domain_column_density_calculation()
		
    call px_py_pz_dx_column_density_calculation ()
    call px_py_pz_dy_column_density_calculation ()
    call px_py_pz_dz_column_density_calculation ()	
    call px_py_nz_dx_column_density_calculation ()
    call px_py_nz_dy_column_density_calculation ()	
    call px_py_nz_dz_column_density_calculation ()		
    call px_ny_pz_dx_column_density_calculation ()	
    call px_ny_pz_dy_column_density_calculation ()	
    call px_ny_pz_dz_column_density_calculation ()	
    call px_ny_nz_dx_column_density_calculation ()
    call px_ny_nz_dy_column_density_calculation ()
    call px_ny_nz_dz_column_density_calculation ()
    call nx_py_pz_dx_column_density_calculation ()
    call nx_py_pz_dy_column_density_calculation ()
    call nx_py_pz_dz_column_density_calculation ()
    call nx_py_nz_dx_column_density_calculation ()	
    call nx_py_nz_dy_column_density_calculation ()	
    call nx_py_nz_dz_column_density_calculation ()	
    call nx_ny_pz_dx_column_density_calculation ()
    call nx_ny_pz_dy_column_density_calculation ()
    call nx_ny_pz_dz_column_density_calculation ()	
    call nx_ny_nz_dx_column_density_calculation ()
    call nx_ny_nz_dy_column_density_calculation ()
    call nx_ny_nz_dz_column_density_calculation ()		
								
  end subroutine pyramid_domain_column_density_calculation	
  
  !!!!!!!!!!
  
  subroutine px_py_pz_dx_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
    integer :: parent_primary,parent_secondary,parent_level_py
    real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
    real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_cell
    real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_cell
    real(kind=dp) :: cellsize
	
    ! i_level_py = 1
    px_py_pz_dx(1)%HI_column_density_in(1,1) = 0
    px_py_pz_dx(1)%HeI_column_density_in(1,1) = 0
    px_py_pz_dx(1)%HeII_column_density_in(1,1) = 0
  
    cellsize = pyramid_grid(1)%cellsize(1,1)

    px_py_pz_dx(1)%HI_column_density_out(1,1) = px_py_pz_dx(1)%HI_density(1,1) * cellsize	
    px_py_pz_dx(1)%HeI_column_density_out(1,1) = px_py_pz_dx(1)%HeI_density(1,1) * cellsize
    px_py_pz_dx(1)%HeII_column_density_out(1,1) = px_py_pz_dx(1)%HeII_density(1,1) * cellsize
    write(*,*) 'column density near source is ', px_py_pz_dx(1)%HI_column_density_out(1,1)				
    ! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)
			
          cellsize = pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  
          if (layer(i_level_py).ne.layer(i_level_py-1)) then
			  
            HI_column_density_in = 0
            HeI_column_density_in = 0
            HeII_column_density_in = 0
			  		  
            do parent_level_py = 1,i_level_py-1
			
              parent_primary = parent_position(i_level_py,i_primary)%with_parent_height(parent_level_py)
              parent_secondary = parent_position(i_level_py,i_secondary)%with_parent_height(parent_level_py)
			
              HI_column_density_in = HI_column_density_in + &
                                px_py_pz_dx(parent_level_py)%HI_density(parent_primary,parent_secondary) * cellsize
				              
              HeI_column_density_in = HeI_column_density_in + &
                                px_py_pz_dx(parent_level_py)%HeI_density(parent_primary,parent_secondary) * cellsize
							  							  
              HeII_column_density_in = HeII_column_density_in + &
                                px_py_pz_dx(parent_level_py)%HeII_density(parent_primary,parent_secondary) * cellsize
							  
            enddo	
		  
            px_py_pz_dx(i_level_py)%HI_column_density_in(i_primary,i_secondary) = HI_column_density_in	
            px_py_pz_dx(i_level_py)%HeI_column_density_in(i_primary,i_secondary) = HeI_column_density_in
            px_py_pz_dx(i_level_py)%HeII_column_density_in(i_primary,i_secondary) = HeII_column_density_in		
		
          else
			  
            px_py_pz_dx(i_level_py)%HI_column_density_in(i_primary,i_secondary) = &
					  px_py_pz_dx(i_level_py-1)%HI_column_density_out(i_primary,i_secondary)
            px_py_pz_dx(i_level_py)%HeI_column_density_in(i_primary,i_secondary) = &
					  px_py_pz_dx(i_level_py-1)%HeI_column_density_out(i_primary,i_secondary)		
            px_py_pz_dx(i_level_py)%HeII_column_density_in(i_primary,i_secondary) = &
					  px_py_pz_dx(i_level_py-1)%HeII_column_density_out(i_primary,i_secondary)
			  
          endif
		  			    			
          px_py_pz_dx(i_level_py)%HI_column_density_out(i_primary,i_secondary) = &
                    px_py_pz_dx(i_level_py)%HI_column_density_in(i_primary,i_secondary) + &
                    px_py_pz_dx(i_level_py)%HI_density(i_primary,i_secondary) * cellsize
          px_py_pz_dx(i_level_py)%HeI_column_density_out(i_primary,i_secondary) = &
                    px_py_pz_dx(i_level_py)%HeI_column_density_in(i_primary,i_secondary) + &
                    px_py_pz_dx(i_level_py)%HeI_density(i_primary,i_secondary) * cellsize	  
          px_py_pz_dx(i_level_py)%HeII_column_density_out(i_primary,i_secondary) = &
                    px_py_pz_dx(i_level_py)%HeII_column_density_in(i_primary,i_secondary) + &
                    px_py_pz_dx(i_level_py)%HeII_density(i_primary,i_secondary) * cellsize		  

        enddo	  
      enddo		    
    enddo	

  end subroutine px_py_pz_dx_column_density_calculation
  
  subroutine px_py_pz_dy_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	integer :: parent_primary,parent_secondary,parent_level_py
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_cell
	real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_cell
    real(kind=dp) :: cellsize
			
	! i_level_py = 1
    px_py_pz_dy(1)%HI_column_density_in(1,1) = 0
    px_py_pz_dy(1)%HeI_column_density_in(1,1) = 0
    px_py_pz_dy(1)%HeII_column_density_in(1,1) = 0

    cellsize = pyramid_grid(1)%cellsize(1,1)
	
    px_py_pz_dy(1)%HI_column_density_out(1,1) = px_py_pz_dy(1)%HI_density(1,1) * cellsize
	px_py_pz_dy(1)%HeI_column_density_out(1,1) = px_py_pz_dy(1)%HeI_density(1,1) * cellsize
	px_py_pz_dy(1)%HeII_column_density_out(1,1) = px_py_pz_dy(1)%HeII_density(1,1) * cellsize

	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)
			
  		  cellsize = pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  
  		  if (layer(i_level_py).ne.layer(i_level_py-1)) then
			  
            HI_column_density_in = 0
            HeI_column_density_in = 0
            HeII_column_density_in = 0
		  		  
            do parent_level_py = 1,i_level_py-1
			
              parent_primary = parent_position(i_level_py,i_primary)%with_parent_height(parent_level_py)
              parent_secondary = parent_position(i_level_py,i_secondary)%with_parent_height(parent_level_py)
			
            HI_column_density_in = HI_column_density_in + &
                              px_py_pz_dy(parent_level_py)%HI_density(parent_primary,parent_secondary) * cellsize

		    HeI_column_density_in = HeI_column_density_in + &
                              px_py_pz_dy(parent_level_py)%HeI_density(parent_primary,parent_secondary) * cellsize

		    HeII_column_density_in = HeII_column_density_in + &
                              px_py_pz_dy(parent_level_py)%HeII_density(parent_primary,parent_secondary) * cellsize

          enddo	
		  
            px_py_pz_dy(i_level_py)%HI_column_density_in(i_primary,i_secondary) = HI_column_density_in	
            px_py_pz_dy(i_level_py)%HeI_column_density_in(i_primary,i_secondary) = HeI_column_density_in
            px_py_pz_dy(i_level_py)%HeII_column_density_in(i_primary,i_secondary) = HeII_column_density_in		
	
      else
		  
        px_py_pz_dy(i_level_py)%HI_column_density_in(i_primary,i_secondary) = &
				  px_py_pz_dy(i_level_py-1)%HI_column_density_out(i_primary,i_secondary)
        px_py_pz_dy(i_level_py)%HeI_column_density_in(i_primary,i_secondary) = &
				  px_py_pz_dy(i_level_py-1)%HeI_column_density_out(i_primary,i_secondary)		
        px_py_pz_dy(i_level_py)%HeII_column_density_in(i_primary,i_secondary) = &
				  px_py_pz_dy(i_level_py-1)%HeII_column_density_out(i_primary,i_secondary)
		  
      endif
	  	  			
	  px_py_pz_dy(i_level_py)%HI_column_density_out(i_primary,i_secondary) = &
                       px_py_pz_dy(i_level_py)%HI_column_density_in(i_primary,i_secondary) + &
                       px_py_pz_dy(i_level_py)%HI_density(i_primary,i_secondary) * cellsize
	  px_py_pz_dy(i_level_py)%HeI_column_density_out(i_primary,i_secondary) = &
                       px_py_pz_dy(i_level_py)%HeI_column_density_in(i_primary,i_secondary) + &
                       px_py_pz_dy(i_level_py)%HeI_density(i_primary,i_secondary) * cellsize	  
	  px_py_pz_dy(i_level_py)%HeII_column_density_out(i_primary,i_secondary) = &
                       px_py_pz_dy(i_level_py)%HeII_column_density_in(i_primary,i_secondary) + &
                       px_py_pz_dy(i_level_py)%HeII_density(i_primary,i_secondary) * cellsize	

			   		enddo	  
			         enddo		    
			       enddo		  
			
  end subroutine px_py_pz_dy_column_density_calculation
  
  subroutine px_py_pz_dz_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	integer :: parent_primary,parent_secondary,parent_level_py
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_cell
	real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_cell
    real(kind=dp) :: cellsize
		
	! i_level_py = 1
    px_py_pz_dz(1)%HI_column_density_in(1,1) = 0
    px_py_pz_dz(1)%HeI_column_density_in(1,1) = 0
    px_py_pz_dz(1)%HeII_column_density_in(1,1) = 0

    cellsize = pyramid_grid(1)%cellsize(1,1)
	
    px_py_pz_dz(1)%HI_column_density_out(1,1) = px_py_pz_dz(1)%HI_density(1,1) * cellsize
	px_py_pz_dz(1)%HeI_column_density_out(1,1) = px_py_pz_dz(1)%HeI_density(1,1) * cellsize
	px_py_pz_dz(1)%HeII_column_density_out(1,1) = px_py_pz_dz(1)%HeII_density(1,1) * cellsize

	! i_level_py > 1
    do i_level_py = 2,level_py
      do i_primary = 1,partition(i_level_py)
        do i_secondary = 1,partition(i_level_py)
			
  		  cellsize = pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
		  
  		  if (layer(i_level_py).ne.layer(i_level_py-1)) then
			  
            HI_column_density_in = 0
            HeI_column_density_in = 0
            HeII_column_density_in = 0
		  		  
            do parent_level_py = 1,i_level_py-1
			
              parent_primary = parent_position(i_level_py,i_primary)%with_parent_height(parent_level_py)
              parent_secondary = parent_position(i_level_py,i_secondary)%with_parent_height(parent_level_py)
			
              HI_column_density_in = HI_column_density_in + &
                              px_py_pz_dz(parent_level_py)%HI_density(parent_primary,parent_secondary) * cellsize

              HeI_column_density_in = HeI_column_density_in + &
                              px_py_pz_dz(parent_level_py)%HeI_density(parent_primary,parent_secondary) * cellsize
  
              HeII_column_density_in = HeII_column_density_in + &
                              px_py_pz_dz(parent_level_py)%HeII_density(parent_primary,parent_secondary) * cellsize
		  
          enddo	
		  
		  px_py_pz_dz(i_level_py)%HI_column_density_in(i_primary,i_secondary) = HI_column_density_in	
		  px_py_pz_dz(i_level_py)%HeI_column_density_in(i_primary,i_secondary) = HeI_column_density_in
		  px_py_pz_dz(i_level_py)%HeII_column_density_in(i_primary,i_secondary) = HeII_column_density_in		
	
      else
		  
        px_py_pz_dz(i_level_py)%HI_column_density_in(i_primary,i_secondary) = &
				  px_py_pz_dz(i_level_py-1)%HI_column_density_out(i_primary,i_secondary)
        px_py_pz_dz(i_level_py)%HeI_column_density_in(i_primary,i_secondary) = &
				  px_py_pz_dz(i_level_py-1)%HeI_column_density_out(i_primary,i_secondary)		
        px_py_pz_dz(i_level_py)%HeII_column_density_in(i_primary,i_secondary) = &
				  px_py_pz_dz(i_level_py-1)%HeII_column_density_out(i_primary,i_secondary)
		  
      endif
	  	  			
	  px_py_pz_dz(i_level_py)%HI_column_density_out(i_primary,i_secondary) = &
                       px_py_pz_dz(i_level_py)%HI_column_density_in(i_primary,i_secondary) + &
                       px_py_pz_dz(i_level_py)%HI_density(i_primary,i_secondary) * cellsize
	  px_py_pz_dz(i_level_py)%HeI_column_density_out(i_primary,i_secondary) = &
                       px_py_pz_dz(i_level_py)%HeI_column_density_in(i_primary,i_secondary) + &
                       px_py_pz_dz(i_level_py)%HeI_density(i_primary,i_secondary) * cellsize	  
	  px_py_pz_dz(i_level_py)%HeII_column_density_out(i_primary,i_secondary) = &
                       px_py_pz_dz(i_level_py)%HeII_column_density_in(i_primary,i_secondary) + &
                       px_py_pz_dz(i_level_py)%HeII_density(i_primary,i_secondary) * cellsize		  
		  	
		enddo	  
      enddo		    
    enddo	  
	  
  end subroutine px_py_pz_dz_column_density_calculation
     
  subroutine px_py_nz_dx_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	integer :: parent_primary,parent_secondary,parent_level_py
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_cell
	real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_cell
    real(kind=dp) :: cellsize
		
	! i_level_py = 1
    px_py_nz_dx(1)%HI_column_density_in(1,-1) = 0
    px_py_nz_dx(1)%HeI_column_density_in(1,-1) = 0
    px_py_nz_dx(1)%HeII_column_density_in(1,-1) = 0
  
    cellsize = pyramid_grid(1)%cellsize(1,1)

    px_py_nz_dx(1)%HI_column_density_out(1,-1) = px_py_nz_dx(1)%HI_density(1,-1) * cellsize	
	px_py_nz_dx(1)%HeI_column_density_out(1,-1) = px_py_nz_dx(1)%HeI_density(1,-1) * cellsize
	px_py_nz_dx(1)%HeII_column_density_out(1,-1) = px_py_nz_dx(1)%HeII_density(1,-1) * cellsize
		
	! i_level_py > 1
	do i_level_py = 2,level_py
	  do i_primary = 1,partition(i_level_py)
	    do i_secondary = 1,partition(i_level_py)
		
		  cellsize = pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
	  
		  if (layer(i_level_py).ne.layer(i_level_py-1)) then
		  
	        HI_column_density_in = 0
	        HeI_column_density_in = 0
	        HeII_column_density_in = 0
		  		  
	        do parent_level_py = 1,i_level_py-1
		
	          parent_primary = parent_position(i_level_py,i_primary)%with_parent_height(parent_level_py)
	          parent_secondary = parent_position(i_level_py,i_secondary)%with_parent_height(parent_level_py)
		
	          HI_column_density_in = HI_column_density_in + &
	                            px_py_nz_dx(parent_level_py)%HI_density(parent_primary,-parent_secondary) * cellsize
			              
	          HeI_column_density_in = HeI_column_density_in + &
	                            px_py_nz_dx(parent_level_py)%HeI_density(parent_primary,-parent_secondary) * cellsize
						  							  
	          HeII_column_density_in = HeII_column_density_in + &
	                            px_py_nz_dx(parent_level_py)%HeII_density(parent_primary,-parent_secondary) * cellsize
						  
	        enddo	
	  
	        px_py_nz_dx(i_level_py)%HI_column_density_in(i_primary,-i_secondary) = HI_column_density_in	
	        px_py_nz_dx(i_level_py)%HeI_column_density_in(i_primary,-i_secondary) = HeI_column_density_in
	        px_py_nz_dx(i_level_py)%HeII_column_density_in(i_primary,-i_secondary) = HeII_column_density_in		
	
	      else
		  
	        px_py_nz_dx(i_level_py)%HI_column_density_in(i_primary,-i_secondary) = &
					  px_py_nz_dx(i_level_py-1)%HI_column_density_out(i_primary,-i_secondary)
	        px_py_nz_dx(i_level_py)%HeI_column_density_in(i_primary,-i_secondary) = &
					  px_py_nz_dx(i_level_py-1)%HeI_column_density_out(i_primary,-i_secondary)		
	        px_py_nz_dx(i_level_py)%HeII_column_density_in(i_primary,-i_secondary) = &
					  px_py_nz_dx(i_level_py-1)%HeII_column_density_out(i_primary,-i_secondary)
		  
	      endif
	  			    			
		  px_py_nz_dx(i_level_py)%HI_column_density_out(i_primary,-i_secondary) = &
	                       px_py_nz_dx(i_level_py)%HI_column_density_in(i_primary,-i_secondary) + &
	                       px_py_nz_dx(i_level_py)%HI_density(i_primary,-i_secondary) * cellsize
		  px_py_nz_dx(i_level_py)%HeI_column_density_out(i_primary,-i_secondary) = &
	                       px_py_nz_dx(i_level_py)%HeI_column_density_in(i_primary,-i_secondary) + &
	                       px_py_nz_dx(i_level_py)%HeI_density(i_primary,-i_secondary) * cellsize	  
		  px_py_nz_dx(i_level_py)%HeII_column_density_out(i_primary,-i_secondary) = &
	                       px_py_nz_dx(i_level_py)%HeII_column_density_in(i_primary,-i_secondary) + &
	                       px_py_nz_dx(i_level_py)%HeII_density(i_primary,-i_secondary) * cellsize		  
	  	
		enddo	  
	  enddo	  	  
	enddo  
	  
  end subroutine px_py_nz_dx_column_density_calculation
  
  subroutine px_py_nz_dy_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	integer :: parent_primary,parent_secondary,parent_level_py
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_cell
	real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_cell
    real(kind=dp) :: cellsize

	! i_level_py = 1
    px_py_nz_dy(1)%HI_column_density_in(-1,1) = 0
    px_py_nz_dy(1)%HeI_column_density_in(-1,1) = 0
    px_py_nz_dy(1)%HeII_column_density_in(-1,1) = 0
	
    cellsize = pyramid_grid(1)%cellsize(1,1)

    px_py_nz_dy(1)%HI_column_density_out(-1,1) = px_py_nz_dy(1)%HI_density(-1,1) * cellsize	
	px_py_nz_dy(1)%HeI_column_density_out(-1,1) = px_py_nz_dy(1)%HeI_density(-1,1) * cellsize
	px_py_nz_dy(1)%HeII_column_density_out(-1,1) = px_py_nz_dy(1)%HeII_density(-1,1) * cellsize
		
	! i_level_py > 1
	do i_level_py = 2,level_py
	  do i_primary = 1,partition(i_level_py)
	    do i_secondary = 1,partition(i_level_py)
		
		  cellsize = pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
	  
		  if (layer(i_level_py).ne.layer(i_level_py-1)) then
		  
	        HI_column_density_in = 0
	        HeI_column_density_in = 0
	        HeII_column_density_in = 0
		  		  
	        do parent_level_py = 1,i_level_py-1
		
	          parent_primary = parent_position(i_level_py,i_primary)%with_parent_height(parent_level_py)
	          parent_secondary = parent_position(i_level_py,i_secondary)%with_parent_height(parent_level_py)
		
	          HI_column_density_in = HI_column_density_in + &
	                            px_py_nz_dy(parent_level_py)%HI_density(-parent_primary,parent_secondary) * cellsize
			              
	          HeI_column_density_in = HeI_column_density_in + &
	                            px_py_nz_dy(parent_level_py)%HeI_density(-parent_primary,parent_secondary) * cellsize
						  							  
	          HeII_column_density_in = HeII_column_density_in + &
	                            px_py_nz_dy(parent_level_py)%HeII_density(-parent_primary,parent_secondary) * cellsize 
						  
	        enddo	
	  
	        px_py_nz_dy(i_level_py)%HI_column_density_in(-i_primary,i_secondary) = HI_column_density_in	
	        px_py_nz_dy(i_level_py)%HeI_column_density_in(-i_primary,i_secondary) = HeI_column_density_in
	        px_py_nz_dy(i_level_py)%HeII_column_density_in(-i_primary,i_secondary) = HeII_column_density_in		
	
	      else
		  
	        px_py_nz_dy(i_level_py)%HI_column_density_in(-i_primary,i_secondary) = &
					  px_py_nz_dy(i_level_py-1)%HI_column_density_out(-i_primary,i_secondary)
	        px_py_nz_dy(i_level_py)%HeI_column_density_in(-i_primary,i_secondary) = &
					  px_py_nz_dy(i_level_py-1)%HeI_column_density_out(-i_primary,i_secondary)		
	        px_py_nz_dy(i_level_py)%HeII_column_density_in(-i_primary,i_secondary) = &
					  px_py_nz_dy(i_level_py-1)%HeII_column_density_out(-i_primary,i_secondary)
		  
	      endif
	  			    			
		  px_py_nz_dy(i_level_py)%HI_column_density_out(-i_primary,i_secondary) = &
	                       px_py_nz_dy(i_level_py)%HI_column_density_in(-i_primary,i_secondary) + &
	                       px_py_nz_dy(i_level_py)%HI_density(-i_primary,i_secondary) * cellsize
		  px_py_nz_dy(i_level_py)%HeI_column_density_out(-i_primary,i_secondary) = &
	                       px_py_nz_dy(i_level_py)%HeI_column_density_in(-i_primary,i_secondary) + &
	                       px_py_nz_dy(i_level_py)%HeI_density(-i_primary,i_secondary) * cellsize	  
		  px_py_nz_dy(i_level_py)%HeII_column_density_out(-i_primary,i_secondary) = &
	                       px_py_nz_dy(i_level_py)%HeII_column_density_in(-i_primary,i_secondary) + &
	                       px_py_nz_dy(i_level_py)%HeII_density(-i_primary,i_secondary) * cellsize		  
	  	
		enddo	  
	  enddo	  	 	  	     
    enddo	  
	  
  end subroutine px_py_nz_dy_column_density_calculation
  
  subroutine px_py_nz_dz_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	integer :: parent_primary,parent_secondary,parent_level_py
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_cell
	real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_cell
    real(kind=dp) :: cellsize

	! i_level_py = 1
    px_py_nz_dz(-1)%HI_column_density_in(1,1) = 0
    px_py_nz_dz(-1)%HeI_column_density_in(1,1) = 0
    px_py_nz_dz(-1)%HeII_column_density_in(1,1) = 0
  
    cellsize = pyramid_grid(1)%cellsize(1,1)

    px_py_nz_dz(-1)%HI_column_density_out(1,1) = px_py_nz_dz(-1)%HI_density(1,1) * cellsize	
	px_py_nz_dz(-1)%HeI_column_density_out(1,1) = px_py_nz_dz(-1)%HeI_density(1,1) * cellsize
	px_py_nz_dz(-1)%HeII_column_density_out(1,1) = px_py_nz_dz(-1)%HeII_density(1,1) * cellsize
		
	! i_level_py > 1
	do i_level_py = 2,level_py
	  do i_primary = 1,partition(i_level_py)
	    do i_secondary = 1,partition(i_level_py)
		
		  cellsize = pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
	  
		  if (layer(i_level_py).ne.layer(i_level_py-1)) then
		  
	        HI_column_density_in = 0
	        HeI_column_density_in = 0
	        HeII_column_density_in = 0
		  		  
	        do parent_level_py = 1,i_level_py-1
		
	          parent_primary = parent_position(i_level_py,i_primary)%with_parent_height(parent_level_py)
	          parent_secondary = parent_position(i_level_py,i_secondary)%with_parent_height(parent_level_py)
		
	          HI_column_density_in = HI_column_density_in + &
	                            px_py_nz_dz(-parent_level_py)%HI_density(parent_primary,parent_secondary) * cellsize
			              
	          HeI_column_density_in = HeI_column_density_in + &
	                            px_py_nz_dz(-parent_level_py)%HeI_density(parent_primary,parent_secondary) * cellsize
						  							  
	          HeII_column_density_in = HeII_column_density_in + &
	                            px_py_nz_dz(-parent_level_py)%HeII_density(parent_primary,parent_secondary) * cellsize 
						  
	        enddo	
	  
	        px_py_nz_dz(-i_level_py)%HI_column_density_in(i_primary,i_secondary) = HI_column_density_in	
	        px_py_nz_dz(-i_level_py)%HeI_column_density_in(i_primary,i_secondary) = HeI_column_density_in
	        px_py_nz_dz(-i_level_py)%HeII_column_density_in(i_primary,i_secondary) = HeII_column_density_in		
	
	      else
		  
	        px_py_nz_dz(-i_level_py)%HI_column_density_in(i_primary,i_secondary) = &
					  px_py_nz_dz(-i_level_py+1)%HI_column_density_out(i_primary,i_secondary)
	        px_py_nz_dz(-i_level_py)%HeI_column_density_in(i_primary,i_secondary) = &
					  px_py_nz_dz(-i_level_py+1)%HeI_column_density_out(i_primary,i_secondary)		
	        px_py_nz_dz(-i_level_py)%HeII_column_density_in(i_primary,i_secondary) = &
					  px_py_nz_dz(-i_level_py+1)%HeII_column_density_out(i_primary,i_secondary)
		  
	      endif
	  			    			
		  px_py_nz_dz(-i_level_py)%HI_column_density_out(i_primary,i_secondary) = &
	                       px_py_nz_dz(-i_level_py)%HI_column_density_in(i_primary,i_secondary) + &
	                       px_py_nz_dz(-i_level_py)%HI_density(i_primary,i_secondary) * cellsize
		  px_py_nz_dz(-i_level_py)%HeI_column_density_out(i_primary,i_secondary) = &
	                       px_py_nz_dz(-i_level_py)%HeI_column_density_in(i_primary,i_secondary) + &
	                       px_py_nz_dz(-i_level_py)%HeI_density(i_primary,i_secondary) * cellsize	  
		  px_py_nz_dz(-i_level_py)%HeII_column_density_out(i_primary,i_secondary) = &
	                       px_py_nz_dz(-i_level_py)%HeII_column_density_in(i_primary,i_secondary) + &
	                       px_py_nz_dz(-i_level_py)%HeII_density(i_primary,i_secondary) * cellsize		  
	  	
		enddo	  
	  enddo	   	  	   
    enddo	  
	  
  end subroutine px_py_nz_dz_column_density_calculation
  	 
  subroutine px_ny_pz_dx_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	integer :: parent_primary,parent_secondary,parent_level_py
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_cell
	real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_cell
    real(kind=dp) :: cellsize

	! i_level_py = 1
	px_ny_pz_dx(1)%HI_column_density_in(-1,1) = 0
	px_ny_pz_dx(1)%HeI_column_density_in(-1,1) = 0
	px_ny_pz_dx(1)%HeII_column_density_in(-1,1) = 0

	cellsize = pyramid_grid(1)%cellsize(1,1)

	px_ny_pz_dx(1)%HI_column_density_out(-1,1) = px_ny_pz_dx(1)%HI_density(-1,1) * cellsize	
	px_ny_pz_dx(1)%HeI_column_density_out(-1,1) = px_ny_pz_dx(1)%HeI_density(-1,1) * cellsize
	px_ny_pz_dx(1)%HeII_column_density_out(-1,1) = px_ny_pz_dx(1)%HeII_density(-1,1) * cellsize
		
	! i_level_py > 1
	do i_level_py = 2,level_py
	  do i_primary = 1,partition(i_level_py)
	    do i_secondary = 1,partition(i_level_py)
		
		  cellsize = pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
	  
		  if (layer(i_level_py).ne.layer(i_level_py-1)) then
		  
	        HI_column_density_in = 0
	        HeI_column_density_in = 0
	        HeII_column_density_in = 0
		  		  
	        do parent_level_py = 1,i_level_py-1
		
	          parent_primary = parent_position(i_level_py,i_primary)%with_parent_height(parent_level_py)
	          parent_secondary = parent_position(i_level_py,i_secondary)%with_parent_height(parent_level_py)
		
	          HI_column_density_in = HI_column_density_in + &
	                            px_ny_pz_dx(parent_level_py)%HI_density(-parent_primary,parent_secondary) * cellsize
			              
	          HeI_column_density_in = HeI_column_density_in + &
	                            px_ny_pz_dx(parent_level_py)%HeI_density(-parent_primary,parent_secondary) * cellsize
						  							  
	          HeII_column_density_in = HeII_column_density_in + &
	                            px_ny_pz_dx(parent_level_py)%HeII_density(-parent_primary,parent_secondary) * cellsize 
						  
	        enddo	
	  
	        px_ny_pz_dx(i_level_py)%HI_column_density_in(-i_primary,i_secondary) = HI_column_density_in	
	        px_ny_pz_dx(i_level_py)%HeI_column_density_in(-i_primary,i_secondary) = HeI_column_density_in
	        px_ny_pz_dx(i_level_py)%HeII_column_density_in(-i_primary,i_secondary) = HeII_column_density_in		
	
	      else
		  
	        px_ny_pz_dx(i_level_py)%HI_column_density_in(-i_primary,i_secondary) = &
					  px_ny_pz_dx(i_level_py-1)%HI_column_density_out(-i_primary,i_secondary)
	        px_ny_pz_dx(i_level_py)%HeI_column_density_in(-i_primary,i_secondary) = &
					  px_ny_pz_dx(i_level_py-1)%HeI_column_density_out(-i_primary,i_secondary)		
	        px_ny_pz_dx(i_level_py)%HeII_column_density_in(-i_primary,i_secondary) = &
					  px_ny_pz_dx(i_level_py-1)%HeII_column_density_out(-i_primary,i_secondary)
		  
	      endif
	  			    			
		  px_ny_pz_dx(i_level_py)%HI_column_density_out(-i_primary,i_secondary) = &
	                       px_ny_pz_dx(i_level_py)%HI_column_density_in(-i_primary,i_secondary) + &
	                       px_ny_pz_dx(i_level_py)%HI_density(-i_primary,i_secondary) * cellsize
		  px_ny_pz_dx(i_level_py)%HeI_column_density_out(-i_primary,i_secondary) = &
	                       px_ny_pz_dx(i_level_py)%HeI_column_density_in(-i_primary,i_secondary) + &
	                       px_ny_pz_dx(i_level_py)%HeI_density(-i_primary,i_secondary) * cellsize	  
		  px_ny_pz_dx(i_level_py)%HeII_column_density_out(-i_primary,i_secondary) = &
	                       px_ny_pz_dx(i_level_py)%HeII_column_density_in(-i_primary,i_secondary) + &
	                       px_ny_pz_dx(i_level_py)%HeII_density(-i_primary,i_secondary) * cellsize		  
	  	
		enddo	  
	  enddo	  
	enddo
	  
  end subroutine px_ny_pz_dx_column_density_calculation
  
  subroutine px_ny_pz_dy_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	integer :: parent_primary,parent_secondary,parent_level_py
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_cell
	real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_cell
    real(kind=dp) :: cellsize

	! i_level_py = 1
	px_ny_pz_dy(-1)%HI_column_density_in(1,1) = 0
	px_ny_pz_dy(-1)%HeI_column_density_in(1,1) = 0
	px_ny_pz_dy(-1)%HeII_column_density_in(1,1) = 0

	cellsize = pyramid_grid(1)%cellsize(1,1)

	px_ny_pz_dy(-1)%HI_column_density_out(1,1) = px_ny_pz_dy(-1)%HI_density(1,1) * cellsize	
	px_ny_pz_dy(-1)%HeI_column_density_out(1,1) = px_ny_pz_dy(-1)%HeI_density(1,1) * cellsize
	px_ny_pz_dy(-1)%HeII_column_density_out(1,1) = px_ny_pz_dy(-1)%HeII_density(1,1) * cellsize
			
	! i_level_py > 1
	do i_level_py = 2,level_py
	  do i_primary = 1,partition(i_level_py)
	    do i_secondary = 1,partition(i_level_py)
		
		  cellsize = pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
	  
		  if (layer(i_level_py).ne.layer(i_level_py-1)) then
		  
	        HI_column_density_in = 0
	        HeI_column_density_in = 0
	        HeII_column_density_in = 0
		  		  
	        do parent_level_py = 1,i_level_py-1
		
	          parent_primary = parent_position(i_level_py,i_primary)%with_parent_height(parent_level_py)
	          parent_secondary = parent_position(i_level_py,i_secondary)%with_parent_height(parent_level_py)
		
	          HI_column_density_in = HI_column_density_in + &
	                            px_ny_pz_dy(-parent_level_py)%HI_density(parent_primary,parent_secondary) * cellsize
			              
	          HeI_column_density_in = HeI_column_density_in + &
	                            px_ny_pz_dy(-parent_level_py)%HeI_density(parent_primary,parent_secondary) * cellsize
						  							  
	          HeII_column_density_in = HeII_column_density_in + &
	                            px_ny_pz_dy(-parent_level_py)%HeII_density(parent_primary,parent_secondary) * cellsize 
						  
	        enddo	
	  
	        px_ny_pz_dy(-i_level_py)%HI_column_density_in(i_primary,i_secondary) = HI_column_density_in	
	        px_ny_pz_dy(-i_level_py)%HeI_column_density_in(i_primary,i_secondary) = HeI_column_density_in
	        px_ny_pz_dy(-i_level_py)%HeII_column_density_in(i_primary,i_secondary) = HeII_column_density_in		
	
	      else
		  
	        px_ny_pz_dy(-i_level_py)%HI_column_density_in(i_primary,i_secondary) = &
					  px_ny_pz_dy(-i_level_py+1)%HI_column_density_out(i_primary,i_secondary)
	        px_ny_pz_dy(-i_level_py)%HeI_column_density_in(i_primary,i_secondary) = &
					  px_ny_pz_dy(-i_level_py+1)%HeI_column_density_out(i_primary,i_secondary)		
	        px_ny_pz_dy(-i_level_py)%HeII_column_density_in(i_primary,i_secondary) = &
					  px_ny_pz_dy(-i_level_py+1)%HeII_column_density_out(i_primary,i_secondary)
		  
	      endif
	  			    			
		  px_ny_pz_dy(-i_level_py)%HI_column_density_out(i_primary,i_secondary) = &
	                       px_ny_pz_dy(-i_level_py)%HI_column_density_in(i_primary,i_secondary) + &
	                       px_ny_pz_dy(-i_level_py)%HI_density(i_primary,i_secondary) * cellsize
		  px_ny_pz_dy(-i_level_py)%HeI_column_density_out(i_primary,i_secondary) = &
	                       px_ny_pz_dy(-i_level_py)%HeI_column_density_in(i_primary,i_secondary) + &
	                       px_ny_pz_dy(-i_level_py)%HeI_density(i_primary,i_secondary) * cellsize	  
		  px_ny_pz_dy(-i_level_py)%HeII_column_density_out(i_primary,i_secondary) = &
	                       px_ny_pz_dy(-i_level_py)%HeII_column_density_in(i_primary,i_secondary) + &
	                       px_ny_pz_dy(-i_level_py)%HeII_density(i_primary,i_secondary) * cellsize		  
	  	
		enddo	  
	  enddo	  
	enddo
	  
  end subroutine px_ny_pz_dy_column_density_calculation
  
  subroutine px_ny_pz_dz_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	integer :: parent_primary,parent_secondary,parent_level_py
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_cell
	real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_cell
    real(kind=dp) :: cellsize

	! i_level_py = 1
	px_ny_pz_dz(1)%HI_column_density_in(1,-1) = 0
	px_ny_pz_dz(1)%HeI_column_density_in(1,-1) = 0
	px_ny_pz_dz(1)%HeII_column_density_in(1,-1) = 0

	cellsize = pyramid_grid(1)%cellsize(1,1)

	px_ny_pz_dz(1)%HI_column_density_out(1,-1) = px_ny_pz_dz(1)%HI_density(1,-1) * cellsize	
	px_ny_pz_dz(1)%HeI_column_density_out(1,-1) = px_ny_pz_dz(1)%HeI_density(1,-1) * cellsize
	px_ny_pz_dz(1)%HeII_column_density_out(1,-1) = px_ny_pz_dz(1)%HeII_density(1,-1) * cellsize
		
	! i_level_py > 1
	do i_level_py = 2,level_py
	  do i_primary = 1,partition(i_level_py)
	    do i_secondary = 1,partition(i_level_py)
		
		  cellsize = pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
	  
		  if (layer(i_level_py).ne.layer(i_level_py-1)) then
		  
	        HI_column_density_in = 0
	        HeI_column_density_in = 0
	        HeII_column_density_in = 0
		  		  
	        do parent_level_py = 1,i_level_py-1
		
	          parent_primary = parent_position(i_level_py,i_primary)%with_parent_height(parent_level_py)
	          parent_secondary = parent_position(i_level_py,i_secondary)%with_parent_height(parent_level_py)
		
	          HI_column_density_in = HI_column_density_in + &
	                            px_ny_pz_dz(parent_level_py)%HI_density(parent_primary,-parent_secondary) * cellsize
			              
	          HeI_column_density_in = HeI_column_density_in + &
	                            px_ny_pz_dz(parent_level_py)%HeI_density(parent_primary,-parent_secondary) * cellsize
						  							  
	          HeII_column_density_in = HeII_column_density_in + &
	                            px_ny_pz_dz(parent_level_py)%HeII_density(parent_primary,-parent_secondary) * cellsize 
						  
	        enddo	
	  
	        px_ny_pz_dz(i_level_py)%HI_column_density_in(i_primary,-i_secondary) = HI_column_density_in	
	        px_ny_pz_dz(i_level_py)%HeI_column_density_in(i_primary,-i_secondary) = HeI_column_density_in
	        px_ny_pz_dz(i_level_py)%HeII_column_density_in(i_primary,-i_secondary) = HeII_column_density_in		
	
	      else
		  
	        px_ny_pz_dz(i_level_py)%HI_column_density_in(i_primary,-i_secondary) = &
					  px_ny_pz_dz(i_level_py-1)%HI_column_density_out(i_primary,-i_secondary)
	        px_ny_pz_dz(i_level_py)%HeI_column_density_in(i_primary,-i_secondary) = &
					  px_ny_pz_dz(i_level_py-1)%HeI_column_density_out(i_primary,-i_secondary)		
	        px_ny_pz_dz(i_level_py)%HeII_column_density_in(i_primary,-i_secondary) = &
					  px_ny_pz_dz(i_level_py-1)%HeII_column_density_out(i_primary,-i_secondary)
		  
	      endif
	  			    			
		  px_ny_pz_dz(i_level_py)%HI_column_density_out(i_primary,-i_secondary) = &
	                       px_ny_pz_dz(i_level_py)%HI_column_density_in(i_primary,-i_secondary) + &
	                       px_ny_pz_dz(i_level_py)%HI_density(i_primary,-i_secondary) * cellsize
		  px_ny_pz_dz(i_level_py)%HeI_column_density_out(i_primary,-i_secondary) = &
	                       px_ny_pz_dz(i_level_py)%HeI_column_density_in(i_primary,-i_secondary) + &
	                       px_ny_pz_dz(i_level_py)%HeI_density(i_primary,-i_secondary) * cellsize	  
		  px_ny_pz_dz(i_level_py)%HeII_column_density_out(i_primary,-i_secondary) = &
	                       px_ny_pz_dz(i_level_py)%HeII_column_density_in(i_primary,-i_secondary) + &
	                       px_ny_pz_dz(i_level_py)%HeII_density(i_primary,-i_secondary) * cellsize		  
	  	
		enddo	  
	  enddo	  
	enddo  
	  
  end subroutine px_ny_pz_dz_column_density_calculation
     
  subroutine px_ny_nz_dx_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	integer :: parent_primary,parent_secondary,parent_level_py
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_cell
	real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_cell
    real(kind=dp) :: cellsize

	! i_level_py = 1
	px_ny_nz_dx(1)%HI_column_density_in(-1,-1) = 0
	px_ny_nz_dx(1)%HeI_column_density_in(-1,-1) = 0
	px_ny_nz_dx(1)%HeII_column_density_in(-1,-1) = 0

	cellsize = pyramid_grid(1)%cellsize(1,1)

	px_ny_nz_dx(1)%HI_column_density_out(-1,-1) = px_ny_nz_dx(1)%HI_density(-1,-1) * cellsize	
	px_ny_nz_dx(1)%HeI_column_density_out(-1,-1) = px_ny_nz_dx(1)%HeI_density(-1,-1) * cellsize
	px_ny_nz_dx(1)%HeII_column_density_out(-1,-1) = px_ny_nz_dx(1)%HeII_density(-1,-1) * cellsize
		
	! i_level_py > 1
	do i_level_py = 2,level_py
	  do i_primary = 1,partition(i_level_py)
	    do i_secondary = 1,partition(i_level_py)
		
		  cellsize = pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
	  
		  if (layer(i_level_py).ne.layer(i_level_py-1)) then
		  
	        HI_column_density_in = 0
	        HeI_column_density_in = 0
	        HeII_column_density_in = 0
		  		  
	        do parent_level_py = 1,i_level_py-1
		
	          parent_primary = parent_position(i_level_py,i_primary)%with_parent_height(parent_level_py)
	          parent_secondary = parent_position(i_level_py,i_secondary)%with_parent_height(parent_level_py)
		
	          HI_column_density_in = HI_column_density_in + &
	                            px_ny_nz_dx(parent_level_py)%HI_density(-parent_primary,-parent_secondary) * cellsize
			              
	          HeI_column_density_in = HeI_column_density_in + &
	                            px_ny_nz_dx(parent_level_py)%HeI_density(-parent_primary,-parent_secondary) * cellsize
						  							  
	          HeII_column_density_in = HeII_column_density_in + &
	                            px_ny_nz_dx(parent_level_py)%HeII_density(-parent_primary,-parent_secondary) * cellsize 
						  
	        enddo	
	  
	        px_ny_nz_dx(i_level_py)%HI_column_density_in(-i_primary,-i_secondary) = HI_column_density_in	
	        px_ny_nz_dx(i_level_py)%HeI_column_density_in(-i_primary,-i_secondary) = HeI_column_density_in
	        px_ny_nz_dx(i_level_py)%HeII_column_density_in(-i_primary,-i_secondary) = HeII_column_density_in		
	
	      else
		  
	        px_ny_nz_dx(i_level_py)%HI_column_density_in(-i_primary,-i_secondary) = &
					  px_ny_nz_dx(i_level_py-1)%HI_column_density_out(-i_primary,-i_secondary)
	        px_ny_nz_dx(i_level_py)%HeI_column_density_in(-i_primary,-i_secondary) = &
					  px_ny_nz_dx(i_level_py-1)%HeI_column_density_out(-i_primary,-i_secondary)		
	        px_ny_nz_dx(i_level_py)%HeII_column_density_in(-i_primary,-i_secondary) = &
					  px_ny_nz_dx(i_level_py-1)%HeII_column_density_out(-i_primary,-i_secondary)
		  
	      endif
	  			    			
		  px_ny_nz_dx(i_level_py)%HI_column_density_out(-i_primary,-i_secondary) = &
	                       px_ny_nz_dx(i_level_py)%HI_column_density_in(-i_primary,-i_secondary) + &
	                       px_ny_nz_dx(i_level_py)%HI_density(-i_primary,-i_secondary) * cellsize
		  px_ny_nz_dx(i_level_py)%HeI_column_density_out(-i_primary,-i_secondary) = &
	                       px_ny_nz_dx(i_level_py)%HeI_column_density_in(-i_primary,-i_secondary) + &
	                       px_ny_nz_dx(i_level_py)%HeI_density(-i_primary,-i_secondary) * cellsize	  
		  px_ny_nz_dx(i_level_py)%HeII_column_density_out(-i_primary,-i_secondary) = &
	                       px_ny_nz_dx(i_level_py)%HeII_column_density_in(-i_primary,-i_secondary) + &
	                       px_ny_nz_dx(i_level_py)%HeII_density(-i_primary,-i_secondary) * cellsize		  
	  	
		enddo	  
	  enddo	  
	enddo
	  
  end subroutine px_ny_nz_dx_column_density_calculation
  
  subroutine px_ny_nz_dy_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	integer :: parent_primary,parent_secondary,parent_level_py
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_cell
	real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_cell
    real(kind=dp) :: cellsize
	
	! i_level_py = 1
	px_ny_nz_dy(-1)%HI_column_density_in(-1,1) = 0
	px_ny_nz_dy(-1)%HeI_column_density_in(-1,1) = 0
	px_ny_nz_dy(-1)%HeII_column_density_in(-1,1) = 0

	cellsize = pyramid_grid(1)%cellsize(1,1)

	px_ny_nz_dy(-1)%HI_column_density_out(-1,1) = px_ny_nz_dy(-1)%HI_density(-1,1) * cellsize	
	px_ny_nz_dy(-1)%HeI_column_density_out(-1,1) = px_ny_nz_dy(-1)%HeI_density(-1,1) * cellsize
	px_ny_nz_dy(-1)%HeII_column_density_out(-1,1) = px_ny_nz_dy(-1)%HeII_density(-1,1) * cellsize
	
	! i_level_py > 1
	do i_level_py = 2,level_py
	  do i_primary = 1,partition(i_level_py)
	    do i_secondary = 1,partition(i_level_py)
		
		  cellsize = pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
	  
		  if (layer(i_level_py).ne.layer(i_level_py-1)) then
		  
	        HI_column_density_in = 0
	        HeI_column_density_in = 0
	        HeII_column_density_in = 0
		  		  
	        do parent_level_py = 1,i_level_py-1
		
	          parent_primary = parent_position(i_level_py,i_primary)%with_parent_height(parent_level_py)
	          parent_secondary = parent_position(i_level_py,i_secondary)%with_parent_height(parent_level_py)
		
	          HI_column_density_in = HI_column_density_in + &
	                            px_ny_nz_dy(-parent_level_py)%HI_density(-parent_primary,parent_secondary) * cellsize
			              
	          HeI_column_density_in = HeI_column_density_in + &
	                            px_ny_nz_dy(-parent_level_py)%HeI_density(-parent_primary,parent_secondary) * cellsize
						  							  
	          HeII_column_density_in = HeII_column_density_in + &
	                            px_ny_nz_dy(-parent_level_py)%HeII_density(-parent_primary,parent_secondary) * cellsize 
						  
	        enddo	
	  
	        px_ny_nz_dy(-i_level_py)%HI_column_density_in(-i_primary,i_secondary) = HI_column_density_in	
	        px_ny_nz_dy(-i_level_py)%HeI_column_density_in(-i_primary,i_secondary) = HeI_column_density_in
	        px_ny_nz_dy(-i_level_py)%HeII_column_density_in(-i_primary,i_secondary) = HeII_column_density_in		
	
	      else
		  
	        px_ny_nz_dy(-i_level_py)%HI_column_density_in(-i_primary,i_secondary) = &
					  px_ny_nz_dy(-i_level_py+1)%HI_column_density_out(-i_primary,i_secondary)
	        px_ny_nz_dy(-i_level_py)%HeI_column_density_in(-i_primary,i_secondary) = &
					  px_ny_nz_dy(-i_level_py+1)%HeI_column_density_out(-i_primary,i_secondary)		
	        px_ny_nz_dy(-i_level_py)%HeII_column_density_in(-i_primary,i_secondary) = &
					  px_ny_nz_dy(-i_level_py+1)%HeII_column_density_out(-i_primary,i_secondary)
		  
	      endif
	  			    			
		  px_ny_nz_dy(-i_level_py)%HI_column_density_out(-i_primary,i_secondary) = &
	                       px_ny_nz_dy(-i_level_py)%HI_column_density_in(-i_primary,i_secondary) + &
	                       px_ny_nz_dy(-i_level_py)%HI_density(-i_primary,i_secondary) * cellsize
		  px_ny_nz_dy(-i_level_py)%HeI_column_density_out(-i_primary,i_secondary) = &
	                       px_ny_nz_dy(-i_level_py)%HeI_column_density_in(-i_primary,i_secondary) + &
	                       px_ny_nz_dy(-i_level_py)%HeI_density(-i_primary,i_secondary) * cellsize	  
		  px_ny_nz_dy(-i_level_py)%HeII_column_density_out(-i_primary,i_secondary) = &
	                       px_ny_nz_dy(-i_level_py)%HeII_column_density_in(-i_primary,i_secondary) + &
	                       px_ny_nz_dy(-i_level_py)%HeII_density(-i_primary,i_secondary) * cellsize		  
	  	
		enddo	  
	  enddo	  
	enddo 
	  
  end subroutine px_ny_nz_dy_column_density_calculation
  
  subroutine px_ny_nz_dz_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	integer :: parent_primary,parent_secondary,parent_level_py
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_cell
	real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_cell
    real(kind=dp) :: cellsize

	! i_level_py = 1
	px_ny_nz_dz(-1)%HI_column_density_in(1,-1) = 0
	px_ny_nz_dz(-1)%HeI_column_density_in(1,-1) = 0
	px_ny_nz_dz(-1)%HeII_column_density_in(1,-1) = 0

	cellsize = pyramid_grid(1)%cellsize(1,1)

	px_ny_nz_dz(-1)%HI_column_density_out(1,-1) = px_ny_nz_dz(-1)%HI_density(1,-1) * cellsize	
	px_ny_nz_dz(-1)%HeI_column_density_out(1,-1) = px_ny_nz_dz(-1)%HeI_density(1,-1) * cellsize
	px_ny_nz_dz(-1)%HeII_column_density_out(1,-1) = px_ny_nz_dz(-1)%HeII_density(1,-1) * cellsize
		
	! i_level_py > 1
	do i_level_py = 2,level_py
	  do i_primary = 1,partition(i_level_py)
	    do i_secondary = 1,partition(i_level_py)
		
		  cellsize = pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
	  
		  if (layer(i_level_py).ne.layer(i_level_py-1)) then
		  
	        HI_column_density_in = 0
	        HeI_column_density_in = 0
	        HeII_column_density_in = 0
		  		  
	        do parent_level_py = 1,i_level_py-1
		
	          parent_primary = parent_position(i_level_py,i_primary)%with_parent_height(parent_level_py)
	          parent_secondary = parent_position(i_level_py,i_secondary)%with_parent_height(parent_level_py)
		
	          HI_column_density_in = HI_column_density_in + &
	                            px_ny_nz_dz(-parent_level_py)%HI_density(parent_primary,-parent_secondary) * cellsize
			              
	          HeI_column_density_in = HeI_column_density_in + &
	                            px_ny_nz_dz(-parent_level_py)%HeI_density(parent_primary,-parent_secondary) * cellsize
						  							  
	          HeII_column_density_in = HeII_column_density_in + &
	                            px_ny_nz_dz(-parent_level_py)%HeII_density(parent_primary,-parent_secondary) * cellsize 
						  
	        enddo	
	  
	        px_ny_nz_dz(-i_level_py)%HI_column_density_in(i_primary,-i_secondary) = HI_column_density_in	
	        px_ny_nz_dz(-i_level_py)%HeI_column_density_in(i_primary,-i_secondary) = HeI_column_density_in
	        px_ny_nz_dz(-i_level_py)%HeII_column_density_in(i_primary,-i_secondary) = HeII_column_density_in		
	
	      else
		  
	        px_ny_nz_dz(-i_level_py)%HI_column_density_in(i_primary,-i_secondary) = &
					  px_ny_nz_dz(-i_level_py+1)%HI_column_density_out(i_primary,-i_secondary)
	        px_ny_nz_dz(-i_level_py)%HeI_column_density_in(i_primary,-i_secondary) = &
					  px_ny_nz_dz(-i_level_py+1)%HeI_column_density_out(i_primary,-i_secondary)		
	        px_ny_nz_dz(-i_level_py)%HeII_column_density_in(i_primary,-i_secondary) = &
					  px_ny_nz_dz(-i_level_py+1)%HeII_column_density_out(i_primary,-i_secondary)
		  
	      endif
	  			    			
		  px_ny_nz_dz(-i_level_py)%HI_column_density_out(i_primary,-i_secondary) = &
	                       px_ny_nz_dz(-i_level_py)%HI_column_density_in(i_primary,-i_secondary) + &
	                       px_ny_nz_dz(-i_level_py)%HI_density(i_primary,-i_secondary) * cellsize
		  px_ny_nz_dz(-i_level_py)%HeI_column_density_out(i_primary,-i_secondary) = &
	                       px_ny_nz_dz(-i_level_py)%HeI_column_density_in(i_primary,-i_secondary) + &
	                       px_ny_nz_dz(-i_level_py)%HeI_density(i_primary,-i_secondary) * cellsize	  
		  px_ny_nz_dz(-i_level_py)%HeII_column_density_out(i_primary,-i_secondary) = &
	                       px_ny_nz_dz(-i_level_py)%HeII_column_density_in(i_primary,-i_secondary) + &
	                       px_ny_nz_dz(-i_level_py)%HeII_density(i_primary,-i_secondary) * cellsize		  
	  	
		enddo	  
	  enddo	  
	enddo
	  
  end subroutine px_ny_nz_dz_column_density_calculation

  subroutine nx_py_pz_dx_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	integer :: parent_primary,parent_secondary,parent_level_py
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_cell
	real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_cell
    real(kind=dp) :: cellsize

	! i_level_py = 1
	nx_py_pz_dx(-1)%HI_column_density_in(1,1) = 0
	nx_py_pz_dx(-1)%HeI_column_density_in(1,1) = 0
	nx_py_pz_dx(-1)%HeII_column_density_in(1,1) = 0

	cellsize = pyramid_grid(1)%cellsize(1,1)

	nx_py_pz_dx(-1)%HI_column_density_out(1,1) = nx_py_pz_dx(-1)%HI_density(1,1) * cellsize	
	nx_py_pz_dx(-1)%HeI_column_density_out(1,1) = nx_py_pz_dx(-1)%HeI_density(1,1) * cellsize
	nx_py_pz_dx(-1)%HeII_column_density_out(1,1) = nx_py_pz_dx(-1)%HeII_density(1,1) * cellsize
		
	! i_level_py > 1
	do i_level_py = 2,level_py
	  do i_primary = 1,partition(i_level_py)
	    do i_secondary = 1,partition(i_level_py)
		
		  cellsize = pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
	  
		  if (layer(i_level_py).ne.layer(i_level_py-1)) then
		  
	        HI_column_density_in = 0
	        HeI_column_density_in = 0
	        HeII_column_density_in = 0
		  		  
	        do parent_level_py = 1,i_level_py-1
		
	          parent_primary = parent_position(i_level_py,i_primary)%with_parent_height(parent_level_py)
	          parent_secondary = parent_position(i_level_py,i_secondary)%with_parent_height(parent_level_py)
		
	          HI_column_density_in = HI_column_density_in + &
	                            nx_py_pz_dx(-parent_level_py)%HI_density(parent_primary,parent_secondary) * cellsize
			              
	          HeI_column_density_in = HeI_column_density_in + &
	                            nx_py_pz_dx(-parent_level_py)%HeI_density(parent_primary,parent_secondary) * cellsize
						  							  
	          HeII_column_density_in = HeII_column_density_in + &
	                            nx_py_pz_dx(-parent_level_py)%HeII_density(parent_primary,parent_secondary) * cellsize 
						  
	        enddo	
	  
	        nx_py_pz_dx(-i_level_py)%HI_column_density_in(i_primary,i_secondary) = HI_column_density_in	
	        nx_py_pz_dx(-i_level_py)%HeI_column_density_in(i_primary,i_secondary) = HeI_column_density_in
	        nx_py_pz_dx(-i_level_py)%HeII_column_density_in(i_primary,i_secondary) = HeII_column_density_in		
	
	      else
		  
	        nx_py_pz_dx(-i_level_py)%HI_column_density_in(i_primary,i_secondary) = &
					  nx_py_pz_dx(-i_level_py+1)%HI_column_density_out(i_primary,i_secondary)
	        nx_py_pz_dx(-i_level_py)%HeI_column_density_in(i_primary,i_secondary) = &
					  nx_py_pz_dx(-i_level_py+1)%HeI_column_density_out(i_primary,i_secondary)		
	        nx_py_pz_dx(-i_level_py)%HeII_column_density_in(i_primary,i_secondary) = &
					  nx_py_pz_dx(-i_level_py+1)%HeII_column_density_out(i_primary,i_secondary)
		  
	      endif
	  			    			
		  nx_py_pz_dx(-i_level_py)%HI_column_density_out(i_primary,i_secondary) = &
	                       nx_py_pz_dx(-i_level_py)%HI_column_density_in(i_primary,i_secondary) + &
	                       nx_py_pz_dx(-i_level_py)%HI_density(i_primary,i_secondary) * cellsize
		  nx_py_pz_dx(-i_level_py)%HeI_column_density_out(i_primary,i_secondary) = &
	                       nx_py_pz_dx(-i_level_py)%HeI_column_density_in(i_primary,i_secondary) + &
	                       nx_py_pz_dx(-i_level_py)%HeI_density(i_primary,i_secondary) * cellsize	  
		  nx_py_pz_dx(-i_level_py)%HeII_column_density_out(i_primary,i_secondary) = &
	                       nx_py_pz_dx(-i_level_py)%HeII_column_density_in(i_primary,i_secondary) + &
	                       nx_py_pz_dx(-i_level_py)%HeII_density(i_primary,i_secondary) * cellsize		  
	  	
		enddo	  
	  enddo	  
	enddo
	  
  end subroutine nx_py_pz_dx_column_density_calculation
  
  subroutine nx_py_pz_dy_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	integer :: parent_primary,parent_secondary,parent_level_py
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_cell
	real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_cell
    real(kind=dp) :: cellsize

	! i_level_py = 1
	nx_py_pz_dy(1)%HI_column_density_in(1,-1) = 0
	nx_py_pz_dy(1)%HeI_column_density_in(1,-1) = 0
	nx_py_pz_dy(1)%HeII_column_density_in(1,-1) = 0

	cellsize = pyramid_grid(1)%cellsize(1,1)

	nx_py_pz_dy(1)%HI_column_density_out(1,-1) = nx_py_pz_dy(1)%HI_density(1,-1) * cellsize	
	nx_py_pz_dy(1)%HeI_column_density_out(1,-1) = nx_py_pz_dy(1)%HeI_density(1,-1) * cellsize
	nx_py_pz_dy(1)%HeII_column_density_out(1,-1) = nx_py_pz_dy(1)%HeII_density(1,-1) * cellsize
		
	! i_level_py > 1
	do i_level_py = 2,level_py
	  do i_primary = 1,partition(i_level_py)
	    do i_secondary = 1,partition(i_level_py)
		
		  cellsize = pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
	  
		  if (layer(i_level_py).ne.layer(i_level_py-1)) then
		  
	        HI_column_density_in = 0
	        HeI_column_density_in = 0
	        HeII_column_density_in = 0
		  		  
	        do parent_level_py = 1,i_level_py-1
		
	          parent_primary = parent_position(i_level_py,i_primary)%with_parent_height(parent_level_py)
	          parent_secondary = parent_position(i_level_py,i_secondary)%with_parent_height(parent_level_py)
		
	          HI_column_density_in = HI_column_density_in + &
	                            nx_py_pz_dy(parent_level_py)%HI_density(parent_primary,-parent_secondary) * cellsize
			              
	          HeI_column_density_in = HeI_column_density_in + &
	                            nx_py_pz_dy(parent_level_py)%HeI_density(parent_primary,-parent_secondary) * cellsize
						  							  
	          HeII_column_density_in = HeII_column_density_in + &
	                            nx_py_pz_dy(parent_level_py)%HeII_density(parent_primary,-parent_secondary) * cellsize 
						  
	        enddo	
	  
	        nx_py_pz_dy(i_level_py)%HI_column_density_in(i_primary,-i_secondary) = HI_column_density_in	
	        nx_py_pz_dy(i_level_py)%HeI_column_density_in(i_primary,-i_secondary) = HeI_column_density_in
	        nx_py_pz_dy(i_level_py)%HeII_column_density_in(i_primary,-i_secondary) = HeII_column_density_in		
	
	      else
		  
	        nx_py_pz_dy(i_level_py)%HI_column_density_in(i_primary,-i_secondary) = &
					  nx_py_pz_dy(i_level_py-1)%HI_column_density_out(i_primary,-i_secondary)
	        nx_py_pz_dy(i_level_py)%HeI_column_density_in(i_primary,-i_secondary) = &
					  nx_py_pz_dy(i_level_py-1)%HeI_column_density_out(i_primary,-i_secondary)		
	        nx_py_pz_dy(i_level_py)%HeII_column_density_in(i_primary,-i_secondary) = &
					  nx_py_pz_dy(i_level_py-1)%HeII_column_density_out(i_primary,-i_secondary)
		  
	      endif
	  			    			
		  nx_py_pz_dy(i_level_py)%HI_column_density_out(i_primary,-i_secondary) = &
	                       nx_py_pz_dy(i_level_py)%HI_column_density_in(i_primary,-i_secondary) + &
	                       nx_py_pz_dy(i_level_py)%HI_density(i_primary,-i_secondary) * cellsize
		  nx_py_pz_dy(i_level_py)%HeI_column_density_out(i_primary,-i_secondary) = &
	                       nx_py_pz_dy(i_level_py)%HeI_column_density_in(i_primary,-i_secondary) + &
	                       nx_py_pz_dy(i_level_py)%HeI_density(i_primary,-i_secondary) * cellsize	  
		  nx_py_pz_dy(i_level_py)%HeII_column_density_out(i_primary,-i_secondary) = &
	                       nx_py_pz_dy(i_level_py)%HeII_column_density_in(i_primary,-i_secondary) + &
	                       nx_py_pz_dy(i_level_py)%HeII_density(i_primary,-i_secondary) * cellsize		  
	  	
		enddo	  
	  enddo	  
	enddo
	  
  end subroutine nx_py_pz_dy_column_density_calculation
  
  subroutine nx_py_pz_dz_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	integer :: parent_primary,parent_secondary,parent_level_py
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_cell
	real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_cell
    real(kind=dp) :: cellsize

	! i_level_py = 1
	nx_py_pz_dz(1)%HI_column_density_in(-1,1) = 0
	nx_py_pz_dz(1)%HeI_column_density_in(-1,1) = 0
	nx_py_pz_dz(1)%HeII_column_density_in(-1,1) = 0

	cellsize = pyramid_grid(1)%cellsize(1,1)

	nx_py_pz_dz(1)%HI_column_density_out(-1,1) = nx_py_pz_dz(1)%HI_density(-1,1) * cellsize	
	nx_py_pz_dz(1)%HeI_column_density_out(-1,1) = nx_py_pz_dz(1)%HeI_density(-1,1) * cellsize
	nx_py_pz_dz(1)%HeII_column_density_out(-1,1) = nx_py_pz_dz(1)%HeII_density(-1,1) * cellsize
		
	! i_level_py > 1
	do i_level_py = 2,level_py
	  do i_primary = 1,partition(i_level_py)
	    do i_secondary = 1,partition(i_level_py)
		
		  cellsize = pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
	  
		  if (layer(i_level_py).ne.layer(i_level_py-1)) then
		  
	        HI_column_density_in = 0
	        HeI_column_density_in = 0
	        HeII_column_density_in = 0
		  		  
	        do parent_level_py = 1,i_level_py-1
		
	          parent_primary = parent_position(i_level_py,i_primary)%with_parent_height(parent_level_py)
	          parent_secondary = parent_position(i_level_py,i_secondary)%with_parent_height(parent_level_py)
		
	          HI_column_density_in = HI_column_density_in + &
	                            nx_py_pz_dz(parent_level_py)%HI_density(-parent_primary,parent_secondary) * cellsize
			              
	          HeI_column_density_in = HeI_column_density_in + &
	                            nx_py_pz_dz(parent_level_py)%HeI_density(-parent_primary,parent_secondary) * cellsize
						  							  
	          HeII_column_density_in = HeII_column_density_in + &
	                            nx_py_pz_dz(parent_level_py)%HeII_density(-parent_primary,parent_secondary) * cellsize 
						  
	        enddo	
	  
	        nx_py_pz_dz(i_level_py)%HI_column_density_in(-i_primary,i_secondary) = HI_column_density_in	
	        nx_py_pz_dz(i_level_py)%HeI_column_density_in(-i_primary,i_secondary) = HeI_column_density_in
	        nx_py_pz_dz(i_level_py)%HeII_column_density_in(-i_primary,i_secondary) = HeII_column_density_in		
	
	      else
		  
	        nx_py_pz_dz(i_level_py)%HI_column_density_in(-i_primary,i_secondary) = &
					  nx_py_pz_dz(i_level_py-1)%HI_column_density_out(-i_primary,i_secondary)
	        nx_py_pz_dz(i_level_py)%HeI_column_density_in(-i_primary,i_secondary) = &
					  nx_py_pz_dz(i_level_py-1)%HeI_column_density_out(-i_primary,i_secondary)		
	        nx_py_pz_dz(i_level_py)%HeII_column_density_in(-i_primary,i_secondary) = &
					  nx_py_pz_dz(i_level_py-1)%HeII_column_density_out(-i_primary,i_secondary)
		  
	      endif
	  			    			
		  nx_py_pz_dz(i_level_py)%HI_column_density_out(-i_primary,i_secondary) = &
	                       nx_py_pz_dz(i_level_py)%HI_column_density_in(-i_primary,i_secondary) + &
	                       nx_py_pz_dz(i_level_py)%HI_density(-i_primary,i_secondary) * cellsize
		  nx_py_pz_dz(i_level_py)%HeI_column_density_out(-i_primary,i_secondary) = &
	                       nx_py_pz_dz(i_level_py)%HeI_column_density_in(-i_primary,i_secondary) + &
	                       nx_py_pz_dz(i_level_py)%HeI_density(-i_primary,i_secondary) * cellsize	  
		  nx_py_pz_dz(i_level_py)%HeII_column_density_out(-i_primary,i_secondary) = &
	                       nx_py_pz_dz(i_level_py)%HeII_column_density_in(-i_primary,i_secondary) + &
	                       nx_py_pz_dz(i_level_py)%HeII_density(-i_primary,i_secondary) * cellsize		  
	  	
		enddo	  
	  enddo	  
	enddo
	  
  end subroutine nx_py_pz_dz_column_density_calculation
     
  subroutine nx_py_nz_dx_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	integer :: parent_primary,parent_secondary,parent_level_py
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_cell
	real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_cell
    real(kind=dp) :: cellsize

	! i_level_py = 1
	nx_py_nz_dx(-1)%HI_column_density_in(1,-1) = 0
	nx_py_nz_dx(-1)%HeI_column_density_in(1,-1) = 0
	nx_py_nz_dx(-1)%HeII_column_density_in(1,-1) = 0

	cellsize = pyramid_grid(1)%cellsize(1,1)

	nx_py_nz_dx(-1)%HI_column_density_out(1,-1) = nx_py_nz_dx(-1)%HI_density(1,-1) * cellsize	
	nx_py_nz_dx(-1)%HeI_column_density_out(1,-1) = nx_py_nz_dx(-1)%HeI_density(1,-1) * cellsize
	nx_py_nz_dx(-1)%HeII_column_density_out(1,-1) = nx_py_nz_dx(-1)%HeII_density(1,-1) * cellsize
		
	! i_level_py > 1
	do i_level_py = 2,level_py
	  do i_primary = 1,partition(i_level_py)
	    do i_secondary = 1,partition(i_level_py)
		
		  cellsize = pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
	  
		  if (layer(i_level_py).ne.layer(i_level_py-1)) then
		  
	        HI_column_density_in = 0
	        HeI_column_density_in = 0
	        HeII_column_density_in = 0
		  		  
	        do parent_level_py = 1,i_level_py-1
		
	          parent_primary = parent_position(i_level_py,i_primary)%with_parent_height(parent_level_py)
	          parent_secondary = parent_position(i_level_py,i_secondary)%with_parent_height(parent_level_py)
		
	          HI_column_density_in = HI_column_density_in + &
	                            nx_py_nz_dx(-parent_level_py)%HI_density(parent_primary,-parent_secondary) * cellsize
			              
	          HeI_column_density_in = HeI_column_density_in + &
	                            nx_py_nz_dx(-parent_level_py)%HeI_density(parent_primary,-parent_secondary) * cellsize
						  							  
	          HeII_column_density_in = HeII_column_density_in + &
	                            nx_py_nz_dx(-parent_level_py)%HeII_density(parent_primary,-parent_secondary) * cellsize 
						  
	        enddo	
	  
	        nx_py_nz_dx(-i_level_py)%HI_column_density_in(i_primary,-i_secondary) = HI_column_density_in	
	        nx_py_nz_dx(-i_level_py)%HeI_column_density_in(i_primary,-i_secondary) = HeI_column_density_in
	        nx_py_nz_dx(-i_level_py)%HeII_column_density_in(i_primary,-i_secondary) = HeII_column_density_in		
	
	      else
		  
	        nx_py_nz_dx(-i_level_py)%HI_column_density_in(i_primary,-i_secondary) = &
					  nx_py_nz_dx(-i_level_py+1)%HI_column_density_out(i_primary,-i_secondary)
	        nx_py_nz_dx(-i_level_py)%HeI_column_density_in(i_primary,-i_secondary) = &
					  nx_py_nz_dx(-i_level_py+1)%HeI_column_density_out(i_primary,-i_secondary)		
	        nx_py_nz_dx(-i_level_py)%HeII_column_density_in(i_primary,-i_secondary) = &
					  nx_py_nz_dx(-i_level_py+1)%HeII_column_density_out(i_primary,-i_secondary)
		  
	      endif
	  			    			
		  nx_py_nz_dx(-i_level_py)%HI_column_density_out(i_primary,-i_secondary) = &
	                       nx_py_nz_dx(-i_level_py)%HI_column_density_in(i_primary,-i_secondary) + &
	                       nx_py_nz_dx(-i_level_py)%HI_density(i_primary,-i_secondary) * cellsize
		  nx_py_nz_dx(-i_level_py)%HeI_column_density_out(i_primary,-i_secondary) = &
	                       nx_py_nz_dx(-i_level_py)%HeI_column_density_in(i_primary,-i_secondary) + &
	                       nx_py_nz_dx(-i_level_py)%HeI_density(i_primary,-i_secondary) * cellsize	  
		  nx_py_nz_dx(-i_level_py)%HeII_column_density_out(i_primary,-i_secondary) = &
	                       nx_py_nz_dx(-i_level_py)%HeII_column_density_in(i_primary,-i_secondary) + &
	                       nx_py_nz_dx(-i_level_py)%HeII_density(i_primary,-i_secondary) * cellsize		  
	  	
		enddo	  
	  enddo	  
	enddo
	  
  end subroutine nx_py_nz_dx_column_density_calculation
  
  subroutine nx_py_nz_dy_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	integer :: parent_primary,parent_secondary,parent_level_py
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_cell
	real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_cell
    real(kind=dp) :: cellsize

	! i_level_py = 1
	nx_py_nz_dy(1)%HI_column_density_in(-1,-1) = 0
	nx_py_nz_dy(1)%HeI_column_density_in(-1,-1) = 0
	nx_py_nz_dy(1)%HeII_column_density_in(-1,-1) = 0

	cellsize = pyramid_grid(1)%cellsize(1,1)

	nx_py_nz_dy(1)%HI_column_density_out(-1,-1) = nx_py_nz_dy(1)%HI_density(-1,-1) * cellsize	
	nx_py_nz_dy(1)%HeI_column_density_out(-1,-1) = nx_py_nz_dy(1)%HeI_density(-1,-1) * cellsize
	nx_py_nz_dy(1)%HeII_column_density_out(-1,-1) = nx_py_nz_dy(1)%HeII_density(-1,-1) * cellsize
		
	! i_level_py > 1
	do i_level_py = 2,level_py
	  do i_primary = 1,partition(i_level_py)
	    do i_secondary = 1,partition(i_level_py)
		
		  cellsize = pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
	  
		  if (layer(i_level_py).ne.layer(i_level_py-1)) then
		  
	        HI_column_density_in = 0
	        HeI_column_density_in = 0
	        HeII_column_density_in = 0
		  		  
	        do parent_level_py = 1,i_level_py-1
		
	          parent_primary = parent_position(i_level_py,i_primary)%with_parent_height(parent_level_py)
	          parent_secondary = parent_position(i_level_py,i_secondary)%with_parent_height(parent_level_py)
		
	          HI_column_density_in = HI_column_density_in + &
	                            nx_py_nz_dy(parent_level_py)%HI_density(-parent_primary,-parent_secondary) * cellsize
			              
	          HeI_column_density_in = HeI_column_density_in + &
	                            nx_py_nz_dy(parent_level_py)%HeI_density(-parent_primary,-parent_secondary) * cellsize
						  							  
	          HeII_column_density_in = HeII_column_density_in + &
	                            nx_py_nz_dy(parent_level_py)%HeII_density(-parent_primary,-parent_secondary) * cellsize 
						  
	        enddo	
	  
	        nx_py_nz_dy(i_level_py)%HI_column_density_in(-i_primary,-i_secondary) = HI_column_density_in	
	        nx_py_nz_dy(i_level_py)%HeI_column_density_in(-i_primary,-i_secondary) = HeI_column_density_in
	        nx_py_nz_dy(i_level_py)%HeII_column_density_in(-i_primary,-i_secondary) = HeII_column_density_in		
	
	      else
		  
	        nx_py_nz_dy(i_level_py)%HI_column_density_in(-i_primary,-i_secondary) = &
					  nx_py_nz_dy(i_level_py-1)%HI_column_density_out(-i_primary,-i_secondary)
	        nx_py_nz_dy(i_level_py)%HeI_column_density_in(-i_primary,-i_secondary) = &
					  nx_py_nz_dy(i_level_py-1)%HeI_column_density_out(-i_primary,-i_secondary)		
	        nx_py_nz_dy(i_level_py)%HeII_column_density_in(-i_primary,-i_secondary) = &
					  nx_py_nz_dy(i_level_py-1)%HeII_column_density_out(-i_primary,-i_secondary)
		  
	      endif
	  			    			
		  nx_py_nz_dy(i_level_py)%HI_column_density_out(-i_primary,-i_secondary) = &
	                       nx_py_nz_dy(i_level_py)%HI_column_density_in(-i_primary,-i_secondary) + &
	                       nx_py_nz_dy(i_level_py)%HI_density(-i_primary,-i_secondary) * cellsize
		  nx_py_nz_dy(i_level_py)%HeI_column_density_out(-i_primary,-i_secondary) = &
	                       nx_py_nz_dy(i_level_py)%HeI_column_density_in(-i_primary,-i_secondary) + &
	                       nx_py_nz_dy(i_level_py)%HeI_density(-i_primary,-i_secondary) * cellsize	  
		  nx_py_nz_dy(i_level_py)%HeII_column_density_out(-i_primary,-i_secondary) = &
	                       nx_py_nz_dy(i_level_py)%HeII_column_density_in(-i_primary,-i_secondary) + &
	                       nx_py_nz_dy(i_level_py)%HeII_density(-i_primary,-i_secondary) * cellsize		  
	  	
		enddo	  
	  enddo	  
	enddo
	  
  end subroutine nx_py_nz_dy_column_density_calculation
  
  subroutine nx_py_nz_dz_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	integer :: parent_primary,parent_secondary,parent_level_py
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_cell
	real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_cell
    real(kind=dp) :: cellsize
	
	! i_level_py = 1
	nx_py_nz_dz(-1)%HI_column_density_in(-1,1) = 0
	nx_py_nz_dz(-1)%HeI_column_density_in(-1,1) = 0
	nx_py_nz_dz(-1)%HeII_column_density_in(-1,1) = 0

	cellsize = pyramid_grid(1)%cellsize(1,1)

	nx_py_nz_dz(-1)%HI_column_density_out(-1,1) = nx_py_nz_dz(-1)%HI_density(-1,1) * cellsize	
	nx_py_nz_dz(-1)%HeI_column_density_out(-1,1) = nx_py_nz_dz(-1)%HeI_density(-1,1) * cellsize
	nx_py_nz_dz(-1)%HeII_column_density_out(-1,1) = nx_py_nz_dz(-1)%HeII_density(-1,1) * cellsize
		
	! i_level_py > 1
	do i_level_py = 2,level_py
	  do i_primary = 1,partition(i_level_py)
	    do i_secondary = 1,partition(i_level_py)
		
		  cellsize = pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
	  
		  if (layer(i_level_py).ne.layer(i_level_py-1)) then
		  
	        HI_column_density_in = 0
	        HeI_column_density_in = 0
	        HeII_column_density_in = 0
		  		  
	        do parent_level_py = 1,i_level_py-1
		
	          parent_primary = parent_position(i_level_py,i_primary)%with_parent_height(parent_level_py)
	          parent_secondary = parent_position(i_level_py,i_secondary)%with_parent_height(parent_level_py)
		
	          HI_column_density_in = HI_column_density_in + &
	                            nx_py_nz_dz(-parent_level_py)%HI_density(-parent_primary,parent_secondary) * cellsize
			              
	          HeI_column_density_in = HeI_column_density_in + &
	                            nx_py_nz_dz(-parent_level_py)%HeI_density(-parent_primary,parent_secondary) * cellsize
						  							  
	          HeII_column_density_in = HeII_column_density_in + &
	                            nx_py_nz_dz(-parent_level_py)%HeII_density(-parent_primary,parent_secondary) * cellsize 
						  
	        enddo	
	  
	        nx_py_nz_dz(-i_level_py)%HI_column_density_in(-i_primary,i_secondary) = HI_column_density_in	
	        nx_py_nz_dz(-i_level_py)%HeI_column_density_in(-i_primary,i_secondary) = HeI_column_density_in
	        nx_py_nz_dz(-i_level_py)%HeII_column_density_in(-i_primary,i_secondary) = HeII_column_density_in		
	
	      else
		  
	        nx_py_nz_dz(-i_level_py)%HI_column_density_in(-i_primary,i_secondary) = &
					  nx_py_nz_dz(-i_level_py+1)%HI_column_density_out(-i_primary,i_secondary)
	        nx_py_nz_dz(-i_level_py)%HeI_column_density_in(-i_primary,i_secondary) = &
					  nx_py_nz_dz(-i_level_py+1)%HeI_column_density_out(-i_primary,i_secondary)		
	        nx_py_nz_dz(-i_level_py)%HeII_column_density_in(-i_primary,i_secondary) = &
					  nx_py_nz_dz(-i_level_py+1)%HeII_column_density_out(-i_primary,i_secondary)
		  
	      endif
	  			    			
		  nx_py_nz_dz(-i_level_py)%HI_column_density_out(-i_primary,i_secondary) = &
	                       nx_py_nz_dz(-i_level_py)%HI_column_density_in(-i_primary,i_secondary) + &
	                       nx_py_nz_dz(-i_level_py)%HI_density(-i_primary,i_secondary) * cellsize
		  nx_py_nz_dz(-i_level_py)%HeI_column_density_out(-i_primary,i_secondary) = &
	                       nx_py_nz_dz(-i_level_py)%HeI_column_density_in(-i_primary,i_secondary) + &
	                       nx_py_nz_dz(-i_level_py)%HeI_density(-i_primary,i_secondary) * cellsize	  
		  nx_py_nz_dz(-i_level_py)%HeII_column_density_out(-i_primary,i_secondary) = &
	                       nx_py_nz_dz(-i_level_py)%HeII_column_density_in(-i_primary,i_secondary) + &
	                       nx_py_nz_dz(-i_level_py)%HeII_density(-i_primary,i_secondary) * cellsize		  
	  	
		enddo	  
	  enddo	  
	enddo
	  
  end subroutine nx_py_nz_dz_column_density_calculation
  	 
  subroutine nx_ny_pz_dx_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	integer :: parent_primary,parent_secondary,parent_level_py
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_cell
	real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_cell
    real(kind=dp) :: cellsize

	! i_level_py = 1
	nx_ny_pz_dx(-1)%HI_column_density_in(-1,1) = 0
	nx_ny_pz_dx(-1)%HeI_column_density_in(-1,1) = 0
	nx_ny_pz_dx(-1)%HeII_column_density_in(-1,1) = 0

	cellsize = pyramid_grid(1)%cellsize(1,1)

	nx_ny_pz_dx(-1)%HI_column_density_out(-1,1) = nx_ny_pz_dx(-1)%HI_density(-1,1) * cellsize	
	nx_ny_pz_dx(-1)%HeI_column_density_out(-1,1) = nx_ny_pz_dx(-1)%HeI_density(-1,1) * cellsize
	nx_ny_pz_dx(-1)%HeII_column_density_out(-1,1) = nx_ny_pz_dx(-1)%HeII_density(-1,1) * cellsize
		
	! i_level_py > 1
	do i_level_py = 2,level_py
	  do i_primary = 1,partition(i_level_py)
	    do i_secondary = 1,partition(i_level_py)
		
		  cellsize = pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
	  
		  if (layer(i_level_py).ne.layer(i_level_py-1)) then
		  
	        HI_column_density_in = 0
	        HeI_column_density_in = 0
	        HeII_column_density_in = 0
		  		  
	        do parent_level_py = 1,i_level_py-1
		
	          parent_primary = parent_position(i_level_py,i_primary)%with_parent_height(parent_level_py)
	          parent_secondary = parent_position(i_level_py,i_secondary)%with_parent_height(parent_level_py)
		
	          HI_column_density_in = HI_column_density_in + &
	                            nx_ny_pz_dx(-parent_level_py)%HI_density(-parent_primary,parent_secondary) * cellsize
			              
	          HeI_column_density_in = HeI_column_density_in + &
	                            nx_ny_pz_dx(-parent_level_py)%HeI_density(-parent_primary,parent_secondary) * cellsize
						  							  
	          HeII_column_density_in = HeII_column_density_in + &
	                            nx_ny_pz_dx(-parent_level_py)%HeII_density(-parent_primary,parent_secondary) * cellsize 
						  
	        enddo	
	  
	        nx_ny_pz_dx(-i_level_py)%HI_column_density_in(-i_primary,i_secondary) = HI_column_density_in	
	        nx_ny_pz_dx(-i_level_py)%HeI_column_density_in(-i_primary,i_secondary) = HeI_column_density_in
	        nx_ny_pz_dx(-i_level_py)%HeII_column_density_in(-i_primary,i_secondary) = HeII_column_density_in		
	
	      else
		  
	        nx_ny_pz_dx(-i_level_py)%HI_column_density_in(-i_primary,i_secondary) = &
					  nx_ny_pz_dx(-i_level_py+1)%HI_column_density_out(-i_primary,i_secondary)
	        nx_ny_pz_dx(-i_level_py)%HeI_column_density_in(-i_primary,i_secondary) = &
					  nx_ny_pz_dx(-i_level_py+1)%HeI_column_density_out(-i_primary,i_secondary)		
	        nx_ny_pz_dx(-i_level_py)%HeII_column_density_in(-i_primary,i_secondary) = &
					  nx_ny_pz_dx(-i_level_py+1)%HeII_column_density_out(-i_primary,i_secondary)
		  
	      endif
	  			    			
		  nx_ny_pz_dx(-i_level_py)%HI_column_density_out(-i_primary,i_secondary) = &
	                       nx_ny_pz_dx(-i_level_py)%HI_column_density_in(-i_primary,i_secondary) + &
	                       nx_ny_pz_dx(-i_level_py)%HI_density(-i_primary,i_secondary) * cellsize
		  nx_ny_pz_dx(-i_level_py)%HeI_column_density_out(-i_primary,i_secondary) = &
	                       nx_ny_pz_dx(-i_level_py)%HeI_column_density_in(-i_primary,i_secondary) + &
	                       nx_ny_pz_dx(-i_level_py)%HeI_density(-i_primary,i_secondary) * cellsize	  
		  nx_ny_pz_dx(-i_level_py)%HeII_column_density_out(-i_primary,i_secondary) = &
	                       nx_ny_pz_dx(-i_level_py)%HeII_column_density_in(-i_primary,i_secondary) + &
	                       nx_ny_pz_dx(-i_level_py)%HeII_density(-i_primary,i_secondary) * cellsize		  
	  	
		enddo	  
	  enddo	  
	enddo  
	  
  end subroutine nx_ny_pz_dx_column_density_calculation
  
  subroutine nx_ny_pz_dy_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	integer :: parent_primary,parent_secondary,parent_level_py
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_cell
	real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_cell
    real(kind=dp) :: cellsize

	! i_level_py = 1
	nx_ny_pz_dy(-1)%HI_column_density_in(1,-1) = 0
	nx_ny_pz_dy(-1)%HeI_column_density_in(1,-1) = 0
	nx_ny_pz_dy(-1)%HeII_column_density_in(1,-1) = 0

	cellsize = pyramid_grid(1)%cellsize(1,1)

	nx_ny_pz_dy(-1)%HI_column_density_out(1,-1) = nx_ny_pz_dy(-1)%HI_density(1,-1) * cellsize	
	nx_ny_pz_dy(-1)%HeI_column_density_out(1,-1) = nx_ny_pz_dy(-1)%HeI_density(1,-1) * cellsize
	nx_ny_pz_dy(-1)%HeII_column_density_out(1,-1) = nx_ny_pz_dy(-1)%HeII_density(1,-1) * cellsize
		
	! i_level_py > 1
	do i_level_py = 2,level_py
	  do i_primary = 1,partition(i_level_py)
	    do i_secondary = 1,partition(i_level_py)
		
		  cellsize = pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
	  
		  if (layer(i_level_py).ne.layer(i_level_py-1)) then
		  
	        HI_column_density_in = 0
	        HeI_column_density_in = 0
	        HeII_column_density_in = 0
		  		  
	        do parent_level_py = 1,i_level_py-1
		
	          parent_primary = parent_position(i_level_py,i_primary)%with_parent_height(parent_level_py)
	          parent_secondary = parent_position(i_level_py,i_secondary)%with_parent_height(parent_level_py)
		
	          HI_column_density_in = HI_column_density_in + &
	                            nx_ny_pz_dy(-parent_level_py)%HI_density(parent_primary,-parent_secondary) * cellsize
			              
	          HeI_column_density_in = HeI_column_density_in + &
	                            nx_ny_pz_dy(-parent_level_py)%HeI_density(parent_primary,-parent_secondary) * cellsize
						  							  
	          HeII_column_density_in = HeII_column_density_in + &
	                            nx_ny_pz_dy(-parent_level_py)%HeII_density(parent_primary,-parent_secondary) * cellsize 
						  
	        enddo	
	  
	        nx_ny_pz_dy(-i_level_py)%HI_column_density_in(i_primary,-i_secondary) = HI_column_density_in	
	        nx_ny_pz_dy(-i_level_py)%HeI_column_density_in(i_primary,-i_secondary) = HeI_column_density_in
	        nx_ny_pz_dy(-i_level_py)%HeII_column_density_in(i_primary,-i_secondary) = HeII_column_density_in		
	
	      else
		  
	        nx_ny_pz_dy(-i_level_py)%HI_column_density_in(i_primary,-i_secondary) = &
					  nx_ny_pz_dy(-i_level_py+1)%HI_column_density_out(i_primary,-i_secondary)
	        nx_ny_pz_dy(-i_level_py)%HeI_column_density_in(i_primary,-i_secondary) = &
					  nx_ny_pz_dy(-i_level_py+1)%HeI_column_density_out(i_primary,-i_secondary)		
	        nx_ny_pz_dy(-i_level_py)%HeII_column_density_in(i_primary,-i_secondary) = &
					  nx_ny_pz_dy(-i_level_py+1)%HeII_column_density_out(i_primary,-i_secondary)
		  
	      endif
	  			    			
		  nx_ny_pz_dy(-i_level_py)%HI_column_density_out(i_primary,-i_secondary) = &
	                       nx_ny_pz_dy(-i_level_py)%HI_column_density_in(i_primary,-i_secondary) + &
	                       nx_ny_pz_dy(-i_level_py)%HI_density(i_primary,-i_secondary) * cellsize
		  nx_ny_pz_dy(-i_level_py)%HeI_column_density_out(i_primary,-i_secondary) = &
	                       nx_ny_pz_dy(-i_level_py)%HeI_column_density_in(i_primary,-i_secondary) + &
	                       nx_ny_pz_dy(-i_level_py)%HeI_density(i_primary,-i_secondary) * cellsize	  
		  nx_ny_pz_dy(-i_level_py)%HeII_column_density_out(i_primary,-i_secondary) = &
	                       nx_ny_pz_dy(-i_level_py)%HeII_column_density_in(i_primary,-i_secondary) + &
	                       nx_ny_pz_dy(-i_level_py)%HeII_density(i_primary,-i_secondary) * cellsize		  
	  	
		enddo	  
	  enddo	  
	enddo
	  
  end subroutine nx_ny_pz_dy_column_density_calculation
  
  subroutine nx_ny_pz_dz_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	integer :: parent_primary,parent_secondary,parent_level_py
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_cell
	real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_cell
    real(kind=dp) :: cellsize
	
	! i_level_py = 1
	nx_ny_pz_dz(1)%HI_column_density_in(-1,-1) = 0
	nx_ny_pz_dz(1)%HeI_column_density_in(-1,-1) = 0
	nx_ny_pz_dz(1)%HeII_column_density_in(-1,-1) = 0

	cellsize = pyramid_grid(1)%cellsize(1,1)

	nx_ny_pz_dz(1)%HI_column_density_out(-1,-1) = nx_ny_pz_dz(1)%HI_density(-1,-1) * cellsize	
	nx_ny_pz_dz(1)%HeI_column_density_out(-1,-1) = nx_ny_pz_dz(1)%HeI_density(-1,-1) * cellsize
	nx_ny_pz_dz(1)%HeII_column_density_out(-1,-1) = nx_ny_pz_dz(1)%HeII_density(-1,-1) * cellsize

	! i_level_py > 1
	do i_level_py = 2,level_py
	  do i_primary = 1,partition(i_level_py)
	    do i_secondary = 1,partition(i_level_py)
		
		  cellsize = pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
	  
		  if (layer(i_level_py).ne.layer(i_level_py-1)) then
		  
	        HI_column_density_in = 0
	        HeI_column_density_in = 0
	        HeII_column_density_in = 0
		  		  
	        do parent_level_py = 1,i_level_py-1
		
	          parent_primary = parent_position(i_level_py,i_primary)%with_parent_height(parent_level_py)
	          parent_secondary = parent_position(i_level_py,i_secondary)%with_parent_height(parent_level_py)
		
	          HI_column_density_in = HI_column_density_in + &
	                            nx_ny_pz_dz(parent_level_py)%HI_density(-parent_primary,-parent_secondary) * cellsize
			              
	          HeI_column_density_in = HeI_column_density_in + &
	                            nx_ny_pz_dz(parent_level_py)%HeI_density(-parent_primary,-parent_secondary) * cellsize
						  							  
	          HeII_column_density_in = HeII_column_density_in + &
	                            nx_ny_pz_dz(parent_level_py)%HeII_density(-parent_primary,-parent_secondary) * cellsize 
						  
	        enddo	
	  
	        nx_ny_pz_dz(i_level_py)%HI_column_density_in(-i_primary,-i_secondary) = HI_column_density_in	
	        nx_ny_pz_dz(i_level_py)%HeI_column_density_in(-i_primary,-i_secondary) = HeI_column_density_in
	        nx_ny_pz_dz(i_level_py)%HeII_column_density_in(-i_primary,-i_secondary) = HeII_column_density_in		
	
	      else
		  
	        nx_ny_pz_dz(i_level_py)%HI_column_density_in(-i_primary,-i_secondary) = &
					  nx_ny_pz_dz(i_level_py-1)%HI_column_density_out(-i_primary,-i_secondary)
	        nx_ny_pz_dz(i_level_py)%HeI_column_density_in(-i_primary,-i_secondary) = &
					  nx_ny_pz_dz(i_level_py-1)%HeI_column_density_out(-i_primary,-i_secondary)		
	        nx_ny_pz_dz(i_level_py)%HeII_column_density_in(-i_primary,-i_secondary) = &
					  nx_ny_pz_dz(i_level_py-1)%HeII_column_density_out(-i_primary,-i_secondary)
		  
	      endif
	  			    			
		  nx_ny_pz_dz(i_level_py)%HI_column_density_out(-i_primary,-i_secondary) = &
	                       nx_ny_pz_dz(i_level_py)%HI_column_density_in(-i_primary,-i_secondary) + &
	                       nx_ny_pz_dz(i_level_py)%HI_density(-i_primary,-i_secondary) * cellsize
		  nx_ny_pz_dz(i_level_py)%HeI_column_density_out(-i_primary,-i_secondary) = &
	                       nx_ny_pz_dz(i_level_py)%HeI_column_density_in(-i_primary,-i_secondary) + &
	                       nx_ny_pz_dz(i_level_py)%HeI_density(-i_primary,-i_secondary) * cellsize	  
		  nx_ny_pz_dz(i_level_py)%HeII_column_density_out(-i_primary,-i_secondary) = &
	                       nx_ny_pz_dz(i_level_py)%HeII_column_density_in(-i_primary,-i_secondary) + &
	                       nx_ny_pz_dz(i_level_py)%HeII_density(-i_primary,-i_secondary) * cellsize		  
	  	
		enddo	  
	  enddo	  
	enddo
	  
  end subroutine nx_ny_pz_dz_column_density_calculation
     
  subroutine nx_ny_nz_dx_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	integer :: parent_primary,parent_secondary,parent_level_py
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_cell
	real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_cell
    real(kind=dp) :: cellsize

	! i_level_py = 1
	nx_ny_nz_dx(-1)%HI_column_density_in(-1,-1) = 0
	nx_ny_nz_dx(-1)%HeI_column_density_in(-1,-1) = 0
	nx_ny_nz_dx(-1)%HeII_column_density_in(-1,-1) = 0

	cellsize = pyramid_grid(1)%cellsize(1,1)

	nx_ny_nz_dx(-1)%HI_column_density_out(-1,-1) = nx_ny_nz_dx(-1)%HI_density(-1,-1) * cellsize	
	nx_ny_nz_dx(-1)%HeI_column_density_out(-1,-1) = nx_ny_nz_dx(-1)%HeI_density(-1,-1) * cellsize
	nx_ny_nz_dx(-1)%HeII_column_density_out(-1,-1) = nx_ny_nz_dx(-1)%HeII_density(-1,-1) * cellsize	
	
	! i_level_py > 1
	do i_level_py = 2,level_py
	  do i_primary = 1,partition(i_level_py)
	    do i_secondary = 1,partition(i_level_py)
		
		  cellsize = pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
	  
		  if (layer(i_level_py).ne.layer(i_level_py-1)) then
		  
	        HI_column_density_in = 0
	        HeI_column_density_in = 0
	        HeII_column_density_in = 0
		  		  
	        do parent_level_py = 1,i_level_py-1
		
	          parent_primary = parent_position(i_level_py,i_primary)%with_parent_height(parent_level_py)
	          parent_secondary = parent_position(i_level_py,i_secondary)%with_parent_height(parent_level_py)
		
	          HI_column_density_in = HI_column_density_in + &
	                            nx_ny_nz_dx(-parent_level_py)%HI_density(-parent_primary,-parent_secondary) * cellsize
			              
	          HeI_column_density_in = HeI_column_density_in + &
	                            nx_ny_nz_dx(-parent_level_py)%HeI_density(-parent_primary,-parent_secondary) * cellsize
						  							  
	          HeII_column_density_in = HeII_column_density_in + &
	                            nx_ny_nz_dx(-parent_level_py)%HeII_density(-parent_primary,-parent_secondary) * cellsize 
						  
	        enddo	
	  
	        nx_ny_nz_dx(-i_level_py)%HI_column_density_in(-i_primary,-i_secondary) = HI_column_density_in	
	        nx_ny_nz_dx(-i_level_py)%HeI_column_density_in(-i_primary,-i_secondary) = HeI_column_density_in
	        nx_ny_nz_dx(-i_level_py)%HeII_column_density_in(-i_primary,-i_secondary) = HeII_column_density_in		
	
	      else
		  
	        nx_ny_nz_dx(-i_level_py)%HI_column_density_in(-i_primary,-i_secondary) = &
					  nx_ny_nz_dx(-i_level_py+1)%HI_column_density_out(-i_primary,-i_secondary)
	        nx_ny_nz_dx(-i_level_py)%HeI_column_density_in(-i_primary,-i_secondary) = &
					  nx_ny_nz_dx(-i_level_py+1)%HeI_column_density_out(-i_primary,-i_secondary)		
	        nx_ny_nz_dx(-i_level_py)%HeII_column_density_in(-i_primary,-i_secondary) = &
					  nx_ny_nz_dx(-i_level_py+1)%HeII_column_density_out(-i_primary,-i_secondary)
		  
	      endif
	  			    			
		  nx_ny_nz_dx(-i_level_py)%HI_column_density_out(-i_primary,-i_secondary) = &
	                       nx_ny_nz_dx(-i_level_py)%HI_column_density_in(-i_primary,-i_secondary) + &
	                       nx_ny_nz_dx(-i_level_py)%HI_density(-i_primary,-i_secondary) * cellsize
		  nx_ny_nz_dx(-i_level_py)%HeI_column_density_out(-i_primary,-i_secondary) = &
	                       nx_ny_nz_dx(-i_level_py)%HeI_column_density_in(-i_primary,-i_secondary) + &
	                       nx_ny_nz_dx(-i_level_py)%HeI_density(-i_primary,-i_secondary) * cellsize	  
		  nx_ny_nz_dx(-i_level_py)%HeII_column_density_out(-i_primary,-i_secondary) = &
	                       nx_ny_nz_dx(-i_level_py)%HeII_column_density_in(-i_primary,-i_secondary) + &
	                       nx_ny_nz_dx(-i_level_py)%HeII_density(-i_primary,-i_secondary) * cellsize		  
	  	
		enddo	  
	  enddo	  
	enddo
	  
  end subroutine nx_ny_nz_dx_column_density_calculation
  
  subroutine nx_ny_nz_dy_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	integer :: parent_primary,parent_secondary,parent_level_py
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_cell
	real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_cell
    real(kind=dp) :: cellsize
	
	! i_level_py = 1
	nx_ny_nz_dy(-1)%HI_column_density_in(-1,-1) = 0
	nx_ny_nz_dy(-1)%HeI_column_density_in(-1,-1) = 0
	nx_ny_nz_dy(-1)%HeII_column_density_in(-1,-1) = 0

	cellsize = pyramid_grid(1)%cellsize(1,1)

	nx_ny_nz_dy(-1)%HI_column_density_out(-1,-1) = nx_ny_nz_dy(-1)%HI_density(-1,-1) * cellsize	
	nx_ny_nz_dy(-1)%HeI_column_density_out(-1,-1) = nx_ny_nz_dy(-1)%HeI_density(-1,-1) * cellsize
	nx_ny_nz_dy(-1)%HeII_column_density_out(-1,-1) = nx_ny_nz_dy(-1)%HeII_density(-1,-1) * cellsize
	
	! i_level_py > 1
	do i_level_py = 2,level_py
	  do i_primary = 1,partition(i_level_py)
	    do i_secondary = 1,partition(i_level_py)
	
		  cellsize = pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
  
		  if (layer(i_level_py).ne.layer(i_level_py-1)) then
	  
	        HI_column_density_in = 0
	        HeI_column_density_in = 0
	        HeII_column_density_in = 0
	  		  
	        do parent_level_py = 1,i_level_py-1
	
	          parent_primary = parent_position(i_level_py,i_primary)%with_parent_height(parent_level_py)
	          parent_secondary = parent_position(i_level_py,i_secondary)%with_parent_height(parent_level_py)
	
	          HI_column_density_in = HI_column_density_in + &
	                            nx_ny_nz_dy(-parent_level_py)%HI_density(-parent_primary,-parent_secondary) * cellsize
		              
	          HeI_column_density_in = HeI_column_density_in + &
	                            nx_ny_nz_dy(-parent_level_py)%HeI_density(-parent_primary,-parent_secondary) * cellsize
					  							  
	          HeII_column_density_in = HeII_column_density_in + &
	                            nx_ny_nz_dy(-parent_level_py)%HeII_density(-parent_primary,-parent_secondary) * cellsize 
					  
	        enddo	
  
	        nx_ny_nz_dy(-i_level_py)%HI_column_density_in(-i_primary,-i_secondary) = HI_column_density_in	
	        nx_ny_nz_dy(-i_level_py)%HeI_column_density_in(-i_primary,-i_secondary) = HeI_column_density_in
	        nx_ny_nz_dy(-i_level_py)%HeII_column_density_in(-i_primary,-i_secondary) = HeII_column_density_in		

	      else
	  
	        nx_ny_nz_dy(-i_level_py)%HI_column_density_in(-i_primary,-i_secondary) = &
					  nx_ny_nz_dy(-i_level_py+1)%HI_column_density_out(-i_primary,-i_secondary)
	        nx_ny_nz_dy(-i_level_py)%HeI_column_density_in(-i_primary,-i_secondary) = &
					  nx_ny_nz_dy(-i_level_py+1)%HeI_column_density_out(-i_primary,-i_secondary)		
	        nx_ny_nz_dy(-i_level_py)%HeII_column_density_in(-i_primary,-i_secondary) = &
					  nx_ny_nz_dy(-i_level_py+1)%HeII_column_density_out(-i_primary,-i_secondary)
	  
	      endif
  			    			
		  nx_ny_nz_dy(-i_level_py)%HI_column_density_out(-i_primary,-i_secondary) = &
	                       nx_ny_nz_dy(-i_level_py)%HI_column_density_in(-i_primary,-i_secondary) + &
	                       nx_ny_nz_dy(-i_level_py)%HI_density(-i_primary,-i_secondary) * cellsize
		  nx_ny_nz_dy(-i_level_py)%HeI_column_density_out(-i_primary,-i_secondary) = &
	                       nx_ny_nz_dy(-i_level_py)%HeI_column_density_in(-i_primary,-i_secondary) + &
	                       nx_ny_nz_dy(-i_level_py)%HeI_density(-i_primary,-i_secondary) * cellsize	  
		  nx_ny_nz_dy(-i_level_py)%HeII_column_density_out(-i_primary,-i_secondary) = &
	                       nx_ny_nz_dy(-i_level_py)%HeII_column_density_in(-i_primary,-i_secondary) + &
	                       nx_ny_nz_dy(-i_level_py)%HeII_density(-i_primary,-i_secondary) * cellsize		  
  	
		enddo	  
	  enddo	  
	enddo 
	  
  end subroutine nx_ny_nz_dy_column_density_calculation
  
  subroutine nx_ny_nz_dz_column_density_calculation ()
	  
    implicit none 
	
    integer :: i_primary,i_secondary,i_level_py
    integer :: o_primary,o_secondary
	integer :: parent_primary,parent_secondary,parent_level_py
	real(kind=dp) :: HI_column_density_in
    real(kind=dp) :: HI_column_density_cell
	real(kind=dp) :: HeI_column_density_in
    real(kind=dp) :: HeI_column_density_cell
	real(kind=dp) :: HeII_column_density_in
    real(kind=dp) :: HeII_column_density_cell
    real(kind=dp) :: cellsize
	
	! i_level_py = 1
	nx_ny_nz_dz(-1)%HI_column_density_in(-1,-1) = 0
	nx_ny_nz_dz(-1)%HeI_column_density_in(-1,-1) = 0
	nx_ny_nz_dz(-1)%HeII_column_density_in(-1,-1) = 0

	cellsize = pyramid_grid(1)%cellsize(1,1)

	nx_ny_nz_dz(-1)%HI_column_density_out(-1,-1) = nx_ny_nz_dz(-1)%HI_density(-1,-1) * cellsize	
	nx_ny_nz_dz(-1)%HeI_column_density_out(-1,-1) = nx_ny_nz_dz(-1)%HeI_density(-1,-1) * cellsize
	nx_ny_nz_dz(-1)%HeII_column_density_out(-1,-1) = nx_ny_nz_dz(-1)%HeII_density(-1,-1) * cellsize
		
	! i_level_py > 1
	do i_level_py = 2,level_py
	  do i_primary = 1,partition(i_level_py)
	    do i_secondary = 1,partition(i_level_py)
	
		  cellsize = pyramid_grid(i_level_py)%cellsize(i_primary,i_secondary)
  
		  if (layer(i_level_py).ne.layer(i_level_py-1)) then
	  
	        HI_column_density_in = 0
	        HeI_column_density_in = 0
	        HeII_column_density_in = 0
	  		  
	        do parent_level_py = 1,i_level_py-1
	
	          parent_primary = parent_position(i_level_py,i_primary)%with_parent_height(parent_level_py)
	          parent_secondary = parent_position(i_level_py,i_secondary)%with_parent_height(parent_level_py)
	
	          HI_column_density_in = HI_column_density_in + &
	                            nx_ny_nz_dz(-parent_level_py)%HI_density(-parent_primary,-parent_secondary) * cellsize
		              
	          HeI_column_density_in = HeI_column_density_in + &
	                            nx_ny_nz_dz(-parent_level_py)%HeI_density(-parent_primary,-parent_secondary) * cellsize
					  							  
	          HeII_column_density_in = HeII_column_density_in + &
	                            nx_ny_nz_dz(-parent_level_py)%HeII_density(-parent_primary,-parent_secondary) * cellsize 
					  
	        enddo	
  
	        nx_ny_nz_dz(-i_level_py)%HI_column_density_in(-i_primary,-i_secondary) = HI_column_density_in	
	        nx_ny_nz_dz(-i_level_py)%HeI_column_density_in(-i_primary,-i_secondary) = HeI_column_density_in
	        nx_ny_nz_dz(-i_level_py)%HeII_column_density_in(-i_primary,-i_secondary) = HeII_column_density_in		

	      else
	  
	        nx_ny_nz_dz(-i_level_py)%HI_column_density_in(-i_primary,-i_secondary) = &
					  nx_ny_nz_dz(-i_level_py+1)%HI_column_density_out(-i_primary,-i_secondary)
	        nx_ny_nz_dz(-i_level_py)%HeI_column_density_in(-i_primary,-i_secondary) = &
					  nx_ny_nz_dz(-i_level_py+1)%HeI_column_density_out(-i_primary,-i_secondary)		
	        nx_ny_nz_dz(-i_level_py)%HeII_column_density_in(-i_primary,-i_secondary) = &
					  nx_ny_nz_dz(-i_level_py+1)%HeII_column_density_out(-i_primary,-i_secondary)
	  
	      endif
  			    			
		  nx_ny_nz_dz(-i_level_py)%HI_column_density_out(-i_primary,-i_secondary) = &
	                       nx_ny_nz_dz(-i_level_py)%HI_column_density_in(-i_primary,-i_secondary) + &
	                       nx_ny_nz_dz(-i_level_py)%HI_density(-i_primary,-i_secondary) * cellsize
		  nx_ny_nz_dz(-i_level_py)%HeI_column_density_out(-i_primary,-i_secondary) = &
	                       nx_ny_nz_dz(-i_level_py)%HeI_column_density_in(-i_primary,-i_secondary) + &
	                       nx_ny_nz_dz(-i_level_py)%HeI_density(-i_primary,-i_secondary) * cellsize	  
		  nx_ny_nz_dz(-i_level_py)%HeII_column_density_out(-i_primary,-i_secondary) = &
	                       nx_ny_nz_dz(-i_level_py)%HeII_column_density_in(-i_primary,-i_secondary) + &
	                       nx_ny_nz_dz(-i_level_py)%HeII_density(-i_primary,-i_secondary) * cellsize		  
  	
		enddo	  
	  enddo	  
	enddo	  
	  
  end subroutine nx_ny_nz_dz_column_density_calculation
  
end module pyramid_column_density
