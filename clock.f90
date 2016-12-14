module clocks

  implicit none

  real :: cputime1_programme
  real :: cputime2_programme 
  real :: cputime1_pyramid
  real :: cputime2_pyramid  
  real :: cputime1_short
  real :: cputime2_short
  real :: cputime1_long
  real :: cputime2_long        
  integer :: cpu_hours = 0 
  integer :: cpu_minutes = 0 
  real :: cpu_seconds = 0.0 
 
  integer(kind=8) :: cntr1_programme 
  integer(kind=8) :: cntr2_programme 
  integer(kind=8) :: cntr1_pyramid
  integer(kind=8) :: cntr2_pyramid  
  integer(kind=8) :: cntr1_short
  integer(kind=8) :: cntr2_short
  integer(kind=8) :: cntr1_long
  integer(kind=8) :: cntr2_long      
  integer(kind=8) :: countspersec 
  integer :: clock_hours = 0
  integer :: clock_minutes = 0 
  real :: clock_seconds = 0.0
  
contains

  subroutine setup_clocks
    
    call cpu_time(cputime1_programme)
    call system_clock(cntr1_programme)
    
  end subroutine setup_clocks
  
  subroutine report_clocks
    
    call report_cpuclock
    call report_wallclock
    
  end subroutine report_clocks
 
  subroutine report_cpuclock
    
    cpu_hours = 0
    cpu_minutes = 0 
    cpu_seconds = 0.0 

    call cpu_time(cputime2_programme)
    cpu_seconds = cpu_seconds+real(cputime2_programme-cputime1_programme)
    cpu_minutes = cpu_minutes + int(cpu_seconds) / 60
    cpu_seconds = MOD ( cpu_seconds , 60.0 )
    cpu_hours = cpu_hours + cpu_minutes / 60
    cpu_minutes = MOD ( cpu_minutes , 60 )

    
!write(*,*) "CPU time: ", cpu_hours, ' hours', cpu_minutes, ' minutes', &
 !              cpu_seconds, ' seconds.'

  end subroutine report_cpuclock
  

  subroutine report_wallclock
    
    clock_hours = 0
    clock_minutes = 0
    clock_seconds = 0.0  
	   
    call system_clock(cntr2_programme,countspersec)
    clock_seconds = clock_seconds+real(cntr2_programme-cntr1_programme)/real(countspersec)
    clock_minutes = clock_minutes + int(clock_seconds) / 60
    clock_seconds = MOD ( clock_seconds , 60.0 )
    clock_hours = clock_hours + clock_minutes / 60
    clock_minutes = MOD ( clock_minutes , 60 )

    !write(*,*) "Wall clock time: ", clock_hours, ' hours', clock_minutes, ' minutes', &
     !          clock_seconds, ' seconds.'

  end subroutine report_wallclock
   
  subroutine start_pyramid_clock()
	  
    call cpu_time(cputime1_pyramid)	
    call system_clock(cntr1_pyramid) 
		  
  end subroutine start_pyramid_clock
  
  subroutine report_pyramid_clock()
	 
    call report_cpuclock_pyramid()
    call report_wallclock_pyramid()
	  
  end subroutine report_pyramid_clock

  subroutine report_cpuclock_pyramid() 
	 
    cpu_hours = 0 
    cpu_minutes = 0 
    cpu_seconds = 0.0 
	
    call cpu_time(cputime2_pyramid)
    cpu_seconds = cpu_seconds+real(cputime2_pyramid-cputime1_pyramid)
    cpu_minutes = cpu_minutes + int(cpu_seconds) / 60
    cpu_seconds = MOD ( cpu_seconds , 60.0 )
    cpu_hours = cpu_hours + cpu_minutes / 60
    cpu_minutes = MOD ( cpu_minutes , 60 )

    write(*,*) "Pyramid CPU time: ",cpu_seconds, ' seconds.'
	  
  end subroutine report_cpuclock_pyramid
 
  subroutine report_wallclock_pyramid()
	  
    clock_hours = 0 
    clock_minutes = 0
    clock_seconds = 0.0 	
	  
    call system_clock(cntr2_pyramid,countspersec)
    clock_seconds = clock_seconds+real(cntr2_pyramid-cntr1_pyramid)/real(countspersec)
    clock_minutes = clock_minutes + int(clock_seconds) / 60
    clock_seconds = MOD ( clock_seconds , 60.0 )
    clock_hours = clock_hours + clock_minutes / 60
    clock_minutes = MOD ( clock_minutes , 60 )

    !write(*,*) "Pyramid Wall clock time: ", clock_hours,' hours', clock_minutes, ' minutes', &
     !          clock_seconds, ' seconds.'
			   
  end subroutine report_wallclock_pyramid  
   
  subroutine start_short_clock()
	  
    call cpu_time(cputime1_short)
    call system_clock(cntr1_short) 
	  	  
  end subroutine start_short_clock
   
  subroutine report_short_clock()
	  
    call report_cpuclock_short()
    call report_wallclock_short()  
	
  end subroutine report_short_clock

  subroutine report_cpuclock_short() 
	  
    cpu_hours = 0 
    cpu_minutes = 0 
    cpu_seconds = 0.0 
   
    call cpu_time(cputime2_short)
    cpu_seconds = cpu_seconds+real(cputime2_short-cputime1_short)
    cpu_minutes = cpu_minutes + int(cpu_seconds) / 60
    cpu_seconds = MOD ( cpu_seconds , 60.0 )
    cpu_hours = cpu_hours + cpu_minutes / 60
    cpu_minutes = MOD ( cpu_minutes , 60 )

    write(*,*) "short CPU time: ", cpu_seconds, ' seconds.'
				  
  end subroutine report_cpuclock_short
 
  subroutine report_wallclock_short() 
	  
    clock_hours = 0
    clock_minutes = 0
    clock_seconds = 0.0 
	
    call system_clock(cntr2_short,countspersec)
    clock_seconds = clock_seconds+real(cntr2_short-cntr1_short)/real(countspersec)
    clock_minutes = clock_minutes + int(clock_seconds) / 60
    clock_seconds = MOD ( clock_seconds , 60.0 )
    clock_hours = clock_hours + clock_minutes / 60
    clock_minutes = MOD ( clock_minutes , 60 )

    !write(*,*) "short Wall clock time: ", clock_hours, ' hours', clock_minutes, ' minutes', &
     !          clock_seconds, ' seconds.'
			   
  end subroutine report_wallclock_short  

  subroutine start_long_clock()
	  
    call cpu_time(cputime1_long)
    call system_clock(cntr1_long)
	 	  
  end subroutine start_long_clock
  
  subroutine report_long_clock()
	  
    call report_cpuclock_long()
    call report_wallclock_long()  
	 
  end subroutine report_long_clock

  subroutine report_cpuclock_long()
	   
    cpu_hours = 0 
    cpu_minutes = 0 
    cpu_seconds = 0.0 
	
    call cpu_time(cputime2_long)
    cpu_seconds = cpu_seconds+real(cputime2_long-cputime1_long)
    cpu_minutes = cpu_minutes + int(cpu_seconds) / 60
    cpu_seconds = MOD ( cpu_seconds , 60.0 )
    cpu_hours = cpu_hours + cpu_minutes / 60
    cpu_minutes = MOD ( cpu_minutes , 60 )

     write(*,*) "long CPU time: ", cpu_seconds,' seconds.' 
				 
  end subroutine report_cpuclock_long
 
  subroutine report_wallclock_long()
	   
    clock_hours = 0 
    clock_minutes = 0 
    clock_seconds = 0.0 
		  
    call system_clock(cntr2_long,countspersec)
    clock_seconds = clock_seconds+real(cntr2_long-cntr1_long)/real(countspersec)
    clock_minutes = clock_minutes + int(clock_seconds) / 60
    clock_seconds = MOD ( clock_seconds , 60.0 )
    clock_hours = clock_hours + clock_minutes / 60
    clock_minutes = MOD ( clock_minutes , 60 )

    ! write(*,*) "long Wall clock time: ", clock_hours, ' hours', &
     !           clock_minutes, ' minutes', clock_seconds, ' seconds.'
				
  end subroutine report_wallclock_long  
  
end module clocks
