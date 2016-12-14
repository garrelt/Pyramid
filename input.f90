module input

  use precision, only: dp

implicit none

  ! number of levels in one pyramid
  real(kind=dp) :: box_size
  real(kind=dp),parameter :: nominal_box_size = 2.0
  integer,parameter :: level_py = 50
  integer,parameter :: level_sl = 2*level_py+1
  real(kind=dp) :: cellsize_py
  real(kind=dp) :: cellsize_py_cube 
  real(kind=dp) :: cellsize_sl
  real(kind=dp) :: volume_factor
  real(kind=dp),parameter :: abu_he=0.0795
  real(kind=dp),parameter :: abu_h=1.0-abu_he
  real(kind=dp),parameter :: number_density = 7.0e-5 ! Hydrogen density
  real(kind=dp),parameter :: xHI = 1.0e-8_dp ! HI fraction
  real(kind=dp),parameter :: xHeI = 1.0e-18_dp ! HI fraction
  real(kind=dp),parameter :: xHeII = 1.0e-18_dp ! HI fraction 
  real(kind=dp),parameter :: clumping = 1.0 ! clumping factor
  real(kind=dp) :: bb_source_temperature
  real(kind=dp),parameter :: bb_nominal_source_temperature = 50000.0
  real(kind=dp),parameter :: S_star = 1.0e55_dp
  real(kind=dp),parameter :: R_solar = 6.9599e10 !< Solar radius
  real(kind=dp),parameter :: L_solar = 3.826e33 !< Solar luminosity  
  real(kind=dp),parameter :: pl_input_flux = 1.0e55_dp
  real(kind=dp) :: pl_source_index
  real(kind=dp),parameter :: pl_nominal_source_index = 2.5
  real(kind=dp),parameter :: pl_minfreq = 13.6
  real(kind=dp),parameter :: pl_maxfreq = 5440.0
  real(kind=dp),parameter :: pi = 3.14159265359
  real(kind=dp),parameter :: four_pi = 4.0*pi 
  real(kind=dp),parameter :: four_over_three_pi = 4.0*pi/3.0  
  real(kind=dp),parameter :: epsilon = 1e-10
  real(kind=dp) :: luminosity
  real(kind=dp),parameter :: nominal_luminosity = 5.0e50_dp
  character :: use_which_source
  character,parameter :: nominal_use_which_source = 'B'
  character :: use_which_field
  character,parameter :: nominal_use_which_field = 'C'

  !> HI cross section at its ionzing frequency
  real(kind=dp), parameter :: sigma_HI_at_ion_freq=6.346e-18
  !> HeI cross section at its ionzing frequency
  real(kind=dp), parameter :: sigma_HeI_at_ion_freq=7.430e-18
  !> HeII cross section at its ionzing frequency
  real(kind=dp), parameter :: sigma_HeII_at_ion_freq=1.589e-18
  !> HI ionization energy in frequency
  real(kind=dp), parameter :: ion_freq_HI=0.241838e15*13.598
  !> HeI ionization energy in frequency
  real(kind=dp), parameter :: ion_freq_HeI=0.241838e15*24.587
  !> HeII ionization energy in frequency
  real(kind=dp), parameter :: ion_freq_HeII=0.241838e15*54.416
  
  integer,parameter :: number_of_source = 1
  integer,dimension(1:3,1:number_of_source) :: source_position 
  real(kind=dp),dimension(1:number_of_source) :: source_luminosity
  integer,dimension(1:level_py) :: partition
  integer,dimension(1:level_py) :: layer
  
  logical, parameter :: ionization_weighted_by_atom = .true.
  !logical, parameter :: ionization_weighted_by_atom = .false.

  logical, parameter :: heating_weighted_by_atom = .true.
  !logical, parameter :: heating_weighted_by_atom = .false.

contains

  subroutine setup_input()
	
    implicit none

    character(len=512) :: inputfile
    character(len=10) :: multiple_string !< temporary
    real(kind=dp) :: multiple !< temporary
    character(len=1) :: answer
    logical :: file_input
    integer, parameter :: stdinput = 10

    if (COMMAND_ARGUMENT_COUNT () .gt. 0) then
      call GET_COMMAND_ARGUMENT(1,inputfile)
      open(unit=stdinput,file=inputfile)
      file_input = .true.
    else
      file_input = .false.
    endif

    ! Ask for black-body source temperature
    if (.not.file_input) then
      write(*,"(A,$)") "Enter black-body source temperature (K): "
    endif
    read(stdinput,*) bb_source_temperature
    if (bb_source_temperature .lt. 2000 .or. bb_source_temperature .gt. 200000) then
      write(*,*) "Black-body source temperature is out of range (2000K < T < 200000K)"
      write(*,*) "Black-body source temperature will be set to the nominal value"
      bb_source_temperature = bb_nominal_source_temperature
    endif

      ! Ask for power law source index
      if (.not.file_input) then
        write(*,"(A,$)") "Enter power-law source index (number of photons): "
      endif
      read(stdinput,*) pl_source_index
      if (pl_source_index .lt. 0) then
        write(*,*) "Power-law source index is out of range (index >= 0)"
        write(*,*) "Power-law source index will be set to the nominal value"
        pl_source_index = pl_nominal_source_index
      endif

      ! Ask for box size in Mpc
      if (.not.file_input) then
        write(*,"(A,$)") "Enter box size (Mpc): "
      endif
      read(stdinput,*) box_size
      if (box_size .lt. 0) then
        write(*,*) "box_size is out of range (d >= 0)"
        write(*,*) "box_size will be set to the nominal value"
        box_size = nominal_box_size
      endif
      box_size = box_size * 3.086e24

      ! Ask for source luminosity
      if (.not.file_input) then
        write(*,"(A,$)") "Enter source luminosity: "
      endif
      read(stdinput,*) luminosity
      if (box_size .lt. 0) then
        write(*,*) "luminosity is out of range (s >= 0)"
        write(*,*) "luminosity will be set to the nominal value"
        luminosity = nominal_luminosity
      endif

      ! use which source? 
      if (.not.file_input) then
        write(*,"(A,$)") "Use which source? "
      endif
      read(stdinput,*) answer
      if (answer.eq.'B' .or. answer.eq.'b') then
        use_which_source = 'B'
      elseif (answer.eq.'P' .or. answer.eq.'p') then
        use_which_source = 'P'
      else
        write(*,*) "use_which_source accepts only B, b, P or p"
        write(*,*) "use_which_source will be set to the nominal value"
        use_which_source = nominal_use_which_source
      endif

      ! use which field? 
      if (.not.file_input) then
        write(*,"(A,$)") "Use which field? "
      endif
      read(stdinput,*) answer
      if (answer.eq.'C' .or. answer.eq.'c') then
        use_which_field = 'C'
      elseif (answer.eq.'O' .or. answer.eq.'o') then
        use_which_field = 'O'
      else
        write(*,*) "use_which_field accepts only B, b, P or p"
        write(*,*) "use_which_field will be set to the nominal value"
        use_which_field = nominal_use_which_field
      endif

  cellsize_py = box_size/real(2.0*level_py)
  cellsize_py_cube = cellsize_py*cellsize_py*cellsize_py  
  cellsize_sl = box_size/level_sl
  volume_factor = (0.5*box_size/level_py)**3 

write(*,*) 'level_py',level_py
write(*,*) 'level_sl',level_sl
  end subroutine setup_input

  subroutine source_information()
	
    implicit none

    source_position(1,1) = level_py+1
    source_position(2,1) = level_py+1
    source_position(3,1) = level_py+1
    source_luminosity(1) = luminosity
	
  end subroutine source_information

  subroutine partition_setup()

    implicit none
	  
    integer :: count

    ! find out the partition of each level
    do count = 1,level_py
      if (int(log(real(count))/log(real(2))) .ne. log(real(count))/log(real(2))) then
        partition(count) = 2**(int(log(real(count))/log(real(2)))+1)
      else
        partition(count) = 2**(int(log(real(count))/log(real(2))))
      endif

    end do
	
  end subroutine partition_setup

  subroutine layer_setup()

    implicit none
	
    integer :: count

    ! find out the partition of each level
    do count = 1,level_py
      if (int(log(real(count))/log(real(2))) .ne. log(real(count))/log(real(2))) then
        layer(count) = int(log(real(count))/log(real(2)))+2
      else
        layer(count) = int(log(real(count))/log(real(2)))+1
      endif
    end do
	
  end subroutine layer_setup
  
end module input
