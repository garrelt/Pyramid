!     This module data and routines which deal with radiative
!     effects. 
!     Its main part deal with photo-ionizing radiation, but it
!     also initializes other radiative properties, such as cooling (which
!     are contained in different modules).
!     It can be used in hydrodynamic or stand-alone radiative transfer 
!     calculations.

module radiation

  use precision, only: dp
  use my_mpi
  use file_admin, only: logf
  use mathconstants, only: pi
  use cgsconstants, only: sigma_SB, &                     ! Stefan-Boltzmann constant
                          hplanck, &                      ! Planck constant
                          k_B, &                          ! Boltzmann constant
                          two_pi_over_c_square            ! two times pi over c aquare
  use cgsphotoconstants, only: ion_freq_HI, &                  ! HI ionization energy in frequency
                               T_low, &                        ! Minimum temperature of blackbody
                               T_high, &                       ! Maximum temperature of blackbody
                               pl_index_cross_section_HI, &    ! Power law index of cross section of HI
                               sigma_HI_at_ion_freq            !HI cross section at its ionzing frequency
  use astroconstants, only: R_SOLAR, &                    ! Solar radius
                            L_SOLAR                       ! Solar luminosity
  use romberg, only: scalar_romberg, &                    ! 1D integration function
                     vector_romberg, &                    ! 1D integration subroutine
                     romberg_initialisation               ! Romberg initialisation procedure
  use c2ray_parameters, only: T_eff_nominal, &            ! Black body  effective temperature for for nominal SED 
                              S_star_nominal              ! Ionizing photon rate for for nominal SED
  use material, only: isothermal

  implicit none

  integer,parameter :: NumFreq = 128       ! Number of integration points in each of the frequency bins
  integer,parameter :: NumTau = 2000       ! Number of table points for the optical depth
  
  ! Optical depths at the entrance of the grid.
  ! It can be used if radiation enters the simulation volume from the outside.
  real(kind=dp) :: boundary_tau = 0.0

  ! Parameters defining the optical depth entries in the table.
  real(kind=dp),parameter :: minlogtau = -20.0                             ! Table position starts at log10(minlogtau) 
  real(kind=dp),parameter :: maxlogtau = 4.0                               ! Table position ends at log10(maxlogtau)  
  real(kind=dp),parameter :: dlogtau = (maxlogtau-minlogtau)/real(NumTau)  ! dlogtau is the step size in log10(tau)

  ! Logical that determines the use of grey opacities
  logical,parameter :: grey = .false.


  ! max frequency for integration 
  real(kind=dp) :: freq_max1 
  real(kind=dp) :: freq_max2

  ! Stellar properties
  real(kind=dp) :: T_eff        ! Black body effective temperature
  real(kind=dp) :: R_star       ! Black body radius
  real(kind=dp) :: L_star       ! Black body luminosity
  real(kind=dp) :: S_star       ! Black body ionizing photons rate

  ! Frequency step width
  real(kind=dp) :: freq_step  

  real(kind=dp), dimension(:,:), allocatable  :: bb_photo_thick_integrand 
  real(kind=dp), dimension(:,:), allocatable  :: bb_heat_thick_integrand
  real(kind=dp), dimension(:,:), allocatable  :: bb_photo_thin_integrand
  real(kind=dp), dimension(:,:), allocatable  :: bb_heat_thin_integrand

  real(kind=dp), dimension(:), allocatable  :: bb_photo_thick_table
  real(kind=dp), dimension(:), allocatable  :: bb_heat_thick_table
  real(kind=dp), dimension(:), allocatable  :: bb_photo_thin_table
  real(kind=dp), dimension(:), allocatable  :: bb_heat_thin_table
  
  ! photrates contains all the photo-ionization rates and heating rates
  type photrates
     real(kind=dp) :: photo_cell       ! HI photoionization rate of the cell
     real(kind=dp) :: heat_cell        ! HI heating rate of the cell       
     real(kind=dp) :: photo_in         ! HI photoionization rate incoming to the cell    
     real(kind=dp) :: heat_in          ! HI heating rate incoming to the cell
     real(kind=dp) :: photo_out        ! HI photoionization rate outgoing from the cell
     real(kind=dp) :: heat_out         ! HI heating rate outgoing from the cell
  end type photrates

  ! tablepos helps to locate correct position of the photoionization and heating tables
  type tablepos
    real(kind=dp) :: tau            
    real(kind=dp) :: odpos          
    real(kind=dp) :: residual       
    integer       :: ipos           
    integer       :: ipos_p1        
  end type tablepos 

#ifdef MPI       
    integer,private :: ierror
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                                             !
  subroutine rad_ini ()                                                                      !
                                                                                             !
    ! Ask for the parameters of the spectrum                                                 !
    call spectrum_parms ()                                                                   !
                                                                                             !
    ! Initialize integration routines                                                        !
    call romberg_initialisation(NumFreq)                                                     !
                                                                                             !
    ! Determine spectrum diagnostics                                                         !
    call spec_diag ()                                                                        !
                                                                                             !
    ! Generate photoionization tables and heating tables                                     !
    call spec_integration ()                                                                 !
                                                                                             !
  end subroutine rad_ini                                                                     !
                                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Ask for the parameters of the spectrum
  subroutine spectrum_parms

    use file_admin, only: stdinput, file_input
    
    integer :: i_choice       ! option number
    real(kind=dp) :: total_flux  ! black body total flux

    ! Ask for the input if you are processor 0 and the
    ! spectral parameters are not set in the c2ray_parameters
    ! Note that it is assumed that if teff_nominal is set, 
    ! S_star_nominal is ALSO set.
    if (rank == 0 .and. T_eff_nominal == 0.0) then

      if (.not.file_input) write(*,'(A)') ' '
      T_eff = 0.0
      do while (T_eff < 2000.0 .or. T_eff > 200000.) 
        if (.not.file_input) write(*,'(A,$)') 'Give black body effective temperature: '
        read(stdinput,*) T_eff
        if (.not.file_input) write(*,*)
        if (T_eff < 2000.0 .or. T_eff > 200000.) then
          write(*,*) 'Error: Effective temperature out of range. Try again'
          write(*,*) 'Valid range: 2000 to 200,000'
        endif
      enddo

      ! Find total flux (Stefan-Boltzmann law)
      total_flux = sigma_SB*T_eff**4
       
      ! Ask for radius, luminosity, ionizing luminosity or ionizing photon rate?
      if (.not.file_input) then
        write(*,'(A)') ' '
        write(*,'(A)') 'You can specify' 
        write(*,'(A)') ' 1) a stellar radius'
        write(*,'(A)') ' 2) a luminosity'
        write(*,'(A)') ' 3) Total number of ionizing photons'
      endif

      i_choice = 0

      do while (i_choice <= 0 .or. i_choice > 3)
        if (.not.file_input) write(*,'(A,$)') 'Preferred option (1, 2 or 3): '
        read(stdinput,*) i_choice
        if (i_choice <= 0 .or. i_choice > 3) then
          write(*,*) 'Error: Choose between 1 2 or 3'
        endif
      enddo

      if (i_choice .eq. 1) then
        if (.not.file_input) write(*,'(A,$)') 'Give radius in solar radii: '
        read(stdinput,*) R_star
        R_star = R_star*r_solar
        L_star = R_star*R_star*(4.0d0*pi*total_flux)
        S_star = 0.0  ! Number of photo-ionizing photons is set to zero here but will be reset in spec_diag routine

      elseif (i_choice .eq. 2) then
        if (.not.file_input) write(*,'(A,$)') 'Give luminosity in solar luminosities: '
        read(stdinput,*) L_star
        L_star = L_star*l_solar
        R_star = dsqrt(L_star/(4.0d0*pi*total_flux))
        S_star = 0.0  ! Number of photo-ionizing photons is set to zero here but will be reset in spec_diag routine

      else
        if (.not.file_input) write(*,'(A,$)') 'Give S_* (ionizing photons s^-1): '
        read(stdinput,*) S_star
        ! Assign some fiducial values for R_star and L_star, 
        ! these are scaled to correspond to S_star in routine spec_diag
        R_star = r_solar
        L_star = R_star*R_star*(4.0d0*pi*total_flux)
      endif
    else
      ! teff and S_star are assumed to have been set in the c2ray_parameter module
      T_eff = T_eff_nominal
      S_star = S_star_nominal
      total_flux = sigma_SB*T_eff**4
      ! Assign some fiducial values for R_star and L_star, 
      ! these are scaled to correspond to S_star in routine spec_diag
      R_star = r_solar
      L_star = R_star*R_star*(4.0d0*pi*total_flux)

    endif

#ifdef MPI       
    ! Distribute the input parameters to the other nodes
    call MPI_BCAST(T_eff,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         ierror)
    call MPI_BCAST(R_star,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         ierror)
    call MPI_BCAST(L_star,1,MPI_DOUBLE_PRECISION,0, & 
         MPI_COMM_NEW,ierror)
    call MPI_BCAST(S_star,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_NEW, &
         ierror)
#endif
    
  end subroutine spectrum_parms

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! Calculates properties of the black body spectrum
  subroutine spec_diag ()

    integer :: i_freq
    real(kind=dp) :: h_over_kT, freq_max
    real(kind=dp) :: flux, S_star_unscaled,scaling
    real(kind=dp), dimension(0:NumFreq) :: frequency
    real(kind=dp), dimension(0:NumFreq) :: weight
    real(kind=dp), dimension(0:NumFreq) :: bb

    ! This is h/kT (unit 1/Hz, or sec)
    h_over_kT = hplanck/(k_B*T_eff)
    
    ! Two upper limit of frequency integration
    freq_max1 = 700.0*T_low*k_B/hplanck 
    freq_max2 = 5.88e-05*T_high*1e15

    ! Upper limit of frequency integration
    freq_max = min(freq_max1,10.0*freq_max2)
    
    ! Frequency step
    freq_step = (freq_max-ion_freq_HI)/real(NumFreq)

    ! Fill the arrays (frequency, weight, spectrum)
    do i_freq=0,NumFreq
      frequency(i_freq) = ion_freq_HI+freq_step*real(i_freq)
      weight(i_freq) = freq_step
      bb(i_freq) = two_pi_over_c_square*frequency(i_freq)*frequency(i_freq)/(exp(frequency(i_freq)*h_over_kT)-1.0)
    enddo

    ! Find flux by integrating
    flux = scalar_romberg(bb,weight,NumFreq,NumFreq,0)

    ! Find out what is the S_star for the radius supplied.
    S_star_unscaled = 4.0*pi*R_star*R_star*flux

    ! If S_star is zero, it is set here.
    if (S_star .eq. 0.0) then
       S_star = S_star_unscaled
    else
       ! Find out the factor by which to change the radius
       ! and luminosity to get the required S_star.
       scaling = S_star/S_star_unscaled
       R_star = sqrt(scaling)*R_star
       L_star = scaling*L_star
    endif
    
    ! Report back
    if (rank == 0) then
       write(logf,'(/a)')            ' Using a black body with'
       write(logf,'(a,1pe10.3,a)')   ' Teff=       ', T_eff, ' K'
       write(logf,'(a,1pe10.3,a)')   ' Radius=     ', R_star/r_solar, ' R_solar'
       write(logf,'(a,1pe10.3,a)')   ' Luminosity= ', L_star/l_solar, ' L_solar'
       write(logf,'(A,1PE10.3,A//)') ' Number of H ionizing photons: ', S_star,' s^-1'
    endif

  end subroutine spec_diag
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Calculates spectral integration cores
  subroutine spec_integration ()

    integer :: i_freq, i_tau
    real(kind=dp) :: freq_max
    real(kind=dp) :: R_star2, h_over_kT
    real(kind=dp), dimension(0:NumTau) :: tau
    real(kind=dp), dimension(0:NumFreq) :: frequency
    real(kind=dp), dimension(0:NumFreq) :: exp_factor
    real(kind=dp), dimension(0:NumTau) :: answer
    real(kind=dp), dimension(0:NumFreq,0:NumTau) :: weight
real(kind=dp), dimension(0:NumFreq) :: weight2, bb,bb2
real(kind=dp) :: mean_energy

    ! Photoionization integrand as a function of frequency and tau
    allocate(bb_photo_thick_integrand(0:NumFreq,0:NumTau))
    allocate(bb_photo_thin_integrand(0:NumFreq,0:NumTau))
  
    ! Heating integrand as a function of frequency and tau
    if (.not.isothermal) then
      allocate(bb_heat_thick_integrand(0:NumFreq,0:NumTau))
      allocate(bb_heat_thin_integrand(0:NumFreq,0:NumTau))
    endif

    ! Photoionization table as a function of tau
    allocate(bb_photo_thick_table(0:NumTau))
    allocate(bb_photo_thin_table(0:NumTau))
   
    ! Heating table as a function of heating tau
    if (.not.isothermal) then
      allocate(bb_heat_thick_table(0:NumTau))
      allocate(bb_heat_thin_table(0:NumTau))
    endif

    ! This is h/kT
    h_over_kT = hplanck/(k_B*T_eff)
    R_star2 = R_star*R_star

    ! fill the optical depth array used to fill the tables 
    ! it is filled in NumTau logarithmic steps 
    ! from minlogtau to maxlogtau
    do i_tau=1,NumTau
       tau(i_tau) = 10.0**(minlogtau+dlogtau*real(i_tau-1))
    enddo

    ! Position zero corresponds to zero optical depth
    tau(0)=0.0

    ! Warn about grey opacities:
    if (grey .and. rank == 0) write(logf,*) 'WARNING: Using grey opacities'

    if (ion_freq_HI.lt.freq_max1) then
       
      ! Upper limit of frequency integration
      freq_max = min(freq_max1,10.0*freq_max2)

! Find the mean photon energy in this frequency band
do i_freq=0,NumFreq
frequency(i_freq) = ion_freq_HI+freq_step*real(i_freq)
weight2(i_freq)=freq_step
bb(i_freq)=two_pi_over_c_square *frequency(i_freq)*frequency(i_freq)/(exp(frequency(i_freq)*h_over_kT)-1.0)
bb2(i_freq)=(frequency(i_freq)-ion_freq_HI)*bb(i_freq)
enddo
! Find flux by integrating
mean_energy=scalar_romberg(bb2,weight2,NumFreq,NumFreq,0)/scalar_romberg(bb,weight2,NumFreq,NumFreq,0) 
write(*,*) mean_energy

      do i_freq=0,NumFreq
         frequency(i_freq) = ion_freq_HI+freq_step*real(i_freq)
         ! Frequency dependence of the absorption cross section:
         if (grey) then
           exp_factor(i_freq) = 1.0
         else
           exp_factor(i_freq) = (frequency(i_freq)/ion_freq_HI)**(-pl_index_cross_section_HI)
         endif
      enddo

      ! Loop through the tau partition
      do i_tau=0,NumTau

        ! Loop through the frequency partition
        do i_freq=0,NumFreq
          weight(i_freq,i_tau)=freq_step
            
          ! Assign values to the photo integrands
          if (tau(i_tau)*exp_factor(i_freq) < 700.0) then
            bb_photo_thick_integrand(i_freq,i_tau) = 4.0*pi*R_star2*two_pi_over_c_square*frequency(i_freq)* &
                                                     frequency(i_freq)*exp(-tau(i_tau)*exp_factor(i_freq))/ &
                                                     (exp(frequency(i_freq)*h_over_kT)-1.0)
            bb_photo_thin_integrand(i_freq,i_tau) = 4.0*pi*R_star2*two_pi_over_c_square*frequency(i_freq)* &
                                                    frequency(i_freq)*exp_factor(i_freq)*exp(-tau(i_tau)* &
                                                    exp_factor(i_freq))/(exp(frequency(i_freq)*h_over_kT)-1.0)
          else
            bb_photo_thick_integrand(i_freq,i_tau) = 0.0
            bb_photo_thin_integrand(i_freq,i_tau) = 0.0
          endif

          ! Assign values to the heating integrands
          if (.not.isothermal) then
            !bb_heat_thick_integrand(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HI)* &
            !                                        bb_photo_thick_integrand(i_freq,i_tau)
            !bb_heat_thin_integrand(i_freq,i_tau) = hplanck*(frequency(i_freq)-ion_freq_HI)* &
            !                                       bb_photo_thin_integrand(i_freq,i_tau)

!bb_heat_thick_integrand(i_freq,i_tau) = hplanck*(1.1*ion_freq_HI)* &
!                                        bb_photo_thick_integrand(i_freq,i_tau)
!bb_heat_thin_integrand(i_freq,i_tau) = hplanck*(1.1*ion_freq_HI)* &
!                                       bb_photo_thin_integrand(i_freq,i_tau)
bb_heat_thick_integrand(i_freq,i_tau) = hplanck*mean_energy* &
                                        bb_photo_thick_integrand(i_freq,i_tau)
bb_heat_thin_integrand(i_freq,i_tau) = hplanck*mean_energy* &
                                       bb_photo_thin_integrand(i_freq,i_tau)
          endif

        enddo

      enddo

    endif
    
    ! Make photo tables
    call vector_romberg (bb_photo_thick_integrand,weight,NumFreq,NumFreq,NumTau,answer)
       bb_photo_thick_table = answer
    call vector_romberg (bb_photo_thin_integrand,weight,NumFreq,NumFreq,NumTau,answer)
       bb_photo_thin_table = answer
    
    ! Make heating tables
    if (.not.isothermal) then
       call vector_romberg (bb_heat_thick_integrand,weight,NumFreq,NumFreq,NumTau,answer)
          bb_heat_thick_table = answer
       call vector_romberg (bb_heat_thin_integrand,weight,NumFreq,NumFreq,NumTau,answer)
          bb_heat_thin_table = answer
    endif

    ! deallocate the useless photo integrand
    deallocate(bb_photo_thick_integrand)
    deallocate(bb_photo_thin_integrand)

    ! deallocate the useless heating integrand
    if (.not.isothermal) then
      deallocate(bb_heat_thick_integrand)
      deallocate(bb_heat_thin_integrand)
    endif

  end subroutine spec_integration

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! this subroutine calculates photo-ionization rates at a particular sets of column density
  subroutine photoion (phi,colum_in_HI,colum_out_HI,vol)
    
    ! result of the routine
    type(photrates), intent(out) :: phi

    ! Incoming HI column density
    real(kind=dp), intent(in) :: colum_in_HI 

    ! Outgoing HI column density
    real(kind=dp), intent(in) :: colum_out_HI 

    ! Volume of shell cell
    real(kind=dp), intent(in) :: vol 

    real(kind=dp) :: tau_in_HI, tau_out_HI
    type(tablepos) :: tau_pos_in, tau_pos_out
    
    ! find the optical depths (in and outgoing)
    tau_in_HI = sigma_HI_at_ion_freq*colum_in_HI
    tau_out_HI = sigma_HI_at_ion_freq*colum_out_HI
    
    ! find the table positions for the optical depth (ingoing)
    tau_pos_in%tau =  log10(max(1.0e-20_dp,tau_in_HI))
    tau_pos_in%odpos = min(real(NumTau,dp),max(0.0_dp,1.0+ &
                       (tau_pos_in%tau-minlogtau)/dlogtau))
    tau_pos_in%ipos = int(tau_pos_in%odpos)
    tau_pos_in%residual = tau_pos_in%odpos-real(tau_pos_in%ipos,dp)
    tau_pos_in%ipos_p1 = min(NumTau,tau_pos_in%ipos+1)
    
    ! find the table positions for the optical depth (outgoing)
    tau_pos_out%tau = log10(max(1.0e-20_dp,tau_out_HI))
    tau_pos_out%odpos = min(real(NumTau,dp),max(0.0_dp,1.0+ &
                        (tau_pos_out%tau-minlogtau)/dlogtau))
    tau_pos_out%ipos = int(tau_pos_out%odpos)
    tau_pos_out%residual = tau_pos_out%odpos-real(tau_pos_out%ipos)
    tau_pos_out%ipos_p1 = min(NumTau,tau_pos_out%ipos+1)

    ! Incoming, outcoming, current cell total photoionization rate
    phi%photo_in = (bb_photo_thick_table(tau_pos_in%ipos)+(bb_photo_thick_table(tau_pos_in%ipos_p1)- &
                   bb_photo_thick_table(tau_pos_in%ipos))*tau_pos_in%residual)

    ! Incoming, outcoming, current cell HI heating rate
    if (.not.isothermal) phi%heat_in = (bb_heat_thick_table(tau_pos_in%ipos)+(bb_heat_thick_table(tau_pos_in%ipos_p1)- &
                                       bb_heat_thick_table(tau_pos_in%ipos))*tau_pos_in%residual)

    ! When current cell is optically thick
    if (abs(tau_out_HI-tau_in_HI) .gt. 1e-2) then 
         
      phi%photo_out = (bb_photo_thick_table(tau_pos_out%ipos)+(bb_photo_thick_table(tau_pos_out%ipos_p1)- &
                      bb_photo_thick_table(tau_pos_out%ipos))*tau_pos_out%residual)
      phi%photo_cell = (phi%photo_in-phi%photo_out)/vol
       
      if (.not.isothermal) then
        phi%heat_out = (bb_heat_thick_table(tau_pos_out%ipos)+(bb_heat_thick_table(tau_pos_out%ipos_p1)- &
                       bb_heat_thick_table(tau_pos_out%ipos))*tau_pos_out%residual)
        phi%heat_cell = (phi%heat_in-phi%heat_out)/vol
      endif
       
    ! When current cell is optically thin
    else 
       
      phi%photo_cell = (tau_out_HI-tau_in_HI)*(bb_photo_thin_table(tau_pos_in%ipos)+ &
                       (bb_photo_thin_table(tau_pos_in%ipos_p1)- &
                       bb_photo_thin_table(tau_pos_in%ipos))*tau_pos_in%residual)/vol
      phi%photo_out = phi%photo_in-phi%photo_cell*vol

      if (.not.isothermal) then
        phi%heat_cell = (tau_out_HI-tau_in_HI)*(bb_heat_thin_table(tau_pos_in%ipos)+ &
                        (bb_heat_thin_table(tau_pos_in%ipos_p1)- &
                        bb_heat_thin_table(tau_pos_in%ipos))*tau_pos_in%residual)/vol
        phi%heat_out = phi%heat_in-phi%heat_cell*vol
      endif

    endif

  end subroutine photoion

end module radiation

