!>
!! \brief This module contains radiation-related physical constants
!!
!! Module for Capreole / C2-Ray (f90)
!!
!! \b Author: Garrelt Mellema
!!
!! \b Date:
!!
!! \b Version: cgs units
!!

module cgsphotoconstants

  ! A collection of physical constants and conversion factors for 
  ! photo-ionization calculations
  ! Units: cgs
  
  use precision, only: dp
  use cgsconstants

  !	!> Helium ionization potentials (eV)
  !	real(kind=dp), dimension(0:1),parameter :: ethe=(/24.587,54.416/)
  !> Hydrogen cross section
  real(kind=dp), parameter :: sigma_HI_at_ion_freq = 6.30e-18
  !> H ionization energy in frequency
  real(kind=dp), parameter :: ion_freq_HI=ev2fr*eth0
  !> Frequency dependence of H cross section parameter
  real(kind=dp), parameter :: pl_index_cross_section_HI=2.8
  !> Frequency dependence of He cross section parameter
  real(kind=dp), parameter :: pl_index_cross_section_HeI=1.7
  !> Frequency dependence of He+ cross section parameter
  real(kind=dp), parameter :: pl_index_cross_section_HeII=2.8

  !> maximum T_eff for black body
  real(kind=dp), parameter :: T_high = 200000.0
  !> minimum T_eff for black body
  real(kind=dp), parameter :: T_low = 2000.0

  !> H optical depth fit parameter for frequency range 2
  real(kind=dp) :: tf2h
  !> H optical depth fit parameter for frequency range 3
  real(kind=dp) :: tf3h
  !> He^0 optical depth fit parameter for frequency range 3
  real(kind=dp) :: tf3he0


end module cgsphotoconstants




