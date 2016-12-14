module trilinear

  use precision, only: dp
  use input, only : level_py, level_sl
  use array, only : interpolated_HI_photoionization_rate_array
  use array, only : interpolated_HeI_photoionization_rate_array  
  use array, only : interpolated_HeII_photoionization_rate_array  
  use array, only : interpolated_photoheating_rate_array
  use array, only : long_global_HI_photoionization_rate_array
  use array, only : long_global_HeI_photoionization_rate_array  
  use array, only : long_global_HeII_photoionization_rate_array  
  use array, only : long_global_photoheating_rate_array 

  contains

  subroutine trilinear_interpolation ()

    implicit none
  
    integer :: i,j,k
    real(kind=dp) :: i_center, j_center, k_center ! tuned to the sl frame
    integer :: i_small, j_small, k_small
    integer :: i_large, j_large, k_large
    real(kind=dp) :: vol_sss, vol_ssl, vol_sls, vol_sll
    real(kind=dp) :: vol_lss, vol_lsl, vol_lls, vol_lll
    real(kind=dp) :: vol

    do i = 1, 2*level_py
      do j = 1, 2*level_py
        do k = 1, 2*level_py
   
          i_center = (i-0.5)*level_sl/(2.0*level_py)+0.5
          j_center = (j-0.5)*level_sl/(2.0*level_py)+0.5
          k_center = (k-0.5)*level_sl/(2.0*level_py)+0.5
          i_small = int(i_center)
          j_small = int(j_center)
          k_small = int(k_center)
          i_large = i_small+1
          j_large = j_small+1
          k_large = k_small+1

          vol_sss = abs(i_center-i_small) * abs(j_center-j_small) * abs(k_center-k_small)
          vol_ssl = abs(i_center-i_small) * abs(j_center-j_small) * abs(k_center-k_large)
          vol_sls = abs(i_center-i_small) * abs(j_center-j_large) * abs(k_center-k_small)
          vol_sll = abs(i_center-i_small) * abs(j_center-j_large) * abs(k_center-k_large)
          vol_lss = abs(i_center-i_large) * abs(j_center-j_small) * abs(k_center-k_small)
          vol_lsl = abs(i_center-i_large) * abs(j_center-j_small) * abs(k_center-k_large)
          vol_lls = abs(i_center-i_large) * abs(j_center-j_large) * abs(k_center-k_small)
          vol_lll = abs(i_center-i_large) * abs(j_center-j_large) * abs(k_center-k_large)

          interpolated_HI_photoionization_rate_array (i,j,k) = &
            long_global_HI_photoionization_rate_array (i_small, j_small, k_small) * vol_sss + &
            long_global_HI_photoionization_rate_array (i_small, j_small, k_large) * vol_ssl + &
            long_global_HI_photoionization_rate_array (i_small, j_large, k_small) * vol_sls + &
            long_global_HI_photoionization_rate_array (i_small, j_large, k_large) * vol_sll + &
            long_global_HI_photoionization_rate_array (i_large, j_small, k_small) * vol_lss + &
            long_global_HI_photoionization_rate_array (i_large, j_small, k_large) * vol_lsl + &
            long_global_HI_photoionization_rate_array (i_large, j_large, k_small) * vol_lls + &
            long_global_HI_photoionization_rate_array (i_large, j_large, k_large) * vol_lll

          interpolated_HeI_photoionization_rate_array (i,j,k) = &
            long_global_HeI_photoionization_rate_array (i_small, j_small, k_small) * vol_sss + &
            long_global_HeI_photoionization_rate_array (i_small, j_small, k_large) * vol_ssl + &
            long_global_HeI_photoionization_rate_array (i_small, j_large, k_small) * vol_sls + &
            long_global_HeI_photoionization_rate_array (i_small, j_large, k_large) * vol_sll + &
            long_global_HeI_photoionization_rate_array (i_large, j_small, k_small) * vol_lss + &
            long_global_HeI_photoionization_rate_array (i_large, j_small, k_large) * vol_lsl + &
            long_global_HeI_photoionization_rate_array (i_large, j_large, k_small) * vol_lls + &
            long_global_HeI_photoionization_rate_array (i_large, j_large, k_large) * vol_lll

          interpolated_HeII_photoionization_rate_array (i,j,k) = &
            long_global_HeII_photoionization_rate_array (i_small, j_small, k_small) * vol_sss + &
            long_global_HeII_photoionization_rate_array (i_small, j_small, k_large) * vol_ssl + &  
            long_global_HeII_photoionization_rate_array (i_small, j_large, k_small) * vol_sls + &  
            long_global_HeII_photoionization_rate_array (i_small, j_large, k_large) * vol_sll + &  
            long_global_HeII_photoionization_rate_array (i_large, j_small, k_small) * vol_lss + &  
            long_global_HeII_photoionization_rate_array (i_large, j_small, k_large) * vol_lsl + &  
            long_global_HeII_photoionization_rate_array (i_large, j_large, k_small) * vol_lls + &  
            long_global_HeII_photoionization_rate_array (i_large, j_large, k_large) * vol_lll  

          interpolated_photoheating_rate_array (i,j,k) = &
            long_global_photoheating_rate_array (i_small, j_small, k_small) * vol_sss + &
            long_global_photoheating_rate_array (i_small, j_small, k_large) * vol_ssl + & 
            long_global_photoheating_rate_array (i_small, j_large, k_small) * vol_sls + & 
            long_global_photoheating_rate_array (i_small, j_large, k_large) * vol_sll + & 
            long_global_photoheating_rate_array (i_large, j_small, k_small) * vol_lss + & 
            long_global_photoheating_rate_array (i_large, j_small, k_large) * vol_lsl + & 
            long_global_photoheating_rate_array (i_large, j_large, k_small) * vol_lls + & 
            long_global_photoheating_rate_array (i_large, j_large, k_large) * vol_lll 

        enddo
      enddo
    enddo

  end subroutine trilinear_interpolation	
	
end module trilinear
