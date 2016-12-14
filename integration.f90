module integration
	
  use precision, only: dp

contains
	
  real function volume_integration(ax,bx,cx,dx,ay,by,cy,dy,az,cz,example,case_number)

    real(kind=dp), intent(in) :: ax
    real(kind=dp), intent(in) :: bx
    real(kind=dp), intent(in) :: cx 
    real(kind=dp), intent(in) :: dx
    real(kind=dp), intent(in) :: ay
    real(kind=dp), intent(in) :: by
    real(kind=dp), intent(in) :: cy
    real(kind=dp), intent(in) :: dy
    real(kind=dp), intent(in) :: az
    real(kind=dp), intent(in) :: cz	  
    integer, intent(in) :: example
    integer, intent(in) :: case_number
    volume_integration = (cx-ax)*(cy-ay)*(cz-az) + &
                         (1.0/2.0)*((cx-ax)*(dy-by)+(cy-ay)*(dx-bx))*(cz*cz-az*az) + &
                         (1.0/3.0)*(dx-bx)*(dy-by)*(cz*cz*cz-az*az*az)	     
  !if (volume_integration.lt.-1e-2) then
  !  write(*,*) 'example',example,'case_number',case_number
  !  write(*,*) ax,bx,cx,dx,ay,by,cy,dy,az,cz,'and volume is',volume_integration 
  !endif
  
  end function volume_integration
	
end module integration
