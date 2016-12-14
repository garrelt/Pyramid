module solid_angle

  use precision, only: dp

  implicit none

  real(kind=dp), parameter :: pi = 3.14159265359

contains

  real(kind=dp) function cc (xx,yy,zz,nn)

    implicit none

    integer, intent(in) :: xx,yy,zz,nn
    real :: x,y,z,n

    x = real(xx)
    y = real(yy)
    z = real(zz)
    n = real(nn)

    if (zz.ne.0) then 
      cc = acos(sqrt((n*n)/(x*x+n*n))*cos(atan(y/z)))
    else
      cc = acos(sqrt((n*n)/(x*x+n*n))*cos(pi/2))
    endif 

  end function cc


  real(kind=dp) function ss (xx,yy,zz,nn)

    implicit none

    integer, intent(in) :: xx,yy,zz,nn
    real :: x,y,z,n

    x = real(xx)
    y = real(yy)
    z = real(zz)
    n = real(nn)
  
    if (zz.ne.0) then   
      ss = asin(sqrt((n*n)/(x*x+n*n))*sin(atan(y/z)))
    else
      ss = asin(sqrt((n*n)/(x*x+n*n))*sin(pi/2))
    endif

  end function ss

end module solid_angle
