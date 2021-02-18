module mod_math

  use iso_fortran_env, only: int32, real32,real64
  implicit none

  private
  public :: erfc, asinh, Thot, Tcold

contains

  pure function erfc(x) result(res)
    real(real64), intent(in) :: x
    real(real64) :: res
    res = 1.0 - erf(x)
  end function erfc


 pure function asinh(x) result(res)
     real(real64), intent(in) :: x
     real(real64) :: res
     res=log(x+sqrt(x**2+1))
 end function asinh

 pure function Thot(time) result(res)
    real(real64), intent (in):: time
    real(real64) :: res

    if(time.le.9.9) then
       res=1.0
    else
        res=1.0
     endif
 end function Thot
 

pure function Tcold(time) result(res)
   real(real64), intent (in):: time
   real(real64) :: res

   if(time.le.9.9) then
      res=1.0
   else
       res=1.0
    endif
end function Tcold

  pure subroutine Gauss(nmax,ncell,A,B,C,F,bx,fx,phi)

     real(real64), intent(in out) :: phi(0:nmax)
     real(real64), intent(in out) :: A(1:nmax)
     real(real64), intent(in out) :: B(0:nmax)
     real(real64), intent(in out) :: C(0:nmax-1)
     real(real64), intent(in out) :: F(0:nmax)
     real(real64), intent(in out) :: bx(0:nmax)
     real(real64), intent(in out) :: fx(0:nmax)
     integer(int32), intent(in out) :: nmax
     integer(int32), intent(in out) :: ncell
     integer(int32) :: j

      j=0
      Bx(j) = B(j)
      Fx(j) = F(j)
    
      do j=1,ncell
         Bx(j) = B(j) - A(j)*C(j-1) /Bx(j-1)
         Fx(j) = F(j) - A(j)*Fx(j-1)/Bx(j-1)
      end do
    
      j = ncell
      phi(j) = Fx(j)/Bx(j)
    
      do j=ncell-1,0,-1
         phi(j) = ( Fx(j)-C(j)*phi(j+1) )/Bx(j)
      end do
  end subroutine Gauss

end module mod_math 

