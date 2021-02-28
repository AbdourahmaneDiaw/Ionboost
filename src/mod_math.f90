module mod_math

  use iso_fortran_env, only: int32, real32,real64
  implicit none

  private
  public :: Thot_, Tcold_, gauss

contains

 pure function Thot_(time) result(res)
    real(real64), intent (in):: time
    real(real64) :: res

    if(time.le.9.9) then
        res=1.0
    else
        res=1.0
     endif
 end function Thot_
 

pure function Tcold_(time) result(res)
   real(real64), intent (in):: time
   real(real64) :: res

   if(time.le.9.9) then
      res=1.0
   else
      res=1.0
    endif
end function Tcold_

  pure subroutine gauss(nmax,ncell,A,B,C,F,bx,fx,phi)

     real(real64), intent(in out) :: phi(0:nmax)
     real(real64), intent(in) :: A(1:nmax), C(0:nmax-1)
     real(real64), intent(in) :: B(0:nmax), F(0:nmax)
     real(real64), intent(in out) :: bx(0:nmax)
     real(real64), intent(in out) :: fx(0:nmax)
     integer(int32), intent(in) :: nmax
     integer(int32), intent(in) :: ncell
     integer(int32) :: j

      Bx(0) = B(0)
      Fx(0) = F(0)

      do concurrent(j = 1:ncell)
         Bx(j) = B(j) - A(j)*C(j-1) /Bx(j-1)
         Fx(j) = F(j) - A(j)*Fx(j-1)/Bx(j-1)
      end do
    
      phi(ncell) = Fx(ncell)/Bx(ncell)
    
     do concurrent(j=ncell-1:0:-1)
         phi(j) = ( Fx(j)-C(j)*phi(j+1) )/Bx(j)
      end do
  end subroutine gauss



end module mod_math 

