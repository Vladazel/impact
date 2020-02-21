!
! Contient deux fonctions :
!
! deriv : dérivée première dy/dx
! deriv2 : dérivée seconde d2y/dx2
!

function deriv(dx, yref, yplus)
    ! calculates differentiation in x
    ! y(x) = yref
    ! y(x+h) = yplus
    ! with formula y'=(yplus-y)/dx
    implicit none
    real(kind=8), intent(in) :: dx, yref, yplus
    real(kind=8) :: deriv

    deriv = (yplus - yref) / dx
end function deriv

function deriv2(dx, yref, yplus, ymoins)
    ! calculates second differentiation in x
    ! y(x) = yref
    ! y(x+h) = yplus
    ! y(x-h) = ymoins
    implicit none
    real(kind=8), intent(in) :: dx, yref, yplus, ymoins
    real(kind=8) :: deriv2
    
    deriv2 = (yplus + ymoins - 2.0*yref) / (dx**2)
end function deriv2
