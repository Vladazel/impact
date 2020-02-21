!
! Définition de l'équa diff de la forme y'' = f(t, y, y')
! avec f = f_equa
! 
! Équation diff est : (M+Ma) y'' + Cs y' + k w = 0
! avec y = ycom + w
!

function f_equa(t, w, wp, M, k, param_geom, param_com)
    !param geom est une liste def le type de géométrie
    !Ref à wetted_correction.f90
    implicit none
    real(kind=8), dimension(3), intent(in) :: t !tn-2, tn-1, tn
    real(kind=8), dimension(2), intent(in) :: w !wn-1, wn
    real(kind=8), intent(in) :: wp, M, k
    real(kind=8), dimension(2), intent(in) :: param_geom
    real(kind=8), dimension(2), intent(in) :: param_com
    real(kind=8) :: f_equa
    real(kind=8) :: y_commande
    real(kind=8) :: yref, yplus, ymoins
    real(kind=8) :: wref, wplus
    real(kind=8) :: dt, dycomdt, dydt, d2ycomdt2
    real(kind=8) :: deriv, deriv2
    real(kind=8) :: Ma, Cs, added_mass, slamming_coef
    
    yref   = y_commande(t(2), param_com) 
    yplus  = y_commande(t(3), param_com) 
    ymoins = y_commande(t(1), param_com) 

    dt = t(3) - t(2)

    dycomdt = deriv(dt, yref, yplus)
    dydt    = dycomdt + wp

    d2ycomdt2 = deriv2(dt, yref, yplus, ymoins)

    Ma =    added_mass(t(3), w(2), param_geom, param_com)
    Cs = slamming_coef(t(2:3), w, param_geom, param_com)

    f_equa = 1/(M+Ma) * (Cs * dydt**2 - k*w(2)) - d2ycomdt2
end function f_equa

function added_mass(t, w, param_c, param_ycom)
    !Calcule la masse ajoutée
    implicit none
    real(kind=8), parameter :: RHO = 999
    real(kind=8), parameter :: PI = 4*atan(1.d0)
    real(kind=8), intent(in) :: t, w
    real(kind=8), dimension(2), intent(in) :: param_c
    real(kind=8), dimension(2), intent(in) :: param_ycom
    real(kind=8) :: added_mass
    real(kind=8) :: correction

    added_mass = 0.5 * RHO * PI * correction(t, w, param_c, param_ycom) ** 2
end function added_mass

function slamming_coef(t, w, param_c, param_ycom)
    !Calcule le coef de slamming Cs
    implicit none
    real(kind=8), parameter :: RHO = 999
    real(kind=8), parameter :: PI = 4*atan(1.d0)
    !on a besoin de deux pas de temps pour la dérivée première de c
    real(kind=8), dimension(2), intent(in) :: t, w 
    real(kind=8), dimension(2), intent(in) :: param_c
    real(kind=8), dimension(2), intent(in) :: param_ycom
    real(kind=8) :: slamming_coef, correction
    real(kind=8) :: y_commande
    real(kind=8) :: yref, ymoins, cref, cmoins
    real(kind=8) :: dcdy, deriv

    yref   = y_commande(t(2), param_ycom) + w(2) !valeur de y pour pas n
    ymoins = y_commande(t(1), param_ycom) + w(1) !valeur de y pour pas n-1

    cref   = correction(t(2), w(2), param_c, param_ycom) !valeur de c pour pas n
    cmoins = correction(t(1), w(1), param_c, param_ycom) !valeur de c pour pas n-1

    dcdy = deriv(yref-ymoins, cmoins, cref)

    !Ajout d'un moins pour avoir un coefficient positif
    slamming_coef = - RHO * PI * cref * dcdy
end function
