!
! Calcul de la correction mouillée selon le modèle de Wagner
!


function correction(t, w, list_param_c, list_param_y)
    implicit none
    real(kind=8), parameter :: PI = 4*atan(1.d0)
    !1ere valeur de list_param_c indique le type de géométrie
    !0 pour wedge
    !1 pour parabole
    ! list_param_y nécessaire pour y_commande
    real(kind=8), intent(in) :: t,w
    real(kind=8), dimension(2), intent(in) :: list_param_c, list_param_y
    real(kind=8) :: c_wedge, c_parabola
    real(kind=8) :: correction, y_commande

    if (list_param_c(1) == 0) then
        ! cor = pi h / (2 tan beta)
        c_wedge = PI * (-y_commande(t, list_param_y)-w) / (2.d0 * tan(list_param_c(2)))

    else if (list_param_c(1) == 1) then
        ! cor = sqrt( 4Rh )
        c_parabola = sqrt(4 * list_param_c(2) * (-y_commande(t,list_param_y)-w))

    end if
end function correction


function y_commande(t, list_param)
    implicit none
    !1ere valeur de list_param indique le type de loi de commande
    !0 pour vitesse constante
    real(kind=8), intent(in) :: t
    real(kind=8), dimension(2), intent(in) :: list_param
    real(kind=8) :: y_commande

    if(list_param(1) == 0) then
        y_commande = -list_param(2) * t !vitesse constante
    end if
end function y_commande 
