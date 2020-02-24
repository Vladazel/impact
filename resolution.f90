!
! Programme principale pour la résolution du mouvement
!

program Resolution_equadiff
    ! Résout le problème aux dérivées partielles
    implicit none
    real(kind=8), parameter :: PI = 4*atan(1.d0)
    real(kind=8) :: M, k, U, beta, rayon
    real(kind=8) :: wpn
    real(kind=8), dimension(2) :: wn
    real(kind=8), dimension(3) :: tn
    integer :: i, nsteps
    real(kind=8) :: dt, deriv
    real(kind=8), dimension(2) :: param_geom0, param_geom1, param_com
    real(kind=8) :: slamming_coef, cs, added_mass, ma
    real(kind=8) :: t_adim, w_adim

    nsteps = 5e3
   
    U = 5.d0
    M = 100.d0
    k = 10000000.d0
    beta = 0.5d0
    rayon = 1.d0

    wn(1) = U/(2.d0*PI)*(cos(2.d0*PI*dt)-1)
    wn(2) = U/(2.d0*PI)*(cos(2.d0*PI*2*dt)-1)
    wpn = deriv(dt, wn(1), wn(2))

    param_geom0 = (/0.d0, beta/)
    param_geom1 = (/1.d0, rayon/)
    param_com  = (/0.d0, U/)

    ! résolution  pour vitesse constante et wedge:

    open(1, file='./resultats_parabole.dat')
    
    ! Calculs pour la parabole
    dt = 1e-5
    tn = (/0.d0, dt, 2.d0*dt/)

    do i = 1,nsteps
        call rk4(dt, tn, wn, wpn, M, k, param_geom1, param_com)
        ma = added_mass(tn(3), wn(2),   param_geom1, param_com)
        cs = slamming_coef(tn(2:3), wn, param_geom1, param_com)
        t_adim = tn(3)*sqrt(k/M)
        w_adim = wn(2)*U*sqrt(M/k)
        write(1,*) t_adim, w_adim, cs/k, ma/k
    end do

    close(1)

    ! Calculs pour le dièdre
    dt = 1e-5
    tn = (/0.d0, dt, 2.d0*dt/)
    
    open(2, file='./resultats_wedge.dat')
    
    do i = 1,nsteps
        call rk4(dt, tn, wn, wpn, M, k, param_geom0, param_com)
        ma = added_mass(tn(3), wn(2),   param_geom0, param_com)
        cs = slamming_coef(tn(2:3), wn, param_geom0, param_com)
        t_adim = tn(3)*sqrt(k/M)
        w_adim = wn(2)*U*sqrt(M/k)
        write(2,*) t_adim, w_adim, cs/k, ma/k
    end do

    close(2)
end program Resolution_equadiff

subroutine rk4(dt, tn, wn, wpn, M, k, param_geom, param_com)
    ! Implémentation de Runge-Kutta d'ordre 4
    implicit none
    real(kind=8), intent(in)  :: dt
    real(kind=8), dimension(3), intent(out) :: tn
    real(kind=8), dimension(2), intent(out) :: wn
    real(kind=8), intent(out) :: wpn
    real(kind=8), intent(in) :: M, k
    real(kind=8), dimension(2), intent(in) :: param_geom, param_com
    real(kind=8) :: k1, k2, k3, k4
    real(kind=8) :: f_equa

    k1 = f_equa(tn,        wn,                       wpn,           M, k, param_geom, param_com)
    k2 = f_equa(tn+dt/2.0, wn+dt/2.0*wpn,            wpn+dt/2.0*k1, M, k, param_geom, param_com)
    k3 = f_equa(tn+dt/2.0, wn+dt/2.0*wpn+dt**2/4*k1, wpn+dt/2.0*k2, M, k, param_geom, param_com)
    k4 = f_equa(tn+dt,     wn+dt*wpn+dt**2/2.0*k2,   wpn+dt*k3,     M, k, param_geom, param_com)

    wn(1) = wn(2)
    wn(2) = wn(2) + dt*wpn + dt**2/6.0*(k1+k2+k3)

    wpn = wpn + dt/6.0*(k1+2.0*k2+2.0*k3+k4)

    tn(1) = tn(2)
    tn(2) = tn(3)
    tn(3) = tn(3) + dt
end subroutine rk4

function y_commande(t, list_param)
    ! Définit la loi de commande
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
