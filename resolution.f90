!
! Programme principale pour la résolution du mouvement
!

program Resolution_equadiff
    ! Résout le problème aux dérivées partielles
    implicit none
    real(kind=8), parameter :: PI = 4*atan(1.d0)
    real(kind=8) :: dt, tn, wn, wpn
    integer :: i, nsteps
    real(kind=8), dimension(2) :: param_geom, param_com
    real(kind=8) :: t_adim, w_adim

    open(unit = 1, file = './params.inp', status = 'old')
    read(1,*) dt
    read(1,*) nsteps
    close(1)

    wn = 0
    wpn = 0

    ! Résolution
    dt = 1e-5
    tn = 0.d0

    open(unit = 2, file = './resultats.dat', status = 'replace', position = 'append', action = 'write')

    do i = 1,nsteps
        call rk4(dt, tn, wn, wpn)
        tn = tn+dt
    end do

    close(2)
end program Resolution_equadiff

subroutine rk4(dt, tn, wn, wpn)
    ! Implémentation de Runge-Kutta d'ordre 4
    implicit none
    real(kind=8), intent(in)  :: dt, tn
    real(kind=8), intent(out) :: wn, wpn
    real(kind=8) :: k1, k2, k3, k4
    real(kind=8) :: f_equa

    k1 = f_equa(tn,        wn,                       wpn)
    k2 = f_equa(tn+dt/2.0, wn+dt/2.0*wpn,            wpn+dt/2.0*k1)
    k3 = f_equa(tn+dt/2.0, wn+dt/2.0*wpn+dt**2/4*k1, wpn+dt/2.0*k2)
    k4 = f_equa(tn+dt,     wn+dt*wpn+dt**2/2.0*k2,   wpn+dt*k3)

    wn = wn + dt*wpn + dt**2/6.0*(k1+k2+k3)
    wpn = wpn + dt/6.0*(k1+2.0*k2+2.0*k3+k4)
end subroutine rk4

subroutine y_commande(t, ycom, dycomdt, d2ycomdt2, list_param)
    ! Définit la loi de commande
    implicit none
    !1ere valeur de list_param indique le type de loi de commande
    !0 pour vitesse constante
    real(kind=8), intent(in) :: t
    real(kind=8), dimension(2), intent(in) :: list_param
    real(kind=8), intent(out) :: ycom, dycomdt, d2ycomdt2

    if(list_param(1) == 0) then
        ycom = -list_param(2) * t !vitesse constante
        dycomdt = -list_param(2)
        d2ycomdt2 = 0
    else
        error stop
    end if
end subroutine y_commande 
