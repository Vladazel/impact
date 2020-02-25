!
! Programme principale pour la résolution du mouvement
!

program Resolution_equadiff
    ! Résout le problème aux dérivées partielles
    implicit none
    real(kind=8), parameter :: PI = 4*atan(1.d0)
    real(kind=8) :: dt, t
    real(kind=8), dimension(2) :: w, wp
    real(kind=8), external :: f_equa
    integer :: i, nsteps

    open(unit = 1, file = './params.inp', status = 'old')
    read(1,*) dt
    read(1,*) nsteps
    close(1)

    w = (/0.d0, 0.d0/)
    wp = (/0.d0, 0.d0/)

    t = 0.d0

    open(unit = 2, file = './resultats.dat', status = 'replace', position = 'append', action = 'write')

    do i = 1,nsteps
        call rk4(f_equa, dt, t, w, wp)
        t = t+dt
    end do

    close(2)
end program Resolution_equadiff

subroutine rk4(f, dt, t, w, wp)
    ! Implémentation de Runge-Kutta d'ordre 4
    implicit none
    real(kind=8), intent(in)  :: dt, t
    real(kind=8), dimension(2), intent(out) :: w, wp
    real(kind=8) :: k1, k2, k3, k4
    real(kind=8), external :: f

    k1 = f(t,        w,                      wp)
    k2 = f(t+dt/2.0, w+dt/2.0*wp,            wp+dt/2.0*k1)
    k3 = f(t+dt/2.0, w+dt/2.0*wp+dt**2/4*k1, wp+dt/2.0*k2)
    k4 = f(t+dt,     w+dt*wp+dt**2/2.0*k2,   wp+dt*k3)

    w(1) = w(1) + dt*w(2) + dt**2/6.0*(k1+k2+k3)
    w(2) = w(2) + dt/6.0*(k1+2.0*k2+2.0*k3+k4)
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
