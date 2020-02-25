!
! Programme principale pour la résolution du mouvement
!

program Resolution_equadiff
    ! Résout le problème aux dérivées partielles
    implicit none
    real(kind=8), parameter :: PI = 4*atan(1.d0)
    real(kind=8) :: dt, t
    ! w0 : solution au pas de temps actuel
    ! w : solution au prochain pas de temps
    real(kind=8), dimension(2) :: w0, w
    external :: f_equa
    integer :: i, nsteps

    open(unit = 1, file = './params.inp', status = 'old')
    read(1,*) dt
    read(1,*) nsteps
    close(1)

    w0 = (/0.d0, 0.d0/)

    t = 0.d0

    open(unit = 2, file = './resultats.dat')

    do i = 1,nsteps
        call rk4(f_equa, dt, t, w0, w)
        write(2,*) t, w(1), w(2)
        t = t+dt
    end do

    close(2)
end program Resolution_equadiff

subroutine rk4(f, dt, t0, w0, w)
    ! Implémentation de Runge-Kutta d'ordre 4
    implicit none
    real(kind=8), intent(in) :: dt, t0
    real(kind=8), dimension(2), intent(in) :: w0
    ! w0 : solution au pas de temps actuel
    ! w  : solution au prochain pas de temps
    real(kind=8), dimension(2), intent(out) :: w
    real(kind=8) :: t1, t2, t3
    real(kind=8), dimension(2) :: k0, k1, k2, k3
    real(kind=8), dimension(2) :: w1, w2, w3
    external :: f

    call f(t0, w0, k0)

    t1 = t0 + dt / 2.0d0
    w1 = w0 + dt / 2.0d0 * k0
    call f(t1, w1, k1)

    t2 = t0 + dt / 2.0d0
    w2 = w0 + dt / 2.0d0 * k1
    call f(t2, w2, k2)

    t3 = t0 + dt / 2.0d0
    w3 = w0 + dt / 2.0d0 * k2
    call f(t3, w3, k3)

    w = w0 + dt / 6.0d0 * (k0 + 2.0d0*k1 + 2.0d0*k2 + k3)
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
