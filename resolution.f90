!
! Programme principale pour la résolution du mouvement
!

program Resolution_equadiff
    ! Résout le problème aux dérivées partielles
    implicit none
    real(kind=8), parameter :: PI = 4*atan(1.d0)
    real(kind=8) :: dt, tstart, tstop, t, ycom, y
    ! w0 : solution au pas de temps actuel
    ! w : solution au prochain pas de temps
    real(kind=8), dimension(2) :: w0, w
    real(kind=8), dimension(2) :: param_com, param_geom
    real(kind=8) :: Ma, Cs, added_mass, slamming_coef
    real(kind=8) :: M, k
    integer :: i, nsteps, schema
    external :: f_equa
    !variable inutile qui sert juste pour le calcul de ycom car on veut pas les autres
    !valeurs que peut sortir la subroutine y_commande
    real(kind=8) :: tmp

    open(unit = 1, file = './param_num.inp', status = 'old')
    read(1,*) w0(1)
    read(1,*) w0(2)
    read(1,*) tstart
    read(1,*) tstop
    read(1,*) nsteps
    read(1,*) schema
    close(1)

    open(unit = 2, file = './param_phy.inp', status = 'old')
    read(2,*) M
    read(2,*) k
    read(2,*) param_com(1)
    read(2,*) param_com(2)
    read(2,*) param_geom(1)
    read(2,*) param_geom(2)
    close(2)

    t = tstart
    dt = (tstop - tstart)/nsteps
    print*, 'Pas de temps =', dt

    open(unit = 3, file = './resultats.dat')

    do i = 1,nsteps
        !On récupère la valeur de la commande pour calculer Ma et Cs
        !On utilise une variable inutile tmp car ici on n'utilise pas les valeurs des dérivées de y
        call y_commande(t, ycom, tmp, tmp, param_com)
        y = ycom + w0(1)
        Ma = added_mass(t, y, param_geom)
        Cs = slamming_coef(t, y, param_geom)
        write(3,*) t*sqrt(k/M), w0(1)/param_com(2)*sqrt(k/M), Cs*param_com(2)/sqrt(k*M), Ma/M

        call integration(f_equa, dt, t, w0, w, schema)
        t = t + dt
        w0 = w
    end do

    close(3)
end program Resolution_equadiff

subroutine integration(f, dt, t, w0, w, choix)
    real(kind=8) :: dt, t
    real(kind=8), dimension(2) :: w0, w
    integer, intent(in) :: choix
    external :: f, rk4

    if (choix == 0) then
        call rk4(f, dt, t, w0, w)
    else if (choix == 1) then
        call rkf45(f, dt, t, w0, w)
    else
        print*, "Erreur choix du schéma d'intégration"
        stop
    end if
end subroutine integration

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
