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
    integer :: i, nsteps
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
        call y_commande(t, ycom, tmp, tmp, param_com) 
        y = ycom + w0(1)
        Ma = added_mass(t, y, param_geom)
        Cs = slamming_coef(t, y, param_geom)
        write(3,*) t*sqrt(k/M), w0(1)/param_com(2)*sqrt(k/M), Cs*param_com(2)/sqrt(k*M), Ma/M

        call rk4(f_equa, dt, t, w0, w)
        t = t + dt
        w0 = w
    end do

    close(3)
end program Resolution_equadiff


