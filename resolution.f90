!
! Programme principale pour la résolution du mouvement
!

program Resolution_equadiff
    ! Résout le problème aux dérivées partielles
    implicit none
    real(kind=8), parameter :: PI = 4*atan(1.d0)
    real(kind=8) :: dt, tstart, tstop, t,t_out, ycom, y,relerr, abserr
    ! w0 : solution au pas de temps actuel
    ! w : solution au prochain pas de temps
    real(kind=8), dimension(2) :: w0, w, wp
    real(kind=8), dimension(2) :: param_com, param_geom
    real(kind=8) :: Ma, Cs, added_mass, slamming_coef
    real(kind=8) :: M, k
    integer(kind=4) :: i, nstep, schema, flag,neqn
    external :: f
    !variable inutile qui sert juste pour le calcul de ycom car on veut pas les autres
    !valeurs que peut sortir la subroutine y_commande
    real(kind=8) :: tmp, dycomdt

    open(unit = 1, file = './param_num.inp', status = 'old')
    read(1,*) w0(1)
    read(1,*) w0(2)
    read(1,*) tstart
    read(1,*) tstop
    read(1,*) nstep
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

    open(unit = 10, file = './resultats.dat')
    open(unit = 11, file = './energie.dat')
    open(unit = 12, file = './oscillateur.dat')

    if (schema == 0) then !rk4
        
        t = tstart
        dt = (tstop - tstart)/nstep
        print*, 'rk4'
        print*, 'Pas de temps =', dt
        
        call y_commande(t, ycom, dycomdt, tmp, param_com)
        y = ycom + w0(1)
        Ma = added_mass(t, y, param_geom)
        Cs = slamming_coef(t, y, param_geom)

        call ecrire(t, w, k, M, Ma, Cs, param_com, param_geom, ycom, dycomdt)

        do i = 1, nstep        

            !On récupère la valeur de la commande pour calculer Ma et Cs
            !On utilise une variable inutile tmp car ici on n'utilise pas les valeurs des dérivées de y
            call y_commande(t, ycom, dycomdt, tmp, param_com)
            y = ycom + w(1)
            Ma = added_mass(t, y, param_geom)
            Cs = slamming_coef(t, y, param_geom)

            call rk4(f, dt, t, w0, w)

            call ecrire(t, w, k, M, Ma, Cs, param_com, param_geom, ycom, dycomdt)

            t = t + dt
            w0 = w
        
        end do
    else if (schema == 1) then !rkf45
        print*, 'rkf45'
        abserr = sqrt ( epsilon ( abserr ) )
        relerr = sqrt ( epsilon ( relerr ) )

        flag = 1
       

        t = 0.0D+00
        t_out = 0.0D+00

        w(1) = w0(1)
        w(2) = w0(2)
        
        call f( t, w, wp )
        
        call y_commande(t, ycom, dycomdt, tmp, param_com)
        y = ycom + w(1)
        Ma = added_mass(t, y, param_geom)
        Cs = slamming_coef(t, y, param_geom)
        
        call ecrire(t, w, k, M, Ma, Cs, param_com, param_geom, ycom, dycomdt)

        do i = 1, nstep
            !définition de l'intervalle de temps pour l'intégration 
            t = ( real ( nstep - i + 1, kind = 8 ) * tstart &
                + real (         i - 1, kind = 8 ) * tstop ) &
                / real ( nstep,         kind = 8 )

            t_out = ( real ( nstep - i, kind = 8 ) * tstart &
                    + real (         i, kind = 8 ) * tstop ) &
                    / real ( nstep,     kind = 8 )
            neqn=2 !equation vectorielle de dimension 2
            call r8_rkf45 ( f, neqn, w, wp, t, t_out, relerr, abserr, flag )
            call y_commande(t, ycom, dycomdt, tmp, param_com)
            y = ycom + w(1)
            Ma = added_mass(t, y, param_geom)
            Cs = slamming_coef(t, y, param_geom)
        
            call ecrire(t, w, k, M, Ma, Cs, param_com, param_geom, ycom, dycomdt)
                    
        end do

    else
        print*, "Schéma d'intégration non défini"
        stop
        
    close(10)
    close(11)

    end if
        
end program Resolution_equadiff

subroutine carac_oscil(w, k, M, Ma, Cs, param_com, param_geom, w0, zeta)
    implicit none
    real(kind=8), parameter :: RHO = 999.d0
    real(kind=8), parameter :: PI = 4.d0*atan(1.d0)
    real(kind=8), intent(in) :: k, M, Ma, Cs
    real(kind=8), dimension(2), intent(in) :: w, param_com, param_geom
    real(kind=8) :: dCsdw
    real(kind=8), intent(out) :: w0, zeta

    if (param_geom(1) == 0.d0) then !dièdre
        dCsdw = 2 * RHO * PI * param_geom(2) * (param_com(1)/w(2) - 1)
    else if (param_geom(1) == 1.d0) then !parabole
        dCsdw = RHO * PI**3 / (4.d0*tan(param_geom(2))**2) * (param_com(1)/w(2) - 1)
    else
        print*, "Géométrie non définie"
        stop
    end if

    w0 = sqrt( (k+dCsdw*param_com(2)**2) / (M+Ma) )
    zeta = - (Cs*param_com(2)) / (M+Ma) / w0
end subroutine carac_oscil

subroutine ecrire(t, w, k, M, Ma, Cs, param_com, param_geom, ycom, dycomdt)
    implicit none
    real(kind=8), parameter :: RHO = 999.d0
    real(kind=8), intent(in) :: t, k, M, Ma, Cs, ycom, dycomdt
    real(kind=8), dimension(2), intent(in) :: w, param_com, param_geom
    external :: carac_oscil
    real(kind=8) :: w0, zeta

    write(10,*) t*param_com(2)*(RHO/M)**(1.d0/3.d0), &
                w(1)*(RHO/M)**(1.d0/3.d0), &
                Cs/M*(M/RHO)**(1.d0/3.d0), &
                Ma/M, &
                dycomdt/param_com(2), &
                w(2)/param_com(2)
                
    write(11,*) t*param_com(2)*(RHO/M)**(1.d0/3.d0), &
                0.5d0*M*(dycomdt+w(2))**2 , & !Ec(t)-Ec(0)
                0.5d0*Ma*(dycomdt+w(2))**2, &  !Ec Masse ajoutée
                0.5d0*k*w(1)**2, & !Energie potentielle élastique
                Cs*(dycomdt+w(2))**2*(w(1)+ycom) + 0.5d0*M*param_com(2) !Travail opérateur pour vitesse constante 
    w0 = 0.d0
    zeta = 0.d0
    call carac_oscil(w, k, M, Ma, Cs, param_com, param_geom, w0, zeta)

    write(12,*) t*param_com(2)*(RHO/M)**(1.d0/3.d0), &
                w0, &
                zeta

end subroutine ecrire
