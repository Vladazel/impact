!
! Définition de l'équa diff de la forme y'' = f(t, y, y')
! avec f = f_equa
! 
! Équation diff est : (M+Ma) y'' + Cs y' + k w = 0
! avec y = ycom + w
!

function f_equa(t, w, wp)
    !fonctionnement vectoriel
    !w = (w, w')T et wp = (w', w'')T
    implicit none
    !Paramètres de la fonction
    real(kind=8) :: t
    real(kind=8), dimension(2), intent(in) :: w
    real(kind=8), dimension(2) :: wp
    !Paramètres de la simulation
    real(kind=8) :: M, k
    real(kind=8), dimension(2) :: param_geom, param_com
    !Variables de calcul
    real(kind=8) :: f_equa
    real(kind=8) :: ycom, y, dycomdt, dydt, d2ycomdt2
    real(kind=8) :: Ma, Cs, added_mass, slamming_coef

    !Lecture des paramètres de la simulation
    open(unit = 1, file = './params.inp', status = 'old')
    read(1,*)
    read(1,*)
    read(1,*) M
    read(1,*) k
    read(1,*) param_com(1)
    read(1,*) param_com(2)
    read(1,*) param_geom(1)
    read(1,*) param_geom(2)
    close(1)
    
    !Calcul des dérivées successives de y_com
    call y_commande(t, ycom, dycomdt, d2ycomdt2, param_com) 
    y = ycom + w(1)
    dydt = dycomdt + w(2)
    
    ! Vérification que nous sommes toujours dans le cadre
    ! du modèle de Wagner
    if (dydt > 0) then
        print*, 'dydt > 0'
        stop 
    end if

    Ma = added_mass(t, y, param_geom)
    Cs = slamming_coef(t, y, param_geom)

    !écrit des variables adimensionnelles
    write(2,*) t*sqrt(k/M), w(1)/param_com(2)*sqrt(k/M), Ma/k, Cs/k

    f_equa = 1/(M+Ma) * (Cs * dydt**2 - k*w(1)) - d2ycomdt2
    wp(2) = f_equa
    wp(1) = w(2)
end function f_equa

function added_mass(t, y, param_c)
    !Calcule le coef de slamming Cs
    implicit none
    real(kind=8), parameter :: RHO = 999
    real(kind=8), parameter :: PI = 4*atan(1.d0)
    real(kind=8), intent(in) :: t, y 
    real(kind=8), dimension(2), intent(in) :: param_c
    real(kind=8) :: added_mass

    !Ajout d'un moins pour avoir un coefficient positif
    if (param_c(1) == 0.d0) then !dièdre
        added_mass = 1.d0/8.d0 * RHO * PI**3 * y**2 / (tan(param_c(2))**2)
    else if (param_c(1) == 1.d0) then !parabole
        added_mass = -2.d0 * RHO * PI * param_c(2) * y
    else
        print*, 'Masse ajoutée non définie'
        stop
    end if
end function

function slamming_coef(t, y, param_c)
    !Calcule la masse ajoutée
    implicit none
    real(kind=8), parameter :: RHO = 999
    real(kind=8), parameter :: PI = 4*atan(1.d0)
    real(kind=8), intent(in) :: t, y
    real(kind=8), dimension(2), intent(in) :: param_c
    real(kind=8) :: slamming_coef

    if (param_c(1) == 0.d0) then !dièdre
        slamming_coef = - RHO * PI**2 * y / (4 * tan(param_c(2))**2)
    else if (param_c(1) == 1.d0) then !parabole
        slamming_coef = - RHO * PI**2 * param_c(2)
    else
        print*, 'Coefficient de slamming non défini'
        stop
    end if
end function slamming_coef
