!
! Définition de l'équa diff de la forme y'' = f(t, y, y')
! avec f = f_equa
! 
! Équation diff est : (M+Ma) y'' + Cs y' + k w = 0
! avec y = ycom + w
!

subroutine f(t, w, wp)
    !fonctionnement vectoriel
    !w = (w, w')T et wp = (w', w'')T
    implicit none
    !Paramètres de la fonction
    real(kind=8), intent(in) :: t
    real(kind=8), dimension(2), intent(in) :: w
    real(kind=8), dimension(2), intent(out) :: wp
    !Paramètres de la simulation
    real(kind=8) :: k

    open(unit = 4, file = './param_phy.inp', status = 'old')
    read(4,*) 
    read(4,*) k !raideur 
    close(4)

    if (k == 0.d0) then
        call f_equa_no_k(t, w, wp)
    else
        call f_equa(t, w, wp)
    end if
end subroutine f

subroutine f_equa(t, w, wp)
    !fonctionnement vectoriel
    !w = (w, w')T et wp = (w', w'')T
    implicit none
    !Paramètres de la fonction
    real(kind=8), intent(in) :: t
    real(kind=8), dimension(2), intent(in) :: w
    real(kind=8), dimension(2), intent(out) :: wp
    !Paramètres de la simulation
    real(kind=8) :: M, k
    real(kind=8), dimension(2) :: param_geom
    real(kind=8), dimension(3) :: param_com
    !Variables de calcul
    real(kind=8) :: ycom, y, dycomdt, dydt, d2ycomdt2
    real(kind=8) :: Ma, Cs, added_mass, slamming_coef

    !Lecture des paramètres de la simulation
    open(unit = 4, file = './param_phy.inp', status = 'old')
    read(4,*) M !masse 
    read(4,*) k !raideur 
    read(4,*) param_com(1) !type de commande
    read(4,*) param_com(2) !vitesse initiale
    read(4,*) param_com(3) !accélération
    read(4,*) param_geom(1) !type de géométrie (parabole ou dièdre)
    read(4,*) param_geom(2) !paramètre géomètrique (rayon de courbure ou angle) 
    close(4)
    
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

    wp(1) = w(2)
    wp(2) = -1/(M+Ma) * (-Cs * dydt**2 + k*w(1)) - d2ycomdt2
end subroutine f_equa

subroutine f_equa_no_k(t, w, wp)
    !!!!!!!!!!!! 
    !fonction pour equation sans ressort
    !!!!!!!!!!!! 
    !fonctionnement vectoriel
    !w = (w, w')T et wp = (w', w'')T
    implicit none
    !Paramètres de la fonction
    real(kind=8), intent(in) :: t
    real(kind=8), dimension(2), intent(in) :: w
    real(kind=8), dimension(2), intent(out) :: wp
    !Paramètres de la simulation
    real(kind=8) :: M
    real(kind=8), dimension(2) :: param_geom
    real(kind=8), dimension(3) :: param_com
    !Variables de calcul
    real(kind=8) :: ycom, y, dycomdt, dydt, d2ycomdt2
    real(kind=8) :: Ma, Cs, added_mass, slamming_coef

    !Lecture des paramètres de la simulation
    open(unit = 4, file = './param_phy.inp', status = 'old')
    read(4,*) M !masse 
    read(4,*)  !raideur 
    read(4,*) param_com(1) !type de commande
    read(4,*) param_com(2) !vitesse initiale
    read(4,*) param_com(3) !accélération
    read(4,*) param_geom(1) !type de géométrie
    read(4,*) param_geom(2) !paramètre géométrique 
    close(4)
    
    !Calcul des dérivées successives de y_com
    call y_commande(t, ycom, dycomdt, d2ycomdt2, param_com) 
    y = ycom
    dydt = dycomdt
    
    ! Vérification que nous sommes toujours dans le cadre
    ! du modèle de Wagner
    if (dydt > 0) then
        print*, 'dydt > 0'
        stop 
    end if

    Ma = added_mass(t, w(1), param_geom)
    Cs = slamming_coef(t, w(1), param_geom)

    wp(1) = w(2)
    wp(2) = -1/(M+Ma) * (-Cs * w(2)**2) - d2ycomdt2
end subroutine f_equa_no_k

function added_mass(t, y, param_c)
    !Calcule le coef de slamming Cs
    implicit none
    real(kind=8), parameter :: RHO = 999.d0
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
    real(kind=8), parameter :: RHO = 999.d0
    real(kind=8), parameter :: PI = 4*atan(1.d0)
    real(kind=8), intent(in) :: t, y
    real(kind=8), dimension(2), intent(in) :: param_c
    real(kind=8) :: slamming_coef

    if (param_c(1) == 0.d0) then !dièdre
        slamming_coef = - RHO * PI**3 * y / (4 * tan(param_c(2))**2)
    else if (param_c(1) == 1.d0) then !parabole
        slamming_coef = 2.d0 * RHO * PI * param_c(2)
    else
        print*, 'Coefficient de slamming non défini'
        stop
    end if
end function slamming_coef

subroutine y_commande(t, ycom, dycomdt, d2ycomdt2, list_param)
    ! Définit la loi de commande
    implicit none
    !1ere valeur de list_param indique le type de loi de commande
    !0 pour vitesse constante
    real(kind=8), intent(in) :: t
    real(kind=8), dimension(3), intent(in) :: list_param
    real(kind=8), intent(out) :: ycom, dycomdt, d2ycomdt2

    if(list_param(1) == 0) then    
        !Loi commande vitesse constante
        ycom = -list_param(2) * t 
        dycomdt = -list_param(2)
        d2ycomdt2 = 0
        
    else if(list_param(1) == 1  ) then
        !Loi commande accélération constante
        ycom = -list_param(3) * t**2/2.d0 - list_param(2)
        dycomdt = -list_param(3) * t - list_param(2)
        d2ycomdt2 = -list_param(3)
    else
        print*, 'Commande non définie'
        stop
        
    end if
end subroutine y_commande 
