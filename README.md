# Cinématique d'un corps rigide lors d'un impact hydrodynamique

## Problème : 

Nous nous intéressons à l'impact hydrodynamique d'une forme 2D avec prise en
compte des vibrations par le biais d'un ressort dont le mouvement est imposé
par une loi de commande. Le problème est 2D avec un degré de liberté selon y.
Le modèle physique choisi est le modèle de Wagner.

```
            ____
            |  |        y_commande déplace le bloc
            |__|
             <
              >         ressort raideur k
             <
           ___>__
           \    /
            \  /        masse M
             \/
       
      ----------------  surface libre
```

## But : 

Résolution de l'équation aux dérivées partielles avec deux formes disponibles
(dièdre, parabole) et une méthodes d'intégration au choix (RK4 ou RKF45).

## Mise en équations

L'équation différentielle est la suivante :

Soit y le mouvement total : y = y_commande + w

Notations : 
- yp  = y' 
- ypp = y''
- wp  = w' 
- wpp = w''

```
(M+Ma)*ypp - Cs*yp**2 + k*w = 0
```
On résout ainsi wpp = f(t, w, wp)

## Organisation

Les fichiers de calculs ont chacun un rôle :

- equadiff.f90 définit la fonction f tel que wpp = f(t, w, wp), calcule la
  masse ajoutée ainsi que le coefficient de slamming et définit la loi de
  commande.

- resolution.f90 contient le programme de résolution qui fait appel aux
  différents schémas d'intégrations et écrit les résultats dans un fichier
  resultats.dat

- schema_int.f90 contient les schemas d'intégration utilisés: RK4 et RKF45
