# TODO list

## Premier temps
- [x] passage à un trou
- [x] diviser et simplifier le code
- [x] modifier les noms afin de respecter les conventions d'écriture (cf. *python_usage.md*)
- [ ] faire la classe **Particule** (5 types) - classe part (nom, charge, x, y vx, vy)
- [ ] Ajouter le type de rebond en param (part : oui / non)
- [ ] Collision élastique (oui/non), avec ou sans changement de charge (oui/non)
- [ ] Params initiaux : Débit de particules, nb de particules (N et p1,p2,p3,p4,p5), mode de fonctionnement (suivi de (t,x(t),y(t),nom(t))+v(tf)+nbimpact(part)+nbimpact et plot trajectoires / résultat)
- [ ] Etudier l'équation de l'équipotentiel d'entrée (une parabole suffira) 
- [ ] Fonction pour traduire la distribution sur la parabole en une distribution horizontale des particules charées (neutre = isotrope).
- [ ] Class grille à créer
- [ ] Factoriser RK pour plusieurs particules en même temps (gérer *dt* pour l'injection de beaucoup de particules) :
    - [ ] y ajouter les contacts inter particules
    - [x] changer les rebonds en téléportation sur certains bords
    - [ ] y ajouter les probas de changement de la charge en cas d'impact
    - [ ] faire une fonction collision boite noire sur les électrodes
- [ ] Plotter les données de sortie nécessaires : on sort les proportions de neutres, la distribution de alpha, V(norme), le volume cumulé sorti.
- [ ] Potentiellement étudier plusieurs trous en parallèle en changeant les distributions d'entrée

## Deuxième temps 
- [ ] Fonction boîte noire sur les électrodes à réaliser (cf. articles ou cours de Benjamin)
- [ ] Modélisation du plasma dans le moteur (au dessus) d'un point de vue particule fluide (! champ EB dans le vide)
