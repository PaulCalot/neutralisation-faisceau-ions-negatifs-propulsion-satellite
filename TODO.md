# TODO list

## Premier temps
- [x] passage à un trou
- [ ] diviser et simplifier le code
- [ ] modifier les noms afin de respecter les conventions d'écriture (cf. *python_usage.md*)
- [ ] faire la classe **Particule**
- [ ] Etudier l'équation de l'équipotentiel d'entrée (une parabole suffira)
- [ ] Fonction pour traduire la distribution sur la parabole en une distribution horizontale des particules charées (neutre = isotrope).
- [ ] Factoriser RK pour plusieurs particules en même temps (gérer *dt* pour l'injection de beaucoup de particules) :
    - [ ] y ajouter les contacts inter particules
    - [ ] changer les rebonds en téléportation sur certains bords
    - [ ] y ajouter les probas de changement de la charge en cas d'impact
    - [ ] faire une fonction collision boite noire sur les électrodes
- [ ] Plotter les données de sortie nécessaires
- [ ] Potentiellement étudier plusieurs trous en parallèle en changeant les distributions d'entrée

## Deuxième temps 
- [ ] Fonction boîte noire sur les électrodes à réaliser (cf. articles ou cours de Benjamin)
- [ ] Modélisation du plasma dans le moteur (au dessus) d'un point de vue particule fluide (! champ EB dans le vide)