# Bonnes pratiques git : 

~ Inspiré de JC Layoun, légèrement modifié !

## I - Initialisation a git :

### 1) Cloner
```shell
git clone https://github.com/PaulCalot/Neutralisation_Faisceau_ions_negatifs_propulsion_satellite.git
git fetch origin le_nom_de_ma_branche
git checkout le_nom_de_ma_branche
```

A l'issue, on se trouve sur sa branche, prêt à réaliser des modifications.

### 2) Faire les premières modifs !

### 3) ajouter ses modifs
```shell
git status
git add -A
git commit -m "commentaires pertinents."
git push origin le_nom_de_ma_branche
```
Plus d'info sur l'[Editing de commits]( https://github.com/k88hudson/git-flight-rules#editing-commits).
### 4) Faire une merge request:

Se fait sur la page du projet sur github directement. Attention à l'option "supprimer branche après fusion".

## II - Bonnes pratiques :

### 1) Toujours checker sur quelle branche on travaille avant de faire des modiffs :
```shell
git branch -a
```
Ou pour plus d'info :
```shell
git status
```

### 2) A chaque fin de séance de prog faire :

git checkout le_nom_de_ma_branche
git add -A
git commit -m "Avec un commentaire pertinant qui décrit simplement ce que vous avez fait."
git push origin branch_Name

### 3) A chaque merge request acceptée faire, pull la branche master avec les nouvelles modifications :
```shell
git checkout master
git pull origin master
git checkout le_nom_de_ma_branche
git merge master
git push origin le_nom_de_ma_branche
```

### 4) Il ne faut pas écrire sur les mêmes lignes que quelqu'un d'autre dans un même fichier !
  
C.f. :
https://git-scm.com/book/en/v2/Git-Branching-Basic-Branching-and-Merging

### 5) Une façon d'utiliser git :
- branche *master* : contient une version fonctionnelle du code ;
- Différentes branches *issXX* (XX = 12 par exemple)  peuvent être créées pour régler des problèmes ou ajouter des fonctionnalités. Ici, *iss* signifie *issue* (mais notations existent). Dans ce cas, on clonera d'abord la branche master en utilisant : 
```shell
git checkout -b issXX
```
- Une fois les modifications, ajouts et tests réalisés, se reporter à I.3) et I.4) pour le *push* et la *merge request*, puis supprimer cette branche en utilisant :
```shell
git branch -d nom_de_la_branche
```

## Liens utiles pour comprendre git : 
- [tuto](https://rachelcarmena.github.io/2018/12/12/how-to-teach-git.html)
- [liste exhaustives commandes git](https://github.com/k88hudson/git-flight-rules)