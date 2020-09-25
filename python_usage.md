# Bonnes pratiques Python
*inspiré de [ce cours](https://openclassrooms.com/fr/courses/235344-apprenez-a-programmer-en-python/235263-de-bonnes-pratiques)*

par *Paul Calot*


## Conventions de nommage
On utilise la norme **PEP 8** disponible [ici](https://www.python.org/dev/peps/pep-0008/).

```python
MaClass 	
MyError 		# Nom d'exception (une exception est une classe)
nom_de_fonction
NOM_DE_MA_CONSTANTE
```

Citation :
>De plus, les modules et packages doivent avoir des noms courts, constitués de lettres minuscules. Les noms de modules peuvent contenir des signes  **_** (souligné). Bien que les noms de packages puissent également en contenir, la PEP 8 nous le déconseille.

## *docstring* (chaîne de documentation)

* Placée juste après la définition d'un module, d'une classe, fonction ou méthode.
* Correspond à l'argument *__doc__* de l'objet.
* un package peut être documenté par une *docstring* placé dans le fichier *__init__.py*.

**Exemple**:
```python
def sinc(x):
    """ Returns the sinus cardinal of x

    Args:
        x (float): the real number from which the sinus cardinal is computed.

    Returns:
        float: the computed result
    """
    return sin(x)/x if x != 0 else 1
```

## Bonnes pratiques :
* ne pas dépasser un 79 (voire 72) caractères par ligne avec césure de ligne après opérateur (s'il y a lieu) en utilisant :
```python
un_long_calcul = variable + \
		taux * 100
```