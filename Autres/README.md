# Negative ions ray neutralization for satellite propulsion

This on-going project takes place in Ecole polytechnique (Palaiseau, France) engineering school third year and started in mid-September 2020 until mid-March 2021.

This project is proposed and tutored by [*Laboratoire de Physique des Plasmas*](https://www.lpp.polytechnique.fr/?lang=fr) (Plasma Physics Laboratory) from Ecole polytechnique. Names of the tutors are given [below](#tutors).

## Installation

We strongly recommend you create a specific *conda* environment to run this project to avoid any compatibility problem.

For clarity purpose, we use the environment name : **NIRN** (for *Negative Ions Ray Neutralization*). However, you are completely free to use your own.


*Linux* and *Mac* users (without *Docker*):
<br>
TODO : test for *Windows* user

Be careful that this command installs a pre-build version of *fenics* (made by the devs) which will not necessarily work on your system. Indeed, this pre-build version does not contain every dependency that is required to use *fenics*. In addition, and since we plan on using *mshr*, you should make sure that you have *dolfin* installed (otherwise you will ge an error).
```shell
conda create --name NIRN -c conda-forge fenics
conda activate NIRN
conda install -c conda-forge matplotlib=3.3.2
conda install -c anaconda scipy=1.5.2 
#conda install pandas=0.20.3
conda install -c conda-forge mshr
```

If you want to use *jupyter notebook* or *jupyter lab*, please install : 
```shell
conda install -c conda-forge notebook
conda install -c conda-forge jupyterlab
```

If it does not work for some reason, you can try installing *mshr* first : 
```
conda create -n NIRN -c conda-forge mshr=2019.1.0
conda install -c conda-forge matplotlib=3.3.2
```
It will automatically download *python*, *fenics* and *numpy* amongst other as they are part of its dependencies.

**Note : *fenicsy* already contains *numpy* which is consequently installed whith *fenicsy*.**

*Sources* : 
* [*fenics* documentation](https://fenicsproject.org/documentation/)
* [*Docker* documentation](https://www.docker.com/)
* [Understanding *Docker*](http://www.science.smith.edu/dftwiki/index.php/Tutorial:_Docker_Anaconda_Python_--_1)

## Authors 

[Edouard Roger](https://www.linkedin.com/in/edouard-roger-a03536194/)
<br>
[Paul Calot](https://www.linkedin.com/in/paul-calot-43549814b/)

## Tutors <a name="tutors"></a>

[Pascal Chabert](https://www.lpp.polytechnique.fr/-Pascal-Chabert-128-?lang=fr)
<br>
[Benjamin Esteves](https://www.linkedin.com/in/benjamin-esteves-9a1234157/?originalSubdomain=fr)
