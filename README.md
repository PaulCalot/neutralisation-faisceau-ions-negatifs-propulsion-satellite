# Negative ions ray neutralization for satellite propulsion

This on-going project takes place in Ecole polytechnique (Palaiseau, France) engineering school third year and started in mid-September 2020 until mid-March 2021.

This project is proposed and tutored by [*Laboratoire de Physique des Plasmas*](https://www.lpp.polytechnique.fr/?lang=fr) (Plasma Physics Laboratory) from Ecole polytechnique. Names of the tutors are given [below](#tutors).

## Installation

We strongly recommend you create a specific *conda* environment to run this project to avoid any compatibility problem.

For clarity purpose, we use the environment name : **NIRN** (for *Negative Ions Ray Neutralization*). However, you are completely free to use your own.

*Linux* and *Mac* users (without *Docker*):
<br>

```shell
<<<<<<< HEAD
conda create --name NIRN -c conda-forge fenics=2018 mshr=2018
conda activate NIRN
conda install -c conda-forge matplotlib=3.3.2
conda install -c anaconda scipy=1.5.2 
conda install -c conda-forge tqdm # progress bar
conda install pandas=1.0.5 
=======
conda create --name NIRN -c conda-forge fenics=2018
conda activate NIRN
conda install -c conda-forge matplotlib=3.3.2
conda install -c anaconda scipy=1.5.2 
#conda install pandas=0.20.3
conda install -c conda-forge mshr=2018
>>>>>>> origin/Paul
```

If you want to use *jupyter notebook* or *jupyter lab*, please install : 
```shell
conda install -c conda-forge notebook
conda install -c conda-forge jupyterlab
```

To handle slide show in the jupyter notebook :
```shell
conda install -c conda-forge rise
```

To handle animation saving :
```shell
conda install -c conda-forge ffmpeg
```

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

