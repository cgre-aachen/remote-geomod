# remote-geomod

![Python 3](https://img.shields.io/badge/Python-3-blue.svg)
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)

> Coupling remote geological mapping in Google Earth with direct 3D geological modeling 
in GemPy.

## Introduction

## Example

The best way to get started is actually to go through the example in the book chapter - this is currently in review and more information will (hopefully) be provided here, soon!

## Installation

Next to remote-geomod, you will require two additional major software components: (1) 
[Google Earth Pro](https://www.google.com/earth/desktop/) for remote mapping and (2) a Python 
distribution with specific software libraries. The software libraries remote-geomod relies on can be difficult to
install on certain systems. Because of this we provide you with two possible installation paths that will take care
of installing the dependecies:

 * Using the open-source **Anaconda** Python distribution which allows convenient, automated installation of the 
 dependencies. 
 * A local installation inside of a **Docker** environment, for wich we can directly ship all dependencies insode
 of a Docker image. We recommend this installation path for users familiar with using Dokcer of if the above one fails.

### Installation using Anaconda

As a first step you need to download the Anaconda Python Distribution on your system. Make sure to select the correct 
Python 3.6 version for your operating system and install it (a comprehensive user guide with installation 
instructions can be found on the [Conda Documentation website](https://conda.io/docs/user-guide/install/)).
 
The next step is to either download and unpack the remote-geomod repository, or to directly clone it using your 
command-line tool of choice:
 
    git clone https://github.com/cgre-aachen/remote-geomod.git

Once you have a local copy of the repository on your computer, you have to create a new Conda environment. You can 
either do this using your command-line tool (1) or the Anaconda Navigator (2): 

1. Open your command-line tool in the downloaded or cloned remote-geomod folder and run 
``conda env create -n rgeomod -f environment.yml``. Afterwards you can activate the environment with the command 
``activate rgeomod`` (Windows) or ``source activate rgeomod`` (macOS, Linux) or in your Anaconda Navigator.
2. For the latter option start the Anaconda Navigator, select ``Environments`` in the navigation bar on the left-hand 
side, then click ``Import`` at the bottom of the window and browse for the ``environment.yml`` file located in the 
remote-geomod folder and give the environment the same name. This will create a separate Conda environment and install 
all necessary dependencies automatically. Once the installation finished, make sure you have selected the newly 
created environment.

### Dependencies

remote-geomod uses Python 3 and has several dependencies for numerical operations, geographical data operations, 
geological modeling and visualization. All dependencies can be found in the  `environment.yml` file:

````
dependencies:
  - python=3.6
  - numpy
  - conda-forge::gdal
  - clinicalgraphics::vtk
  - theano=1.0.1
  - jupyter
  - nb_conda
  - matplotlib
  - pandas
  - seaborn
  - scikit-image
  - tqdm
  - pip:
    - mplstereonet
    - gempy
````

We also provide precompiled Docker images hosted on Docker Hub with all necessary dependencies to get remote-geomod up
and running. 

### Installation using Docker

Docker is an operating-system-level-visualization software,
meaning that we can package a tiny operating system with pre-installed
software into a Docker image. This Docker image can then be shared
with and run by others, enabling them to use intricate dependencies
with just a few commands. For this to work the user needs to have a
working [Docker](https://www.docker.com/) installation.

#### (a) Pull Docker image from DockerHub

The easiest way to get remote-geomod running is by running the pre-compiled Docker image (containing everything you
need) directly from the cloud service Docker Hub to get a locally running Docker container. Make sure to set your 
Docker daemon to Linux containers in Docker's context menu.

    docker run -it -p 8888:8888 cgreaachen/rgeomod_complete
    
This will automatically pull the Docker image from Docker Hub and run it, opening a command line shell inside of the
running Docker container. There you have access to the file system inside of the container. Note that this pre-compiled
Docker image already contains the rgeomod repository. As rgeomod is undergoing active development, we can not 
guarantee that this Docker image always contains the latest release version.

If you just want the dependencies:

    docker run -it -p 8888:8888 cgreaachen/rgeomod_dependencies

Alternatively, you can also build the Docker image yourself from the Dockerfile provided with the repository by running
``docker build . -t <image-tag>`` in its root directory. Once the Docker image is built you can look up it's image-id 
using ``docker images``. Then, instead of using the DockerHub repository name, you can run the Docker image by using
its id: ``docker run -it -p 888:8888 <image-id>``.

## Getting started

A good way to get started is to follow the steps in the publication (see below) and to run the jupyter notebooks, supplied with this package.

## References

Wellmann, F., A. Schaaf and C. von Hagke: "From GoogleEarth to 3-D Geology Problem 2: Seeing below the surface of the Digital Earth". Submitted as book chapter for "Structural Geology and Tectonics: Problems and Solutions".

## Contact

This library is developed by LuFG Computational Geoscience and Reservoir Engineering (CGRE) at RWTH Aachen University, Germany.
