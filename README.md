# DETO: Discrete Element Topology Optimization

This is a simple code for topology optimization of systems of discrete elements written in MATLAB.

## Download
To download use:
```
$ cd <path> 
$ git clone https://github.com/Connor-OS/DETO
$ cd DETO
```
<path> should be the location you wish to store DETO 

## Folder Layout
After downloading the DETO folder will contain the following folders:

| Folder | Description                        |
|--------|------------------------------------|
| src    | source files                       |
| tests  | collection of test cases           |
| dump   | repository for dump files          |
| rocket | files for running on rocket HPC    |

## Running
### Local machine 
To run the code on your own machine it can be called from the src folder with the function:

$ top(nelx,nely,mass,Diam,rmin,kspr)

Where the input parameters are:

| Input  | Description                        |
|--------|------------------------------------|
| nelx   | Elements in the x dimension        |
| nely   | Elements in the x dimension        |
| mass   | Constraint on total mass fraction  |
| Diam   | Diameter of particles              |
| rmin   | Filter length                      |
| kspr   | Unmodified spring stiffness        |

a simple benchmark simulation to run could be:

$ top(45,15,0.6,1,1.2,100)

respectively
### Tests
A number of example simulations are included in the tests folder to showcase the capabilities of the code and some adapted boundry conditions without the user having to input material parameters. These can be run from the tests folder with the matlab commands:

$ top
a simply supported beam exaple.
sim time: approx 15 mins

$ topCant
A cantillever beam example.
sim time: approx 10 mins

$ topCentA
$ topCentB
two centrally supported beam examples designed to be compared.
In each case all input and boundry conditions are identical except sim A applies roller supports while sim B applies pinned. Entended to demonstrate the geometric non-linearity in the code as described in the related research paper.
sim time: approx 10 mins each

### Running on Newcastle University Rocket cluster
The code can also be run on the Newcastle University Rocket cluster or a simmilar HPC. To do this you should upload the contents of the rocket folder to your chosen a directory on the cluster. The code can then be submited as a batch job using the command:
```
$ sbatch topR.sh
```
To alter the dimensions and and properties of this simulation you will need to eddit the values defiened at the top of the script topR.m

## Output files

The code will produce a matlab figure displaying the position and relative mass of all particles at each step of the optimisation. It will also produce two dump files deposited in the dump folder labeled XY and Itt followed by the input parameters of that simulation. These files contain usefull information on the simulation.

XY: containes x,y coordinate and relative mass data on the system that can be opened in Ovito for analysis.
Itt: Contains information on the simulation time an performance indicators at each iteration.

## Alternate energy minimization
A second energy minimsation algorithm by stepest descent has been included in the src folder labeled SDmin.m it is possible to swap out the conventional quick min algorithm by making changes to the top.m file to call this instead. In almost all cases this will slow down the simulation but has been included here so the interested user may compare the results of different optimisation techniques.

