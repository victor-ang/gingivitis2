# Delaunay Mass Spring Damper  
  
This folder contains the body and plugin classes that could be used to simulate a 3D mass spring damper physics for cells using *MecaCell*  

## Prerequisites

- work on a 3D *MecaCell* simulation
- have a general body class that will inherit the **BodyDelaunayMassSpringDamper** class
- have a general plugin class or struct  that will use the **PluginDelaunayMassSpringDamper** class


## How to use it

- include the folder in your MecaCell project
- import the **BodyDelaunayMassSpringDamper.hpp** file in the file containing your simulation general body class
- make your general body class inherit from the **BodyDelaunayMassSpringDamper** class
```cpp 
template <class cell_t> class Bodies3D : public BodyDelaunayMassSpringDamper
```
- import the **PluginDelaunayMassSpringDamper.hpp** file in the file containing your general plugin class or struct
- add a **PluginDelaunayMassSpringDamper** attribute to your  general plugin class or struct, this class has to be templated by a cell class
```cpp
PluginDelaunayMassSpringDamper<cell_t> physicsPlugin;
```
- in the onRegister method, add the new attribute
```cpp
template <typename world_t>  
void onRegister(world_t* w){  
    w->registerPlugin(physicsPlugin);  
}
```
Your simulation will now implement cells that are subjected to a mass spring damper physics

## How to calibrate your simulation

### For your cells
If you want the physics to be coherent, you will need to set a few parameters.

#### Radius
If your cells' radius has to change during the simulation, you can modify it by calling the setRadius method of the body class. 
For example within the cell class :
```cpp
this->getBody().setRadius(50.);
```
The default radius at construction is set to **40µm**.
Setting the radius will automatically update the mass knowing the current density.

#### Density
If your cells' density has to change during the simulation, you can modify it by calling the setDensity method of the body class. 
For example within the cell class :
```cpp
this->getBody().setDensity(0.0013);
```
The default density at construction is set to **0.001 ng/µm³** = 1000kg/m³, the density of water.
Setting the density will automatically update the mass knowing the current radius.

#### Adhesion
The adhesion attribute of the body class represents the rate of interpenetration of a cell and is settable.
For example within the cell class :
```cpp
this->getBody().setAdhesion(0.8);
```
The default value at construction is set to **0.75**, it means that the stable position between 2 cells is equal to the sum of their radius multiplied by their adhesion.

### Physics parameters

You can modify a few parameters to adapt the accuracy of the calculations and set some cellular limitations.

#### Physics dt
You can set the dt used to compute. Having a smaller dt will make the simulation more accurate but will also increase the calculation time and vice versa.
You'll have to access it from the world in your scenario at initialization.
```cpp
w.cellPlugin.physicsPlugin.setPhysicsDt(0.001);
```
The default value is set to **0.01 s**.

#### Coherence Coefficient
After each calculation of the forces and new positions of each cells, some cells might have moved too much that their real position is not coherent anymore with their position in the Delaunay triangulation.
This coefficient is used as an indicator to know if the Delaunay has to be recomputed or if a cell has to be replace in the triangulation.
You'll have to access it from the world in your scenario at initialization.
```cpp
w.cellPlugin.physicsPlugin.setCoherenceCoeff(0.25);
```
The default value is set to **0.2**, it means that a cell is considered incoherent in the triangulation if its real position its position in the triangulation are distant from 0.2 radius away.
It prevents from useless calculation of the triangulation.

#### Maximal Speed
Maximal movement speed of the cells can be set.
You'll have to access it from the world in your scenario at initialization.
```cpp
w.cellPlugin.physicsPlugin.setMaxSpeed(20.);
```
The default value is set to **10 µm/s**, it means that a cell can't move faster than 10µm/s during a physics step. It prevents the cells to keep a coherent speed during simulation.

#### Minimal Force
The minimal force to be considered can be set.
You'll have to access it from the world in your scenario at initialization.
```cpp
w.cellPlugin.physicsPlugin.setMinForce(0.5);
```
The default value is set to **1 pN** = 1 ng.µm/s², it means that if all forces are weaker than 1 pN, it is not necessary anymore to compute the forces on each cell.
It prevents from useless calculation of the physics.

#### Neighbour Coefficient
 Each step of the simulation, the connections between cells are recomputed. Each connection is represented by a spring. A spring is created if 2 cells are connected in the triangulation and if their distance from one to another is smaller than the sum of their radius multiplied by the neighbour coefficient.
 You'll have to access it from the world in your scenario at initialization.
```cpp
w.cellPlugin.physicsPlugin.setNeighbourCoeff(1.25);
```
The default values is set to **1.1**, it means that 2 connected cells in the triangulation are considered as attached to each other if the distance between them is smaller than the sum of their 2 radius multiplied by 1.1. 

