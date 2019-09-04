# Diffusion

This folder contains the body and plugin classes that could be used to simulate diffusion of molecules between cells using *MecaCell*.

It was tested on Fedora 29.

## Prerequisites

- CGAL 4.13
- work on a 3D *MecaCell* simulation
- have a general body class that will inherit the **BodyDiffusion** class
- have a general plugin class or struct that will use the **PluginDiffusion** class

## How to use it

- include the folder in your *MecaCell* project
- import the **BodyDiffusion.hpp** file in the file containing your simulation general body class
- make your general body class inherit from the **BodyDiffusion** class
```cpp 
template <class cell_t> class Bodies3D : public BodyDiffusion 
```
- in the onRegister method, add the new attribute 
```cpp 
template <typename world_t>  
void onRegister(world_t* w){      
w->registerPlugin(diffusionPlugin);  } 
``` 
Your simulation will now be able to add and diffuse molecules in the environment and make your cells interact with them.

## How to calibrate your simulation

### Instantiating the plugin

In the general plugin class, you'll have to instantiate this plugin with the other ones.
The constructor needs 2 arguments :
- dx : represents the spatial precision of your diffusion
- accuracy : represents the minimal changes in quantities between 2 diffusion steps to determine if this state is stable or not
```cpp 
PluginDiffusion diffusionPlugin = PluginDiffusion(5., 0.1);
``` 
This calibration means that the unit distance will be set to 5µm and the stability will be considered when the total quantity difference for each molecule between 2 steps of diffusion will be lesser than 0.1 mmHg.

### Adding new molecules in the simulation

You can add molecules and their characteristics in your simulation. You will have to do this through the world in your scenario at initialization.
```cpp
Molecule oxygen(2000.0, 37.50319, 1.428, 0.1);
Molecule carbon(4000.0, 4.2, 1.977, 0.01);

w.cellPlugin.diffusionPlugin.addMolecule(oxygen);
w.cellPlugin.diffusionPlugin.addMolecule(carbon);
```
The first thing to do is to instantiate your molecules.
The constructor needs 4 parameters :
- diffusion constant in µm²/s
- default quantity in the environment in mmHg
- density in kg/m³
- default absorption of the environment in µm^3/ng/s

Setting the default quantity with a negative value will make the external environment act as a sponge for this molecule.
Setting the default consumption with a negative value will make the molecule appear from nowhere

The indexes of each molecule will be in the same order as the order you added them.
In order to make it easier to access a molecule from its index, you can create an enum like this one :
```cpp
enum MOLECULE{
    OXYGEN = 0,
    CARBON = 1
};
```
It can prevent you from accessing the wrong molecule.

### Setting consumptions and accessing quantities from your cells

You can access and modify the quantity of every molecule from every cell thanks to the **BodyDiffusion** class.
```cpp
double oxygen = this->getBody().getQuantities()[MOLECULE::OXYGEN];
double carbon = this->getBody().getQuantities()[MOLECULE::CARBON];
```

You can do the same with the consumptions.
```cpp
this->getBody().setConsumption(MOLECULE::OXYGEN, 0.5634224);
this->getBody().setConsumption(MOLECULE::CARBON, -0.5);
```
A positive consumption will make the cell absorb the molecule while a negative one will make it produce the molecule.

## Other features

### Spheroids simulations

If your cells tend to form a spheroid, you can acces some interesting information from the world.

#### Spheroid's radius

You can compute and get the spheroid radius in µm like this :
```cpp
 w->cellPlugin.diffusionPlugin.getSpheroidRadius(w);
```
where w is of type World *.

#### Spheroid's center

You can also compute and get the spheroid center as a MecaCell::Vec this way :
```cpp
w->cellPlugin.diffusionPlugin.getCentroid();
```
where w is of type World *.