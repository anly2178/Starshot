# Starshot

A Python library for the Starshot initiative.

## Installation

### Clone

* Clone this repo to your local machine using [https://github.com/anly2178/Starshot.git](https://github.com/anly2178/Starshot.git)

### Setup

The directory tree should appear like this:

```bash
root
├── test.py
├── saved_materials
└── Starshot
    ├── README.md
    ├── __pycache__
    ├── __init__.py
    ├── sail.py
    ├── multilayer_sail.py
    ├── diffractive_sail.py
    ├── motion.py
    ├── gaussbeam.py
    ├── results.py
    ├── gaussbeam.py
    ├── motion.py
    ├── materials
    │   ├── README.txt
    │   ├── __pycache__
    │   ├── material.py
    │   └── save_load_mat.py
    └── tmm
        ├── __pycache__
        ├── tmm.py
        └── make_transfer_matrix.py
```
* Starshot directory is downloaded from github.

* test.py is to be created by the user. It is the script that the user runs.
The user may name the script something else. It is important that the script
is located parallel to the Starshot directory.

* saved_materials is automatically created when the user initialises a Material object.

## Usage

**To initialise a new** ```Material```:

```python
from Starshot.materials.material import Material

new_material = Material(name=insert_name, density=insert_density, n_list=insert_n_list, k_list=insert_k_list)
```
* For more detail, see the Material section.

**To initialise a new** ```Sail```:

```python
from Starshot.sail import Sail

new_sail = Sail(name=insert_name, mass=insert_mass, area=insert_area, reflectance=insert_reflectance,
  target=insert_target, power=insert_power, wavelength=insert_wavelength)
```
* For more detail, see the Sail section.

**To initialise a new** ```MultilayerSail```:

```python
from Starshot.multilayer_sail import MultilayerSail

new_multi = MultilayerSail(name=insert_name, materials=insert_materials, mass=insert_mass,
  thickness=insert_thickness, area=insert_area, reflectance=insert_reflectance, abs_coeff=insert_abs_coeff,
  target=insert_target, max_temp=insert_max_temp, power=insert_power, wavelength=insert_wavelength)
```
* For more detail, see the Multilayer Sail section.

**To calculate mission scenario**:

```python
sail_name.calculate_mission()
```

*Note*: The ```DiffractiveSail``` class will be included when it is in a useful state.

## Sail

* The ```Sail``` class is the superclass for all subclasses of sails, such as ```MultilayerSail``` and ```DiffractiveSail```. Therefore, these subclasses inherit the ```Sail``` attributes and methods.
* Sail is flat and circular. Generalising the shape is to be completed.
* Gaussian beam is produced by circular laser array.

### Attributes

* name (str) - a unique name or code that identifies the sail. Defaults to None.
* mass (float) [kg] - mass of lightsail (excluding payload). It is assumed that payload mass equals to lightsail mass as per optimal mass condition (Kulkarni 2018). Defaults to None.
* area (float) [m^2] - surface area of lightsail on one side. Defaults to None.
* radius (float) [m] - radius of lightsail.
* s_density (float) [kg/m^2] - surface density of lightsail.
* reflectance (float) - fraction of incident power that is reflected by lighsail.
* transmittace (float) - fraction of incident power that is transmitted through lightsail.
* target (float) - target speed as fraction of speed of light. Defaults to 0.2c.
* power (float) [W] - power of laser array. Defaults to None.
* wavelength (float) [m] - laser wavelength, not Doppler-shifted. Defaults to 1.064e-6 m.
* W (float) [sqrt(g)/m] - square root of 'reflectivity-adjusted-area-density' as defined by Ilic et al. (2018).
* diameter (float) [m] - diameter of circular laser array.
* angles_coeffs (list of tuples of three floats) - angle [degrees], reflection efficiency and transmission efficiency of each order.

### Methods

```python
__init__(   self, name=None, mass=None, area=None, reflectance=None,
                  target=0.2, power=None, wavelength=1.064e-6)
```
* constructor for Sail class
      The constructor for Sail class

```python
calculate_mission()
```
* Calculates the mission scenario, including distance, speed and time.
* A folder is created with 2 txt files and 1 png file. ```trajectory.txt``` file includes distance, speed and time results. ```variables.txt``` file includes the variables of the mission. ```plots.png``` file includes speed vs distance and speed vs time graphs.

## Multilayer Sail

### Attributes

### Methods

## Material

### Attributes

### Methods

## Creators
**Andrew Ly**
* [anly2178@uni.sydney.edu.au](anly2178@uni.sydney.edu.au)
* [https://github.com/anly2178](https://github.com/anly2178)

**Justin Widjaja**
* [jwid8259@uni.sydney.edu.au](jwid8259@uni.sydney.edu.au)
* [https://github.com/jwid8259](https://github.com/jwid8259)
