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
sail.calculate_mission()
```

*Note*: The ```DiffractiveSail``` class will be included when it is in a useful state.

##Creators
**Andrew Ly**
* [anly2178@uni.sydney.edu.au](anly2178@uni.sydney.edu.au)
* [https://github.com/anly2178](https://github.com/anly2178)

**Justin Widjaja**
* [jwid8259@uni.sydney.edu.au](jwid8259@uni.sydney.edu.au)
* [https://github.com/jwid8259](https://github.com/jwid8259)
