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
