from pathlib import Path
import dill as pickle
from numpy import loadtxt, pi

def mkmatdir():
    """Make saved_materials directory if it does not exist. Return Path object for directory."""
    matdir = Path('saved_materials')
    matdir.mkdir(exist_ok=True)
    return matdir

def save_material(material):
    """Save material."""
    matdir = mkmatdir()
    with matdir.joinpath(material.get_name() + '.pkl').open(mode='wb') as f:
        pickle.dump(material, f, pickle.HIGHEST_PROTOCOL)

def del_material(name):
    """Delete material according to its name."""
    matdir = mkmatdir()
    matdir.joinpath(name + '.pkl').unlink()

def material_exists(name):
    """Check if material already exists."""
    matdir = mkmatdir()
    for mat in matdir.glob('*.pkl'):
        if mat.stem == name:
            return True
    return False

def load_material(name):
    """Load material from pkl file."""
    matdir = mkmatdir()
    try:
        with matdir.joinpath(name + '.pkl').open(mode='rb') as f:
            return pickle.load(f)
    except FileNotFoundError:
        raise ValueError(f"ValueError: '{name}' does not exist. Initialise Material object for '{name}'.")

def make_list_from_file(path_flag):
    """ Takes in the path of a CSV or txt or dat file with each entry organised as
        [wavelength],[n/k] and forms a list of tuples. If any tuples are a
        different size to any of the others, will raise a ValueError.
        If the file is not of valid extension, will also raise a ValueError.

        Perhaps future implementation for tab separated files?
    """
    filepath, flag = path_flag
    try:
        ls = loadtxt(filepath, delimiter=',') #comma-separated
    except ValueError:
        ls = loadtxt(filepath) #space-separated
    if flag == 4:
        for r, wl in enumerate(ls[:,0]):
            ls[r,0] = 2*pi/wl
    print(ls)
    ls = ls[ls[:,0].argsort()]
    if flag == 2:
        ls[:, 0] *= 1e-6
    if flag == 3:
        ls[:, 0] *= 1e-9
    return ls
