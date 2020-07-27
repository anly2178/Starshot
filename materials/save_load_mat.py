from pathlib import Path
import pickle
import os

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



#Create saved_materials folder if it does not exist.
#Create files for each material according to name in saved_materials
#Navigate by stem name.
