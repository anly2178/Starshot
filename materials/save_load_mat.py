import pickle
import os

def load_materials():
    """Loads a list of saved materials"""
    materials = []
    if not os.path.exists('material_data.pkl'):
        with open('material_data.pkl', 'w'): pass
    with open('material_data.pkl', 'rb') as input:
        while True:
            try:
                materials.append(pickle.load(input))
            except EOFError:
                break
    return materials

def has_saved(material):
    """Check if material has previously been saved."""
    materials = load_materials()
    for mat in materials:
        if material.get_name() == mat.get_name():
            return True
    return False

def update_material(material, old_name):
    """Update material in pkl file with the same name."""
    old_materials = load_materials()
    new_materials = []
    for mat in old_materials:
        if mat.get_name() == old_name:
            new_materials.append(material)
        elif mat.get_name() == material.get_name():
            continue
        else:
            new_materials.append(mat)
    with open('material_data.pkl', 'wb') as f:
        for mat in new_materials:
            pickle.dump(mat, f, pickle.HIGHEST_PROTOCOL)

def delete_material_data():
    """Deletes pkl file containing materials."""
    if os.path.exists('material_data.pkl'):
        os.remove('material_data.pkl')


#To allow users to update a material e.g. if defined incorrectly,
#must reload a list of all materials, delete the incorrect material,
#save a new list of materials.
#Also, need to delete old pkl file.
