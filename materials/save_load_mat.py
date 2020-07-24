import pickle

def load_materials():
    """Loads a list of saved materials"""
    materials = []
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
