import pickle

def load_materials():
    materials = []
    with open('material_data.pkl', 'rb') as input:
        while True:
            try:
                materials.append(pickle.load(input))
            except EOFError:
                break
    return materials
