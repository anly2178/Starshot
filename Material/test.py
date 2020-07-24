from Material import Material
from load_materials import load_materials

a = Material(name='bob', density='123', n_list=[4,5,6], k_list=[7,8,9],
    has_equations_for_k=True, has_equations_for_n=True)
