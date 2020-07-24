from material import Material
from save_load_mat import load_materials

a = Material(name='bob', density='123', n_list=[4,5,6], k_list=[7,8,9],
    has_equations_for_k=True, has_equations_for_n=True)

b = Material(name='woo', density='hello', n_list=[4,5,6], k_list=[7,8,9],
    has_equations_for_k=True, has_equations_for_n=True)

c = load_materials()
print(f'There are {len(c)} materials saved.')
print(c[1].get_density())
