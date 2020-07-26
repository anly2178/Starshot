from material import Material
from save_load_mat import load_materials

a = Material(name='bob', density='123', n_list=[4,5,6], k_list=[7,8,9],
    has_equations_for_k=True, has_equations_for_n=True)

b = Material(name='woo', density='456', n_list=[4,5,6], k_list=[7,8,9],
    has_equations_for_k=True, has_equations_for_n=True)

for mat in load_materials():
    print(f'The density of {mat.get_name()} is {mat.get_density()}')

b.set_density('soo')

for mat in load_materials():
    print(mat.get_name())

a.set_name('12345')

for mat in load_materials():
    print(mat.get_name())

f = Material(name='loo', density='hello', n_list=[4,5,6], k_list=[7,8,9],
    has_equations_for_k=True, has_equations_for_n=True)

for mat in load_materials():
    print(mat.get_name())
