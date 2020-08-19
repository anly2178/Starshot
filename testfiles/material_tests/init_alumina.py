from Starshot.materials.material import Material

alumina = Material(name="Al2O3", density=3.986e3, max_temp=2072, abs_coeff=1e-3, n_list_path = ('n_alumina.txt',3), k_list_path = ('k_alumina.txt',3))

alumina.print_variables()
