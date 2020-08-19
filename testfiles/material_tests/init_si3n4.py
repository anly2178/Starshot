from Starshot.materials.material import Material

si3n4 = Material(name="Si3N4", density=3.17e3, max_temp=1900, abs_coeff=1e-3, n_list_path = ('n_si3n4.txt',2), k_list_path = ('k_si3n4.txt',2))

si3n4.print_variables()
