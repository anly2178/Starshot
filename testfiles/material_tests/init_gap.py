from Starshot.materials.material import Material

#Initialise gap. Saves gap.pkl into saved_materials directory.
gap = Material(name="gap", density=0, max_temp=float('inf'), abs_coeff=0, n_list_path = ('n_gap.txt',2), k_list_path = ('k_gap.txt',2))

#Print variables.
gap.print_variables()
