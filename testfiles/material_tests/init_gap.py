from Starshot.materials.material import Material

#Initialise gap. Saves gap.pkl into saved_materials directory.
gap = Material(name="gap", density=0, max_temp=1000, abs_coeff=0, n_list_path = 'n_gap.txt', k_list_path = 'k_gap.txt')

#Print variables.
gap.print_variables()
