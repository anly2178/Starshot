#Before running this script, ensure that 'SiO2', 'gap' and 'GeO2' have been initialised; see material_tests.

#Import MultilayerSail class
from Starshot.multilayer_sail import MultilayerSail

#First, let's try a single layer sail.
#Initialise the single layer sail - S1 as defined in Ilic et al. (2018).
#Note area is not defined as it will be automatically calculated. It does not even have to be
#written in the constructor.
single_test = MultilayerSail(name='S1', materials=['SiO2'], mass=0.001,
    thickness=[206e-9], area=None, target=0.2, max_Starchip_temp=1000, power=1e11,
    wavelength=1.2e-6)
#Note that initialising a multilayer sail will print out its variables in terminal.

#Run calculate_mission() method. Description of output is given in README.md.
#Importantly, a directory named 'S1' will be produced in current working directory.
single_test.calculate_mission()

#Initialise multilayer sail - three layer sail as defined in Ilic et al. (2018).
multi_test = MultilayerSail(name='S3', materials=['SiO2','gap','SiO2'], mass=0.001,
    thickness=[197e-9,399e-9,197e-9], area=None, target=0.2, max_Starchip_temp=1000,
    power=1e11, wavelength=1.2e-6)

#Run calculate_mission() method.
multi_test.calculate_mission()

#Initialise the same multilayer sail, but this time instead of mass given, area is given.
multi_test2 = MultilayerSail(name='S3', materials=['SiO2','gap','SiO2'], mass=None,
    thickness=[197e-9,399e-9,197e-9], area=1.1557700664798942, target=0.2, max_Starchip_temp=1000,
    power=1e11, wavelength=1.2e-6)

#Run calculate_mission() again. Note that another directory is produced.
#However, since the names of the sails are identical, the new directory has (1).
#If more directories are produced with the same name, then (2), (3), (4), ... will be appended.
multi_test2.calculate_mission()
