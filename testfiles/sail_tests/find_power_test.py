#Using the same sail as in multi_test (S3), instead of providing 100GW power as input
# we find the highest power the sail can be subject to, based on the max temps
# given for each material and the Starchip. Uses Newton's method.

#All the user has to do is not provide a power into the constructor.

#Import MultilayerSail class
from Starshot.multilayer_sail import MultilayerSail

#Initialise sail without power.
multi_test = MultilayerSail(name='S3_nopower', materials=['SiO2','gap','SiO2'], mass=0.001,
    thickness=[197e-9,399e-9,197e-9], area=None, target=0.2, max_Starchip_temp=1000,
    power=None, wavelength=1.2e-6)
#Run calculate_mission()
#It will automatically calculate power, it may take a few minutes.
multi_test.calculate_mission()
