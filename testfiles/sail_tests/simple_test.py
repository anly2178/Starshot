#Import Sail class
#It is rare that Sail class will be used. It is probably only used when the main variables of the sail are known,
#e.g. reflectance, mass, area, laser power and wavelength. AND the user does not wish to calculate temperature.
from Starshot.sail import Sail

#Initialise simple sail.
simple_sail = Sail(name='simp', mass=0.001, area=10, reflectance=1,
                target=0.2, power=100e9, wavelength=1.2e-6)

#Run calculate_mission. Creates a directory named after the sail, in the current working directory.
simple_sail.calculate_mission()

#Note that the sail object itself is not saved, only the results.
#The sail can easily be reconstructed by following the variables.txt file.
