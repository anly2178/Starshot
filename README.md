#starshot_lib
Python library relating to the Starshot Breakthrough project.

As of 22/04/2020...

This directory should be placed in the same directory as the file using it.

To use the functions, one must define parameters for the spacecraft in this format:
	
	params = {"m_sail": <value (kg)>, "thickness": <value (m)>, "density": <value (kgm^-3)>, "reflectivity": <value>, "absorptance": <value>,  "k": <value>, "power": <value (W)>, "laser_size": <value (m)>, "wavelength": <value (m)>, "alpha": <value>}
	
Definitions:
