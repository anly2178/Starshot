#starshot_lib
Python library relating to the Starshot Breakthrough project.

As of 22/04/2020...

###

This directory should be placed in the same directory as the file using it.

###

To use the functions, one must define parameters for the spacecraft in this format:
	
	params = { "material": 'SiO2',
			 "m_sail": 1e-3,
			 "thickness": None,
			 "area": 10,
			 "density": 1400,
			 "reflectivity": 1,
			 "abs_coeff": 1e-3,
			 "absorptance": None,
			 "k": 1,
			 "power": 1e11,
			 "laser_size": 1e4, 
			 "wavelength": 1064e-9,
			 "alpha": 1}

========================================================================================================
Parameter							Defintion
--------------------------------------------------------------------------------------------------------
"material"  	---------------------------- 	sail material; see the list below
								Input: 	string
								List of materials and their corresponding strings
       			 						Material              String
        			 						Silicon  ------------ 'Si'
        			 						Silica    ------------ 'SiO2'
        			 						Will be updated
									
"m_sail"     	---------------------------- 	mass of the sail
								Input:	int, float
								Units:		kg
								
"thickness"	----------------------------		thickness of the sail
								Input:	int, float
								Units:		m
								
"area"		----------------------------		surface area of one side of the sail
								Input:	int, float
								Units:		m^2
								
"density"	----------------------------		volume density of sail material (not including the payload)
								Input:	int, float
								Units:		kgm^-3
								
"reflectivity"	----------------------------		fraction of radiant flux reflected off the sail
								Input:	int, float

"abs_coeff"	----------------------------		absorption coefficient of the material
								Input:	int, float
								
"absorptance"	----------------------------		fraction of radiant flux absorbed by the sail
								Input:	int, float
								
"k"		----------------------------		constant relating to the shape of the sail; see the list below
								Input:	int, float
        								List of k constants
        									Shape		Constant
									Square ---------- 1
									Circle ---------- Pi/4
								
"power"		----------------------------		laser power from the DE (directed energy) system
								Input:	int, float
								Units:	W

"laser_size"	----------------------------		width of the DE system
								Input:	int, float
								Units:	m

"wavelength"	----------------------------		wavelength of the laser, not including relativistic doppler shift
								Input:	int, float
								Units:	m
								
"alpha"		----------------------------		constant relating to the shape of the DE system; see the list below
								Input:	int, float
								List of alpha constants
									Shape		Constant
									Square ---------- 1
									Circle ---------- 1.22
========================================================================================================

