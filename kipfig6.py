import frac_power

mass_tot = 1e-3 #kg
sail_area = 16 #m^2 for one side
beta_targ = 0.2 #fraction of c
max_temp = 300 + 273.15 #K
wavelength = 650e-9 #m
emissivity = 1
transmitter_D = [1, 1e1, 1e2, 1e3] #m
absorptivity = np.logspace(-12,0,100) #absolute absorptivity
reflectivity = 1 - absorptivity #absolute reflectivity

fig, ax1 = plt.subplots(constrained_layout=True)

i = 0
while i < len(transmitter_D):
    diameter = transmitter_D[i]
    r_tot = find_total_relative_energy(beta_targ, reflectivity)
    t = find_firing_time(mass_tot, sail_area, absorptivity, r_tot,\
    max_temp, emissivity)
    L = find_acceleration_dist(beta_targ, t)
    W = find_beam_width(L, wavelength, diameter)
    F = find_fraction_incident(sail_area, W)
    ax1.loglog(absorptivity, F)
    i += 1

ax1.set_ylim(1e-10,1.5)
ax1.set_ylabel('laser power fraction on the sail, $F_P$')

ax1.set_xlim(1e-12, 1)
ax1.set_xlabel('absorptivity of sail')

tick_times = np.array([1e-3, 1, 60, 3600, 3600 * 24, 3600 * 24 * 30])
tick_location = np.array([])
A_guesses = np.array([1e-10, 1e-7, 1e-6,1e-4,1e-3,1e-2])

surface_density = mass_tot / sail_area
e = emissivity
beta = beta_targ

i = 0
while i < len(tick_times):
    t = tick_times[i]
    eq = lambda A : (surface_density * A * (c**2) * (-(1+1-A)*(1-beta)+\
    np.sqrt(8*beta*(1-A)*(1-beta)+((1-beta)**2) * ((1+1-A)**2)))/\
    (4*(1-A)*(1-beta))) / (2 * e * stefan_boltzmann * max_temp**4) - t
    A_guess = A_guesses[i]
    solution = fsolve(eq, A_guess)
    tick_location = np.append(tick_location, solution)
    i += 1

ax2 = ax1.twiny()
ax2.set_xscale('log')
ax2.set_xlim(ax1.get_xlim())
ax2.set_xticks(tick_location)
ax2.set_xticklabels(["1 ms", "1 s", "1 min", "1 hr", "1 d", "1 mo"])
ax2.set_xlabel("continuous firing time")

plt.show()
