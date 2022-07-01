import propagator

# Tests the propagator class

orbit = propagator.Propagator(0, 40000, 0.1, 3.986004e14, 6370e3)
orbit.prop([7e6, 0, 0, 0, 10000, 3000])
orbit.plot_orbit()
