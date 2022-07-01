from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

# Propagator class, 2 body problem + plotting

class Propagator:

    def __init__(self, t0, tf, dt, mu, R):
        self.t_init = t0
        self.t_final = tf
        self.t_step = dt
        self.mu = mu
        self.radius = R

        return

    def __ode(self, r, t):
        mu = self.mu   # Earth's gravity constant
        rmod = np.sqrt(r[1]**2 + r[0]**2)

        dr0 = r[3]
        dr1 = r[4]
        dr2 = r[5]
        dr3 = -(mu/(rmod**3)) * r[0]
        dr4 = -(mu/(rmod**3)) * r[1]
        dr5 = -(mu/(rmod**3)) * r[2]

        return [dr0, dr1, dr2, dr3, dr4, dr5]

    def prop(self, r0):
        t = np.arange(self.t_init, self.t_final, self.t_step)
        self.r = odeint(self.__ode, r0, t)

        return self.r

    def plot_orbit(self):
        r = self.r
        fig = plt.figure()
        ax = plt.axes(projection='3d')

        ax.plot3D(r[:,0],r[:,1], r[:,2])

        R_earth = self.radius

        # plot sphere
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        x = R_earth*np.cos(u)*np.sin(v)
        y = R_earth*np.sin(u)*np.sin(v)
        z = R_earth*np.cos(v)
        ax.plot_surface(x, y, z)

        # plot graphic properties
        s = 30 # scale parameter for axis
        ax.set_xlim(-s*1e6,s*1e6)
        ax.set_ylim(-s*1e6,s*1e6)
        ax.set_zlim(-s*1e6,s*1e6)
        ax.set_xticklabels([]) # disables tick labels
        ax.set_yticklabels([])
        ax.set_zticklabels([])
        ax.set_box_aspect((np.ptp(x), np.ptp(y), np.ptp(z)))

        plt.show()