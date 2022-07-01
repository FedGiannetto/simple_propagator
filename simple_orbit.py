from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt
import numerical_tools as nt

CD = 2.2

def f(r,t):
    mu = 3.986004e14    # Earth's gravity constant
    RE = 6370
    rmod = np.sqrt(r[2]**2 + r[1]**2 + r[0]**2)
    # J2 terms
    J2 = 1.08262668e-3
    J2coeff = (3*J2*mu*RE**2)/(2*rmod**5)

    # Central body acceleration plus J2 perturbation (Cowell's method)

    dr0 = r[3]
    dr1 = r[4]
    dr2 = r[5]
    dr3 = -(mu/(rmod**3)) * r[0] 
    dr4 = -(mu/(rmod**3)) * r[1] 
    dr5 = -(mu/(rmod**3)) * r[2]

    return [dr0, dr1, dr2, dr3, dr4, dr5]


def SMA(r):
    rmod = np.sqrt(r[:,2]**2 + r[:,1]**2 + r[:,0]**2)
    vmod = np.sqrt(r[:,5]**2 + r[:,4]**2 + r[:,3]**2)
    mu = 3.986004e14
    a = 1 / (2/rmod - vmod**2/mu)
    return a

def periapsis(r,a):
    rmod = np.sqrt(r[:,2]**2 + r[:,1]**2 + r[:,0]**2)
    vmod = np.sqrt(r[:,5]**2 + r[:,4]**2 + r[:,3]**2)


# Starting conditions
r0 = [7e6, 0, 0, 0, 0, 8000]  # x, y, z, vx, vy, vz

# Time
t = np.arange(0, 48000, 10)

# Solve ODE
r = odeint(f, r0, t)

# --- plotting ---

fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot3D(r[:,0],r[:,1], r[:,2])

R_earth = 6370e3

# plot sphere
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x = R_earth*np.cos(u)*np.sin(v)
y = R_earth*np.sin(u)*np.sin(v)
z = R_earth*np.cos(v)
ax.plot_surface(x, y, z, alpha=0.4)

# plot graphic properties
s = 30 # scale parameter for axis
ax.set_xlim(-s*1e6,s*1e6)
ax.set_ylim(-s*1e6,s*1e6)
ax.set_zlim(-s*1e6,s*1e6)
#ax.set_xticklabels([]) # disables tick labels
#ax.set_yticklabels([])
#ax.set_zticklabels([])
ax.set_box_aspect((np.ptp(x), np.ptp(y), np.ptp(z)))


a = SMA(r)

#fig2 = plt.figure()
#ax2 = plt.subplots(2)

#ax2[0].plot(t,a)
#ax2[1].plot(t,)

print(r)

plt.show()