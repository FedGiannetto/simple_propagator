from prop_perturbation_class import prop_perturbation
import numpy as np
import numerical_tools as nt
import matplotlib.pyplot as plt

# Uses the prop_perturbation class to propagate an orbit accounting for J2 effect
# Note that Python is not the best at plotting surfaces, so the figure is not going to be super clear

prop = prop_perturbation(0, 10*80e3, 100, 3.986004e14, 6370e3)
x = prop.prop([7e6, 0, 0, 0, 7000, 9000])
prop.plot_orbit()

c = 0
inc = []

for i in range(1,len(x)):
    if x[i,1] > 0 and x[i-1,1] < 0:
        inc.append(np.arcsin(x[i,2]/nt.norm(x[i,1:3])))
        c += 1

print(inc)

fig, ax = plt.subplots()
ax = plt.plot(inc)
plt.show()