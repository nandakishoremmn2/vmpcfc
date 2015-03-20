'''
	Generates surface plot of input file
	Usage:
		python plot.py [input file]

	'input fie' default = 'out.dat'
'''

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from numpy import loadtxt, linspace, meshgrid, cos, sin
import sys

fig = plt.figure()
ax = fig.add_subplot(111)
# ax = fig.add_subplot(111, projection='3d')

xi = loadtxt('xi.dat' if len(sys.argv) < 2 else sys.argv[1])[1:-1,:-85]

coords = open('coords.dat')
t = coords.readline()
r = coords.readline()
t = [float(tt) for tt in t.strip().split()][1:-1]
r = [float(rr) for rr in r.strip().split()][:-85]

r, t = meshgrid(r, t)

x = r*cos(t)
y = r*sin(t)

fi = r*cos(t) + xi

# stride = max([(n-1)/64, 1])

# ax.plot_surface(x, y, z, 
# 	rstride=n, 
# 	cstride=m, 
# 	cmap=cm.coolwarm, 
# 	# linewidth=0
# )
# CS = ax.plot_wireframe(x, y, z)
# CS = ax.contourf(x, y, xi, 100)
CS = ax.contour(x, y, fi, 100)
# CS = ax.contourf(r, t, z)
CB = plt.colorbar(CS, shrink=0.8, extend='both')
ax.axis('image')

# print z.max()

plt.show()

