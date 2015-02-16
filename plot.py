'''
	Generates surface plot of input file
	Usage:
		python plot.py [input file]

	'input fie' default = 'out.dat'
'''

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from numpy import loadtxt, linspace, meshgrid
import sys

fig = plt.figure()
ax = fig.add_subplot(111)
# ax = fig.add_subplot(111, projection='3d')

z = loadtxt('xi.dat' if len(sys.argv) < 2 else sys.argv[1])

m = z.shape[0]
n = z.shape[1]

x, y = meshgrid(linspace(1,10,n), linspace(0,1,m))

# stride = max([(n-1)/64, 1])

# ax.plot_surface(x, y, z, 
# 	rstride=n, 
# 	cstride=m, 
# 	cmap=cm.coolwarm, 
# 	# linewidth=0
# )
CS = ax.contourf(x, y, z)
CB = plt.colorbar(CS, shrink=0.8, extend='both')
# ax.axis('image')

print z.max()

plt.show()

