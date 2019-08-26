from mpl_toolkits.mplot3d import Axes3D    

import numpy
from matplotlib import pyplot, cm

nx = 81
ny = 81
dx = 2.0 / (nx - 1)
dy = 2.0 / (ny - 1)
c = 1
nt = 120
sigma = 0.2
dt = sigma * dx

x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 2, ny)

un = numpy.ones( (ny, nx) )
u = numpy.ones( (nx, ny) )
u[int(0.5 / dy) : int (1.0 / dy + 1) , int(0.5 / dx) : int (1.0 / dx + 1)] = 2

X, Y = numpy.meshgrid(x, y)

fig = pyplot.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, u[:], cmap = cm.viridis)


for n in range(nt + 1):
    un = u.copy()
    row, col = u.shape
    for i in range(1, row):
        for j in range(1, col):
            u[i,j] = un[i,j] - (c * dt / dy * (un[i,j] - un[i-1,j])) - (c * dt / dx * (un[i,j] - un[i,j-1]))

fig = pyplot.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
surf2 = ax.plot_surface(X, Y, u[:], cmap = cm.viridis)

pyplot.show()