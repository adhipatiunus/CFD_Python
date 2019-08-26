from mpl_toolkits.mplot3d import Axes3D    

import numpy
from matplotlib import pyplot, cm

## Generate computational mesh
nx = 81
ny = 81
dx = 2.0 / (nx - 1)
dy = 2.0 / (ny - 1)
c = 1
nt = 100
sigma = 0.2
dt = sigma * dx

## Generate axes
x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 2, ny)
X, Y = numpy.meshgrid(x, y)

## Generate initial velocity along x-axis
un = numpy.ones( (ny, nx) )
u = numpy.ones( (nx, ny) )
u[int(0.5 / dy) : int (1.0 / dy + 1) , int(0.5 / dx) : int (1.0 / dx + 1)] = 2

## Generate initial velocity along y-axis
vn = numpy.ones( (ny, nx) )
v = numpy.ones( (nx, ny) )
v[int(0.5 / dy) : int (1.0 / dy + 1) , int(0.5 / dx) : int (1.0 / dx + 1)] = 2

## Calculate continuity equation for given initial and boundary condition
for n in range(nt + 1):
    un = u.copy()
    vn = v.copy()

    row, col = u.shape
    for i in range(1, row):
        for j in range(1, col):
            u[i,j] = un[i,j] - (un[i,j] * dt / dy * (un[i,j] - un[i-1,j])) - (un[i,j] * dt / dx * (un[i,j] - un[i,j-1]))
            v[i,j] = vn[i,j] - (vn[i,j] * dt / dy * (vn[i,j] - vn[i-1,j])) - (vn[i,j] * dt / dx * (un[i,j] - un[i,j-1]))

            u[-1,:] = 1
            u[:,-1] = 1
            u[0,:] = 1
            u[:,0] = 1

            v[-1,:] = 1
            v[:,-1] = 1
            v[0,:] = 1
            v[:,0] = 1

fig = pyplot.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, u[:], cmap = cm.viridis)

fig = pyplot.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
surf2 = ax.plot_surface(X, Y, v[:], cmap = cm.viridis)

pyplot.show()