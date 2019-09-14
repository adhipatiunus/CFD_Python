#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D


# In[2]:


nx = 31
ny = 31
nt = 100
nit = 100
c = 1
dx = 2.0 / (nx - 1)
dy = 2.0 / (ny - 1)
x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 2, ny)
X, Y = numpy.meshgrid(x, y)

rho = 1
niu = .1
dt = .001

u = numpy.zeros((ny, nx))
v = numpy.zeros((ny, nx))
p = numpy.zeros((ny, nx)) 
b = numpy.zeros((ny, nx))

u[0, :]  = 0
u[:, 0]  = 0
u[:, -1] = 0
u[-1, :] = 1.0    # set velocity on cavity lid equal to 1
v[0, :]  = 0
v[-1, :] = 0
v[:, 0]  = 0
v[:, -1] = 0


# In[3]:


def generate_input_function(b, u, v, dx, dy, dt):
    row, col = u.shape

    for i in range(1,row-1):
        for j in range(1, col-1):
            b[i,j] = ((u[i,j+1] - u[i,j-1]) / (2 * dx *dt) 
            + (v[i+1,j] - v[i-1,j]) / (2 * dy * dt)
            + ((u[i,j+1] - u[i,j-1])/(2 * dx))**2
            + 2 * (u[i+1,j] - u[i-1,j]) / (2 * dy) 
            * (v[i,j+1] - v[i, j-1]) / (2 * dx)
            + ((v[i+1,j] - u[i-1,j])/(2 * dy))**2)
    
    return b


# In[4]:


def generate_pressure_poisson(p, dx, dy, b, nit):
    row, col = p.shape

    for it in range(nit):
        pn = p.copy()
        
        for i in range(1,row-1):
            for j in range(1,col-1):
                p[i,j] = (((pn[i,j+1] + pn[i,j-1]) * dx**2 
                            + (pn[i+1,j] + pn[i-1,j]) * dy**2 
                            - rho * b[i,j] * (dx * dy)**2)
                            / (2*(dx**2 + dy**2)))
        
        p[:, -1] = p[:, -2] # dp/dx = 0 at x = 2
        p[0, :] = p[1, :]   # dp/dy = 0 at y = 0
        p[:, 0] = p[:, 1]   # dp/dx = 0 at x = 0
        p[-1, :] = 0        # p = 0 at y = 2

    return p


# In[5]:


def generate_cavity_flow(u, v, b, p, dx, dy, dt, nt, rho, niu):
    row, col = u.shape
    un = numpy.empty_like(u)
    vn = numpy.empty_like(v)

    for t in range(nt):
        b = generate_input_function(b, u, v, dx, dy, dt)
        p = generate_pressure_poisson(p, dx, dy, b, nit)
        
        un = u.copy()
        vn = v.copy()
        pn = p.copy()
        
        u[1:-1, 1:-1] = (un[1:-1, 1:-1]-
                         un[1:-1, 1:-1] * dt / dx *
                        (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
                         vn[1:-1, 1:-1] * dt / dy *
                        (un[1:-1, 1:-1] - un[0:-2, 1:-1]) -
                         dt / (2 * rho * dx) * (p[1:-1, 2:] - p[1:-1, 0:-2]) +
                         niu * (dt / dx**2 *
                        (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
                         dt / dy**2 *
                        (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])))

        v[1:-1,1:-1] = (vn[1:-1, 1:-1] -
                        un[1:-1, 1:-1] * dt / dx *
                       (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
                        vn[1:-1, 1:-1] * dt / dy *
                       (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) -
                        dt / (2 * rho * dy) * (p[2:, 1:-1] - p[0:-2, 1:-1]) +
                        niu * (dt / dx**2 *
                       (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +
                        dt / dy**2 *
                       (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))
        
        u[0, :]  = 0
        u[:, 0]  = 0
        u[:, -1] = 0
        u[-1, :] = 1.0    # set velocity on cavity lid equal to 1
        v[0, :]  = 0
        v[-1, :] = 0
        v[:, 0]  = 0
        v[:, -1] = 0

    return u, v, p


# In[6]:


def plot_streamline(u, v, p, X, Y):
    fig = pyplot.figure(figsize=(11, 7), dpi=100)
    pyplot.contourf(X, Y, p, alpha=0.5, cmap=cm.viridis)
    pyplot.colorbar()
    pyplot.contour(X, Y, p, cmap=cm.viridis)
    pyplot.streamplot(X, Y, u, v)
    pyplot.xlabel('X')
    pyplot.ylabel('Y')

    pyplot.show()


# In[7]:


u, v, p = generate_cavity_flow(u, v, b, p, dx, dy, dt, nt, rho, niu)


# In[8]:


plot_streamline(u, v, p, X, Y)




