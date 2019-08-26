import numpy
import sympy
from sympy import *
from matplotlib import pyplot
import matplotlib.pyplot as plt

nx = 101
nt = 100
niu = 0.07
dx = 2 * sympy.pi / (nx - 1)
dt = dx * niu

x, nu, t = sympy.symbols('x nu t')

phi = sympy.exp(-(x - 4 * t)**2 / (4 * nu * (t + 1))) + sympy.exp(-(x - 4 * t - 2 * sympy.pi)**2 / (4 * nu * (t + 1)))
dphi = phi.diff(x)

u = - 2 * nu * dphi / phi + 4

ufunc = lambdify((x, nu, t), u)

x = numpy.linspace(0, 2 * numpy.pi, nx)

un = numpy.empty(nx)

t0 = 0

un = numpy.asarray([ufunc(x0, niu, t0) for x0 in x])

for t in range(nt):
    u = un.copy()
    for i in range (1, nx - 1):
        u[i] = un[i] - un[i] * dt * (un[i] - un[i-1]) / dx + niu * dt * (un[i+1] + un[i-1] - 2 * un[i]) / dx**2

u[0] =  un[0] - un[0] * dt * (un[0] - un[-1]) / dx + niu * dt * (un[1] + un[-1] - 2 * un[0]) / dx**2
u[-1] = u[0]

plt.figure(figsize=(11, 7), dpi=100)
plt.plot(x,u, marker='o', lw=2, label='Computational')
plt.xlim([0, 2 * numpy.pi])
plt.ylim([0, 10])
plt.legend()
plt.show()
