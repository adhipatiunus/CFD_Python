#!/usr/bin/env

import numpy
from matplotlib import pyplot
import matplotlib.pyplot as plt
import time, sys

nx = 41
dx = 2.0 / ( nx - 1 )
nt = 25
nu = 0.3
sigma = .2
dt = sigma * dx ** 2 / nu
c = 1

u = numpy.ones(nx)
u[ int( .5 / dx ) : int( 1 / dx + 1 ) ] = 2
print(u)

un = numpy.ones(nx)

for n in range(nt):
    un = u.copy()
    for i in range( 1 , nx - 1):
        u[i] = un[i] + nu * dt / dx**2 * (un[i+1] - 2 * un[i] + un[i-1])

pyplot.plot(u)
pyplot.show()

