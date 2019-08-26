#!/usr/bin/env

import numpy
from matplotlib import pyplot
import matplotlib.pyplot as plt
import time, sys

nx = 41
dx = 2.0 / ( nx - 1 )
nt = 25
dt = 0.025
c = 1

u = numpy.ones(nx)
u[ int( .5 / dx ) : int( 1 / dx + 1 ) ] = 2
print(u)

un = numpy.ones(nx)

for n in range(nt):
    un = u.copy()
    for i in range( 1 , nx ):
        u[i] = un[i] - c * dt / dx * ( un[i] - un[i - 1] )

pyplot.plot(numpy.linspace(0, 2, nx), u)
pyplot.show()



