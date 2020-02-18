from visual import *
from time import clock
from random import random

# Stars interacting gravitationally
# Program uses numpy arrays for high speed computations

win = 600

# change this to have more or fewer stars
Nstars = 20  

# Universal gravitational constant
G = 6.7e-11  

# Typical values
Msun = 2E30
Rsun = 2E9
L = 4e10
vsun = 0.8 * sqrt(G * Msun / Rsun)

scene = canvas(title="Stars", width=win, height=win, range=2 * L, forward=vector(-1, -1, -1))

xaxis = curve(pos=[(0, 0, 0), (L, 0, 0)], color=vector(0.5, 0.5, 0.5))
yaxis = curve(pos=[(0, 0, 0), (0, L, 0)], color=vector(0.5, 0.5, 0.5))
zaxis = curve(pos=[(0, 0, 0), (0, 0, L)], color=vector(0.5, 0.5, 0.5))

Stars = []
colors = [color.red, color.green, color.blue,
          color.yellow, color.cyan, color.magenta]
poslist = []
plist = []
mlist = []
rlist = []

for i in range(Nstars):
    x = -L + 2 * L * random()
    y = -L + 2 * L * random()
    z = -L + 2 * L * random()
    r = Rsun / 2 + Rsun * random()
    Stars = Stars + [sphere(pos=vector(x, y, z), radius=r, color=colors[i % 6], make_trail=True, interval=10)]
    mass = Msun * r**3 / Rsun**3
    px = mass * (-vsun + 2 * vsun * random())
    py = mass * (-vsun + 2 * vsun * random())
    pz = mass * (-vsun + 2 * vsun * random())
    poslist.append((x, y, z))
    plist.append((px, py, pz))
    mlist.append(mass)
    rlist.append(r)

pos = array(poslist)
p = array(plist)
m = array(mlist)
# Numeric Python: (1 by Nstars) vs. (Nstars by 1)
m.shape = vector(Nstars, 1,0)  
radius = array(rlist)

# velocity of center of mass
vcm = sum(p) / sum(m)  
# make total initial momentum equal zero
p = p - m * vcm  

dt = 1000.0
Nsteps = 0
# initial half-step
pos = pos + (p / m) * (dt / 2.)  
time = clock()
Nhits = 0

while True:
    rate(100)

    # Compute all forces on all stars
# all pairs of star-to-star vectors
    r = pos - pos[:, newaxis]  
    for n in range(Nstars):
# otherwise the self-forces are infinite
        r[n, n] = 1e6  
# star-to-star scalar distances
    rmag = sqrt(sum(square(r), -1))  
    hit = less_equal(rmag, radius + radius[:, newaxis]) - identity(Nstars)
# 1,2 encoded as 1*Nstars+2
    hitlist = sort(nonzero(hit.flat)[0]).tolist()  
# all force pairs
    F = G * m * m[:, newaxis] * r / rmag[:, :, newaxis]**3  

    for n in range(Nstars):
# no self-forces
        F[n, n] = 0  
    p = p + sum(F, 1) * dt

    # Having updated all momenta, now update all positions
    pos = pos + (p / m) * dt

    # Update positions of display objects; add trail
    for i in range(Nstars):
        Stars[i].pos = pos[i]

    # If any collisions took place, merge those stars
    for ij in hitlist:
# decode star pair
        i, j = divmod(ij, Nstars)  
        if not Stars[i].visible:
            continue
        if not Stars[j].visible:
            continue
        # m[i] is a one-element list, e.g. [6e30]
        # m[i,0] is an ordinary number, e.g. 6e30
        newpos = vector(pos[i] * m[i, 0] + pos[j] * m[j, 0]) / (m[i, 0] + m[j, 0])
        newmass = m[i, 0] + m[j, 0]
        newp = p[i] + p[j]
        newradius = Rsun * ((newmass / Msun)**(1. / 3.))
        iset, jset = i, j
        if radius[j] > radius[i]:
            iset, jset = j, i
        Stars[iset].radius = newradius
        m[iset, 0] = newmass
        pos[iset] = newpos
        p[iset] = newp
        Stars[jset].trail_object.visible = False
        Stars[jset].visible = 0
        p[jset] = vector(0, 0, 0)
# give it a tiny mass
        m[jset, 0] = Msun * 1E-30  
        Nhits = Nhits + 1
# put it far away
        pos[jset] = vector(10. * L * Nhits, 0, 0)  

    Nsteps += 1
