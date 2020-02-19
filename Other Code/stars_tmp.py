from vpython import *
import math
import numpy as np
# from array import array

# _____________
# Stars interacting gravitationally
# Program uses numpy arrays for high speed computations

wid = 1300      # width
heit = 800      # height

Nplanets = 100  # <- change this to have more or fewer stars

G = 6.7e-11     # universal gravitational constant

# Typical values
Msun = 2E30
Rsun = 2E9
Rtrail = 1e8
L = 4e10
vsun = 0.8*sqrt(G*Msun/Rsun)

scene = display(title="Solar System", width=wid, height=heit,
                range=0.5*L, forward=(-1,-1,-1))

xaxis = curve(pos=[(0,0,0), (L,0,0)], color=vector(0.5,0.5,0.5))
yaxis = curve(pos=[(0,0,0), (0,L,0)], color=vector(0.5,0.5,0.5))
zaxis = curve(pos=[(0,0,0), (0,0,L)], color=vector(0.5,0.5,0.5))

Stars = []
colors = [color.red, color.green, color.blue,
          color.yellow, color.cyan, color.magenta]
poslist = []
plist = []
mlist = []
rlist = []

# Set up Sun
x=0.0
y=0.0
z=0.0
r = Rsun
Stars = Stars+[sphere(pos=vector(x,y,z), radius=r, color=color.yellow)]
Stars[-1].trail = curve(pos=[Stars[-1].pos], color=color.yellow, radius=Rtrail)
mass = Msun
px = 0.0
py = 0.0
pz = 0.0
poslist.append((x,y,z))
plist.append((px,py,pz))
mlist.append(mass)
rlist.append(r)

for i in range(Nplanets):
    rad = L*uniform(0.1,0.5)              # Radius
    theta = 2.0*math.pi*random()          # Theta ()
    x = rad*math.sin(theta)
    y = rad*0.03*random()
    z = rad*math.cos(theta)
    r = Rsun*uniform(0.1,0.13)

    Stars = Stars+[sphere(pos=vector(x,y,z), radius=r, color=colors[i % 6])]
    Stars[-1].trail = curve(pos=[Stars[-1].pos], color=colors[i % 6], radius=Rtrail)
#     print(Stars)

    mass = Msun*r**3/Rsun**3
    v0 = math.sqrt(G*Msun/rad)
    px = -1.0*mass*math.cos(theta)*v0
    py = 0.0
    pz = mass*math.sin(theta)*v0
    poslist.append((x,y,z))
    plist.append((px,py,pz))
    mlist.append(mass)
    rlist.append(r)

pos = np.array(poslist)        #
p = np.array(plist)            #
m = np.array(mlist)            #
m.shape = ( Nplanets+1 ,1 )      # Numeric Python: (1 by Nstars) vs. (Nstars by 1)
radius = np.array(rlist)

vcm = sum(p) / sum(m)         # velocity of center of mass
p = p - m * vcm               # make total initial momentum equal zero

dt = 300.0
Nsteps = 0
pos = pos+(p/m)*(dt/2.)       # initial half-step
time = clock()
Nhits = 0

while 1:
    rate(100)
    # _____________
    # Compute all forces on all stars

    try:  # numpy
        r = pos-pos[:, np.newaxis]                                                  # all pairs of star-to-star vectors
        # print(r)                                                                    # ie: [20, 21, 45]

        for n in range(Nplanets+1):
            r[n,n] = 1e6                                                            # otherwise the self-forces are infinite

#         print(list(r*r))
#         z = np.log10(np.abs(r*r))
#         z = np.add.reduce(r*r)
#         z = (r*r).shape
#         print(z)

        rmag = np.sqrt(np.add.reduce(r*r, -1))                                      # star-to-star scalar distances
        # print(rmag)                                                               # ie: 14.6567

        hit = np.less_equal(rmag, radius + radius[:, np.newaxis]) - np.identity( Nplanets + 1 )

        hitlist = np.sort(np.nonzero(hit.flat)[0]).tolist()                         # 1,2 encoded as 1 * Nstars + 2

        F = G*m*m[:, np.newaxis] * r / rmag[:,:, np.newaxis]**3                     # all force pairs

    except: # old Numeric
        r = pos-pos[:, np.newaxis]                                                  # all pairs of star-to-star vectors

        for n in range(Nplanets+1):
            r[n,n] = 1e6                                                            # otherwise the self-forces are infinite

        rmag = np.sqrt(np.add.reduce(r*r,-1))                                       # star-to-star scalar distances

        hit = np.less_equal(rmag, radius + radius[:, np.newaxis]) - np.identity( Nplanets + 1)

        hitlist = sort(nonzero(hit.flat))                                           # 1,2 encoded as 1*Nstars+2

        F = G*m*m[:, np.newaxis] * r / rmag[:, :, np.newaxis]**3                    # all force pairs

    for n in range( Nplanets + 1):
        F[n,n] = 0                                                                  # no self-forces

    p = p + sum(F,1) * dt
    # print(p)

    # Having updated all momenta, now update all positions
    pos = pos + (p/m) * dt
    # print(pos)

    # Update positions of display objects; add trail
    for i in range( Nplanets + 1 ):

#         print(i)
        print(Stars[i].pos)
        print("\n")
        print(pos[i])
        print("\n")

        Stars[i].pos = pos[i]

        print(Stars[i])
        print("\n")

#         if Nsteps % 10 == 0:
#             Stars[i].trail.append(pos=pos[i])

    # If any collisions took place, merge those stars
    for ij in hitlist:
        i, j = divmod(ij, Nplanets+1) # decode star pair

        if not Stars[i].visible: continue

        if not Stars[j].visible: continue

        # m[i] is a one-element list, e.g. [6e30]

        # m[i,0] is an ordinary number, e.g. 6e30

        newpos = (pos[i]*m[i,0]+pos[j]*m[j,0])/(m[i,0]+m[j,0])

        newmass = m[i,0]+m[j,0]

        newp = p[i]+p[j]

        newradius = Rsun*((newmass/Msun)**(1./3.))

        iset, jset = i, j

        if radius[j] > radius[i]:
            iset, jset = j, i

        Stars[iset].radius = newradius
        m[iset,0] = newmass
        pos[iset] = newpos
        p[iset] = newp
        Stars[jset].trail.visible = 0
        Stars[jset].visible = 0
        p[jset] = vector(0,0,0)
        m[jset,0] = Msun*1E-30  # give it a tiny mass
        Nhits = Nhits+1
        pos[jset] = (10.*L*Nhits, 0, 0) # put it far away

    Nsteps += 1
