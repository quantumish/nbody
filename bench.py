from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import time
import nbody

import math
import random

def time_bench(fmethod, fmethodstr):
    sim = nbody.Sim(1, fmethod, nbody.Euler)
    times = []
    for i in range(100):
        before = time.time()
        for j in range(i+10):
            sim.add_body(random.randrange(10**5, 10**25, 10**3), [random.randrange(-10**10, 10**10, 10**3),random.randrange(-10**10, 10**10, 10**3), random.randrange(-10**10, 10**10, 10**3)], [random.randrange(-10**4, 10**4, 10),random.randrange(-10**4, 10**4, 10),random.randrange(-10**4, 10**4, 10)], [0,0,0])
        for j in range(1):
            sim.update()
        after = time.time()
        times.append(after-before)
        print("Loop %s completed in %s seconds" % (i, after-before))
    plt.xlabel("Number of bodies in N-body simulation")
    plt.ylabel("Time per tick (s)")
    plt.plot(times, label=fmethodstr)

def acc_bench(tmethod, tmethodstr):
    errors = []
    before = time.time()
    bodies = []
    for j in range(3):
        bodies.append(nbody.Body(random.randrange(10**5, 10**25, 10**3), [random.randrange(-10**10, 10**10, 10**3),random.randrange(-10**10, 10**10, 10**3), random.randrange(-10**10, 10**10, 10**3)], [random.randrange(-10**4, 10**4, 10),random.randrange(-10**4, 10**4, 10),random.randrange(-10**4, 10**4, 10)], [0,0,0]))
    for i in range(100):
        sim1 = nbody.Sim(1, nbody.Direct, nbody.Euler)
        sim1.set_bodies(bodies)
        sim2 = nbody.Sim(i, nbody.Direct, tmethod)
        sim2.set_bodies(bodies)
        for j in range (1000):
            sim1.update()
            sim2.update()
        errors.append(math.sqrt(pow(sim2.bodies[0].position[0] - sim1.bodies[0].position[0],2)+pow(sim2.bodies[0].position[1] - sim1.bodies[0].position[1],2)+pow(sim2.bodies[0].position[2] - sim1.bodies[0].position[2],2)))
        after = time.time()
        print("Loop %s completed in %s seconds" % (i, after-before))
    plt.plot(errors, label=tmethodstr)
    plt.xlabel("Size of timestep (s)")
    plt.ylabel("Approximate position error (m)")

def sample_orbit():
    sim = nbody.Sim(200, nbody.Tree, nbody.Leapfrog)
    sim.add_body(10**9, [0,10**5, 10], [1,0,0], [0,0,0])
    sim.add_body(10**15, [0,0,0], [0,0,0], [0,0,0])
    
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    x1,x2=[],[]
    y1,y2=[],[]
    z1,z2=[],[]
    for j in range(15000):
        if (j % 10 == 0):
            x1.append(sim.bodies[0].position[0])
            y1.append(sim.bodies[0].position[1])
            z1.append(sim.bodies[0].position[2])
            x2.append(sim.bodies[1].position[0])
            y2.append(sim.bodies[1].position[1])
            z2.append(sim.bodies[1].position[2])
        sim.update()

    #print (sim.bodies[0].net_force)
    ax.scatter3D(x1, y1, z1, label="Planet");
    ax.scatter3D(x2, y2, z2, label="Star");
    plt.legend()
    plt.show()

# acc_bench(nbody.Euler, "Euler")
# acc_bench(nbody.Leapfrog, "Leapfrog")
# plt.legend()
# plt.show()

time_bench(nbody.Direct, "Direct")
plt.show()
#sample_orbit()
