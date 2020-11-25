from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import time
import traceback 
import nbody

import math
import random

def time_bench(fmethod, fmethodstr):
    times = []
    vals = []
    i = 2
    total = 0
    while (total < 60):
        sim = nbody.Sim(1, fmethod, nbody.Euler)
        for j in range(10**i):
            sim.add_body(random.randrange(10**5, 10**25, 10**3), [random.randrange(-10**10, 10**10, 10**3),random.randrange(-10**10, 10**10, 10**3), random.randrange(-10**10, 10**10, 10**3)], [random.randrange(-10**4, 10**4, 10),random.randrange(-10**4, 10**4, 10),random.randrange(-10**4, 10**4, 10)], [0,0,0])
        before = time.time()
        sim.update()
        after = time.time()
        i+=1
        total=after-before
        times.append(after-before)
        vals.append(10**i)
        print("%s-body tick completed in %s seconds" % (10**i, after-before))
    plt.xlabel("Number of bodies in N-body simulation")
    plt.ylabel("Time per tick (s)")
    plt.plot(vals, times, label=fmethodstr)

def acc_bench(fmethod, tmethod, methodstr):
    errors = []
    before = time.time()
    bodies = []
    for j in range(3):
        bodies.append(nbody.Body(random.randrange(10**5, 10**25, 10**3), [random.randrange(-10**10, 10**10, 10**3),random.randrange(-10**10, 10**10, 10**3), random.randrange(-10**10, 10**10, 10**3)], [random.randrange(-10**4, 10**4, 10),random.randrange(-10**4, 10**4, 10),random.randrange(-10**4, 10**4, 10)], [0,0,0]))
    for i in range(100):
        sim1 = nbody.Sim(1, nbody.Direct, nbody.Euler)
        sim1.set_bodies(bodies)
        sim2 = nbody.Sim(i, fmethod, tmethod)
        sim2.set_bodies(bodies)
        for j in range (1000):
            sim1.update()
            sim2.update()
        errors.append(math.sqrt(pow(sim2.bodies[0].position[0] - sim1.bodies[0].position[0],2)+pow(sim2.bodies[0].position[1] - sim1.bodies[0].position[1],2)+pow(sim2.bodies[0].position[2] - sim1.bodies[0].position[2],2)))
        after = time.time()
        print("Loop %s completed in %s seconds" % (i, after-before))
    plt.plot(errors, label=methodstr)
    plt.xlabel("Size of timestep (s)")
    plt.ylabel("Approximate position error (m)")
import matplotlib.animation as animation

def sample_orbit():
    sim = nbody.Sim(200, nbody.TreePM, nbody.Leapfrog)
    for j in range(100):
        sim.add_body(random.randrange(10**5, 10**25, 10**3), [random.randrange(-10**6, 10**6, 10**2),random.randrange(-10**6, 10**6, 10**2), random.randrange(-10**6, 10**6, 10**2)], [0,0,0], [0,0,0])
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    def animate(test):
        ax.clear()
        ax.set_zlim(-10**6, 10**6)
        ax.set_ylim(-10**6, 10**6)
        ax.set_xlim(-10**6, 10**6)
        for j in range(100):
            sim.update()            

            #print (sim.bodies[0].net_force)
        for j in range(len(sim.bodies)):
            ax.scatter3D([sim.bodies[j].position[0]], [sim.bodies[j].position[1]], [sim.bodies[j].position[2]]);
        
    ani = animation.FuncAnimation(fig, animate, interval=1) 
    plt.show()
    
def max_bench(fmethod, fmethodstr):
    i=2
    while (True):
        sim = nbody.Sim(1, fmethod, nbody.Euler)
        for j in range(i):
            sim.add_body(random.randrange(10**5, 10**25, 10**3), [random.randrange(-10**10, 10**10, 10**3),random.randrange(-10**10, 10**10, 10**3), random.randrange(-10**10, 10**10, 10**3)], [random.randrange(-10**4, 10**4, 10),random.randrange(-10**4, 10**4, 10),random.randrange(-10**4, 10**4, 10)], [0,0,0])
        before = time.time()
        for j in range(1):
            sim.update()
        after = time.time()
        print(i, after-before)
        if (after-before > 0.5): break
        i+=10000
    return i
# acc_bench(nbody.Direct, nbody.Euler, "Direct")
# acc_bench(nbody.TreePM, nbody.Euler, "TreePM")
# plt.legend()
# plt.show()


#time_bench(nbody.Direct, "Direct")
#time_bench(nbody.TreePM, "TreePM")
#plt.legend();
#plt.show()
sample_orbit()
#print (max_bench(nbody.TreePM, "TreePM"))
#print (max_bench(nbody.Direct, "Direct"))

