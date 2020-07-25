from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import nbody

sim = nbody.Sim(10000)
sim.add_body(10**24, [0,150 * 10**9, 0], [1600,0,0], [0,0,0])
sim.add_body(10**30, [0,0,0], [0,0,0], [0,0,0])

fig = plt.figure()
ax = plt.axes(projection='3d')

x1,x2,x3=[],[],[]
y1,y2,y3=[],[],[]
z1,z2,z3=[],[],[]
for j in range(500):
   x1.append(sim.bodies[0].position[0])
   y1.append(sim.bodies[0].position[1])
   z1.append(sim.bodies[0].position[2])
   x2.append(sim.bodies[1].position[0])
   y2.append(sim.bodies[1].position[1])
   z2.append(sim.bodies[1].position[2])
   sim.update()

print (sim.bodies[0].net_force)
ax.scatter3D(x1, y1, z1, label="Planet");
ax.scatter3D(x2, y2, z2, label="Star");
plt.legend()
plt.show()
