from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import nbody

sim = nbody.Sim(1)
sim.add_body(10**10, [10*9,0,0], [0,0,0], [0,0,0])
sim.add_body(10**12, [0,0,0], [0,0,0], [0,0,0])

fig = plt.figure()
ax = plt.axes(projection='3d')

x1,x2,x3=[],[],[]
y1,y2,y3=[],[],[]
z1,z2,z3=[],[],[]
for j in range(10000):
   x1.append(scene.matter[0].position[0])
   y1.append(scene.matter[0].position[1])
   z1.append(scene.matter[0].position[2])
   x2.append(scene.stars[0].position[0])
   y2.append(scene.stars[0].position[1])
   z2.append(scene.stars[0].position[2])

ax.scatter3D(fixedx1, fixedy1, fixedz1, label="Planet");
ax.scatter3D(fixedx2, fixedy2, fixedz2, label="Star", s=s);
plt.legend()
plt.show()
