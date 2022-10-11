import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcolors
import newtonian_panel_code as npc
import unstruct_mesh_gen as mesh

# generate mesh
points, tri = mesh.triangulate_sphere_cone()

# run panel code
aoa_vec = np.linspace(0, np.radians(20), 20)
cd_vec, cl_vec, cp_vec = npc.newtonian_driver(points, tri, aoa_vec)

# plotting
plt.plot(np.degrees(aoa_vec), cd_vec, label='Newtonian', color='black')
#plt.scatter([7.15, 8.5, 10.5, 12.2, 17.1, 24.7], [1.506, 1.510, 1.515, 1.4992, 1.4939, 1.4816], label='LAURA Axisymm. CFD', color='black')
plt.scatter([0, 5, 10], [1.4939, 1.4806, 1.4208], label='LAURA CFD $M_{a}=17$', color='black', marker='s')
plt.scatter([0, 5, 10], [1.4816, 1.4708, 1.4303], label='LAURA CFD $M_{a}=25$', color='black', marker='^')
plt.xlabel(r'$\alpha$ $(^{\circ})$')
plt.ylabel('$C_{A}$')
plt.title(r'Stardust Hypersonic Axial Force Coefficient')
plt.ylim((1,2))
plt.grid()
plt.legend()
plt.show()

# plot the geometry
fig, ax =  plt.subplots(subplot_kw={"projection": "3d"})
#ax.quiver(centroid[0], centroid[1], centroid[2], 0.5*normal[0], 0.5*normal[1], 0.5*normal[2])
collec = ax.plot_trisurf(points[:,0].flatten(), points[:,1].flatten(), points[:,2].flatten(), triangles=tri.simplices, cmap=cm.viridis)
collec.set_array(np.array(cp_vec[0]))
fig.colorbar(collec, shrink=0.5, aspect=5)
ax.set_box_aspect((np.ptp(points[:,0]), np.ptp(points[:,1]), np.ptp(points[:,2])))
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax.set_title(r'Pressure Coefficient at $\alpha=0^{\circ}$')
plt.show()


