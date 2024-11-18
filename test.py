from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
from intersectiongraph import IntersectionGraph
import simulator
import numpy as np
from simulation_parameters import *


surface_list=[[[(0, 0, 0), (0, 0.1, 0.1),(0.1, 0.2, 0.1), (0.1, 0.1, 0)]],
            [[(0.1, 0.1, 0.1), (0.1, 0.2, 0.2),(0.2, 0.3, 0.2), (0.2, 0.2, 0.1)]] ]

grains, surfaces, graph, substrate = simulator.simulate_nucleation(rows=3,cols=3)
grains, event_handler = simulator.find_K_events(grains,graph, substrate)
history = simulator.simulate_vertical_growth(grains, graph, event_handler, time_limit=40, max_steps=10)

# coord=[graph.get_node_position(node,event_handler.get_time()) for node in graph.nodes]
coord=[graph.get_node_position(node,10) for node in graph.nodes]


# get node coordinates, and plot nodes in 3D  

# change surface set to node set 

# get ordered vertices, plot vertices in 3D 

# plot surfaces using ordered vertices, with time and color specified 



fig = plt.figure()
ax = Axes3D(fig, auto_add_to_figure=False)
fig.add_axes(ax)


for surface in surfaces: 

    if surface.id == -1:
        continue

    surfaceNodes_ordered_list = []
    surfaceNodes = [node for node in graph.nodes if surface in node]
    # print(len(surfaceNodes))
    
    while len(surfaceNodes)>0:

        surfaceNodes_ordered = []

        node1=surfaceNodes[0]

        surfaceNodes_ordered.append(node1)

        # neighbors=graph.neighbors(node1) 

        node_prev = node1

        node_prev_prev = None

        surfaceNodes.remove(node1)

        complete = False

        while not complete: 
            # complete=False

            for node_i in graph.neighbors(node_prev): 

                complete = True
                if surface in node_i and node_i in surfaceNodes: 
                    # print(node_i in surfaceNodes_ordered)
                    # print(node_i in surfaceNodes)

                    surfaceNodes_ordered.append(node_i)

                    surfaceNodes.remove(node_i)

                    node_prev = node_i
                    complete = False
                    break

                # complete=True
            # print(surfaceNodes_ordered)
            # print(len(surfaceNodes))
        
            # if complete:
            #     break

        # print(surfaceNodes_ordered)
        # print(len(surfaceNodes))
        # surfaceNodes_ordered.append(surfaceNodes_ordered[0])
        surfaceNodes_ordered_list.append(surfaceNodes_ordered)

    surface_vertices_coor= []

    normal = surface.v_mu[:3]

    intensity = 0.7+0.3*np.dot(normal, np.array([0,0,1]))
    
    alpha=0.5
    linewidths=1

    if surface.id in [0,1,2,3,4,5]:
        color = [(intensity*0.5647,intensity*0.9333,intensity*0.5647)]
    elif surface.id in [6,7]:
        color = [(intensity,intensity*0.7137,intensity*0.7569)]
    else:
        color = [(0,0,0)]
        alpha=0
        linewidths=0

    # print(color)

    for i in surfaceNodes_ordered_list: 

        surface_vertices_coord = [graph.get_node_position(node,event_handler.get_time()+1) for node in i]

        # print(surface_vertices_coord)

        # check surface id to decide color         

        surface_vertices_coord = [list( tuple(i) for i in surface_vertices_coord)]

        # print(surface_vertices_coord)

        ax.add_collection3d(Poly3DCollection(surface_vertices_coord,facecolor=color, alpha=alpha, linewidths=linewidths, edgecolors='dimgray'))


xy_margin=1
zlim=50
ax.set_xlim([-xy_margin, SIZE[0]+xy_margin])
ax.set_ylim([-xy_margin, SIZE[1]+xy_margin])
ax.set_zlim([0, zlim])

ax.set_box_aspect([SIZE[0],SIZE[1],zlim])


plt.show()