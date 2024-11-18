# -*- coding: utf-8 -*-


from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
from intersectiongraph import IntersectionGraph
import simulator


surface_list=[[[(0, 0, 0), (0, 0.1, 0.1),(0.1, 0.2, 0.1), (0.1, 0.1, 0)]],
            [[(0.1, 0.1, 0.1), (0.1, 0.2, 0.2),(0.2, 0.3, 0.2), (0.2, 0.2, 0.1)]] ]



grains, surfaces, graph, substrate = simulator.simulate_nucleation(rows=2,cols=2)
# grains, event_handler = simulator.find_K_events(grains,graph, substrate)
# history = simulator.simulate_vertical_growth(grains, graph, event_handler, time_limit=1000, max_steps=10)

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

    surfaceNodes_ordered_list = []

    surfaceNodes = [node for node in graph.nodes if surface in node]

    print(len(surfaceNodes))
    
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
                    # print("node added:",node_i)
                    # print("previous:",node_prev_prev)
                    # try:
                    surfaceNodes.remove(node_i)
                    # except:
                    #     print(surfaceNodes)

                    # neighbors=graph.neighbors(node_i)
                    
                    # node_prev_prev=node_prev
                    node_prev = node_i
                    complete = False
                    break

                # complete=True
            # print(surfaceNodes_ordered)
            # print(len(surfaceNodes))
        
            # if complete:
            #     break

        print(surfaceNodes_ordered)
        print(len(surfaceNodes))
        surfaceNodes_ordered_list.append(surfaceNodes_ordered)

    surface_vertices_coor= []

    normal = surface.v_mu[:3]

    for i in surfaceNodes_ordered_list: 

        surface_vertices_coord = [graph.get_node_position(node,event_handler.get_time()) for node in i]

        print(surface_vertices_coord)

        # check surface id to decide color 

        color = 'b'

        

        surface_vertices_coord = [list( tuple(i) for i in surface_vertices_coord)]

        print(surface_vertices_coord)



        ax.add_collection3d(Poly3DCollection(surface_vertices_coord,facecolor=color))




ax.set_xlim([-1, 100])
ax.set_ylim([-1, 100])
ax.set_zlim([-1, 100])

ax.set_box_aspect([1,1,1])


plt.show()



