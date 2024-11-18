from eventhandler import Event
from intersectiongraph import IntersectionGraph
from surface import Surface, Grain
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection



class Record:
    def __init__(self, time: float, event: Event):
        self.time = time
        self.event = event

class Historian: # historian keeps track of the history of the system, and writes a chronological record of the system
    
    timestamps = []
    history = []
    surfaces = []
    grains = []
    
    def __init__(self):
        self.history = []

    def add_record(self, record):
        self.history.append(record)

    def sort_history(self):
        self.history.sort(key=lambda x: x.get_time())

    def write_surface_area_history(self):
        pass

    def write_volume_history(self):
        pass

    def get_surface_area_polynomial(self, surface: Surface, interval: tuple[float, float]): # get the polynomial parameters of the surface area of a surface in a given time interval
        def surface_area(t):
            pass
        return surface_area

    def get_volume_polynomial(self, grain: Grain, interval: tuple[float, float]): # get the polynomial parameters of the volume of a grain in a given time interval
        pass

    def save(self, filename):
        pass

    def visualize(self, graph: IntersectionGraph, surfaces: list[Surface], time_interval: tuple[float, float] = None, time_step: float = 0.1):
        fig = plt.figure()
        ax = Axes3D(fig, auto_add_to_figure=False)
        fig.add_axes(ax)

        for surface in surfaces: 

            surfaceNodes_ordered_list = []
            surfaceNodes = [node for node in graph.nodes if surface in node]
            # print(len(surfaceNodes))
            
            while len(surfaceNodes)>0:
                surfaceNodes_ordered = []
                node1=surfaceNodes[0]
                surfaceNodes_ordered.append(node1)
                node_prev = node1
                surfaceNodes.remove(node1)
                complete = False

                while not complete: 
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

            intensity = 0.7+0.3*np.dot(normal, np.array([0,0,1]))
            
            alpha=1

            if surface.id in [0,1,2,3,4,5]:
                color = [(intensity*0.5647,intensity*0.9333,intensity*0.5647)]
            elif surface.id in [6,7]:
                color = [(intensity,intensity*0.7137,intensity*0.7569)]
            else:
                color = [(0,0,0)]
                alpha=0.1

            print(color)

            for i in surfaceNodes_ordered_list: 

                surface_vertices_coord = [graph.get_node_position(node,event_handler.get_time()) for node in i]

                print(surface_vertices_coord)

                # check surface id to decide color         

                surface_vertices_coord = [list( tuple(i) for i in surface_vertices_coord)]

                print(surface_vertices_coord)

                ax.add_collection3d(Poly3DCollection(surface_vertices_coord,facecolor=color, alpha=alpha, linewidths=1, edgecolors='dimgray'))


        ax.set_xlim([-1, 100])
        ax.set_ylim([-1, 100])
        ax.set_zlim([-1, 100])

        ax.set_box_aspect([1,1,1])


        plt.show()


    def plot_surface_area(self):
        for interval in (self.timestamps, np.roll(self.timestamps, -1)):
            for surface in self.surfaces:
                poly=self.get_surface_area_polynomial(surface, interval)
                x = np.linspace(interval[0], interval[1], 100)
                y = poly(x)

        
    def plot_volume(self):
        pass


