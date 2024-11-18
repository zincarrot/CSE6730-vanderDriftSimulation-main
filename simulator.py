#################################################################################
#       This file is part of the Surface Growth Simulation package.             #
#       This file contains the classes for the simulation process.              #
#################################################################################

# simulation process:
# 1. initialize substrate and grains (simulator, surface)
# 2. initialize surfaces (surface)
# 3. initialize intersection graph (intersection graph)
# 4. use linear programming to find collision events (event handler)

#### repeat until time limit is reached ####
# 5. get the earliest event (event handler)
# 6. update the intersection graph (intersection graph)
# 7. update history (historian)
# 8. update event handler (event handler)
#### repeat until time limit is reached ####

# 9. rescale the time (historian)
# 10. visualize/output history (historian)


import numpy as np
# import networkx as nx
from intersectiongraph import IntersectionGraph
from surface import Surface, Grain, Substrate
from eventhandler import Event, Event_Handler
from historian import Historian, Record
from itertools import combinations
from simulation_parameters import *


def generate_jittered_points(size, rows, cols, jitter_factor=0.5): # use jittered points to ensure some evenness in the distribution of grains

    width, height = size
    cell_width, cell_height = width / cols, height / rows

    points = []
    for i in range(rows):
        for j in range(cols):
            x = j * cell_width + cell_width * 0.5 + np.random.uniform(-jitter_factor, jitter_factor) * cell_width
            y = i * cell_height + cell_height * 0.5 + np.random.uniform(-jitter_factor, jitter_factor) * cell_height
            points.append((x, y))

    return points

def simulate_nucleation(rows=ROWS,cols=COLS, size: tuple[float, float] = SIZE, nucleation_type='simultaneous', jitter_factor=0.5):
    print("Simulating nucleation...")
    np.random.seed(SEED)

    # initialize substrate and grains
    substrate = Substrate(size=size)
    grains:list[Grain] = []
    surfaces:list[Surface] = []
    surfaces.append(substrate)

    # simultaneous nucleation
    if nucleation_type == "simultaneous":
        points = generate_jittered_points(size, rows, cols, jitter_factor)
        for (x,y) in points:
            r0 = (x, y, 0, 0)
            orientation = (np.random.uniform(0, 360), np.random.uniform(0, 180), np.random.uniform(0, 360))
            grain = Grain(r0=r0, orientation=orientation)
            grain.resolve_substrate(substrate)
            grains.append(grain)
            surfaces.extend(grain.surfaces)
    else:
        raise NotImplementedError("Only simultaneous nucleation is implemented.")
    # elif nucleation_type == "sequential":
    #     sequential_nucleation(kwargs['interval'])
    # elif nucleation_type == "random_time":
    #     random_time_nucleation(kwargs['interval'])
    # else:
    #     raise ValueError("Unknown nucleation type")

    # initialize intersection graph

    graph = IntersectionGraph(surfaces=surfaces, grains=grains)

    return grains, surfaces, graph, substrate

# def simultaneous_nucleation(n_grains: int, size: tuple[float, float], eps: float = 1e-6):
#         for i in range(n_grains):
#             r0 = (np.random.uniform(0, size[0]), np.random.uniform(0, size[1]), 0, 0)
#             orientation = (np.random.uniform(0, 360), np.random.uniform(0, 180), np.random.uniform(0, 360))
#             grains.append(Grain(r0=r0, orientation=orientation))
#             surfaces.extend(grains[i].surfaces)

# def sequential_nucleation(n_grains: int, size: tuple[float, float], interval: float = 1, eps: float = 1e-6):
#     for i in range(n_grains):
#         r0 = (np.random.uniform(0, size[0]), np.random.uniform(0, size[1]), 0, interval*i)
#         orientation = (np.random.uniform(0, 360), np.random.uniform(0, 180), np.random.uniform(0, 360))
#         grains.append(Grain(r0=r0, orientation=orientation))
#         surfaces.extend(grains[i].surfaces)

# def random_time_nucleation(n_grains: int, size: tuple[float, float], interval: float = 1, eps: float = 1e-6):
#     for i in range(n_grains):
#         r0 = (np.random.uniform(0, size[0]), np.random.uniform(0, size[1]), 0, np.random.uniform(0, interval*n_grains))
#         orientation = (np.random.uniform(0, 360), np.random.uniform(0, 180), np.random.uniform(0, 360))
#         grains.append(Grain(r0=r0, orientation=orientation))
#         surfaces.extend(grains[i].surfaces)


# objective of this function is to find all K events on the substrate. the events are stored in the event handler, but not handled yet.
def find_K_events(grains: list[Grain], graph:IntersectionGraph, substrate:Substrate):
    print("finding initial collision events")
    event_handler = Event_Handler()
    ub = np.sqrt(4/3*VN_100**2+VN_002**2)
    lb = min(VN_002, VN_100)
    ratio = ub/lb
    # calculate the upper bound of the distance between two grains that can collide
    len_ub = np.sqrt((SIZE[0]/COLS)**2+(SIZE[1]/ROWS)**2)*ratio
    for grain1, grain2 in combinations(grains, 2):
        if np.linalg.norm(grain1.r0-grain2.r0) < len_ub:
            event_handler.add_K_event(grain1, grain2, graph, substrate)        

    return grains, event_handler

# TODO: implement vertical growth
# the objective of this function is to run the simulation starting from the stored K events
def simulate_lateral_growth(grains: list[Grain], graph: IntersectionGraph, event_handler:Event_Handler, time_limit=TIME_LIMIT, max_steps=MAX_STEPS, eps: float = EPS):
    print("simulating crystal growth")
    time = 0
    history = Historian()
    # history.add_record(grains, graph, time)
    steps=0
    while steps < max_steps:
        if event_handler.handle_next(graph,lateral_growth=True):
            time = event_handler.get_time()
            print(time)
            # history.add_record(grains, graph, time)
        else:
            time = event_handler.get_time()
            print(time)
            # history.add_record(grains, graph, time)
            break
        steps+=1

        
    # history.sort_history()

    return history, grains, graph, event_handler

def simulate_vertical_growth(grains: list[Grain], graph: IntersectionGraph, event_handler:Event_Handler, time_limit=TIME_LIMIT, max_steps=MAX_STEPS, eps: float = EPS):
    print("simulating crystal growth")
    time = 0
    history = Historian()
    # history.add_record(grains, graph, time)
    steps=0
    while time < time_limit and steps < max_steps:
        if event_handler.handle_next(graph):
            time = event_handler.get_time()
            print(time)
            # history.add_record(grains, graph, time)
        else:
            time = event_handler.get_time()
            print(time)
            # history.add_record(grains, graph, time)
            break
        steps+=1

        
    # history.sort_history()

    return history


def simulate_growth(n_grains: int, size: tuple[float, float], nucleation_type='simultaneous'):
    grains, surfaces, graph = simulate_nucleation(n_grains, size, nucleation_type)
    grains, surfaces, graph, event_handler = find_K_events(grains, surfaces, graph)
    history = simulate_vertical_growth(grains, surfaces, graph, event_handler)
    
    # save or visualize the history 
    history.save(filename)
    history.visualize()
    
        



if __name__ == "__main__":
    simulate_growth(N_GRAINS, SIZE, nucleation_type='simultaneous')