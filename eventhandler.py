from queue import PriorityQueue
from surface import Surface, Grain
from intersectiongraph import IntersectionGraph, IntersectionEdge, IntersectionNode
from scipy.optimize import linprog
import numpy as np
from simulation_parameters import *

class Event:
    # event is a turning point. The turning point can either happen in future or in past.
    # events are found using linear programming in the lateral growth stage
    # an event can be found for each edge in the vertical growth stage
    def __init__(self, surfaces: list[Surface], nodes_in: set[IntersectionNode]):
        try:
            assert len(surfaces) == 4
        except AssertionError:
            raise AssertionError(nodes_in)

        assert len(nodes_in) <= 4 # 0: X-event 1: K-event, 2: H-event, 3: A-event, 4: Y-event
        self.surfaces = surfaces
        self.nodes = set([frozenset([surfaces[i], surfaces[j], surfaces[k]]) for i, j, k in [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]])
        self.nodes_in = nodes_in
        self.nodes_out = self.nodes - nodes_in
        self.trn_pt = Surface.intersect(surfaces)

    def __repr__(self) -> str:
        return f"Event: {self.trn_pt}, Surfaces: {self.surfaces}"
        
    def get_time(self):
        return self.trn_pt[3]
    
    def get_z(self):
        return self.trn_pt[2]
    
    def get_x(self):
        return self.trn_pt[0]
    
    def get_y(self):
        return self.trn_pt[1]
    
    def get_nodes_in(self):
        return self.nodes_in
    
    def get_nodes_out(self):
        return self.nodes_out
    
    def get_surfaces(self):
        return self.surfaces
        

class Event_Handler(PriorityQueue): # event handler using priority queue
    
    current_time = 0
    event_id = 0 # event id is used to break ties in the priority queue
    
    def __init__(self):

        super().__init__()
        self.event_buffer:list[Event] = []
        self.edges:set[IntersectionEdge] = set()
        # self.initialize_from_graph(graph)


############# init functions ####################

    def load_from_graph(self, graph: IntersectionGraph):
        for edge in graph.edges:
            self.add_event(Event(edge))

    def add_K_event(self, grain1: Grain, grain2: Grain, graph:IntersectionGraph, substrate: Surface, sub_event=False): # find K event using linear programming
        res=self.find_K_event(grain1, grain2, substrate, sub_event=sub_event)
        
        #TODO: find the surfaces that are involved in the event
        # use slack variables to find the surfaces that are involved in the event
        c = 1
        surfaces = []
        for surface in grain1.get_surfaces():
            if res.slack[c] == 0:
                surfaces.append(surface)
            c += 1
        for surface in grain2.get_surfaces():
            if res.slack[c] == 0:
                surfaces.append(surface)
            c += 1
        if res.slack[c] == 0:
            surfaces.append(substrate)
        # else: # if the event is aerial, then add the substrate event too
        #     self.add_K_event(grain2, grain1, graph, substrate, sub_event=True)
        try:
            assert len(surfaces)==4
        except:
            print(len(surfaces))
            print(res.x, res.fun)
            raise ValueError

        nodes = [frozenset([surfaces[i], surfaces[j], surfaces[k]]) for i, j, k in [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]]
        nodes_in = set([node for node in nodes if node in list(graph.nodes)])
        event=Event(surfaces, nodes_in)
        self.add_event(event)

    def find_K_event(self, grain1: Grain, grain2: Grain, substrate: Surface, sub_event=False):
        c=np.array([0,0,0, 1]) # objective function: minimize time
        A_ub=np.array([[0,0,0, -1]]) # time >= 0
        b_ub=np.array([0])
        for surface in grain1.get_surfaces():
            A_ub=np.concatenate((A_ub, np.array([surface.v_mu])))
            b_ub=np.concatenate((b_ub, np.array([surface.b])))
        for surface in grain2.get_surfaces():
            A_ub=np.concatenate((A_ub, np.array([surface.v_mu])))
            b_ub=np.concatenate((b_ub, np.array([surface.b])))
        
        A_ub=np.concatenate((A_ub, -np.array([substrate.v_mu])))
        b_ub=np.concatenate((b_ub, -np.array([substrate.b])))
        
        if sub_event:
            return linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=[(None,None), (None,None), (0,0), (0,None)])
        else:
            return linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=[(None,None), (None,None), (0,None), (0,None)])


########## update functions ##########

    def new_event_from_edge(self, edge: IntersectionEdge):
        surfaces = set()
        for node in edge:
            surfaces = surfaces | node
        # print(list(surfaces))
        try:
            assert len(surfaces)==4
        except:
            raise ValueError(f"Edge {edge} does not have 4 surfaces")
        return Event(list(surfaces), set(edge))
    
    def add_event_from_edge(self, edge: IntersectionEdge):
        self.add_event(self.new_event_from_edge(edge))

    def new_event_from_spear_node(self, spear_node: IntersectionNode, surface: Surface, graph: IntersectionGraph, time: float):
        
        node_pos = np.append(Surface.nodeposition(spear_node, time), time)
        adj_surfaces = list(graph.get_adjacent_surfaces(surface).difference(spear_node))
        print(surface)
        print(spear_node)
        print(adj_surfaces)
        
        c=np.array([0,0,0, -1]) # objective function: maximize time
        A_ub=None
        b_ub=None
        A_eq=None
        b_eq=None

        for spear_surface in spear_node:
            if A_eq is None:
                A_eq=np.array([spear_surface.v_mu])
                b_eq=np.array([spear_surface.b])
            else:
                A_eq=np.concatenate((A_eq, np.array([spear_surface.v_mu])))
                b_eq=np.concatenate((b_eq, np.array([spear_surface.b])))

        for adj_surface in adj_surfaces:
            det = np.dot(node_pos, adj_surface.v_mu) - adj_surface.b
            if A_ub is None:
                if det<-EPS:
                    A_ub=np.array([adj_surface.v_mu])
                    b_ub=np.array([adj_surface.b])
                elif det>EPS:
                    A_ub=-np.array([adj_surface.v_mu])
                    b_ub=-np.array([adj_surface.b])
                else:
                    A_ub=np.array([[0.,0.,0.,0.]])
                    b_ub=np.array([1.])

            else:
                if det < -EPS:
                    A_ub=np.concatenate((A_ub, np.array([adj_surface.v_mu])))
                    b_ub=np.concatenate((b_ub, np.array([adj_surface.b])))
                elif det > EPS:
                    A_ub=np.concatenate((A_ub, -np.array([adj_surface.v_mu])))
                    b_ub=np.concatenate((b_ub, -np.array([adj_surface.b])))
                else:
                    A_ub=np.concatenate((A_ub, np.array([[0.,0.,0.,0.]])))
                    b_ub=np.concatenate((b_ub, np.array([1.])))

        res = linprog(c, A_eq=A_eq, b_eq=b_eq, A_ub=A_ub, b_ub=b_ub, bounds=[(None,None), (None,None), (0,None), (time, TIME_LIMIT)])

        if res.success:
            if res.fun == -TIME_LIMIT:
                return None
            
            count = 0
            surfaces = []
            for adj_surface in adj_surfaces:
                if res.slack[count] == 0:
                    surfaces.append(adj_surface)
                    print(adj_surface, -res.fun)
                count += 1

            if len(surfaces)==0:
                if res.x[2] == 0:
                    for adj_surface in adj_surfaces:
                        if adj_surface.id ==-1:
                            surfaces.append(adj_surface)

            surfaces.extend(list(spear_node))

            try:
                assert len(surfaces)==4
            except:
                print(len(surfaces))
                print(res.x, res.fun)
                # raise ValueError
                return None

            return Event(surfaces, set([spear_node]))
        
        else:
            raise ValueError("no valid event found")


    def handle_next(self, graph: IntersectionGraph, lateral_growth: bool = False):
        if self.empty():
            return False
        empty=False
        while self.is_same_event_all(self.event_buffer):
            if not self.empty():
                event=self.pop_event()
            else:
                empty=True
                break

            if (not self.is_lateral_growth(event)) and lateral_growth:
                self.add_event(event)
                print("lateral growth completed")
                return False
            # print(event)
            if self.is_valid_event(event, graph, check_nodes_out=False):
                self.event_buffer.append(event) # pop the next event from the queue and add it to the buffer
            else:
                # print("invalid event")
                continue # if the event is not valid, discard it and pop the next event from the queue
            
        # add back the last event in buffer to the queue (if not empty)
        if not empty:
            self.add_event(self.event_buffer.pop(-1), suppress_print=True)

        # update graph
        # this_event:Event = None
        
        this_event=self.event_buffer[0]
        nodes = self.event_buffer[0].get_nodes_in() | self.event_buffer[0].get_nodes_out()
        nodes_in = set()
        for node in nodes:
            if graph.has_node(node):
                nodes_in.add(node)
        nodes_out = nodes - nodes_in
        new_edges, spear_nodes = graph.update_graph(nodes_in, nodes_out)
        self.event_buffer=[]
        # if len(self.event_buffer) == 1: # initial K or H events
        #     print("initial K or H events")
        #     this_event = self.event_buffer.pop(0)
        #     # print (len(this_event.get_nodes_in()))
        #     # try:
        #     new_edges, spear_nodes = graph.update_graph(this_event.get_nodes_in(), this_event.get_nodes_out())
        #     # except:
        #     #     print("# of nodes:",len(list(graph.nodes)))

        # elif len(self.event_buffer) == 2: # 2 events should not happen. If it does re-check nodes_in
        #     this_event=self.event_buffer.pop(0)
        #     nodes = self.event_buffer[0].get_nodes_in() | self.event_buffer[0].get_nodes_out()
        #     nodes_in = set()
        #     for node in nodes:
        #         if graph.has_node(node):
        #             nodes_in.add(node)
        #     nodes_out = nodes - nodes_in
        #     new_edges, spear_nodes = graph.update_graph(nodes_in, nodes_out)
        #     self.event_buffer=[]

        # elif len(self.event_buffer) == 3: # 3 events corresponds to A (or K if time is reverted) type event. 3 nodes disappear and 1 node appears
        #     print("A (or K if time is reverted) type event")
        #     nodes_in = self.event_buffer[0].get_nodes_in() | self.event_buffer[1].get_nodes_in() | self.event_buffer[2].get_nodes_in()
        #     nodes_out= self.event_buffer[0].get_nodes_out() & self.event_buffer[1].get_nodes_out() & self.event_buffer[2].get_nodes_out()
        #     new_edges, spear_nodes = graph.update_graph(nodes_in, nodes_out)
        #     this_event=self.event_buffer.pop(0)
        #     self.event_buffer=[]

        # elif len(self.event_buffer) >= 4: # 4 events corresponds to Y (or X if time is reverted) type event. 4 nodes disappear and no node appears
        #     print("Y (or X if time is reverted) type event")
        #     nodes_in = self.event_buffer[0].get_nodes_in() | self.event_buffer[1].get_nodes_in() | self.event_buffer[2].get_nodes_in() | self.event_buffer[3].get_nodes_in()
        #     nodes_out= self.event_buffer[0].get_nodes_out() & self.event_buffer[1].get_nodes_out() & self.event_buffer[2].get_nodes_out() & self.event_buffer[3].get_nodes_out()
        #     new_edges, spear_nodes =graph.update_graph(nodes_in, nodes_out)
        #     this_event=self.event_buffer.pop(0)
        #     self.event_buffer=[]

        # else:
        #     raise ValueError("unknown event type:", len(self.event_buffer), "events in buffer")
        
        # update time
        
        self.update_time(this_event.get_time())
        print(this_event)

        # add new events
        # from edges

        # newedges=set(graph.edges).difference(self.edges)
        # self.edges=set(graph.edges)

        for edge in new_edges:
            new_event = self.new_event_from_edge(edge)
            if not self.is_same_event(new_event, this_event):
                # print("Potential new event:", new_event)
                # try:
                # if self.is_valid_event(new_event, graph, check_nodes_out=False):
                self.add_event(new_event)
                # except:
                #     print(new_event.surfaces)
                #     print(this_event.surfaces)
                #     print(Surface.intersect(new_event.surfaces))
                #     print(Surface.intersect(this_event.surfaces))
                # print("new event added")

        # from spear nodes
        for node in spear_nodes:
            surface=spear_nodes[node]
            new_event = self.new_event_from_spear_node(node, surface, graph, self.current_time)
            if new_event is None:
                print("no new event from spear node")
                continue
            if not self.is_same_event(new_event, this_event):
                # print("Potential new event:", new_event)
                # if self.is_valid_event(new_event, graph, check_nodes_out=False):
                self.add_event(new_event)
                # print("new event added")

        return True
            
############## key utils: event ##############

    def add_event(self, event: Event, suppress_print=False):
        if event.trn_pt[3] < self.current_time:
            if event.trn_pt[3] <= EPS:
                return
            event.type='past'
            print("past event!", event.trn_pt)
            # TODO: handle past event right away
            # self.handle_past_event(event)
            return
            
        else:
            event.type='future'
            # print("future event")
        try:
            if not suppress_print:
                print("Add Event:",len(event.get_nodes_in()),(event.trn_pt[3], self.event_id, event))
            
            # if self.has_event(event):
            #     print("Event already in queue!")

            #     if self.has_event_exact(event):
            #         print("Event already in queue!")
            #         return
            self.put((event.trn_pt[3], self.event_id, event)) # put the event in the queue. The priority is z of the event
            self.event_id += 1
        except:
            while not self.empty():
                print(self.get())
            raise ValueError("add_event failed")

    def has_event(self, event: Event):
        for e in self.queue:
            if self.is_same_event(event, e[2]):
                print("same event", event.get_nodes_in(), e[2].get_nodes_in())
                return True
        return False
    
    def has_event_exact(self, event: Event):
        for e in self.queue:
            if self.is_same_event_exact(event, e[2]):
                return True
        return False

    def pop_event(self) -> Event:
        # try:
        return self.get_nowait()[2]
        # except:
        #     while not self.empty():
        #         print(self.get_nowait()[1])
    
    def is_same_event_exact(self, event1: Event, event2: Event):
        if self.is_same_event(event1, event2):
            if event1.get_nodes_in() == event2.get_nodes_in():
                return True
        return False

    def is_same_event(self, event1: Event, event2: Event):
        # print(event1,event2)
        for surface in event1.surfaces:
            if surface not in event2.surfaces:
                return False
        return True
    
    def is_same_event_all(self, list_of_events: list[Event]):
        if len(list_of_events) <= 1:
            return True
        e1 = list_of_events[0]
        for e2 in list_of_events[1:]:
            if not self.is_same_event(e1, e2):
                return False
        return True

    def is_valid_event(self, event: Event, graph: IntersectionGraph, check_nodes_out=True): # check if the event is valid by checking if the nodes in the event are in the graph
        # nodes=list(graph.nodes)

        if event.get_z() < 0:
            return False
        # if event.get_x() < 0 or event.get_x() > SIZE[0]: # check if the event is out of bounds
        #     return False
        # if event.get_y() < 0 or event.get_y() > SIZE[1]:
        #     return False

        nodes_in= event.get_nodes_in()
        nodes_in_check = all([graph.has_node(node) for node in nodes_in])

        if not nodes_in_check:
            print("nodes in check failed", nodes_in)

        if check_nodes_out:
            nodes_out= event.get_nodes_out()
            nodes_out_check = not any([graph.has_node(node) for node in nodes_out])
            if not nodes_out_check:
                print("nodes out check failed", nodes_out)
        else:
            nodes_out_check = True

        # print(nodes_in_check,nodes_out_check)
        # if len(event.get_nodes_in()) == 1:
        #     IncidentSurface = set(event.surfaces) - nodes_in
        #     if len(graph.get_adjacent_surfaces(IncidentSurface)) == 0:
        #         return False
        return nodes_in_check and nodes_out_check
    
    def is_buffer_event(self, event: Event): # avoid adding the same event to the buffer
        return event.get_nodes_in() in [ev.get_nodes_in() for ev in self.event_buffer]

    def is_lateral_growth(self, event: Event):
        return event.get_z() == 0


############# utils: misc #############
    
    # def handle_past_event(self, event: Event, graph: IntersectionGraph):
    #     # handle past event
    #     if not self.is_valid_event(event, self.graph, checkbuffer=False):
    #         return False

    def update_time(self, time: float):
        self.current_time = max(time, self.current_time)

    def get_time(self):
        return self.current_time
