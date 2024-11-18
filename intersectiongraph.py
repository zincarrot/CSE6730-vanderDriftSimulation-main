#####################################################################################################
#       This file is part of the Surface Growth Simulation package.                                 #
#       This file contains the classes for the intersection graph, including nodes and edges.       #
#####################################################################################################

import networkx as nx
from surface import Surface, Grain
import itertools as it
import numpy as np
from simulation_parameters import *

IntersectionNode = frozenset[Surface] # 3-intersection of 3 surfaces. One-dimensional line in 4D spacetime
IntersectionEdge = frozenset[IntersectionNode] # 2-intersection of 2 surfaces. The key to finding turning points


class IntersectionGraph(nx.Graph):
    def __init__(self, surfaces: list[Surface], grains: list[Grain]):
        super().__init__()
        self.surfaces = surfaces
        self.grains = grains
        # self.nodes:set[IntersectionNode] = set()
        # self.edges:list[IntersectionEdge] = []
        self.build_graph()

    def detect_spear_node(self, node: IntersectionNode, time):
        time+=STEP_TIME
        neighbors = list(self.neighbors(node))
        pos=np.array([0.0,0.0,0.0])
        for neighbor in neighbors:
            pos+=self.get_node_position(neighbor,time)
        pos/=3 #len(neighbors)
        pos=np.append(pos,time)
        above=[]
        below=[]
        for surface in node:
            det = np.dot(surface.v_mu,pos)-surface.b
            if det>EPS:
                above.append(surface)
            elif det<-EPS:
                below.append(surface)
            else:
                return None

        if len(above)==1 and len(below)==2:
            return above[0]
        elif len(above)==2 and len(below)==1:
            return below[0]
        else:
            return None

    def edge_to_surfaces_all(self, edge: IntersectionEdge):
        return frozenset([surface for node in edge for surface in node])
    
    def edge_to_surfaces(self, edge: IntersectionEdge):
        edge_list=list(edge)
        return frozenset(edge_list[0].intersection(edge_list[1]))
    
    def get_adjacent_surfaces(self, surface: Surface):
        adjacent_surfaces = set()
        for edge in self.edges:
            if surface in self.edge_to_surfaces(edge):
                adjacent_surfaces= adjacent_surfaces.union(self.edge_to_surfaces(edge))

        if len(adjacent_surfaces) != 0:
            adjacent_surfaces.remove(surface)

        return adjacent_surfaces
    
    def get_boundary_edges(self, surface: Surface):
        boundary_edges = set()
        for edge in self.edges:
            if surface in self.edge_to_surfaces(edge):
                boundary_edges.add(edge)
        return boundary_edges
    
    def resolve_incident_edge(self, IncidentEdge: IntersectionEdge, TempEdge: IntersectionEdge, time, return_only=False):
        new_edges = []
        TempEdge = list(TempEdge)
        IncidentEdge = list(IncidentEdge)
        refLocation = self.get_node_position(list(IncidentEdge[0]),time+STEP_TIME)
        potLocation1 = self.get_node_position(list(TempEdge[0]),time+STEP_TIME)
        potLocation2 = self.get_node_position(list(TempEdge[1]),time+STEP_TIME)
        dist1 = np.linalg.norm(refLocation-potLocation1)
        dist2 = np.linalg.norm(refLocation-potLocation2)
        if dist1 < dist2:
            if not return_only:
                self.add_edge(TempEdge[0],IncidentEdge[0])
                self.add_edge(TempEdge[1],IncidentEdge[1])
            new_edges.append(frozenset((TempEdge[0],IncidentEdge[0])))
            new_edges.append(frozenset((TempEdge[1],IncidentEdge[1])))
        elif dist2 < dist1:
            if not return_only:
                self.add_edge(TempEdge[0],IncidentEdge[1])
                self.add_edge(TempEdge[1],IncidentEdge[0])
            new_edges.append(frozenset((TempEdge[0],IncidentEdge[1])))
            new_edges.append(frozenset((TempEdge[1],IncidentEdge[0])))
        else:
            raise ValueError("Nodes are static; can't determine new edge")
        if self.has_edge(IncidentEdge[0],IncidentEdge[1]):
            self.remove_edge(IncidentEdge[0],IncidentEdge[1])
        if self.has_edge(TempEdge[0],TempEdge[1]):
            self.remove_edge(TempEdge[0],TempEdge[1])
        # print("resolve:",new_edges)
        return new_edges
    
    def get_spear_node(self, node: IntersectionNode, edge: IntersectionEdge):
        surface = list(node.difference(self.edge_to_surfaces(edge))).pop()
        return {node: surface}

    def get_node_position(self, node: IntersectionNode, t: float):
        # calculate the position of the node at time t
        return Surface.nodeposition(list(node), t)

    def build_graph(self): # graph initialization
        visited = set()
        #Acquire node list from surfaces
        for grain in self.grains:
            for i, surface in enumerate(grain.surfaces):
                visited.add(surface)
                for j, adj_surface in enumerate(surface.adjacent_surfaces):
                    neighbor1 = adj_surface
                    Index1 = self.surfaces.index(neighbor1)
                    if neighbor1 not in visited:
                        for k in range(j+1, len(surface.adjacent_surfaces)):                        
                            neighbor2 = surface.adjacent_surfaces[k]
                            if neighbor2 not in visited and neighbor2 in self.surfaces[Index1].adjacent_surfaces:
                                NodeSurface = frozenset((surface, neighbor1, neighbor2))
                                pos=np.append(self.get_node_position(NodeSurface,0.0001),0.0001)
                                killnode=False
                                for nsurface in grain.surfaces:
                                    if nsurface.id==-1:
                                        continue
                                    if np.dot(nsurface.v_mu,pos)-nsurface.b>EPS:
                                        killnode=True
                                        break
                                if not killnode:
                                    self.add_node(NodeSurface)

                                # self.nodes.add(NodeSurface)
                            

        #Acquire edges list from nodes
        for node1, node2 in it.combinations(self.nodes,2):
                matches = len(node1.intersection(node2))
                # if matches >0:
                #     print(matches)
                if matches == 2:
                    # edge = frozenset([node1,node2])
                    # self.edges.append(edge)
                    self.add_edge(node1,node2)
                    # print([(surface.grain_id, surface.id) for surface in node1])
                    # print([(surface.grain_id, surface.id) for surface in node2])
                    # print("add:",[(surface.grain_id, surface.id) for surface in self.edge_to_surfaces(frozenset([node1,node2]))])

    def update_graph(self, nodes_in: set[IntersectionNode], nodes_out: set[IntersectionNode]) -> list[IntersectionEdge]:
        # edges to be added and removed are determined by the nodes_in and nodes_out.
       
        # nodes_in and nodes_out are set(frozenset(surfaces)); format to list(surfaces)
        new_edges = []
        spear_nodes = dict() # store the nodes that can cause collision events
        nodes = nodes_in.union(nodes_out)
        surfaces= set()
        for node in nodes:
            surfaces.update(node)
        time = Surface.intersect(list(surfaces))[3]  

        # H type:2 nodes disappear and 2 nodes appear
        if len(nodes_in) == len(nodes_out):
            print('H type')
        # Replace nodes
            self.add_nodes_from(nodes_out)
            # print(len(self.edges))

            tempEdge = frozenset((list(nodes_out)[0],list(nodes_out)[1]))
            if self.has_edge(list(nodes_out)[0],list(nodes_out)[1]):
                self.remove_edge(list(nodes_out)[0],list(nodes_out)[1])

            new_resolved_edge = []
            resolve= False
            for edge in self.edges:
                if self.edge_to_surfaces(edge)==(self.edge_to_surfaces(tempEdge)):
                    if edge == tempEdge:
                        continue
                    new_resolved_edge.extend(self.resolve_incident_edge(edge, tempEdge, time))
                    # new_edges.remove(tempEdge)
                    resolve = True
                    print("resolve:",new_resolved_edge)
                    break
            if not resolve:
                self.add_edge(list(nodes_out)[0],list(nodes_out)[1])
                new_edges.append(tempEdge)
            new_edges.extend(new_resolved_edge)
            # self.add_edges_from(new_resolved_edge)

            #replace edges
            for node_in in nodes_in:
                NeighborNodes = self.neighbors(node_in)
                for NeighborNode in NeighborNodes:
                    if NeighborNode in nodes_in:
                        continue
                    for node_out in nodes_out:
                        Match = len(NeighborNode.intersection(node_out))
                        if Match == 2:
                            self.add_edge(NeighborNode,node_out)
                            new_edges.append(frozenset((NeighborNode,node_out)))
            

            
            # print(len(new_edges))
            # print(len(self.edges))

            # check incident edges
            # EdgeSurfaces = list(list(nodes_in)[0].union(list(nodes_in)[1]))
            # time  = Surface.intersect(EdgeSurfaces)[3] 
            
            # print(len(new_edges))
            
            self.remove_nodes_from(nodes_in)
            # print(len(self.edges))

        # A/K type: 3 nodes disappear and 1 node appears
        elif len(nodes_in) -2 == len(nodes_out):
            print('A type')

            self.add_nodes_from(nodes_out)

            #replace edges
            for node_in in nodes_in:
                NeighborNodes = self.neighbors(node_in)
                for NeighborNode in NeighborNodes:
                    for node_out in nodes_out:
                        Match = len(NeighborNode.intersection(node_out))
                        if Match == 2:
                            self.add_edge(NeighborNode,node_out)
                            new_edges.append(frozenset((NeighborNode,node_out)))
            
            self.remove_nodes_from(nodes_in)

        # Y event: 4 nodes disappear and no node appears
        elif len(nodes_in) -4 == len(nodes_out):
            # self.nodes = self.nodes.difference(nodes_in)
            print('Y type')
            self.remove_nodes_from(nodes_in)
            

        elif len(nodes_in) + 4 == len(nodes_out):
            print('X type')
            self.add_nodes_from(nodes_out)

            #find existing edge and remove
            NodeList = list(nodes_out)
            EdgeSurfaces = list(NodeList[0].union(NodeList[1],NodeList[2],NodeList[3]))

            # add all possible edges in TempEdges
            TempEdges = []
            # OldEdgeSurfaces = []
            for node1, node2 in it.combinations(nodes_out,2):
                TempEdges.append(frozenset((node1,node2)))

            # time = Surface.intersect(EdgeSurfaces)[3]     

            # resolve incident edges
            for edge in self.edges:
                edgeSurfaces=self.edge_to_surfaces(edge)
                if edgeSurfaces.issubset(EdgeSurfaces):
                    for tempEdge in TempEdges:
                        if edgeSurfaces==(self.edge_to_surfaces(tempEdge)):
                            if edge == tempEdge:
                                continue
                            new_edges.extend(self.resolve_incident_edge(edge,tempEdge,time, return_only=False))
                            TempEdges.remove(tempEdge)
                            for node in tempEdge:
                                spear_nodes.update(self.get_spear_node(node, tempEdge))
            
            if len(new_edges) != 4:
                raise ValueError('Unexpected number of external edges')
            # self.add_edges_from(new_edges)
            self.add_edges_from(TempEdges)
            

        # Initial collision event: 1 nodes disappear and 3 node appears
        elif len(nodes_in) + 2 == len(nodes_out):
            print('K type')
            #replace nodes
            self.add_nodes_from(nodes_out)
            
            node_in = list(nodes_in)[0]
            NodeList = list(nodes_out)
            EdgeSurfaces = NodeList[0].union(NodeList[1],NodeList[2])
            IncidentSurface = list(EdgeSurfaces.difference(node_in))[0]
            TempEdges = [frozenset(comb) for comb in it.combinations(nodes_out,2)]
            # time = Surface.intersect(list(EdgeSurfaces))[3]

            for node in nodes_out:
                spear_nodes.update({node:IncidentSurface})


            for edge in self.get_boundary_edges(IncidentSurface):
                # print("edge",edge)
                if self.edge_to_surfaces(edge).issubset(EdgeSurfaces):
                    for tempEdge in TempEdges:
                        # print("compare",self.edge_to_surfaces(edge),(self.edge_to_surfaces(tempEdge)))
                        if self.edge_to_surfaces(edge)==(self.edge_to_surfaces(tempEdge)):
                            if edge == tempEdge:
                                continue
                            new_edges.extend(self.resolve_incident_edge(edge,tempEdge,time))
                            TempEdges.remove(tempEdge)
                            print("remove",tempEdge)
                            print("keep",edge)
                            for node in tempEdge:
                                spear_nodes.pop(node)
                            break
                    break
                self.add_edges_from(TempEdges)
                # new_edges.extend(TempEdges)
                # self.add_edges_from(new_edges)

            #replace edges
            NeighborNodes = self.neighbors(node_in)
            for NeighborNode in NeighborNodes:
                for node_out in nodes_out:
                    Match = len(NeighborNode.intersection(node_out))
                    if Match == 2:
                        self.add_edge(NeighborNode,node_out)
                        new_edges.append(frozenset((NeighborNode,node_out)))
                        break
            
            self.remove_nodes_from(nodes_in)

        else:
            raise ValueError("Nodes not matching any event")
        # print("new edges:",new_edges)

        for node_out in nodes_out:
            if node_out in spear_nodes:
                continue
            spear_surface = self.detect_spear_node(node_out, time)
            if spear_surface is not None and spear_surface.id != -1:
                spear_nodes.update({node_out:spear_surface})

        return new_edges, spear_nodes

