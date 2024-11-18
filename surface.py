#################################################################################
#       This file is part of the Surface Growth Simulation package.             #
#       This file contains the classes for the surfaces and grains.             #
#################################################################################

#### important notes ####
# 1. the surfaces are essentially STATIC. The 4-vector for each surface is fixed (velocity does not change).
# 2. the adjacency of the surfaces are DYANMIC. The adjacency of a surface is determined by the intersection graph.
# 3.



import numpy as np
from scipy.linalg import null_space
from scipy.spatial.transform import Rotation as R
from simulation_parameters import *


class Surface:

    def __init__(self, v_mu, r0, start_time=0, eps=1e-6):

        self.v_mu = np.array(v_mu) # 4-vector of the surface normal and the velocity (nx,ny,nz,-vn)
        self.b = np.dot(r0,v_mu) # dot product of the surface normal and the grain origin n*r0

        self.start_time = start_time
        self.eps = eps

        self.adjacent_surfaces = []
        # growing = True
        self.floating = False
        self.id = None
        self.grain_id = None
        self.killnode = None

    def __repr__(self) -> str:
        # return f'Surface: {self.v_mu}, {self.b}'
        return f'Surface: {self.grain_id}, {self.id}'

    # def intersect(self, surfaces: list):
    #     # calculate the intersection of a surface and a list of surfaces
    #     pass

    def intersect(surfaces: list):
        # calculate the intersection of a list of surfaces
        if len(surfaces) == 4:
            A=np.vstack([surface.v_mu for surface in surfaces])
            b=np.array([surface.b for surface in surfaces])
            x = np.linalg.solve(A,b)
            return x
        elif len(surfaces) == 0:
            return None
        elif len(surfaces) == 1:
            return None
        else:
            raise NotImplementedError('Intersection of less than 4 surfaces is not implemented yet.')

    def nodeposition(surfaces:list, t):
        # calculate the position of the node at time t
        assert len(surfaces) == 3
        A=np.vstack([surface.v_mu for surface in surfaces])
        A=np.vstack([A, np.array([0,0,0,1])])
        b=np.array([surface.b for surface in surfaces])
        b=np.append(b, t)
        x = np.linalg.solve(A,b)
        return x[:3]

    def area(self, adjacent_surfaces: list):
        # calculate the area of the surface 
        pass

    def stop_growing(self, stop_time):
        for surface in self.adjacent_surfaces:
            surface.adjacent_surfaces.remove(self)
        self.adjacent_surfaces = []
        self.growing = False
        self.stop_time = stop_time

    # def remove_surface(self):
    #     self.grain.surfaces.remove(self)
    #     for surface in self.adjacent_surfaces:
    #         surface.adjacent_surfaces.remove(self)
    #     self.adjacent_surfaces = []

 
class Substrate(Surface):
    def __init__(self, size=(1,1), eps=1e-6):
        super().__init__(v_mu=(0,0,1,0), r0=(0,0,0,0), eps=eps)

        self.growing = False
        self.size = size
        self.id = -1

class Surface_Hexagonal(Surface):
    def __init__(self, r0, orientation, id, grain_id, eps=EPS):
        
        # rotate the normal vector based on grain orientation
        r = R.from_euler('xyz', orientation, degrees=True)
        v_mu = np.append(r.apply(ORIENTATION[id]), -VN[TYPE[id]])
        # print(f"{self.id}: {v_mu}")
        super().__init__(v_mu=v_mu, r0=r0, eps=eps)
        self.id = id
        self.grain_id = grain_id


class Grain:

    id=0

    def __init__(self, orientation=(0,0,0), r0=(0,0,0,0), eps=EPS, grain_type='hexagonal'):
        '''
        orientation: (phi,theta,psi) in degrees
        r0: (x0,y0,z0,t0) in microns
        eps: tolerance for floating surfaces
        grain_type: 'hexagonal' only for now
        '''
        self.surfaces:list[Surface]=[]
        self.orientation = np.array(orientation) # orientation of the grain (phi,theta,psi)
        self.r0 = np.array(r0) # origin of the grain (x0,y0,z0=0,t0)
        self.eps = eps
        self.id = Grain.id
        Grain.id += 1
        
        if grain_type == 'hexagonal': # add the 8 surfaces of a hexagonal grain
            for i in range(8):
                self.surfaces.append(Surface_Hexagonal(r0, orientation, id=i, grain_id=self.id))
        else:
            raise ValueError('Grain type not implemented')

        if all(self.r0==0):
            print('Warning: grain origin is at (0,0,0,0)')

        if all(self.orientation==0):
            print('Warning: grain orientation is (0,0,0)')

        
        for i, surface in enumerate(self.surfaces): # set the adjacent surfaces in the grain
            for j in ADJACENCY[i]:
                surface.adjacent_surfaces.append(self.surfaces[j])
            
    def __repr__(self) -> str:
        return f'Grain: {self.r0}, {self.orientation}'

    def get_surface(self, index):
        return self.surfaces[index]
    
    def get_surfaces(self):
        return self.surfaces

    def resolve_substrate(self, substrate: Substrate): # resolve the initial intersection with the substrate
        
        # if the surface is below the substrate, it is not growing

        surfaces_to_stop = []
        surfaces_to_resolve = []
        triangular_surfaces = []
        adjacency_to_remove = set()
        for surface in self.surfaces:
            above=[]
            for adj_surface in surface.adjacent_surfaces:
                adj_adj_surfaces=set(adj_surface.adjacent_surfaces).intersection(set(surface.adjacent_surfaces))
                node_above = 0
                assert len(adj_adj_surfaces) == 2
                for adj_adj_surface in adj_adj_surfaces:
                    mat=np.vstack((adj_surface.v_mu,adj_adj_surface.v_mu,surface.v_mu, (0,0,0,1))) 
                    sol=np.linalg.solve(mat, (adj_surface.b,adj_adj_surface.b,surface.b,1))

                    if sol[2]>0: # if z increases with time (z>0 when t=1)
                        node_above += 1
                        # print(f"above: {surface.id}, {adj_surface.id}, {adj_adj_surface.id}")
                
                above.append(node_above)
                if node_above == 0:
                    adjacency_to_remove.add(frozenset((surface, adj_surface)))
                    # print(f"remove: {surface.id}, {adj_surface.id}")

            if max(above)==0:
                surfaces_to_stop.append(surface)
            elif min(above)==2:
                surface.floating = True
                # print(f"floating{surface.id}")
            elif max(above)==1:
                triangular_surfaces.append(surface)
                surfaces_to_resolve.append(surface)
            else:
                surfaces_to_resolve.append(surface)

        for surface, adj_surface in adjacency_to_remove:
            surface.adjacent_surfaces.remove(adj_surface)
            adj_surface.adjacent_surfaces.remove(surface)
        
        for surface_r in surfaces_to_resolve:
            # add the substrate as an adjacent surface
            # print("resolve:",surface_r.grain_id,surface_r.id)
            surface_r.adjacent_surfaces.append(substrate)
            substrate.adjacent_surfaces.append(surface_r)

        for surface in triangular_surfaces:
            assert len(surface.adjacent_surfaces) == 3
            surface.killnode = frozenset(surface.adjacent_surfaces)

        for surface_s in surfaces_to_stop:
            surface_s.stop_growing(0) # remove the surface not growing
            # self.surfaces.remove(surface_s)
            # print(f"stop growing!: {surface_s.id}")
            # self.surfaces.remove(surface) # remove the surface not growing
        # print(len(substrate.adjacent_surfaces))


                

    def volume(self):
        # calculate the volume of the grain
        pass



