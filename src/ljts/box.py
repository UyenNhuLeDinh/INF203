import random
import math
from .molecule import Molecule

class Box:
    def __init__(self, Lx: float, Ly: float, rho_liquid: float, rho_vapor: float):
        # Lengths of the box:
        self._Lx = Lx
        self._Ly = Ly
        self._Lz = Lx
        
        # Density at different phases:
        self._rho_liquid = rho_liquid
        self._rho_vapor = rho_vapor
        
        # Volume of the box:
        self._V = self._Lx * self._Ly * self._Lz
        
        # Dictionary stores volumes of three compartments:
        self._compartments = {}
        self.compartment()
        
        # Dictionary stores densities of three compartments:
        self._densities = {}
        
        # List stores all molecules' coordinates in the box:
        self._molecules = []
        
        self.compartment()
        
        
    def compartment(self):
        self._lower = 0.4 * self._V
        self._middle = 0.2 * self._V
        self._upper = 0.4 * self._V
        
        self._compartments = {
        "lower": self._lower,
        "middle": self._middle,
        "upper": self._upper
        }
        
        
    def compartment_densities(self):
        self._densities = {
            "lower": self._rho_vapor,
            "middle": self._rho_liquid,
            "upper": self._rho_vapor
            }
        
        
    def populate_box(self):
        # Call assigning density method:
        self.compartment_densities()
        
        # Loop through each items in self._densities (region and density):
        for region, density in self._densities.items():
            volume = self._compartments[region]
            # Compute total number of molecules in each compartment:
            num_molecules = int(density * volume)
            
            # Populate molecules in specific region:
            if region == "lower":
                y_min, y_max = 0, 0.4 * self._Ly
            elif region  == "middle":
                y_min, y_max = 0.4 * self._Ly, 0.6 * self._Ly
            else: # "upper" region
                y_min, y_max = 0.6 * self._Ly, self._Ly
            
            # Random configuration generation on positions of molecules:
            for _ in range(num_molecules):
                x = random.uniform(0, self._Lx)
                y = random.uniform(y_min, y_max)
                z = random.uniform(0, self._Lz)
                self._molecules.append(Molecule(x,y,z))
                
    def compute_potential(self):
        # Different molecule positions on each run by random:
        random.seed()
        
        self.populate_box()
        # Get positions of all molecules:
        for mol in self._molecules:
            positions = [(mol._x, mol._y, mol._z)]
            
            # Number of molecules:
            N = len(positions)
            
        def LJ_potential(r_ij):
            if r_ij < 2.5:
                inv_r6 = (1 / r_ij**6)
                inv_r12 = inv_r6 ** 2
                return 4 * (inv_r12 - inv_r6) + 0.01631689
            else:
                return 0
            
        def distance():
            pass
            
        E_pot = 0
        for i in range(N):
            xi, yi, zi = positions[i]
            for j in range (i+1, N):
                xj, yj, zj = positions[j]
                
                dx = xj - xi
                dy = yj - yi
                dz = zj - zi
                
                r_ij = math.sqrt(dx**2 + dy**2 + dz**2)
                E_pot += LJ_potential(r_ij)
        
        print('The total potential energy:', E_pot)
        return E_pot    
                

        

            