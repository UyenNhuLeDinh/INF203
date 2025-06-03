import random
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
            
            for _ in range(num_molecules):
                x = random.uniform(0, self._Lx)
                y = random.uniform(y_min, y_max)
                z = random.uniform(0, self._Lz)
                self._molecules.append(Molecule(x,y,z))
        

            