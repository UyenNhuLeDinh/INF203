import random
import numpy as np
from .molecule import Molecule

class Box:
    def __init__(self, Lx: float, Ly: float, rho_liquid: float, rho_vapor: float):
        """Create the simulation cuboid with given dimensions and phase densities.

        Args:
            Lx (float): Length of the box in the x and z directions (box is square in xz).
            Ly (float): Length of the box in the y direction.
            rho_liquid (float): Initial density of molecules in the liquid compartment.
            rho_vapor (float): Initial density of molecules in the vapor compartment.
            
        Raises:
            ValueError: If any of dimensions and densities are non-positive.
        """
        # Lengths of the box:
        self._Lx = Lx
        self._Ly = Ly
        self._Lz = Lx
        
        # Density at different phases:
        self._rho_liquid = rho_liquid
        self._rho_vapor = rho_vapor
        
        if Lx <= 0 or Ly <= 0 or rho_liquid < 0 or rho_vapor < 0:
            raise ValueError("Box dimensions and densities must be positive.")
        
        # Volume of the box:
        self._V = self._Lx * self._Ly * self._Lz
        
        # Cut-off radius for Lennard-Jones (LJ) Potential:
        self._cutoff = 2.5
        
        # Dictionary stores volumes of three compartments:
        self._compartments = {}
        self.compartment()
        
        # Dictionary stores densities of three compartments:
        self._densities = {}
        
        # List stores all molecules' coordinates in the box:
        self._molecules = []
        
        
    def compartment(self):
        """Define compartment volumes and store them in self._compartments.
        Note: 
            The box is divided along the y-axis into three regions: 
            lower (40%), middle (20%), and upper (40%).
        """
        self._compartments = {
        "lower": 0.4 * self._V,
        "middle": 0.2 * self._V,
        "upper": 0.4 * self._V
        }
        
        
    def compartment_densities(self):
        """Assign densities to compartments and store them in self._densities.
        """
        self._densities = {
            "lower": self._rho_vapor,
            "middle": self._rho_liquid,
            "upper": self._rho_vapor
            }
        
        
    def populate_box(self):
        """Populate the simulation box with molecules randomly placed in each compartment.
            The number of molecules is the product of densities and compartment volumes.
        Note:
            The function prints the total number of molecules in the box.
        """
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
            
            # Random generation of the molecules' positions:
            for _ in range(num_molecules):
                x = random.uniform(0, self._Lx)
                y = random.uniform(y_min, y_max)
                z = random.uniform(0, self._Lz)
                self._molecules.append(Molecule(x,y,z))
                
        print('The total number of molecules: ', len(self._molecules))       
                
                
    def compute_potential(self):
        """Compute the total potential energy of the box using the pairwise 
        Lennard-Jones truncated-shifted (LJTS) energy.
        
        Returns:
            float: The total potential energy of the box.
            
        Raises:
            ValueError: If no molecules are initialized in the box.     
        
        Notes:
            - Assumes all molecules have unit mass of 1.
            - Only pairs within a cutoff radius of 2.5 are considered.
            - The LJ potential is shifted to ensure continuity at the cutoff.     
        """
        # Get positions of all molecules:
        positions = [(mol._x, mol._y, mol._z) for mol in self._molecules]
            
        # Number of molecules:
        N = len(positions)
        
            
        def LJ_potential(r_ij: float):
            """Calculate Lennard-Jones potential energy between a pair of molecules.
            
            Args:
                r_ij (float): The distance between molecule i and j.
                
            Returns:
                float: Pairwise Lennard-Jones energy with truncation and shift.
            """
            shift = 4 * ((1 / self._cutoff ** 12) - (1 / self._cutoff** 6))
            
            if r_ij ** 2 < self._cutoff ** 2:
                inv_r6 = r_ij ** -6
                inv_r12 = inv_r6 ** 2
                return 4 * (inv_r12 - inv_r6) - shift
            else:
                return 0
            
            
        def pbc_distance(dx: float, L: float):
            """Apply minimum image convention for periodic boundary conditions.
            
            Args:
                dx (float): Displacement along one coordinate axis (xj - xi).
                L (float): Length of the periodic box in that direction.
                
            Returns:
                float: The shortest distance between two molecules' image on a particular direction.
                
            Notes:
                round(dx/L) gives the nearest image of the second molecule.
            """
            return dx - round(dx / L) * L
        
        
        # Iterations to cummulate the potential energy:
        E_pot = 0
        for i in range(N):
            xi, yi, zi = positions[i]
            
            for j in range (i+1, N):
                xj, yj, zj = positions[j]
                
                dx = pbc_distance(xj - xi, self._Lx) 
                dy = pbc_distance(yj - yi, self._Ly)
                dz = pbc_distance(zj - zi, self._Lz)
                
                r_ij = np.sqrt(dx**2 + dy**2 + dz**2)
                E_pot += LJ_potential(r_ij)
        
        print('The total potential energy:', E_pot)
        return E_pot    
                

        

            