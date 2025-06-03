from .molecule import Molecule

class Box:
    def __init__(self, Lx, Ly):
        self._Lx = Lx
        self._Ly = Ly
        self._Lz = Lx
        self._V = self._Ly * self._Lx **2
        self._molecules = []
        self._compartments = {}
        self._densities = {}
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
        
    def density_to_compartment(self, rho_liquid, rho_vapour):
        self._densities = {
            "lower": rho_vapour,
            "middle": rho_liquid,
            "upper": rho_vapour
            }
        
    def compute_nr_molecule(self):
        total_molecules = 0
        for region, density in self._densities.items():
            volume = self._compartments[region]
            nr_molecules = density * volume
            total_molecules += nr_molecules
        print("The total number of molecules in the box:", total_molecules)

            