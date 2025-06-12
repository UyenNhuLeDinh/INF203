from .box import Box

class DistortedBox(Box):
    def __init__(self, Lx, Ly, Lz, rho_liquid, rho_vapor, s_x, s_y, s_z):
        super().__init__(rho_liquid, rho_vapor)
        
        # Save original and scale factors
        self._Lx_orig = Lx
        self._Ly_orig = Ly
        self._Lz_orig = Lz

        self._s_x = s_x
        self._s_y = s_y
        self._s_z = s_z

        # Distorted box dimensions
        self._Lx = Lx * s_x
        self._Ly = Ly * s_y
        self._Lz = Lz * s_z

        self._dimensions = (self._Lx, self._Ly, self._Lz)
        self._scaled_molecules = [mol.clone() for mol in self.molecules]
        
        
    def compute_distorted_energy(self):
        # Rescale position:
        for mol in self._scaled_molecules:
            x, y, z = mol.position()
            mol.set_position((self._s_x * x, self._s_y * y , self._s_z * z))
            
        # Compute the distorted energy:
        distorted_energy = self.compute_potential()
        
        return distorted_energy
    
    
    def internal_energy_chaneg(self):
        undistorted_energy = Box.compute_potential()
        distorted_energy = self.compute_distorted_energy()
        
        return distorted_energy - undistorted_energy
    
    
    def surface_area_chaneg(self):
        distorted_A = 2 * self._Lx * self._Lf 
        undistorted_A = 2 * Box._Lx * Box._Lz
        
        return distorted_A - undistorted_A