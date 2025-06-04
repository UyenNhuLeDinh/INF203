
class Molecule:
    def __init__(self, x: float, y: float, z: float):
        """Initialize a molecule as a point of mass.

        Args:
            x, y, z (float): Coordinates in a 3D-volume.
        
        Attributes:
            _m (float): Mass of the molecule, fixed unit mass of 1.
            _x (float): x-coordinate.
            _y (float): y-coordinate.
            _z (float): z-coordinate.
        """
        self._m = 1.0
        self._x = x
        self._y = y
        self._z = z
        
    def position(self):
        " Return the 3D coordinates of the molecule as a tuple."
        self._position = (self._x, self._y, self._z)
    
        
    