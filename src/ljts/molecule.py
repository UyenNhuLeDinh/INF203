
class Molecule:
    def __init__(self, x: float, y: float, z: float):
        self._m = 1.0
        self._x = x
        self._y = y
        self._z = z
        
    def position(self):
        "Method to return the position, updated when x,y,z coordinates change."
        self._position = (self._x, self._y, self._z)
    
        
    