class Molecule:
    def __init__(self, position):
        """_summary_

        Args:
            position (_tuple_): Coordinates of a point mass in 3 dimensional space.
        """
        self._mass = 1
        self.position = position