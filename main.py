import random
import time
from src.ljts.box import Box

def main(): 
    random.seed(time.time())
    
    Lx = 5.0
    Ly = 40
    rho_liquid = 0.73
    rho_vapor = 0.02
    
    our_box = Box(Lx, Ly, rho_liquid, rho_vapor)
    
    # Calculate the number of molecules in the box:
    our_box.populate_box()
    
    # Calculate the random generated Potential Energy of the box:
    our_box.compute_potential()
    
    
if __name__ == "__main__":
    main() 
    
    
