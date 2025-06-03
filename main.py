from src.ljts.box import Box

def main(): 
    Lx = 5
    Ly = 40
    rho_liquid = 0.73
    rho_vapour = 0.02
    
    our_box = Box(Lx, Ly)
    our_box.compartment
    our_box.density_to_compartment(rho_liquid, rho_vapour)
    
    
    
    
if __name__ == "__main__" :
    main() 
    
    
