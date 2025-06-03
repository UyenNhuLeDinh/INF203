from src.ljts.box import Box

def main(): 
    Lx = 5.0
    Ly = 40
    rho_liquid = 0.73
    rho_vapor = 0.02
    
    our_box = Box(Lx, Ly, rho_liquid, rho_vapor)
    print(our_box._compartments)
    print('The total number of molecules: ', len(our_box._molecules))
    
    
    
if __name__ == "__main__":
    main() 
    
    
