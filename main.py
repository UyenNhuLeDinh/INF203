import random
import time
from src.ljts.box import Box
from src.ljts.montecarlo import MonteCarloSimulator

def main(): 
    random.seed(time.time())
    
    Lx = 5.0
    Ly = 40
    rho_liquid = 0.02
    rho_vapor = 0.02
    
    our_box = Box(Lx, Ly, rho_liquid, rho_vapor)
    
    # Calculate the number of molecules in the box:
    our_box.populate_box()
    print("Number of molecules:", our_box.num_molecules)
    
    # Calculate the random generated Potential Energy of the box:
    our_box.compute_potential()
    
    # Monte Carlo simulator:
    sim = MonteCarloSimulator(T = 0.8, b = 1/8, box=our_box, steps = 10000)
    sim.run_simulation()
    
    #print(f"Acceptance ratio: {sim._accepted_moves / sim._total_moves}")
    
    # Show results:
    sim.energy_analysis()
    
    
if __name__ == "__main__":
    main() 
    
    
