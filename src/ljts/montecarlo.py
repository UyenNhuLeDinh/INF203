import random
import numpy as np
import math
import matplotlib.pyplot as plt
from .box import pbc_distance

class MonteCarloSimulator:
    def __init__(self, T, b, box, steps):
        self._T = T
        self._b = b
        self._box = box
        self._steps = steps
        self._accepted_moves = 0 
        self._total_moves = 0
        self._energy_log = []
        
        
    def metropolis_algorithm(self):
        """One trial move to a random molecule i
        """
        
        # Select a random molecule i:
        i = random.randint(0, self._box.num_molecules - 1)
        mol = self._box.molecules[i]

        
        old_pos = mol.position()
        
        # Compute its potential energy before the test action:
        old_energy = self._box.molecular_energy(i)
        
        # Random displacement in [-b, b] for each dimension with PBC conditions:
        old_x, old_y, old_z = old_pos
        Lx, Ly, Lz = self._box.dimensions
        
        new_x = pbc_distance(old_x + random.uniform(-self._b, self._b), Lx)
        new_y = pbc_distance(old_y + random.uniform(-self._b, self._b), Ly)
        new_z = pbc_distance(old_z + random.uniform(-self._b, self._b), Lz)
        
        mol.set_position((new_x, new_y, new_z))
        
        # Compute new potential energy:
        new_energy = self._box.molecular_energy(i)
        
        delta_E = new_energy - old_energy
        
        # Evaluation for accepting or rejecting the move:
        if delta_E < 0 : 
            # Accept if energy decrease
            self._accepted_moves += 1
        
        else:
            # Accept with probability exp(-ΔE / T) (Boltzman factor) with increasing energy:
            if random.random() < math.exp(- delta_E / self._T):
                self._accepted_moves += 1
            else:
                mol.set_position(old_pos)
        self._total_moves += 1
        
        
    def MC_step(self):
        """Perform one MC step containing N test actions (a sweep)
        """
        for _ in range(self._box.num_molecules):
            self.metropolis_algorithm()
            
            
    def run_simulation(self):
        for step in range(self._steps):
            self.MC_step()
            
            # Collect output every 20 steps:
            if step % 10 == 0:
                total_energy = self._box.compute_potential()
                print(f"Step {step}: Potential energy of the box = {total_energy}")
                self._energy_log.append(total_energy)
                
                
    def energy_analysis(self):
        # Skip the first few lines:
        sample_data = self._energy_log[2:]
        avg_energy = sum(sample_data) / len(sample_data)
        print(f"The average energy over {self._steps} steps : {avg_energy}")
        
        plt.plot(sample_data)
        plt.axhline(avg_energy, color='red', linestyle='--', label='Avg. Energy')
        plt.xlabel("Every 10 steps")
        plt.ylabel("Potential Energy")
        plt.title("Energy of the Box during Monte Carlo Simulation")
        plt.legend()
        plt.show()
        
                
            
            
    
            
        