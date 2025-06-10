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
        i = random.randint(0, len(self._box._molecules) - 1)
        mol = self._box._molecules[i]
        
        old_pos = mol.position()
        
        # Compute its potential energy before the test action:
        old_energy = self._box.molecular_energy(i)
        
        # Random displacement in [-b, b] for each dimension with PBC conditions:
        old_x, old_y, old_z = old_pos
        Lx, Ly, Lz = self._box.dimensions()
        
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
        
        
    def compute_total_energy(self):
        energy = 0
        for i in range(len(self._box._molecules)):
            energy += self._box.molecular_energy(i)
        return energy / 2  #Avoid counting pair-wise interactions twice
        
        
    def MC_step(self):
        """Perform one MC step containing N test actions (a sweep)
        """
        for _ in range(len(self._box._molecules)):
            self.metropolis_algorithm()
            
    def run_simulation(self):
        for step in range(self._steps):
            self.MC_step()
            
            # Collect output every 20 steps:
            if step % 20 == 0:
                total_energy = self.compute_total_energy()
                print(f"Step {step}: Potential energy of the box = {total_energy}")
                self._energy_log.append(total_energy)
                
    def energy_analysis(self):
        avg_energy = sum(self._energy_log) / len(self._energy_log)
        print(f"The average energy over {self._steps} steps : {avg_energy}")
        
        plt.plot(self._energy_log)
        plt.xlabel("Every 10 steps")
        plt.ylabel("Potential Energy")
        plt.title("Energy of the Box during Monte Carlo Simulation")
        plt.show()
        
                
            
            
    
            
        