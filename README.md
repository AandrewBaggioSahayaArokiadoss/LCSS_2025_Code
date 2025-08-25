# Network Synchronization of Lorenz Oscillators

MATLAB and Python code for simulating and visualizing the synchronization of a network of coupled Lorenz oscillators.  
Uses **graph theory** (weighted Laplacian matrices, vertex imbalances) to achieve and analyze network synchronization.

---

## Overview
- **MATLAB**: Implements the coupled Lorenz network, synchronization conditions, and simulation (`LCSS_synchronization_plot.m`).  
- **Python**: Jupyter Notebook for visualizing results (`synchronization_plot.ipynb`).  

---

## Files
- `LCSS_synchronization_plot.m` – Main simulation script (saves results to `sync_data.xlsx`)  
- `CoupledDynamics.m`, `SimulateCoupledSystems.m`, `LorenzOscillator.m` – Dynamics + solver  
- `SyncCouplingAssign.m`, `NegativeImbalanceVectorSCC.m`, `CycleBasisVector.m` – Coupling/imbalance assignment functions  
- `VertexImbalancePlot.m` – Visualization of graph imbalances  
- `synchronization_plot.ipynb` – Python visualization of synchronization error  

---

## How to Run
1. **MATLAB**  
   ```matlab
   run('LCSS_synchronization_plot.m')
