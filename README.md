# Cell Pattern Prediction using Theoretical and Simulation Models

This repository contains tools to predict the preferred patterns of cell types in an array based on the adjacency matrix and the minimal eigenvalue's eigenvector. The computations are done using Matlab.

### Overview

1. **Adjacency Matrix Calculation**: To use our method, we first compute the adjacency matrix for an array of cells.
   
2. **Theoretical Prediction**: After computing the adjacency matrix, we calculate the eigenvector corresponding to the minimal eigenvalue. This allows us to predict the preferred pattern of cell types in the array.

3. **Verification with Simulations**: The theoretical predictions are verified through simulations, using quasistatic variation of parameters.

### Key Files

- **nn_theoretical_predictions.m**: This script makes theoretical predictions for a square array of cells with nearest and next-nearest neighbor couplings.
- **nn_quasi_sim.m**: This script runs simulations to verify the predictions made in `nn_theoretical_predictions.m`.

- **filo_theoretical_predictions.m**: This script makes theoretical predictions for a square array of cells with long-range signaling, potentially via filopodia or paracrine signaling.
- **filo_quasi_sim.m**: This script runs simulations to verify the predictions made in `filo_theoretical_predictions.m`.
