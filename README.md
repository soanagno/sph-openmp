# SPH
Smoothed Particle Hydrodynamics (SPH) using OpenMP

![](sph_coarse_dense.png)

Figs: Snapshot of 20x10 m tank with 2000 particles (left), 20000 particles (right).


**Compilation**

with OpenMP: g++ sph.cpp -O3 -fopenmp

without OpenMP: g++ sph.cpp -O3


**Post-processing**

python3 post_sph.py


Further implementations:
