# SPH on OpenMP
Smoothed Particle Hydrodynamics (SPH) using OpenMP

![](sph_results.png)

Figs: Indicative snapshots of a 20x10 m tank with 2000 particles (left), 20000 particles (right). Timings: 2000 particles ~ 1 min for 4 s of real time, 20000 particles ~ 15 min for 4 s of real time, (time-step: 0.0001 s ran on an 8-core i5 CPU).



***Compilation***

with OpenMP: 
```
g++ sph.cpp -O3 -fopenmp
```

without OpenMP:
```
g++ sph.cpp -O3
```


***Post-processing***
```
python3 post_sph.py
```


***Further implementations***

- Pre-computed Kernel and Grad functions.
- 2nd order time integration method (Predictor-Corrector), for stability and larger time-step.
- Stencil for faster neighbor calculations (only once per particle).
- GPU enabled computations.
- Print results to binary.
