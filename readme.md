# N-body

## About
This is intended to be a backend for other astronomical/cosmological code, ensuring an accurate and fast gravitational simulation is present without needing a new implementation each time. This library has Python bindings and is intended to be used from both C++ and Python.

## Planned Features
- More methods for time integration to allow for accuracy in large timesteps.
- More methods for faster particle simulation to allow for larger simulations.
- Parallelization of particle simulation to allow for effective use of multi-core machines and clusters.
- Specialized variations for different usecases, i.e. relativistic simulations, celestial mechanics (high precision simulations with N < 20), or simulations with collisions.
- Non-particle-based simulations for higher fidelity simulations.

And more!
