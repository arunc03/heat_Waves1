Heat Wave Propagation through a Cuboidal Block
Arun Chonkaria

A 3D simulation propagation of different input thermal waves through a cuboid block. 
The block can be configured to different dimensional sizes, thermal conductivity, density and holes can be drilled into the faces. 

CONTENTS OF ZIP FILE:
- README.MD
- Simulation.py - contains code associated with the Simulation class, i.e. running the animation and plotting temperature variational plots.
- Fourier.py - contains code associated with the FourierAnalysis class, i.e. computing thermal diffusivity using Fourier methods.
- Input_Heat_Waves.py - contains associated code with the HeatingFunction class, i.e. to generate different waveforms on one face.
- Calculations.py - contains code associated code with the HeatEquationSolver class, i.e. to solve the wave equation using the Finite Differences method.
- Main.py - Main code file for exempler simulations.
- Tests.py - code used for testing through development.
- block_three_holes.pkl AND block_one_hole.pkl - saved objects run for 100s, as shown in the Main.py file.

BASIC REQUIREMENTS:
The packages used are Numpy, Scipy, Pickle and Matplotlib.

HOW TO RUN:
To run the simulation, use the Main.py file and run each block.

It is adviced that all animations and plots should be run in the 'Inline' graphical backend (using Spyder).

It is important that the filenames are not changed.