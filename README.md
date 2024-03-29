# Intrinsic-Interface-and-Fractal-Dimension
Manuscript title: The Intrinsic Fragility of the Liquid-Vapor Interface: A Stress Network Perspective
https://doi.org/10.1021/acs.langmuir.2c00201

Codes to reproduce the data of this manuscript: 

<b>Flowmol MD solver</b> was used for the simulations. This is an open source platform developed by one of the authors (E. R. Smith) and can be found here: https://github.com/edwardsmith999/flowmol#flowmol

<b>intrinsic.in</b> : input file to reproduce the data.
Steps: 
1. Initially, run in NVT ensemble with INTRINSIC_INTERFACE turned OFF (to speed up the process).
2. Using the final state of step 1 as restart file, run NVE ensembles with INTRINSIC_INTERFACE turned OFF (to speed up the process) for sufficient amount of steps.
3. Using the final state of step 2 as restart file, run NVE ensemble with INTRINSIC_INTERFACE turned ON. When this option is turned on, the MPI topology should be 1x1x1.Data generated with INTRINSIC_INTERFACE turned ON were used for the surface stress analysis.

<b>NetworkMaps.py</b> : generates the stress network maps at the intrinsic surface with moelcular positions overlaid.
