My code is split up into sections which handle the sub-question of the project

Analysis_Methods contains general functions like the NLL function that i call from other modules

1. Parabolic routines
Parabolic_v (validation)
Parabolic_f (final) for application to data
2. Univariate routines
Univariate_v (validation)
Univariate_f (final) for application to data
3. Simulatous 2D Newton routines
Sim_2D_v (validation)
Sim_2D_f (final) for application to data
4. Simultaneous 3D Newton routine
Newton_3D (final) for application to data
5. Simuluated amnealing routine
Annealing_v (validation)
Annealing_f (final) for application to data

One simply needs to unpack the contents of the zip into a single folder and then run 
the routines in any order to view the results from that routine, or the validation exercises performed on each method

->Histograms.py plots some histograms that were used in the report but are not essential to view

->Annealing routines take a while to run without GPU, but all other routines will run ~instantly