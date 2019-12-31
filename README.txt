The different simulations are all found under directory CODE.
In order to run the different models:
original model:
run producer_vs_resistant.m when defining type_mutation = 'gamma'
multiple draws model:
run producer_vs_resistant.m when defining type_mutation = 'multiple_draws'
no overshoots model:
producer_vs_resistant_no_overshoots.m
reduced pressure model:
producer_vs_resistant_reduced_evolution.m
    
Other importanat functions : 
plot_fitness_map.m- creates a fitness map
single_droplet.m - The function that calculate the fitness of a single
pairwise competition


The directory "analayzeing results" contains some of the analyzing scripts used.
overshootsEffects.m - creates figure 9 and 10, looking at the probability
of antibiotic production collapse due to overshoots.
DrawAverageProduction.m - creates figure 12, comparing the average 
production of the original model to that of the mulatiple draws model.
similair codes were used for figures 6 and 16.
distribution_resistance.m - calculates the averaging effect in figure 11 
for different production windows.
plot_phen_on_fitness_mat.m - plot a single burst on a fitness map to see
the selection in each fixation event. Usually uses an arrow function taken from the matlab file exchange:
Erik Johnson (2019). arrow (https://www.mathworks.com/matlabcentral/fileexchange/278-arrow), MATLAB Central File Exchange. 
but can be modified to use standard functions.