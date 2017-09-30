# collective-motion-swarm
A group of programs simulating collective motion, crowd dynamics and various agent-based models. These include the classic Vicsek model, a topological Vicsek model, a repulsion-based Vicsek model, and crowd dynamics as modeled by a steepest-descent equation. 

#simulation.c
An implementation of the classic vicsek model. The program takes four command line arguments: (1) Board size (2) number of agents
(3) Noise coefficient and (4) iterations.  They should be provided in that order, as illustrated: 
./simulation BOARD_SIZE AGENT_NUMBER NOISE ITERATIONS
If no arguments are provided, then the program will use default values of 100, 10000, 0.1 and 1000, respectively. 
Outpoot is printed to standard output in CSV format. Each line represents the information for a single agent at a 
particular iteration and is given by: 
x_coord, y_coord, x_dir, y_dir, global order
global order is the same for all agents within a single iteration, and iterations are denoted by 
printing our an empty line.

#video-maker-classic-vicsek.ipynb: 
Jupyter notebook to visualize order parameter plots and create videos of simulations for simulation.c