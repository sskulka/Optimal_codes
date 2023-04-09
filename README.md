# Optimal_codes
ME822_project: Hierarchical Path Planner for Off-road Vehicles

The attached matlab script 'hierarchical_planner_code_final.mlx' contains the implementation for the hierarchical path planner.

The planner is characterized into 2 algorithmic layers- 
Firstly, the global planner which is achieved through an implementation of adapting Dijkstra's planner as a Dynamic Programming algorithm.
The subsequent local planner is realized by implementing a finite receding-horizon Model Predictive Control algorithm.

The input to the planner is an excerpt of a geospatial map data tile available through the United States Geological Survey TopoView (https://www.usgs.gov/programs/national-geospatial-program/topographic-maps)
From the map data layers, the elevation and soil data layers have been extracted which are form a part of consideration for the cost function of the planners. 
Additionally, to simulate a stealth prone trajectory, a visibility layer consisting of binary information has been spoofed which is also incorporated into the cost function.

The Dijkstra's based DP computes a global plan by factoring in the elevation cost, soil texture cost & visibility cost and then iteratively calculating the cost-to-go 
of each of the adjacent 8 vertices of a particular vertex while going in the direction of the end goal.

For the purpose of this project, we did not segregate the input map data as low resolution for DP and high resolution for MPC, but instead did an interpolation on the 
low-resolution map data of DP for utilising it towards the MPC's finite-horizon. Based on the DP path co-ordinates, obstacles were spoofed into the elevation data 
along those co-ordiantes to make the MPC consider such points as dynamic obstacles that would be detected by sensors in a real-life application. This was done to check whether 
the obstacle's implication on cost function would make the MPC diverge from the global path which did happen (as it should).

The MPC optimization was done using fmincon for which inequality constraints were modeled based on kinematic considerations of the vehicle. The MPC cost function was similar 
to the DP cost function with the addition of an obstacle-based cost and a terminal cost in the form of the DP cost of that co-ordinate which served as a heuristic.
