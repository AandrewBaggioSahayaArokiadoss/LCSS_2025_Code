The matlab codes used for checking the computations of my first submission to LCS-2025
1) SCC_vertices_with_inedges – Adds a variable "is_inedge" to each vertex of the input digraph, which indicates whether a vertex has incoming edges from other SCCs.
2) mergeUpstreamSCCs - For each SCC in the given digraph, it merges all its upstream SCCs into a single digraph and returns all the adjacency matrices of these digraphs into a matrix of size maxbinsize X maxbinsize*no_of_SCCs
3) coupledDynamics - Given the state vector for the entire dynamical network, gives dX/dt for the entire network
4) SimulateCoupledSystems -  Gives a list of state vectors for the entire time span
5) NegativeImbalanceVectorSCC - Assigns a negative imbalance vector to a strongly connected digraph and plots it with the edge weights and imbalances
6) NegativeImbalanceVector - Uses NegativeImbalanceVectorSCC to assign vertex imbalance vector for digraphs that are either strongly connected or that have two SCCs with one of them being a vertex devoid of incoming edges
7) LCSS_synchronization_plot - Creates the plot for the synchronization of a dynamical network of Lorentz oscillators, plots the image and stores the image data in an excel sheet of the name "sync_data.xlsx"
8) CycleBasisVector - finds the sum of cycle basis with positive entries of the input strongly connected digraph
9) synchronization_plot - Uses sync_data to create the plot
