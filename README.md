The matlab codes used for checking the computations of my first submission to LCS-2025
1) SCC_vertices_with_inedges – Adds a variable "is_inedge" to each vertex of the input digraph, which indicates whether a vertex has incoming edges from other SCCs.
2) random_strong_digraph – Creates a strongly connected digraph with n vertices
3) mergeUpstreamSCCs - For each SCC in the given digraph, it merges all its upstream SCCs into a single digraph and returns all the adjacency matrices of these digraphs in block diagonal format
