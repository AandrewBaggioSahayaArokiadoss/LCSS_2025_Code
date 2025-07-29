% This function takes in a digraph and for each SCC in it,
% one of the followingactions are performed : 
% 1) If the SCC has no incoming arcs from other SCCs, it is directly input to the algorithm.
% 2) Otherwise, a subgraph containing the SCC and its upstream SCCs is formed. These upstream SCCs are 
%    merged into a single vertex (combining parallel arcs),
%    and the simplified subgraph is input to the algorithm

function A = mergeUpstreamSCCs(G)
A = null;
A_condense = adjacency(condensation(G));
if sum(~(sum(A_condense,2)>0))>1
    disp("More than one source (initial) SCC");
else
    n = G.numnodes;
    [bins,binsize] = conncomp(G);
    G = SCC_vertices_with_inedges(G);
    bin_num = length(binsize);
    A = zeros(2*n-1,2*n-1);
    for i = 1:bin_num
        % indices of the vertices in bin i
        indx = abs(bins-i)<1;
        G_temp = subgraph(G,find(indx));
        A_temp = adjacency(G_temp).';
        % Contains which vertices of the SCC (G_temp) had incoming edges
        % from other SCCs
        in_edge_temp = G_temp.Nodes.is_inedge;
        if any(X(indx))
            inc_vertices = X(indx)>0
            A=[zeros(1,binsize+1); A_temp]
        else

        end
    end
end
