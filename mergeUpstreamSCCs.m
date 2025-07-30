% This function takes in a digraph and for each SCC in it,
% one of the followingactions are performed : 
% 1) If the SCC has no incoming arcs from other SCCs, it is directly input to the algorithm.
% 2) Otherwise, a subgraph containing the SCC and its upstream SCCs is formed. These upstream SCCs are 
%    merged into a single vertex (combining parallel arcs),
%    The adjacency matrices of such simplified subgraphs are combined in
%    block diagonal format and returned

function A = mergeUpstreamSCCs(G)
A_condense = adjacency(condensation(G));
if sum(~(sum(A_condense,2)>0))>1
    % A = [];
    error("More than one source (initial) SCC");
else
    % n = G.numnodes;
    [bins,binsize] = conncomp(G);
    A = zeros(max(binsize),max(bins)*max(binsize));
    G = SCC_vertices_with_inedges(G);
    bin_num = length(binsize);
    for i = 1:bin_num
        % indices of the vertices in the ith bin
        idx = abs(bins-i)<1;
        G_temp = subgraph(G,find(idx));
        A_temp = adjacency(G_temp).';
        % Adjacency vector for vertices of the SCC (G_temp) with incoming edges
        % from other SCCs
        in_edge_temp = G_temp.Nodes.is_inedge;
        if i<2
            indx = 1;
        else
            indx = find(sum(A)>0,1,'last')+1;
        end
        disp([indx:indx+binsize(i) indx:indx+binsize(i)])
        if any(in_edge_temp)
            A(indx:indx+binsize(i),indx:indx+binsize(i))=[zeros(1,binsize(i)+1);in_edge_temp A_temp].';
        else
            A(indx:indx-1+binsize(i),indx:indx-1+binsize(i))=A_temp.';
        end
    end
end
