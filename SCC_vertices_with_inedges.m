% Description : This function takes in a digraph, 
% Performs SCC decomposition
% For each SCC, it finds the vertices that have incoming edges
% from the other components

function G = SCC_vertices_with_inedges(G)
n = G.numnodes;

% bins is the label for each vertex as to which SCC it belongs to
[bins,~] = conncomp(G);

% This is the new property of each node as to whether
% it has incoming edges from the other components or not
G.Nodes.is_inedge = zeros(n,1);
for i = 1:n
    [~,inc_nodes] = inedges(G,i);
    if ~isempty(inc_nodes)
        if any(bins(i)~=bins(inc_nodes))
            G.Nodes.is_inedge(i)=1;
        end
    end    
end
end