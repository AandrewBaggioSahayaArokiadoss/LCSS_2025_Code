% Description : This function takes in a digraph, 
% Performs SCC decomposition
% For each SCC, it finds the vertices that have incoming edges
% from the other components

function X = SCC_vertices_with_inedges(G)

n = G.numnodes;
X = zeros(n,1);

% bins is the label for each vertex as to which SCC it belongs to
[bins,~] = conncomp(G);

% This vector says whether th ith node
% has incoming edges from the other components or not
for i = 1:n
    [~,inc_nodes] = inedges(G,i);
    if ~isempty(inc_nodes)
        if any(bins(i)~=bins(inc_nodes))
            X(i)=1;
        end
    end    
end
end
