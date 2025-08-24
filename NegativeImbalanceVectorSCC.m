% NegativeImbalanceVectorSCC takes in a digraph and a vertex id and returns
% a the same digraph with its edge weights updated such that the vertex
% imbalances of all the verteices are negative
%   Inputs:
%     G   – strongly‑connected digraph
%     src – source vertex index or name
%   Outputs:
%     G   – updated digraph, with the edge variable edge_weight updated

function G = NegativeImbalanceVectorSCC(G,src)

n = numnodes(G);

for t = 1:n    
    if t==src
        continue;
    else
        % Compute one path — for example shortest path
        [~,~,edgePath] = shortestpath(G, src, t);
        L = numel(edgePath);
        % for each edge on the path, starting from source
        for i = 1:L
            % weight contribution = (path_len − i + 1)
            w = (L-i+1);
            G.Edges.edge_weight(edgePath(i)) = G.Edges.edge_weight(edgePath(i)) + w;
        end
    end
end

% Compute imbalance for each node: outgoing minus incoming weight sums
for v = 1:n
    outgoing = sum(G.Edges.edge_weight(G.outedges(v)));
    incoming = sum(G.Edges.edge_weight(G.inedges(v)));
    G.Nodes.imbalance(v) = outgoing - incoming;
end

end
