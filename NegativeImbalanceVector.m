function G = NegativeImbalanceVector(G, src, a)
% assignWeightsAndImbalance  — Add edge/node variables and compute weights & imbalance
%   G = assignWeightsAndImbalance(G, src, a)
%   Inputs:
%     G   – strongly‑connected digraph
%     src – source vertex index or name
%     a   – scaling constant
%   Outputs:
%     G   – updated digraph, with:
%           G.Edges.weight   – numeric weight per edge
%           G.Nodes.imbalance – numeric imbalance per node
%   And a plot of the graph with edge‑labels showing weight and node‑labels showing imbalance.

% Add custom edge and node attributes, initialized to zero
G.Edges.weight = zeros(numedges(G),1);
G.Nodes.imbalance = zeros(numnodes(G),1);


% For each target node (except src), find one directed path
n = numnodes(G);

G.Nodes.names = (1:n).';

for t = 1:n
    % if isequal(t,src) || isequal(findnode(G,src), findnode(G,t)), continue; end
    
    if t==src
        continue;
    else
        % Compute one path — for example shortest path
        [~,L,edgePath] = shortestpath(G, src, t);
        disp(edgePath)
        % L = numel(edgePath);
        % for each edge on the path, starting from source
        for i = 1:L
            % weight contribution = (path_len − i + 1)
            w = (L-i+1);
            G.Edges.weight(edgePath(i)) = G.Edges.weight(edgePath(i)) + w;
        end
    end
end
% disp(G.Edges.weight.')
% Scale edge weights by a
G.Edges.weight = a * G.Edges.weight;

% Compute imbalance for each node: outgoing minus incoming weight sums
for v = 1:n
    outgoing = sum(G.Edges.weight(G.outedges(v)));
    incoming = sum(G.Edges.weight(G.inedges(v)));
    G.Nodes.imbalance(v) = outgoing - incoming;
end

% Plot the graph
plot(G, 'EdgeLabel', G.Edges.weight, 'NodeLabel', strcat(string(G.Nodes.names),"(",string(G.Nodes.imbalance),")"));
title(sprintf('Digraph with weight-scaled by %g and node imbalance', a));

end
