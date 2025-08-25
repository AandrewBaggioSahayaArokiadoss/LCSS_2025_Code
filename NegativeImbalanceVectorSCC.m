% NegativeImbalanceVectorSCC assigns edge weights in a strongly connected
% component so that all vertex imbalances are negative.
%
%   Inputs:
%     G   – strongly connected digraph
%     src – source vertex index (or name)
%
%   Outputs:
%     G – updated digraph with edge_weight property modified
%
function G = NegativeImbalanceVectorSCC(G,src)

    %% Setup
    n = numnodes(G);

    %% Assign edge weights along paths from source
    for t = 1:n
        if t == src
            continue;   % skip source node
        end

        % Compute a path (use shortest path)
        [~,~,edgePath] = shortestpath(G,src,t);
        L = numel(edgePath);

        % For each edge along the path, increment weight
        for i = 1:L
            % Contribution decreases along the path
            w = (L - i + 1);
            G.Edges.edge_weight(edgePath(i)) = ...
                G.Edges.edge_weight(edgePath(i)) + w;
        end
    end

    %% Compute vertex imbalances (outgoing - incoming)
    for v = 1:n
        outgoing = sum(G.Edges.edge_weight(G.outedges(v)));
        incoming = sum(G.Edges.edge_weight(G.inedges(v)));
        G.Nodes.imbalance(v) = outgoing - incoming;
    end
end
