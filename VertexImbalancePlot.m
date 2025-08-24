% VertexImbalancePlot takes in a digraph and returns
% a the vertex imbalance vector of the digraph and
% plots the digraph with its edge weights and the respective 
% verteximbalances
%   Inputs:
%     G   – strongly‑connected digraph
%   Outputs:
%     vertImb   – vertex imbalance vector

function vertImb = VertexImbalancePlot(G)
    % Validate input
    assert(isa(G, 'digraph'), 'Input must be a digraph object.');

    % Number of nodes
    n = numnodes(G);

    % Extract edge list and weights
    s = G.Edges.EndNodes(:,1);  % sources
    t = G.Edges.EndNodes(:,2);  % targets
    w = G.Edges.edge_weight;         % edge weights

    % Initialize imbalance vector
    vertImb = zeros(n,1);

    % Accumulate outgoing and incoming weights
    for k = 1:numel(w)
        vertImb(s(k)) = vertImb(s(k)) + w(k);  % outgoing
        vertImb(t(k)) = vertImb(t(k)) - w(k);  % incoming
    end

    % Create node labels
    nodeLabels = arrayfun(@(x) num2str(x), vertImb, 'UniformOutput', false);

    % Plot the graph
    h = plot(G, 'EdgeLabel', G.Edges.Weight, 'NodeLabel', nodeLabels);

    % Improve readability
    layout(h, 'force');                % force-directed layout
    title('Digraph with Vertex Imbalance and Edge Weights');
end
