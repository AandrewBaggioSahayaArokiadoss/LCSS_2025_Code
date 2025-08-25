% VertexImbalancePlot computes the vertex imbalance vector of a digraph
% and plots the graph with edge weights and vertex imbalances.
%
%   Inputs:
%     G – strongly connected digraph
%
%   Outputs:
%     vertImb – vertex imbalance vector
%
function vertImb = VertexImbalancePlot(G)

    %% Validate input
    assert(isa(G,'digraph'), 'Input must be a digraph object.');

    %% Setup
    n = numnodes(G);                         % number of nodes
    s = G.Edges.EndNodes(:,1);               % source nodes of edges
    t = G.Edges.EndNodes(:,2);               % target nodes of edges
    w = G.Edges.edge_weight;                 % edge weights

    %% Compute vertex imbalance
    vertImb = zeros(n,1);
    for k = 1:numel(w)
        vertImb(s(k)) = vertImb(s(k)) + w(k);   % outgoing weight
        vertImb(t(k)) = vertImb(t(k)) - w(k);   % incoming weight
    end

    %% Plot graph with labels
    % Node labels: imbalance values
    nodeLabels = arrayfun(@num2str, vertImb, 'UniformOutput', false);

    % Edge labels: rounded weights
    % h = plot(G, ...
    %     'EdgeLabel', round(G.Edges.edge_weight,3), ...
    %     'NodeLabel', G.Nodes.node_id);

    h = plot(G,'NodeLabel',strcat(string(G.Nodes.node_id),'(',nodeLabels,')'));

    % Layout and styling
    layout(h,'force');   % force-directed layout
    title('Digraph with Vertex Imbalance and Edge Weights');
end
