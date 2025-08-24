% CycleBasisVector takes in a digraph returns the same digraph with its edge weights updated
% such that the vertex imbalances are unaltered but the edge weights are non zero
%   Inputs:
%     G   – strongly‑connected digraph
%   Outputs:
%     G   – updated digraph, with the edge variable edge_weight updated

function G = CycleBasisVector(G)

m = G.numedges;
null_weight = zeros(m,1);

% Keep track until all null_weight > 0
while any(null_weight == 0)
    % Step 1: pick a random edge index among those still zero
    zeroIdx = find(null_weight == 0);
    eidx = zeroIdx(randi(numel(zeroIdx)));
    
    % Step 2: get source and target nodes of that edge
    uv = G.Edges.EndNodes(eidx,:);
    t = uv(1);  % tail/source
    h = uv(2);  % head/target
    
    % Step 3: find a directed path in edges from h to t
    % use shortestpath to get edge indices directly
    [~,~,edgePath] = shortestpath(G, h, t);

    % Build cycle-edge list C
    if isempty(edgePath)
        % no path, skip this edge this iteration
        continue;
    else
        C = [eidx, edgePath(:)'];
        % Step 4: increment null_weight on all edges in C
        null_weight(C) = null_weight(C)+1;
    end  
end
% Update the edge weights of G
G.Edges.edge_weight = G.Edges.edge_weight + null_weight;
end
