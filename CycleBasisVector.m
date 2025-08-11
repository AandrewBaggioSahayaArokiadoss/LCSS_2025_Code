% CycleBasisVector  Increment 'null_weight' along random cycles until every edge has positive weight
%   G = CycleBasisVector(G)
%   Input: digraph G
%   Returns: updated digraph G with edge weight G.Edges.null_weight > 0

function G = CycleBasisVector(G)

% Step 1: initialize null_weight to zero on all edges
G.Edges.null_weight = zeros(G.numedges,1);
G.Edges.id = (1:G.numedges).';
% Keep track until all null_weight > 0
disp(G.Edges)

while any(G.Edges.null_weight == 0)
    % Step 2: pick a random edge index among those still zero
    zeroIdx = find(G.Edges.null_weight == 0);
    eidx = zeroIdx(randi(numel(zeroIdx)));
    
    % Step 3: get source and target nodes of that edge
    uv = G.Edges.EndNodes(eidx,:);
    t = uv(1);  % tail/source
    h = uv(2);  % head/target
    
    % Step 4: find a directed path in edges from h to t
    % use shortestpath to get edge indices directly
    [~,~,edgePath] = shortestpath(G, h, t);
    disp(G.Edges.id(edgePath).')
    % Build cycle-edge list C
    if isempty(edgePath)
        % no path, skip this edge this iteration
        continue;
    end
    C = [eidx, edgePath(:)'];
    % Step 5: increment null_weight on all edges in C
    G.Edges.null_weight(C) = G.Edges.null_weight(C) + 1;
end
end
