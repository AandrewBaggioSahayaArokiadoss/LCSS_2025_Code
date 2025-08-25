% CycleBasisVector updates edge weights of a strongly connected digraph
% so that all edges have nonzero weight while preserving vertex imbalances.
%
%   Inputs:
%     G – strongly connected digraph
%
%   Outputs:
%     G – updated digraph with edge_weight property modified
%
function G = CycleBasisVector(G)

    %% Initialize temporary weights
    null_weight = zeros(G.numedges,1);

    % Continue until all edges have been assigned nonzero weights
    while any(null_weight == 0)

        % --- Step 1: Pick a random edge among those still unweighted ---
        zeroIdx = find(null_weight == 0);
        eidx    = zeroIdx(randi(numel(zeroIdx)));

        % --- Step 2: Get source (tail) and target (head) nodes of edge ---
        uv = G.Edges.EndNodes(eidx,:);
        t  = uv(1);   % tail/source
        h  = uv(2);   % head/target

        % --- Step 3: Find a directed path from head back to tail ---
        % shortestpath can return the sequence of edge indices directly
        [~,~,edgePath] = shortestpath(G,h,t);

        % --- Step 4: Build cycle and update weights ---
        if ~isempty(edgePath)
            C = [eidx, edgePath(:)'];            % cycle edges
            null_weight(C) = null_weight(C) + 1; % increment weights
        end
    end

    %% Update edge weights of G
    G.Edges.edge_weight = G.Edges.edge_weight + null_weight;
end
