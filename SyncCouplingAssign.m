% SyncCouplingAssign updates edge weights in a strongly connected digraph
% to enforce negative vertex imbalances and ensure synchronization
% conditions are satisfied.
%
%   Inputs:
%     G – strongly connected digraph
%     a – scalar parameter
%
%   Outputs:
%     G – updated digraph with edge_weight property assigned
%
function G = SyncCouplingAssign(G,a)

    %% Initialize graph properties
    m = G.numedges;     % Number of edges
    n = G.numnodes;     % Number of nodes

    % Initialize edge properties
    G.Edges.edge_weight = zeros(m,1);
    G.Edges.SCC_start   = zeros(m,1);
    G.Edges.SCC_end     = zeros(m,1);
    G.Edges.edge_id     = (1:m).';

    % Initialize node properties
    G.Nodes.node_id   = (1:n).';

    %% Strongly connected components
    bins = conncomp(G,'Type','strong');

    % Assign SCC_start and SCC_end for each edge
    for e = 1:m
        u = G.Edges.EndNodes(e,1);
        v = G.Edges.EndNodes(e,2);
        G.Edges.SCC_start(e) = bins(u);
        G.Edges.SCC_end(e)   = bins(v);
    end

    % Mark nodes that have incoming edges from different SCCs
    has_in_edge = zeros(n,1);
    for e = 1:m
        if G.Edges.SCC_start(e) ~= G.Edges.SCC_end(e)
            v = G.Edges.EndNodes(e,2);
            has_in_edge(v) = 1;
        end
    end
    G.Nodes.has_in_edge = has_in_edge;

    %% Process each strongly connected component (SCC)
    uniqueSCCs = unique(bins);

    for s = uniqueSCCs

        % Extract subgraph for the current SCC
        nodesInSCC = find(bins == s);
        G_SCC      = subgraph(G,nodesInSCC);

        % Root SCC selection:
        %   - If no incoming edges, choose any node
        %   - Otherwise, choose among nodes with incoming inter-SCC edges

        if (sum(has_in_edge(nodesInSCC)) < 1)
            nid = randsample(nodesInSCC,1);
        else
            nid = randsample(find(and(bins == s,has_in_edge.'>0)),1);
        end

        src = find(G_SCC.Nodes.node_id == nid);

        %% Compute negative imbalance vector for this SCC
        G_temp = NegativeImbalanceVectorSCC(G_SCC,src);

        % Path length sum 
        P_sum = max(-incidence(G_temp)*G_temp.Edges.edge_weight);

        % Scale negative imbalance edge weights according 
        % to the theorem in my paper
        G_temp.Edges.edge_weight = 2*a*G_temp.Edges.edge_weight;

        % Update global graph edge weights
        for ei = 1:G_temp.numedges
            eid = G_temp.Edges.edge_id(ei);
            G.Edges.edge_weight(eid) = G.Edges.edge_weight(eid) + G_temp.Edges.edge_weight(ei);
        end

        %% Cycle basis sum
        G_temp = CycleBasisVector(G_SCC);

        % Scale cycle basis vector edge weights

        G_temp.Edges.edge_weight = ...
            (2*a/numel(nodesInSCC))*(1+P_sum)*P_sum*G_temp.Edges.edge_weight;

        % Update global edge weights with cycle basis contributions
        for ei = 1:G_temp.numedges
            eid = G_temp.Edges.edge_id(ei);
            G.Edges.edge_weight(eid) = G.Edges.edge_weight(eid) + G_temp.Edges.edge_weight(ei);
        end

        %% Assign imbalance correction
        outSum = sum(G_temp.Edges.edge_weight(G_temp.Edges.EndNodes(:,1) == src));
        inSum  = sum(G_temp.Edges.edge_weight(G_temp.Edges.EndNodes(:,2) == src));
        diff   = outSum - inSum;

        % Nodes in current SCC with incoming inter-SCC edges
        hi_nodes = nodesInSCC(has_in_edge(nodesInSCC) == 1);

        if ~isempty(hi_nodes)
            % Distribute imbalance equally among incoming inter-SCC edges
            for u = hi_nodes
                incomingE = find(G.Edges.EndNodes(:,2) == u & ...
                                 G.Edges.SCC_start ~= G.Edges.SCC_end);
                if u ~= nid
                    inedge_weight = a;
                else
                    inedge_weight = a + diff/2;
                end

                share = inedge_weight / numel(incomingE);
                for e = incomingE
                    G.Edges.edge_weight(e) = G.Edges.edge_weight(e) + share;
                end
            end
        end
    end
end
