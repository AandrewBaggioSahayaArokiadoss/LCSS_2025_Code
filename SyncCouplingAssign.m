% SyncCouplingAssign takes in a digraph and a scalar and returns
% the same digraph with its edge weights updated such that the vertex
% imbalances of all the vertices are negative and the synchroniztion
% conditions are satisfied

%   Inputs:
%     G – strongly‑connected digraph
%     a – scalar parameter a
%   Outputs:
%     G   – updated digraph, with the edge variable edge_weight updated

function G = SyncCouplingAssign(G,a)

    % Number of edges
    m = G.numedges;
    % Number of nodes
    n = G.numnodes;

    % Initialize edge and node properties
    G.Edges.edge_weight = zeros(m,1);
    G.Edges.SCC_start   = zeros(m,1);
    G.Edges.SCC_end     = zeros(m,1);
    G.Edges.edge_id     = (1:m).';
    G.Nodes.node_id     = (1:n).';
    G.Nodes.imblance    = zeros(n,1);

    % Compute strongly connected components
    bins = conncomp(G, 'Type', 'strong');

    % Assign SCC_start and SCC_end for each edge
    for e = 1:m
        u = G.Edges.EndNodes(e, 1);
        v = G.Edges.EndNodes(e, 2);
        G.Edges.SCC_start(e) = bins(u);
        G.Edges.SCC_end(e)   = bins(v);
    end

    % Initialize node property has_in_edge
    has_in_edge = zeros(n, 1);

    % Mark nodes that have incoming edges from different SCCs
    for e = 1:m
        if G.Edges.SCC_start(e) ~= G.Edges.SCC_end(e)
            v = G.Edges.EndNodes(e, 2);
            has_in_edge(v) = 1;
        end
    end

    G.Nodes.has_in_edge = has_in_edge;

    % Set the imbalances of nodes with edges from other SCCs to '-a'
    G.Nodes.imbalance(find(has_in_edge==1)) = -a;

    % Number of SCCs
    uniqueSCCs = unique(bins);

    % Process each SCC
    for s = uniqueSCCs

        % Extract the subgraph of the current SCC
        nodesInSCC = find(bins == s);
        G_SCC = subgraph(G,nodesInSCC);

        src = randsample(intersect(nodesInSCC,find(has_in_edge==1)),1);

        % Call NegativeVertexImbalance and update edge weights in G
        G_temp = NegativeImbalanceVectorSCC(G_SCC,src,a);
        % Add G2 edge weights to G for matching edge_id
        for ei = 1:G_temp.numedges
            eid = G_temp.Edges.edge_id(ei);
            G.Edges.edge_weight(eid) = G.Edges.edge_weight(eid) + G_temp.Edges.edge_weight(ei);
        end

        % Call CycleBasisVector and update edge weights in G
        G_temp = CycleBasisVector(G_SCC);
        for ei = 1:numedges(G_temp)
            eid = G_temp.Edges.edge_id(ei);
            G.Edges.edge_weight(eid) = G.Edges.edge_weight(eid) + G_temp.Edges.edge_weight(ei);
        end

        % (ii) Among those nodes, find one with positive (out − in) sum
        hi_nodes = nodesInSCC(has_in_edge(nodesInSCC) == 1);
        bestNode = [];
        if ~isempty(hi_nodes)
            for u = hi_nodes.'
                outSum = sum(G.Edges.edge_weight(G.Edges.EndNodes(:,1) == u & ...
                    G.Edges.SCC_start == s & G.Edges.SCC_end == s));
                inSum  = sum(G.Edges.edge_weight(G.Edges.EndNodes(:,2) == u & ...
                    G.Edges.SCC_start == s & G.Edges.SCC_end == s));
                diff = outSum - inSum;
                if diff > 0
                    bestNode = u;
                    imbalanceIncrement = diff / 2;
                    G.Nodes.imbalance(u) = G.Nodes.imbalance(u) + imbalanceIncrement;
                    break;
                end
            end

        % (iii) Distribute imbalance equally over incoming inter-SCC edges
        for u = hi_nodes.'
            imb = G.Nodes.imbalance(u);
            incomingE = find(G.Edges.EndNodes(:,2) == u & ...
                             G.Edges.SCC_start ~= G.Edges.SCC_end);
            if ~isempty(incomingE)
                share = abs(imb / numel(incomingE));
                for e = incomingE.'
                    G.Edges.edge_weight(e) = G.Edges.edge_weight(e) + share;
                end
            end
        end
        end
    end
end
