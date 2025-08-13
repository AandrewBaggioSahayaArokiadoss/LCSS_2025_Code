function G = SyncCouplingAssign (G, a)
    % processGraph — Process directed graph G with parameter a
    %
    % G : digraph
    % a : scalar parameter
    %
    % For each edge, initialize edge_weight, SCC_start, SCC_end, edge_id.
    % Identify SCCs, mark nodes with incoming inter-SCC edges.
    % For each SCC, adjust imbalance according to provided rules.

    % Number of edges
    m = numedges(G);

    % Initialize edge properties
    G.Edges.edge_weight = zeros(m, 1);
    G.Edges.SCC_start    = zeros(m, 1);
    G.Edges.SCC_end      = zeros(m, 1);
    G.Edges.edge_id      = (1:m).';

    % Compute strongly connected components
    bins = conncomp(G, 'Type', 'strong');  % Using MATLAB's conncomp
    % bins(i) gives SCC index of node i :contentReference[oaicite:1]{index=1}

    % Assign SCC_start and SCC_end for each edge
    for e = 1:m
        u = G.Edges.EndNodes(e, 1);
        v = G.Edges.EndNodes(e, 2);
        G.Edges.SCC_start(e) = bins(u);
        G.Edges.SCC_end(e)   = bins(v);
    end

    % Initialize node property has_in_edge
    n = numnodes(G);
    has_in_edge = zeros(n, 1);

    % Mark nodes that have incoming edges from different SCCs
    for e = 1:m
        if G.Edges.SCC_start(e) ~= G.Edges.SCC_end(e)
            v = G.Edges.EndNodes(e, 2);
            has_in_edge(v) = 1;
        end
    end
    G.Nodes.has_in_edge = has_in_edge;

    % Process each SCC
    uniqueSCCs = unique(bins);
    for s = uniqueSCCs
        % Extract subgraph G1 of this SCC
        nodesInSCC = find(bins == s);
        edgeMask = ismember(G.Edges.SCC_start, s) & ismember(G.Edges.SCC_end, s);
        G1 = subgraph(G, nodesInSCC);
        % Adjust edge_id correspondence in G1
        G1.Edges.edge_id = G.Edges.edge_id(edgeMask);

        % Call NegativeVertexImbalance and update edge weights in G
        G2 = NegativeVertexImbalance(G1, a);  %#ok<NASGU>
        % Add G2 edge weights to G for matching edge_id
        for ei = 1:numedges(G2)
            eid = G2.Edges.edge_id(ei);
            G.Edges.edge_weight(eid) = G.Edges.edge_weight(eid) + G2.Edges.edge_weight(ei);
        end

        % Call CycleBasisVector and update edge weights in G
        G2 = CycleBasisVector(G1, a);  %#ok<NASGU>
        for ei = 1:numedges(G2)
            eid = G2.Edges.edge_id(ei);
            G.Edges.edge_weight(eid) = G.Edges.edge_weight(eid) + G2.Edges.edge_weight(ei);
        end

        % (i) Set imbalance of has_in_edge nodes in SCC to a
        G.Nodes.imbalance(nodesInSCC & has_in_edge == 1) = a;

        % (ii) Among those nodes, find one with positive (out − in) sum
        hi_nodes = nodesInSCC(has_in_edge(nodesInSCC) == 1);
        bestNode = NaN;
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
                share = imb / numel(incomingE);
                for e = incomingE.'
                    G.Edges.edge_weight(e) = G.Edges.edge_weight(e) + share;
                end
            end
        end
    end
end
