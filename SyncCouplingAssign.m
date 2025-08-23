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

        % Find the sum of all path lengths from src in G_SCC 
        % by finding the largest value of vector imbalance
        P_sum = max(-incidence(G_temp)*G_temp.edge_weight);

        % Scaling the negative vertex imbalance vector edge weights
        % by 2a
        G_temp.edge_weight = 2*a*G_temp.edge_weight;

        % Call CycleBasisVector and update edge weights in G
        G_temp = CycleBasisVector(G_SCC);

        % Scaling the cycle basis vector edge weights to meet
        % the synchronization condition

        G_temp.edge_weight = (2*a/nodesInSCC)*(1+P_sum)*P_sum*G_temp.edge_weight;

        % Adding this to the edge weights of G by comparing edge ids
        for ei = 1:numedges(G_temp)
            eid = G_temp.Edges.edge_id(ei);
            G.Edges.edge_weight(eid) = G.Edges.edge_weight(eid) + G_temp.Edges.edge_weight(ei);
        end

        % Compute the vertex imbalance of src and assign it to the
        % imbalance variable of src node

        outSum = sum(G.Edges.edge_weight(G.Edges.EndNodes(:,1) == src & ...
                    G.Edges.SCC_start == s & G.Edges.SCC_end == s));
        inSum  = sum(G.Edges.edge_weight(G.Edges.EndNodes(:,2) == src & ...
                    G.Edges.SCC_start == s & G.Edges.SCC_end == s));
        diff = outSum - inSum;
        imbalanceIncrement = diff / 2;
        G.Nodes.imbalance(src) = a + imbalanceIncrement;

        % Find the nodes among the nodes of the current SCC with
        % incoming edges from other SCCs
        hi_nodes = nodesInSCC(has_in_edge(nodesInSCC) == 1);

        % Distribute imbalance equally over incoming inter-SCC edges
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
