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
        u = G.Edges.EndNodes(e,1);
        v = G.Edges.EndNodes(e,2);
        G.Edges.SCC_start(e) = bins(u);
        G.Edges.SCC_end(e)   = bins(v);
    end

    % Initialize node property has_in_edge
    has_in_edge = zeros(n, 1);

    % Mark nodes that have incoming edges from different SCCs
    for e = 1:m
        if G.Edges.SCC_start(e) ~= G.Edges.SCC_end(e)
            v = G.Edges.EndNodes(e,2);
            has_in_edge(v) = 1;
        end
    end

    G.Nodes.has_in_edge = has_in_edge;

    % Number of SCCs
    uniqueSCCs = unique(bins);

    % Process each SCC
 for s = uniqueSCCs

        % Extract the subgraph of the current SCC
        nodesInSCC = find(bins == s);
        G_SCC = subgraph(G,nodesInSCC);
        
        % In case G_SCC is a root SCC, we can select
        % src arbitrarily otherwise, it is selected
        % from among the nodes with incoming
        % edges from other SCCs
     
        if(sum(has_in_edge(nodesInSCC))<1)
            src = randi([1,numel(nodesInSCC)]);
            nid = G_SCC.Nodes.node_id(src);
        else
            nid = randsample(intersect(nodesInSCC,find(has_in_edge==1)),1);
            src = find(G_SCC.Nodes.node_id==nid);
        end

        disp(nid)

        % Call NegativeVertexImbalanceSCC with G_SCC
        G_temp = NegativeImbalanceVectorSCC(G_SCC,src);

        % Find the sum of path lengths from src in G_SCC 
        % by finding the largest value of vector imbalance
        P_sum = max(-incidence(G_temp)*G_temp.Edges.edge_weight);

        % Scaling the negative vertex imbalance vector edge weights
        % by 2a
        G_temp.Edges.edge_weight = 2*a*G_temp.Edges.edge_weight;

        % Add G_temp edge weights to G for matching edge_id
        for ei = 1:G_temp.numedges
            eid = G_temp.Edges.edge_id(ei);
            G.Edges.edge_weight(eid) = G.Edges.edge_weight(eid) + G_temp.Edges.edge_weight(ei);
        end

        % Call CycleBasisVector and update edge weights in G
        G_temp = CycleBasisVector(G_SCC);

        % Scaling the cycle basis vector edge weights to meet
        % the synchronization condition

        G_temp.Edges.edge_weight = (2*a/numel(nodesInSCC))*(1+P_sum)*P_sum*G_temp.Edges.edge_weight;

        % Adding this to the edge weights of G by comparing edge ids
        for ei = 1:numedges(G_temp)
            eid = G_temp.Edges.edge_id(ei);
            G.Edges.edge_weight(eid) = G.Edges.edge_weight(eid) + G_temp.Edges.edge_weight(ei);
        end

        figure
        VertexImbalancePlot(G)

        % Compute the vertex imbalance of src and assign it to the
        % imbalance variable of src node

        outSum = sum(G.Edges.edge_weight(G.Edges.EndNodes(:,1) == nid & ...
                    G.Edges.SCC_start == s & G.Edges.SCC_end == s));
        inSum  = sum(G.Edges.edge_weight(G.Edges.EndNodes(:,2) == nid & ...
                    G.Edges.SCC_start == s & G.Edges.SCC_end == s));
        diff = outSum - inSum;

        % Find the nodes among the nodes of the current SCC with
        % incoming edges from other SCCs
        hi_nodes = nodesInSCC(has_in_edge(nodesInSCC) == 1);

        if (~isempty(hi_nodes))
            % Distribute imbalance equally over incoming inter-SCC edges
            for u = hi_nodes
                incomingE = find(G.Edges.EndNodes(:,2) == u & ...
                                G.Edges.SCC_start ~= G.Edges.SCC_end);
                if(u==nid)
                    inedge_weight = a + (diff/2);
                else
                    inedge_weight = a;
                end
                if ~isempty(incomingE)
                    share = inedge_weight/numel(incomingE);
                    for e = incomingE.'
                        G.Edges.edge_weight(e) = G.Edges.edge_weight(e) + share;
                    end
                end
            end
        end
        figure
        VertexImbalancePlot(G)
end
