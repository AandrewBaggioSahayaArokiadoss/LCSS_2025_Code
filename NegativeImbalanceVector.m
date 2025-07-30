% NegativeImbalanceVector - computes a vector with positive entries that
% can make all but one vertex imbalance to have a negative values

function G = NegativeImbalanceVector(G, a)

    m = numedges(G);                   % total number of edges
    G.Edges.edge_id = (1:m)';          % assign unique ID 1..m
    
    % find strongly connected components
    comp = conncomp(G, 'Type', 'strong');
    numComp = max(comp);
    compSizes = histcounts(comp, 1:numComp+1);
    
    if numComp == 1
        % Case 1: strongly connected
        src = randi(numnodes(G));      % random node index
        G = NegativeImbalanceVectorSCC(G, src, a);
        
    elseif numComp == 2 && any(compSizes == 1)
        % Case 2: exactly two SCCs, one of size 1
        % identify isolated singleton component
        singletonComp = find(compSizes == 1);
        v_r = find(comp == singletonComp);  % node index of that singleton
        % verify no incoming edges
        if indegree(G, v_r) > 0
            error('Expected the singleton to have no incoming edges.');
        end
        
        % subgraph G1 by removing v_r
        nodesToKeep = setdiff(1:numnodes(G), v_r);
        G1 = subgraph(G, nodesToKeep);
        
        % find an edge from v_r: find successors
        succ = successors(G, v_r);
        if isempty(succ)
            error('Expected v_r to have outgoing edges.');
        end
        src = succ(randi(numel(succ)));
        
        % map src into index in G1
        % In subgraph, node indices are preserved order so find index:
        srcInG1 = find(nodesToKeep == src);
        G1 = NegativeImbalanceVectorSCC(G1, srcInG1, a);
        
        % initialize weights in G
        G.Edges.Weight = zeros(m,1);
        % copy weights by matching edge_id
        % First, ensure G1 has edge_id and Weight
        % match edge‑ids:
        [tf, loc] = ismember(G1.Edges.edge_id, G.Edges.edge_id);
        G.Edges.Weight(loc(tf)) = G1.Edges.Weight(tf);
        
        % now handle edges out of v_r
        outEdges = outedges(G, v_r);
        for k = 1:numel(outEdges)
            eIdx = outEdges(k);
            e = G.Edges(eIdx, :);
            tgt = e.EndNodes(2);
            if tgt == src
                imb = indegree(G1, srcInG1) - outdegree(G1, srcInG1);
                G.Edges.Weight(eIdx) = a + 0.5 * imb;
            else
                G.Edges.Weight(eIdx) = a;
            end
        end
        
    else
        % Case 3: failure structure
        disp('The digraph is not in a proper structure');
        G = digraph();  % return empty digraph
    end
end
