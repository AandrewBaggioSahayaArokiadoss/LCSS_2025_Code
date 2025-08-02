% Graph plot with edge weight labels

syms a;

w11 = [12*a+2 6*a+2 1 2*a+1 1].';
w21 = [7*a a a/2 a/2 a].';

w1 = [w11;w21];

w2 = a*[6+30 3+30 1+60 0+30 0+30 3+60 3+30 1+60].';

w = [w2;w1];

Edge_label = string(w);

disp(Edge_label)

G = digraph([1 2 3 4 3 5 8 8 7 9 8 6 7 9 10 5 5 7],[2 3 4 1 1 1 2 3 3 4 6 7 8 10 5 9 8 9]);

G.Edges.edge_id = (1:G.numedges).';
G.Edges.label = string(zeros(G.numedges,1));

for i = 1:G.numedges
    idx = G.Edges.edge_id(i);
    G.Edges.label(i)=Edge_label(i);
end

G.Nodes.imbalance = string(zeros(G.numnodes,1));

plot(G,"EdgeLabel",G.Edges.label)