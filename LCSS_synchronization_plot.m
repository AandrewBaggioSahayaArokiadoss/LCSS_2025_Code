% Graph plot with edge weight labels

syms a;

sigma = 10;
rho = 25;
beta = 8/3;

% syms a
a_num = -sigma+(beta*(beta+1)*(rho+sigma)^2)/(16*(beta-1)); %Lorentz Osillator parameter "a" parameter


w11 = [12*a+2 6*a+2 1 2*a+1 1].';
w21 = [7*a a a/2 a/2 a].';

w1 = [w11;w21];

w2 = a*[3+60 1+60 0+30 1+60 0+60 3+30 6+30 0+30].';

w = [w1;w2];

Edge_label = string(w);

w = double(subs(w,a,a_num));

tail = [1 2 3 3 4 5 8 8 7 9 8 6 7 9 10 5 5 7];
head = [2 3 1 4 1 1 2 3 3 4 6 7 8 10 5 9 8 9];

G = digraph(tail,head,w);

N = G.numnodes; % Number of vertices
numStates = 3; % Number of states

data_length = 20;
t_end = 3;
tSpan = linspace(0,t_end,data_length); % Time span

edge_weights = G.Edges.Weight;

EdgeLabel = round(edge_weights,2);

plot(G,'EdgeLabel',EdgeLabel)

mu = 5;
sigma = 10;


% x = -10;
% y = 10;

P = diag([1,0,0]);

X0 = mu + sigma*rand(1,numStates*N);
[X,t] = SimulateCoupledSystems(@LorenzOscillator,tSpan,X0,G,P);

state_index = 1:numStates;
state_index_rem = 1:numStates*N;

colors = lines(N);
linestyle = {'-','--','-.',':','-','--','-.',':','-','--'};
lw = [2*ones(1,4) 1.5*ones(1,4) 1 1];
mark = {'none','none','none','none','*','*','o','o','.','.'};

E = zeros(N,length(t));

figure
grid on
hold on
for i = 1:N
    slice_i = (i-1)*numStates+state_index;
    slice_rem = setdiff(state_index_rem,slice_i);
    e = vecnorm(repmat(X(:,slice_i),1,N-1)-X(:,slice_rem),2,2).';
    E(i,:) = e;
    plot(0:length(t)-1,e,'Color',colors(i,:),'LineWidth',lw(i),'LineStyle',linestyle{i}(:),'Marker',mark{i}(:),'MarkerFaceColor','none','DisplayName', sprintf('System %d',i))
end
xlabel('Time Steps')
ylabel('Sum of square of pairwise distances')
title('Synchronization of 10 Lorentz Oscillators');
legend show
hold off