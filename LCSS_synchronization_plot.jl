using Plots

# Close all open GR windows
try
    Plots.closeall()
catch
end

ENV["GKSwstype"] = "gksqt"

using Plots
using Graphs
using MetaGraphs
using OrdinaryDiffEq
using CSV
using DataFrames
using Random
using LinearAlgebra
using Statistics

gr()

###############################################################
# OUTPUT DIRECTORY
###############################################################
output_dir = raw"C:\Users\Aandrew Baggio\Julia files"
isdir(output_dir) || mkpath(output_dir)

println("🚀 Random Digraph with Exactly One Root SCC")
println("============================================")

###############################################################
# 1️⃣ Lorenz Oscillator
###############################################################
function LorenzOscillator!(du, u, p, t)
    σ, ρ, β = 10.0, 28.0, 8/3
    du[1] = σ*(u[2] - u[1])
    du[2] = u[1]*(ρ - u[3]) - u[2]
    du[3] = u[1]*u[2] - β*u[3]
end

###############################################################
# 2️⃣ Coupled Network Dynamics (In-degree Laplacian)
###############################################################
function CoupledDynamics!(du, u, p, t)

    systemDynamics!, L, P, a, N, stateDim = p

    U  = reshape(u, stateDim, N)
    dU = reshape(du, stateDim, N)

    temp = zeros(stateDim)

    for i in 1:N
        systemDynamics!(temp, view(U,:,i), nothing, t)
        dU[:, i] .= temp

        for j in 1:N
            Lij = L[i,j]
            if Lij != 0.0
                dU[:,i] -= Lij * (P * view(U,:,j))
            end
        end

        dU[:,i] -= a * (P * view(U,:,i))
    end
end

###############################################################
# 3️⃣ Simulation Wrapper
###############################################################
function SimulateCoupledSystems(systemDynamics!, tspan, X0, G, P, a)

    n = nv(G)
    A = zeros(n,n)

    # In-degree adjacency convention:
    # edge i → j  ⇒  A[j,i]
    for e in edges(G)
        A[dst(e), src(e)] = get_prop(G, e, :weight, 1.0)
    end

    D = Diagonal(vec(sum(A, dims=2)))
    L = D - A

    stateDim = size(P,1)
    params = (systemDynamics!, L, P, a, n, stateDim)

    prob = ODEProblem(
        CoupledDynamics!,
        X0,
        (tspan[1], tspan[end]),
        params
    )

    sol = solve(prob, Tsit5(); saveat=tspan)

    return Array(sol)', sol.t
end

###############################################################
# 4️⃣ Generate Directed Graph with EXACTLY ONE Root SCC
###############################################################
function generateOneRootSCC(n::Int)

    G = MetaDiGraph(n)

    # Root SCC size
    k = max(2, floor(Int, n÷4))
    root = collect(1:k)

    # Make root strongly connected (cycle)
    for i in 1:k
        add_edge!(G, root[i], root[mod1(i+1,k)])
    end

    # Add extra internal edges
    for i in root, j in root
        if i != j && rand() < 0.3
            add_edge!(G, i, j)
        end
    end

    # Remaining vertices
    for v in k+1:n
        parent = rand(root)  # incoming edge from root SCC
        add_edge!(G, parent, v)
    end

    return G
end

###############################################################
# 5️⃣ Assign Coupling Weights
###############################################################
function SyncCouplingAssign(G, weight)
    for e in edges(G)
        set_prop!(G, e, :weight, weight)
    end
    return G
end

println("✓ All functions ready")

###############################################################
# ===================== MAIN EXECUTION =======================
###############################################################

N = 100

σ, ρ, β = 10.0, 25.0, 8/3
a = -σ + (β*(β+1)*(ρ+σ)^2)/(16*(β-1))

println("🔧 Coupling a = ", round(a,digits=4))

G = generateOneRootSCC(N)
G = SyncCouplingAssign(G, 1.01a)

###############################################################
# FIGURE 1 — Directed Graph with Arrows
###############################################################

θ = range(0,2π,length=N+1)[1:end-1]
x = cos.(θ)
y = sin.(θ)

p1 = scatter(x,y,
             markersize=6,
             markercolor=:lightblue,
             legend=false,
             axis=false,
             size=(600,600))

for e in edges(G)
    plot!(p1,
          [x[src(e)], x[dst(e)]],
          [y[src(e)], y[dst(e)]],
          arrow=:arrow,
          lw=1.2,
          c=:black)
end

display(p1)  # Separate window 1
savefig(p1, joinpath(output_dir,"network.svg"))

println("📊 Saved: network.svg")

###############################################################
# FIGURE 2 — Pairwise Distances Between Consecutive Nodes
###############################################################

P = Diagonal([1.0,0.0,0.0])
tspan = range(0,100,length=10)
X0 = 1e-6 .* randn(3N)

println("⏱️ Simulating...")

X,t = SimulateCoupledSystems(
    LorenzOscillator!,
    collect(tspan),
    X0,
    G,
    P,
    a
)

println("📏 Computing pairwise distances...")

num_time = size(X,1)
pairwise_distances = zeros(num_time, N-1)

for k in 1:num_time
    U = reshape(X[k,:], 3, N)  # reshape row into (stateDim × N)
    
    for i in 1:N-1
        pairwise_distances[k,i] = norm(U[:,i] - U[:,i+1])
    end
end

###############################################################
# Plot all pairwise distances on same graph (no legend)
###############################################################

p2 = plot(size=(900,500),
          grid=true,
          legend=false,
          xlabel="Time",
          ylabel="Pairwise Distances",
          title="Total Pairwise Distance vs Time")

for i in 1:N-1
    plot!(p2, t, pairwise_distances[:,i], lw=1.5)
end

display(p2)
savefig(p2, joinpath(output_dir,"pairwise_distances.svg"))

println("📈 Saved: pairwise_distances.svg")

###############################################################
# Compute average of final-time pairwise distances
###############################################################

final_distances = pairwise_distances[end, :]
avg_final_distance = mean(final_distances)

println("📊 Average final pairwise distance = ",
        round(avg_final_distance, digits=6))

println("✅ COMPLETE")
