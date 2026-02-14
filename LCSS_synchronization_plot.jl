###############################################################
#  FORCE EXTERNAL WINDOWS (IMPORTANT)
###############################################################
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

N = 50

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
savefig(p1, joinpath(output_dir,"network.pdf"))

println("📊 Saved: network.pdf")

###############################################################
# FIGURE 2 — Synchronization Error
###############################################################

P = Diagonal([1.0,0.0,0.0])
tspan = range(0,100,length=50)
X0 = 1e-4 .* randn(3N)

println("⏱️ Simulating...")

X,t = SimulateCoupledSystems(
    LorenzOscillator!,
    collect(tspan),
    X0,
    G,
    P,
    a
)

sync_error = [var(X[i,1:3:3N]) for i in 1:size(X,1)]

p2 = plot(t,sync_error,
          lw=3,
          c=:red,
          yscale=:log10,
          xlabel="Time",
          ylabel="Var(x)",
          legend=false,
          size=(900,400))

display(p2)  # Separate window 2
savefig(p2, joinpath(output_dir,"sync.pdf"))

println("📈 Saved: sync.pdf")
println("✅ COMPLETE")