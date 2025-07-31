\documentclass[journal,twoside,web]{ieeecolor}
\usepackage{lcsys}
\usepackage[labelsep=period,skip=2pt]{caption}
\captionsetup[table]{name=table}
\usepackage{booktabs}
\usepackage{xcolor}
\usepackage[noadjust]{cite}
\usepackage{amsmath,amsfonts}
\usepackage{mathtools}
\usepackage{mathrsfs}
\usepackage{relsize}
\usepackage{bbm}
\usepackage{algorithm}
\usepackage{algorithmic}
\usepackage{graphicx}
\usepackage{caption}
%\usepackage{subfigure}
\usepackage{subcaption}
\usepackage{comment}
\usepackage{csquotes}
\usepackage{textcomp}
\usepackage{tikz}
\usetikzlibrary{positioning, arrows, calc}
\usepackage{ragged2e}
\usepackage{hyperref}
\usepackage{cleveref}

\hypersetup{linkbordercolor=red}

\usetikzlibrary{decorations.markings,decorations.pathreplacing,positioning}
\tikzset{every node/.style={circle}, 
	strike through/.append style={
		decoration={markings, mark=at position 0.5 with {
				\draw[-] ++ (-2pt,-2pt) -- (2pt,2pt);}
		},postaction={decorate}}
}

\markboth{\journalname}{}
\newtheorem{theorem}{Theorem}
\newtheorem{lemma}{Lemma} 
\newtheorem{problem}{Problem} 
\newtheorem{corollary}{Corollary} 
\newtheorem{assumption}{Assumption}
\newtheorem{definition}{Definition} 
\newtheorem{proposition}{Proposition} 
\newtheorem{remark}{Remark}

\newcommand{\G}{\boldsymbol{\mathrm{G}}}
\newcommand{\HH}{\boldsymbol{\mathrm{H}}}
\newcommand{\PP}{\boldsymbol{\mathcal{P}}}
\newcommand{\CC}{\boldsymbol{\mathcal{C}}}
\DeclareMathOperator*{\EE}{\mathtt{E}}
\DeclareMathOperator*{\VV}{\mathtt{V}}
\newcommand{\Break}{\State \textbf{break} }

\DeclareMathOperator*{\argmin}{\text{argmin }}
\DeclareMathOperator*{\argmax}{\text{argmax }}
\DeclareMathOperator*{\D}{\mathcal{D}}
\DeclareMathOperator*{\E}{\mathcal{E}}


\def\BibTeX{{\rm B\kern-.05em{\sc i\kern-.025em b}\kern-.08em
    T\kern-.1667em\lower.7ex\hbox{E}\kern-.125emX}}
\markboth{\journalname}{}
\begin{document}
\title{A coupling strength allocation algorithm for the synchronization of dynamical networks using cycle basis}
\author{Aandrew Baggio Sahaya Arokiadoss and G. Arunkumar
\thanks{This paper was submitted on 23 May, 2025}
\thanks{Aandrew Baggio Sahaya Arokiadoss is now with the Department of Electrical Engineering, Indian Institute of Technology Madras, Chennai, Tamilnadu, India (e-mail: ee20d067@smail.iitm.ac.in/s.aandrewbaggio@protonmail.com).}
\thanks{G. Arunkumar is with the Department of Mathematics, Indian Institute of Technology Madras, Chennai, India (e-mail: garunkumar@iitm.ac.in).}
}

\maketitle

\begin{abstract}
We present an efficient method for computing coupling strengths that guarantee synchronization in time‐varying dynamical networks. Building on prior work that reformulated a Matrix Inequality into a system of linear inequalities via spectral graph theory, our approach addresses two major drawbacks: high computational cost and the absence of guaranteed solutions. We prove that a solution always exists for any digraph containing a directed spanning tree. By exploiting the structure of strongly connected components, we reduce the path count from $\binom{n}{2}$ to just $n-1$. Rather than relying on generic inequality solvers, we use the cycle basis of each component to satisfy the synchronization condition. The result is a scalable, computationally efficient alternative to the traditional lyapunov‐based method that entirely depends upon graph theoretic tools.
\end{abstract}


\begin{IEEEkeywords}
Coupling strength allocation, Cycle basis, Spectral graph theory, Synchronization, Time varying dynamical networks
\end{IEEEkeywords}

\section{Introduction}
\label{sec:introduction}
\IEEEPARstart{S}{ynchronization} in dynamical networks is a fundamental phenomenon studied extensively across disciplines such as neuroscience \cite{breakspear2010generative}, physics \cite{levis2017synchronization}, biology \cite{karakaya2022effective} engineering \cite{dorfler2013synchronization} and various other fields. Achieving synchronization typically involves understanding how the interplay between network topology and node dynamics governs collective behavior \cite{pecora1998master}.

A commonly used technique for analyzing and ensuring synchronization in dynamical networks is the Lyapunov function method. When applied in this context, it frequently yields conditions that can be formulated as a Matrix Inequality (MI) \cite[Theorem 4.4]{wu2007synchronization}. Solving these LMIs provides remark into the required coupling strength between the dynamical systems in the network to achieve synchronization. However, as network size grows, solving such an MI to identify suitable parameters for synchronization becomes increasingly computationally intensive. This issue is especially critical in time-varying networks, where each change in connectivity requires the solution of a new MI. One notable contribution that attempted to simplify this process is the work in \cite{liu2015synchronization}, where the authors proposed a method called the \emph{Generalized Connection Graph Method} (GCGM). This algorithm determines coupling strengths that satisfy a prescribed set of inequalities to guarantee synchronization. It requires computing a set of undirected paths between each pair of vertices of the network’s connectivity digraph and solving the set of the resulting inequalities. This process must be repeated each time the network connections change.

Our approach for achieving network synchronization is based on the GCGM algorithm. It involves extracting two vectors for each strongly connected component (SCC) of the connectivity graph of the dynamical network. The two vectors are derived exclusively through graph-theoretic methods. The scaled sum of these vectors provides the coupling strengths used in the synchronization process. We prove that the star graph can be used for computing the arc weights of the arcs in-between SCCs. We demonstrate that the resulting coupling strengths inherently fulfill the necessary synchronization inequalities.
\section{Preliminaries}
\subsection{Graphs and Digraphs}
A \emph{digraph} (directed graph) is defined as an ordered pair $(\VV,\EE)$. $\VV = \{v_{1}, \dots, v_{n}\}$ is called the set of vertices and $\EE = \{e_{1}, \dots, e_{m}\}$ ($\EE \subseteq \VV \times \VV$) is the set of arcs of the digraph. An \emph{undirected graph} or (simply a graph) is defined as an ordered pair $(\VV, \EE)$, where $\VV = \{v_{1}, \dots, v_{n}\}$ is its set of vertices and $\EE = \{e_{1}, \dots, e_{m}\}$ is its multiset of edges. The frequently used definitions regarding these two mathematical objects have been tabulated in the table \ref{tab:Graph definitions}. In this work, we restrict our attention to simple digraphs i.e. digraphs devoid of self-loops. Note that the graphs are allowed parallel edges. Unless stated otherwise, all digraphs considered henceforth are assumed to be simple. 

\begin{table}[h]
\centering
\caption{Functions, Terms, and Definitions}
\renewcommand{\arraystretch}{1.3}
\begin{tabular}{|p{2.5cm}|p{5cm}|}
\hline
\textbf{Term / Function} & \textbf{Definition / Description} \\
\hline
$e_{k} = (v_{i}, v_{j})$ & Arc from vertex $v_{i}$ to vertex $v_{j}$ \\
\hline
$e_{k} = \{v_{i}, v_{j}\}$ & Edge between $v_{i}$ and $v_{j}$ \\
\hline
End vertices & Vertices that make up an arc (edge)\\
\hline
Incident arc (edge) of a vertex & Arcs (edges) that contain the vertex \\
\hline
Tail & First vertex of an arc (ordered pair) \\
\hline
Head & Second vertex of an arc (ordered pair) \\
\hline
Incoming (Outgoing) arc & An arc is incoming to (outgoing from) a vertex if the vertex is its head (tail)\\
\hline
Indegree & Number of arcs with the vertex as head \\
\hline
Outdegree & Number of arcs with the vertex as tail \\
\hline
Adjacent arcs (edges) & Arcs (edges) sharing at least one end vertex \\
\hline
Source vertex & Vertex with 0 indegree\\
\hline
Parallel arcs (edges) & Arcs (edges) sharing the same end vertices \\
\hline
Underlying graph & Graph formed by replacing each arc in a digraph with an edge between the same vertices\\
\hline
Self-loop & Arc $e_{i}$ such that $e_{i} = (v_{k}, v_{k})$ \\
\hline
Traversal & A sequence of vertices with each consecutive pair is connected by an arc\\
\hline
$\VV(\G)$ & Vertex set of digraph  (graph) $\G$ \\
\hline
$\EE(\G)$ & Arc set  (edge multiset) of digraph (graph) $\G$ \\
\hline
$|\G|$ & Number of vertices in digraph $\G$ \\
\hline
$||\G||$ & Number of arcs in digraph $\G$ \\
\hline
$\mathcal{H}$ & A function that maps an arc to its head\\
\hline
$\mathcal{T}$ & A function that maps an arc to its tail\\
\hline
$\E$ & Vector of arc weights \\
\hline
$\varepsilon_{i}$ & Weight assigned to arc $e_{i}$ \\
\hline
\end{tabular}
\label{tab:Graph definitions}
\end{table}


%Arcs are ordered pairs of vertices i.e. $E \subseteq V \times V$. If $e_{k}$ is an arc i.e. $e_{k}=(v_{i}, v_{j})$, it is said to originate from $v_{i}$ and terminate at $v_{j}$. Each edge $e_{k}\in \EE$ is a two element subset of $\VV$ i.e. $e_{k}=\{v_{i}, v_{j}\} \subseteq \VV$. The vertices that make up an arc (edge) are called its \emph{end vertices}. An arc (edge) is said to be incident on a vertex if the vertex is one of its end vertices. The number of edges incident on a vertex is called its \textit{degree}. The first vertex of the ordered pair of vertices of an arc is its tail vertex and the other one, its head vertex. The \emph{in-degree} (\emph{out-degree}) of a vertex $v$ is the number of arcs that have $v$ as their tail (head). An arc is said to leave its tail vertex and enter its head vertex a vertex. Two arcs (edges) are said to be adjacent if they share at least one end vertex. Arcs (edges) sharing the same end vertices are called \emph{parallel arcs (edges)}. We refer to the \emph{underlying graph} of a digraph as the undirected graph obtained by replacing each arc with an edge connecting the same pair of vertices.\\For a digraph $\G$, we use $\VV(\G)$ and $\EE(\G)$ to represent its vertex set and arc set respectively. The cardinality of the vertex set and arc set are denoted by $|\G|$ and $||\G||$ respectively. For a digraph $\G$, we define two mappings : $\mathcal{H}: \EE(\G) \rightarrow \VV(\G)$ and $\mathcal{T}: \EE(\G) \rightarrow \VV(\G)$, where $\mathcal{H}(e_{k})$ and $\mathcal{T}(e_{k})$ denote the head (terminating vertex) and tail (starting vertex) of the arc $e_{k}$ respectively.\\An arc $e_{i}$ is referred to as a self-loop if its end vertices are identical i.e. $e_{i}=$($v_{k},v_{k}$). \\
A weighted digraph is a digraph in which each arc has an associated weight. It is represented as a triplet ($\VV, \EE, \E$), where $\E$ is the vector of arc weights, given by $\E = \left[\varepsilon_{1}, \dots, \varepsilon_{m}\right]^{\top}$.
\begin{definition}[Vertex imbalance]
The vertex imbalance of a vertex $v_{i}$ of a weighted digraph is defined as the difference between the sum of theweights of the outgoing arcs and the sum of the weights of the incoming arcs of $v_{i}$ . Formally, it is expressed as:
\begin{equation*}
\text{D}_{i}=\smashoperator{\sum_{\{s :\text{$\mathcal{T}$}(e_{s})=v_{i}\}}}\varepsilon_{s}\;\;\;-\;\;\;\;\smashoperator{\sum_{\{ t :\text{$\mathcal{H}$}(e_{t})=v_{i}\}}}\varepsilon_{t}
\end{equation*}    
\end{definition}
\subsection{Subgraphs}
A \emph{subgraph} $\HH$ of a digraph (graph) $\G$ is a digraph (graph) such that $\VV(\HH) \subseteq \VV(\G)$ and $\EE(\HH) \subseteq \EE(\G)$ and is denoted as $\HH \subseteq \G$. A \emph{path} $\PP$ is a graph such that $\VV(\PP) = \{v_{k_{0}}, \dots, v_{k_{\ell}}\}$ and $\EE(\PP) = \{e_{k_1}, \dots, e_{k_{\ell}}\}$ such that $e_{k_i} = \{v_{k_{i-1}}, v_{k_{i}}\}$. graph is \emph{connected} if, for every pair of distinct vertices, there exists a path subgraph that includes them both. A \emph{cycle} is a connected graph in which every vertex has degree 2. A \emph{tree} is a connected graph wherein none of its subgraphs is a cycle. A \emph{spanning tree} is a subgraph of a graph that is a tree and contains all its vertices. A \emph{directed path} $\PP$ is a digraph with $\VV(\PP)=\{v_{k_{0}}, \dots, v_{k_{\ell}}\}$ and $\EE(\PP)=\{e_{k_{1}}, \dots, e_{k_{\ell}}\}$  such that each arc $e_{k_{i}}=(v_{k_{i-1}},v_{k_{i}})$ for all $1 \leq i \leq  \ell $. The tail and head of $\PP$ are defined to be $\mathcal{T}(\PP) = \mathcal{T}(e_{k_{1}})$ and $\mathcal{H}(\PP) = \mathcal{H}(e_{k_{\ell}})$, respectively. As a subgraph, a directed path induces a natural ordering on its edges where the edge incident on the source vertex is the first edge, its adjacent edge is the second edge and so on.  A digraph is connected if its underlying graph is connected. A \emph{directed cycle} is a digraph whose underlying undirected graph forms a cycle. A subgraph of a digraph is considered its \emph{directed spanning tree} if it has exactly one source vertex and its underlying graph is a spanning tree. A directed spanning tree is said to be rooted at the source vertex.\\
The length of a directed path/directed cycle (path/cycle) is given by the number of arcs (edges) they contain. 
%A digraph is said to contain a directed cycle/path/tree if it has a subgraph that is a directed cycle/path/tree and it is the same for their directed counterparts.
A digraph is strongly connected if there exists a directed path between any two ordered pair of distinct vertices.
\subsection{Matrices and Vectors}
Let $\text{X} \in \mathbb{R}^{n}$ be a vector. We use the notation $\text{X} > 0$ to indicate that all the coordinates of $\text{X}$ are strictly positive. The $n \times 1$ zero vector is denoted by $\mathbf{0}_{n}$, while the vector of all ones is denoted by $\mathbf{1}_{n}$. We define the adjacency (A), laplacian (L) and incidence matrix (Q) elementwise for a weighted digraph $\G$ with $n$ vertices and $m$ arcs as :
\begin{equation*}
\text{A}_{ij} =
\begin{cases}
\varepsilon_{k}, & \text{if } e_{k}=(j,i) \in \EE(\G), \\
0, & \text{otherwise},
\end{cases}
\end{equation*}

Let D$_{ii}= \sum_{j=1}^{n}$A$_{ij}$.

\begin{equation*}
    \text{L}_{ij}=\begin{cases}\text{D}_{ii} \text{, if }i=j\\
        -\text{A}_{ij}\text{, otherwise}
    \end{cases}
\end{equation*}
In case of unweighted digraphs/graphs, all the weights assigned are assumed to be equal to $1$. We define the \emph{incidence matrix} of $\G$ as follows :
	\begin{equation*}
		\begin{aligned}
			\text{Q}_{ij}=\begin{cases}
				\phantom{-}1,\text{ if } v_{i}= \text{$\mathcal{T}$(}e_{j}\text{)}\\
				-1,\text{ if } v_{i}= \text{$\mathcal{H}$(}e_{j}\text{)}\\
				\phantom{-}0,\text{ otherwise}
			\end{cases}
		\end{aligned}
	\end{equation*}
The vertex imbalances of a digraph $\G$ relate to the arc weights through the equation Q$\E = \mathcal{D}$, where $\mathcal{D} \in \mathbb{R}^n$ is the \emph{vertex imbalance vector}, and its $i^{\text{th}}$ entry represents the imbalance at vertex $v_i \in \VV(\G)$. Each directed cycle $\HH$ in $\G$ can be associated with a \emph{signed incidence vector} $\text{X} = \left[x_{1}, \dots, x_{m}\right]^{\top} \in \mathbb{R}^{||\G||}$, defined component-wise as:
\begin{equation*}
x_{i} =
\begin{cases}
\phantom{-}1, & \text{if } e_{i} \in \EE(\HH) \text{ traversed from tail to head}, \\
-1, & \text{if } e_{i} \in \EE(\HH) \text{ traversed from head to tail}, \\
\phantom{-}0, & \text{otherwise}.
\end{cases}
\end{equation*}

The sign of $\text{X}$ depends on traversal direction, but this ambiguity is inconsequential in our context, as we consider cycles from directed ear decompositions. Traversal becomes relevant when arc weight vectors with positive entries are required.

The \emph{cycle space} $\mathscr{C}$ of $\G$ is the nullspace of its incidence matrix $Q$, i.e., $\mathscr{C} = \{\text{X} \in \mathbb{R}^m \mid \text{QX} = \boldsymbol{0}_n\}$. A \emph{cycle basis} of $\G$ is a set of linearly independent signed incidence vectors corresponding to distinct directed cycles that span $\mathscr{C}$.

\section{Problem Formulation}
We analyze a dynamical network framework with $n$ dynamical systems, similar to the one presented in \cite{liu2015synchronization}:
\begin{equation}
\dot{\boldsymbol{z}}_{i}=f(\boldsymbol{z}_{i})+\smashoperator{\sum}\varepsilon_{k}(t)\text{P}(\boldsymbol{z}_{j}-\boldsymbol{z}_{i})
\label{eq:Dynamical Network Model}
\end{equation}
where the sum is over the set $\{j:\text{$\mathcal{T}$}(e_{k})=v_{j},\text{$\mathcal{H}$}(e_{k})=v_{i}\}$. here $f$ represents the self-dynamics, $\boldsymbol{z}_{i} \in \mathbb{R}^{d}$ represents the state vector of the $i^{\text{th}}$ dynamical system, $\varepsilon_{k}(t)$ is the positive coupling strength indicating how strongly system $j$ influences system $i$, and $\text{P}$ is a diagonal matrix with entries 0 or 1 that specifies which components of the dynamical system $j$ affect the dynamical system $i$. We have access to vary the coupling strengths each time the connectivity changes. For simplicity in computation, we omit the time variable $t$ under the assumption that the network connection remains static during the interval of interest. The information about network’s connections and coupling strengths is represented using a weighted digraph $\boldsymbol{\G}$ wherein the weight of the arc $(v_{i},v_{j})$ is equal to the coupling strength with which the system $i$ influences the system $j$. So, the terms \emph{coupling strength} and \emph{arc weight} . This digraph is referred to as the \emph{connectivity digraph}. We use the following assumptions :
\begin{enumerate}
    \item $ \forall \; \boldsymbol{x}, \boldsymbol{y} \in \mathbb{R}^d \; \exists$ a$>0$ such that \begin{equation}
(\boldsymbol{x} - \boldsymbol{y})^{\top} \left[f(\boldsymbol{x}) - f(\boldsymbol{y})-a\text{P}(\boldsymbol{x}-\boldsymbol{y})\right]\leq 0
\label{eq:Synchronization Assumption}
\end{equation}

\item $\G$ admits a directed spanning tree
\end{enumerate}

These two assumptions guarantee the global stability of the synchronization manifold i.e. $\{\boldsymbol{z}\big| \boldsymbol{z} = z\mathbf{1}_{d},z\in \mathbb{R}\}$ .
%staaligns with a condition established in \cite{belykh2004connection}, which ensures network synchronization i.e. $\lim_{t \to \infty} ||z_{i}(t) - z_{j}(t)||_{2} = 0 \quad \text{for all } i, j$ via coupling strength allocation (appendix of \cite{liu2013coupling}). Consequently, verifying whether a network can be synchronized reduces to the question of existence of a positive constant $a$.
The goal is to ensure network synchronization by finding a suitable arc weight vector $\E$ for $\G$. The typically used technique is to apply  Lyapunov function method and solve for the resulting matrix inequality which is of the form $\big(\text{U} \otimes \text{V}\big) \big( \text{L}_{\G}\otimes (-\text{P})-\text{I}_{n}\otimes\text{Y}\big)\preceq 0$. The authors of \cite{liu2015synchronization} using clever algebraic manipulations converted this into a simpler inequality comprised of laplacian matrices :
\begin{equation}
\text{L}_{\G}\text{L}_{\G_{0}}-a\text{L}_{\G_{0}}\succeq 0 
\label{eq:Graph inequality}
\end{equation}

This was again further simplified by using a complete graph L$_{\text{K}_{n}}$ in place of L$_{\G_{0}}$ resulting in a set of linear inequalities of the following form :
\begin{equation}
\frac{n}{2}\varepsilon_{k}^{s}\geq \smashoperator{\sum_{\substack{\{\PP \in \mathcal{S}:e_{k}\in \PP\}}}}||\PP|| \chi(\omega(\PP))
\label{eq:Synchronization Condition}
\end{equation}
Here, $\varepsilon_{k}^{s}$ equals $\varepsilon_{k}+\varepsilon_{l}$ only if $(e_{l},e_{k})$,$(e_{k},e_{l})$ $\in \EE(\G)$; otherwise, it is simply $\varepsilon_{k}$. $\omega(\PP)=\text{D}_{\mathcal{T}(\PP)}+\text{D}_{\mathcal{H}(\PP)}+2a$ (where $\chi(z)=z$ if $z>0$ and $0$ otherwise) and $\mathcal{S}$ is the set of undirected paths in the underlying graph of $\G$. A key difficulty in satisfying this synchronization condition lies in the interdependence between the arc weights and the vertex imbalances. Conventionally, one must first fix the vertex imbalance at each vertex before attempting to find a suitable arc weight vector $\E$ that satisfies the inequalities in \eqref{eq:Synchronization Condition}, constrained to the solution space defined by the equation Q$\E=\mathcal{D}$.\\
Our method construction of $\E$ removes the interdependence and lets us increase the arc weights till the synchronization condition in \eqref{eq:Synchronization Condition} is satisfied thus, eliminating the need for inequality solvers. The proposed method leverages two key remarks from the GCGM algorithm :

\begin{remark}
\label{thm:Negative Vertex Imbalance}
Any path $\PP$ whose end vertices have vertex imbalances $\leq -a$ in a weighted digraph $\G$ does not raise the lower bound on arc weights, since $\chi(\omega(\PP)) = 0$.
\end{remark}
\begin{remark}
\label{thm:SCC with 1 or 2 components}
The connectivity graph can be processed by the GCGM algorithm in segments, each segment being either :
\begin{enumerate}
\item A strongly connected component (SCC), or
\item A subgraph with two SCCs, where one of them is a source vertex.
\end{enumerate}
\end{remark}
The first remark implies that we need not compute the paths that have their head and tail vertex imbalances negative. This reduces the number of computations and motivates us to find an arc weight vector that can make most of the vertex imbalances to be less than or equal to $-a$. By the virtue of handshaking lemma for weighted digraphs, we observe that we can make at most $n-1$ among the $n$ vertices of a weighted digraph to have a negative vertex imbalance.
\begin{lemma}[Handshaking lemma for weighted digraphs]
    For a weighted digraph with $n$ vertices, 
    \begin{equation*}
    \sum_{i=1}^{n}\text{D}_{i}=0
    \end{equation*}
    \label{thm:Handshaking Lemma}
\end{lemma}
We call the arc weight vector with positive entries that can achieve maximum negative vertex imbalance for a weighted digraph as its \emph{negative imbalance arc weight vector} X$^{-}$.
The second remark is due to a result in \cite[section IV]{liu2015synchronization} that enables us to allocate arc weights to the connectivity digraph one SCC at a time. It also helps reduce computational complexity since we now need only compute paths within each SCC and not the whole digraph. For the rest of this paper, we refer to the single vertex SCCs in a digraph as \emph{trivial} SCCs and the others as \emph{non-trivial} SCCs. For each non-trivial SCC, the following operation is applied before being fed to the GCGM algorithm :
\begin{enumerate}
    \item If the SCC has no incoming arcs from other SCCs, it is directly input to the algorithm.
    \item Otherwise, a subgraph containing the SCC and its upstream SCCs is formed. These upstream SCCs are merged into a single vertex (combining parallel arcs), and the simplified subgraph is input to the algorithm.
\end{enumerate}
\begin{figure}[h]
\centering
\begin{tikzpicture}[>=stealth,scale=0.75, vertex/.style={circle, draw, fill=blue!20, minimum size=14pt, inner sep=0pt},arrow/.style={ultra thick}]
        \node[vertex] (v1) at (-4.5,-3.8) {1};
        \node[vertex] (v2) at (-1.5,-3.8) {2};
        \node[vertex] (v3) at (-1.5,-1.85) {3};
        \node[vertex] (v4) at (-4.5,-1.85) {4};
        
        \node[vertex] (v0) at (2,1) {0};
        \node[vertex] (v11) at (0.5,-3.8) {1};
        \node[vertex] (v21) at (3.5,-3.8) {2};
        \node[vertex] (v31) at (3.5,-1.85) {3};
        \node[vertex] (v41) at (0.5,-1.85) {4};
    
        \node[vertex] (v5) at (-3,2.5) {5};
        \node[vertex] (v6) at (-1.5,1.75) {6};
        \node[vertex] (v7) at (-1.5,0.25) {7};
        \node[vertex] (v8) at (-3,-0.5) {8};
        \node[vertex] (v9) at (-4.5,0.25) {9};
        \node[vertex] (v10) at (-4.5,1.75) {10};

        %\node[draw=none,above of = v5] (G) {$\G^{\text{SCC}}_{1}$};
        \node[draw=none,above of = v0,xshift = 1mm] (Gscc2tild) { };
 
        % Draw the surrounding circle
        \draw[red, ultra thick, dashed] (-3,1) circle (2.1);

        %\draw[->, thick,double] (-0.5,1) -- (2.5,1);
        
        \draw[->] (v1) -- (v2);
        \draw[->] (v2) -- (v3);
        \draw[->] (v3) -- (v4);
        \draw[->] (v4) -- (v1);
        \draw[->] (v3) -- (v1);

        \draw[->] (v9) -- (v4);
        \draw[->] (v5) -- (v1);
        \draw[->] (v7) -- (v3);
        \draw[->] (v8) -- (v3);
        \draw[->] (v8) -- (v2);
        
        \draw[->] (v6) -- (v7);
        \draw[->] (v8) -- (v6);
        \draw[->] (v7) -- (v8);
        
        \draw[->] (v5) -- (v9);
        \draw[->] (v9) -- (v10);
        \draw[->] (v10) -- (v5);
        
        \draw[->] (v7) -- (v9);
        \draw[->] (v5) -- (v8);

        \draw[->] (v11) -- (v21);
        \draw[->] (v21) -- (v31);
        \draw[->] (v31) -- (v41);
        \draw[->] (v41) -- (v11);
        \draw[->] (v31) -- (v11);
        
        \draw[->] (v0) -- (v11);
        \draw[->] (v0) -- (v21);
        \draw[->] (v0) -- (v31);
        \draw[->] (v0) -- (v41);
\end{tikzpicture}
\caption{The SCC of the digraph in the left (encircled by the dashed circle) is merged into vertex $0$, and the arcs $(7,3)$ and $(8,3)$ are replaced by the single arc $(0,3)$. This results in the digraph given on the right-hand side, which will be fed to the arc weight allocation algorithm.}
\label{fig:Example Condensation}
\end{figure}
From Remark~\ref{thm:Negative Vertex Imbalance}, minimizing the number of required paths is equivalent to ensuring that most vertex imbalances are below $-a$. Leveraging this observation, along with the structural properties of the input digraph used in the GCGM algorithm (as described in Remark~\ref{thm:SCC with 1 or 2 components}), we pose the following problem:

\emph{Given a weighted digraph $\G = (\VV, \EE, \E)$ that is either strongly connected or composed of two SCCs, one of which is a source vertex, the goal is to compute an arc weight vector $\E > 0$ such that:}
\begin{enumerate}
    \item The arc weights satisfy the constraints in \eqref{eq:Synchronization Condition}, and
    \item All but one vertex in $\G$ have imbalance less than or equal to $-a$.
\end{enumerate}

\section{Coupling Strength Allocation}
In this section, we show that any strongly connected weighted digraph admits an arc weight vector lying in the affine space defined by the sum of a negative vertex imbalance arc weight vector and its cycle space both with appropriate scaling. We prove that this vector satisfies the conditions outlined earlier. The result is then extended to digraphs with a source vertex and a non-trivial strongly connected component using a star-graph Laplacian.
\subsection{Arc Weight Assignment for Strongly Connected Digraphs}
First, we prove that any strongly connected weighted digraph admits a negative imbalance arc weight vector. We use the following proposition to help us in constructing one such vector.
\begin{proposition}
\label{thm:maximum negative vertex imbalance vector}
Let $\PP$ be a directed path in the weighted digraph $\G$, and $\E = [\varepsilon_{1}, \dots, \varepsilon_{m}]^{\top}$. Suppose:
\begin{enumerate}
    \item For adjacent arcs $e_i, e_j \in \EE(\PP)$ with $\mathcal{H}(e_j) = \mathcal{T}(e_i)$, we have $0 < \varepsilon_i < \varepsilon_j$.
    \item $\varepsilon_i = 0$ for all $e_i \in \EE(\G) \setminus \EE(\PP)$.
\end{enumerate}
Then the resulting vertex imbalance vector $\mathcal{D}$ satisfies:
\begin{enumerate}
    \item $\text{D}_k > 0$ if $v_k = \mathcal{T}(\PP)$,
    \item $\text{D}_k < 0$ for all $v_k \in \VV(\PP) \setminus \{\mathcal{T}(\PP)\}$,
    \item $\text{D}_k = 0$ for all $v_k \in \VV(\G) \setminus \VV(\PP)$.
\end{enumerate}
\label{thm:negative imbalance arc weight vector for a path}
\end{proposition}

\begin{proof}
Let $\PP$ be a path of length $\ell$ in $\G = (\VV, \EE, \E)$ with $|\G| = n$ and $||\G|| = m$. Relabel the vertices so that $\PP$ consists of arcs $e_1 = (v_0, v_1)$ to $e_{\ell} = (v_{\ell-1}, v_{\ell})$, with $v_0 = \mathcal{T}(\PP)$ and $v_{\ell} = \mathcal{H}(\PP)$. Define $\E$ such that only arcs in $\PP$ have non-zero weights, i.e.$\mathcal{D} = \begin{bmatrix} \varepsilon_1 & \varepsilon_2 - \varepsilon_1 & \cdots & -\varepsilon_{\ell} & 0 & \cdots & 0 \end{bmatrix}^{\top}$. Choosing $\varepsilon_i > 0$ and $\varepsilon_{i+1} > \varepsilon_i$ ensures $\text{D}_1 = \varepsilon_1 > 0$, $\text{D}_i < 0$ for $1 < i \leq \ell$, and $\text{D}_i = 0$ for $i > \ell$.
\end{proof}

\begin{corollary}
For a strongly connected digraph $\G$ with $|\G|=n$, there exists at least one negative imbalance arc weight vector such that only $n-1$ directed paths have positive weights
\label{thm:negative imbalance arc weight vector for strongly connected digraphs}
\end{corollary}

\begin{proof}
Any strongly connected digraph admits a directed spanning tree rooted at any vertex, ensuring directed paths from that vertex to all others. For each such path, we construct a vector $\text{X}^{-}_j$ as in Proposition~\ref{thm:path negative vertex imbalance}, where Q$\text{X}^{-}_{j}$ has a negative $j^\text{th}$ entry, a positive $r^\text{th}$ entry, and non-positive entries elsewhere. Since these paths span all vertices, the sum $\sum_{j \neq r} \text{QX}^{-}_{j}$ yields a vector with a positive $r^\text{th}$ entry and the rest strictly negative. By linearity, X$^{-} = \sum_{j \neq r}$X$_{j}^{-}$ defines a negative imbalance arc weight vector.
\end{proof}
To simplify computation, we revise the first condition in proposition~\ref{thm:negative imbalance arc weight vector for a path} as : for $e_i, e_j \in \EE(\PP)$ with $\mathcal{H}(e_j) = \mathcal{T}(e_i)$, let $0 < \varepsilon_i + 1 = \varepsilon_j$. We scale $\text{X}^{-}$ by $a$ to obtain $\E^{-}$. Using $\E^{-}$ as the arc weight vector for $\G$, makes all paths not starting from $v_r$ to have non-positive weights, i.e., $\omega(\PP_{i,j}) < 0$ $\forall$ $i,j \neq r$, and can be ignored. The synchronization condition in \eqref{eq:Synchronization Condition} then becomes:
For $e_{k} = (v_{i},v_{j})$ and $e_{l}= (v_{j},v_{i})\in \EE(\mathcal{\G}),$
\begin{equation}
\varepsilon_{k}^{s}\geq \frac{2a}{n}\Big(1+\smashoperator{\sum_{\PP \in \mathcal{S}}}||\PP||\Big)\Big(\smashoperator{\sum_{\PP \in \mathcal{S}}}||\PP||\Big)
%\varepsilon_{k}^{s}\geq\begin{cases}
%\frac{2a}{n}\Big(1+\smashoperator{\sum_{\PP\in \mathcal{S}}}||\PP||\Big)^{2},\text{ if }v_{i}\text{ or }v_{j} \in \{v_{r}\}\\
%\frac{2a}{n}\Big(1+\smashoperator{\sum_{\PP\in \mathcal{S}}}||\PP||\Big)\Big(\smashoperator{\sum_{\PP\in \mathcal{S}}}||\PP||\Big) \text{, otherwise}
%\end{cases}    
\label{eq:Modified Synchronization Condition}    
\end{equation}
%In this context, the term $\varepsilon_{k}^{s}$ denotes the sum $\varepsilon_{k} + \varepsilon_{l}$ if there exists an arc $e_{l}$ that is antiparallel to $e_{k}$; otherwise, $\varepsilon_{k}^{s}$ is defined to be $\varepsilon_{k}$ alone.\\
The vector $\E^{-}$ fixes the number of directed paths that need to be computed to be $n-1$, but does not necessarily help in satisfying the synchronization condition as described by the synchronization inequality \eqref{eq:Modified Synchronization Condition}. To accomplish that, we pick up  can pick up a vector from the cycle space of $\G$ so that the arc weights can be increased to meet the inequalities in \eqref{eq:Modified Synchronization Condition} without disturbing the vertex imbalances i.e. Q(X$^{0}$+$\E^{-}$)=Q$\E^{-}$ since X$^{0}\in \mathcal{C}$. We compute one such vector by finding the sum of the vectors of a cycle basis of $\G$ with non-negative entries. Such a cycle basis can be computed by completing the ears of the ear decomposition of $\G$. This basis, called the \emph{directed ear basis} \cite{loebl2001some}, comprises of a set of linearly independent vectors with each vector being the signed incidence vector corresponding to a directed cycle of the non-trivial SCC. Directed cycles from completed ears can be traversed to yield 0-1 signed incidence vectors. The last step is to find the sum of the \emph{directed ear basis} (X$^{0}$) and scale it by $\Delta \varepsilon = \frac{2a}{n}\Big(1+\smashoperator{\sum_{\PP \in \mathcal{S}}}||\PP||\Big)\Big(\smashoperator{\sum_{\PP \in \mathcal{S}}}||\PP||\Big)$ to obtain our desired vector $\E^{0}$. The vector ($\mathcal{E}_{sync}$) formed by summing $\E^{-}$ and $\E^{0}$, when used as the arc weight vector, can help us meet our synchronization condition and also reduce the number of paths to be computed.
\begin{equation}
    \mathcal{E}_{sync} = \mathcal{E}^{-}+\Delta \varepsilon \text{X}^{0}=\mathcal{E}^{-}+\mathcal{E}^{0}
\end{equation}
\subsection{Extension to Digraphs with Trivial and Non-trivial SCCs}
We prove that using a undirected star graph laplacian in the condition $\frac{1}{2} \left( \text{L}_{\G_0} \text{L}_{\G} + \text{L}_{\G}^\top \text{L}_{\G_0} \right) \succeq a \text{L}_{\G_0}$ which is equivalent to the condition \eqref{eq:Graph inequality} can give us simpler conditions on arc weights.
\begin{theorem}
Let $\G$ be a weighted digraph formed by attaching a new vertex with outgoing arcs to a strongly connected weighted digraph, and let $\text{L}_{\G}$ denote its Laplacian matrix. Let $\G_0$ be an undirected star graph on the same vertex set as $\G$, with Laplacian matrix $\text{L}_{\G_0}$. Then, $\G$ admits a negative imbalance arc weight vector such that
\begin{equation}
\frac{1}{2} \left( \text{L}_{\G_0} \text{L}_{\G} + \text{L}_{\G}^\top \text{L}_{\G_0} \right) \succeq a \, \text{L}_{\G_0},
\label{eq:symmetric-form}
\end{equation}
for some scalar $a > 0$.
\end{theorem}

%\begin{align}
%\varepsilon_{i,0} &\geq \frac{\text{D}_i}{2} + a \quad \text{for all } 1 \leq i \leq m_1,
%\label{eq:cond-eps1} \\
%\text{D}_j &\leq -2a \quad \text{for all } m_1 + 1 \leq j \leq n - 1,
%\label{eq:cond-Dj}
%\end{align}
%where D$_i$ denotes the vertex imbalance at $v_i$ in the non-trivial SCC.
\begin{proof}
Consider a weighted digraph $\G$ that has two SCCs - one trivial SCC that is a source vertex and one non-trivial SCC where $|\G|=n$ and $||\G||=m$. Using the corollary \ref{thm:negative imbalance arc weight vector for strongly connected digraphs}, we can find an arc weight vector that can make all except one vertex (arbitrary) of a strongly connected digraph to have vertex imbalances that do not exceed $-2a$. We find one such vector for the non-trivial SCC in $\G$. Let $v_{0}$ be the newly attached vertex, $\{e_{0,1},...,e_{0,m_{1}}\}$ be its incident arcs. Let the heads of these arcs be labeled $v_{1},\dots v_{m_{1}}$ ($m_{1} \leq n-1$). We can express $\frac{1}{2} \left( \text{L}_{\G_0} \text{L}_{\G} + \text{L}_{\G}^\top \text{L}_{\G_0} \right)$ as :
\begin{equation*}
\begin{aligned}
&\begin{bmatrix}
0 & \mathbf{0}_{n-1}^\top \\
\mathbf{0}_{n-1} & \frac{1}{2} \left( \text{L}_{\text{non-triv}} + \text{L}_{\text{non-triv}}^\top + \mathrm{diag}(\mathcal{D}) \right)
\end{bmatrix}&+\\
&\begin{bmatrix}
\text{D}_0 & \frac{1}{2} \mathcal{D}^\top - \E_{\text{triv}}^\top \\
\frac{1}{2} \mathcal{D} - \E_{\text{triv}} & \mathrm{diag}(\E_{\text{triv}} - \frac{1}{2} \mathcal{D})
\end{bmatrix}&
\end{aligned}
\label{eq:lhs-decomposition}
\end{equation*}
The first matrix above is positive semi-definite, being a symmetric Laplacian. Hence, it suffices to show that
\begin{equation*}
\begin{bmatrix}
\text{D}_0 & \frac{1}{2} \mathcal{D}^\top - \E_{\text{triv}}^\top \\
\frac{1}{2} \mathcal{D} - \E_{\text{triv}} & \mathrm{diag}(\E_{\text{triv}} - \frac{1}{2} \mathcal{D})
\end{bmatrix}
\succeq
a \begin{bmatrix}
n-1 & -\mathbf{1}^\top_{n-1} \\
-\mathbf{1}_{n-1} & \text{I}_{n-1}
\end{bmatrix}
\label{eq:reduced-matrix-inequality}
\end{equation*}
Since both sides in the above equation share the same sparsity structure with zero row sum, both of them can be written as a linear combination of the elementary laplacians of an undirected star graph. The coefficients on the left hand side are : 
\begin{align*}
\varepsilon_{0,i}-\dfrac{\text{D}_{i}}{2} \:\: \forall \:\: 1 \leq i \leq m_{1}\\
\text{ and }
-\dfrac{\text{D}_{i}}{2} \:\: \forall \:\: m_{1}+1 \leq i \leq n-1
\end{align*}
and all those on the right hand side are equal to $a$. The inequality \eqref{eq:symmetric-form} holds if the coefficients on the left hand side are greater than or equal to $a$. Since the vertex imbalances except for D$_{1}$ are already lesser than $-2a$ we need only make sure that
\begin{equation}
\varepsilon_{0,1} \geq \frac{\text{D}_1}{2} + a, \text{ and } \varepsilon_{0,i} \geq 0 \quad \forall\, 1<i \leq n-1,
\label{eq:eps-final}
\end{equation}
which are feasible given that $\varepsilon_{i,0}$ are free design parameters.
\end{proof}
\subsection{Example}
We consider a digraph $\G$ depicted in figure \ref{fig:Example Condensation}, which consists of two SCCs. The first SCC is the $\G_{1}$, includes all vertices labeled from 5 to 10, along with the arcs connecting them as shown in the figure \ref{fig:SCC1}. The second SCC includes all vertices labeled from 1 to 4 and their respective incident arcs.
\begin{figure}[h]
\centering
\begin{subfigure}[t]{0.2\textwidth}
\begin{tikzpicture}[>=stealth,scale=0.8, vertex/.style={circle, draw, fill=blue!20, minimum size=15pt, inner sep=0pt}]
\node[vertex] (v5) at (0,2.5) {5};
\node[vertex] (v6) at (1.5,1.75) {6};
\node[vertex] (v7) at (1.5,0.25) {7};
\node[vertex] (v8) at (0,-0.5) {8};
\node[vertex] (v9) at (-1.5,0.25) {9};
\node[vertex] (v10) at (-1.5,1.75) {10};

\draw[->] (v5) -- (v8) node[near start,above right] {$1$};
\draw[->] (v5) -- (v9) node[midway,below right] {$2$};
\draw[->] (v6) -- (v7) node[midway,right] {$3$};
\draw[->] (v7) -- (v8) node[near start,below] {$4$};
\draw[->] (v7) -- (v9) node[near end,below,anchor=north east] {$5$};
\draw[->] (v8) -- (v6) node[midway,above] {$6$};
\draw[->] (v9) -- (v10) node[midway,anchor=east] {$7$};
\draw[->] (v10) -- (v5) node[midway,anchor=south east] {$8$};
\end{tikzpicture}
\caption{$\G_{1}$}
\label{fig:SCC1}
\end{subfigure}\hfill\centering
\begin{subfigure}[t]{0.2\textwidth}
\begin{tikzpicture}[>=stealth,scale=0.8, vertex/.style={circle, draw, fill=blue!20, minimum size=15pt, inner sep=0pt}]
\node[vertex] (v1) at (-0.5,-3.8) {1};
\node[vertex] (v2) at (2.5,-3.8) {2};
\node[vertex] (v3) at (2.5,-1.85) {3};
\node[vertex] (v4) at (-0.5,-1.85) {4};

\draw[->] (v1) -- (v2) node[midway,below] {$5$};
\draw[->] (v2) -- (v3) node[midway,right] {$6$};
\draw[->] (v3) -- (v4) node[midway,yshift=0.075cm,below] {$7$};
\draw[->] (v4) -- (v1) node[midway,left] {$8$};
\draw[->] (v3) -- (v1) node[midway,below,xshift=0.1cm,yshift=0.15cm] {$9$};
\end{tikzpicture}
\caption{$\G_{2}$}
\label{fig:SCC2}
\end{subfigure}
\begin{subfigure}[b]{0.2\textwidth}
\begin{tikzpicture}[>=stealth,scale=0.8, vertex/.style={circle, draw, fill=blue!20, minimum size=15pt, inner sep=0pt}]
\node[vertex] (v0) at (1,0.5) {0};
\node[vertex] (v1) at (-0.5,-3.8) {1};
\node[vertex] (v2) at (2.5,-3.8) {2};
\node[vertex] (v3) at (2.5,-1.85) {3};
\node[vertex] (v4) at (-0.5,-1.85) {4};

\draw[->] (v0) -- (v1) node[near start,below,xshift=0.1cm] {$1$};
\draw[->] (v0) -- (v2) node[near start,below,xshift=-0.1cm] {$2$};
\draw[->] (v0) -- (v3) node[midway,above,xshift=0.15cm,yshift=-0.1cm] {$3$};
\draw[->] (v0) -- (v4) node[midway,above,xshift=-0.15cm,yshift=-0.1cm] {$4$};

\draw[->] (v1) -- (v2) node[midway,below] {$5$};
\draw[->] (v2) -- (v3) node[midway,right] {$6$};
\draw[->] (v3) -- (v4) node[midway,yshift=0.075cm,below] {$7$};
\draw[->] (v4) -- (v1) node[midway,left] {$8$};
\draw[->] (v3) -- (v1) node[midway,below,xshift=0.1cm,yshift=0.15cm] {$9$};
\end{tikzpicture}
\caption{$\tilde{\G}_{2}$}
\label{fig:processed SCC2}
\end{subfigure}
\label{fig:SCCs of G}
\caption{Two SCCs of $\G$: $\G_1$ and $\G_2$, with their respective edge labels (relabelled for processing) are shown in the figures \ref{fig:SCC1} and \ref{fig:SCC2}. Figure \ref{fig:processed SCC2}  shows $\G_2$ with its upstream SCCs condensed into a single vertex $0$}
\end{figure}

\begin{table}[h]
\centering
\begin{tabular}{|c|c|}
\hline
 & \textbf{Negative Imbalance}\\
\textbf{Directed Path}&\textbf{Arc Weight Vector}\\
&(for paths)\\\hline
$5\rightarrow8\rightarrow6$&$\left[\;2\;0\;0\;0\;0\;1\;0\;0\;\right]^{\top}$\\\hline
$5\rightarrow8\rightarrow6\rightarrow7$&$\left[\;3\;0\;1\;0\;0\;2\;0\;0\;\right]^{\top}$ \\ \hline
$5\rightarrow8$&$\left[\;1\;0\;0\;1\;0\;0\;0\;0\;\right]^{\top}$\\\hline
$5\rightarrow9$&$\left[\;0\;1\;0\;0\;0\;0\;0\;0\;\right]^{\top}$\\\hline
$5\rightarrow9\rightarrow10$&$\left[\;0\;2\;0\;0\;0\;0\;1\;0\;\right]^{\top}$\\\hline
\end{tabular}
\label{tab: Negative imbalance arc weight vector for G1}
\caption{Directed paths in $\G_{1}$ and their respective vectors\\for negative imbalance arc weight vector computation}
\end{table}
 We computed directed paths from vertex $5$ to all the other vertices in $\G_{1}$ and computed the negative imbalance arc weight vector to be : $\E^{-}=a[6\:3\:1\:1\:0\:3\:1\:0]^\top$
\begin{table}[h]
\centering
\begin{tabular}{|c|c|}
\hline
\textbf{Directed cycle} & \textbf{Signed Incidence Vector}\\\hline
$5\rightarrow9\rightarrow10\rightarrow5$ & $\left[\;0\;1\;0\;0\;0\;0\;1\;1\;\right]^{\top}$\\\hline
$6\rightarrow7\rightarrow8\rightarrow6$ & $\left[\;0\;0\;1\;1\;0\;1\;0\;0\;\right]^{\top}$\\\hline $5\rightarrow8\rightarrow6\rightarrow7\rightarrow9\rightarrow10\rightarrow5$ & $\left[\;1\;0\;1\;0\;1\;1\;1\;1\;\right]^{\top}$\\\hline
\end{tabular}
\label{tab:Cycle basis sum for G1}
\caption{Sum of cycle basis for $\G_{1}$}
\end{table}
A cycle basis for $\G_{1}$ is listed in the table \ref{tab:Cycle basis sum for G1}, and from this, we determine its summation i.e. X$^{0}$. We computed $\Delta \varepsilon$ to be equal to $\dfrac{20a}{6}$ and scale X$^{0}$ to get  $\E^{0}= \dfrac{20a}{6}[1\:1\:2\:1\:1\:2\:2\:2]^\top$. The second component, $\G_{2}$, consists of all vertices labeled from 1 to 4, along with their respective arcs. We construct a digraph $\tilde{\G}_{2}$ by condensing $\G_{1}$ into a single vertex (vertex $0$) while preserving the arcs from $\G_{1}$ to $\G_{2}$, as illustrated in \cref{fig:Example Condensation}. 
For $\tilde{\G}_{2}$, we find the negative imbalance arc weight vector for its SCC ($G_{2}$) by finding directed paths from vertex "1" and scale the sum by $2a$. This results in the following vertex imbalances : D$_{1}=12a$, D$_{2}=-6a$, D$_{3}=-4a$ and D$_{4}=-2a$. Fixing vertex $0$, X$^{-}$ is computed to be $[\;6\;3\;0\;1\;0]^{\top}$.

\begin{table}[h]
\centering
\begin{tabular}{|c|c|}
\hline
 & \textbf{Negative Imbalance}\\
\textbf{Directed Path}&\textbf{Arc Weight Vector}\\
&(for paths)\\\hline
$1\rightarrow2$&$[\;1\;0\;0\;0\;0\;]^{\top}$\\\hline
$1\rightarrow2\rightarrow3$&$[\;2\;1\;0\;0\;0\;]^{\top}$ \\ \hline
$1\rightarrow2\rightarrow3\rightarrow4$&$[\;3\;2\;0\;1\;0\;]^{\top}$\\\hline
\end{tabular}
\label{tab: Negative imbalance arc weight vector for G2}
\caption{Directed paths in $\G_{2}$ and their respective vectors\\for negative imbalance arc weight vector computation}
\end{table}

\begin{table}[h]
\begin{tabular}{|c|c|}
\hline
\textbf{Directed Cycle} & \textbf{Signed Incidence Vector}\\\hline
$1\rightarrow2\rightarrow3\rightarrow4\rightarrow1$ & $\left[\;0\;0\;0\;0\;1\;1\;1\;1\;0\;\right]^{\top}$\\\hline
$1\rightarrow2\rightarrow3\rightarrow1$ & $\left[\;0\;0\;0\;0\;1\;1\;0\;0\;1\;\right]^{\top}$\\\hline
\end{tabular}
\label{tab:G2 vectors}
%\caption{}
%\label{tab:Circuit Basis of G2} 
\caption{Weighted Incidence Vectors of directed paths from $v_{0}$ in $\tilde{\G}_{2}$ and cycle basis of $\tilde{\G}_{2}$}
\end{table}


\section{Simulation and Results}
\subsection{Simulation}
\begin{figure}[h]
    \centering
    \includegraphics[scale=0.8]{Synchronization_Image.pdf}
    \caption{10 Lorentz oscillators with connectivity digraph given in \cref{fig:Example Condensation} synchronize due to coupling strength allocation by the enhanced algorithm.}
    \label{fig:Error Convergence}
\end{figure}
For each vertex in the digraph $\G$ of \cref{fig:Example Condensation}, we considered the Lorentz attractor described by :
\begin{equation}
\begin{aligned}
  \dot{x} &= \sigma (y - x), &
  \dot{y} &= x (r - z) - y, &
  \dot{z} &= x y - b z,&
\end{aligned}
\end{equation}
with the following parameters $\sigma = 10$, $r = 28$ and $b = 8/3$ (chaotic behaviour) and computed the system parameter $a$ using the formula $a = \dfrac{b(b+1)(r+\sigma)^{2}}{16(b-1)} -\sigma$ as in \cite{liu2013coupling}. The initial conditions for each system were assigned randomly. At each time step, the squared distances between state i and all other states were summed, and this sum was plotted. The sum consistently decreases to zero as seen in \cref{fig:Error Convergence} indicating network synchronization.

\subsection{Complexity}
For each subgraph with $n$ vertices and $m$ arcs formed from an SCC, it involves computing $n-1$ directed paths along with determining a circuit basis. This task has the same computational complexity as finding a cycle basis using a spanning tree, as discussed in \cite{ryan1981comparison}, which is : $\mathcal{O}(n\gamma)$, where $\gamma = m - n + 1$. In the worst-case scenario, where $m = \frac{n(n-1)}{2}$, the time complexity increases to $\mathcal{O}(n^3)$.

\section{Conclusion}
\label{sec:conclusion}
We introduced a novel method to determine coupling strengths that guarantee synchronization in time-varying dynamical networks. By exploiting the characteristics of strongly connected digraphs, we showed that our method assigns arc weights that minimizes the number of directed paths that need to be calculated to $n-1$. Furthermore, we used the cycle basis for each strongly connected component (SCC) to meet the synchronization criterion. However, since the method relies on strong connectivity, it is not applicable to connectivity graphs lacking a cycle basis, such as trees. Despite this, our algorithm enhances existing methodologies and extends their applicability to more intricate and dynamically evolving network structures. In summary, the proposed solution offers an efficient computational framework for achieving synchronization in time-varying networks and sets the stage for future developments in network control and dynamics.

\section*{Acknowledgment}

The authors acknowledge the use of language assistance tools, including ChatGPT, DeepSeek, and Writefull, for paraphrasing suggestions and grammar refinement during manuscript preparation. These tools were utilized solely for linguistic enhancement and not for the generation of any technical or scientific content.

	
\bibliographystyle{ieeetr}
\bibliography{References}
		
\end{document}
%\begin{figure}[h]
%\centering
%\begin{tikzpicture}[>=stealth,scale=0.75, vertex/.style={circle, draw, fill=blue!20, minimum size=15pt, inner sep=0pt},arrow/.style={ultra thick}]
    % Compute rotated coordinates manually (y remains unchanged)
    %\node[vertex] (a1) at (-4.5,-3.8) {$1$};
    %\node[draw=none,xshift=-0.19cm,yshift=-0.75cm] (d1) at (a1.north east) {D$_{1}=+3$};
    %\node[vertex] (a2) at (-1.5,-3.8) {$2$};
    %\node[draw=none,xshift=0.19cm,yshift=-0.75cm] (d2) at (a2.north west) {D$_{2}=-1$};
    %\node[vertex] (a3) at (-1.5,-1.85) {$3$};
    %\node[draw=none,xshift=0.19cm,yshift=0.75cm] (d3) at (a3.south west) {D$_{3}=-1$};
    %\node[vertex] (a4) at (-4.5,-1.85) {$4$};
    %\node[draw=none,xshift=-0.19cm,yshift=0.75cm] (d4) at (a4.south east) {D$_{4}=-1$};

    % Edges for Ga
    %\draw[->] (a1) -- (a2) node[midway,anchor=north] {$3$};
    %\draw[->] (a2) -- (a3) node[midway,anchor=west] {$2$};
    %\draw[->] (a3) -- (a4) node[midway,anchor=south] {$1$};
    %\draw[->] (a4) -- (a1) node[midway,anchor=east] {$0$};
    %\draw[->] (a3) -- (a1) node[midway,anchor= east] {$0$};

    %\node[draw=none,yshift=-0.4cm] (Vc1) at ($(a1) + (0.3*\xshift,-0.2*\xshift)$) {X$=\left[\;3\;2\;1\;0\;0\;\right]^{\top}$};
%\end{tikzpicture}
%\caption{For the directed path $\mathcal P: 1\rightarrow2\rightarrow3\rightarrow4$ in the digraph shown above, the weighted incidence vector X $= [3\,2\,1\,0\,0]^{\top}$ results in all non-tail vertices of $\mathcal P$ having a negative vertex imbalance.}
%\label{fig:Circuit Decompositon of G}
%\end{figure}

%Any change in an arc weight vector might change the vertex imbalance vector through \eqref{eq:vertex imbalance equality} which in turn might result in some of the inequalities of \eqref{eq:Modified Synchronization Condition} to not be satisfied. So, we need increase the arc weights such that the inequalities of \eqref{eq:Modified Synchronization Condition} are met without disturbing the vertex imbalances. This is 

