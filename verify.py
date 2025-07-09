import multiprocessing as mp
import sys
import time

from itertools import combinations
from sage.all import *

def k_orbits(G, k):
    A = G.automorphism_group()

    if k == 1:
        return [tuple((v,) for v in G.vertices())]

    V = list(G.vertices())

    objs   = list(frozenset(c) for c in combinations(V, k))
    action = "OnSets"
    
    seen   = set()
    result = []

    for x in objs:
        if x in seen:
            continue
        orb = tuple(A.orbit(x, action=action))  # one orbit at a time
        result.append(orb)
        seen.update(orb)

    return result


def _check_for_conflicts(edges, mult, vtr, alpha):
    graph = Graph(edges)
    subgraph = graph.subgraph([v for v in graph.vertices() if v not in vtr])
    subgraph_alpha = subgraph.complement().clique_number()
    return (subgraph_alpha <= alpha - mult)

def check_for_conflicts(graph, k, generalized=False, pool=None):
    graph_edges = [(u, v) for u, v, _ in graph.edges()]
    degree_sequence = graph.degree_sequence()

    if set(degree_sequence) != set([k]):
        raise ValueError(f"Degree sequence not k-regular: {degree_sequence}")

    vertices_to_remove = []

    LG = graph.line_graph()
    for mult in range(1, k) if generalized else [1]:
        orbits = k_orbits(LG, mult*(k-1))
        edges_to_remove = [orbit[0] for orbit in orbits]

        for edges in edges_to_remove:
            assert all(graph.has_edge(e) for e in edges), f"{edges}"
            vtr = set([edge[0] for edge in edges] + [edge[1] for edge in edges])
            vertices_to_remove.append((mult, vtr))

    alpha = graph.complement().clique_number()

    args = [(graph_edges, mult, vtr, alpha) for mult, vtr in vertices_to_remove]
    results = [None] * len(args)
    pending = [pool.apply_async(_check_for_conflicts, arg) for arg in args]

    completed = set()
    while len(completed) < len(pending):
        updated = False
        for i, res in enumerate(pending):
            if i in completed:
                continue
            if res.ready():
                results[i] = res.get()
                completed.add(i)
                updated = True
        if not updated:
            time.sleep(0.01)
    
    return not any(results)

def get_generalized_petersen_graph(n, k):
    edges = []

    for i in range(n):
        edges.append(((0, i), (0, (i + 1) % n)))
        edges.append(((1, i), (1, (i + k) % n)))

    for i in range(n):
        edges.append(((0, i), (1, i)))

    return Graph(edges)

if __name__ == "__main__":

    with mp.Pool() as pool:

        print("\nChecking the three base cases for Theorem 1 ...")
        for k in [0, 1, 2]:
            graph = get_generalized_petersen_graph(5*k + 11, 2)
            print(f"  Checking k={k} ", end="")
            print(f"✅" if (min(graph.degree_sequence()) == max(graph.degree_sequence()) == 3) else f"❌", end=" ")
            print(f"✅" if graph.line_graph().chromatic_number() == 3 else f"❌", end=" ")
            print(f"✅" if check_for_conflicts(graph, 3, pool=pool) else f"❌")

        print("\nChecking r=3 counterexamples ...")
        with open("counterexamples_r3.txt", "r") as f:
            for line in f:
                graph = Graph(line.strip())
                print(f"  Checking graph of order {graph.order()} ", end="")
                print(f"✅" if (min(graph.degree_sequence()) == max(graph.degree_sequence()) == 3) else f"❌", end=" ")
                print(f"✅" if graph.line_graph().chromatic_number() == 3 else f"❌", end=" ")
                print(f"✅" if check_for_conflicts(graph, 3, pool=pool) else f"❌")

        print("\nChecking r=4 counterexamples ")
        with open("counterexamples_r4.txt", "r") as f:
            for line in f:
                graph = Graph(line.strip())
                print(f"  Checking graph of order {graph.order()} ", end="")
                print(f"✅" if (min(graph.degree_sequence()) == max(graph.degree_sequence()) == 4) else f"❌", end=" ")
                print(f"✅" if graph.line_graph().chromatic_number() == 4 else f"❌", end=" ")
                print(f"✅" if check_for_conflicts(graph, 4, pool=pool) else f"❌")

