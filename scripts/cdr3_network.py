#!/usr/bin/env python3

import pandas as pd
import networkx as nx
import argparse
import json
from itertools import combinations


def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        return None
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def build_network(df, max_dist=1):

    G = nx.Graph()

    # clone size
    clone_size = df.groupby("cdr3").size().to_dict()

    for cdr3, size in clone_size.items():
        G.add_node(cdr3, size=int(size))

    # group by VJ + length
    df["cdr3_len"] = df["cdr3"].str.len()

    groups = df.groupby(["v_call", "j_call", "cdr3_len"])

    for (v, j, l), group in groups:

        cdr3_list = group["cdr3"].unique()

        for c1, c2 in combinations(cdr3_list, 2):

            d = hamming_distance(c1, c2)

            if d is not None and d <= max_dist:
                G.add_edge(
                    c1,
                    c2,
                    v_call=v,
                    j_call=j,
                    distance=d
                )

    return G


def network_stats(G):

    stats = {}

    stats["nodes"] = G.number_of_nodes()
    stats["edges"] = G.number_of_edges()

    if G.number_of_nodes() > 0:

        comps = list(nx.connected_components(G))
        stats["components"] = len(comps)

        largest = max(comps, key=len)
        stats["largest_component"] = len(largest)

    else:
        stats["components"] = 0
        stats["largest_component"] = 0

    return stats


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--input", required=True)
    parser.add_argument("--edges", required=True)
    parser.add_argument("--nodes", required=True)
    parser.add_argument("--json", required=True)
    parser.add_argument("--stats", required=True)

    parser.add_argument("--max-dist", type=int, default=1)

    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")

    df = df.dropna(subset=["cdr3", "v_call", "j_call"])

    G = build_network(df, args.max_dist)

    # edges
    edges = nx.to_pandas_edgelist(G)
    edges.to_csv(args.edges, sep="\t", index=False)

    # nodes
    nodes = pd.DataFrame([
        {"cdr3": n, **G.nodes[n]}
        for n in G.nodes()
    ])

    nodes.to_csv(args.nodes, sep="\t", index=False)

    # json
    data = nx.node_link_data(G)

    with open(args.json, "w") as f:
        json.dump(data, f)

    # stats
    stats = network_stats(G)

    with open(args.stats, "w") as f:
        json.dump(stats, f, indent=2)


if __name__ == "__main__":
    main()