#!/usr/bin/env python3
"""Choose which samples to probe, using only the profiles.

A contested group is resolved by a sample that contains SOME but not ALL of its
members. If a sample contains every member, the contested k-mers stay assigned
whichever member owns them and the flip leaves no trace in the novel section; if
it contains none, the k-mers are novel either way and rule nothing out. Only the
partial overlap is informative.

Profiles say which genomes each sample contains, so this choice costs no reads
and opens no .sylspr. Samples with the fewest members of a group present are
preferred, since each one rules out more of the group at once.
"""

import argparse
import collections
import csv
import os
import sys


def read_list(path):
    with open(path) as fh:
        return [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]


def stem(path):
    return os.path.basename(path).split(".sylsp")[0]


def load_graph_tsv(path):
    gid_of, groups = {}, {}
    with open(path) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            if row[0] == "genome":
                gid_of[row[2]] = int(row[1])
                gid_of.setdefault(os.path.basename(row[2]), int(row[1]))
            elif row[0] == "group":
                groups[int(row[1])] = [int(x) for x in row[2].split(",")]
    return gid_of, groups


def load_profiles(paths, gid_of):
    contains = collections.defaultdict(set)
    for prof in paths:
        with open(prof) as fh:
            for rec in csv.DictReader(fh, delimiter="\t"):
                g = gid_of.get(rec["Genome_file"])
                if g is None:
                    g = gid_of.get(os.path.basename(rec["Genome_file"]))
                if g is not None:
                    contains[rec["Sample_file"]].add(g)
    return contains


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph-tsv", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--profiles", required=True)
    ap.add_argument("--per-group", type=int, default=4)
    ap.add_argument("--out", required=True)
    ap.add_argument("--report", required=True)
    args = ap.parse_args()

    gid_of, groups = load_graph_tsv(args.graph_tsv)
    contains = load_profiles(read_list(args.profiles), gid_of)

    # Only samples that actually have a .sylspr carry evidence; a collection
    # normally holds plain sketches too, from samples too diverse to have been
    # reference-compressed.
    have = {stem(s): s for s in read_list(args.samples) if s.endswith(".sylspr")}

    picked, per_group = set(), collections.Counter()
    resolvable = 0
    for gi, members in sorted(groups.items()):
        ms = set(members)
        cands = []
        for sample, genomes in contains.items():
            inside = ms & genomes
            if inside and len(inside) < len(ms):
                path = have.get(stem(sample))
                if path:
                    cands.append((len(inside), sample, path))
        if cands:
            resolvable += 1
        cands.sort()
        for _, _, path in cands[: args.per_group]:
            picked.add(path)
            per_group[gi] += 1

    with open(args.out, "w") as fh:
        for p in sorted(picked):
            fh.write(p + "\n")
    with open(args.report, "w") as fh:
        fh.write("groups\tresolvable_from_profiles\tselected_samples\n")
        fh.write(f"{len(groups)}\t{resolvable}\t{len(picked)}\n")

    if resolvable < len(groups):
        sys.stderr.write(
            f"WARNING: {len(groups) - resolvable} of {len(groups)} contested groups have no "
            "sample containing some-but-not-all of their members. Those cannot be resolved "
            "from evidence and are left to `ref-recover search`.\n"
        )


if __name__ == "__main__":
    main()
