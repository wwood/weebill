#!/usr/bin/env python3
"""Turn probe evidence into a `ref-build --genome-order` file.

Each contested group is won by whichever of its members came first in the lost
reference's build order. The probe gives one-directional evidence: a group's
k-mer showing up in a sample's *novel* section proves the winner is a member that
sample did not contain, because a k-mer owned by a genome the sample hit would
have been encoded against that genome instead. Absence of the k-mer from the
novel section proves nothing -- the sample may simply not have had it -- so only
the positive observations are used.

Combining that with the sample's genome content (from its profile) gives, per
observation, `winner in members - contained`. Intersecting those candidate sets
across samples narrows each group, and a singleton is a solved group: a strict
"this member precedes all the others" constraint. Those constraints are
topologically sorted into a total order, and any group left unresolved is
reported rather than guessed at, since the fingerprint check downstream is what
decides whether the answer is right.
"""

import argparse
import collections
import csv
import os
import subprocess
import sys


def read_list(path):
    with open(path) as fh:
        return [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph-tsv", required=True)
    ap.add_argument("--evidence", required=True)
    ap.add_argument("--profiles", required=True)
    ap.add_argument("--graph", required=True)
    ap.add_argument("--weebill", default="weebill")
    ap.add_argument("--out", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument(
        "--solved-out",
        help="Write the settled groups as <group_id><TAB><winner> for "
        "`weebill ref-recover search` to take as given.",
    )
    ap.add_argument(
        "target_fingerprint",
        nargs="?",
        default="",
        help="If given, the recovered order is scored against it analytically "
        "before any reference is built.",
    )
    args = ap.parse_args()

    names = {}
    gid_of = {}
    groups = {}
    with open(args.graph_tsv) as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            if row[0] == "genome":
                gid = int(row[1])
                names[gid] = row[2]
                gid_of[row[2]] = gid
                gid_of.setdefault(os.path.basename(row[2]), gid)
            elif row[0] == "group":
                groups[int(row[1])] = [int(x) for x in row[2].split(",")]

    # sample -> genomes it contains
    contains = collections.defaultdict(set)
    for prof in read_list(args.profiles):
        with open(prof) as fh:
            for rec in csv.DictReader(fh, delimiter="\t"):
                g = gid_of.get(rec["Genome_file"]) or gid_of.get(
                    os.path.basename(rec["Genome_file"])
                )
                if g is not None:
                    contains[rec["Sample_file"]].add(g)

    def contained_for(sylspr_path):
        stem = os.path.basename(sylspr_path).split(".sylsp")[0]
        for sample, genomes in contains.items():
            if os.path.basename(sample).split(".sylsp")[0] == stem:
                return genomes
        return None

    # winner candidates per group, intersected over observations
    candidates = {gi: set(ms) for gi, ms in groups.items()}
    observations = collections.Counter()
    missing_profile = set()
    with open(args.evidence) as fh:
        for rec in csv.DictReader(fh, delimiter="\t"):
            if not rec.get("group_id"):
                continue
            gi = int(rec["group_id"])
            held = contained_for(rec["sample_file"])
            if held is None:
                missing_profile.add(rec["sample_file"])
                continue
            # the k-mer was novel, so its owner is not one of this sample's hits
            narrowed = candidates[gi] - held
            if narrowed:
                candidates[gi] = narrowed
                observations[gi] += 1

    solved = {gi: next(iter(c)) for gi, c in candidates.items() if len(c) == 1}
    unresolved = {gi: sorted(c) for gi, c in candidates.items() if len(c) != 1}

    # winner precedes every other member
    succ = collections.defaultdict(set)
    indeg = collections.Counter({g: 0 for g in names})
    for gi, w in solved.items():
        for m in groups[gi]:
            if m != w and m not in succ[w]:
                succ[w].add(m)
                indeg[m] += 1

    # Kahn, breaking ties by genome id so the output is reproducible
    import heapq

    ready = [g for g in sorted(names) if indeg[g] == 0]
    heapq.heapify(ready)
    order = []
    while ready:
        g = heapq.heappop(ready)
        order.append(g)
        for m in sorted(succ[g]):
            indeg[m] -= 1
            if indeg[m] == 0:
                heapq.heappush(ready, m)

    cyclic = len(order) != len(names)
    if cyclic:
        # Contradictory evidence: keep going with the genomes that did sort, then
        # append the rest by id, and let the fingerprint check be the judge.
        placed = set(order)
        order += [g for g in sorted(names) if g not in placed]

    with open(args.out, "w") as fh:
        for g in order:
            fh.write(names[g] + "\n")

    if args.solved_out:
        with open(args.solved_out, "w") as fh:
            for gi, w in sorted(solved.items()):
                fh.write(f"{gi}\t{w}\n")

    scored = ""
    if args.target_fingerprint:
        got = subprocess.run(
            [
                args.weebill,
                "ref-recover",
                "fingerprint",
                "--graph",
                args.graph,
                "--genome-order",
                args.out,
            ],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        scored = got
        sys.stderr.write(
            f"analytic fingerprint {got}, target {args.target_fingerprint}: "
            + ("MATCH\n" if got == args.target_fingerprint else "MISMATCH\n")
        )

    with open(args.report, "w") as fh:
        fh.write(f"genomes\t{len(names)}\n")
        fh.write(f"contested_groups\t{len(groups)}\n")
        fh.write(f"groups_with_evidence\t{len(observations)}\n")
        fh.write(f"groups_solved\t{len(solved)}\n")
        fh.write(f"groups_unresolved\t{len(unresolved)}\n")
        fh.write(f"cyclic_evidence\t{cyclic}\n")
        if scored:
            fh.write(f"analytic_fingerprint\t{scored}\n")
            fh.write(f"target_fingerprint\t{args.target_fingerprint}\n")
        if missing_profile:
            fh.write(f"samples_without_profile\t{len(missing_profile)}\n")
        for gi, cands in sorted(unresolved.items()):
            fh.write(
                "unresolved_group\t{}\tmembers={}\tcandidates={}\n".format(
                    gi,
                    ",".join(str(m) for m in groups[gi]),
                    ",".join(str(m) for m in cands),
                )
            )

    if unresolved:
        sys.stderr.write(
            f"{len(unresolved)} of {len(groups)} contested groups unresolved; the "
            "order written is a best effort. Each unresolved group is one "
            "remaining choice -- score candidates with "
            "`weebill ref-recover fingerprint` rather than rebuilding each one.\n"
        )


if __name__ == "__main__":
    main()
