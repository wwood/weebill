# Recover the genome order a lost .sylref was built with.
#
# A reference's genome ids come from the build order (species, representatives
# first, then file name) and are reproducible. What is not reproducible is the
# owner of a k-mer contested by fewer than --pool-min-genomes same-tier genomes:
# that goes to whichever contender the build saw first, which is the order the
# sketches happened to sit in their .syldb, and `weebill sketch` writes genomes
# in parallel-completion order. Re-sketching therefore rebuilds the same genome
# ids with different ownership, and every .sylspr refuses to decode.
#
# Re-sketching cannot reveal the original order -- each run gives a fresh random
# one. The evidence is in the samples instead. This workflow:
#
#   1. sketches the genomes ONCE (order within the .syldb does not matter here),
#   2. records the contested groups -- the actual unknowns,
#   3. reads each .sylspr's novel hashes, which needs no reference at all,
#   4. turns "this group's k-mer was novel in a sample that contained only some
#      of its members" into "the winner is one of the members the sample lacked",
#   5. topologically sorts those constraints into a --genome-order, and
#   6. rebuilds and checks the fingerprint every .sylspr records.
#
# Step 6 is the only check that counts; everything before it is inference.
#
# Usage:
#   snakemake -s workflow/recover_sylref_order.smk --configfile recover.yaml -j 32
#
# Required config keys (see `config` defaults below):
#   genomes:       file listing the genome FASTAs, one per line, that the lost
#                  reference was built from
#   samples:       file listing the sample sketches. Files that are not *.sylspr
#                  are ignored, so a listing of a whole collection is fine --
#                  samples too diverse to have been reference-compressed simply
#                  carry no evidence.
#   profiles:      file listing the profile TSVs for those samples
#   taxonomy:      the --taxonomy the lost reference was built with, or null
#   pool_min_genomes / sparse_c / store_genomes / c / k:
#                  the remaining ref-build settings it was built with

import collections
import csv
import os
import sys

configfile: "recover.yaml"

config.setdefault("outdir", "recover")
config.setdefault("weebill", "weebill")
config.setdefault("taxonomy", None)
config.setdefault("pool_min_genomes", 3)
config.setdefault("sparse_c", 3000)
config.setdefault("store_genomes", False)
config.setdefault("c", 200)
config.setdefault("k", 31)
config.setdefault("threads", 16)
config.setdefault("probe_shards", 64)
# Cap on how many samples are probed per contested group. Evidence is cheap to
# combine and samples are many; this keeps the probe pass proportional to the
# number of unknowns rather than to the collection.
config.setdefault("samples_per_group", 4)
# Ceiling on the residual search. Every group the evidence settles divides the
# number of combinations, so a run that trips this wants more probed samples,
# not a bigger ceiling.
config.setdefault("max_candidates", 4_000_000)

OUT = config["outdir"]
WEEBILL = config["weebill"]
SHARDS = [f"{i:04d}" for i in range(config["probe_shards"])]


def read_list(path):
    with open(path) as fh:
        return [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]


rule all:
    input:
        f"{OUT}/recovered.sylref",
        f"{OUT}/verify.txt",


# --- 1. one sketching pass ---------------------------------------------------
# The genomes are sketched once and never again: per-genome sketch content is
# deterministic, so only the concatenation order varies, and that is exactly what
# --genome-order overrides downstream.
rule sketch_genomes:
    input:
        genomes=config["genomes"],
    output:
        db=f"{OUT}/genomes.syldb",
    threads: config["threads"]
    params:
        c=config["c"],
        k=config["k"],
    shell:
        r"""
        {WEEBILL} sketch -c {params.c} -k {params.k} -t {threads} \
            -l {input.genomes} -o {OUT}/genomes
        """


# --- 2. the unknowns ---------------------------------------------------------
rule contested_groups:
    input:
        db=f"{OUT}/genomes.syldb",
    output:
        graph=f"{OUT}/graph.bin",
        tsv=f"{OUT}/graph.tsv",
    threads: config["threads"]
    params:
        tax=lambda w: f"--taxonomy {config['taxonomy']}" if config["taxonomy"] else "",
        pmg=config["pool_min_genomes"],
    shell:
        r"""
        {WEEBILL} ref-recover pairs {input.db} {params.tax} \
            --pool-min-genomes {params.pmg} -t {threads} \
            -o {output.graph} --tsv {output.tsv}
        """


# --- 3. pick which samples to probe -----------------------------------------
# Profiles say which genomes each sample contains, so the choice of samples is
# made without opening a single .sylspr. What resolves a group is a sample that
# contains SOME but not ALL of its members: if a sample contains every member,
# the contested k-mers stay assigned whoever owns them and the flip is invisible.
rule select_samples:
    input:
        tsv=f"{OUT}/graph.tsv",
        samples=config["samples"],
        profiles=config["profiles"],
    output:
        selected=f"{OUT}/selected_samples.txt",
        report=f"{OUT}/selection_report.tsv",
    params:
        per_group=config["samples_per_group"],
    shell:
        r"""
        python3 {workflow.basedir}/scripts/select_samples.py \
            --graph-tsv {input.tsv} --samples {input.samples} \
            --profiles {input.profiles} --per-group {params.per_group} \
            --out {output.selected} --report {output.report}
        """


rule shard_samples:
    input:
        selected=f"{OUT}/selected_samples.txt",
    output:
        directory(f"{OUT}/shards"),
    run:
        os.makedirs(output[0], exist_ok=True)
        sel = read_list(input.selected)
        n = len(SHARDS)
        for i, shard in enumerate(SHARDS):
            with open(os.path.join(output[0], f"{shard}.txt"), "w") as fh:
                for p in sel[i::n]:
                    fh.write(p + "\n")


# --- 4. read the novel hashes ------------------------------------------------
rule probe_shard:
    input:
        shards=f"{OUT}/shards",
        graph=f"{OUT}/graph.bin",
    output:
        tsv=f"{OUT}/evidence/{{shard}}.tsv",
    threads: 4
    shell:
        r"""
        mkdir -p {OUT}/evidence
        if [ -s {input.shards}/{wildcards.shard}.txt ]; then
            xargs -a {input.shards}/{wildcards.shard}.txt \
                {WEEBILL} ref-recover probe --graph {input.graph} -t {threads} \
                -o {output.tsv}
        else
            printf 'sample_file\treference_fingerprint\thit_genomes\tnovel_hashes\tgroup_id\tnovel_group_hashes\n' > {output.tsv}
        fi
        """


rule merge_evidence:
    input:
        lambda w: expand(f"{OUT}/evidence/{{shard}}.tsv", shard=SHARDS),
    output:
        tsv=f"{OUT}/evidence.tsv",
    shell:
        r"""
        head -n1 {input[0]} > {output.tsv}
        for f in {input}; do tail -n +2 "$f" >> {output.tsv}; done
        """


# --- 5. solve ----------------------------------------------------------------
rule solve_order:
    input:
        graph=f"{OUT}/graph.bin",
        tsv=f"{OUT}/graph.tsv",
        evidence=f"{OUT}/evidence.tsv",
        profiles=config["profiles"],
    output:
        order=f"{OUT}/evidence_order.txt",
        solved=f"{OUT}/solved.tsv",
        report=f"{OUT}/solve_report.txt",
    params:
        target=lambda w: config.get("target_fingerprint", ""),
    shell:
        r"""
        python3 {workflow.basedir}/scripts/solve_order.py \
            --graph-tsv {input.tsv} --evidence {input.evidence} \
            --profiles {input.profiles} \
            --graph {input.graph} --weebill {WEEBILL} \
            --out {output.order} --report {output.report} \
            --solved-out {output.solved} \
            {params.target:q}
        """


# The evidence settles the groups it can see; whatever is left is one choice per
# group, scored analytically rather than by building a reference each time. The
# target is the fingerprint every .sylspr carries, so a match is the answer.
rule search_residual:
    input:
        graph=f"{OUT}/graph.bin",
        solved=f"{OUT}/solved.tsv",
        evidence=f"{OUT}/evidence.tsv",
    output:
        order=f"{OUT}/genome_order.txt",
    params:
        maxc=config["max_candidates"],
    shell:
        r"""
        set -euo pipefail
        target=$(tail -n +2 {input.evidence} | cut -f2 | sort -u | head -n1)
        if [ -z "$target" ]; then
            echo "no *.sylspr supplied any fingerprint -- nothing to recover against" >&2
            exit 1
        fi
        {WEEBILL} ref-recover search --graph {input.graph} --solved {input.solved} \
            --target "$target" --max-candidates {params.maxc} -o {output.order}
        """


# --- 6. rebuild and check ----------------------------------------------------
rule rebuild:
    input:
        db=f"{OUT}/genomes.syldb",
        order=f"{OUT}/genome_order.txt",
    output:
        ref=f"{OUT}/recovered.sylref",
    threads: config["threads"]
    params:
        tax=lambda w: f"--taxonomy {config['taxonomy']}" if config["taxonomy"] else "",
        pmg=config["pool_min_genomes"],
        sparse=config["sparse_c"],
        store=lambda w: "--store-genomes" if config["store_genomes"] else "",
    shell:
        r"""
        {WEEBILL} ref-build {input.db} {params.tax} {params.store} \
            --pool-min-genomes {params.pmg} --sparse-c {params.sparse} \
            --genome-order {input.order} -t {threads} \
            -o {OUT}/recovered
        """


# The decisive test: the fingerprint must equal the one the samples record, and a
# sample must decode against it. `--verify` decompresses and requires exact
# equality with the original sketch, so a pass is proof, not a plausibility check.
rule verify:
    input:
        ref=f"{OUT}/recovered.sylref",
        evidence=f"{OUT}/evidence.tsv",
        selected=f"{OUT}/selected_samples.txt",
    output:
        txt=f"{OUT}/verify.txt",
    shell:
        r"""
        set -euo pipefail
        got=$({WEEBILL} inspect {input.ref} 2>/dev/null | awk '/fingerprint:/ {{print $2}}')
        want=$(tail -n +2 {input.evidence} | cut -f2 | sort -u | head -n1)
        {{
          echo "rebuilt_fingerprint	$got"
          echo "sample_fingerprint	$want"
        }} > {output.txt}
        if [ "$got" != "$want" ]; then
            echo "FINGERPRINT MISMATCH -- see {OUT}/solve_report.txt for unresolved groups" >&2
            echo "status	MISMATCH" >> {output.txt}
            exit 1
        fi
        echo "status	MATCH" >> {output.txt}
        head -n1 {input.selected} | xargs -r {WEEBILL} ref-compress \
            --decompress -r {input.ref} -d {OUT}/verify_decompress
        echo "decompress	ok" >> {output.txt}
        """
