#!/usr/bin/env python3
"""
Generate a taxonomic-distribution tree of IPR013902 domain-containing proteins
across fungal reference proteomes.

Input
-----
A tab-separated table (TSV) with one row per UniProt Reference Proteome
("List of IPR013902 domain-containing proteins across N fungal reference
proteomes"), containing at least:

  - Proteome Id          : UniProt proteome accession (e.g. UP000002311)
  - Organism             : organism name (may include strain information)
  - Taxonomic lineage    : NCBI lineage, comma-separated
                           (e.g. "cellular organisms, Eukaryota, ..., Fungi, ...")
  - interpro_IPR013902_count : number of IPR013902 domain-containing proteins
                           detected in that proteome

Column names can be overridden with command-line options.

Output
------
A text file drawing the taxonomy tree with box-drawing characters, in the
style of "Taxonomic distribution of IPR013902 domain-containing proteins
mapped onto a fungal species tree". Every proteome is a leaf; internal-node
counts are the exact sums of their descendant proteomes (no double counting).

By default, proteomes of the same species (different strains / assemblies)
are NOT merged: they are listed as sibling leaves under their species node.
Use --merge-strains to collapse them into the species node instead.

Usage
-----
python generate_ipr_taxonomy_tree.py \
    --input fungal_ref_proteomes_with_interpro.tsv \
    --output IPR013902_taxonomy_tree.txt

Options allow a different root (default: Fungi), a different count column,
and strain merging.

Requires: pandas
"""

import argparse
import re
import sys
from collections import OrderedDict

import pandas as pd


# --------------------------------------------------------------------------
# Rank assignment for internal nodes
# --------------------------------------------------------------------------

# Suffix-based heuristics cover the standard botanical ranks used in NCBI
# fungal lineages. Anything not matched falls back to "lineage".
SUFFIX_RANKS = [
    ("subgen.", "subgenus"),
    ("mycotina", "subphylum"),
    ("mycota", "phylum"),
    ("mycetidae", "subclass"),
    ("mycetes", "class"),
    ("ineae", "suborder"),
    ("ales", "order"),
    ("aceae", "family"),
    ("oideae", "subfamily"),
]

UNRESOLVED_WORDS = {"sp.", "sp", "aff.", "cf.", "var.", "subsp.", "sensu", "nov."}
BINOMIAL_RE = re.compile(r"^[A-Z][a-z]+ [a-z]+$")


def assign_rank(name, has_binomial_children):
    """Return a rank label for an internal tree node."""
    if "incertae sedis" in name:
        return "group"
    if name.endswith("species complex") or name.endswith("species group"):
        return "species complex"
    if BINOMIAL_RE.match(name):
        return "species"
    lower = name.lower()
    for suffix, rank in SUFFIX_RANKS:
        if lower.endswith(suffix):
            return rank
    if "clade" in lower:
        return "clade"
    if has_binomial_children:
        return "genus"
    return "lineage"


# --------------------------------------------------------------------------
# Tree construction
# --------------------------------------------------------------------------

class Node:
    __slots__ = ("children", "proteomes", "count")

    def __init__(self):
        self.children = OrderedDict()   # name -> Node
        self.proteomes = []             # list of (leaf_label, count)
        self.count = 0


def species_binomial(organism_name):
    """Return 'Genus epithet' from an organism name, or None if unresolved."""
    clean = re.sub(r"\s*\(.*$", "", organism_name).strip()
    words = clean.split()
    if len(words) >= 2 and words[1].islower() and words[1] not in UNRESOLVED_WORDS:
        return f"{words[0]} {words[1]}"
    return None


def build_tree(df, lineage_col, count_col, organism_col, proteome_col,
               root_taxon, complete_species=True):
    """
    Build the taxonomic tree.

    Each row of `df` becomes one proteome leaf attached to the last element
    of its lineage. If the lineage ends above species level but the organism
    name is a resolved binomial (and --complete-species), a species node is
    appended so that all proteomes of one species share the species node.
    """
    root = Node()
    for _, row in df.iterrows():
        lineage = [p.strip() for p in str(row[lineage_col]).split(",")]
        if root_taxon not in lineage:
            sys.exit(f"error: root taxon '{root_taxon}' not found in lineage of "
                     f"proteome {row[proteome_col]}")
        path = lineage[lineage.index(root_taxon):]

        if complete_species:
            binom = species_binomial(str(row[organism_col]))
            if (binom and path[-1] != binom and not path[-1].startswith(binom)
                    and path[-1].split()[0] == binom.split()[0]):
                path = path + [binom]

        node = root
        for name in path[1:]:  # skip the root taxon itself
            node = node.children.setdefault(name, Node())
        leaf = f"{row[organism_col]} [{row[proteome_col]}]"
        node.proteomes.append((leaf, int(row[count_col])))

    def compute(node):
        total = sum(c for _, c in node.proteomes)
        for child in node.children.values():
            total += compute(child)
        node.count = total
        return total

    root.count = compute(root)
    return root


# --------------------------------------------------------------------------
# Rendering
# --------------------------------------------------------------------------

def _plural(n):
    return "protein" if n == 1 else "proteins"


def render(node, name, prefix, is_last, is_root, rank_override=None):
    """Yield lines of the tree drawing (depth-first, pre-order)."""
    if is_root:
        yield f"\u2500 {name} ({node.count} IPR013902 {_plural(node.count)})"
        child_prefix = "  "
    else:
        branch = "\u2514\u2500 " if is_last else "\u251c\u2500 "
        rank = rank_override if rank_override else assign_rank(
            name, any(BINOMIAL_RE.match(c) for c in node.children))
        yield (f"{prefix}{branch}{name} ({node.count} IPR013902 "
               f"{_plural(node.count)} in this {rank})")
        child_prefix = prefix + ("   " if is_last else "\u2502  ")

    # taxon children first (largest count first), then proteome leaves
    taxa = sorted(node.children.items(), key=lambda kv: (-kv[1].count, kv[0]))
    items = [(cname, child, None) for cname, child in taxa]
    items += [(None, None, (label, cnt)) for label, cnt in node.proteomes]

    n_items = len(items)
    for idx, (cname, child, leaf) in enumerate(items):
        last = idx == n_items - 1
        if cname is not None:
            yield from render(child, cname, child_prefix, last, False)
        else:
            label, cnt = leaf
            branch = "\u2514\u2500 " if last else "\u251c\u2500 "
            yield (f"{child_prefix}{branch}{label} ({cnt} IPR013902 "
                   f"{_plural(cnt)} in this proteome)")


def merge_strains_into_species(root):
    """
    Collapse proteome leaves into their parent nodes (used only with
    --merge-strains): every node keeps its summed count and the individual
    proteome lines are dropped, reproducing the style of the original
    supplementary tree.
    """
    def walk(node):
        node.proteomes = []
        for child in node.children.values():
            walk(child)
    walk(root)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Build the taxonomic-distribution tree of IPR013902 "
                    "domain-containing proteins across fungal reference proteomes.")
    ap.add_argument("--input", required=True,
                    help="input TSV, one row per reference proteome")
    ap.add_argument("--output", required=True, help="output tree (text)")
    ap.add_argument("--count-col", default="interpro_IPR013902_count",
                    help="column with IPR013902 protein counts per proteome")
    ap.add_argument("--lineage-col", default="Taxonomic lineage",
                    help="column with the NCBI lineage")
    ap.add_argument("--organism-col", default="Organism",
                    help="column with the organism name")
    ap.add_argument("--proteome-col", default="Proteome Id",
                    help="column with the proteome accession")
    ap.add_argument("--root", default="Fungi",
                    help="taxon to root the tree at (default: Fungi)")
    ap.add_argument("--merge-strains", action="store_true",
                    help="collapse proteomes of the same species into the "
                         "species node (default: list them as sibling leaves)")
    ap.add_argument("--no-complete-species", action="store_true",
                    help="do not create species nodes for lineages that end "
                         "above species level")
    args = ap.parse_args(argv)

    df = pd.read_csv(args.input, sep="\t")
    for col in (args.count_col, args.lineage_col, args.organism_col, args.proteome_col):
        if col not in df.columns:
            sys.exit(f"error: required column '{col}' not found in {args.input}")

    root = build_tree(df, args.lineage_col, args.count_col,
                      args.organism_col, args.proteome_col,
                      args.root, complete_species=not args.no_complete_species)

    # ---- summary / self-check (computed before any optional merging) ----
    leaf_counts = []

    def collect(node):
        for _, c in node.proteomes:
            leaf_counts.append(c)
        for child in node.children.values():
            collect(child)

    collect(root)
    arr = pd.Series(leaf_counts)

    if args.merge_strains:
        merge_strains_into_species(root)

    header = (
        "Taxonomy tree illustrating the distribution of IPR013902 domain proteins:\n"
        f"(rooted at {args.root}; internal counts are exact sums of descendant "
        "proteomes)\n"
    )
    lines = [header]
    lines.extend(render(root, args.root, "", True, True))

    with open(args.output, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

    print(f"proteomes (leaves)          : {len(arr)}")
    print(f"total IPR013902 proteins    : {arr.sum()}")
    print(f"proteomes with >=1 protein  : {(arr > 0).sum()}")
    print(f"proteomes with >1 protein   : {(arr > 1).sum()}")
    print(f"proteomes with >5 proteins  : {(arr > 5).sum()}")
    print(f"max proteins in one proteome: {arr.max()}")
    print(f"written: {args.output}")


if __name__ == "__main__":
    main()
