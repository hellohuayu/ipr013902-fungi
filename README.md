# IPR013902 taxonomic-distribution tree across fungal reference proteomes

Generates **"Taxonomic distribution of IPR013902 domain-containing proteins
mapped onto a fungal species tree"** from **"List of IPR013902
domain-containing proteins across 1,430 fungal reference proteomes"**.

The input is a TSV with one row per UniProt Reference Proteome; the output is
a text file that draws the NCBI taxonomy as a box-drawing tree, annotated
with the number of IPR013902 domain-containing proteins at every node.

## Input format

A tab-separated table with (at least) these columns:

| column | description |
|---|---|
| `Proteome Id` | UniProt proteome accession, e.g. `UP000002311` |
| `Organism` | organism name, may include strain information |
| `Taxonomic lineage` | NCBI lineage, comma-separated |
| `interpro_IPR013902_count` | number of IPR013902 (PIRSF-related) domain-containing proteins in that proteome |

Column names can be overridden with `--count-col`, `--lineage-col`,
`--organism-col`, `--proteome-col`.

## Usage

```bash
# default: every proteome is a leaf; strains of one species are listed
# as sibling leaves under their species node (not merged)
python generate_ipr_taxonomy_tree.py \
    --input fungal_ref_proteomes_with_interpro.tsv \
    --output IPR013902_taxonomy_tree.txt

# collapse proteomes of the same species into the species node
python generate_ipr_taxonomy_tree.py \
    --input fungal_ref_proteomes_with_interpro.tsv \
    --output IPR013902_taxonomy_tree_merged.txt \
    --merge-strains
```

Other options: `--root` (default `Fungi`), `--no-complete-species`.

## Output

```
─ Fungi (478 IPR013902 proteins)
  ├─ Dikarya (469 IPR013902 proteins in this lineage)
  │  ├─ Basidiomycota (371 IPR013902 proteins in this phylum)
  │  │  ├─ Agaricomycotina (369 IPR013902 proteins in this subphylum)
  ...
  │  │  │  │  │  ├─ Pichia membranifaciens (11 IPR013902 proteins in this species)
  │  │  │  │  │  │  ├─ Pichia membranifaciens NRRL Y-2026 [UP000094455] (5 IPR013902 proteins in this proteome)
  │  │  │  │  │  │  └─ Pichia membranifaciens [UP000186136] (6 IPR013902 proteins in this proteome)
```

Internal-node counts are exact sums of their descendant proteomes; the script
prints a self-check summary (leaf count, total, positives) on every run.

## Implementation notes

Compared with a naive lineage traversal, the script handles three edge cases
that occur in UniProt reference-proteome lineages:

1. **Strains are not merged by default.** Proteomes of the same species
   (species-level representative + strain proteomes, or several assemblies of
   one organism) are drawn as separate sibling leaves labelled with the
   organism name and proteome accession. Use `--merge-strains` to reproduce
   the collapsed style.
2. **Lineages ending above species level.** Some proteomes have lineages that
   stop at the genus (e.g. `..., Pichia`) even though the organism name is a
   resolved binomial. By default the script appends the binomial as a species
   node so that all proteomes of one species hang under the same node
   (disable with `--no-complete-species`). Unresolved organisms
   (`Penicillium sp.`, `Metschnikowia aff.`, species complexes, etc.) stay at
   their lineage position.
3. **Rank labels.** Internal nodes are labelled with standard botanical ranks
   inferred from suffixes (`-mycota` → phylum, `-ales` → order, `-aceae` →
   family, ...); binomials are labelled *species*; anything else is labelled
   *lineage*.

## Requirements

- Python >= 3.8
- pandas

```bash
pip install -r requirements.txt
```

## Citation / methods

If you use this script, please cite it via the GitHub link and describe it in
the methods as:

> The taxonomic distribution of IPR013902 domain-containing proteins across
> fungal reference proteomes was visualised on the NCBI taxonomy using a
> Python script (https://github.com/hellohuayu/ipr013902-fungi) that maps
> each reference proteome onto its UniProt taxonomic lineage; proteomes of the
> same species are displayed as separate leaves and node counts are exact
> sums of descendant proteomes.

## License

MIT
