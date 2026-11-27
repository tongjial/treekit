# TreeKit
A toolkit for the mass trees manipulation

## Dependency
* [python3](https://www.python.org/downloads/) 3.9 or newer (tested with 3.12)
* [DendroPy](https://dendropy.org/) 5.0 or newer

## Installation
Clone the repository and install it with pip (this provides the `treekit` command):
```
git clone https://github.com/tongjial/treekit.git
cd treekit
pip install .
```
Or use `treekit.py` directly as a stand-alone script with `python3 treekit.py <command> ...`.

### Download an executable (no Python needed)
Every [release](https://github.com/tongjial/treekit/releases) has a single-file executable for
Linux (x86_64, glibc 2.17 or newer), macOS (Apple Silicon and Intel) and Windows (x86_64). Each archive
holds `treekit` (`treekit.exe` on Windows), the license, this README and the notices of the bundled software.
```
tar -xzf treekit-v0.1.0-linux-x86_64.tar.gz        # Windows: unzip the .zip
./treekit-v0.1.0-linux-x86_64/treekit --version
```
Compare the checksum with the `.sha256` file of the release (`sha256sum -c`, `shasum -a 256 -c` on
macOS, or `Get-FileHash` in PowerShell). The executables are not signed: on macOS run
`xattr -d com.apple.quarantine treekit` once, on Windows choose "More info" and "Run anyway" in SmartScreen.
The executable unpacks itself on every start, which takes about 0.2 s.

## Input
Input files must be newick or nexus tree files, UTF-8 encoded and without a byte-order mark.
Each file may hold one tree or many trees, and `-i` accepts several files (of either format),
relative or absolute paths, and shell-style wildcards (for example `-i "trees/*.tre"`).
Taxon labels are case sensitive: `a` and `A` are two taxa. Labels with spaces are written quoted
(`'A b'`). Comments such as `[&posterior=0.98,note="red, blue"]` are kept in the trees that commands
write, exactly as written. A file may hold a single tip with comments (`A[&state=blue];`).

## Commands
usage:
`treekit <command> [<args>]` or
`treekit <command> -h`
```
The options of commands are:
    summary     view trees summary
    taxa        view taxa summary
    remove      remove specific taxa/tips from trees
    keep        pick-up the trees with specific taxa
    mrca        check the monophyly or extract the clade based on given taxa/group
    extract     generate subtrees only consist of the given taxa
    prune       prune trees with specific branchlength
    root        root/reroot/unroot trees
    draw        display trees as textplot
    compare     compare trees between two different sources
    collapse    collapse branches with low support into polytomies
Use treekit <command> -h for help with arguments of the command of interest
```

### Options shared by all commands
| option | meaning |
|---|---|
| `-i FILE [FILE ...]` | input tree files (required) |
| `-o FILE` | write the result to this file instead of stdout; it must not be a file the command reads (`-i`, `-r`, `-f`, `-g`) |
| `-c N` | use up to N processes, one per input file, for `summary` and `compare` (default 1) |
| `-l FILE` | append run details and warnings to this log file; not a file the command reads, nor the `-o` file |

Results go to stdout; warnings and errors go to stderr, so the output can be piped safely.
The exit status is 0 on success, 1 on an error (bad input, nothing to do) and 2 on a usage error.

### summary
`treekit summary -i trees.tre` prints one row per input file; `-e` prints one row per tree.

| column | meaning |
|---|---|
| `num_binary`, `is_binary` | strictly bifurcating (a trifurcating root counts for an unrooted tree) |
| `num_rooted`, `is_rooted` | the root has exactly two children |
| `num_hasbrlen`, `has_brlen` | the total tree length is positive |
| `num_ultrametric`, `is_ultrametric` | rooted, with all tips at the same distance from the root (within 1e-9) |
| `avg_tree_length`, `tree_length_sd` | mean and population standard deviation of the tree lengths |
| `avg_leaf_height`, `leaf_height_sd` | mean and population standard deviation of the tip heights; `NA` for unrooted trees |
| `num_taxa`, `num_common_taxa` | taxa in the file, and taxa present in every tree of the file |

### taxa
`treekit taxa -i trees.tre` lists all taxa, sorted; `-k` lists only the taxa present in every tree.

### remove
`treekit remove -n A B -i trees.tre` removes the taxa A and B from every tree
(or `-f names.txt` with one name per line). `-u nexus` writes nexus instead of newick.

* A taxon that is absent from some trees is normal for sets of gene trees; it is only noted in the log.
* A taxon absent from **all** trees produces a warning that lists them. If none of the requested
  taxa is found, the run stops with an error (exit status 1) and writes no trees.
* `--strict` turns the warning for absent taxa into an error.
* A warning is given when trees are left with a single tip; the run stops if a tree would be left empty.
* Trees that have tips without a label are rejected, because pruning would delete those tips too.

### extract
`treekit extract -n A B C -i trees.tre` keeps only the taxa A, B and C in every tree
(or `-f names.txt`); the other tips are removed and branch lengths are merged. `-u nexus`
writes nexus.

* Taxa missing from some trees are noted in the log; taxa missing from all trees give one
  warning (an error with `--strict`), and the run stops if none of them is found.
* A tree that contains none of the requested taxa stops the run: select the trees first
  with `treekit keep`. A warning is given for trees left with a single tip.

### keep
`treekit keep -n A B -i trees.tre` writes, unchanged and in input order, only the trees that
contain **all** the taxa A and B (or `-f names.txt`). The log lists the kept trees. The run
stops with an error if no tree contains all of them.

### root
Exactly one of:

* `treekit root -n O1 O2 -i trees.tre` roots every tree on the outgroup O1, O2 (or `-f names.txt`),
  in the middle of the outgroup's branch. A tree that lacks some outgroup taxa is rooted on the
  ones it has. A tree already rooted between the outgroup and the rest is kept as it is, so time
  trees stay ultrametric.
* `treekit root --midpoint -i trees.tre` roots every tree halfway along its longest tip-to-tip path
  (needs branch lengths). A tree that is already rooted at its midpoint, such as a time tree, is
  kept as it is, with everything on its nodes.
* `treekit root --unroot -i trees.tre` removes the root, merging the two root branches.

A tree that cannot be rooted (no outgroup taxon, an outgroup that is not monophyletic, no branch
lengths for `--midpoint`) is written unchanged, so the output has the same trees in the same order
as the input; one warning per reason lists the trees. `--strict` makes this an error.
When the root moves, what describes a branch stays on that branch: node labels and the items of the
comments such as `posterior`, `rate` or `length` follow the split. What describes a node itself, the items
`height`, `height_*`, `age` and `age_*` (BEAST, MrBayes), is only valid for the tips below that node: it
stays on a node whose tips are unchanged and is removed from the others, never moved to another node
(a warning lists the trees). The two branches at an old root become one branch: their labels and items
are combined, and when both have a label or an item of the same name with different values, the first
one is kept (a warning lists the trees).

### mrca
`treekit mrca -n A B C -i trees.tre` tests whether the group A, B, C is monophyletic in every
tree; `-g groups.tsv` tests several groups at once (one `taxon<TAB>group` per line, `#` starts a
comment). The table has one row per tree and group:

| column | meaning |
|---|---|
| `source`, `index` | input file and tree number in that file |
| `rooted` | whether the rooted rule was used (the tree has a bifurcating root) |
| `members_present`, `members_total` | group members in this tree, and in the group |
| `monophyletic` | `True`/`False`; `NA` if fewer than two members are present, or no other tip |
| `num_intruders`, `intruders` | for rooted trees, the other tips inside the clade of the group's MRCA |

* Rooted trees: the group is monophyletic if the clade of its most recent common ancestor holds
  only group members. Unrooted trees: if one branch separates the group from all other tips,
  wherever the root is placed. `--unrooted` applies the unrooted rule to every tree.
* Groups are tested on the members present in each tree; the log has the count per group.
* `--extract` writes, instead of the table, the clade of every group in every tree (for rooted
  trees the MRCA clade, intruders included). In nexus output (`-u nexus`) each clade is named
  `<group>_tree<number>`. Clades that cannot be extracted are skipped with a warning
  (an error with `--strict`).

### draw
`treekit draw -i trees.tre` draws every tree as text, each under a title with its number, file
and index in the file. `-t 3 7` draws only trees 3 and 7 (numbers count over all input files),
`-w 120` sets the width (default 80) and `-s` shows internal node labels such as supports.
`-m length` scales the branches, `-m topology` ignores their lengths; the default `auto`
scales them unless lengths are missing or so uneven that most branches would be shorter than
one character (noted in the log).

### compare
`treekit compare -i a.tre -r b.tre` compares every tree of `-i` with every tree of `-r`
(both accept several files) and writes one row per pair:

| column | meaning |
|---|---|
| `source_a`, `index_a`, `source_b`, `index_b` | the two trees |
| `common_taxa`, `only_a`, `only_b` | taxa in both trees, and in only one of them |
| `rf` | Robinson-Foulds distance: branches (splits) found in only one of the two trees |
| `norm_rf` | `rf` divided by the number of branches of both trees (0 = same, 1 = no branch shared) |
| `weighted_rf` | sum of the branch length differences over all branches; `NA` if a branch length is missing |

* Trees with different taxa are both reduced to their common taxa before the comparison.
  Pairs with fewer than 4 common taxa (3 with `--rooted`) have no testable branch: `rf` is 0
  and `norm_rf` is `NA`, with a warning.
* Trees are compared as unrooted trees by default; `--rooted` compares clades instead.
* The number of rows is the product of the two tree counts (846 x 846 trees: about 0.7 million
  rows in about 30 s).

### prune
`treekit prune --max-ratio 10 -i trees.tre` removes the tips whose branch is longer than 10 times
the median branch length of their tree; `--max-length 0.5` removes the tips whose branch is longer
than 0.5. With both, a tip longer than either limit is removed. Only tip branches are tested, in one
pass, on the branch lengths before pruning; the log lists every removed tip.

* A tree without branch lengths (or, with `--max-ratio` alone, with a median branch length of 0)
  is written unchanged with a warning; `--strict` makes this an error.
* The run stops if every tip of a tree would be removed, and warns about trees left with one tip.

### collapse
`treekit collapse -m 70 -i trees.tre` collapses every internal branch whose support is lower than
70 into a polytomy. The length of a collapsed branch is added to the branches below it, so
root-to-tip distances (and ultrametric trees) are kept. Tips and the root are never collapsed.

* Supports are read from the node labels, as written by RAxML or IQ-TREE (`)95:0.1`). For labels
  like `85.3/97` (SH-aLRT/UFBoot), `--support-field 2` uses the second number.
* `--support-key posterior` reads them from annotations such as `[&posterior=0.98]` (BEAST,
  MrBayes); the annotations are kept in the output.
* Branches without a readable support are left as they are. A tree without any support is
  written unchanged with a warning (an error with `--strict`).
* A threshold on the wrong scale gets a warning (e.g. `-m 0.7` with supports up to 100).

## Examples
Sample data is in `examples/`:
```
treekit summary -i examples/test
treekit taxa -k -i examples/concat_6genes.trees
treekit remove -f examples/taxa.txt -i examples/concat_6genes.trees -l remove.log > reduced.tre
treekit keep -n Badula_crassa_Y24864 Coris_monspeliensis_Y925 -i examples/concat_6genes.trees > kept.tre
treekit extract -f examples/taxa.txt -i examples/concat_6genes.trees.all > pair.tre
treekit root --midpoint -i examples/RAxML_bipartitions.all > rooted.tre
treekit mrca -n Badula_crassa_Y24864 Coris_monspeliensis_Y925 -i examples/concat_6genes.trees
treekit draw -t 1 -w 100 -i examples/RAxML_bipartitions.all
treekit compare -i examples/concat_6genes.trees -r examples/concat_6genes.mcc.tre
treekit prune --max-ratio 12 -i examples/RAxML_bipartitions.all -l prune.log > pruned.tre
treekit collapse -m 70 -i examples/RAxML_bipartitions.all > collapsed.tre
treekit collapse -m 0.95 --support-key posterior -i examples/concat_6genes.mcc.tre > collapsed_mcc.tre
```

## Tests
```
pip install ".[test]"
pytest
```
