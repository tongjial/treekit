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

## Input
Input files must be newick or nexus tree files, UTF-8 encoded and without a byte-order mark.
Each file may hold one tree or many trees, and `-i` accepts several files (of either format),
relative or absolute paths, and shell-style wildcards (for example `-i "trees/*.tre"`).
Taxon labels are case sensitive: `a` and `A` are two taxa.

## Commands
usage:
`treekit <command> [<args>]` or
`treekit <command> -h`
```
The options of commands are:
    summary     view trees summary
    taxa        view taxa summary
    remove      remove specific taxa/tips from trees
    keep        (planned) pick-up the trees with specific taxa
    mrca        (planned) check the monophyly or extract the clade based on given taxa/group
    extract     (planned) generate subtrees only consist of the given taxa
    prune       (planned) prune trees with specific branchlength
    root        (planned) root/reroot/unroot trees
    draw        (planned) display trees as textplot
    compare     (planned) compare trees between two different sources
Use treekit <command> -h for help with arguments of the command of interest
```

### Options shared by all commands
| option | meaning |
|---|---|
| `-i FILE [FILE ...]` | input tree files (required) |
| `-o FILE` | write the result to this file instead of stdout; it must not be an input file |
| `-c N` | use up to N processes, one per input file (default 1) |
| `-l FILE` | append run details and warnings to this log file |

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

## Examples
Sample data is in `examples/`:
```
treekit summary -i examples/test
treekit taxa -k -i examples/concat_6genes.trees
treekit remove -f examples/taxa.txt -i examples/concat_6genes.trees -l remove.log > reduced.tre
```

## Tests
```
pip install ".[test]"
pytest
```
