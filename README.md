# TreeKit
A toolkit for the mass trees manipulation

## Dependency
* [DendroPy](https://dendropy.org/)

## Current implementation
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
Use treekit <command> -h for help with arguments of the command of interest
```
