#! /usr/bin/env python3

import argparse, sys, multiprocessing as mp
from pathlib import Path
from collections import Counter
import dendropy
import statistics
import time

class ParsedArgs:

    def __init__(self):

        # initialize the main ArgumentParser object
        parser = argparse.ArgumentParser(
            usage='''treekit <command> [<args>]
The treekit commands are:
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
'''
        )

        # add the main argument "command"
        parser.add_argument(
            "command", 
            help="Subcommand to run"
        )

        # parse the command by excluding the rest of args
        self.args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, self.args.command):
            print('Unrecognized command')
            parser.print_help()
            exit(1)
        
        # use dispatch pattern to invoke method for parsing the rest of args
        getattr(self, self.args.command)()

    def add_common_args(self, parser):

        # define required arguments for all commands
        requiredNamed = parser.add_argument_group(
            'required arguments'
        )
        
        # required arguments
        requiredNamed.add_argument(
            "-i",
            "--in-files",
            nargs = "+",
            dest = "in_files",
            type = str,
            required = True,
            help = "Tree files to be taken as input."
        )

        # given the outfile name, or output to stdout
        parser.add_argument(
            "-o",
            "--out-file",
            dest = "out_file",
            type = str,
            help = "give the outfile name."
        )

        # parallelization is used for file parsing and calculating stats
        parser.add_argument(
            "-c",
            "--cores",
            dest = "cores",
            default = 1,
            type = int,
            help = "Number of cores used. Default: 1"
        )

    def summary(self):

        # summarize the tree features
        parser = argparse.ArgumentParser(
            description="View trees summary",
        )

        # if give option '-e', then output the details
        parser.add_argument(
            "-e",
            "--details-out",
            dest = "details_out",
            action = "store_true",
            default = False,
            help = "Details of each tree summary."
        )

        # add shared arguments
        self.add_common_args(parser)
        args = parser.parse_args(sys.argv[2:])
        return args
    
    def taxa(self):

        # view the taxon namespace in various ways
        parser = argparse.ArgumentParser(
            description="View taxa summary",
        )

        # list common taxa names
        parser.add_argument(
            "-k",
            "--common-taxa",
            dest = "common_taxa",
            action = "store_true",
            default = False,
            help = "only list common taxa names in all trees."
        )

        # add shared arguments
        self.add_common_args(parser)
        args = parser.parse_args(sys.argv[2:])
        return args
    
    def extract(self):
        pass
    
    def remove(self):
        parser = argparse.ArgumentParser(
            description="remove specific taxa from given trees",
        )

        # two options for input taxa names
        group = parser.add_mutually_exclusive_group(
            required=True
        )

        group.add_argument(
            "-n",
            "--taxa-list",
            nargs = "+",
            dest = "taxa_list",
            type = str,
            help = "Input taxa names to be removed (with white-space separated).",
        )
        
        group.add_argument(
            "-f",
            "--taxa-file",
            dest = "taxa_file",
            type = str,
            help = "File with taxa to be removed (one taxa per line)."
        )

        parser.add_argument(
            "-u",
            "--out-format",
            dest = "out_format",
            choices = ["newick", "nexus"],
            default = "newick",
            help = "Format for output trees. Default: newick"
        )

        # add shared arguments
        self.add_common_args(parser)
        args = parser.parse_args(sys.argv[2:])
        return args
    
    def keep(self):
        pass

    def mrca(self):
        pass

    def prune(self):
        pass

    def root(self):
        pass

    def draw(self):
        pass

    def compare(self):
        pass

    def get_args_dict(self):

        # store arguments in a dictionary
        command = self.args.__dict__
        arguments = getattr(self, self.args.command)().__dict__
        argument_dictionary = command.copy()
        argument_dictionary.update(arguments)
        
        return argument_dictionary


class FuncTree:
    """Class for encapsulate treekit functions"""

    # code loading at the initialization phase

    def __init__(self, **kwargs):

        # set defaults and get values from kwargs
        self.command = kwargs.get("command")
        self.in_files = self.get_paths(kwargs.get("in_files"))
        self.in_formats = self.get_formats()
        
        # set the optional parameters
        if kwargs.get("out_file"):
            self.out_file = self.get_path(kwargs.get("out_file"))
        else: 
            self.out_file = None

        if kwargs.get("out_format"):
            self.out_format = kwargs.get("out_format")
        else: 
            self.out_format = None
        
        if kwargs.get("cores"):
            self.cores = self.get_cores(kwargs.get("cores"))
        else: 
            self.cores = 1

        # specific attributes for each command
        if self.command == "summary":
            self.details_out = kwargs.get("details_out")

        if self.command == "taxa":
            self.common_taxa = kwargs.get("common_taxa")
        
        if self.command == "remove":
            
            if kwargs.get("taxa_list"):
                self.taxa_list = kwargs.get("taxa_list")
            else:
                self.taxa_list = None
            
            if kwargs.get("taxa_file"):
                self.taxa_file = self.get_path(kwargs.get("taxa_file"))
            else:
                self.taxa_file = None

        """
        if self.command == "extract":
            pass
        
        if self.command == "remove":
            self.taxa_to_remove = kwargs.get("taxa_to_remove")
            self.taxa_file = kwargs.get("taxa_file")
            self.out_format = kwargs.get("out_format")
        """

    def get_paths(self, files):
        
        # each in_file may contain asterisk wildcard
        out_paths = []
        for file in files:
            paths = Path().glob(file)
            abs_paths = [path.resolve() for path in paths]
            out_paths.extend(abs_paths)
        return out_paths
    
    def get_path(self, file):

        # get outfile absolute path
        out_path = Path(file).resolve()
        return out_path
    
    def get_formats(self):

        # identify tree in nexus or newick
        in_formats=[]
        for in_file in self.in_files:
            with open(in_file, 'r') as f:
                while True:
                    line = f.readline().strip()
                    if line != '': break
            if line[0] == '#':
                in_formats.append("nexus")
            elif line[0] == '(':
                in_formats.append("newick")
            else:
                print("Please input the newick or nexus format tree files")
                sys.exit()
        return in_formats
    
    def get_cores(self, cores):

        # check the cores attribute
        assert type(cores) == int and cores > 0, \
            "-c(--cores) requires int and >1" 

        # if cores > real_cores, use real_cores
        real_cores = mp.cpu_count()
        if cores > real_cores:
            print("the current computer has up to " + \
                str(real_cores) + " cores, instead of " + \
                str(cores) + " cores given by you.\n" + \
                "as a better way, this job will be running in " + \
                str(real_cores) + " cores."
            )
            return real_cores
        else:
            return cores

    # code for Meta Functions on a tree object

    def is_binary(self, tree_object):

        # is the rooted/unrooted tree strictly bipartition 
        multi_nodes = []
        multifurcating = lambda x: True if len(x.child_nodes()) > 2 else False
        for node in tree_object.preorder_node_iter(multifurcating):
            multi_nodes.append(node)
        
        # none of multi nodes included in rooted tree 
        if self.is_rooted(tree_object):
            if len(multi_nodes) > 0:
                return False
            else:
                return True
        
        # at least one multi node included in unrooted tree
        if not self.is_rooted(tree_object):
            if len(multi_nodes) > 1:
                return False
            else:
                return True
    
    def is_rooted(self, tree_object):

        # the tree is rooted when seed node has two children
        if len(tree_object.seed_node.child_nodes()) == 2:
            return True
        
        # the tree is unrooted when seed node has three children
        elif len(tree_object.seed_node.child_nodes()) > 2: 
            return False

        # warning if failed to parse the seed node
        else:
            print("Invalid tree object, please check the import file")

    def has_branchlength(self, tree_object):

        # the tree with branch length
        if tree_object.length() > 0:
            return True
        
        # the tree without branch length
        elif tree_object.length() == 0:
            return False
        
        # warning if failed to parse the tree length
        else:
            print('Invalid tree object, please check the import file')

    def is_ultrametric(self, tree_object):

        # two conditions required in an ultrametric tree
        # rooted tree and the same tip heights
        if self.is_rooted(tree_object) and tree_object.length() > 0:
            tip_heights = tree_object.minmax_leaf_distance_from_root()
            min_tip_heights = decimal_correct(tip_heights[0])
            max_tip_heights = decimal_correct(tip_heights[1])
            if min_tip_heights == max_tip_heights:
                return True
            else:
                return False
        else:
            return False

    def get_leaf_heights(self, tree_object):

        # if the tree is rooted, calculate leaf heights
        if self.is_rooted(tree_object):
            leaf_heights = tree_object.calc_node_root_distances(
                return_leaf_distances_only=True
            )
            return leaf_heights
        
        # or return None
        else:
            return None

    # methods corresponding to summary block
    
    def get_summaries(self):
        
        # get summaries for each of tree objects
        summaries = []
        header = [
            "source",
            "index",
            "is_binary",
            "is_rooted",
            "has_brlen",
            "is_ultrametric",
            "tree_length",
            "avg_leaf_height",
            "leaf_height_variance",
            "num_tips"
        ]

        # generate the tree list object from each of files
        for in_file, in_format in zip(self.in_files, self.in_formats):
            tree_list = dendropy.TreeList.get(
                path = in_file, 
                schema = in_format
            )
            
            # invoke get_summs function to deal with each tree_list object
            summs = self.get_summs(tree_list, in_file.name)
            summaries.extend(summs)

        return header, summaries

    def get_summs(self, tree_list, infile_name):

        # initialize the properties
        summs = []
        source = []
        index = []
        is_binary = []
        is_rooted = []
        has_brlen = []
        is_ultrametric = []
        tree_length = []
        avg_leaf_height = []
        leaf_height_variance =[]
        num_tips = []

        # iterate tree in tree_list
        for tree_idx, tree_object in enumerate(tree_list):

            # get values for each tree object
            source.append(infile_name)
            index.append(tree_idx + 1)
            is_binary.append(self.is_binary(tree_object))
            is_rooted.append(self.is_rooted(tree_object))
            has_brlen.append(self.has_branchlength(tree_object))
            is_ultrametric.append(self.is_ultrametric(tree_object))
            tree_length.append(
                decimal_correct(tree_object.length())
            )
            num_tips.append(len(tree_object.leaf_nodes()))

            # if leaf heights list is not None, stat the leaf heights
            leaf_heights = self.get_leaf_heights(tree_object)
            if leaf_heights:
                avg_leaf_height.append(
                    decimal_correct(statistics.mean(leaf_heights))
                )
                leaf_height_variance.append(
                    decimal_correct(statistics.pstdev(leaf_heights))
                )
            else:
                avg_leaf_height.append('NA')
                leaf_height_variance.append('NA')

        # combine and format the results
        summs = list(zip(
            source,
            map(str, index),
            map(str, is_binary),
            map(str, is_rooted),
            map(str, has_brlen),
            map(str, is_ultrametric),
            map(str, tree_length),
            map(str, avg_leaf_height),
            map(str, leaf_height_variance),
            map(str, num_tips)
        ))

        return summs
    
    def get_summstats(self):

        # get summstat for all tree objects
        header = [
            "file_name",
            "num_tree",
            "num_binary",
            "num_rooted",
            "num_hasbrlen",
            "num_ultrametric",
            "avg_tree_length",
            "tree_length_variance",
            "num_taxa",
            "num_common_taxa"
        ]

        # init the property lists of all input files
        summstats = []
        file_names = []
        num_trees = []
        num_binary = []
        num_rooted = []
        num_hasbrlen = []
        num_ultrametric = []
        avg_tree_length = []
        tree_length_variance = []
        num_taxa = []
        num_common_taxa = []

        # pass the trees to TreeList object
        for in_file, in_format in zip(self.in_files, self.in_formats):
            tree_list = dendropy.TreeList.get(
                path = in_file, 
                schema = in_format
            )

            # get the values of each file
            file_names.append(in_file.name)
            num_trees.append(len(tree_list))
            num_binary.append(self.get_num_binary(tree_list))
            num_rooted.append(self.get_num_rooted(tree_list))
            num_hasbrlen.append(self.get_num_hasbrlen(tree_list))
            num_ultrametric.append(self.get_num_ultrametric(tree_list))
            avg_tree_length.append(self.get_avg_tree_length(tree_list))
            tree_length_variance.append(self.get_tree_length_variance(tree_list))
            num_taxa.append(len(tree_list.taxon_namespace))
            num_common_taxa.append(len(self.get_common_taxa(tree_list)))

        # combine and format the results
        summstats = list(zip(
            file_names, 
            map(str, num_trees), 
            map(str, num_binary),
            map(str, num_rooted),
            map(str, num_hasbrlen),
            map(str, num_ultrametric),
            map(str, avg_tree_length),
            map(str, tree_length_variance),
            map(str, num_taxa),
            map(str, num_common_taxa)
        ))

        return header, summstats

    def get_num_binary(self, tree_list):

        # get number of binary trees in given treelist
        if self.cores == 1:
            binary = [self.is_binary(tree_object) for tree_object in tree_list]
            return Counter(binary)[True]

        elif self.cores > 1:
            pool = mp.Pool(int(self.cores))
            binary = pool.map(self.is_binary, tree_list)
            return Counter(binary)[True]
    
    def get_num_rooted(self, tree_list):

        # get number of rooted trees in given treelist
        if self.cores == 1:
            rooted = [self.is_rooted(tree_object) for tree_object in tree_list]
            return Counter(rooted)[True]

        elif self.cores > 1:
            pool = mp.Pool(int(self.cores))
            rooted = pool.map(self.is_rooted, tree_list)
            return Counter(rooted)[True]

    def get_num_hasbrlen(self, tree_list):

        # get the number of trees with branch length
        if self.cores == 1:
            hasbrlen = [self.has_branchlength(tree_object) for tree_object in tree_list]
            return Counter(hasbrlen)[True]

        elif self.cores > 1:
            pool = mp.Pool(int(self.cores))
            hasbrlen = pool.map(self.has_branchlength, tree_list)
            return Counter(hasbrlen)[True]        

    def get_num_ultrametric(self, tree_list):

        # get the number of trees are ultrametric
        if self.cores == 1:
            ultrametric = [self.is_ultrametric(tree_object) for tree_object in tree_list]
            return Counter(ultrametric)[True]

        elif self.cores > 1:
            pool = mp.Pool(int(self.cores))
            ultrametric = pool.map(self.is_ultrametric, tree_list)
            return Counter(ultrametric)[True]

    def get_avg_tree_length(self, tree_list):

        # calculate the average of tree lengths
        if self.cores == 1:
            tree_lengths = [tree_object.length() for tree_object in tree_list]
            return decimal_correct(statistics.mean(tree_lengths))

        elif self.cores > 1:
            pool = mp.Pool(int(self.cores))
            tree_lengths = pool.map(dendropy.Tree.length, tree_list)
            return decimal_correct(statistics.mean(tree_lengths))

    def get_tree_length_variance(self, tree_list):

        # calculate the pstdev of tree lengths
        if self.cores == 1:
            tree_lengths = [tree_object.length() for tree_object in tree_list]
            return decimal_correct(statistics.pstdev(tree_lengths))

        elif self.cores > 1:
            pool = mp.Pool(int(self.cores))
            tree_lengths = pool.map(dendropy.Tree.length, tree_list)
            return decimal_correct(statistics.pstdev(tree_lengths))

    def write_summaries(self):

        # get the detailed summaries
        if self.details_out:
            header, summaries = self.get_summaries()
            header_out = '\t'.join(header)
            summ_out = ['\t'.join(summ) for summ in summaries]

            # if outfile is given, then write to outfile
            if self.out_file != None:
                self.file_overwrite_warning()
                with open(self.out_file, "w") as summ_file:
                    summ_file.write(header_out + '\n')
                    summ_file.write('\n'.join(summ_out))
                print("Wrote summaries to file '" + str(self.out_file) + "'")

            # otherwise print to stdout
            else:
                print(header_out + '\n' + '\n'.join(summ_out))
        
        # Or get the brief summaries
        else:
            header, summstats = self.get_summstats()
            header_out = '\t'.join(header)
            summ_out = ['\t'.join(summ) for summ in summstats]

            # if outfile is given, then write to outfile
            if self.out_file != None:
                self.file_overwrite_warning()
                with open(self.out_file, "w") as summ_file:
                    summ_file.write(header_out + '\n')
                    summ_file.write('\n'.join(summ_out))
                print("Wrote summaries to file '" + str(self.out_file) + "'")
            
            # otherwise print to stdout
            else:
                print(header_out + '\n' + '\n'.join(summ_out))

    # methods corresponding to taxa code block

    def get_taxa_list(self):
        
        # initial treelist, Namespace and out taxalist
        taxa_space = dendropy.TaxonNamespace()
        tree_list = dendropy.TreeList(taxon_namespace = taxa_space)
        taxa_list = []

        # generate tree objects for all input tree files
        for in_file, in_format in zip(self.in_files, self.in_formats):
            tree_list.read(
                path = in_file, 
                schema = in_format,
                preserve_underscores=True
            )

        # if common names are required, get the intersections of labels
        if self.common_taxa:
            taxa_list = self.get_common_taxa(tree_list)

        # if required all names (default), get namespace labels
        else:
            taxa_list = taxa_space.labels()

        # sort name list and return
        taxa_list.sort()
        return taxa_list

    def get_common_taxa(self, tree_list):

        # get taxa in the first tree
        taxa_tree0 = set(
            leaf.taxon.label for leaf in tree_list[0].leaf_nodes()
        )

        # get taxa list in all trees
        taxa_trees = [
            [leaf.taxon.label for leaf in tree.leaf_nodes()] 
            for tree in tree_list
        ]

        # calculate the intersection between taxa of trees
        taxa_list = list(taxa_tree0.intersection(*taxa_trees))

        return taxa_list

    def write_taxa_out(self):
        
        # get taxa for all trees or get the common taxa names
        taxa_out = self.get_taxa_list()

        # if outfile is given, then write to outfile
        if self.out_file != None:
            self.file_overwrite_warning()
            with open(self.out_file, "w") as taxa_file:
                taxa_file.write('\n'.join(taxa_out))
            print("Wrote taxa list to the file '" + str(self.out_file) + "'")

        # otherwise print to stdout
        else:
            print('\n'.join(taxa_out))

    # methods corresponding to remove code block

    def get_reduced_trees(self):
        
        # initial treelist, namespace
        taxa_space = dendropy.TaxonNamespace()
        tree_list = dendropy.TreeList(taxon_namespace = taxa_space)
        
        # generate tree list object
        for in_file, in_format in zip(self.in_files, self.in_formats):
            tree_list.read(
                path = in_file, 
                schema = in_format,
                preserve_underscores=True
            )

        # parse the list for taxa to be removed 
        if self.taxa_file:
            with open(self.taxa_file, "r") as taxa_file:
                taxa_list = [line.strip(" \r\n") for line in taxa_file.readlines()]
                
        elif self.taxa_list:
            taxa_list = self.taxa_list
            
        # prune the taxa base on labels
        for tree in tree_list:
            tree.prune_taxa_with_labels(taxa_list)

        return tree_list

    def write_reduced_trees(self):

        tree_list = self.get_reduced_trees()
        tree_out = tree_list.as_string(schema=self.out_format).replace("'","")
        
        # if out_file is given
        if self.out_file != None:
            self.file_overwrite_warning()
            with open(self.out_file, 'w') as out_file:
                out_file.write(tree_out)
        
        # or output to stdout
        else:
            print(tree_out)

    # general methods for all common blocks

    def file_overwrite_warning(self):

        # print warning when overwriting a file
        if self.out_file.exists():
            print("WARNING: You are overwriting '" + str(self.out_file) + "'")

# general functions for global invoking

def decimal_correct(m, n=9):
    
    # keep n decimal places behind the point for input float m
    return int(float(m) * (10 ** n)) / (10 ** n)

# main function

def main():
    
    # initialize parsed arguments and FuncTree objects
    kwargs = configs()
    func_tree = FuncTree(**kwargs)

    # run FuncTree methods to get outputs   
    if func_tree.command == "summary":
        func_tree.write_summaries()

    if func_tree.command == "taxa":
        func_tree.write_taxa_out()

    if func_tree.command == "remove":
        func_tree.write_reduced_trees()
    
def configs():

    # initialize parsed arguments
    config = ParsedArgs()

    # get arguments
    config_dict = config.get_args_dict()
    return config_dict
    
if __name__ == '__main__':

    # start = time.time()
    main()
    # end = time.time()
    # print('time consume:', end - start)
