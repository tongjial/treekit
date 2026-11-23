#! /usr/bin/env python3

import argparse, sys, multiprocessing as mp
import glob, logging, math, os, re
from collections import Counter
from pathlib import Path
import dendropy
from dendropy.utility.error import DataParseError
import statistics

# every subcommand: name, one-line help, and the FuncTree method that runs it;
# the order is the order of the usage text. A method of None means "planned".
COMMANDS = (
    ("summary", "view trees summary", "write_summaries"),
    ("taxa", "view taxa summary", "write_taxa_out"),
    ("remove", "remove specific taxa/tips from trees", "write_reduced_trees"),
    ("keep", "pick-up the trees with specific taxa", "write_kept_trees"),
    ("mrca", "check the monophyly or extract the clade based on given taxa/group", "write_mrca"),
    ("extract", "generate subtrees only consist of the given taxa", "write_extracted_trees"),
    ("prune", "prune trees with specific branchlength", "write_pruned_trees"),
    ("root", "root/reroot/unroot trees", "write_rooted_trees"),
    ("draw", "display trees as textplot", "write_drawings"),
    ("compare", "compare trees between two different sources", "write_comparison"),
    ("collapse", "collapse branches with low support into polytomies", "write_collapsed_trees"),
)
COMMAND_RUNNERS = {name: method for name, _, method in COMMANDS if method}
IMPLEMENTED_COMMANDS = tuple(name for name, _, method in COMMANDS if method)
PLANNED_COMMANDS = tuple(name for name, _, method in COMMANDS if not method)

def usage_text():

    # the usage text lists every command, marking the ones not implemented yet
    lines = ["treekit <command> [<args>]", "The treekit commands are:"]
    for name, help_text, method in COMMANDS:
        prefix = "" if method else "(planned) "
        lines.append("    %-12s%s%s" % (name, prefix, help_text))
    lines.append("Use treekit <command> -h for help with arguments of the command of interest")
    return "\n".join(lines) + "\n"

class TreekitError(Exception):
    """User-facing error: reported on stderr with exit status 1."""

# relative/absolute tolerance for comparing leaf heights (ultrametric check)
ULTRAMETRIC_TOL = 1e-9

# a newick tree made of a single leaf, such as "A;" or "'A b':1.5;"
BARE_NEWICK = re.compile(
    r"^('([^']|'')*'|[^\s()\[\]:;,'#]+)(:[-+0-9.eE]+)?\s*;"
)

# warnings and errors go to stderr; with --log, info lines also go to a log file
log = logging.getLogger("treekit")
logging.addLevelName(logging.WARNING, "warning")
logging.addLevelName(logging.ERROR, "error")

def expand_path(file):

    # an existing path is taken literally, so names containing [ ] work;
    # anything else is a wildcard pattern
    if os.path.exists(file):
        return [file]
    return sorted(glob.glob(file))

def setup_logging(log_file=None, in_files=()):

    # the log is appended to, so it must never be one of the input files;
    # this has to be checked before the file is opened
    if log_file:
        log_path = Path(log_file).resolve()
        for pattern in in_files:
            if any(Path(match).resolve() == log_path for match in expand_path(pattern)):
                raise TreekitError(
                    "-l file '%s' is also an input file; choose a different file"
                    % log_path
                )

    # reset handlers so that repeated calls do not duplicate output
    log.handlers.clear()
    log.setLevel(logging.INFO)
    log.propagate = False

    console = logging.StreamHandler(sys.stderr)
    console.setLevel(logging.WARNING)
    console.setFormatter(logging.Formatter("treekit: %(levelname)s: %(message)s"))
    log.addHandler(console)

    if log_file:
        try:
            # append, so that repeated runs keep a record of every run
            file_handler = logging.FileHandler(log_file, mode="a", encoding="utf-8")
        except OSError as err:
            raise TreekitError("cannot write log file '%s': %s" % (log_file, err.strerror))
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(
            logging.Formatter("%(asctime)s %(levelname)s: %(message)s"))
        log.addHandler(file_handler)
    return log

def positive_float(value):

    # argparse type for thresholds: a number > 0
    try:
        number = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError("invalid number: '%s'" % value)
    if not number > 0 or math.isinf(number):
        raise argparse.ArgumentTypeError("must be a number > 0, got '%s'" % value)
    return number

def non_negative_float(value):

    try:
        number = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError("invalid number: '%s'" % value)
    if not number >= 0 or math.isinf(number):
        raise argparse.ArgumentTypeError("must be a number >= 0, got '%s'" % value)
    return number

def positive_int(value):

    # argparse type for options such as -c: an integer >= 1
    try:
        number = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("invalid integer value: '%s'" % value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be an integer >= 1, got '%s'" % value)
    return number

class ParsedArgs:

    def __init__(self):

        # initialize the main ArgumentParser object
        parser = argparse.ArgumentParser(usage=usage_text())

        # add the main argument "command"
        parser.add_argument(
            "command",
            help="Subcommand to run"
        )

        # parse the command by excluding the rest of args
        self.args = parser.parse_args(sys.argv[1:2])
        if self.args.command in PLANNED_COMMANDS:
            parser.error("command '%s' is not implemented yet" % self.args.command)
        if self.args.command not in IMPLEMENTED_COMMANDS:
            parser.error("unrecognized command '%s'" % self.args.command)

        # build the parser of the command and parse the rest of args, once
        command = self.args.command
        help_text = dict((name, text) for name, text, _ in COMMANDS)[command]
        sub_parser = argparse.ArgumentParser(
            prog = "treekit " + command,
            description = help_text[0].upper() + help_text[1:],
        )
        getattr(self, command)(sub_parser)
        self.add_common_args(sub_parser)
        self.parsed = sub_parser.parse_args(sys.argv[2:])

    # helpers for arguments shared by several commands

    def add_taxa_args(self, parser, purpose, group=None):

        # the taxa a command works on: from the command line or from a file
        if group is None:
            group = parser.add_mutually_exclusive_group(required=True)
        group.add_argument(
            "-n",
            "--taxa-list",
            nargs = "+",
            dest = "taxa_list",
            type = str,
            help = "Taxa names %s (white-space separated)." % purpose,
        )
        group.add_argument(
            "-f",
            "--taxa-file",
            dest = "taxa_file",
            type = str,
            help = "File with the taxa names %s (one name per line)." % purpose,
        )
        return group

    def add_tree_out_args(self, parser):

        # commands that write trees can write newick or nexus
        parser.add_argument(
            "-u",
            "--out-format",
            dest = "out_format",
            choices = ["newick", "nexus"],
            default = "newick",
            help = "Format for output trees. Default: newick"
        )

    def add_strict_arg(self, parser, help_text):

        parser.add_argument(
            "--strict",
            dest = "strict",
            action = "store_true",
            default = False,
            help = help_text
        )

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

        # parallelization is over input files: one process parses one file
        parser.add_argument(
            "-c",
            "--cores",
            dest = "cores",
            default = 1,
            type = positive_int,
            help = "Number of cores used, at most one per input file. Default: 1"
        )

        # optional log file for all commands
        parser.add_argument(
            "-l",
            "--log",
            dest = "log_file",
            type = str,
            help = "Append run details and warnings to this log file."
        )

    def summary(self, parser):

        # if give option '-e', then output the details
        parser.add_argument(
            "-e",
            "--details-out",
            dest = "details_out",
            action = "store_true",
            default = False,
            help = "Details of each tree summary."
        )

    def taxa(self, parser):

        # list common taxa names
        parser.add_argument(
            "-k",
            "--common-taxa",
            dest = "common_taxa",
            action = "store_true",
            default = False,
            help = "only list common taxa names in all trees."
        )

    def remove(self, parser):

        self.add_taxa_args(parser, "to be removed")
        self.add_tree_out_args(parser)
        self.add_strict_arg(
            parser, "Exit with an error if a taxon to remove is absent from all trees.")

    def extract(self, parser):

        self.add_taxa_args(parser, "to keep in the subtrees")
        self.add_tree_out_args(parser)
        self.add_strict_arg(
            parser, "Exit with an error if a taxon to extract is absent from all trees.")

    def keep(self, parser):

        self.add_taxa_args(parser, "that a tree must contain")
        self.add_tree_out_args(parser)

    def root(self, parser):

        # exactly one way of rooting
        group = parser.add_mutually_exclusive_group(required=True)
        self.add_taxa_args(parser, "of the outgroup", group=group)
        group.add_argument(
            "--midpoint",
            dest = "midpoint",
            action = "store_true",
            help = "Root every tree at the midpoint of its longest tip-to-tip path."
        )
        group.add_argument(
            "--unroot",
            dest = "unroot",
            action = "store_true",
            help = "Remove the root: write unrooted trees."
        )
        self.add_tree_out_args(parser)
        self.add_strict_arg(
            parser, "Exit with an error if a tree cannot be rooted, or if an outgroup "
                    "taxon is absent from all trees.")

    def mrca(self, parser):

        # one group (-n/-f) or several groups from a table (-g)
        group = parser.add_mutually_exclusive_group(required=True)
        self.add_taxa_args(parser, "of one group", group=group)
        group.add_argument(
            "-g",
            "--groups-file",
            dest = "groups_file",
            type = str,
            help = "Tab-separated file of 'taxon<TAB>group' lines, to check several groups."
        )
        parser.add_argument(
            "--extract",
            dest = "extract_clades",
            action = "store_true",
            help = "Write the clade of every group instead of the monophyly table."
        )
        parser.add_argument(
            "--unrooted",
            dest = "force_unrooted",
            action = "store_true",
            help = "Treat every tree as unrooted, even when its root is bifurcating."
        )
        self.add_tree_out_args(parser)
        self.add_strict_arg(
            parser, "Exit with an error if a group taxon is absent from all trees, or if a "
                    "clade cannot be extracted.")

    def draw(self, parser):

        parser.add_argument(
            "-t",
            "--tree",
            nargs = "+",
            dest = "tree_numbers",
            type = positive_int,
            help = "Draw only these trees (numbers counted over all input files)."
        )
        parser.add_argument(
            "-w",
            "--width",
            dest = "width",
            type = positive_int,
            default = 80,
            help = "Width of the drawing in characters. Default: 80"
        )
        parser.add_argument(
            "-m",
            "--metric",
            dest = "metric",
            choices = ["auto", "length", "topology"],
            default = "auto",
            help = "Scale branches by length, or draw the topology only. Default: auto "
                   "(length when every branch has a positive length)"
        )
        parser.add_argument(
            "-s",
            "--supports",
            dest = "show_supports",
            action = "store_true",
            help = "Show internal node labels, such as branch supports."
        )

    def compare(self, parser):

        parser.add_argument(
            "-r",
            "--ref-files",
            nargs = "+",
            dest = "ref_files",
            type = str,
            required = True,
            help = "Tree files of the second source; every input tree is compared with "
                   "every tree of these files."
        )
        parser.add_argument(
            "--rooted",
            dest = "compare_rooted",
            action = "store_true",
            help = "Compare the trees as rooted trees (clades) instead of unrooted (splits)."
        )

    def prune(self, parser):

        parser.add_argument(
            "--max-length",
            dest = "max_length",
            type = positive_float,
            help = "Remove the tips whose branch is longer than this length."
        )
        parser.add_argument(
            "--max-ratio",
            dest = "max_ratio",
            type = positive_float,
            help = "Remove the tips whose branch is longer than this many times the median "
                   "branch length of their tree."
        )
        self.add_tree_out_args(parser)
        self.add_strict_arg(
            parser, "Exit with an error if a tree cannot be pruned (no branch lengths).")

    def collapse(self, parser):

        parser.add_argument(
            "-m",
            "--min-support",
            dest = "min_support",
            type = non_negative_float,
            required = True,
            help = "Collapse the internal branches whose support is lower than this value."
        )
        source = parser.add_mutually_exclusive_group()
        source.add_argument(
            "--support-key",
            dest = "support_key",
            type = str,
            help = "Read the support from this annotation, e.g. posterior for "
                   "[&posterior=0.98] (BEAST, MrBayes). Default: the node label."
        )
        source.add_argument(
            "--support-field",
            dest = "support_field",
            type = positive_int,
            help = "For node labels like 80.5/97, use this '/'-separated field (1 = first)."
        )
        self.add_tree_out_args(parser)
        self.add_strict_arg(
            parser, "Exit with an error if a tree has no support value that can be read.")

    def get_args_dict(self):

        # store arguments in a dictionary
        argument_dictionary = dict(self.args.__dict__)
        argument_dictionary.update(self.parsed.__dict__)
        return argument_dictionary


class FuncTree:
    """Class for encapsulate treekit functions"""

    # code loading at the initialization phase

    def __init__(self, **kwargs):

        # set defaults and get values from kwargs; every option is also kept in
        # self.options for the commands that need one of their own
        self.options = dict(kwargs)
        self.command = kwargs.get("command")
        self.log_file = kwargs.get("log_file")
        self.strict = bool(kwargs.get("strict"))
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

        cores = kwargs.get("cores")
        self.cores = 1 if cores is None else self.get_cores(cores)

        # the second source of compare
        self.ref_files = self.get_paths(kwargs.get("ref_files") or [])
        self.ref_formats = [self.detect_format(path) for path in self.ref_files]

        # never write the output over one of the input files
        if self.out_file is not None and self.out_file in self.in_files + self.ref_files:
            raise TreekitError(
                "-o file '%s' is also an input file; choose a different file"
                % self.out_file
            )

        # specific attributes for each command
        self.details_out = bool(kwargs.get("details_out"))
        self.common_taxa = bool(kwargs.get("common_taxa"))
        self.taxa_list = kwargs.get("taxa_list")
        if kwargs.get("taxa_file"):
            self.taxa_file = self.get_path(kwargs.get("taxa_file"))
        else:
            self.taxa_file = None

    def get_paths(self, files):

        # each in_file may be a relative or absolute path, with optional wildcards
        out_paths = []
        for file in files:

            matches = expand_path(file)
            if not matches:
                raise TreekitError("input file not found: '%s'" % file)
            for match in matches:
                if not os.path.isfile(match):
                    raise TreekitError("input is not a file: '%s'" % match)
                out_paths.append(Path(match).resolve())
        return out_paths

    def get_path(self, file):

        # get outfile absolute path
        out_path = Path(file).resolve()
        return out_path

    def get_formats(self):

        # identify tree in nexus or newick
        in_formats=[]
        for in_file in self.in_files:
            in_formats.append(self.detect_format(in_file))
        return in_formats

    def detect_format(self, in_file):

        # only the beginning of the file is needed to tell the format
        try:
            with open(in_file, "rb") as f:
                head = f.read(65536)
        except OSError as err:
            raise TreekitError("cannot read '%s': %s" % (in_file, err.strerror))

        # the tree parser cannot read a byte-order mark
        if head.startswith(b"\xef\xbb\xbf"):
            raise TreekitError(
                "'%s' starts with a UTF-8 byte-order mark (BOM), which cannot be "
                "parsed; please save the file without BOM" % in_file
            )

        # skip blank lines and leading [comments], e.g. "[&R]"
        text = head.decode("utf-8", errors="replace")
        while True:
            text = text.lstrip()
            if not text.startswith("["):
                break
            end = text.find("]")
            text = text[end + 1:] if end != -1 else ""

        if text.startswith("#"):
            return "nexus"
        if text.startswith("(") or BARE_NEWICK.match(text):
            return "newick"
        raise TreekitError(
            "'%s' is not a newick or nexus tree file (empty or unrecognized content)"
            % in_file
        )

    def read_trees(self, in_file, in_format, tree_list=None, **kwargs):

        # parse one file; parser errors become treekit errors, not tracebacks.
        # labels are case sensitive: "a" and "A" are two taxa, and removing "a"
        # must not remove "A" (the dendropy default would merge them)
        try:
            if tree_list is None:
                tree_list = new_tree_list()
            tree_list.read(
                path = in_file,
                schema = in_format,
                case_sensitive_taxon_labels = True,
                **kwargs
            )
            return tree_list
        except DataParseError as err:
            raise TreekitError(
                "cannot parse '%s' as %s: %s" % (in_file, in_format, " ".join(str(err).split()))
            )
        except UnicodeDecodeError:
            raise TreekitError("'%s' is not valid UTF-8 text" % in_file)

    def get_cores(self, cores):

        # check the cores attribute
        if not isinstance(cores, int) or cores < 1:
            raise TreekitError("-c/--cores requires an integer >= 1")

        # if cores > real_cores, use real_cores
        real_cores = mp.cpu_count()
        if cores > real_cores:
            log.warning(
                "this computer has %d cores, fewer than the %d requested; "
                "using %d cores" % (real_cores, cores, real_cores)
            )
            return real_cores
        return cores

    def map_files(self, function, files=None, formats=None):

        # apply function(in_file, in_format) to every input file, one process per file
        if files is None:
            files, formats = self.in_files, self.in_formats
        jobs = list(zip(files, formats))

        # a single file is handled in-process (no pool, no overhead)
        if self.cores > 1 and len(jobs) > 1:
            with mp.Pool(min(self.cores, len(jobs))) as pool:
                return pool.starmap(function, jobs)
        return [function(in_file, in_format) for in_file, in_format in jobs]

    # code for Meta Functions on a tree object

    def is_rooted(self, tree_object):

        # the tree is rooted when seed node has two children;
        # three or more children (unrooted) or a single child are not rooted
        return len(tree_object.seed_node.child_nodes()) == 2

    def is_binary(self, tree_object):

        # strictly bifurcating: no node other than the root may have > 2 children
        seed_node = tree_object.seed_node
        for node in tree_object.preorder_node_iter():
            if node is not seed_node and len(node.child_nodes()) > 2:
                return False

        # the root has two children if rooted, or three if unrooted;
        # more children (a star-like root) is not binary
        return len(seed_node.child_nodes()) in (2, 3)

    def has_branchlength(self, tree_object):

        # the tree with branch length has a positive total length
        return tree_object.length() > 0

    def is_ultrametric(self, tree_object):

        # two conditions required in an ultrametric tree
        # rooted tree and the same tip heights
        if self.is_rooted(tree_object) and tree_object.length() > 0:
            min_height, max_height = tree_object.minmax_leaf_distance_from_root()

            # heights are floats: compare with a tolerance, not exact equality
            return math.isclose(min_height, max_height,
                                rel_tol=ULTRAMETRIC_TOL, abs_tol=ULTRAMETRIC_TOL)
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

    def tip_labels(self, tree_object):

        # labels of all tips; tips without a label are skipped
        return set(
            leaf.taxon.label for leaf in tree_object.leaf_node_iter()
            if leaf.taxon is not None
        )

    # methods corresponding to summary block

    def get_summaries(self):

        # get summaries for each of tree objects
        header = [
            "source",
            "index",
            "is_binary",
            "is_rooted",
            "has_brlen",
            "is_ultrametric",
            "tree_length",
            "avg_leaf_height",
            "leaf_height_sd",
            "num_tips"
        ]

        # one list of rows per file, in the order of the input files
        per_file = self.map_files(self.get_file_summaries)
        summaries = [row for rows in per_file for row in rows]
        return header, summaries

    def get_file_summaries(self, in_file, in_format):

        # one row (all values as strings) for each tree of one file
        tree_list = self.read_trees(in_file, in_format)
        if len(tree_list) == 0:
            raise TreekitError("no trees found in '%s'" % in_file)

        rows = []
        for tree_idx, tree_object in enumerate(tree_list, start=1):

            # if leaf heights list is not None, stat the leaf heights
            leaf_heights = self.get_leaf_heights(tree_object)
            if leaf_heights:
                avg_leaf_height = decimal_correct(statistics.mean(leaf_heights))
                leaf_height_sd = decimal_correct(statistics.pstdev(leaf_heights))
            else:
                avg_leaf_height = 'NA'
                leaf_height_sd = 'NA'

            values = (
                in_file.name,
                tree_idx,
                self.is_binary(tree_object),
                self.is_rooted(tree_object),
                self.has_branchlength(tree_object),
                self.is_ultrametric(tree_object),
                decimal_correct(tree_object.length()),
                avg_leaf_height,
                leaf_height_sd,
                len(tree_object.leaf_nodes()),
            )
            rows.append(tuple(str(value) for value in values))

        return rows

    def get_summstats(self):

        # get summstat for all tree objects, one row per input file
        header = [
            "file_name",
            "num_tree",
            "num_binary",
            "num_rooted",
            "num_hasbrlen",
            "num_ultrametric",
            "avg_tree_length",
            "tree_length_sd",
            "num_taxa",
            "num_common_taxa"
        ]

        # each worker parses one file and returns one small row
        summstats = self.map_files(self.get_file_summstat)
        return header, summstats

    def get_file_summstat(self, in_file, in_format):

        # parse one file and compute its summary row (all values as strings)
        tree_list = self.read_trees(in_file, in_format)
        if len(tree_list) == 0:
            raise TreekitError("no trees found in '%s'" % in_file)

        # compute the per-tree properties once, then aggregate
        metrics = [self.tree_metrics(tree_object) for tree_object in tree_list]
        lengths = [m[4] for m in metrics]
        values = (
            len(tree_list),
            sum(1 for m in metrics if m[0]),
            sum(1 for m in metrics if m[1]),
            sum(1 for m in metrics if m[2]),
            sum(1 for m in metrics if m[3]),
            decimal_correct(statistics.mean(lengths)),
            decimal_correct(statistics.pstdev(lengths)),
            len(tree_list.taxon_namespace),
            len(self.get_common_taxa(tree_list)),
        )
        return (in_file.name,) + tuple(str(v) for v in values)

    def tree_metrics(self, tree_object):

        # all per-tree properties needed by get_summstats, computed once per tree
        return (
            self.is_binary(tree_object),
            self.is_rooted(tree_object),
            self.has_branchlength(tree_object),
            self.is_ultrametric(tree_object),
            tree_object.length(),
        )

    def write_summaries(self):

        # detailed (one row per tree) or brief (one row per file) summaries
        if self.details_out:
            header, rows = self.get_summaries()
        else:
            header, rows = self.get_summstats()

        lines = ['\t'.join(header)] + ['\t'.join(row) for row in rows]
        self.emit('\n'.join(lines), "summaries")

    # methods corresponding to taxa code block

    def read_all_trees(self, **kwargs):

        # every file is parsed on its own, because a NEXUS file declares the number
        # of its taxa and cannot be read into a namespace already filled by another
        # file; the trees are then merged into one shared namespace
        tree_list = new_tree_list()
        for in_file, in_format in zip(self.in_files, self.in_formats):
            file_trees = self.read_trees(
                in_file,
                in_format,
                preserve_underscores = True,
                **kwargs
            )
            for tree in file_trees:
                tree_list.append(tree)
        return tree_list

    def get_taxa_list(self):

        # all trees of all input files, and the taxa names in their shared namespace
        tree_list = self.read_all_trees()
        taxa_space = tree_list.taxon_namespace
        taxa_list = []

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

        if len(tree_list) == 0:
            raise TreekitError("no trees found in the input files")

        # intersection of the tip labels of all trees
        taxa_trees = [self.tip_labels(tree) for tree in tree_list]
        return list(set.intersection(*taxa_trees))

    def write_taxa_out(self):

        # get taxa for all trees or get the common taxa names
        taxa_out = self.get_taxa_list()
        self.emit('\n'.join(taxa_out), "taxa list")

    # methods corresponding to remove code block

    def get_taxa_to_remove(self):

        return self.get_requested_taxa("to remove")

    def get_requested_taxa(self, purpose):

        # the taxa names given with -n or with -f (one name per line)
        if self.taxa_file:
            if not self.taxa_file.is_file():
                raise TreekitError("taxa file not found: '%s'" % self.taxa_file)
            try:
                with open(self.taxa_file, "r", encoding="utf-8") as taxa_file:
                    names = [line.strip() for line in taxa_file]
            except UnicodeDecodeError:
                raise TreekitError("'%s' is not valid UTF-8 text" % self.taxa_file)
        else:
            names = [name.strip() for name in self.taxa_list]

        # drop blank lines and duplicates, keeping the given order
        names = list(dict.fromkeys(name for name in names if name))
        if not names:
            raise TreekitError("no taxa given %s" % purpose)
        return names

    def get_reduced_trees(self):

        # read all input trees into one shared taxon namespace
        tree_list = self.read_all_trees()

        # check the request against every tree before anything is changed
        taxa_to_remove = self.get_taxa_to_remove()
        self.check_removal(tree_list, taxa_to_remove)

        # prune the taxa based on labels
        for tree in tree_list:
            tree.prune_taxa_with_labels(taxa_to_remove)
            clear_stem(tree)

        # no tree uses the removed taxa any more, so do not declare them in the output
        taxon_namespace = tree_list.taxon_namespace
        for name in taxa_to_remove:
            taxon = taxon_namespace.get_taxon(name)
            if taxon is not None:
                taxon_namespace.remove_taxon(taxon)

        return tree_list

    def check_removal(self, tree_list, taxa_to_remove):

        # count, for each requested taxon, the trees that contain it
        total = len(tree_list)
        requested = set(taxa_to_remove)
        found = Counter()
        single_tip_trees = []
        for index, tree in enumerate(tree_list, start=1):
            present = set()
            remaining = 0
            for leaf in tree.leaf_node_iter():

                # pruning would delete a tip without a label as well, silently
                if leaf.taxon is None:
                    raise TreekitError(
                        "tree %d has a tip without a label; taxa cannot be removed "
                        "from trees with unlabelled tips" % index
                    )
                present.add(leaf.taxon.label)
                if leaf.taxon.label not in requested:
                    remaining += 1
            found.update(requested & present)

            # a tree with no tips left cannot be written, so stop with an error
            if remaining == 0:
                raise TreekitError(
                    "removing these taxa would leave tree %d with no tips" % index
                )
            if remaining == 1:
                single_tip_trees.append(index)

        # one summary line for all trees, instead of one line per tree
        self.warn_single_tip_trees(single_tip_trees, "removing these taxa")
        self.report_taxa_presence(taxa_to_remove, found, total, "to remove", "removed from")

    def warn_single_tip_trees(self, tree_numbers, action):

        # one warning for all the trees left with a single tip
        if tree_numbers:
            log.warning(
                "%s leaves %d tree(s) with a single tip (tree %s)"
                % (action, len(tree_numbers), format_numbers(tree_numbers))
            )

    def report_taxa_presence(self, names, found, total, purpose, done):

        # taxa absent from every tree: nothing would be done with them
        absent = [name for name in names if found[name] == 0]
        if absent and len(absent) == len(names):
            noun = "taxon" if len(names) == 1 else "taxa"
            raise TreekitError(
                "none of the %d %s %s were found in any input tree: %s"
                % (len(names), noun, purpose, ", ".join(absent))
            )

        # one message for all absent taxa, instead of one line per taxon
        if absent:
            noun, verb = ("taxon", "was") if len(absent) == 1 else ("taxa", "were")
            message = "%d %s %s not found in any input tree: %s" % (
                len(absent), noun, verb, ", ".join(absent))
            if self.strict:
                raise TreekitError(message + " (--strict)")
            log.warning(message)

        for name in names:
            count = found[name]
            if 0 < count < total:
                # partial presence is normal for sets of gene trees: log, do not warn
                log.info("taxon '%s' %s %d of %d trees" % (name, done, count, total))
            elif count == total:
                log.info("taxon '%s' %s all %d trees" % (name, done, total))

    def write_reduced_trees(self):

        self.write_trees(self.get_reduced_trees())

    # methods corresponding to extract code block

    def get_extracted_trees(self):

        # keep only the requested taxa in every tree
        tree_list = self.read_all_trees()
        taxa_to_keep = self.get_requested_taxa("to extract")
        wanted = set(taxa_to_keep)

        # check every tree before anything is changed
        found = Counter()
        empty_trees = []
        single_tip_trees = []
        for index, tree in enumerate(tree_list, start=1):
            present = wanted & self.tip_labels(tree)
            found.update(present)
            if not present:
                empty_trees.append(index)
            elif len(present) == 1:
                single_tip_trees.append(index)

        # taxa missing everywhere are reported first: that is the likely mistake
        self.report_taxa_presence(
            taxa_to_keep, found, len(tree_list), "to extract", "extracted from")
        if empty_trees:
            raise TreekitError(
                "%d tree(s) contain none of the taxa to extract (tree %s); use "
                "'treekit keep' first to select the trees that contain them"
                % (len(empty_trees), format_numbers(empty_trees))
            )
        self.warn_single_tip_trees(single_tip_trees, "extracting these taxa")

        # tips not asked for are removed, including tips without a label
        for tree in tree_list:
            tree.retain_taxa_with_labels(wanted)
            clear_stem(tree)
        self.drop_unused_taxa(tree_list)
        return tree_list

    def write_extracted_trees(self):

        self.write_trees(self.get_extracted_trees())

    # methods corresponding to keep code block

    def get_kept_trees(self):

        # select, unchanged, the trees that contain every requested taxon
        tree_list = self.read_all_trees()
        required = self.get_requested_taxa("that a tree must contain")
        wanted = set(required)

        found = Counter()
        kept = new_tree_list(tree_list.taxon_namespace)
        kept_numbers = []
        for index, tree in enumerate(tree_list, start=1):
            present = wanted & self.tip_labels(tree)
            found.update(present)
            if present == wanted:
                kept.append(tree)
                kept_numbers.append(index)

        total = len(tree_list)
        absent = [name for name in required if found[name] == 0]
        if absent:
            noun, verb = ("taxon", "is") if len(absent) == 1 else ("taxa", "are")
            raise TreekitError(
                "no tree can be kept: %d %s %s not in any input tree: %s"
                % (len(absent), noun, verb, ", ".join(absent))
            )
        if not kept_numbers:
            raise TreekitError(
                "none of the %d input trees contains all %d requested taxa"
                % (total, len(required))
            )

        for name in required:
            log.info("taxon '%s' is in %d of %d trees" % (name, found[name], total))
        log.info("kept %d of %d trees (tree %s)"
                 % (len(kept_numbers), total, format_numbers(kept_numbers)))
        self.drop_unused_taxa(kept)
        return kept

    def write_kept_trees(self):

        self.write_trees(self.get_kept_trees())

    # methods corresponding to root code block

    def get_rooted_trees(self):

        tree_list = self.read_all_trees()
        total = len(tree_list)

        if self.options.get("unroot"):
            for tree in tree_list:
                keep_split_labels(tree, unroot_tree)
            log.info("unrooted %d trees" % total)
            return tree_list

        # outgroup rooting checks the outgroup against all trees first
        if self.options.get("midpoint"):
            mode = "midpoint"
            root_one = root_at_midpoint
        else:
            mode = "outgroup"
            outgroup = self.get_requested_taxa("of the outgroup")
            found = Counter()
            for tree in tree_list:
                found.update(set(outgroup) & self.tip_labels(tree))
            self.report_taxa_presence(outgroup, found, total, "of the outgroup", "used in")
            root_one = lambda tree: root_on_outgroup(tree, outgroup)

        # a tree that cannot be rooted is written unchanged, so the output keeps
        # the number and the order of the input trees
        failed = {}
        for index, tree in enumerate(tree_list, start=1):
            reason = keep_split_labels(tree, root_one)
            if reason:
                failed.setdefault(reason, []).append(index)

        for reason, numbers in failed.items():
            message = "%d tree(s) cannot be rooted (%s): tree %s" % (
                len(numbers), reason, format_numbers(numbers))
            if self.strict:
                raise TreekitError(message + " (--strict)")
            log.warning(message + "; written unchanged")
        not_rooted = sum(len(numbers) for numbers in failed.values())
        log.info("rooted %d of %d trees (%s)" % (total - not_rooted, total, mode))
        return tree_list

    def write_rooted_trees(self):

        self.write_trees(self.get_rooted_trees())

    # methods corresponding to mrca code block

    def get_groups(self):

        # {group name: [taxa]} in the order of the input
        groups_file = self.options.get("groups_file")
        if not groups_file:
            return {"group": self.get_requested_taxa("of the group")}

        path = self.get_path(groups_file)
        if not path.is_file():
            raise TreekitError("groups file not found: '%s'" % path)
        groups = {}
        try:
            with open(path, "r", encoding="utf-8") as handle:
                for number, line in enumerate(handle, start=1):
                    line = line.rstrip("\r\n")
                    if not line.strip() or line.lstrip().startswith("#"):
                        continue
                    fields = line.split("\t")
                    if len(fields) != 2 or not fields[0].strip() or not fields[1].strip():
                        raise TreekitError(
                            "line %d of '%s' is not 'taxon<TAB>group': %r" % (number, path, line))
                    taxon, name = fields[0].strip(), fields[1].strip()
                    members = groups.setdefault(name, [])
                    if taxon not in members:
                        members.append(taxon)
        except UnicodeDecodeError:
            raise TreekitError("'%s' is not valid UTF-8 text" % path)
        if not groups:
            raise TreekitError("no groups in '%s'" % path)
        return groups

    def read_all_trees_with_sources(self):

        # like read_all_trees, with the file name and index of every tree
        tree_list = new_tree_list()
        sources = []
        for in_file, in_format in zip(self.in_files, self.in_formats):
            file_trees = self.read_trees(in_file, in_format, preserve_underscores = True)
            for index, tree in enumerate(file_trees, start=1):
                tree_list.append(tree)
                sources.append((in_file.name, index))
        return tree_list, sources

    def get_mrca_results(self):

        # test every group in every tree
        tree_list, sources = self.read_all_trees_with_sources()
        groups = self.get_groups()
        force_unrooted = bool(self.options.get("force_unrooted"))

        # tips without a label cannot be told apart from group members
        for number, tree in enumerate(tree_list, start=1):
            if any(leaf.taxon is None for leaf in tree.leaf_node_iter()):
                raise TreekitError(
                    "tree %d has a tip without a label; monophyly cannot be tested" % number)

        all_taxa = list(dict.fromkeys(t for members in groups.values() for t in members))
        found = Counter()
        for tree in tree_list:
            found.update(set(all_taxa) & self.tip_labels(tree))
        self.report_taxa_presence(all_taxa, found, len(tree_list), "of the groups", "found in")

        results = []
        for number, (tree, source) in enumerate(zip(tree_list, sources), start=1):
            for name, members in groups.items():
                results.append((number, tree, source, name, len(members),
                                group_monophyly(tree, members, force_unrooted)))
        return tree_list, results

    def write_mrca(self):

        tree_list, results = self.get_mrca_results()
        if self.options.get("extract_clades"):
            self.write_mrca_clades(tree_list, results)
            return

        header = ["source", "index", "group", "rooted", "members_present", "members_total",
                  "monophyletic", "num_intruders", "intruders"]
        rows = []
        tested = Counter()
        monophyletic = Counter()
        for number, tree, source, name, total, status in results:
            mono = status["monophyletic"]
            intruders = status["intruders"]
            rows.append((
                source[0], source[1], name, status["rooted"], status["present"], total,
                "NA" if mono is None else mono,
                "NA" if intruders is None else len(intruders),
                "NA" if intruders is None else ",".join(intruders),
            ))
            if mono is not None:
                tested[name] += 1
                monophyletic[name] += int(mono)
        for name in dict.fromkeys(r[3] for r in results):
            log.info("group '%s': monophyletic in %d of %d testable trees"
                     % (name, monophyletic[name], tested[name]))

        lines = ["\t".join(header)] + ["\t".join(str(v) for v in row) for row in rows]
        self.emit("\n".join(lines), "monophyly table")

    def write_mrca_clades(self, tree_list, results):

        # one clade per tree and group; a group that is not a clade is skipped
        clades = new_tree_list(tree_list.taxon_namespace)
        skipped = {}
        for number, tree, source, name, total, status in results:
            clade = extract_group_clade(tree, status)
            if isinstance(clade, str):
                skipped.setdefault(clade, []).append("%d/%s" % (number, name))
                continue
            clade.label = "%s_tree%d" % (name, number)
            clades.append(clade)
            log.info("extracted the clade of group '%s' from tree %d (%s, tree %d)"
                     % (name, number, source[0], source[1]))

        for reason, items in skipped.items():
            message = "%d clade(s) cannot be extracted (%s): tree/group %s" % (
                len(items), reason, format_numbers(items))
            if self.strict:
                raise TreekitError(message + " (--strict)")
            log.warning(message + "; skipped")
        if len(clades) == 0:
            raise TreekitError("no clade could be extracted")
        self.drop_unused_taxa(clades)
        self.write_trees(clades)

    # methods corresponding to draw code block

    def write_drawings(self):

        tree_list, sources = self.read_all_trees_with_sources()
        total = len(tree_list)
        numbers = self.options.get("tree_numbers") or list(range(1, total + 1))
        too_large = [n for n in numbers if n > total]
        if too_large:
            raise TreekitError(
                "there are %d input trees; no tree %s" % (total, format_numbers(too_large)))
        width = self.options.get("width") or 80
        if width < 20:
            raise TreekitError("-w/--width must be at least 20")

        drawings = []
        topology_only = []
        for number in dict.fromkeys(numbers):
            tree = tree_list[number - 1]
            source, index = sources[number - 1]
            picture, metric = draw_tree(
                tree, width, self.options.get("metric") or "auto",
                bool(self.options.get("show_supports")))
            if metric == "topology" and self.options.get("metric") == "auto":
                topology_only.append(number)
            title = "tree %d (%s, tree %d)" % (number, source, index)
            if tree.label:
                title += " %s" % tree.label
            drawings.append(title + "\n" + picture)
        if topology_only:
            log.info("drawn as topology, because branch lengths are missing, zero, or too "
                     "uneven to draw to scale at this width: tree %s"
                     % format_numbers(topology_only))
        self.emit("\n\n".join(drawings), "drawings")

    # methods corresponding to compare code block

    def get_file_shapes(self, in_file, in_format):

        # a compact description of every tree of a file, cheap to send between
        # processes: tip labels and, per node in preorder, the parent and the length
        shapes = []
        for number, tree in enumerate(self.read_trees(in_file, in_format,
                                                      preserve_underscores=True), start=1):
            nodes = list(tree.preorder_node_iter())
            position = {node: i for i, node in enumerate(nodes)}
            parents = []
            lengths = []
            labels = []
            for node in nodes:
                parents.append(-1 if node is tree.seed_node else position[node.parent_node])
                lengths.append(None if node is tree.seed_node or node.edge.length is None
                               else float(node.edge.length))
                if node.is_leaf():
                    if node.taxon is None:
                        raise TreekitError(
                            "tree %d of '%s' has a tip without a label; trees with unlabelled "
                            "tips cannot be compared" % (number, in_file))
                    labels.append(node.taxon.label)
                else:
                    labels.append(None)
            rooted = sum(1 for p in parents if p == 0) == 2
            shapes.append((in_file.name, number, parents, lengths, labels, rooted))
        if not shapes:
            raise TreekitError("no trees found in '%s'" % in_file)
        return shapes

    def get_comparison(self):

        rooted_mode = bool(self.options.get("compare_rooted"))
        shapes_a = [s for f in self.map_files(self.get_file_shapes) for s in f]
        shapes_b = [s for f in self.map_files(self.get_file_shapes, self.ref_files,
                                              self.ref_formats) for s in f]

        # one bit per taxon over both sources
        bit = {}
        for shape in shapes_a + shapes_b:
            for label in shape[4]:
                if label is not None and label not in bit:
                    bit[label] = 1 << len(bit)
        label_of = {b: label for label, b in bit.items()}
        trees_a = [TreeSplits(shape, bit) for shape in shapes_a]
        trees_b = [TreeSplits(shape, bit) for shape in shapes_b]

        if rooted_mode:
            unrooted = [t.name for t in trees_a + trees_b if not t.rooted]
            if unrooted:
                log.warning("--rooted: %d tree(s) have no bifurcating root and are compared "
                            "as written: %s" % (len(unrooted), format_numbers(unrooted)))

        header = ["source_a", "index_a", "source_b", "index_b", "common_taxa",
                  "only_a", "only_b", "rf", "norm_rf", "weighted_rf"]
        rows = []
        too_small = 0
        pruned = 0
        min_taxa = 3 if rooted_mode else 4
        for a in trees_a:
            for b in trees_b:
                common = a.taxa & b.taxa
                n_common = popcount(common)
                only_a = popcount(a.taxa & ~common)
                only_b = popcount(b.taxa & ~common)
                pruned += bool(only_a or only_b)
                if n_common < min_taxa:
                    too_small += 1
                    rows.append((a.source, a.index, b.source, b.index, n_common,
                                 only_a, only_b, 0, "NA", "NA"))
                    continue
                splits_a, weights_a = a.restricted(common, rooted_mode)
                splits_b, weights_b = b.restricted(common, rooted_mode)
                rf = len(splits_a ^ splits_b)
                denominator = len(splits_a) + len(splits_b)
                norm = decimal_correct(rf / denominator) if denominator else "NA"
                if a.has_lengths and b.has_lengths:
                    weighted = decimal_correct(sum(
                        abs(weights_a.get(key, 0.0) - weights_b.get(key, 0.0))
                        for key in set(weights_a) | set(weights_b)))
                else:
                    weighted = "NA"
                rows.append((a.source, a.index, b.source, b.index, n_common,
                             only_a, only_b, rf, norm, weighted))

        if too_small:
            log.warning("%d pair(s) share fewer than %d taxa: rf is 0 and norm_rf is NA"
                        % (too_small, min_taxa))
        log.info("compared %d x %d trees (%s); %d pair(s) compared on their common taxa only"
                 % (len(trees_a), len(trees_b), "rooted" if rooted_mode else "unrooted", pruned))
        return header, rows

    def write_comparison(self):

        header, rows = self.get_comparison()
        lines = ["\t".join(header)] + ["\t".join(str(v) for v in row) for row in rows]
        self.emit("\n".join(lines), "comparison")

    # methods corresponding to prune code block

    def get_pruned_trees(self):

        max_length = self.options.get("max_length")
        max_ratio = self.options.get("max_ratio")
        if max_length is None and max_ratio is None:
            raise TreekitError("give --max-length, --max-ratio, or both")

        tree_list = self.read_all_trees()
        failed = {}
        plans = []
        single_tip_trees = []
        removed_total = 0
        for index, tree in enumerate(tree_list, start=1):
            edges = [node.edge for node in tree.preorder_node_iter() if node is not tree.seed_node]
            lengths = [float(edge.length) for edge in edges if edge.length is not None]
            if not lengths or not any(lengths):
                failed.setdefault("no branch lengths", []).append(index)
                plans.append([])
                continue

            # the thresholds of this tree, from the branch lengths before pruning
            limits = []
            if max_length is not None:
                limits.append(max_length)
            median = statistics.median(lengths)
            if max_ratio is not None:
                if median > 0:
                    limits.append(max_ratio * median)
                elif max_length is None:
                    failed.setdefault("the median branch length is 0", []).append(index)
                    plans.append([])
                    continue
            limit = min(limits)

            long_tips = [leaf for leaf in tree.leaf_node_iter()
                         if leaf.edge.length is not None and float(leaf.edge.length) > limit]
            n_tips = len(tree.leaf_nodes())
            if long_tips and len(long_tips) == n_tips:
                raise TreekitError(
                    "every tip of tree %d has a branch longer than %s; nothing would be left"
                    % (index, format_number(limit)))
            if long_tips and n_tips - len(long_tips) == 1:
                single_tip_trees.append(index)
            for leaf in long_tips:
                log.info("tree %d: removed %s (branch %s, %s x the median %s)" % (
                    index, leaf.taxon.label if leaf.taxon is not None else "(no label)",
                    format_number(leaf.edge.length),
                    format_number(float(leaf.edge.length) / median) if median > 0 else "NA",
                    format_number(median)))
            removed_total += len(long_tips)
            plans.append(long_tips)

        for reason, numbers in failed.items():
            message = "%d tree(s) cannot be pruned (%s): tree %s" % (
                len(numbers), reason, format_numbers(numbers))
            if self.strict:
                raise TreekitError(message + " (--strict)")
            log.warning(message + "; written unchanged")
        self.warn_single_tip_trees(single_tip_trees, "pruning")

        for tree, long_tips in zip(tree_list, plans):
            for leaf in long_tips:
                tree.prune_subtree(leaf, suppress_unifurcations=True)
            if long_tips:
                clear_stem(tree)

        changed = sum(1 for plan in plans if plan)
        log.info("removed %d tip(s) from %d of %d trees" % (removed_total, changed, len(tree_list)))
        if removed_total == 0:
            log.warning("no tip has a branch longer than the threshold; trees are unchanged")
        self.drop_unused_taxa(tree_list)
        return tree_list

    def write_pruned_trees(self):

        self.write_trees(self.get_pruned_trees())

    # methods corresponding to collapse code block

    def get_collapsed_trees(self):

        threshold = self.options.get("min_support")
        key = self.options.get("support_key")
        field = self.options.get("support_field")

        # annotations such as [&posterior=0.98] are only read when they are used
        tree_list = self.read_all_trees(extract_comment_metadata=bool(key))
        self.keep_annotations = bool(key)

        failed = {}
        values = []
        unreadable = 0
        collapsed_total = 0
        changed = 0
        for index, tree in enumerate(tree_list, start=1):
            supports = {}
            for node in tree.postorder_internal_node_iter():
                if node is tree.seed_node:
                    continue
                raw = node_support_text(node, key)
                if raw is None:
                    continue
                value = parse_support(raw, field)
                if value is None:
                    unreadable += 1
                else:
                    supports[node] = value
            if not supports:
                reason = ("no '%s' annotation" % key) if key else \
                         "no numeric support in the node labels"
                failed.setdefault(reason, []).append(index)
                continue
            values.extend(supports.values())

            # collapse from the tips up; a collapsed branch's length goes to its children
            low = [node for node in supports if supports[node] < threshold]
            for node in low:
                node.edge.collapse(adjust_collapsed_head_children_edge_lengths=True)
            if low:
                changed += 1
                collapsed_total += len(low)

        for reason, numbers in failed.items():
            message = "%d tree(s) have %s: tree %s" % (len(numbers), reason, format_numbers(numbers))
            if not key and reason.startswith("no numeric"):
                message += " (for labels like 80/95 use --support-field)"
            if self.strict:
                raise TreekitError(message + " (--strict)")
            log.warning(message + "; written unchanged")
        if unreadable:
            log.warning("%d node label(s) are not numbers and were left as they are%s"
                        % (unreadable, "" if field else
                           " (for labels like 80/95 use --support-field)"))

        # a threshold on the wrong scale would collapse everything or nothing
        if values:
            largest = max(values)
            if largest <= 1 and threshold > 1:
                log.warning("supports are between 0 and 1 but --min-support is %s: "
                            "every supported branch is collapsed" % format_number(threshold))
            elif largest > 1 and 0 < threshold < 1:
                log.warning("supports go up to %s but --min-support is %s: did you mean %s?"
                            % (format_number(largest), format_number(threshold),
                               format_number(threshold * 100)))
        log.info("collapsed %d branch(es) in %d of %d trees"
                 % (collapsed_total, changed, len(tree_list)))
        return tree_list

    def write_collapsed_trees(self):

        self.write_trees(self.get_collapsed_trees())

    # general methods for commands that write trees

    def drop_unused_taxa(self, tree_list):

        # do not declare, in nexus output, taxa that no output tree uses
        used = set()
        for tree in tree_list:
            used |= self.tip_labels(tree)
        taxon_namespace = tree_list.taxon_namespace
        for taxon in list(taxon_namespace):
            if taxon.label not in used:
                taxon_namespace.remove_taxon(taxon)

    def write_trees(self, tree_list):

        # keep quotes where labels need them, e.g. 'A(1)' or 'B,2'
        tree_out = tree_list.as_string(
            schema=self.out_format,
            unquoted_underscores=True,
            suppress_annotations=not getattr(self, "keep_annotations", False)
        )
        self.emit(tree_out, "trees")

    # general methods for all common blocks

    def emit(self, text, what):

        # write the result to the outfile if given, or to stdout; always end with a newline
        text = text.rstrip("\n") + "\n"
        if self.out_file is None:
            sys.stdout.write(text)
            return

        self.file_overwrite_warning()
        try:
            with open(self.out_file, "w", encoding="utf-8") as out_file:
                out_file.write(text)
        except OSError as err:
            raise TreekitError("cannot write '%s': %s" % (self.out_file, err.strerror))
        log.info("wrote %s to '%s'" % (what, self.out_file))

    def file_overwrite_warning(self):

        # warn when overwriting a file
        if self.out_file.exists():
            log.warning("overwriting '%s'" % self.out_file)

# general functions for global invoking

def new_tree_list(taxon_namespace=None):

    # an empty tree list; a new taxon namespace is case sensitive
    if taxon_namespace is None:
        taxon_namespace = dendropy.TaxonNamespace(is_case_sensitive = True)
    return dendropy.TreeList(taxon_namespace = taxon_namespace)

# ---------------------------------------------------------------- tree algorithms

def leaf_bit_map(tree):

    # each tip gets one bit; tip nodes are kept by the rerooting operations
    return {leaf: 1 << i for i, leaf in enumerate(tree.leaf_node_iter())}

def node_masks(tree, bits):

    # for every node, the bit mask of the tips below it
    masks = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            masks[node] = bits[node]
        else:
            mask = 0
            for child in node.child_node_iter():
                mask |= masks[child]
            masks[node] = mask
    return masks

def split_key(mask, full):

    # a split and its complement are the same split: use the side without tip 0
    return full ^ mask if mask & 1 else mask

def keep_split_labels(tree, operation):

    # node labels (e.g. branch supports) belong to the split of the branch above
    # the node; rerooting moves nodes, so labels are put back by split afterwards
    bits = leaf_bit_map(tree)
    full = (1 << len(bits)) - 1
    labels = {}
    for node, mask in node_masks(tree, bits).items():
        if node.is_internal() and node is not tree.seed_node and node.label is not None:
            labels[split_key(mask, full)] = node.label

    result = operation(tree)

    for node, mask in node_masks(tree, bits).items():
        if node.is_internal():
            node.label = None if node is tree.seed_node else labels.get(split_key(mask, full))
    return result

def unroot_tree(tree):

    # merge the two branches at a bifurcating root into one branch
    if len(tree.seed_node.child_nodes()) == 2 and len(tree.leaf_nodes()) > 2:
        tree.deroot()
    tree.is_rooted = False

def float_or_none(value):

    return None if value is None else float(value)

def root_on_edge(tree, edge, distance_from_tail):

    # put a bifurcating root on edge, at distance_from_tail from its tail node
    length = edge.length
    if length is None or distance_from_tail is None:
        tail_part = head_part = None
    else:
        tail_part = float(distance_from_tail)
        head_part = float(length) - tail_part
    tree.reroot_at_edge(edge, length1=tail_part, length2=head_part)
    tree.is_rooted = True

def root_on_outgroup(tree, outgroup):

    # returns None when rooted, or the reason why the tree cannot be rooted
    wanted = set(outgroup)
    is_out = {leaf: leaf.taxon is not None and leaf.taxon.label in wanted
              for leaf in tree.leaf_node_iter()}
    n_out = sum(is_out.values())
    n_all = len(is_out)
    if n_out == 0:
        return "no outgroup taxon in the tree"
    if n_out == n_all:
        return "every taxon of the tree is in the outgroup"

    # a tree already rooted between the outgroup and the ingroup is kept as it is,
    # so that e.g. a time tree stays ultrametric
    root_children = tree.seed_node.child_nodes()
    if len(root_children) == 2:
        for child in root_children:
            flags = {is_out[leaf] for leaf in child.leaf_iter()}
            if len(flags) == 1 and (n_out if True in flags else n_all - n_out) == len(child.leaf_nodes()):
                tree.is_rooted = True
                return None

    unroot_tree(tree)

    # the outgroup must be one side of a split: a clade, or the complement of one
    out_below = {}
    all_below = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            out_below[node] = int(is_out[node])
            all_below[node] = 1
        else:
            out_below[node] = sum(out_below[c] for c in node.child_node_iter())
            all_below[node] = sum(all_below[c] for c in node.child_node_iter())

    for node in tree.postorder_node_iter():
        if node is tree.seed_node:
            continue
        outgroup_clade = out_below[node] == n_out == all_below[node]
        ingroup_clade = out_below[node] == 0 and all_below[node] == n_all - n_out
        if outgroup_clade or ingroup_clade:
            length = node.edge.length
            root_on_edge(tree, node.edge, None if length is None else float(length) / 2)
            return None
    return "the outgroup is not monophyletic"

def root_at_midpoint(tree):

    # returns None when rooted, or the reason why the tree cannot be rooted
    unroot_tree(tree)
    edges = [node.edge for node in tree.preorder_node_iter() if node is not tree.seed_node]
    if not edges:
        return "the tree has a single tip"
    if any(edge.length is None for edge in edges):
        return "midpoint rooting needs branch lengths"

    # neighbours of every node with the length of the branch between them
    neighbours = {node: [] for node in tree.preorder_node_iter()}
    for node in tree.preorder_node_iter():
        if node is not tree.seed_node:
            length = float(node.edge.length)
            neighbours[node].append((node.parent_node, length))
            neighbours[node.parent_node].append((node, length))

    def farthest_tip(start):
        distance = {start: 0.0}
        previous = {start: None}
        stack = [start]
        while stack:
            node = stack.pop()
            for other, length in neighbours[node]:
                if other not in distance:
                    distance[other] = distance[node] + length
                    previous[other] = node
                    stack.append(other)
        tips = [n for n in distance if n.is_leaf()]
        far = max(tips, key=lambda n: distance[n])
        return far, distance, previous

    # the longest tip-to-tip path: farthest tip from any tip, then farthest from it
    a, _, _ = farthest_tip(next(tree.leaf_node_iter()))
    b, distance, previous = farthest_tip(a)
    diameter = distance[b]
    if diameter <= 0:
        return "all branch lengths are zero"
    half = diameter / 2.0

    # walk from b towards a until the branch that contains the midpoint
    node = b
    while previous[node] is not None:
        parent_on_path = previous[node]
        if distance[parent_on_path] <= half <= distance[node]:
            near, far = parent_on_path, node      # near is on the side of a
            break
        node = parent_on_path

    # the branch lies between near and far; find which one is the child
    if far.parent_node is near:
        edge = far.edge
        distance_from_tail = half - distance[near]
    else:
        edge = near.edge
        distance_from_tail = distance[far] - half
    root_on_edge(tree, edge, min(max(distance_from_tail, 0.0), float(edge.length)))
    return None

def group_monophyly(tree, members, force_unrooted=False):

    # rooted trees: the group is monophyletic if the clade of its MRCA holds only
    # group members. Unrooted trees: if one branch separates the group from all
    # other tips, wherever the root is. Groups are tested on the members present.
    wanted = set(members)
    n_all = 0
    in_group = {}
    tips_below = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            n_all += 1
            in_group[node] = int(node.taxon.label in wanted)
            tips_below[node] = 1
        else:
            in_group[node] = sum(in_group[c] for c in node.child_node_iter())
            tips_below[node] = sum(tips_below[c] for c in node.child_node_iter())
    present = in_group[tree.seed_node]
    rooted = not force_unrooted and len(tree.seed_node.child_nodes()) == 2
    status = {"rooted": rooted, "present": present, "monophyletic": None,
              "intruders": None, "node": None, "side": None}

    # one member, or no tip outside the group: nothing to test
    if present < 2 or present == n_all:
        return status

    if rooted:
        # the first node in postorder that holds every member is their MRCA
        for node in tree.postorder_node_iter():
            if in_group[node] == present:
                break
        intruders = sorted(leaf.taxon.label for leaf in node.leaf_iter()
                           if leaf.taxon.label not in wanted)
        status.update(monophyletic=not intruders, intruders=intruders, node=node)
        return status

    for node in tree.postorder_node_iter():
        if node is tree.seed_node:
            continue
        if in_group[node] == present == tips_below[node]:
            status.update(monophyletic=True, intruders=[], node=node, side="below")
            return status
        if in_group[node] == 0 and tips_below[node] == n_all - present:
            status.update(monophyletic=True, intruders=[], node=node, side="above")
            return status
    status.update(monophyletic=False)
    return status

def extract_group_clade(tree, status):

    # the clade of a group as a new tree, or the reason why there is none
    if status["monophyletic"] is None:
        return "fewer than 2 group members, or no other tip, in the tree"
    if status["rooted"]:
        labels = [leaf.taxon.label for leaf in status["node"].leaf_iter()]
    elif not status["monophyletic"]:
        return "the group is not monophyletic in an unrooted tree"
    elif status["side"] == "below":
        labels = [leaf.taxon.label for leaf in status["node"].leaf_iter()]
    else:
        # the group is the rest of the tree above this node: on a copy, root on the
        # node's branch so that the group becomes one side of the root
        below = set(leaf.taxon.label for leaf in status["node"].leaf_iter())
        copy = tree.clone(depth=1)
        for node in copy.postorder_node_iter():
            if node is not copy.seed_node and \
                    set(leaf.taxon.label for leaf in node.leaf_iter()) == below:
                break
        length = node.edge.length
        root_on_edge(copy, node.edge, None if length is None else float(length) / 2)
        group_side = [child for child in copy.seed_node.child_nodes()
                      if not below & set(leaf.taxon.label for leaf in child.leaf_iter())][0]
        tree = copy
        labels = [leaf.taxon.label for leaf in group_side.leaf_iter()]
    clade = tree.extract_tree_with_taxa_labels(labels)
    clade.seed_node.edge.length = None
    clade.is_rooted = True
    return clade

def lengths_are_drawable(tree, width):

    # drawn to scale, most branches should be at least one character long;
    # one very long branch would otherwise squeeze the rest of the tree
    labels = [len(leaf.taxon.label) if leaf.taxon is not None else 10
              for leaf in tree.leaf_node_iter()]
    columns = max(width - max(labels) - 2, 10)
    depth = max(leaf.distance_from_root() for leaf in tree.leaf_node_iter())
    lengths = sorted(float(node.edge.length) for node in tree.preorder_node_iter()
                     if node is not tree.seed_node)
    if depth <= 0:
        return False
    return statistics.median(lengths) * columns / depth >= 1

def draw_tree(tree, width, metric="auto", show_supports=False):

    # a text drawing of one tree; returns the drawing and the metric used
    leaves = tree.leaf_nodes()
    if len(leaves) == 1:
        leaf = leaves[0]
        return "--- " + (leaf.taxon.label if leaf.taxon else "(no label)"), "topology"

    edges = [node.edge for node in tree.preorder_node_iter() if node is not tree.seed_node]
    has_lengths = all(edge.length is not None and edge.length > 0 for edge in edges)
    if metric == "length" and not has_lengths:
        raise TreekitError(
            "-m length needs a positive length on every branch; use -m topology")
    if metric == "auto":
        metric = "length" if has_lengths and lengths_are_drawable(tree, width) else "topology"

    def node_label(node):
        if node.is_leaf():
            return node.taxon.label if node.taxon is not None else "(no label)"
        return node.label if (show_supports and node.label is not None) else "+"

    picture = tree.as_ascii_plot(
        plot_metric = "length" if metric == "length" else "depth",
        width = width,
        show_internal_node_labels = show_supports,
        node_label_compose_fn = node_label,
    )
    lines = [line.rstrip() for line in picture.splitlines()]
    while lines and not lines[-1]:
        lines.pop()
    return "\n".join(lines), metric

def popcount(value):

    return bin(value).count("1")

class TreeSplits:
    """The branches of one tree as bit masks of the tips below them."""

    def __init__(self, shape, bit):
        self.source, self.index, parents, lengths, labels, self.rooted = shape
        self.name = "%s/%d" % (self.source, self.index)
        masks = [0] * len(parents)
        for i in range(len(parents) - 1, -1, -1):     # children come after parents
            if labels[i] is not None:
                masks[i] |= bit[labels[i]]
            if parents[i] >= 0:
                masks[parents[i]] |= masks[i]
        self.taxa = masks[0]
        self.edges = [(masks[i], lengths[i]) for i in range(1, len(parents))]
        self.has_lengths = any(length for _, length in self.edges)
        self._cache = {}

    def restricted(self, taxa, rooted):

        # the branches of the tree reduced to `taxa`: the non-trivial ones as a set
        # (for the RF distance) and all of them with summed lengths (weighted RF).
        # Reducing a tree to a subset of tips keeps exactly the branches whose tip
        # sets still split that subset, and branches that become one add up.
        key = (taxa, rooted)
        if key in self._cache:
            return self._cache[key]
        lowest = taxa & -taxa
        n_taxa = popcount(taxa)
        weights = {}
        for mask, length in self.edges:
            side = mask & taxa
            if side == 0 or side == taxa:
                continue
            if not rooted and side & lowest:
                side = taxa ^ side        # a split and its complement are one split
            weights[side] = weights.get(side, 0.0) + (length or 0.0)
        if rooted:
            nontrivial = frozenset(s for s in weights if popcount(s) >= 2)
        else:
            nontrivial = frozenset(s for s in weights
                                   if 2 <= popcount(s) <= n_taxa - 2)
        self._cache[key] = (nontrivial, weights)
        return nontrivial, weights

def clear_stem(tree):

    # after pruning, the new root may keep the length of the branch above it;
    # that is not a branch of the tree (it would count in the tree length)
    tree.seed_node.edge.length = None

def node_support_text(node, key=None):

    # the raw support of the branch above node: an annotation, or the node label
    if key is None:
        return node.label
    for annotation in node.annotations:
        if annotation.name == key:
            return annotation.value
    return None

def parse_support(text, field=None):

    # a number, or the n-th '/'-separated number of labels such as "80.5/97"
    text = str(text).strip()
    if field is not None:
        parts = text.split("/")
        if field > len(parts):
            return None
        text = parts[field - 1].strip()
    try:
        value = float(text)
    except ValueError:
        return None
    return value if math.isfinite(value) else None

def format_number(value):

    return "%g" % float(value)

def format_numbers(numbers, limit=20):

    # "1, 4, 7" for a short list; the first `limit` numbers and "..." otherwise
    text = ", ".join(map(str, numbers[:limit]))
    return text + (" ..." if len(numbers) > limit else "")

def decimal_correct(m, n=9):

    # round (not truncate) float m to n decimal places
    return round(float(m), n)

# main function

def main():

    try:
        # initialize parsed arguments and FuncTree objects
        kwargs = configs()
        setup_logging(kwargs.get("log_file"),
                      list(kwargs.get("in_files") or []) + list(kwargs.get("ref_files") or []))
        log.info("run: treekit " + " ".join(sys.argv[1:]))
        func_tree = FuncTree(**kwargs)

        # run the FuncTree method of the command
        getattr(func_tree, COMMAND_RUNNERS[func_tree.command])()
    except TreekitError as err:
        log.error(str(err))
        sys.exit(1)
    except BrokenPipeError:
        # the reader closed the pipe (e.g. "| head"): stop quietly
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        sys.exit(1)

def configs():

    # initialize parsed arguments
    config = ParsedArgs()

    # get arguments
    config_dict = config.get_args_dict()
    return config_dict

if __name__ == '__main__':

    main()
