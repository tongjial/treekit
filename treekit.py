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
    ("mrca", "check the monophyly or extract the clade based on given taxa/group", None),
    ("extract", "generate subtrees only consist of the given taxa", "write_extracted_trees"),
    ("prune", "prune trees with specific branchlength", None),
    ("root", "root/reroot/unroot trees", "write_rooted_trees"),
    ("draw", "display trees as textplot", None),
    ("compare", "compare trees between two different sources", None),
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

        # never write the output over one of the input files
        if self.out_file is not None and self.out_file in self.in_files:
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

    def map_files(self, function):

        # apply function(in_file, in_format) to every input file, one process per file
        jobs = list(zip(self.in_files, self.in_formats))

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

    def read_all_trees(self):

        # every file is parsed on its own, because a NEXUS file declares the number
        # of its taxa and cannot be read into a namespace already filled by another
        # file; the trees are then merged into one shared namespace
        tree_list = new_tree_list()
        for in_file, in_format in zip(self.in_files, self.in_formats):
            file_trees = self.read_trees(
                in_file,
                in_format,
                preserve_underscores = True
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
            unquoted_underscores=True
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
        setup_logging(kwargs.get("log_file"), kwargs.get("in_files"))
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
