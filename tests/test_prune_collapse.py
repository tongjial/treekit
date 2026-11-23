import pathlib
import re

import pytest

import treekit

EXAMPLES = pathlib.Path(treekit.__file__).parent / "examples"
# D's branch (5) is 50 times the median branch length (0.1)
LONG_D = "((A:0.1,B:0.1):0.1,(C:0.1,D:5):0.1,E:0.1);\n"
EVEN = "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1,E:0.1);\n"
NO_LENGTHS = "((A,B),(C,D),E);\n"


def read(newick, **kwargs):
    tl = treekit.new_tree_list()
    tl.read(data=newick, schema="newick", case_sensitive_taxon_labels=True,
            preserve_underscores=True, **kwargs)
    return tl[0]


def tips(newick):
    return sorted(l.taxon.label for l in read(newick).leaf_node_iter() if l.taxon)


def tree_lines(stdout):
    return [line for line in stdout.splitlines() if line.strip()]


def root_to_tip(newick):
    tree = read(newick)
    return {l.taxon.label: round(l.distance_from_root(), 9) for l in tree.leaf_node_iter()}


# ---------------------------------------------------------------- prune

class TestPrune:

    def test_ratio_removes_the_long_tip(self, run, write):
        path = write("p.tre", LONG_D)
        result = run("prune", "--max-ratio", "10", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stderr == ""
        assert result.stdout.strip() == "((A:0.1,B:0.1):0.1,C:0.2,E:0.1);"

    def test_ratio_uses_the_median_of_each_tree(self, run, write):
        # in tree 2 every branch is 10, so 10 is not long there
        path = write("p.tre", LONG_D + "((A:10,B:10):10,(C:10,D:10):10,E:10);\n")
        result = run("prune", "--max-ratio", "10", "-i", str(path))
        first, second = tree_lines(result.stdout)
        assert tips(first) == ["A", "B", "C", "E"]
        assert tips(second) == ["A", "B", "C", "D", "E"]

    def test_absolute_length(self, run, write):
        path = write("p.tre", LONG_D + "((A:2,B:0.1):0.1,(C:0.1,D:0.1):0.1,E:0.1);\n")
        result = run("prune", "--max-length", "1", "-i", str(path))
        first, second = tree_lines(result.stdout)
        assert tips(first) == ["A", "B", "C", "E"]
        assert tips(second) == ["B", "C", "D", "E"]

    def test_both_thresholds_remove_tips_longer_than_either(self, run, write):
        # B (0.5) is over the ratio limit 5 x 0.1; D (5) is over both
        path = write("p.tre", "((A:0.1,B:0.5):0.1,(C:0.1,D:5):0.1,E:0.1,F:0.1);\n")
        only_length = run("prune", "--max-length", "1", "-i", str(path))
        both = run("prune", "--max-length", "1", "--max-ratio", "4", "-i", str(path))
        assert tips(only_length.stdout.strip()) == ["A", "B", "C", "E", "F"]
        assert tips(both.stdout.strip()) == ["A", "C", "E", "F"]

    def test_internal_branches_are_not_pruned(self, run, write):
        # the long branch above (C, D) is internal: its tips are kept
        path = write("p.tre", "((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):9,E:0.1);\n")
        result = run("prune", "--max-length", "1", "-i", str(path))
        assert tips(result.stdout.strip()) == ["A", "B", "C", "D", "E"]

    def test_log_lists_every_removed_tip(self, run, write, tmp_path):
        path = write("p.tre", LONG_D + EVEN)
        log = tmp_path / "run.log"
        run("prune", "--max-ratio", "10", "-l", str(log), "-i", str(path))
        text = log.read_text()
        assert "tree 1: removed D (branch 5, 50 x the median 0.1)" in text
        assert "removed 1 tip(s) from 1 of 2 trees" in text

    def test_nothing_to_prune_warns(self, run, write):
        path = write("p.tre", EVEN)
        result = run("prune", "--max-ratio", "10", "-i", str(path))
        assert result.returncode == 0
        assert "no tip has a branch longer than the threshold" in result.stderr
        assert tips(result.stdout.strip()) == ["A", "B", "C", "D", "E"]

    def test_tree_without_lengths_is_written_unchanged(self, run, write):
        path = write("p.tre", LONG_D + NO_LENGTHS)
        result = run("prune", "--max-ratio", "10", "-i", str(path))
        assert result.returncode == 0
        assert "1 tree(s) cannot be pruned (no branch lengths): tree 2; written unchanged" in result.stderr
        assert len(tree_lines(result.stdout)) == 2

    def test_strict(self, run, write):
        path = write("p.tre", LONG_D + NO_LENGTHS)
        result = run("prune", "--max-ratio", "10", "--strict", "-i", str(path))
        assert result.returncode == 1
        assert result.stdout == ""

    def test_zero_median_with_ratio_only(self, run, write):
        path = write("z.tre", "((A:0,B:0):0,(C:0,D:3):0,E:0);\n")
        result = run("prune", "--max-ratio", "10", "-i", str(path))
        assert "cannot be pruned (the median branch length is 0)" in result.stderr
        # with an absolute threshold as well, that one is used
        both = run("prune", "--max-ratio", "10", "--max-length", "1", "-i", str(path))
        assert tips(both.stdout.strip()) == ["A", "B", "C", "E"]

    def test_every_tip_too_long_is_an_error(self, run, write):
        path = write("p.tre", "((A:5,B:5):1,C:5);\n")
        result = run("prune", "--max-length", "1", "-i", str(path))
        assert result.returncode == 1
        assert "every tip of tree 1 has a branch longer than 1" in result.stderr

    def test_single_tip_left_warns(self, run, write):
        path = write("p.tre", "((A:5,B:5):1,C:0.5);\n")
        result = run("prune", "--max-length", "1", "-i", str(path))
        assert result.returncode == 0
        assert "pruning leaves 1 tree(s) with a single tip (tree 1)" in result.stderr

    def test_unlabelled_tips_are_kept_unless_long(self, run, write):
        path = write("u.tre", "((A:0.1,:0.1):0.1,(C:0.1,D:5):0.1,E:0.1);\n")
        result = run("prune", "--max-length", "1", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert len(read(result.stdout.strip()).leaf_nodes()) == 4

    def test_tip_on_the_root_of_a_rooted_tree(self, run, write):
        path = write("r.tre", "[&R] ((A:1,B:1):1,C:9);\n")
        result = run("prune", "--max-length", "5", "-i", str(path))
        assert result.stdout.strip().endswith("(A:1.0,B:1.0);")

    @pytest.mark.parametrize("args", [(), ("--max-length", "0"), ("--max-ratio", "-1"),
                                      ("--max-length", "x")])
    def test_invalid_thresholds(self, run, write, args):
        path = write("p.tre", LONG_D)
        result = run("prune", *args, "-i", str(path))
        assert result.returncode in (1, 2)
        assert "Traceback" not in result.stderr

    def test_nexus_output_drops_removed_taxa(self, run, write):
        path = write("p.tre", LONG_D)
        result = run("prune", "--max-ratio", "10", "-u", "nexus", "-i", str(path))
        assert "NTAX=4" in result.stdout


# ---------------------------------------------------------------- collapse

class TestCollapse:

    LABELLED = "((A:1,B:1)40:1,(C:1,D:1)90:1,E:1);\n"

    def test_low_support_branch_is_collapsed(self, run, write):
        path = write("s.tre", self.LABELLED)
        result = run("collapse", "-m", "50", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stderr == ""
        assert result.stdout.strip() == "(A:2.0,B:2.0,(C:1.0,D:1.0)90:1.0,E:1.0);"

    def test_root_to_tip_distances_are_kept(self, run, write):
        path = write("s.tre", "[&R] (((A:1,B:1)40:1,C:2)30:1,(D:1,E:1)90:2);\n")
        before = root_to_tip(path.read_text().strip())
        result = run("collapse", "-m", "50", "-i", str(path))
        assert root_to_tip(result.stdout.strip()) == before

    def test_threshold_is_strict_lower_bound(self, run, write):
        # a support equal to the threshold is kept
        path = write("s.tre", self.LABELLED)
        result = run("collapse", "-m", "40", "-i", str(path))
        assert "40" in result.stdout

    def test_nested_low_branches(self, run, write):
        path = write("s.tre", "(((A:1,B:1)10:1,C:1)20:1,(D:1,E:1)90:1,F:1);\n")
        result = run("collapse", "-m", "50", "-i", str(path))
        tree = read(result.stdout.strip())
        assert len(tree.seed_node.child_nodes()) == 5      # A, B, C, (D,E), F

    def test_tips_and_root_are_never_collapsed(self, run, write):
        path = write("s.tre", "((A:1,B:1)90:1,(C:1,D:1)90:1)5;\n")
        result = run("collapse", "-m", "50", "-i", str(path))
        assert tips(result.stdout.strip()) == ["A", "B", "C", "D"]
        assert result.stdout.count("90") == 2

    def test_proportions(self, run, write):
        path = write("s.tre", "((A:1,B:1)0.4:1,(C:1,D:1)0.99:1,E:1);\n")
        result = run("collapse", "-m", "0.5", "-i", str(path))
        assert result.stderr == ""
        assert len(read(result.stdout.strip()).seed_node.child_nodes()) == 4

    def test_scale_mismatch_warnings(self, run, write):
        percent = write("p.tre", self.LABELLED)
        proportion = write("q.tre", "((A:1,B:1)0.4:1,(C:1,D:1)0.99:1,E:1);\n")
        assert "did you mean 50?" in run("collapse", "-m", "0.5", "-i", str(percent)).stderr
        assert "supports are between 0 and 1 but --min-support is 50" in \
            run("collapse", "-m", "50", "-i", str(proportion)).stderr

    def test_support_field(self, run, write):
        # IQ-TREE style SH-aLRT/UFBoot labels
        path = write("iq.tre", "((A:1,B:1)85.3/40:1,(C:1,D:1)20/99:1,E:1);\n")
        first = run("collapse", "-m", "50", "--support-field", "1", "-i", str(path))
        second = run("collapse", "-m", "50", "--support-field", "2", "-i", str(path))
        assert "85.3/40" in first.stdout and "20/99" not in first.stdout
        assert "20/99" in second.stdout and "85.3/40" not in second.stdout

    def test_slash_labels_without_field(self, run, write):
        path = write("iq.tre", "((A:1,B:1)85.3/40:1,(C:1,D:1)20/99:1,E:1);\n")
        result = run("collapse", "-m", "50", "-i", str(path))
        assert result.returncode == 0
        assert "use --support-field" in result.stderr
        assert "85.3/40" in result.stdout and "20/99" in result.stdout

    def test_support_key_reads_annotations_and_keeps_them(self, run, write):
        path = write("b.tre", "((A:1,B:1)[&posterior=0.4]:1,(C:1,D:1)[&posterior=0.99]:1,E:1);\n")
        result = run("collapse", "-m", "0.95", "--support-key", "posterior", "-i", str(path))
        assert result.returncode == 0, result.stderr
        out = result.stdout.strip()
        assert "posterior=0.4" not in out
        assert "posterior=0.99" in out
        assert len(read(out).seed_node.child_nodes()) == 4

    def test_missing_annotation(self, run, write):
        path = write("s.tre", self.LABELLED)
        result = run("collapse", "-m", "50", "--support-key", "posterior", "-i", str(path))
        assert "1 tree(s) have no 'posterior' annotation: tree 1; written unchanged" in result.stderr
        strict = run("collapse", "-m", "50", "--support-key", "posterior", "--strict", "-i", str(path))
        assert strict.returncode == 1

    def test_tree_without_supports_is_written_unchanged(self, run, write):
        path = write("s.tre", self.LABELLED + "((A:1,B:1):1,(C:1,D:1):1,E:1);\n")
        result = run("collapse", "-m", "50", "-i", str(path))
        assert "1 tree(s) have no numeric support in the node labels: tree 2" in result.stderr
        assert len(tree_lines(result.stdout)) == 2

    def test_raxml_example(self, run, tmp_path):
        # one branch of the example has a support below 70 (58)
        out = tmp_path / "c.tre"
        log = tmp_path / "run.log"
        source = EXAMPLES / "RAxML_bipartitions.all"
        result = run("collapse", "-m", "70", "-o", str(out), "-l", str(log), "-i", str(source))
        assert result.returncode == 0, result.stderr
        assert "collapsed 1 branch(es) in 1 of 1 trees" in log.read_text()
        assert "58" not in re.findall(r"\)([0-9.]+):", out.read_text())

    def test_beast_mcc_example_stays_ultrametric(self, run, tmp_path):
        out = tmp_path / "m.tre"
        log = tmp_path / "run.log"
        source = EXAMPLES / "concat_6genes.mcc.tre"
        result = run("collapse", "-m", "0.95", "--support-key", "posterior",
                     "-o", str(out), "-l", str(log), "-i", str(source))
        assert result.returncode == 0, result.stderr
        assert "collapsed 4 branch(es)" in log.read_text()
        summary = run("summary", "-i", str(out)).stdout.splitlines()[1].split("\t")
        assert summary[2] == "0"          # no longer binary
        assert summary[5] == "1"          # still ultrametric

    @pytest.mark.parametrize("args", [(), ("-m", "x"), ("-m", "-1"),
                                      ("-m", "50", "--support-key", "p", "--support-field", "1")])
    def test_invalid_options(self, run, write, args):
        path = write("s.tre", self.LABELLED)
        result = run("collapse", *args, "-i", str(path))
        assert result.returncode == 2
        assert "Traceback" not in result.stderr
