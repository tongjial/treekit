import random

import dendropy
import pytest

# D is in two of the three trees
MULTI = (
    "((A:1,B:1):1,(C:1,D:1):1);\n"
    "((A:1,C:1):1,(B:1,E:1):1);\n"
    "((B:1,D:1):1,(C:1,E:1):1);\n"
)
SINGLE = "((A:1,B:1):1,(C:1,D:1):1);\n"


def tree_lines(stdout):
    return [line for line in stdout.splitlines() if line.strip()]


def leaf_set(newick):
    tree = dendropy.Tree.get(data=newick, schema="newick", preserve_underscores=True)
    return {leaf.taxon.label for leaf in tree.leaf_node_iter()}


def random_tree(labels, rng):
    nodes = ["%s:%.3f" % (label, rng.uniform(0.1, 2)) for label in labels]
    while len(nodes) > 2:
        a = nodes.pop(rng.randrange(len(nodes)))
        b = nodes.pop(rng.randrange(len(nodes)))
        nodes.append("(%s,%s):%.3f" % (a, b, rng.uniform(0.1, 2)))
    return "(%s,%s);" % tuple(nodes)


# ---------------------------------------------------------------- extract

class TestExtract:

    def test_single_tree(self, run, write):
        path = write("s.tre", SINGLE)
        result = run("extract", "-n", "A", "B", "C", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stderr == ""
        assert result.stdout.strip() == "((A:1.0,B:1.0):1.0,C:2.0);"

    def test_branch_lengths_are_merged(self, run, write):
        path = write("s.tre", SINGLE)
        result = run("extract", "-n", "A", "C", "-i", str(path))
        assert result.stdout.strip() == "(A:2.0,C:2.0);"

    def test_multi_tree_partial_presence_is_logged_not_warned(self, run, write, tmp_path):
        path = write("m.tre", MULTI)
        log = tmp_path / "run.log"
        result = run("extract", "-n", "B", "C", "D", "-l", str(log), "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stderr == ""
        lines = tree_lines(result.stdout)
        assert len(lines) == 3
        assert [leaf_set(line) for line in lines] == [{"B", "C", "D"}, {"B", "C"}, {"B", "C", "D"}]
        assert "taxon 'D' extracted from 2 of 3 trees" in log.read_text()

    def test_absent_taxa_warn_once(self, run, write):
        path = write("m.tre", MULTI)
        result = run("extract", "-n", "B", "C", "ZZZ", "YYY", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "warning: 2 taxa were not found in any input tree: ZZZ, YYY" in result.stderr
        assert result.stderr.count("warning") == 1

    def test_strict(self, run, write):
        path = write("m.tre", MULTI)
        result = run("extract", "-n", "B", "C", "ZZZ", "--strict", "-i", str(path))
        assert result.returncode == 1
        assert "1 taxon was not found in any input tree: ZZZ (--strict)" in result.stderr
        assert result.stdout == ""

    def test_all_absent_is_an_error(self, run, write):
        path = write("s.tre", SINGLE)
        result = run("extract", "-n", "ZZZ", "-i", str(path))
        assert result.returncode == 1
        assert "none of the 1 taxon to extract were found in any input tree" in result.stderr
        assert "treekit keep" not in result.stderr

    def test_tree_without_any_requested_taxon_is_an_error(self, run, write):
        path = write("m.tre", "((A:1,B:1):1,C:2);\n((X:1,Y:1):1,Z:2);\n")
        result = run("extract", "-n", "A", "B", "-i", str(path))
        assert result.returncode == 1
        assert "1 tree(s) contain none of the taxa to extract (tree 2)" in result.stderr
        assert "treekit keep" in result.stderr
        assert result.stdout == ""

    def test_single_tip_warning(self, run, write):
        path = write("m.tre", "((A:1,B:1):1,C:2);\n((A:1,X:1):1,Y:2);\n")
        result = run("extract", "-n", "A", "B", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "extracting these taxa leaves 1 tree(s) with a single tip (tree 2)" in result.stderr

    def test_unlabelled_tips_are_dropped(self, run, write):
        path = write("u.tre", "((A:1,(:1,B:1):1):1,(C:1,D:1):1);\n")
        result = run("extract", "-n", "A", "B", "C", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert leaf_set(result.stdout.strip()) == {"A", "B", "C"}

    def test_taxa_file(self, run, write):
        path = write("s.tre", SINGLE)
        names = write("keep.txt", "A\n\nB\n")
        result = run("extract", "-f", str(names), "-i", str(path))
        # no stem: the branch above the new root is not a branch of the tree
        assert result.stdout.strip() == "(A:1.0,B:1.0);"

    def test_nexus_output_declares_only_extracted_taxa(self, run, write):
        path = write("m.tre", MULTI)
        result = run("extract", "-n", "B", "C", "D", "-u", "nexus", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "NTAX=3" in result.stdout
        taxlabels = result.stdout.split("TAXLABELS")[1].split(";")[0].split()
        assert sorted(taxlabels) == ["B", "C", "D"]

    def test_output_can_be_read_back(self, run, write, tmp_path):
        path = write("m.tre", MULTI)
        out = tmp_path / "out.tre"
        run("extract", "-n", "B", "C", "D", "-o", str(out), "-i", str(path))
        again = run("summary", "-i", str(out))
        assert again.returncode == 0, again.stderr

    def test_same_result_as_removing_the_other_taxa(self, run, write):
        # extract S and remove (all taxa - S) must give the same trees
        rng = random.Random(5)
        labels = ["t%d" % i for i in range(12)]
        trees = "".join(random_tree(labels, rng) + "\n" for _ in range(20))
        path = write("r.tre", trees)
        keep = labels[:7]
        extracted = run("extract", "-n", *keep, "-i", str(path))
        removed = run("remove", "-n", *labels[7:], "-i", str(path))
        assert extracted.returncode == 0 and removed.returncode == 0
        assert extracted.stdout == removed.stdout

    def test_input_files_are_not_changed(self, run, write):
        path = write("s.tre", SINGLE)
        run("extract", "-n", "A", "B", "-i", str(path))
        assert path.read_text() == SINGLE


# ---------------------------------------------------------------- keep

class TestKeep:

    def test_keeps_trees_with_all_taxa_in_input_order(self, run, write):
        # trees 1 and 3 contain both B and D
        path = write("m.tre", MULTI)
        result = run("keep", "-n", "B", "D", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stderr == ""
        assert tree_lines(result.stdout) == [
            "((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0);",
            "((B:1.0,D:1.0):1.0,(C:1.0,E:1.0):1.0);",
        ]

    def test_kept_trees_are_unchanged(self, run, write):
        # every tree has A..E? no: only keep the ones with all; their tips stay
        path = write("m.tre", MULTI)
        result = run("keep", "-n", "B", "-i", str(path))
        lines = tree_lines(result.stdout)
        assert len(lines) == 3
        assert [leaf_set(line) for line in lines] == [
            {"A", "B", "C", "D"}, {"A", "B", "C", "E"}, {"B", "C", "D", "E"}]

    def test_single_tree_input(self, run, write):
        path = write("s.tre", SINGLE)
        result = run("keep", "-n", "A", "D", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert len(tree_lines(result.stdout)) == 1

    def test_no_tree_contains_all_taxa_is_an_error(self, run, write):
        path = write("m.tre", MULTI)
        result = run("keep", "-n", "A", "D", "E", "-i", str(path))
        assert result.returncode == 1
        assert "none of the 3 input trees contains all 3 requested taxa" in result.stderr
        assert result.stdout == ""

    def test_taxon_absent_from_all_trees_is_an_error(self, run, write):
        path = write("m.tre", MULTI)
        result = run("keep", "-n", "B", "ZZZ", "YYY", "-i", str(path))
        assert result.returncode == 1
        assert "no tree can be kept: 2 taxa are not in any input tree: ZZZ, YYY" in result.stderr

    def test_log_records_the_kept_trees(self, run, write, tmp_path):
        path = write("m.tre", MULTI)
        log = tmp_path / "run.log"
        run("keep", "-n", "B", "D", "-l", str(log), "-i", str(path))
        text = log.read_text()
        assert "kept 2 of 3 trees (tree 1, 3)" in text
        assert "taxon 'D' is in 2 of 3 trees" in text

    def test_several_files(self, run, write):
        write("a.tre", "((A:1,B:1):1,C:2);\n")
        write("b.tre", "((A:1,X:1):1,Y:2);\n((A:1,B:1):1,Z:2);\n")
        result = run("keep", "-n", "A", "B", "-i", "a.tre", "b.tre")
        assert result.returncode == 0, result.stderr
        assert len(tree_lines(result.stdout)) == 2

    def test_nexus_output_declares_only_taxa_of_kept_trees(self, run, write):
        path = write("m.tre", MULTI)
        result = run("keep", "-n", "A", "C", "E", "-u", "nexus", "-i", str(path))
        assert result.returncode == 0, result.stderr
        taxlabels = result.stdout.split("TAXLABELS")[1].split(";")[0].split()
        assert sorted(taxlabels) == ["A", "B", "C", "E"]

    def test_taxa_file(self, run, write):
        path = write("m.tre", MULTI)
        names = write("need.txt", "B\nD\n")
        result = run("keep", "-f", str(names), "-i", str(path))
        assert len(tree_lines(result.stdout)) == 2

    def test_keep_then_extract(self, run, write, tmp_path):
        # the error message of extract points to keep; the two work together
        path = write("m.tre", "((A:1,B:1):1,C:2);\n((X:1,Y:1):1,Z:2);\n")
        kept = tmp_path / "kept.tre"
        run("keep", "-n", "A", "B", "-o", str(kept), "-i", str(path))
        result = run("extract", "-n", "A", "B", "-i", str(kept))
        assert result.returncode == 0, result.stderr
        assert result.stdout.strip() == "(A:1.0,B:1.0);"


class TestNoStemAfterPruning:

    @pytest.mark.parametrize("args", [("remove", "-n", "C"), ("extract", "-n", "A", "B")])
    def test_tree_length_is_not_inflated(self, run, write, tmp_path, args):
        # removing the whole side C of the root used to leave "(A:1,B:1):1;",
        # and summary then counted that stem: a tree length of 3 instead of 2
        path = write("s.tre", "((A:1,B:1):1,C:2);\n")
        out = tmp_path / "out.tre"
        run(*args, "-o", str(out), "-i", str(path))
        assert out.read_text().strip() == "(A:1.0,B:1.0);"
        detail = run("summary", "-e", "-i", str(out)).stdout.splitlines()[1].split("\t")
        assert detail[6] == "2.0"
