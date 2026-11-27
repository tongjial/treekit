"""Regression tests for the issues found in the review of pull request #2."""
import pathlib
import random

import pytest

import treekit

EXAMPLES = pathlib.Path(treekit.__file__).parent / "examples"


def read(newick, **kwargs):
    tl = treekit.new_tree_list()
    tl.read(data=newick, schema="newick", case_sensitive_taxon_labels=True,
            preserve_underscores=True, **kwargs)
    return tl[0]


def labels_of(newick):
    return sorted(l.taxon.label for l in read(newick).leaf_node_iter())


# ---------------------------------------------------------------- 1. labels with spaces

class TestLabelsWithSpaces:

    @pytest.mark.parametrize("args", [
        ("remove", "-n", "D"), ("extract", "-n", "A B", "A_B", "C"), ("keep", "-n", "C"),
        ("root", "-n", "C"), ("prune", "--max-length", "100"),
    ])
    def test_a_label_with_a_space_is_not_renamed(self, run, write, args):
        path = write("sp.tre", "(('A B':1,A_B:1):1,C:1,D:1);\n")
        result = run(*args, "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "'A B'" in result.stdout
        # the output is readable, and the two taxa are still two taxa
        assert {"A B", "A_B"} <= set(labels_of(result.stdout.strip().split("\n")[-1]))

    def test_nexus_output(self, run, write):
        path = write("sp.tre", "(('A B':1,A_B:1):1,C:1,D:1);\n")
        result = run("remove", "-n", "D", "-u", "nexus", "-i", str(path))
        assert "'A B'" in result.stdout
        out = write("out.nex", result.stdout)
        assert run("taxa", "-i", str(out)).stdout.splitlines() == ["A B", "A_B", "C"]


# ---------------------------------------------------------------- 2. failed rooting

class TestFailedRootingLeavesTheTree:

    @pytest.mark.parametrize("args, newick", [
        (("-n", "A", "C"), "[&R] ((A:1,B:1):1,(C:1,D:1):2);"),
        (("-n", "Z"), "[&R] ((A:1,B:1):1,(C:1,D:1):2);"),
        (("--midpoint",), "[&R] ((A,B),(C,D));"),
        (("--midpoint",), "[&R] ((A:0,B:0):0,(C:0,D:0):0);"),
    ])
    def test_written_exactly_as_read(self, run, write, args, newick):
        path = write("r.tre", newick + "\n((A:1,B:1):1,(C:1,D:1):1,Z:1);\n")
        result = run("root", *args, "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "cannot be rooted" in result.stderr
        first = result.stdout.splitlines()[0]
        original = read(newick)
        assert first == original.as_string(schema="newick").strip()


# ---------------------------------------------------------------- 3. weighted RF and missing lengths

class TestWeightedRfMissingLengths:

    def columns(self, run, a, b):
        line = run("compare", "-i", str(a), "-r", str(b)).stdout.splitlines()[1].split("\t")
        return line[7], line[9]

    def test_one_missing_length_gives_na(self, run, write):
        a = write("a.tre", "((A:1,B:1),C:1,D:1,E:1);\n")
        b = write("b.tre", "((A:1,B:1):5,C:1,D:1,E:1);\n")
        assert self.columns(run, a, b) == ("0", "NA")

    def test_explicit_zero_lengths_are_measured(self, run, write):
        a = write("a.tre", "((A:0,B:0):0,C:0,D:0,E:0);\n")
        b = write("b.tre", "((A:1,B:1):5,C:1,D:1,E:1);\n")
        assert self.columns(run, a, b) == ("0", "10.0")

    def test_missing_length_on_a_branch_removed_by_pruning_is_ignored(self, run, write):
        # F's branch has no length, but F is not in the other tree
        a = write("a.tre", "((A:1,B:1):1,C:1,(D:1,F):1,E:1);\n")
        b = write("b.tre", "((A:1,B:1):1,C:1,D:2,E:1);\n")
        assert self.columns(run, a, b) == ("0", "0.0")


# ---------------------------------------------------------------- 4. annotations are kept

class TestAnnotationsAreKept:

    MCC = EXAMPLES / "concat_6genes.mcc.tre"

    @pytest.mark.parametrize("args", [
        ("keep", "-n", "Coris_monspeliensis_Y925"),
        ("remove", "-n", "Badula_crassa_Y24864"),
        ("extract", "-n", "Coris_monspeliensis_Y925", "Badula_crassa_Y24864",
         "Ardisia_solanacea_Y17988", "Embelia_vestita_YG205"),
        ("prune", "--max-length", "1000"),
    ])
    @pytest.mark.parametrize("schema", ["newick", "nexus"])
    def test_posterior_survives_and_feeds_collapse(self, run, tmp_path, args, schema):
        out = tmp_path / "out.tre"
        result = run(*args, "-u", schema, "-o", str(out), "-i", str(self.MCC))
        assert result.returncode == 0, result.stderr
        assert "posterior=" in out.read_text()
        collapsed = run("collapse", "-m", "0.95", "--support-key", "posterior", "-i", str(out))
        assert collapsed.returncode == 0
        assert "no 'posterior' annotation" not in collapsed.stderr

    def test_rooting_moves_annotations_with_their_split(self, run, write):
        path = write("a.tre", "((A:1,B:1)[&posterior=0.9]:1,(C:1,D:1)[&posterior=0.7]:1,E:1);\n")
        result = run("root", "-n", "A", "-i", str(path))
        assert result.returncode == 0, result.stderr
        tree = read(result.stdout.strip(), extract_comment_metadata=True)
        by_clade = {}
        for node in tree.postorder_internal_node_iter():
            for annotation in node.annotations:
                by_clade[tuple(sorted(l.taxon.label for l in node.leaf_iter()))] = annotation.value
        # rooted on A, the split {A, B} | {C, D, E} is the clade (C, D, E): it keeps 0.9;
        # the split {C, D} | {A, B, E} is the clade (C, D): it keeps 0.7
        assert by_clade.get(("C", "D", "E"), by_clade.get(("A", "B"))) == "0.9"
        assert by_clade.get(("C", "D")) == "0.7"
        assert sorted(by_clade.values()) == ["0.7", "0.9"]


# ---------------------------------------------------------------- 5. compare cache stays small

class TestCompareCache:

    def test_only_full_restrictions_are_kept(self):
        rng = random.Random(1)
        labels = ["t%d" % i for i in range(40)]

        def tree():
            members = [l for l in labels if rng.random() < 0.9]
            nodes = ["%s:1" % l for l in members]
            while len(nodes) > 3:
                a = nodes.pop(rng.randrange(len(nodes)))
                b = nodes.pop(rng.randrange(len(nodes)))
                nodes.append("(%s,%s):1" % (a, b))
            return read("(%s);" % ",".join(nodes))

        bit = {l: 1 << i for i, l in enumerate(labels)}
        splits = []
        for t in [tree() for _ in range(30)]:
            nodes = list(t.preorder_node_iter())
            pos = {n: i for i, n in enumerate(nodes)}
            shape = ("x", 1, [-1 if n is t.seed_node else pos[n.parent_node] for n in nodes],
                     [None if n is t.seed_node else 1.0 for n in nodes],
                     [n.taxon.label if n.is_leaf() else None for n in nodes], False)
            splits.append(treekit.TreeSplits(shape, bit))
        for a in splits:
            for b in splits:
                common = a.taxa & b.taxa
                a.restricted(common, False)
                b.restricted(common, False)
        # at most one stored restriction per tree, whatever the number of pairs
        assert all(not hasattr(s, "_cache") for s in splits)
        assert sum(s._full is not None for s in splits) <= len(splits)


# ---------------------------------------------------------------- 6. a root that does not move

class TestUnmovedRootKeepsLabels:

    def test_outgroup_already_at_the_root(self, run, write):
        newick = "[&R] ((A:1,B:1)90:1,(C:1,D:1)70:2);"
        path = write("l.tre", newick + "\n")
        result = run("root", "-n", "C", "D", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stderr == ""
        assert result.stdout.strip() == read(newick).as_string(schema="newick").strip()

    def test_unrooted_input_with_unroot(self, run, write):
        newick = "((A:1,B:1)90:1,(C:1,D:1)70:1,E:1);"
        path = write("u.tre", newick + "\n")
        result = run("root", "--unroot", "-i", str(path))
        assert "90" in result.stdout and "70" in result.stdout

    def test_different_labels_merged_at_the_old_root_are_reported(self, run, write):
        path = write("l.tre", "[&R] ((A:1,B:1)90:1,(C:1,D:1)70:2);\n")
        result = run("root", "-n", "A", "-i", str(path))
        assert result.returncode == 0
        assert "different labels or annotations on the two branches at the old root" in result.stderr
        unrooted = run("root", "--unroot", "-i", str(path))
        assert "different labels or annotations on the two branches at the old root" in unrooted.stderr

    def test_equal_labels_at_the_old_root_are_not_reported(self, run, write):
        path = write("l.tre", "[&R] ((A:1,B:1)90:1,(C:1,D:1)90:2);\n")
        result = run("root", "-n", "A", "-i", str(path))
        assert result.stderr == ""


# ---------------------------------------------------------------- 7 + 8. output and log destinations

class TestWriteTargets:

    TREE = "((A:1,B:1):1,C:1,D:1);\n"

    def test_output_and_log_cannot_be_the_same_file(self, run, write, tmp_path):
        path = write("t.tre", self.TREE)
        report = tmp_path / "report.tsv"
        result = run("summary", "-i", str(path), "-o", str(report), "-l", "./report.tsv")
        assert result.returncode == 1
        assert "treekit: error: -o and -l name the same file" in result.stderr
        assert not report.exists()

    @pytest.mark.parametrize("flag", ["-o", "-l"])
    def test_groups_file_is_protected(self, run, write, flag):
        path = write("t.tre", self.TREE)
        groups = write("groups.tsv", "A\tg\nB\tg\n")
        result = run("mrca", "-g", "groups.tsv", "-i", str(path), flag, str(groups))
        assert result.returncode == 1
        assert "treekit: error: %s file" % flag in result.stderr
        assert groups.read_text() == "A\tg\nB\tg\n"

    @pytest.mark.parametrize("command", ["remove", "extract", "keep", "root"])
    @pytest.mark.parametrize("flag", ["-o", "-l"])
    def test_taxa_file_is_protected(self, run, write, command, flag):
        path = write("t.tre", self.TREE)
        taxa = write("taxa.txt", "A\n")
        result = run(command, "-f", "taxa.txt", "-i", str(path), flag, "./taxa.txt")
        assert result.returncode == 1
        assert taxa.read_text() == "A\n"
