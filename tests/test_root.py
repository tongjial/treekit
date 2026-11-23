import random

import dendropy
import pytest

import treekit

UNROOTED = "((A:1,B:1):1,(C:1,D:1):1,E:4);\n"


def read(newick):
    tl = treekit.new_tree_list()
    tl.read(data=newick, schema="newick", case_sensitive_taxon_labels=True,
            preserve_underscores=True)
    return tl[0]


def tree_lines(stdout):
    return [line for line in stdout.splitlines() if line.strip()]


def root_sides(newick):
    """The tip sets on the two sides of the root, and the two root branch lengths."""
    tree = read(newick)
    kids = tree.seed_node.child_nodes()
    sides = sorted(tuple(sorted(l.taxon.label for l in k.leaf_iter())) for k in kids)
    lengths = [None if k.edge.length is None else round(float(k.edge.length), 9) for k in kids]
    lengths.sort(key=lambda value: (value is None, value or 0))
    return sides, lengths


def split_labels(newick):
    """Node labels by split, independent of where the root is."""
    tree = read(newick)
    full = frozenset(l.taxon.label for l in tree.leaf_node_iter())
    labels = {}
    for node in tree.postorder_internal_node_iter():
        if node is tree.seed_node or node.label is None:
            continue
        side = frozenset(l.taxon.label for l in node.leaf_iter())
        key = min(side, full - side, key=lambda s: sorted(s))
        labels[tuple(sorted(key))] = node.label
    return labels


def random_unrooted(rng, n):
    labels = ["t%d" % i for i in range(n)]
    nodes = ["%s:%.4f" % (l, rng.uniform(0.1, 3)) for l in labels]
    while len(nodes) > 3:
        a = nodes.pop(rng.randrange(len(nodes)))
        b = nodes.pop(rng.randrange(len(nodes)))
        nodes.append("(%s,%s):%.4f" % (a, b, rng.uniform(0.1, 3)))
    return "(%s,%s,%s);" % tuple(nodes)


# ---------------------------------------------------------------- outgroup rooting

class TestOutgroup:

    def test_single_taxon_outgroup(self, run, write):
        path = write("u.tre", UNROOTED)
        result = run("root", "-n", "E", "-i", str(path))
        assert result.returncode == 0, result.stderr
        out = result.stdout.strip()
        assert out.startswith("[&R]")
        # the root is in the middle of E's branch (length 4)
        assert root_sides(out) == ([("A", "B", "C", "D"), ("E",)], [2.0, 2.0])

    def test_clade_outgroup(self, run, write):
        path = write("u.tre", UNROOTED)
        result = run("root", "-n", "C", "D", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert root_sides(result.stdout.strip()) == ([("A", "B", "E"), ("C", "D")], [0.5, 0.5])

    def test_outgroup_that_is_the_complement_of_a_clade(self, run, write):
        # {A, B, C, D} is not a clade below the seed node here, but {E} is its complement
        path = write("u.tre", UNROOTED)
        result = run("root", "-n", "A", "B", "C", "D", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert root_sides(result.stdout.strip())[0] == [("A", "B", "C", "D"), ("E",)]

    def test_rerooting_a_rooted_tree_merges_the_old_root_branches(self, run, write):
        # the old root joined two branches of lengths 1 and 2 into one branch of length 3
        path = write("r.tre", "[&R] ((A:1,B:1):1,(C:1,D:1):2);\n")
        result = run("root", "-n", "C", "-i", str(path))
        assert result.returncode == 0, result.stderr
        out = result.stdout.strip()
        assert root_sides(out) == ([("A", "B", "D"), ("C",)], [0.5, 0.5])
        tree = read(out)
        ab = tree.mrca(taxon_labels=["A", "B"])
        assert float(ab.edge.length) == 3.0

    def test_tree_already_rooted_on_the_outgroup_is_unchanged(self, run, write):
        # an ultrametric time tree must stay ultrametric
        newick = "[&R] ((A:1,B:1):2,(C:2,D:2):1);"
        path = write("t.tre", newick + "\n")
        result = run("root", "-n", "C", "D", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert root_sides(result.stdout.strip()) == root_sides(newick)
        summary = run("summary", "-i", str(path))
        assert summary.stdout.splitlines()[1].split("\t")[5] == "1"     # num_ultrametric

    def test_no_branch_lengths(self, run, write):
        path = write("n.tre", "((A,B),(C,D),E);\n")
        result = run("root", "-n", "E", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert root_sides(result.stdout.strip()) == ([("A", "B", "C", "D"), ("E",)], [None, None])

    def test_partial_outgroup_uses_the_members_present(self, run, write, tmp_path):
        # tree 2 has only D of the outgroup {D, E}
        path = write("m.tre", "((A:1,B:1):1,(D:1,E:1):1,C:1);\n((A:1,B:1):1,(C:1,D:1):1,X:1);\n")
        log = tmp_path / "run.log"
        result = run("root", "-n", "D", "E", "-l", str(log), "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stderr == ""
        first, second = tree_lines(result.stdout)
        assert root_sides(first)[0] == [("A", "B", "C"), ("D", "E")]
        assert root_sides(second)[0] == [("A", "B", "C", "X"), ("D",)]
        text = log.read_text()
        assert "taxon 'E' used in 1 of 2 trees" in text
        assert "rooted 2 of 2 trees (outgroup)" in text

    def test_tree_without_outgroup_is_written_unchanged(self, run, write):
        path = write("m.tre", "((A:1,B:1):1,(C:1,D:1):1,E:1);\n((A:1,B:1):1,(C:1,X:1):1,Y:1);\n")
        result = run("root", "-n", "E", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert ("1 tree(s) cannot be rooted (no outgroup taxon in the tree): tree 2; "
                "written unchanged") in result.stderr
        lines = tree_lines(result.stdout)
        assert len(lines) == 2                      # same number and order as the input
        assert lines[0].startswith("[&R]")
        assert root_sides(lines[1])[0] == root_sides("((A:1,B:1):1,(C:1,X:1):1,Y:1);")[0]

    def test_non_monophyletic_outgroup_is_written_unchanged(self, run, write):
        path = write("u.tre", UNROOTED)
        result = run("root", "-n", "A", "C", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "cannot be rooted (the outgroup is not monophyletic): tree 1" in result.stderr

    def test_strict_stops_on_a_tree_that_cannot_be_rooted(self, run, write):
        path = write("u.tre", UNROOTED)
        result = run("root", "-n", "A", "C", "--strict", "-i", str(path))
        assert result.returncode == 1
        assert "(the outgroup is not monophyletic): tree 1 (--strict)" in result.stderr
        assert result.stdout == ""

    def test_reasons_are_reported_separately(self, run, write):
        # trees 1 and 3: {A, C, E} is not a clade; tree 2 has none of them
        path = write("m.tre", UNROOTED + "((X:1,Y:1):1,(Z:1,W:1):1,V:1);\n" + UNROOTED)
        result = run("root", "-n", "A", "C", "E", "-i", str(path))
        assert "(the outgroup is not monophyletic): tree 1, 3;" in result.stderr
        assert "(no outgroup taxon in the tree): tree 2;" in result.stderr
        assert result.stderr.count("cannot be rooted") == 2
        assert len(tree_lines(result.stdout)) == 3

    def test_outgroup_with_every_taxon(self, run, write):
        path = write("s.tre", "((A:1,B:1):1,C:1);\n")
        result = run("root", "-n", "A", "B", "C", "-i", str(path))
        assert "every taxon of the tree is in the outgroup" in result.stderr

    def test_outgroup_absent_from_all_trees_is_an_error(self, run, write):
        path = write("u.tre", UNROOTED)
        result = run("root", "-n", "ZZZ", "-i", str(path))
        assert result.returncode == 1
        assert "none of the 1 taxon of the outgroup were found in any input tree" in result.stderr

    def test_outgroup_from_file(self, run, write):
        path = write("u.tre", UNROOTED)
        names = write("out.txt", "C\nD\n")
        result = run("root", "-f", str(names), "-i", str(path))
        assert root_sides(result.stdout.strip())[0] == [("A", "B", "E"), ("C", "D")]

    def test_agrees_with_dendropy_on_random_trees(self):
        rng = random.Random(2)
        for _ in range(100):
            newick = random_unrooted(rng, rng.randint(5, 25))
            tree = read(newick)
            nodes = [n for n in tree.postorder_node_iter() if n is not tree.seed_node]
            clade = sorted(l.taxon.label for l in rng.choice(nodes).leaf_iter())
            assert treekit.root_on_outgroup(tree, clade) is None
            sides = [sorted(l.taxon.label for l in k.leaf_iter())
                     for k in tree.seed_node.child_nodes()]
            assert clade in sides
            # the total length is unchanged by rooting
            assert tree.length() == pytest.approx(read(newick).length())


# ---------------------------------------------------------------- midpoint rooting

class TestMidpoint:

    def test_known_tree(self, run, write):
        # A-B and C-B are the longest paths (5); the midpoint is 2.5 from B
        path = write("m.tre", "(A:1,B:4,C:1);\n")
        result = run("root", "--midpoint", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert root_sides(result.stdout.strip()) == ([("A", "C"), ("B",)], [1.5, 2.5])

    def test_ultrametric_tree_keeps_its_root(self, run, write):
        newick = "[&R] (((A:1,B:1):1,C:2):1,(D:2.5,E:2.5):0.5);"
        path = write("t.tre", newick + "\n")
        result = run("root", "--midpoint", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert root_sides(result.stdout.strip()) == root_sides(newick)

    def test_root_is_halfway_between_the_farthest_tips(self, run, write):
        rng = random.Random(4)
        trees = "".join(random_unrooted(rng, rng.randint(4, 30)) + "\n" for _ in range(50))
        path = write("r.tre", trees)
        result = run("root", "--midpoint", "-i", str(path))
        assert result.returncode == 0, result.stderr
        for line in tree_lines(result.stdout):
            tree = read(line)
            depths = [max(l.distance_from_root() for l in k.leaf_iter())
                      for k in tree.seed_node.child_nodes()]
            assert depths[0] == pytest.approx(depths[1])

    def test_agrees_with_dendropy(self):
        rng = random.Random(9)
        for _ in range(100):
            newick = random_unrooted(rng, rng.randint(4, 30))
            mine = read(newick)
            assert treekit.root_at_midpoint(mine) is None
            reference = read(newick)
            reference.reroot_at_midpoint()
            assert root_sides(mine.as_string(schema="newick")) == \
                root_sides(reference.as_string(schema="newick"))

    def test_no_branch_lengths_cannot_be_rooted(self, run, write):
        path = write("n.tre", "((A,B),(C,D),E);\n")
        result = run("root", "--midpoint", "-i", str(path))
        assert result.returncode == 0
        assert "cannot be rooted (midpoint rooting needs branch lengths): tree 1" in result.stderr

    def test_zero_lengths_cannot_be_rooted(self, run, write):
        path = write("z.tre", "((A:0,B:0):0,C:0,D:0);\n")
        result = run("root", "--midpoint", "--strict", "-i", str(path))
        assert result.returncode == 1
        assert "all branch lengths are zero" in result.stderr

    def test_log(self, run, write, tmp_path):
        path = write("m.tre", "(A:1,B:4,C:1);\n")
        log = tmp_path / "run.log"
        run("root", "--midpoint", "-l", str(log), "-i", str(path))
        assert "rooted 1 of 1 trees (midpoint)" in log.read_text()


# ---------------------------------------------------------------- unrooting

class TestUnroot:

    def test_unroot_merges_the_root_branches(self, run, write):
        path = write("r.tre", "[&R] ((A:1,B:1):1,(C:1,D:1):2);\n")
        result = run("root", "--unroot", "-i", str(path))
        assert result.returncode == 0, result.stderr
        out = result.stdout.strip()
        assert out.startswith("[&U]")
        tree = read(out)
        assert len(tree.seed_node.child_nodes()) == 3
        assert tree.length() == pytest.approx(7.0)

    def test_unrooted_input_stays_unrooted(self, run, write):
        path = write("u.tre", UNROOTED)
        result = run("root", "--unroot", "-i", str(path))
        assert len(read(result.stdout.strip()).seed_node.child_nodes()) == 3

    def test_summary_sees_the_change(self, run, write, tmp_path):
        path = write("r.tre", "[&R] ((A:1,B:1):1,(C:1,D:1):2);\n")
        out = tmp_path / "u.tre"
        run("root", "--unroot", "-o", str(out), "-i", str(path))
        assert run("summary", "-i", str(out)).stdout.splitlines()[1].split("\t")[3] == "0"

    def test_two_tip_tree(self, run, write):
        path = write("two.tre", "(A:1,B:1);\n")
        result = run("root", "--unroot", "-i", str(path))
        assert result.returncode == 0, result.stderr


# ---------------------------------------------------------------- branch supports follow their split

class TestNodeLabels:

    LABELLED = "((A:1,B:1)90:1,(C:1,D:1)70:1,E:1);\n"

    @pytest.mark.parametrize("args", [
        ("-n", "A"), ("-n", "E"), ("-n", "C", "D"), ("--midpoint",), ("--unroot",),
    ])
    def test_labels_stay_on_their_split(self, run, write, args):
        path = write("l.tre", self.LABELLED)
        result = run("root", *args, "-i", str(path))
        assert result.returncode == 0, result.stderr
        labels = split_labels(result.stdout.strip())
        assert labels.get(("A", "B")) == "90"
        # the split {C, D} | {A, B, E} is written {A, B, E} or {C, D}
        assert labels.get(("C", "D"), labels.get(("A", "B", "E"))) == "70"
        # no label ends up on a split that had none
        assert set(labels.values()) <= {"90", "70"}

    def test_dendropy_alone_would_move_the_label(self):
        # documents why treekit re-assigns labels: plain dendropy puts 90 on {A}
        tree = read(self.LABELLED)
        a = [l for l in tree.leaf_node_iter() if l.taxon.label == "A"][0]
        tree.reroot_at_edge(a.edge, length1=0.5, length2=0.5)
        assert split_labels(tree.as_string(schema="newick")).get(("A", "B")) != "90"


# ---------------------------------------------------------------- options

class TestRootOptions:

    @pytest.mark.parametrize("args", [
        (), ("--midpoint", "--unroot"), ("-n", "A", "--midpoint"), ("-n", "A", "-f", "x.txt"),
    ])
    def test_exactly_one_mode(self, run, write, args):
        path = write("u.tre", UNROOTED)
        result = run("root", *args, "-i", str(path))
        assert result.returncode == 2

    def test_nexus_output(self, run, write):
        path = write("u.tre", UNROOTED)
        result = run("root", "-n", "E", "-u", "nexus", "-i", str(path))
        assert result.stdout.startswith("#NEXUS")
        assert "[&R]" in result.stdout

    def test_several_files_keep_their_order(self, run, write):
        write("a.tre", UNROOTED)
        write("b.tre", "((A:1,E:1):1,(C:1,D:1):1,B:1);\n")
        result = run("root", "-n", "E", "-i", "a.tre", "b.tre")
        first, second = tree_lines(result.stdout)
        assert root_sides(first)[0] == [("A", "B", "C", "D"), ("E",)]
        assert root_sides(second)[0] == [("A", "B", "C", "D"), ("E",)]

    def test_beast_example_midpoint_keeps_every_root(self, run):
        # the BEAST time trees in examples/ are ultrametric: midpoint = original root
        import pathlib
        example = pathlib.Path(treekit.__file__).parent / "examples" / "concat_6genes.trees"
        before = run("summary", "-i", str(example))
        result = run("root", "--midpoint", "-i", str(example))
        assert result.returncode == 0, result.stderr
        original = treekit.new_tree_list()
        original.read(path=str(example), schema="nexus", case_sensitive_taxon_labels=True,
                      preserve_underscores=True)
        for line, tree in zip(tree_lines(result.stdout), original):
            assert root_sides(line)[0] == root_sides(tree.as_string(schema="newick"))[0]
