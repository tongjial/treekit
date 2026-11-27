"""Regression tests for the five issues left open after the review of pull request #2."""
import pathlib

import pytest

import treekit

EXAMPLES = pathlib.Path(treekit.__file__).parent / "examples"
MCC = EXAMPLES / "concat_6genes.mcc.tre"


def read(source, schema="newick"):
    """A tree list with comments kept as text, from a string or a path."""
    tl = treekit.new_tree_list()
    where = {"path": str(source)} if isinstance(source, pathlib.Path) else {"data": source}
    tl.read(schema=schema, case_sensitive_taxon_labels=True, preserve_underscores=True,
            extract_comment_metadata=False, **where)
    return tl


def tips_below(node):
    return frozenset(leaf.taxon.label for leaf in node.leaf_iter())


def comments_of(tree_list):
    return sorted(c for tree in tree_list for node in tree.postorder_node_iter()
                  for c in node.comments)


def items_of(tree_list):
    return {item for comment in comments_of(tree_list)
            for item in (treekit.split_comment(comment) or [])}


def node_scoped(tree):
    """{(tips below the node, key): item text} for the internal nodes, root excluded."""
    found = {}
    for node in tree.postorder_internal_node_iter():
        if node is tree.seed_node:
            continue
        for comment in node.comments:
            for item in treekit.split_comment(comment) or []:
                if treekit.is_node_scoped(treekit.comment_key(item)):
                    found[(tips_below(node), treekit.comment_key(item))] = item
    return found


def branch_data(tree, clade):
    """(label, items) of the node whose tips below it are `clade`."""
    for node in tree.postorder_internal_node_iter():
        if node is not tree.seed_node and tips_below(node) == frozenset(clade):
            items = [i for c in node.comments for i in (treekit.split_comment(c) or [])]
            return node.label, items
    raise AssertionError("no node with the clade %s" % sorted(clade))


# ---------------------------------------------------------------- 1. an annotated single tip

class TestAnnotatedSingleTip:

    @pytest.mark.parametrize("newick", [
        "A[&state=blue];", "A:1[&state=blue];", "A[&state=blue]:1;",
        "A[&a=1][&b=2]:3 [&c=4];", "'A[1] b'[&x=1]:3;",
    ])
    def test_a_file_with_one_annotated_tip_is_a_tree(self, run, write, newick):
        path = write("one.tre", newick + "\n")
        result = run("taxa", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert len(result.stdout.splitlines()) == 1

    @pytest.mark.parametrize("text", ["A[&x=1]B;", "A[&x=1", "A:1:2;", "hello world",
                                      "'A b'[&x=\"1]2\"]:3;"])
    def test_other_text_is_still_not_a_tree(self, run, write, text):
        path = write("bad.tre", text + "\n")
        result = run("taxa", "-i", str(path))
        assert result.returncode == 1
        assert "not a newick or nexus tree file" in result.stderr

    @pytest.mark.parametrize("command", [("extract", "-n", "A"), ("keep", "-n", "A", "B", "C"),
                                         ("remove", "-n", "B", "C")])
    def test_what_treekit_writes_it_can_read(self, run, write, tmp_path, command):
        path = write("t.tre", "(A[&state=blue]:1,B:1,C:1);\n")
        out = tmp_path / "out.tre"
        result = run(*command, "-o", str(out), "-i", str(path))
        assert result.returncode == 0, result.stderr
        again = run("taxa", "-i", str(out))
        assert again.returncode == 0, again.stderr
        assert "&state=blue" in out.read_text()


# ---------------------------------------------------------------- 2. labels and annotations merge

class TestMergedBranchAtTheOldRoot:

    def reroot(self, run, write, newick):
        path = write("m.tre", newick + "\n")
        result = run("root", "-n", "A", "-i", str(path))
        assert result.returncode == 0, result.stderr
        return result

    def test_label_on_one_branch_and_annotation_on_the_other(self, run, write):
        result = self.reroot(run, write, "[&R] ((A:1,B:1)[&posterior=0.9]:1,(C:1,D:1)70:2);")
        assert result.stderr == ""
        label, items = branch_data(read(result.stdout.strip())[0], "CD")
        assert label == "70" and items == ["posterior=0.9"]

    def test_the_same_in_the_other_order(self, run, write):
        result = self.reroot(run, write, "[&R] ((A:1,B:1)70:1,(C:1,D:1)[&posterior=0.9]:2);")
        assert result.stderr == ""
        label, items = branch_data(read(result.stdout.strip())[0], "CD")
        assert label == "70" and items == ["posterior=0.9"]

    def test_annotations_with_different_keys_are_both_kept(self, run, write):
        result = self.reroot(
            run, write, "[&R] ((A:1,B:1)[&posterior=0.9]:1,(C:1,D:1)[&rate=1.1]:2);")
        assert result.stderr == ""
        _, items = branch_data(read(result.stdout.strip())[0], "CD")
        assert sorted(items) == ["posterior=0.9", "rate=1.1"]

    def test_the_same_value_on_both_branches_is_no_conflict(self, run, write):
        result = self.reroot(
            run, write, "[&R] ((A:1,B:1)[&posterior=0.9]:1,(C:1,D:1)[&posterior=0.9]:2);")
        assert result.stderr == ""
        assert branch_data(read(result.stdout.strip())[0], "CD")[1] == ["posterior=0.9"]

    def test_the_same_key_with_other_values_is_reported_and_the_first_kept(self, run, write):
        result = self.reroot(
            run, write, "[&R] ((A:1,B:1)[&posterior=0.9]:1,(C:1,D:1)[&posterior=0.7]:2);")
        assert "different labels or annotations on the two branches at the old root" in result.stderr
        assert branch_data(read(result.stdout.strip())[0], "CD")[1] == ["posterior=0.9"]

    def test_two_different_labels_are_reported(self, run, write):
        result = self.reroot(run, write, "[&R] ((A:1,B:1)90:1,(C:1,D:1)70:2);")
        assert "different labels or annotations" in result.stderr
        assert branch_data(read(result.stdout.strip())[0], "CD")[0] == "90"

    def test_unrooting_merges_the_same_way(self, run, write):
        path = write("m.tre", "[&R] ((A:1,B:1)[&posterior=0.9]:1,(C:1,D:1)70:2);\n")
        result = run("root", "--unroot", "-i", str(path))
        assert result.returncode == 0 and result.stderr == ""
        tree = read(result.stdout.strip())[0]
        found = {node.label: [i for c in node.comments for i in (treekit.split_comment(c) or [])]
                 for node in tree.postorder_internal_node_iter() if node is not tree.seed_node}
        assert found == {"70": ["posterior=0.9"]}

    def test_collapse_finds_the_numeric_label_after_rooting(self, run, write, tmp_path):
        # the case of the review: 70 is the only numeric label of the input
        path = write("m.tre", "[&R] ((A:1,B:1)[&posterior=0.9]:1,(C:1,D:1)70:2);\n")
        rooted = tmp_path / "rooted.tre"
        run("root", "-n", "A", "-o", str(rooted), "-i", str(path))
        result = run("collapse", "-m", "80", "-i", str(rooted))
        assert result.returncode == 0, result.stderr
        assert "no numeric support" not in result.stderr


# ---------------------------------------------------------------- 3. an existing root stays

class TestRootAlreadyOnTheSplit:

    CHAIN = "[&R] (((A:1,B:1):1):2,(C:1,D:1):3);"

    @pytest.mark.parametrize("outgroup", [("A", "B"), ("C", "D")])
    def test_a_single_child_chain_does_not_move_the_root(self, run, write, outgroup):
        path = write("c.tre", self.CHAIN + "\n")
        result = run("root", "-n", *outgroup, "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stderr == ""
        tree = read(result.stdout.strip())[0]
        distances = {l.taxon.label: l.distance_from_root() for l in tree.leaf_node_iter()}
        assert distances == {"A": 4.0, "B": 4.0, "C": 4.0, "D": 4.0}
        assert result.stdout.strip() == read(self.CHAIN)[0].as_string(schema="newick").strip()

    def test_an_unrooted_tree_with_a_chain_is_rooted_on_the_split_and_keeps_its_distances(
            self, run, write):
        # two branches of the chain have the same split; whichever is used, the distances
        # between the tips do not change (DendroPy merges the chain while rerooting)
        newick = "(((A:1,B:1):1):2,C:1,D:1);"
        path = write("u.tre", newick + "\n")
        result = run("root", "-n", "A", "B", "-i", str(path))
        assert result.returncode == 0, result.stderr
        tree = read(result.stdout.strip())[0]
        sides = sorted(sorted(tips_below(c)) for c in tree.seed_node.child_nodes())
        assert sides == [["A", "B"], ["C", "D"]]

        def distances(t):
            matrix = t.phylogenetic_distance_matrix()
            taxa = {x.label: x for x in t.taxon_namespace}
            return {(a, b): round(matrix.distance(taxa[a], taxa[b]), 9)
                    for a in "ABCD" for b in "ABCD" if a < b}
        assert distances(tree) == distances(read(newick)[0])

    def test_a_normal_rooted_tree_is_still_rerooted_when_the_split_is_elsewhere(self, run, write):
        path = write("r.tre", "[&R] (((A:1,B:1):1,C:2):1,(D:1,E:1):2);\n")
        result = run("root", "-n", "A", "B", "-i", str(path))
        tree = read(result.stdout.strip())[0]
        sides = sorted(sorted(tips_below(c)) for c in tree.seed_node.child_nodes())
        assert sides == [["A", "B"], ["C", "D", "E"]]


# ---------------------------------------------------------------- 4. heights are not copied to other nodes

class TestNodeHeightsWhenTheRootMoves:

    def test_a_height_is_not_moved_to_another_node(self, run, write):
        path = write("h.tre", "((A:1,B:1)[&height=1]:1,C:2,D:2);\n")
        result = run("root", "-n", "A", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "height" not in result.stdout
        assert "node annotations of heights or ages" in result.stderr
        assert "tree 1" in result.stderr

    def test_clades_that_do_not_change_keep_their_heights(self, run, write):
        path = write("h.tre", "(((A:1,B:1)[&height=1,posterior=0.9]:1,C:2)[&height=3]:1,D:4,E:4);\n")
        result = run("root", "-n", "E", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stderr == ""
        tree = read(result.stdout.strip())[0]
        assert branch_data(tree, "AB")[1] == ["posterior=0.9", "height=1"]
        assert branch_data(tree, "ABC")[1] == ["height=3"]

    def test_other_annotations_still_follow_the_split(self, run, write):
        path = write("h.tre", "((A:1,B:1)[&height=1,posterior=0.9]:1,C:2,D:2);\n")
        result = run("root", "-n", "A", "-i", str(path))
        # the split {A, B} | {C, D} is now the clade {C, D}: it has the posterior, not the height
        _, items = branch_data(read(result.stdout.strip())[0], "CD")
        assert items == ["posterior=0.9"]

    @pytest.mark.parametrize("item", ["height=1", "height_95%_HPD={0.5,1.5}", "age_mean=2", "age=2"])
    def test_every_height_or_age_item_is_node_scoped(self, run, write, item):
        path = write("h.tre", "((A:1,B:1)[&%s]:1,C:2,D:2);\n" % item)
        result = run("root", "-n", "A", "-i", str(path))
        assert item.split("=")[0] not in result.stdout

    @pytest.mark.parametrize("item", ["posterior=0.9", "rate=1.1", "length=1", "heights=1",
                                      "length_95%_HPD={0.5,1.5}"])
    def test_other_items_are_branch_scoped(self, run, write, item):
        path = write("h.tre", "((A:1,B:1)[&%s]:1,C:2,D:2);\n" % item)
        result = run("root", "-n", "A", "-i", str(path))
        assert result.stderr == ""
        assert item in result.stdout

    def test_unrooting_removes_the_heights_of_the_root(self, run, write):
        path = write("r.tre", "[&R] ((A:1,B:1)[&height=1]:1,(C:1,D:1)[&height=1]:2)[&height=3];\n")
        result = run("root", "--unroot", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "height=3" not in result.stdout
        assert "node annotations of heights or ages" in result.stderr

    def test_in_a_real_time_tree_no_height_ends_on_another_clade(self, run, tmp_path):
        before = node_scoped(read(MCC, "nexus")[0])
        assert len(before) > 100
        source = read(MCC, "nexus")[0]
        clade = next(sorted(tips_below(n)) for n in source.preorder_node_iter()
                     if n.is_internal() and n.parent_node is not None
                     and n.parent_node is not source.seed_node and 3 <= len(tips_below(n)) <= 20)
        for outgroup in (["Coris_monspeliensis_Y925"], clade):
            out = tmp_path / "rooted.tre"
            result = run("root", "-n", *outgroup, "-o", str(out), "-i", str(MCC))
            assert result.returncode == 0, result.stderr
            after = node_scoped(read(out)[0])
            # every remaining item is on the clade it belonged to, with its own text
            assert after and set(after.items()) <= set(before.items())
            assert len(after) < len(before)
            assert "node annotations of heights or ages" in result.stderr

    def test_the_midpoint_of_a_time_tree_changes_nothing(self, run, tmp_path):
        out = tmp_path / "mid.tre"
        result = run("root", "--midpoint", "-o", str(out), "-i", str(MCC))
        assert result.returncode == 0, result.stderr
        assert "node annotations" not in result.stderr
        assert comments_of(read(out)) == comments_of(read(MCC, "nexus"))


# ---------------------------------------------------------------- 5. comments are kept as written

class TestCommentsAreKeptAsWritten:

    QUOTED = '(A[&note="red, blue"]:1,B:1,C:1);\n'

    @pytest.mark.parametrize("command", [
        ("keep", "-n", "A"), ("remove", "-n", "B"), ("extract", "-n", "A", "B"),
        ("prune", "--max-length", "100"), ("root", "-n", "B"), ("root", "--midpoint"),
        ("collapse", "-m", "0"),
    ])
    def test_a_quoted_value_with_a_comma_is_not_cut(self, run, write, command):
        path = write("q.tre", self.QUOTED)
        result = run(*command, "-i", str(path))
        if command[0] == "collapse":
            assert result.returncode == 0       # no supports: written unchanged, with a warning
        assert result.returncode == 0, result.stderr
        assert '[&note="red, blue"]' in result.stdout

    def test_the_same_in_nexus(self, run, write):
        path = write("q.tre", self.QUOTED)
        result = run("keep", "-n", "A", "-u", "nexus", "-i", str(path))
        assert '[&note="red, blue"]' in result.stdout

    def test_keep_gives_every_comment_of_a_beast_tree_back(self, run, tmp_path):
        out = tmp_path / "kept.tre"
        result = run("keep", "-n", "Coris_monspeliensis_Y925", "-o", str(out), "-i", str(MCC))
        assert result.returncode == 0, result.stderr
        before, after = comments_of(read(MCC, "nexus")), comments_of(read(out))
        assert len(before) > 100
        assert after == before                       # the same comments, character for character

    @pytest.mark.parametrize("command", [
        ("remove", "-n", "Coris_monspeliensis_Y925"),
        ("extract", "-n", "Coris_monspeliensis_Y925", "Badula_crassa_Y24864",
         "Ardisia_solanacea_Y17988", "Embelia_vestita_YG205"),
        ("prune", "--max-length", "1000"),
        ("collapse", "-m", "0.95", "--support-key", "posterior"),
    ])
    def test_no_command_writes_a_comment_that_was_not_in_the_input(self, run, tmp_path, command):
        out = tmp_path / "out.tre"
        result = run(*command, "-o", str(out), "-i", str(MCC))
        assert result.returncode == 0, result.stderr
        written = set(comments_of(read(out)))
        assert written and written <= set(comments_of(read(MCC, "nexus")))

    def test_rooting_writes_only_items_that_were_in_the_input(self, run, tmp_path):
        out = tmp_path / "rooted.tre"
        result = run("root", "-n", "Coris_monspeliensis_Y925", "-o", str(out), "-i", str(MCC))
        assert result.returncode == 0, result.stderr
        assert items_of(read(out)) <= items_of(read(MCC, "nexus"))

    def test_collapse_reads_a_support_that_follows_a_quoted_comma(self, run, write):
        path = write("c.tre", '((A:1,B:1)[&note="red, blue",posterior=0.4]:1,'
                              '(C:1,D:1)[&posterior=0.99]:1,E:1);\n')
        result = run("collapse", "-m", "0.95", "--support-key", "posterior", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "posterior=0.4" not in result.stdout      # collapsed with its node
        assert "posterior=0.99" in result.stdout

    def test_collapse_reads_a_quoted_support(self, run, write):
        path = write("c.tre", '((A:1,B:1)[&posterior="0.4"]:1,(C:1,D:1)[&posterior="0.99"]:1,E:1);\n')
        result = run("collapse", "-m", "0.95", "--support-key", "posterior", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert '"0.99"' in result.stdout and '"0.4"' not in result.stdout


class TestCommentHelpers:

    def test_split_comment(self):
        assert treekit.split_comment('&a=1,b="x, y",c={1,2},d') == ['a=1', 'b="x, y"', 'c={1,2}', 'd']
        assert treekit.split_comment("&height_95%_HPD={0.0,2.8E-14},height=1.3") == \
            ["height_95%_HPD={0.0,2.8E-14}", "height=1.3"]
        assert treekit.split_comment("&") == []
        assert treekit.split_comment("free text") is None

    def test_comment_value(self):
        comments = ['&note="red, blue",posterior=0.4', "free", '&q="0.5"']
        assert treekit.comment_value(comments, "posterior") == "0.4"
        assert treekit.comment_value(comments, "q") == "0.5"
        assert treekit.comment_value(comments, "note") == "red, blue"
        assert treekit.comment_value(comments, "missing") is None

    @pytest.mark.parametrize("key, expected", [
        ("height", True), ("height_95%_HPD", True), ("height_median", True), ("age", True),
        ("age_mean", True), ("age_95%HPD", True),
        ("posterior", False), ("rate", False), ("length", False), ("heights", False),
        ("ageing", False), ("prob", False),
    ])
    def test_is_node_scoped(self, key, expected):
        assert treekit.is_node_scoped(key) is expected


# ---------------------------------------------------------------- 4 (again). nodes with the same tips

class TestHeightsOnNodesWithTheSameTips:

    CHAIN = "(((A:1,B:1)[&height=1]:1)[&height=2]:1,C:3,D:3,E:3);\n"

    def test_a_surviving_node_keeps_its_own_height(self, run, write):
        # a node with one child and its child have the same tips; DendroPy removes the
        # upper one while rerooting, and the node that stays must keep height=1
        path = write("c.tre", self.CHAIN)
        result = run("root", "-n", "E", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "height=1" in result.stdout and "height=2" not in result.stdout
        tree = read(result.stdout.strip())[0]
        assert branch_data(tree, "AB")[1] == ["height=1"]
        # the height of the node that is gone is reported as removed
        assert "node annotations of heights or ages" in result.stderr

    def test_unrooting_keeps_each_height_on_its_own_node(self, run, write):
        # the chain stays after unrooting: both nodes keep their own height. The node of
        # (C, D) is merged into the root, so its height is removed.
        path = write("c.tre", "[&R] (((A:1,B:1)[&height=1]:1)[&height=2]:1,(C:3,D:3)[&height=3]:1);\n")
        result = run("root", "--unroot", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "(A:1.0,B:1.0):1.0[&height=1]" in result.stdout     # the lower node of the chain
        assert "[&height=2]" in result.stdout                      # the upper node
        assert "height=3" not in result.stdout
        assert "node annotations of heights or ages" in result.stderr

    def test_two_nodes_with_the_same_tips_and_a_third_that_stays(self, run, write):
        path = write("c.tre", "(((A:1,B:1)[&height=1]:1)[&height=2]:1,(C:3,D:3)[&height=3]:1,E:5);\n")
        result = run("root", "-n", "E", "-i", str(path))
        assert result.returncode == 0, result.stderr
        tree = read(result.stdout.strip())[0]
        assert branch_data(tree, "AB")[1] == ["height=1"]
        assert branch_data(tree, "CD")[1] == ["height=3"]
        assert "height=2" not in result.stdout

    def test_on_a_time_tree_with_a_unary_node_above_every_second_node(self, run, tmp_path):
        # every unary node gets a height of its own that no other node has (negative):
        # rooting removes the unary nodes, and no height may end up on another node
        source = read(MCC, "nexus")[0]
        original = node_scoped(read(MCC, "nexus")[0])
        inner = [n for n in source.postorder_internal_node_iter()
                 if n is not source.seed_node and n.edge.length]
        for index, node in enumerate(inner[::2]):
            parent = node.parent_node
            unary = type(node)()
            unary.edge.length = node.edge.length / 2
            node.edge.length = node.edge.length / 2
            unary.comments.append("&height=%s" % (-1.0 - index))
            position = parent.child_nodes().index(node)
            parent.remove_child(node)
            parent.insert_child(position, unary)
            unary.add_child(node)
        path = tmp_path / "unary.tre"
        path.write_text(source.as_string(schema="newick", suppress_item_comments=False,
                                         unquoted_underscores=True) + "\n")
        assert path.read_text().count("height=-") == len(inner[::2]) > 30   # the unary nodes are there

        out = tmp_path / "rooted.tre"
        result = run("root", "-n", "Coris_monspeliensis_Y925", "-o", str(out), "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "height=-" not in out.read_text()
        after = node_scoped(read(out)[0])
        assert len(after) > 100
        # every height left is a height of the original tree, on the clade it belonged to
        assert set(after.items()) <= set(original.items())


# ---------------------------------------------------------------- 6. NHX comments

class TestNhxComments:

    NHX = "((A:1,B:1)[&&NHX:B=40]:1,(C:1,D:1)[&&NHX:B=99]:1,E:1);\n"
    NHX2 = ("((A:1,B:1)[&&NHX:S=human:B=40:E=1.1]:1,(C:1,D:1)[&&NHX:S=mouse:B=99]:1,E:1);\n")

    def test_collapse_reads_the_support_of_an_nhx_comment(self, run, write):
        path = write("n.tre", self.NHX)
        result = run("collapse", "-m", "95", "--support-key", "B", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stderr == ""
        assert result.stdout.strip() == "(A:2.0,B:2.0,(C:1.0,D:1.0):1.0[&&NHX:B=99],E:1.0);"

    def test_collapse_with_several_fields(self, run, write):
        path = write("n.tre", self.NHX2)
        result = run("collapse", "-m", "95", "--support-key", "B", "-i", str(path))
        assert result.stderr == ""
        assert "&&NHX:S=mouse:B=99" in result.stdout and "human" not in result.stdout

    def test_other_keys_of_an_nhx_comment_can_be_the_support(self, run, write):
        path = write("n.tre", self.NHX2)
        result = run("collapse", "-m", "1.5", "--support-key", "E", "-i", str(path))
        # E=1.1 is the only value of its key: collapsed; the other branch has no E
        assert "human" not in result.stdout

    def test_nhx_comments_are_written_as_read(self, run, write):
        path = write("n.tre", self.NHX2)
        result = run("keep", "-n", "A", "-i", str(path))
        assert "[&&NHX:S=human:B=40:E=1.1]" in result.stdout
        assert "[&&NHX:S=mouse:B=99]" in result.stdout

    @pytest.mark.parametrize("outgroup", ["A", "E", "C"])
    def test_rerooting_keeps_the_nhx_format_and_all_fields(self, run, write, outgroup):
        path = write("n.tre", self.NHX2)
        result = run("root", "-n", outgroup, "-i", str(path))
        assert result.returncode == 0, result.stderr
        tree = read(result.stdout.strip())[0]
        comments = sorted(c for n in tree.postorder_internal_node_iter() for c in n.comments)
        assert comments == ["&&NHX:S=human:B=40:E=1.1", "&&NHX:S=mouse:B=99"]

    def test_both_styles_on_one_node_stay_two_comments(self, run, write):
        path = write("n.tre", "((A:1,B:1)[&&NHX:B=40][&posterior=0.9]:1,(C:1,D:1):1,E:1);\n")
        result = run("root", "-n", "E", "-i", str(path))
        assert result.returncode == 0, result.stderr
        tree = read(result.stdout.strip())[0]
        comments = sorted(c for n in tree.postorder_internal_node_iter() for c in n.comments)
        assert comments == ["&&NHX:B=40", "&posterior=0.9"]

    def test_helpers(self):
        assert treekit.split_comment("&&NHX:S=human:B=40") == ["S=human", "B=40"]
        assert treekit.comment_style("&&NHX:B=1") == "nhx" and treekit.comment_style("&B=1") == "beast"
        assert treekit.comment_value(["&&NHX:S=human:B=40"], "B") == "40"
        # items stay in their order; a comment ends where the style changes
        assert treekit.join_items([("nhx", "B=40"), ("nhx", "S=x"), ("beast", "p=1"), ("beast", "q=2")]) == \
            ["&&NHX:B=40:S=x", "&p=1,q=2"]
        assert treekit.join_items([("nhx", "B=40"), ("beast", "p=1"), ("nhx", "S=x")]) == \
            ["&&NHX:B=40", "&p=1", "&&NHX:S=x"]
        assert treekit.join_items([]) == []


# ---------------------------------------------------------------- 7. one key, two comment styles

def b_values(text):
    """The values of every item with the key B, in the order they are written."""
    tree = read(text.strip())[0]
    values = []
    for node in tree.preorder_node_iter():
        for comment in node.comments:
            for item in treekit.split_comment(comment) or []:
                if treekit.comment_key(item) == "B":
                    values.append(item.partition("=")[2])
    return values


class TestSameKeyInBothStyles:

    MIX = "((A:1,B:1)[&S=mouse][&&NHX:B=99][&B=40]:1,(C:1,D:1):1,E:1);\n"

    @pytest.mark.parametrize("outgroup", ["C", "E", "A"])
    def test_rerooting_keeps_the_order_of_the_items(self, run, write, outgroup):
        path = write("m.tre", self.MIX)
        result = run("root", "-n", outgroup, "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert b_values(result.stdout) == ["99", "40"]
        assert "[&S=mouse][&&NHX:B=99][&B=40]" in result.stdout      # as written

    def test_collapse_gives_the_same_answer_before_and_after_rerooting(self, run, write, tmp_path):
        path = write("m.tre", self.MIX)
        rooted = tmp_path / "rooted.tre"
        run("root", "-n", "C", "-o", str(rooted), "-i", str(path))
        before = run("collapse", "-m", "95", "--support-key", "B", "-i", str(path))
        after = run("collapse", "-m", "95", "--support-key", "B", "-i", str(rooted))
        # the first value of the key (99) is the support: the branch is kept in both
        assert before.returncode == after.returncode == 0, (before.stderr, after.stderr)
        assert before.stderr == "" == after.stderr
        assert b_values(before.stdout) == b_values(after.stdout) == ["99", "40"]

    def test_collapse_does_remove_the_branch_when_the_first_value_is_low(self, run, write, tmp_path):
        # the same tree with the order of the two B items swapped: now 40 is the support
        path = write("m.tre", self.MIX.replace("[&&NHX:B=99][&B=40]", "[&B=40][&&NHX:B=99]"))
        rooted = tmp_path / "rooted.tre"
        run("root", "-n", "C", "-o", str(rooted), "-i", str(path))
        for source in (path, rooted):
            result = run("collapse", "-m", "95", "--support-key", "B", "-i", str(source))
            assert b_values(result.stdout) == [], source             # collapsed, both times

    CLASH = "[&R] (((A:1,B:1):1,C:2)[&&NHX:B=40]:1,(D:1,E:1)[&B=99]:1);\n"

    def test_different_values_on_the_two_old_root_branches_are_reported(self, run, write):
        path = write("c.tre", self.CLASH)
        result = run("root", "--unroot", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "different labels or annotations on the two branches at the old root" in result.stderr
        assert b_values(result.stdout) == ["40"]                   # the first value, in its own style
        assert "[&&NHX:B=40]" in result.stdout

    def test_the_same_value_in_both_styles_is_no_conflict(self, run, write):
        path = write("c.tre", self.CLASH.replace("B=99", "B=40"))
        result = run("root", "--unroot", "-i", str(path))
        assert result.stderr == ""
        assert b_values(result.stdout) == ["40"]

    def test_rooting_on_an_outgroup_reports_it_too(self, run, write):
        path = write("c.tre", self.CLASH)
        result = run("root", "-n", "A", "-i", str(path))
        assert "different labels or annotations" in result.stderr
        assert b_values(result.stdout) == ["40"]

    def test_the_value_that_collapse_reads_is_the_one_that_was_kept(self, run, write, tmp_path):
        path = write("c.tre", self.CLASH)
        unrooted = tmp_path / "u.tre"
        run("root", "--unroot", "-o", str(unrooted), "-i", str(path))
        result = run("collapse", "-m", "50", "--support-key", "B", "-i", str(unrooted))
        # B=40 is below 50: the merged branch is collapsed, and no 99 is left to read
        assert "B=99" not in result.stdout and "B=40" not in result.stdout

    def test_different_keys_in_the_two_styles_are_both_kept_in_order(self, run, write):
        path = write("t.tre", "((A:1,B:1)[&posterior=0.9][&&NHX:S=human:B=40]:1,(C:1,D:1):1,E:1);\n")
        result = run("root", "-n", "E", "-i", str(path))
        assert result.stderr == ""
        assert "[&posterior=0.9][&&NHX:S=human:B=40]" in result.stdout
