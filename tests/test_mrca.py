import itertools
import random

import pytest

import treekit

# rooted: ((A,(B,C)),(D,E))
ROOTED = "[&R] ((A:1,(B:1,C:1):1):1,(D:1,E:1):1);\n"
UNROOTED = "((A:1,B:1):1,(C:1,D:1):1,E:1);\n"
HEADER = ["source", "index", "group", "rooted", "members_present", "members_total",
          "monophyletic", "num_intruders", "intruders"]


def read(newick):
    tl = treekit.new_tree_list()
    tl.read(data=newick, schema="newick", case_sensitive_taxon_labels=True,
            preserve_underscores=True)
    return tl[0]


def rows(stdout):
    lines = stdout.splitlines()
    assert lines[0].split("\t") == HEADER
    return [dict(zip(HEADER, line.split("\t"))) for line in lines[1:]]


def tree_lines(stdout):
    return [line for line in stdout.splitlines() if line.strip()]


def leaves(newick):
    return sorted(l.taxon.label for l in read(newick).leaf_node_iter())


def random_tree(rng, n, rooted):
    labels = ["t%d" % i for i in range(n)]
    nodes = ["%s:1" % l for l in labels]
    stop = 2 if rooted else 3
    while len(nodes) > stop:
        a = nodes.pop(rng.randrange(len(nodes)))
        b = nodes.pop(rng.randrange(len(nodes)))
        nodes.append("(%s,%s):1" % (a, b))
    return "(%s);" % ",".join(nodes), labels


def splits_by_brute_force(tree):
    """Every tip set separated by one branch, both sides, for the unrooted check."""
    full = frozenset(l.taxon.label for l in tree.leaf_node_iter())
    sides = set()
    for node in tree.preorder_node_iter():
        if node is tree.seed_node:
            continue
        below = frozenset(l.taxon.label for l in node.leaf_iter())
        sides.add(below)
        sides.add(full - below)
    return sides


def clades_by_brute_force(tree):
    return {frozenset(l.taxon.label for l in node.leaf_iter())
            for node in tree.preorder_node_iter()}


# ---------------------------------------------------------------- the monophyly table

class TestTable:

    def test_rooted_monophyletic(self, run, write):
        path = write("r.tre", ROOTED)
        result = run("mrca", "-n", "B", "C", "-i", str(path))
        assert result.returncode == 0, result.stderr
        [row] = rows(result.stdout)
        assert row == {"source": "r.tre", "index": "1", "group": "group", "rooted": "True",
                       "members_present": "2", "members_total": "2", "monophyletic": "True",
                       "num_intruders": "0", "intruders": ""}

    def test_rooted_not_monophyletic_lists_intruders(self, run, write):
        path = write("r.tre", ROOTED)
        [row] = rows(run("mrca", "-n", "A", "D", "-i", str(path)).stdout)
        assert row["monophyletic"] == "False"
        assert row["num_intruders"] == "3"
        assert row["intruders"] == "B,C,E"

    def test_rooted_and_unrooted_rules_differ_when_the_root_is_inside_the_group(self, run, write):
        # {A, D, E}: its MRCA is the root, so B and C intrude; unrooted, one branch
        # separates {B, C} from {A, D, E}
        path = write("r.tre", ROOTED)
        [rooted] = rows(run("mrca", "-n", "A", "D", "E", "-i", str(path)).stdout)
        [unrooted] = rows(run("mrca", "-n", "A", "D", "E", "--unrooted", "-i", str(path)).stdout)
        assert (rooted["rooted"], rooted["monophyletic"], rooted["intruders"]) == ("True", "False", "B,C")
        assert (unrooted["rooted"], unrooted["monophyletic"], unrooted["intruders"]) == ("False", "True", "")

    def test_unrooted_tree_uses_the_split_rule(self, run, write):
        path = write("u.tre", UNROOTED)
        result = run("mrca", "-n", "C", "D", "E", "-i", str(path))
        [row] = rows(result.stdout)
        assert row["rooted"] == "False"
        assert row["monophyletic"] == "True"

    def test_unrooted_not_monophyletic_has_no_intruder_list(self, run, write):
        path = write("u.tre", UNROOTED)
        [row] = rows(run("mrca", "-n", "A", "C", "-i", str(path)).stdout)
        assert row["monophyletic"] == "False"
        assert row["num_intruders"] == "NA" and row["intruders"] == "NA"

    @pytest.mark.parametrize("members", [["A"], ["A", "B", "C", "D", "E"]])
    def test_untestable_groups_are_na(self, run, write, members):
        path = write("r.tre", ROOTED)
        [row] = rows(run("mrca", "-n", *members, "-i", str(path)).stdout)
        assert row["monophyletic"] == "NA"

    def test_partial_group_is_tested_on_the_members_present(self, run, write, tmp_path):
        path = write("m.tre", ROOTED + "[&R] ((A:1,X:1):1,(D:1,E:1):1);\n")
        log = tmp_path / "run.log"
        result = run("mrca", "-n", "B", "C", "X", "-l", str(log), "-i", str(path))
        assert result.returncode == 0, result.stderr
        first, second = rows(result.stdout)
        assert (first["members_present"], first["members_total"], first["monophyletic"]) == ("2", "3", "True")
        assert (second["members_present"], second["monophyletic"]) == ("1", "NA")
        assert "group 'group': monophyletic in 1 of 1 testable trees" in log.read_text()

    def test_several_groups_from_a_file(self, run, write):
        path = write("r.tre", ROOTED)
        groups = write("groups.tsv", "# taxon\tgroup\nB\tbc\nC\tbc\n\nA\tad\nD\tad\n")
        result = run("mrca", "-g", str(groups), "-i", str(path))
        assert result.returncode == 0, result.stderr
        table = rows(result.stdout)
        assert [(r["group"], r["monophyletic"]) for r in table] == [("bc", "True"), ("ad", "False")]

    def test_groups_file_errors(self, run, write):
        path = write("r.tre", ROOTED)
        bad = write("bad.tsv", "B bc\n")
        result = run("mrca", "-g", str(bad), "-i", str(path))
        assert result.returncode == 1
        assert "line 1" in result.stderr and "taxon<TAB>group" in result.stderr
        missing = run("mrca", "-g", "no.tsv", "-i", str(path))
        assert missing.returncode == 1 and "groups file not found" in missing.stderr

    def test_source_and_index_for_several_files(self, run, write):
        write("a.tre", ROOTED)
        write("b.tre", ROOTED + ROOTED)
        table = rows(run("mrca", "-n", "B", "C", "-i", "a.tre", "b.tre").stdout)
        assert [(r["source"], r["index"]) for r in table] == [("a.tre", "1"), ("b.tre", "1"), ("b.tre", "2")]

    def test_absent_taxa(self, run, write):
        path = write("r.tre", ROOTED)
        result = run("mrca", "-n", "B", "C", "ZZZ", "-i", str(path))
        assert "1 taxon was not found in any input tree: ZZZ" in result.stderr
        strict = run("mrca", "-n", "B", "C", "ZZZ", "--strict", "-i", str(path))
        assert strict.returncode == 1
        none = run("mrca", "-n", "ZZZ", "-i", str(path))
        assert none.returncode == 1 and "none of the 1 taxon of the groups" in none.stderr

    def test_unlabelled_tips_are_rejected(self, run, write):
        path = write("u.tre", "((A:1,:1):1,(C:1,D:1):1);\n")
        result = run("mrca", "-n", "C", "D", "-i", str(path))
        assert result.returncode == 1
        assert "tree 1 has a tip without a label" in result.stderr

    def test_agrees_with_brute_force_on_random_trees(self):
        rng = random.Random(8)
        for rooted in (True, False):
            for _ in range(60):
                newick, labels = random_tree(rng, rng.randint(4, 9), rooted)
                tree = read(newick)
                truth = clades_by_brute_force(tree) if rooted else splits_by_brute_force(tree)
                for size in range(2, len(labels)):
                    for group in itertools.combinations(labels, size):
                        status = treekit.group_monophyly(tree, group)
                        assert status["rooted"] is rooted
                        assert status["monophyletic"] == (frozenset(group) in truth), (newick, group)
                        if rooted:
                            tree.is_rooted = True
                            mrca = tree.mrca(taxon_labels=list(group))
                            expected = sorted(l.taxon.label for l in mrca.leaf_iter()
                                              if l.taxon.label not in group)
                            assert status["intruders"] == expected


# ---------------------------------------------------------------- extracting clades

class TestExtractClades:

    def test_rooted_clade_includes_intruders(self, run, write):
        # the clade of {A, B} is the MRCA's clade: A, B and C
        path = write("r.tre", ROOTED)
        result = run("mrca", "-n", "A", "B", "--extract", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stdout.strip() == "[&R] (A:1.0,(B:1.0,C:1.0):1.0);"

    def test_unrooted_monophyletic_group(self, run, write):
        path = write("u.tre", UNROOTED)
        result = run("mrca", "-n", "C", "D", "--extract", "-i", str(path))
        assert result.stdout.strip() == "[&R] (C:1.0,D:1.0);"

    def test_unrooted_group_above_a_node_is_rooted_on_its_branch(self, run, write):
        # {W, Q, R} is everything except V: its clade hangs from V's branch
        path = write("u.tre", "((W:1,V:1):1,Q:1,R:1);\n")
        result = run("mrca", "-n", "W", "Q", "R", "--extract", "-i", str(path))
        assert result.stdout.strip() == "[&R] (W:1.0,(Q:1.0,R:1.0):1.0);"

    def test_unrooted_non_monophyletic_is_skipped(self, run, write):
        path = write("m.tre", UNROOTED + UNROOTED.replace("(C:1,D:1)", "(C:1,X:1)").replace("E:1", "D:1"))
        result = run("mrca", "-n", "A", "B", "--extract", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert len(tree_lines(result.stdout)) == 2
        result = run("mrca", "-n", "A", "C", "--extract", "-i", str(path))
        assert result.returncode == 1
        assert "cannot be extracted (the group is not monophyletic in an unrooted tree)" in result.stderr
        assert "no clade could be extracted" in result.stderr

    def test_skipped_clades_warn_and_strict_fails(self, run, write):
        path = write("m.tre", ROOTED + "[&R] ((A:1,X:1):1,(D:1,E:1):1);\n")
        result = run("mrca", "-n", "B", "C", "--extract", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "1 clade(s) cannot be extracted" in result.stderr
        assert "tree/group 2/group" in result.stderr
        assert len(tree_lines(result.stdout)) == 1
        strict = run("mrca", "-n", "B", "C", "--extract", "--strict", "-i", str(path))
        assert strict.returncode == 1

    def test_nexus_names_every_clade(self, run, write):
        path = write("r.tre", ROOTED)
        groups = write("groups.tsv", "B\tbc\nC\tbc\nD\tde\nE\tde\n")
        result = run("mrca", "-g", str(groups), "--extract", "-u", "nexus", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "TREE bc_tree1 =" in result.stdout
        assert "TREE de_tree1 =" in result.stdout

    def test_extracted_clades_have_no_stem_branch(self, run, write):
        path = write("r.tre", ROOTED)
        out = run("mrca", "-n", "D", "E", "--extract", "-i", str(path)).stdout.strip()
        assert out == "[&R] (D:1.0,E:1.0);"
        assert read(out).seed_node.edge.length is None

    def test_input_tree_is_not_changed_by_extraction(self, run, write):
        path = write("u.tre", "((W:1,V:1):1,Q:1,R:1);\n")
        run("mrca", "-n", "W", "Q", "R", "--extract", "-i", str(path))
        assert path.read_text() == "((W:1,V:1):1,Q:1,R:1);\n"

    def test_clade_tips_match_the_table(self, run, write):
        rng = random.Random(3)
        trees = "".join(random_tree(rng, 12, True)[0] + "\n" for _ in range(20))
        path = write("r.tre", trees)
        table = rows(run("mrca", "-n", "t1", "t2", "t3", "-i", str(path)).stdout)
        clades = tree_lines(run("mrca", "-n", "t1", "t2", "t3", "--extract", "-i", str(path)).stdout)
        assert len(clades) == len(table)
        for row, clade in zip(table, clades):
            intruders = row["intruders"].split(",") if row["intruders"] else []
            assert leaves(clade) == sorted(["t1", "t2", "t3"] + intruders)
