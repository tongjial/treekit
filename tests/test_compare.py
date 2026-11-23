import random

import dendropy
import pytest
from dendropy.calculate import treecompare

import treekit

HEADER = ["source_a", "index_a", "source_b", "index_b", "common_taxa",
          "only_a", "only_b", "rf", "norm_rf", "weighted_rf"]


def rows(stdout):
    lines = stdout.splitlines()
    assert lines[0].split("\t") == HEADER
    return [dict(zip(HEADER, line.split("\t"))) for line in lines[1:]]


def random_tree(rng, labels, rooted, lengths=True):
    nodes = [("%s:%.4f" % (l, rng.uniform(0.1, 3))) if lengths else l for l in labels]
    stop = 2 if rooted else 3
    while len(nodes) > stop:
        a = nodes.pop(rng.randrange(len(nodes)))
        b = nodes.pop(rng.randrange(len(nodes)))
        nodes.append(("(%s,%s):%.4f" % (a, b, rng.uniform(0.1, 3))) if lengths
                     else "(%s,%s)" % (a, b))
    return "(%s);" % ",".join(nodes)


def dendropy_distances(newick_a, newick_b, rooted):
    """Reference values: both trees reduced to their common taxa, then dendropy."""
    ns = dendropy.TaxonNamespace(is_case_sensitive=True)
    ta = dendropy.Tree.get(data=newick_a, schema="newick", taxon_namespace=ns,
                           case_sensitive_taxon_labels=True)
    tb = dendropy.Tree.get(data=newick_b, schema="newick", taxon_namespace=ns,
                           case_sensitive_taxon_labels=True)
    common = ({l.taxon.label for l in ta.leaf_node_iter()}
              & {l.taxon.label for l in tb.leaf_node_iter()})
    ta = ta.extract_tree_with_taxa_labels(common)
    tb = tb.extract_tree_with_taxa_labels(common)
    for t in (ta, tb):
        t.is_rooted = rooted
        # extraction leaves a "stem" above the new root; it is not a branch of the
        # reduced tree, but dendropy would count its length in the weighted RF
        t.seed_node.edge.length = None
    rf = treecompare.symmetric_difference(ta, tb)
    weighted = treecompare.weighted_robinson_foulds_distance(ta, tb)
    return rf, weighted


class TestCompare:

    def test_identical_trees(self, run, write):
        a = write("a.tre", "((A:1,B:1):1,(C:1,D:1):1,E:1);\n")
        result = run("compare", "-i", str(a), "-r", str(a).replace("a.tre", "a.tre"))
        # comparing a file with itself is allowed (-r is only read)
        assert result.returncode == 0, result.stderr
        [row] = rows(result.stdout)
        assert (row["rf"], row["norm_rf"], row["weighted_rf"]) == ("0", "0.0", "0.0")

    def test_all_pairs(self, run, write):
        write("a.tre", "((A:1,B:1):1,(C:1,D:1):1,E:1);\n((A:1,C:1):1,(B:1,D:1):1,E:1);\n")
        write("b.tre", "((A:1,B:1):1,(C:1,D:1):1,E:1);\n((A:1,D:1):1,(B:1,C:1):1,E:1);\n"
                       "((A:1,B:1):1,(C:1,E:1):1,D:1);\n")
        table = rows(run("compare", "-i", "a.tre", "-r", "b.tre").stdout)
        assert [(r["index_a"], r["index_b"]) for r in table] == [
            (a, b) for a in "12" for b in "123"]

    def test_different_taxa_are_pruned_to_the_common_ones(self, run, write):
        # without F, b is the same topology as a; E's branch becomes 1 + 1
        write("a.tre", "((A:1,B:1):1,(C:1,D:1):1,E:1);\n")
        write("b.tre", "((A:1,B:1):1,(C:1,D:1):1,(E:1,F:1):1);\n")
        [row] = rows(run("compare", "-i", "a.tre", "-r", "b.tre").stdout)
        assert (row["common_taxa"], row["only_a"], row["only_b"]) == ("5", "0", "1")
        assert (row["rf"], row["norm_rf"], row["weighted_rf"]) == ("0", "0.0", "1.0")

    def test_known_distances(self, run, write):
        write("a.tre", "((A:1,C:1):1,(B:1,D:1):1,E:1);\n")
        write("b.tre", "((A:1,B:1):1,(C:1,D:1):1,E:2);\n")
        [row] = rows(run("compare", "-i", "a.tre", "-r", "b.tre").stdout)
        assert (row["rf"], row["norm_rf"], row["weighted_rf"]) == ("4", "1.0", "5.0")

    def test_unrooted_by_default_rooted_on_request(self, run, write):
        # the same unrooted tree, rooted in two places
        write("a.tre", "[&R] (((A:1,B:1):1,C:1):1,(D:1,E:1):1);\n")
        write("b.tre", "[&R] ((A:1,B:1):1,(C:1,(D:1,E:1):1):1);\n")
        [unrooted] = rows(run("compare", "-i", "a.tre", "-r", "b.tre").stdout)
        [rooted] = rows(run("compare", "--rooted", "-i", "a.tre", "-r", "b.tre").stdout)
        assert unrooted["rf"] == "0"
        assert rooted["rf"] == "2"

    def test_rooted_mode_warns_about_unrooted_trees(self, run, write):
        write("a.tre", "((A:1,B:1):1,(C:1,D:1):1,E:1);\n")
        result = run("compare", "--rooted", "-i", "a.tre", "-r", "a.tre")
        assert "have no bifurcating root" in result.stderr

    def test_no_branch_lengths_gives_na_weighted(self, run, write):
        write("a.tre", "((A,B),(C,D),E);\n")
        write("b.tre", "((A:1,B:1):1,(C:1,D:1):1,E:1);\n")
        [row] = rows(run("compare", "-i", "a.tre", "-r", "b.tre").stdout)
        assert row["rf"] == "0" and row["weighted_rf"] == "NA"

    def test_too_few_common_taxa(self, run, write):
        write("a.tre", "((A:1,B:1):1,(C:1,X:1):1,Y:1);\n")
        write("b.tre", "((A:1,B:1):1,(C:1,P:1):1,Q:1);\n")
        result = run("compare", "-i", "a.tre", "-r", "b.tre")
        [row] = rows(result.stdout)
        assert row["common_taxa"] == "3" and row["norm_rf"] == "NA"
        assert "1 pair(s) share fewer than 4 taxa" in result.stderr

    def test_polytomies(self, run, write):
        # a star tree has no split: rf = number of splits of the other tree
        write("a.tre", "(A:1,B:1,C:1,D:1,E:1);\n")
        write("b.tre", "((A:1,B:1):1,(C:1,D:1):1,E:1);\n")
        [row] = rows(run("compare", "-i", "a.tre", "-r", "b.tre").stdout)
        assert (row["rf"], row["norm_rf"]) == ("2", "1.0")

    @pytest.mark.parametrize("rooted", [False, True])
    def test_agrees_with_dendropy_on_random_trees(self, rooted):
        rng = random.Random(21 + rooted)
        for _ in range(150):
            labels = ["t%d" % i for i in range(rng.randint(5, 25))]
            # different, overlapping taxon sets in the two trees
            la = [l for l in labels if rng.random() < 0.85] or labels[:5]
            lb = [l for l in labels if rng.random() < 0.85] or labels[:5]
            if len(set(la) & set(lb)) < 4:
                continue
            na, nb = random_tree(rng, la, rooted), random_tree(rng, lb, rooted)
            if rng.random() < 0.3:
                nb = na                                   # also test equal trees
            rf_ref, weighted_ref = dendropy_distances(na, nb, rooted)

            bit = {}
            shapes = []
            for newick in (na, nb):
                tree = dendropy.Tree.get(data=newick, schema="newick")
                nodes = list(tree.preorder_node_iter())
                pos = {n: i for i, n in enumerate(nodes)}
                parents = [-1 if n is tree.seed_node else pos[n.parent_node] for n in nodes]
                lengths = [None if n is tree.seed_node else float(n.edge.length) for n in nodes]
                names = [n.taxon.label if n.is_leaf() else None for n in nodes]
                for name in names:
                    if name is not None and name not in bit:
                        bit[name] = 1 << len(bit)
                shapes.append(("x", 1, parents, lengths, names, rooted))
            a, b = (treekit.TreeSplits(s, bit) for s in shapes)
            common = a.taxa & b.taxa
            sa, wa = a.restricted(common, rooted)
            sb, wb = b.restricted(common, rooted)
            assert len(sa ^ sb) == rf_ref, (na, nb)
            weighted = sum(abs(wa.get(k, 0) - wb.get(k, 0)) for k in set(wa) | set(wb))
            assert weighted == pytest.approx(weighted_ref), (na, nb)

    def test_mixed_formats_and_several_files(self, run, write):
        write("a.tre", "((A:1,B:1):1,(C:1,D:1):1,E:1);\n")
        write("n.nex", "#NEXUS\nBegin trees;\nTree t1 = ((A:1,B:1):1,(C:1,D:1):1,E:1);\n"
                       "Tree t2 = ((A:1,C:1):1,(B:1,D:1):1,E:1);\nEnd;\n")
        write("c.tre", "((A:1,B:1):1,(C:1,D:1):1,E:1);\n")
        table = rows(run("compare", "-i", "a.tre", "c.tre", "-r", "n.nex").stdout)
        assert [(r["source_a"], r["index_b"], r["rf"]) for r in table] == [
            ("a.tre", "1", "0"), ("a.tre", "2", "4"), ("c.tre", "1", "0"), ("c.tre", "2", "4")]

    def test_parallel_gives_the_same_table(self, run, write):
        rng = random.Random(1)
        labels = ["t%d" % i for i in range(10)]
        write("a1.tre", "".join(random_tree(rng, labels, False) + "\n" for _ in range(5)))
        write("a2.tre", "".join(random_tree(rng, labels, False) + "\n" for _ in range(5)))
        write("b.tre", "".join(random_tree(rng, labels, False) + "\n" for _ in range(4)))
        serial = run("compare", "-i", "a1.tre", "a2.tre", "-r", "b.tre")
        parallel = run("compare", "-c", "2", "-i", "a1.tre", "a2.tre", "-r", "b.tre")
        assert parallel.returncode == 0, parallel.stderr
        assert parallel.stdout == serial.stdout
        assert len(rows(serial.stdout)) == 40

    def test_reference_is_required(self, run, write):
        write("a.tre", "((A:1,B:1):1,(C:1,D:1):1,E:1);\n")
        assert run("compare", "-i", "a.tre").returncode == 2

    def test_output_cannot_overwrite_the_reference(self, run, write):
        write("a.tre", "((A:1,B:1):1,(C:1,D:1):1,E:1);\n")
        b = write("b.tre", "((A:1,B:1):1,(C:1,D:1):1,E:1);\n")
        result = run("compare", "-i", "a.tre", "-r", "b.tre", "-o", "b.tre")
        assert result.returncode == 1
        assert "is also an input file" in result.stderr
        assert b.read_text() == "((A:1,B:1):1,(C:1,D:1):1,E:1);\n"

    def test_log_cannot_be_the_reference(self, run, write):
        write("a.tre", "((A:1,B:1):1,(C:1,D:1):1,E:1);\n")
        b = write("b.tre", "((A:1,B:1):1,(C:1,D:1):1,E:1);\n")
        result = run("compare", "-i", "a.tre", "-r", "b.tre", "-l", "b.tre")
        assert result.returncode == 1
        assert b.read_text() == "((A:1,B:1):1,(C:1,D:1):1,E:1);\n"

    def test_unlabelled_tips_are_rejected(self, run, write):
        write("a.tre", "((A:1,:1):1,(C:1,D:1):1,E:1);\n")
        write("b.tre", "((A:1,B:1):1,(C:1,D:1):1,E:1);\n")
        result = run("compare", "-i", "a.tre", "-r", "b.tre")
        assert result.returncode == 1
        assert "tip without a label" in result.stderr
