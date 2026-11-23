import math
import statistics

import dendropy
import pytest

import treekit
from treekit import FuncTree


def make_ft(**overrides):
    kwargs = dict(command="summary", in_files=[], out_file=None,
                  out_format=None, cores=1, details_out=False)
    kwargs.update(overrides)
    return FuncTree(**kwargs)


def parse_tree(newick):
    return dendropy.Tree.get(data=newick, schema="newick")


def columns(stdout):
    lines = stdout.splitlines()
    return dict(zip(lines[0].split("\t"), lines[1].split("\t")))


# ---------------------------------------------------------------- A1: *_sd columns

class TestStandardDeviationColumns:

    def test_brief_header_names(self, run, write):
        path = write("a.tre", "((A:1,B:1):1,C:2);\n")
        header = run("summary", "-i", str(path)).stdout.splitlines()[0].split("\t")
        assert "tree_length_sd" in header
        assert "tree_length_variance" not in header

    def test_detail_header_names(self, run, write):
        path = write("a.tre", "((A:1,B:1):1,C:2);\n")
        header = run("summary", "-e", "-i", str(path)).stdout.splitlines()[0].split("\t")
        assert "leaf_height_sd" in header
        assert "leaf_height_variance" not in header

    def test_tree_length_sd_is_a_standard_deviation(self, run, write):
        # tree lengths 5 and 4.5: population sd 0.25 (the variance would be 0.0625)
        path = write("two.tre", "((A:1,B:1):1,C:2);\n((A:1,B:1):1,C:1.5);\n")
        result = run("summary", "-i", str(path))
        assert columns(result.stdout)["tree_length_sd"] == "0.25"

    def test_leaf_height_sd_is_a_standard_deviation(self, run, write):
        # leaf heights 2, 2, 1.5 -> population sd, not variance
        path = write("a.tre", "((A:1,B:1):1,C:1.5);\n")
        result = run("summary", "-e", "-i", str(path))
        expected = round(statistics.pstdev([2.0, 2.0, 1.5]), 9)
        assert float(columns(result.stdout)["leaf_height_sd"]) == pytest.approx(expected)


# ---------------------------------------------------------------- A2: is_binary

class TestIsBinary:

    @pytest.mark.parametrize("newick, expected", [
        ("((A:1,B:1):1,(C:1,D:1):1);", True),             # rooted, bifurcating
        ("(A:1,(B:1,(C:1,D:1):1):1);", True),             # rooted, caterpillar
        ("((A:1,B:1):1,C:1,D:1);", True),                 # unrooted, 4 tips, one internal edge
        ("(A:1,B:1,(C:1,D:1):1);", True),                 # unrooted, bifurcating below the trifurcating root
        ("(A:1,B:1,C:1);", True),                         # unrooted, smallest binary tree
        ("(A:1,B:1,C:1,D:1);", False),                    # star tree with 4 tips
        ("(A:1,B:1,C:1,D:1,E:1);", False),                # star tree with 5 tips
        ("((A:1,B:1,C:1):1,(D:1,E:1,F:1):1,G:1);", False),  # several multifurcations
        ("((A:1,B:1,C:1):1,D:1);", False),                # rooted with a multifurcation
    ])
    def test_cases(self, newick, expected):
        assert make_ft().is_binary(parse_tree(newick)) is expected

    def test_star_trees_are_not_binary_in_summary(self, run, write):
        path = write("star.tre", "(A:1,B:1,C:1,D:1);\n(A:1,B:1,C:1);\n")
        result = run("summary", "-i", str(path))
        assert columns(result.stdout)["num_binary"] == "1"

    def test_returns_a_bool_for_degenerate_trees(self):
        # a root with a single child is neither binary nor rooted, and nothing is printed
        tree = parse_tree("((A:1,B:1));")
        assert make_ft().is_binary(tree) is False


class TestIsRooted:

    def test_values(self):
        ft = make_ft()
        assert ft.is_rooted(parse_tree("((A:1,B:1):1,C:2);")) is True
        assert ft.is_rooted(parse_tree("(A:1,B:1,C:1);")) is False
        assert ft.is_rooted(parse_tree("((A:1,B:1));")) is False

    def test_single_child_root_prints_nothing_and_is_not_none(self, run, write):
        path = write("chain.tre", "((A:1,B:1));\n")
        result = run("summary", "-e", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "Invalid tree object" not in result.stdout
        row = result.stdout.splitlines()[1].split("\t")
        assert row[3] == "False"       # is_rooted
        assert row[7] == "NA"          # avg_leaf_height


# ---------------------------------------------------------------- A3: rounding, not truncation

class TestDecimalCorrect:

    @pytest.mark.parametrize("value", [520.938418, 539.617448, 518.678284, 0.57, 4.35])
    def test_short_decimals_are_unchanged(self, value):
        # truncation used to turn 520.938418 into 520.938417999
        assert treekit.decimal_correct(value) == value

    def test_rounds_to_nine_places(self):
        assert treekit.decimal_correct(1.23456789049) == 1.23456789
        assert treekit.decimal_correct(1.23456789051) == 1.234567891

    def test_non_finite_values_do_not_raise(self):
        assert math.isnan(treekit.decimal_correct(float("nan")))
        assert treekit.decimal_correct(float("inf")) == float("inf")

    def test_many_six_decimal_values_survive(self):
        import random
        rng = random.Random(1)
        for _ in range(20000):
            value = round(rng.uniform(0, 1000), 6)
            assert treekit.decimal_correct(value) == value

    def test_summary_prints_the_value_as_written(self, run, write):
        path = write("a.tre", "((A:260.469209,B:260.469209):260.469209,C:520.938418);\n")
        result = run("summary", "-e", "-i", str(path))
        assert columns(result.stdout)["avg_leaf_height"] == "520.938418"
