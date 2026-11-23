import re

import pytest

TREES = "((A:1,B:3)90:1,(C:0.5,D:0.5):2,E:1);\n((A,B),(C,D),E);\nA:1;\n"


def blocks(stdout):
    """Drawings by title line."""
    out = {}
    current = None
    for line in stdout.splitlines():
        if line.startswith("tree "):
            current = line
            out[current] = []
        elif current is not None and line:
            out[current].append(line)
    return out


class TestDraw:

    def test_every_tree_has_a_title(self, run, write):
        path = write("d.tre", TREES)
        result = run("draw", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert list(blocks(result.stdout)) == [
            "tree 1 (d.tre, tree 1)", "tree 2 (d.tre, tree 2)", "tree 3 (d.tre, tree 3)"]

    def test_every_tip_is_drawn(self, run, write):
        path = write("d.tre", TREES)
        drawing = blocks(run("draw", "-t", "1", "-i", str(path)).stdout)["tree 1 (d.tre, tree 1)"]
        for name in "ABCDE":
            assert sum(line.endswith(" " + name) for line in drawing) == 1

    def test_select_trees(self, run, write):
        path = write("d.tre", TREES)
        result = run("draw", "-t", "3", "1", "-i", str(path))
        assert list(blocks(result.stdout)) == ["tree 3 (d.tre, tree 3)", "tree 1 (d.tre, tree 1)"]

    def test_tree_numbers_count_over_all_files(self, run, write):
        write("a.tre", "((A:1,B:1):1,C:1);\n")
        write("b.tre", "((X:1,Y:1):1,Z:1);\n((P:1,Q:1):1,R:1);\n")
        result = run("draw", "-t", "3", "-i", "a.tre", "b.tre")
        assert list(blocks(result.stdout)) == ["tree 3 (b.tre, tree 2)"]

    def test_tree_number_out_of_range(self, run, write):
        path = write("d.tre", TREES)
        result = run("draw", "-t", "2", "9", "-i", str(path))
        assert result.returncode == 1
        assert "there are 3 input trees; no tree 9" in result.stderr

    def test_width(self, run, write):
        path = write("d.tre", TREES)
        result = run("draw", "-w", "50", "-t", "1", "2", "-i", str(path))
        lines = [l for l in result.stdout.splitlines() if not l.startswith("tree ")]
        assert max(len(l) for l in lines) <= 50

    def test_width_too_small(self, run, write):
        path = write("d.tre", TREES)
        assert run("draw", "-w", "5", "-i", str(path)).returncode == 1

    def test_no_trailing_spaces_or_blank_tail(self, run, write):
        path = write("d.tre", TREES)
        out = run("draw", "-i", str(path)).stdout
        assert all(line == line.rstrip() for line in out.splitlines())
        assert out.endswith("\n") and not out.endswith("\n\n")

    def test_supports_are_shown_only_on_request(self, run, write):
        path = write("d.tre", TREES)
        plain = run("draw", "-t", "1", "-i", str(path)).stdout
        with_supports = run("draw", "-t", "1", "-s", "-i", str(path)).stdout
        assert "90" not in plain
        assert "90" in with_supports
        assert "@" not in with_supports

    def test_length_metric_scales_branches(self, run, write):
        # B's branch (3) is three times A's (1)
        path = write("s.tre", "((A:1,B:3):1,C:1);\n")
        drawing = run("draw", "-m", "length", "-w", "60", "-i", str(path)).stdout
        dashes = {name: len(re.search(r"(-+) " + name + "$", line).group(1))
                  for line in drawing.splitlines() for name in "AB"
                  if line.endswith(" " + name)}
        assert dashes["B"] > 2 * dashes["A"]

    def test_topology_metric_ignores_lengths(self, run, write):
        path = write("s.tre", "((A:1,B:3):1,C:1);\n")
        drawing = run("draw", "-m", "topology", "-w", "60", "-i", str(path)).stdout
        ends = [line.index(" " + name) for line in drawing.splitlines() for name in "AB"
                if line.endswith(" " + name)]
        assert ends[0] == ends[1]

    def test_auto_draws_trees_without_lengths_as_topology(self, run, write, tmp_path):
        path = write("d.tre", TREES)
        log = tmp_path / "run.log"
        result = run("draw", "-t", "2", "-l", str(log), "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "drawn as topology" in log.read_text()

    def test_length_metric_without_lengths_is_an_error(self, run, write):
        path = write("d.tre", TREES)
        result = run("draw", "-m", "length", "-t", "2", "-i", str(path))
        assert result.returncode == 1
        assert "-m length needs a positive length on every branch" in result.stderr

    def test_auto_uses_topology_for_very_uneven_lengths(self, run, write, tmp_path):
        # one branch is 10000 times the others: to scale, the rest would collapse
        path = write("u.tre", "((A:0.001,B:0.001):0.001,(C:0.001,D:0.001):10);\n")
        log = tmp_path / "run.log"
        run("draw", "-l", str(log), "-i", str(path))
        assert "too uneven" in log.read_text()

    def test_single_tip_tree(self, run, write):
        path = write("d.tre", TREES)
        drawing = blocks(run("draw", "-t", "3", "-i", str(path)).stdout)
        assert drawing == {"tree 3 (d.tre, tree 3)": ["--- A"]}

    def test_nexus_tree_names_are_shown(self, run, write):
        path = write("n.nex", "#NEXUS\nBegin trees;\nTree best = ((A:1,B:1):1,C:1);\nEnd;\n")
        assert "tree 1 (n.nex, tree 1) best" in run("draw", "-i", str(path)).stdout

    def test_output_file(self, run, write, tmp_path):
        path = write("d.tre", TREES)
        out = tmp_path / "drawing.txt"
        result = run("draw", "-o", str(out), "-i", str(path))
        assert result.stdout == ""
        assert out.read_text().startswith("tree 1 (d.tre, tree 1)")
