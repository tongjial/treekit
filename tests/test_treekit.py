import multiprocessing as mp

import dendropy
import pytest

import treekit
from treekit import FuncTree

# trees used across tests
ULTRA = "((A:1,B:1):1,C:2);"                 # tips all at depth 2
NON_ULTRA = "((A:1,B:1):1,C:1.5);"           # tips at depth 2 and 1.5
NOISE = "((A:1,B:1):1,C:1.9999999999995);"   # tips differ only by float noise
UNROOTED = "(A:1,B:1,C:1);"
NEXUS_NO_TREES = (
    "#NEXUS\n"
    "Begin taxa;\n"
    "\tDimensions ntax=3;\n"
    "\tTaxlabels A B C;\n"
    "End;\n"
)


def make_ft(**overrides):
    """FuncTree with no input files, for calling methods directly."""
    kwargs = dict(command="summary", in_files=[], out_file=None,
                  out_format=None, cores=1, details_out=False)
    kwargs.update(overrides)
    return FuncTree(**kwargs)


def parse_tree(newick):
    return dendropy.Tree.get(data=newick, schema="newick")


# ---------------------------------------------------------------- fix 1: empty input

class TestEmptyInput:

    def test_empty_file_is_an_error_not_a_hang(self, run, write):
        path = write("empty.tre", "")
        # timeout makes a regression (infinite loop) fail the test instead of hanging it
        result = run("summary", "-i", str(path), timeout=20)
        assert result.returncode == 1
        assert "not a newick or nexus tree file" in result.stderr
        assert "Traceback" not in result.stderr

    def test_blank_lines_only_is_an_error(self, run, write):
        path = write("blank.tre", "\n\n   \n\t\n")
        result = run("summary", "-i", str(path), timeout=20)
        assert result.returncode == 1
        assert "not a newick or nexus tree file" in result.stderr

    def test_unrecognized_content_is_an_error(self, run, write):
        path = write("junk.tre", "hello world\n")
        result = run("summary", "-i", str(path), timeout=20)
        assert result.returncode == 1
        assert "not a newick or nexus tree file" in result.stderr

    def test_nexus_without_trees_is_an_error(self, run, write):
        path = write("notrees.nex", NEXUS_NO_TREES)
        result = run("summary", "-i", str(path), timeout=20)
        assert result.returncode == 1
        assert "no trees found" in result.stderr
        assert "Traceback" not in result.stderr

    def test_taxa_common_on_no_trees_is_an_error(self, run, write):
        path = write("notrees.nex", NEXUS_NO_TREES)
        result = run("taxa", "-k", "-i", str(path), timeout=20)
        assert result.returncode == 1
        assert "no trees found" in result.stderr


# ---------------------------------------------------------------- fix 2 and 3: input paths

class TestInputPaths:

    def test_absolute_path_works_from_any_directory(self, run, write, tmp_path):
        path = write("ok.tre", ULTRA + "\n")
        elsewhere = tmp_path / "elsewhere"
        elsewhere.mkdir()
        result = run("summary", "-i", str(path), cwd=elsewhere)
        assert result.returncode == 0, result.stderr
        header, row = result.stdout.splitlines()[:2]
        assert header.startswith("file_name\t")
        assert row.startswith("ok.tre\t1\t")

    def test_relative_wildcard_matches_sorted_files(self, run, write):
        write("b.tre", ULTRA + "\n")
        write("a.tre", NON_ULTRA + "\n")
        result = run("summary", "-i", "*.tre")
        assert result.returncode == 0, result.stderr
        rows = result.stdout.splitlines()[1:]
        assert [row.split("\t")[0] for row in rows] == ["a.tre", "b.tre"]

    def test_missing_file_is_an_error(self, run):
        result = run("summary", "-i", "does_not_exist.tre")
        assert result.returncode == 1
        assert "input file not found: 'does_not_exist.tre'" in result.stderr
        assert "Traceback" not in result.stderr

    def test_wildcard_with_no_match_is_an_error(self, run, write):
        write("a.tre", ULTRA + "\n")
        result = run("summary", "-i", "*.nope")
        assert result.returncode == 1
        assert "input file not found" in result.stderr

    def test_directory_is_rejected(self, run, tmp_path):
        (tmp_path / "subdir").mkdir()
        result = run("summary", "-i", "subdir")
        assert result.returncode == 1
        assert "input is not a file" in result.stderr


# ---------------------------------------------------------------- fix 4: unimplemented commands

class TestCommandDispatch:

    @pytest.mark.parametrize("command", treekit.PLANNED_COMMANDS)
    def test_planned_command_reports_not_implemented(self, run, write, command):
        path = write("ok.tre", ULTRA + "\n")
        result = run(command, "-i", str(path))
        assert result.returncode == 2
        assert "not implemented yet" in result.stderr
        assert "Traceback" not in result.stderr

    def test_unknown_command_is_rejected(self, run):
        result = run("bogus")
        assert result.returncode == 2
        assert "unrecognized command 'bogus'" in result.stderr

    @pytest.mark.parametrize("name", ["get_paths", "write_summaries", "__init__"])
    def test_method_names_are_not_commands(self, run, name):
        # previously any method on ParsedArgs was dispatched as a subcommand
        result = run(name)
        assert result.returncode == 2
        assert "unrecognized command" in result.stderr

    def test_implemented_commands_are_registered(self):
        assert set(treekit.IMPLEMENTED_COMMANDS) == {"summary", "taxa", "remove"}
        assert not set(treekit.IMPLEMENTED_COMMANDS) & set(treekit.PLANNED_COMMANDS)


# ---------------------------------------------------------------- fix 5: ultrametric check

class TestUltrametric:

    def test_exactly_ultrametric_tree(self):
        assert make_ft().is_ultrametric(parse_tree(ULTRA)) is True

    def test_non_ultrametric_tree(self):
        assert make_ft().is_ultrametric(parse_tree(NON_ULTRA)) is False

    def test_float_noise_is_still_ultrametric(self):
        # tip depths differ by 5e-13; the old check truncated to 9 decimals and said False
        assert make_ft().is_ultrametric(parse_tree(NOISE)) is True

    def test_unrooted_tree_is_not_ultrametric(self):
        assert make_ft().is_ultrametric(parse_tree(UNROOTED)) is False

    def test_summary_detail_is_consistent(self, run, write):
        # the sd rounds to 0.0 (all tips equal) so ultrametric must also be True
        path = write("noise.tre", NOISE + "\n")
        result = run("summary", "-e", "-i", str(path))
        assert result.returncode == 0, result.stderr
        columns = dict(zip(result.stdout.splitlines()[0].split("\t"),
                           result.stdout.splitlines()[1].split("\t")))
        assert columns["leaf_height_sd"] == "0.0"
        assert columns["is_ultrametric"] == "True"


# ---------------------------------------------------------------- fix 6: parallel summary

class TestSummaryStats:

    def test_brief_summary_values(self, write):
        path = write("two.tre", ULTRA + "\n" + NON_ULTRA + "\n")
        header, rows = make_ft(in_files=[str(path)]).get_summstats()
        assert header[0] == "file_name"
        # lengths are 5 and 4.5 -> mean 4.75, population stdev 0.25
        assert rows == [("two.tre", "2", "2", "2", "2", "1", "4.75", "0.25", "3", "3")]

    def test_parallel_matches_serial_over_multiple_files(self, run, write):
        write("a.tre", ULTRA + "\n" + NON_ULTRA + "\n")
        write("b.tre", NOISE + "\n")
        write("c.tre", UNROOTED + "\n")
        serial = run("summary", "-c", "1", "-i", "*.tre")
        parallel = run("summary", "-c", "2", "-i", "*.tre")
        assert serial.returncode == 0, serial.stderr
        assert parallel.returncode == 0, parallel.stderr
        assert parallel.stdout == serial.stdout

    def test_multiple_files_use_a_pool(self, write, monkeypatch):
        a = write("a.tre", ULTRA + "\n")
        b = write("b.tre", NON_ULTRA + "\n")
        calls = []
        real_pool = mp.Pool

        def spy_pool(*args, **kwargs):
            calls.append(args)
            return real_pool(*args, **kwargs)

        monkeypatch.setattr(treekit.mp, "Pool", spy_pool)
        ft = make_ft(in_files=[str(a), str(b)])
        ft.cores = 2  # set directly: get_cores caps this at the machine's core count
        _, rows = ft.get_summstats()
        assert len(calls) == 1
        assert [row[0] for row in rows] == ["a.tre", "b.tre"]

    def test_single_file_never_creates_a_pool(self, write, monkeypatch):
        a = write("a.tre", ULTRA + "\n")

        def no_pool(*args, **kwargs):
            raise AssertionError("a pool must not be created for a single file")

        monkeypatch.setattr(treekit.mp, "Pool", no_pool)
        ft = make_ft(in_files=[str(a)])
        ft.cores = 4
        _, rows = ft.get_summstats()
        assert rows[0][0] == "a.tre"
