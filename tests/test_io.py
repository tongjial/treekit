import subprocess
import sys

import pytest

from conftest import TREEKIT_CMD

TREE = "((A:1,B:1):1,C:2);\n"
MULTI = "((A:1,B:1):1,C:2);\n((A:1,C:1):1,D:2);\n((B:1,C:1):1,D:2);\n"
NEXUS_NO_TREES = "#NEXUS\nBegin taxa;\nDimensions ntax=2;\nTaxlabels A B;\nEnd;\n"
NEXUS_TREES = (
    "#NEXUS\nBegin taxa;\nDimensions ntax=3;\nTaxlabels A B C;\nEnd;\n"
    "Begin trees;\nTree t1 = ((A:1,B:1):1,C:2);\nEnd;\n"
)


def no_traceback(result):
    assert "Traceback" not in result.stderr, result.stderr


# ---------------------------------------------------------------- B1 / B2: output streams

class TestStreams:

    def test_core_count_notice_goes_to_stderr(self, run, write):
        # asking for far more cores than any machine has
        path = write("a.tre", TREE)
        result = run("summary", "-c", "9999", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "warning: this computer has" in result.stderr
        assert "fewer than the 9999 requested" in result.stderr
        # stdout holds only the table: header + one row
        lines = result.stdout.splitlines()
        assert len(lines) == 2 and lines[0].startswith("file_name\t")

    def test_overwrite_warning_goes_to_stderr(self, run, write, tmp_path):
        path = write("a.tre", TREE)
        out = tmp_path / "out.tsv"
        first = run("summary", "-o", str(out), "-i", str(path))
        second = run("summary", "-o", str(out), "-i", str(path))
        assert first.stderr == "" and first.stdout == ""
        assert "warning: overwriting" in second.stderr
        assert second.stdout == ""

    def test_output_file_ends_with_a_newline(self, run, write, tmp_path):
        path = write("a.tre", TREE)
        out = tmp_path / "out.tsv"
        run("summary", "-o", str(out), "-i", str(path))
        text = out.read_text()
        assert text.endswith("\n") and not text.endswith("\n\n")

    def test_stdout_ends_with_exactly_one_newline(self, run, write):
        path = write("a.tre", TREE)
        result = run("remove", "-n", "A", "-i", str(path))
        assert result.stdout.endswith("\n") and not result.stdout.endswith("\n\n")

    def test_closed_pipe_exits_quietly(self, write, tmp_path):
        # "treekit ... | head -n 1" must not print a BrokenPipeError traceback
        path = write("many.tre", "((A:1,B:1):1,C:2);\n" * 6000)
        proc = subprocess.Popen(
            [*TREEKIT_CMD, "summary", "-e", "-i", str(path)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )
        proc.stdout.readline()
        proc.stdout.close()
        stderr = proc.stderr.read()
        proc.wait(timeout=60)
        assert "Traceback" not in stderr and "BrokenPipe" not in stderr


# ---------------------------------------------------------------- B3: output over input

class TestOutputOverInput:

    def test_out_file_equal_to_input_is_refused(self, run, write):
        path = write("a.tre", TREE)
        result = run("remove", "-n", "A", "-o", str(path), "-i", str(path))
        assert result.returncode == 1
        assert "-o file" in result.stderr and "is also an input file" in result.stderr
        assert path.read_text() == TREE          # the input is untouched

    def test_out_file_equal_to_input_through_a_relative_path(self, run, write, tmp_path):
        path = write("a.tre", TREE)
        result = run("remove", "-n", "A", "-o", "./a.tre", "-i", "a.tre")
        assert result.returncode == 1
        assert path.read_text() == TREE

    def test_log_file_equal_to_input_is_refused(self, run, write):
        path = write("a.tre", TREE)
        result = run("summary", "-l", str(path), "-i", str(path))
        assert result.returncode == 1
        assert "-l file" in result.stderr
        # the check must come before the log is opened, which would append to the input
        assert path.read_text() == TREE

    def test_log_file_equal_to_input_through_a_wildcard(self, run, write):
        path = write("a.tre", TREE)
        result = run("summary", "-l", "a.tre", "-i", "a.*")
        assert result.returncode == 1
        assert path.read_text() == TREE


# ---------------------------------------------------------------- B4: unwritable output

class TestUnwritableOutput:

    @pytest.mark.parametrize("command", [
        ("remove", "-n", "A"), ("summary",), ("taxa",),
    ])
    def test_unwritable_out_file_is_an_error(self, run, write, tmp_path, command):
        path = write("a.tre", TREE)
        # a folder that does not exist; the message shows the resolved path, which is
        # not the literal argument on every system (/x is C:\\x on Windows)
        target = tmp_path / "no" / "such" / "dir" / "out.txt"
        result = run(*command, "-o", str(target), "-i", str(path))
        assert result.returncode == 1
        assert "treekit: error: cannot write '" in result.stderr
        assert "out.txt" in result.stderr
        assert not target.exists()
        no_traceback(result)


# ---------------------------------------------------------------- B5 / B6: bad input files

class TestBadInput:

    def test_malformed_newick_is_an_error(self, run, write):
        path = write("bad.tre", "((A:1,B:1):1,C:2;\n")
        result = run("summary", "-i", str(path))
        assert result.returncode == 1
        assert "cannot parse" in result.stderr and "as newick" in result.stderr
        no_traceback(result)

    @pytest.mark.parametrize("command", [
        ("summary",), ("summary", "-e"), ("taxa",), ("remove", "-n", "A"),
    ])
    def test_malformed_input_in_every_command(self, run, write, command):
        path = write("bad.tre", "((A:1,B:1):1,C:2;\n")
        result = run(*command, "-i", str(path))
        assert result.returncode == 1
        assert "cannot parse" in result.stderr
        no_traceback(result)

    def test_non_utf8_is_an_error(self, run, tmp_path):
        path = tmp_path / "latin.tre"
        path.write_bytes(b"((A\xe9:1,B:1):1,C:2);\n")
        result = run("summary", "-i", str(path))
        assert result.returncode == 1
        assert "not valid UTF-8" in result.stderr
        no_traceback(result)

    def test_bom_gets_a_clear_message(self, run, tmp_path):
        path = tmp_path / "bom.tre"
        path.write_bytes(b"\xef\xbb\xbf" + TREE.encode())
        result = run("summary", "-i", str(path))
        assert result.returncode == 1
        assert "byte-order mark" in result.stderr
        no_traceback(result)

    def test_parse_error_inside_a_worker_is_reported(self, run, write):
        good = write("a.tre", TREE)
        bad = write("b.tre", "((A:1,B:1):1,C:2;\n")
        result = run("summary", "-c", "2", "-i", str(good), str(bad))
        assert result.returncode == 1
        assert "cannot parse" in result.stderr
        no_traceback(result)

    def test_summary_details_with_no_trees_is_an_error(self, run, write):
        path = write("nt.nex", NEXUS_NO_TREES)
        result = run("summary", "-e", "-i", str(path))
        assert result.returncode == 1
        assert "no trees found" in result.stderr
        assert result.stdout == ""


# ---------------------------------------------------------------- B7: -c validation

class TestCores:

    @pytest.mark.parametrize("value", ["0", "-1", "abc", "1.5"])
    def test_invalid_core_counts_are_rejected(self, run, write, value):
        path = write("a.tre", TREE)
        result = run("summary", "-c", value, "-i", str(path))
        assert result.returncode == 2
        assert "--cores" in result.stderr or "-c/--cores" in result.stderr
        no_traceback(result)

    def test_details_mode_uses_the_cores_too(self, run, write):
        write("a.tre", MULTI)
        write("b.tre", TREE)
        serial = run("summary", "-e", "-c", "1", "-i", "a.tre", "b.tre")
        parallel = run("summary", "-e", "-c", "2", "-i", "a.tre", "b.tre")
        assert parallel.returncode == 0, parallel.stderr
        assert parallel.stdout == serial.stdout
        assert len(serial.stdout.splitlines()) == 1 + 3 + 1   # header + 3 trees + 1 tree


# ---------------------------------------------------------------- C: input formats and paths

class TestInputFormats:

    def test_leading_comment_before_newick(self, run, write):
        # remove writes "[&R] (...)" for some inputs; that must be readable again
        path = write("c.tre", "[&R] ((A:1,B:1):1,C:2);\n")
        result = run("summary", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert columns_of(result)["num_tree"] == "1"

    def test_multiline_comment_before_newick(self, run, write):
        path = write("c.tre", "\n[a comment\nover two lines]\n((A:1,B:1):1,C:2);\n")
        result = run("summary", "-i", str(path))
        assert result.returncode == 0, result.stderr

    def test_single_leaf_newick(self, run, write):
        path = write("one.tre", "A;\n")
        result = run("summary", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert columns_of(result)["num_tree"] == "1"

    def test_unterminated_comment_is_an_error(self, run, write):
        path = write("c.tre", "[never closed ((A:1,B:1):1,C:2);\n")
        result = run("summary", "-i", str(path))
        assert result.returncode == 1
        assert "not a newick or nexus" in result.stderr

    def test_bracket_in_file_name(self, run, write):
        path = write("trees[1].tre", TREE)
        result = run("summary", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stdout.splitlines()[1].startswith("trees[1].tre\t")

    def test_bracket_in_relative_file_name(self, run, write):
        write("trees[1].tre", TREE)
        result = run("summary", "-i", "trees[1].tre")
        assert result.returncode == 0, result.stderr

    def test_wildcards_still_work(self, run, write):
        write("a1.tre", TREE)
        write("a2.tre", TREE)
        result = run("summary", "-i", "a?.tre")
        assert [line.split("\t")[0] for line in result.stdout.splitlines()[1:]] == ["a1.tre", "a2.tre"]

    def test_remove_output_can_be_read_back_for_nexus_input(self, run, write, tmp_path):
        path = write("n.nex", NEXUS_TREES)
        out = tmp_path / "out.tre"
        result = run("remove", "-n", "A", "-o", str(out), "-i", str(path))
        assert result.returncode == 0, result.stderr
        again = run("summary", "-i", str(out))
        assert again.returncode == 0, again.stderr
        assert columns_of(again)["num_tree"] == "1"

    def test_nexus_output_format(self, run, write, tmp_path):
        path = write("m.tre", MULTI)
        out = tmp_path / "out.nex"
        result = run("remove", "-n", "D", "-u", "nexus", "-o", str(out), "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert out.read_text().startswith("#NEXUS")
        again = run("summary", "-i", str(out))
        assert again.returncode == 0, again.stderr
        assert columns_of(again)["num_tree"] == "3"

    def test_nexus_output_to_stdout(self, run, write):
        path = write("m.tre", MULTI)
        result = run("remove", "-n", "D", "-u", "nexus", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stdout.startswith("#NEXUS")


# ---------------------------------------------------------------- the log option on every command

class TestLogOnAllCommands:

    @pytest.mark.parametrize("command", [("summary",), ("summary", "-e"), ("taxa",)])
    def test_log_option(self, run, write, tmp_path, command):
        path = write("a.tre", TREE)
        log = tmp_path / "run.log"
        result = run(*command, "-l", str(log), "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "run: treekit " + command[0] in log.read_text()

    def test_wrote_message_goes_to_the_log_only(self, run, write, tmp_path):
        path = write("a.tre", TREE)
        out = tmp_path / "out.tsv"
        log = tmp_path / "run.log"
        result = run("summary", "-o", str(out), "-l", str(log), "-i", str(path))
        assert result.stdout == "" and result.stderr == ""
        assert "wrote summaries to" in log.read_text()


def columns_of(result):
    lines = result.stdout.splitlines()
    return dict(zip(lines[0].split("\t"), lines[1].split("\t")))


# ---------------------------------------------------------------- non-ASCII names, any default encoding

# the default encoding of text files differs between systems (cp1252 on Windows);
# these variables ask for an ASCII default on every system that honours them
ASCII_DEFAULT = {"LC_ALL": "C", "PYTHONCOERCECLOCALE": "0", "PYTHONUTF8": "0",
                 "PYTHONIOENCODING": ""}


class TestNonAsciiNames:

    NEWICK = "((Café_é:1,'Zoë b':1):1,C:1,D:1);\n"
    NEXUS = ("#NEXUS\nBegin taxa;\nDimensions ntax=3;\nTaxlabels Café 'Zoë b' C;\nEnd;\n"
             "Begin trees;\nTree t = ((Café:1,'Zoë b':1):1,C:1);\nEnd;\n")

    @pytest.mark.parametrize("env", [None, ASCII_DEFAULT], ids=["default", "ascii-default"])
    def test_newick_names_are_read_as_utf8(self, run, tmp_path, env):
        path = tmp_path / "u.tre"
        path.write_bytes(self.NEWICK.encode("utf-8"))
        result = run("taxa", "-i", str(path), env=env)
        assert result.returncode == 0, result.stderr
        assert result.stdout.splitlines() == ["C", "Café_é", "D", "Zoë b"]

    @pytest.mark.parametrize("env", [None, ASCII_DEFAULT], ids=["default", "ascii-default"])
    def test_nexus_names_are_read_as_utf8(self, run, tmp_path, env):
        path = tmp_path / "u.nex"
        path.write_bytes(self.NEXUS.encode("utf-8"))
        result = run("taxa", "-i", str(path), env=env)
        assert result.returncode == 0, result.stderr
        assert result.stdout.splitlines() == ["C", "Café", "Zoë b"]

    @pytest.mark.parametrize("env", [None, ASCII_DEFAULT], ids=["default", "ascii-default"])
    def test_names_survive_a_remove_round_trip(self, run, tmp_path, env):
        path = tmp_path / "u.tre"
        path.write_bytes(self.NEWICK.encode("utf-8"))
        out = tmp_path / "out.tre"
        result = run("remove", "-n", "D", "-o", str(out), "-i", str(path), env=env)
        assert result.returncode == 0, result.stderr
        assert out.read_bytes().decode("utf-8").strip() == \
            "((Café_é:1.0,'Zoë b':1.0):1.0,C:1.0);"

    def test_a_name_given_on_the_command_line(self, run, tmp_path):
        path = tmp_path / "u.tre"
        path.write_bytes(self.NEWICK.encode("utf-8"))
        result = run("remove", "-n", "Café_é", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "Café" not in result.stdout

    def test_other_encodings_are_rejected_not_garbled(self, run, tmp_path):
        path = tmp_path / "l1.tre"
        path.write_bytes("((Café:1,B:1):1,C:1);\n".encode("latin-1"))
        result = run("taxa", "-i", str(path), env=ASCII_DEFAULT)
        assert result.returncode == 1
        assert "is not valid UTF-8 text" in result.stderr


@pytest.mark.skipif(sys.platform == "win32" or sys.platform == "darwin",
                    reason="needs file names that are not valid text (Linux)")
class TestFileNamesThatAreNotText:

    def make(self, tmp_path):
        import os
        name = os.fsdecode(b"trees_\xe9.tre")          # surrogate escape for the byte e9
        path = tmp_path / name
        path.write_text("((A:1,B:1):1,C:1,D:1);\n")
        return path

    def test_summary_prints_the_name_instead_of_crashing(self, run, tmp_path):
        path = self.make(tmp_path)
        result = run("summary", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "Traceback" not in result.stderr
        assert "trees_\\udce9.tre" in result.stdout

    def test_the_output_is_valid_utf8(self, tmp_path):
        import subprocess
        path = self.make(tmp_path)
        out = subprocess.run([*TREEKIT_CMD, "summary", "-i", str(path)], capture_output=True)
        assert out.returncode == 0
        out.stdout.decode("utf-8")                     # raises if the bytes are not UTF-8

    def test_an_error_message_with_such_a_name_does_not_crash(self, run, tmp_path):
        result = run("remove", "-n", "A", "-i", str(tmp_path / "missing_\udce9.tre"))
        assert result.returncode == 1
        assert "Traceback" not in result.stderr
