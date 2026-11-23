import pytest

# three trees sharing most taxa; D appears in only two of them
MULTI = (
    "((A:1,B:1):1,C:2);\n"
    "((A:1,C:1):1,D:2);\n"
    "((B:1,C:1):1,D:2);\n"
)
SINGLE = "((A:1,B:1):1,C:2);\n"
# X is a tip that the quoting tests remove, so they exercise a real removal
QUOTED = "(('A(1)':1,'B,2':1):1,'C:3':2,Badula_crassa_Y24864:1,X:1);\n"


def tree_lines(stdout):
    """Output lines that hold trees (the warnings go to stderr, not here)."""
    return [line for line in stdout.splitlines() if line.strip()]


# ---------------------------------------------------------------- absent taxa: warnings

class TestMissingTaxa:

    def test_absent_taxon_on_single_tree_warns(self, run, write):
        path = write("s.tre", SINGLE)
        result = run("remove", "-n", "A", "ZZZ", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "treekit: warning: 1 taxon was not found in any input tree: ZZZ" in result.stderr
        assert "(B:2.0,C:2.0);" in result.stdout

    def test_absent_taxa_are_summarised_in_one_warning(self, run, write):
        path = write("m.tre", MULTI)
        result = run("remove", "-n", "D", "ZZZ", "YYY", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert ("treekit: warning: 2 taxa were not found in any input tree: ZZZ, YYY"
                in result.stderr)
        assert result.stderr.count("warning") == 1

    def test_absent_taxon_on_multi_tree_warns_once(self, run, write):
        path = write("m.tre", MULTI)
        result = run("remove", "-n", "D", "ZZZ", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stderr.count("not found in any input tree") == 1

    def test_partial_presence_does_not_warn(self, run, write):
        # D is in 2 of 3 trees: normal for gene-tree sets, so no warning
        path = write("m.tre", MULTI)
        result = run("remove", "-n", "D", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stderr == ""
        assert len(tree_lines(result.stdout)) == 3

    def test_partial_presence_is_recorded_in_log(self, run, write, tmp_path):
        path = write("m.tre", MULTI)
        log = tmp_path / "run.log"
        result = run("remove", "-n", "D", "-l", str(log), "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "taxon 'D' removed from 2 of 3 trees" in log.read_text()

    def test_warning_is_also_written_to_log(self, run, write, tmp_path):
        path = write("m.tre", MULTI)
        log = tmp_path / "run.log"
        run("remove", "-n", "D", "ZZZ", "-l", str(log), "-i", str(path))
        assert "warning: 1 taxon was not found in any input tree: ZZZ" in log.read_text()

    def test_duplicate_names_are_handled_once(self, run, write):
        path = write("s.tre", SINGLE)
        result = run("remove", "-n", "A", "ZZZ", "ZZZ", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stderr.count("ZZZ") == 1


# ---------------------------------------------------------------- all taxa absent: error

class TestAllAbsent:

    def test_all_absent_on_multi_tree_is_an_error(self, run, write):
        path = write("m.tre", MULTI)
        result = run("remove", "-n", "ZZZ", "YYY", "-i", str(path))
        assert result.returncode == 1
        assert ("none of the 2 taxa to remove were found in any input tree: ZZZ, YYY"
                in result.stderr)
        assert "warning" not in result.stderr
        assert result.stdout == ""

    def test_all_absent_on_single_tree_is_an_error(self, run, write):
        path = write("s.tre", SINGLE)
        result = run("remove", "-n", "ZZZ", "-i", str(path))
        assert result.returncode == 1
        assert "none of the 1 taxon to remove were found in any input tree: ZZZ" in result.stderr
        assert result.stdout == ""

    def test_all_absent_from_taxa_file_is_an_error(self, run, write):
        path = write("s.tre", SINGLE)
        taxa = write("absent.txt", "ZZZ\nYYY\n")
        result = run("remove", "-f", str(taxa), "-i", str(path))
        assert result.returncode == 1
        assert "none of the 2 taxa" in result.stderr

    def test_all_absent_is_logged_as_error(self, run, write, tmp_path):
        path = write("s.tre", SINGLE)
        log = tmp_path / "run.log"
        run("remove", "-n", "ZZZ", "-l", str(log), "-i", str(path))
        assert "error: none of the 1 taxon to remove were found" in log.read_text()


# ---------------------------------------------------------------- --strict

class TestStrict:

    def test_strict_fails_on_absent_taxon(self, run, write):
        path = write("m.tre", MULTI)
        result = run("remove", "-n", "D", "ZZZ", "--strict", "-i", str(path))
        assert result.returncode == 1
        assert "1 taxon was not found in any input tree: ZZZ (--strict)" in result.stderr
        # nothing is written to stdout when the run fails
        assert result.stdout == ""

    def test_strict_reports_all_absent_taxa_together(self, run, write):
        path = write("m.tre", MULTI)
        result = run("remove", "-n", "D", "ZZZ", "YYY", "--strict", "-i", str(path))
        assert result.returncode == 1
        assert "2 taxa were not found in any input tree: ZZZ, YYY (--strict)" in result.stderr

    def test_strict_allows_partial_presence(self, run, write):
        path = write("m.tre", MULTI)
        result = run("remove", "-n", "D", "--strict", "-i", str(path))
        assert result.returncode == 0, result.stderr

    def test_without_strict_absent_taxon_is_not_fatal(self, run, write):
        path = write("m.tre", MULTI)
        result = run("remove", "-n", "D", "ZZZ", "-i", str(path))
        assert result.returncode == 0


# ---------------------------------------------------------------- trees left without tips

class TestDegenerateTrees:

    def test_removing_all_tips_of_a_tree_is_an_error(self, run, write):
        path = write("two.tre", "((A:1,B:1):1,C:2);\n(A:1,B:1);\n")
        result = run("remove", "-n", "A", "B", "-i", str(path))
        assert result.returncode == 1
        assert "would leave tree 2 with no tips" in result.stderr
        assert "Traceback" not in result.stderr
        assert result.stdout == ""

    def test_removing_all_taxa_from_one_tree_is_an_error(self, run, write):
        path = write("s.tre", SINGLE)
        result = run("remove", "-n", "A", "B", "C", "-i", str(path))
        assert result.returncode == 1
        assert "would leave tree 1 with no tips" in result.stderr
        assert "Traceback" not in result.stderr

    def test_leaving_a_single_tip_warns(self, run, write):
        path = write("s.tre", SINGLE)
        result = run("remove", "-n", "A", "B", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "leaves 1 tree(s) with a single tip (tree 1)" in result.stderr
        assert "C:2.0;" in result.stdout

    def test_single_tip_warning_is_one_line_for_many_trees(self, run, write):
        # every tree loses A and B, leaving one tip each: one summary line, not one per tree
        trees = "".join("((A:1,B:1):1,C%d:2);\n" % i for i in range(30))
        path = write("many.tre", trees)
        result = run("remove", "-n", "A", "B", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stderr.count("single tip") == 1
        assert "leaves 30 tree(s)" in result.stderr


# ---------------------------------------------------------------- taxa file

class TestTaxaFile:

    def test_blank_lines_in_taxa_file_are_ignored(self, run, write):
        # four tips, so removing two still leaves two and no single-tip warning
        path = write("four.tre", "((A:1,B:1):1,(C:1,D:1):1);\n")
        taxa = write("taxa.txt", "A\n\n   \nB\n")
        result = run("remove", "-f", str(taxa), "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stderr == ""
        assert "C:1.0" in result.stdout and "D:1.0" in result.stdout
        assert "A:" not in result.stdout and "B:" not in result.stdout

    def test_missing_taxa_file_is_an_error(self, run, write):
        path = write("s.tre", SINGLE)
        result = run("remove", "-f", "no_such_taxa.txt", "-i", str(path))
        assert result.returncode == 1
        assert "taxa file not found" in result.stderr
        assert "Traceback" not in result.stderr

    def test_empty_taxa_file_is_an_error(self, run, write):
        path = write("s.tre", SINGLE)
        taxa = write("empty_taxa.txt", "\n\n")
        result = run("remove", "-f", str(taxa), "-i", str(path))
        assert result.returncode == 1
        assert "no taxa given to remove" in result.stderr


# ---------------------------------------------------------------- label quoting in output

class TestQuoting:

    def test_special_labels_keep_their_quotes(self, run, write):
        path = write("q.tre", QUOTED)
        result = run("remove", "-n", "X", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "'A(1)'" in result.stdout
        assert "'B,2'" in result.stdout
        assert "'C:3'" in result.stdout

    def test_plain_labels_stay_unquoted(self, run, write):
        path = write("q.tre", QUOTED)
        result = run("remove", "-n", "X", "-i", str(path))
        assert "Badula_crassa_Y24864:1.0" in result.stdout
        assert "'Badula_crassa_Y24864'" not in result.stdout

    def test_output_can_be_read_back(self, run, write, tmp_path):
        # the output of remove is valid input for summary (round trip)
        path = write("q.tre", QUOTED)
        out = tmp_path / "out.tre"
        result = run("remove", "-n", "X", "-o", str(out), "-i", str(path))
        assert result.returncode == 0, result.stderr
        summary = run("summary", "-i", str(out))
        assert summary.returncode == 0, summary.stderr
        assert "out.tre\t1\t1\t" in summary.stdout


# ---------------------------------------------------------------- log file

class TestLogFile:

    def test_no_log_file_without_option(self, run, write, tmp_path):
        path = write("s.tre", SINGLE)
        run("remove", "-n", "A", "-i", str(path))
        assert not any(p.suffix == ".log" for p in tmp_path.iterdir())

    def test_log_is_appended_across_runs(self, run, write, tmp_path):
        path = write("m.tre", MULTI)
        log = tmp_path / "run.log"
        run("remove", "-n", "D", "-l", str(log), "-i", str(path))
        run("remove", "-n", "C", "-l", str(log), "-i", str(path))
        runs = [line for line in log.read_text().splitlines() if "run: treekit" in line]
        assert len(runs) == 2

    def test_log_records_the_command(self, run, write, tmp_path):
        path = write("s.tre", SINGLE)
        log = tmp_path / "run.log"
        run("remove", "-n", "A", "-l", str(log), "-i", str(path))
        assert "run: treekit remove -n A -l" in log.read_text()

    def test_unwritable_log_path_is_an_error(self, run, write):
        path = write("s.tre", SINGLE)
        result = run("remove", "-n", "A", "-l", "/no/such/dir/x.log", "-i", str(path))
        assert result.returncode == 1
        assert "cannot write log file" in result.stderr
        assert "Traceback" not in result.stderr
