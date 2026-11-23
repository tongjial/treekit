NEXUS_NO_TREES = "#NEXUS\nBegin taxa;\nDimensions ntax=2;\nTaxlabels A B;\nEnd;\n"
MULTI = "((A:1,B:1):1,C:2);\n((A:1,C:1):1,D:2);\n((B:1,C:1):1,D:2);\n"


class TestTaxa:

    def test_all_taxa_are_listed_sorted(self, run, write):
        path = write("m.tre", MULTI)
        result = run("taxa", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stdout.splitlines() == ["A", "B", "C", "D"]

    def test_common_taxa_only(self, run, write):
        # only C is in all three trees
        path = write("m.tre", MULTI)
        result = run("taxa", "-k", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stdout.splitlines() == ["C"]

    def test_taxa_from_several_files_are_merged(self, run, write):
        write("a.tre", "((A:1,B:1):1,C:2);\n")
        write("b.tre", "((X:1,Y:1):1,C:2);\n")
        result = run("taxa", "-i", "a.tre", "b.tre")
        assert result.stdout.splitlines() == ["A", "B", "C", "X", "Y"]

    def test_common_taxa_across_several_files(self, run, write):
        write("a.tre", "((A:1,B:1):1,C:2);\n")
        write("b.tre", "((X:1,Y:1):1,C:2);\n")
        result = run("taxa", "-k", "-i", "a.tre", "b.tre")
        assert result.stdout.splitlines() == ["C"]

    def test_no_common_taxa_gives_an_empty_list(self, run, write):
        write("a.tre", "((A:1,B:1):1,C:2);\n")
        write("b.tre", "((X:1,Y:1):1,Z:2);\n")
        result = run("taxa", "-k", "-i", "a.tre", "b.tre")
        assert result.returncode == 0, result.stderr
        assert result.stdout.strip() == ""

    def test_underscores_are_preserved(self, run, write):
        path = write("u.tre", "((Genus_species_1:1,B:1):1,C:2);\n")
        result = run("taxa", "-i", str(path))
        assert "Genus_species_1" in result.stdout.splitlines()

    def test_quoted_labels_with_spaces(self, run, write):
        path = write("q.tre", "(('A b':1,B:1):1,C:2);\n")
        result = run("taxa", "-i", str(path))
        assert "A b" in result.stdout.splitlines()

    def test_out_file(self, run, write, tmp_path):
        path = write("m.tre", MULTI)
        out = tmp_path / "taxa.txt"
        result = run("taxa", "-o", str(out), "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stdout == ""
        assert out.read_text() == "A\nB\nC\nD\n"

    def test_taxa_file_can_feed_remove(self, run, write, tmp_path):
        # the taxa command writes one name per line, which is what remove -f reads
        path = write("m.tre", MULTI)
        out = tmp_path / "taxa.txt"
        run("taxa", "-o", str(out), "-i", str(path))
        write("four.tre", "((A:1,B:1):1,(C:1,D:1):1);\n")
        result = run("remove", "-f", str(out), "-i", "four.tre")
        # every taxon is requested, so no tips would be left
        assert result.returncode == 1
        assert "no tips" in result.stderr

    def test_tips_without_a_label_are_skipped_in_common_taxa(self, run, write):
        path = write("nolabel.tre", "((:1,B:1):1,C:2);\n")
        result = run("taxa", "-k", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stdout.splitlines() == ["B", "C"]
        assert "Traceback" not in result.stderr

    def test_tips_without_a_label_are_skipped_in_all_taxa(self, run, write):
        path = write("nolabel.tre", "((:1,B:1):1,C:2);\n")
        result = run("taxa", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stdout.splitlines() == ["B", "C"]

    def test_no_trees_with_common_taxa_is_an_error(self, run, write):
        path = write("nt.nex", NEXUS_NO_TREES)
        result = run("taxa", "-k", "-i", str(path))
        assert result.returncode == 1
        assert "no trees found" in result.stderr

    def test_mixed_newick_and_nexus_input(self, run, write):
        write("a.tre", "((A:1,B:1):1,C:2);\n")
        write("n.nex", "#NEXUS\nBegin taxa;\nDimensions ntax=2;\nTaxlabels C D;\nEnd;\n"
                       "Begin trees;\nTree t = (C:1,D:1);\nEnd;\n")
        result = run("taxa", "-i", "a.tre", "n.nex")
        assert result.returncode == 0, result.stderr
        assert result.stdout.splitlines() == ["A", "B", "C", "D"]


class TestRemoveWithUnlabelledTips:

    def test_unlabelled_tips_are_rejected(self, run, write):
        # dendropy would delete a tip without a label together with the removed taxa,
        # so the run stops instead of losing data silently
        path = write("ul.tre", "((:1,B:1):1,C:2);\n")
        result = run("remove", "-n", "C", "-i", str(path))
        assert result.returncode == 1
        assert "tree 1 has a tip without a label" in result.stderr
        assert result.stdout == ""
        assert "Traceback" not in result.stderr

    def test_the_tree_number_is_reported(self, run, write):
        path = write("ul.tre", "((A:1,B:1):1,C:2);\n((:1,B:1):1,C:2);\n")
        result = run("remove", "-n", "C", "-i", str(path))
        assert result.returncode == 1
        assert "tree 2 has a tip without a label" in result.stderr


class TestSeveralInputFiles:

    NEXUS_CD = ("#NEXUS\nBegin taxa;\nDimensions ntax=2;\nTaxlabels C D;\nEnd;\n"
                "Begin trees;\nTree t = (C:1,D:1);\nEnd;\n")

    def test_remove_over_newick_and_nexus_files(self, run, write):
        # the NEXUS file declares two taxa, fewer than the newick file already holds
        write("a.tre", "((A:1,B:1):1,C:2);\n")
        write("n.nex", self.NEXUS_CD)
        result = run("remove", "-n", "C", "-i", "a.tre", "n.nex")
        assert result.returncode == 0, result.stderr
        lines = [line for line in result.stdout.splitlines() if line.strip()]
        assert lines == ["(A:1.0,B:1.0):1.0;", "D:1.0;"]
        assert "single tip" in result.stderr

    def test_nexus_output_covers_the_taxa_of_all_files(self, run, write):
        write("a.tre", "((A:1,B:1):1,C:2);\n")
        write("n.nex", self.NEXUS_CD)
        result = run("remove", "-n", "A", "-u", "nexus", "-i", "a.tre", "n.nex")
        assert result.returncode == 0, result.stderr
        # B, C and D are used by the trees; the removed A is not declared any more
        assert "NTAX=3" in result.stdout
        taxlabels = result.stdout.split("TAXLABELS")[1].split(";")[0].split()
        assert sorted(taxlabels) == ["B", "C", "D"]

    def test_two_nexus_files_with_different_taxa(self, run, write):
        write("n1.nex", "#NEXUS\nBegin taxa;\nDimensions ntax=3;\nTaxlabels A B C;\nEnd;\n"
                        "Begin trees;\nTree t = ((A:1,B:1):1,C:2);\nEnd;\n")
        write("n2.nex", self.NEXUS_CD)
        result = run("taxa", "-i", "n1.nex", "n2.nex")
        assert result.returncode == 0, result.stderr
        assert result.stdout.splitlines() == ["A", "B", "C", "D"]


class TestNexusTaxaBlock:

    def test_removed_taxa_are_not_declared(self, run, write):
        path = write("m.tre", "((A:1,B:1):1,(C:1,D:1):1);\n((A:1,C:1):1,(B:1,D:1):1);\n")
        result = run("remove", "-n", "A", "-u", "nexus", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "NTAX=3" in result.stdout
        taxlabels = result.stdout.split("TAXLABELS")[1].split(";")[0].split()
        assert sorted(taxlabels) == ["B", "C", "D"]

    def test_partially_present_taxon_is_removed_from_the_declaration(self, run, write):
        # D is in 2 of 3 trees; after the removal no tree has it
        path = write("m.tre", "((A:1,B:1):1,C:2);\n((A:1,C:1):1,D:2);\n((B:1,C:1):1,D:2);\n")
        result = run("remove", "-n", "D", "-u", "nexus", "-i", str(path))
        assert "NTAX=3" in result.stdout


class TestCaseSensitiveLabels:

    TREE = "((a:1,A:1):1,C:2);\n"

    def test_labels_differing_in_case_are_different_taxa(self, run, write):
        path = write("case.tre", self.TREE)
        result = run("taxa", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stdout.splitlines() == ["A", "C", "a"]

    def test_summary_counts_both(self, run, write):
        path = write("case.tre", self.TREE)
        result = run("summary", "-i", str(path))
        row = result.stdout.splitlines()[1].split("\t")
        assert row[8] == "3"            # num_taxa
        assert row[1] == "1"

    def test_removing_a_lowercase_label_keeps_the_uppercase_one(self, run, write):
        path = write("case.tre", self.TREE)
        result = run("remove", "-n", "a", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stdout.strip() == "(A:2.0,C:2.0);"

    def test_a_label_in_the_wrong_case_is_not_found(self, run, write):
        # the tree has "A" only; "a" is another taxon and must be reported as absent,
        # not silently remove "A"
        path = write("upper.tre", "((A:1,B:1):1,C:2);\n")
        result = run("remove", "-n", "a", "-i", str(path))
        assert result.returncode == 1
        assert "none of the 1 taxon to remove were found" in result.stderr
        assert result.stdout == ""

    def test_wrong_case_among_present_taxa_warns_and_keeps_the_tip(self, run, write):
        path = write("upper.tre", "((A:1,B:1):1,(C:1,D:1):1);\n")
        result = run("remove", "-n", "a", "C", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert "1 taxon was not found in any input tree: a" in result.stderr
        assert "A:1.0" in result.stdout        # still there

    def test_nexus_input_keeps_case(self, run, write):
        path = write("c.nex", "#NEXUS\nBegin taxa;\nDimensions ntax=3;\nTaxlabels a A C;\nEnd;\n"
                              "Begin trees;\nTree t = ((a:1,A:1):1,C:2);\nEnd;\n")
        result = run("taxa", "-i", str(path))
        assert result.returncode == 0, result.stderr
        assert result.stdout.splitlines() == ["A", "C", "a"]
