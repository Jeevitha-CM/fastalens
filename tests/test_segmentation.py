import pandas as pd

from segmentation import classify, count_groups, selection_mask


def small_table():
    return pd.DataFrame({
        "Source file": ["a.fasta"] * 4,
        "ID": ["TRINITY_DN1_c0_g1_i1", "TRINITY_DN1_c0_g1_i2", "TRINITY_DN2_c0_g1_i1", "TRINITY_DN3_c0_g1_i1"],
        "Gene": ["TRINITY_DN1_c0_g1", "TRINITY_DN1_c0_g1", "TRINITY_DN2_c0_g1", "TRINITY_DN3_c0_g1"],
        "Format": ["Trinity"] * 4,
        "Type": ["DNA/RNA"] * 4,
        "Length": [150, 199, 200, 5000],
        "GC%": [30.0, 35.0, 50.0, 70.0],
        "Flags": ["short", "short, no ORF", "", ""],
    })


def test_length_class_boundaries():
    # lower bound included, upper bound excluded: 199 -> "< 200", 200 -> "200–499"
    classes = classify(small_table())["Length classes"].astype(str).tolist()
    assert classes == ["< 200", "< 200", "200–499", "≥ 5,000"]


def test_gc_classes():
    classes = classify(small_table())["GC content classes"].astype(str).tolist()
    assert classes == ["AT-rich (< 35%)", "35–44.9%", "45–54.9%", "GC-rich (≥ 65%)"]


def test_isoforms_per_gene():
    classes = classify(small_table())["Isoforms per gene (Trinity)"].astype(str).tolist()
    assert classes == ["2 isoforms", "2 isoforms", "1 isoform", "1 isoform"]


def test_quality_flag_counts_allow_several_flags():
    table = small_table()
    counts = count_groups(classify(table), table["Flags"], "Quality flags").set_index("Group")["Sequences"]
    assert counts["short"] == 2
    assert counts["Clean"] == 2
    assert counts["no ORF"] == 1


def test_click_selection_mask():
    table = small_table()
    mask = selection_mask(classify(table), table["Flags"], "Length classes", ["< 200"])
    assert mask.tolist() == [True, True, False, False]


def test_custom_classes_equal_width_and_range():
    from segmentation import custom_classes
    values = pd.Series([0, 10, 25, 60, 100, 150])
    classes = custom_classes(values, 0, 100, 4, whole_numbers=True)
    # 100 (the maximum of the range) is kept in the last class; 150 is outside the range
    assert classes.iloc[:5].astype(str).tolist() == ["0–25", "0–25", "25–50", "50–75", "75–100"]
    assert pd.isna(classes.iloc[5])