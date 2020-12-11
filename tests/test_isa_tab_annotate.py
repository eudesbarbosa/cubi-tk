"""Tests for ``cubi_tk.isa_tab``.

"""
import os
import glob
import filecmp

from shutil import copytree

from cubi_tk.isa_tpl import run_cookiecutter, TEMPLATES
from cubi_tk.__main__ import setup_argparse, main


def test_run_isatab_annotate_case1(tmp_path):
    path_input = os.path.join(os.path.dirname(__file__), "data", "isa_tab_annotate")
    path_input = copytree(path_input, os.path.join(tmp_path, "isa_tab_annotate"))
    #fs.add_real_directory(path_input, read_only=False)

    path_input_annotation = os.path.join(os.path.dirname(__file__), "data", "isa_tab_annotate",
                                         "isa_tab_annotation.csv")

    argv = [
        "isa-tab",
        "annotate",
        "--yes",
        os.path.join(path_input, "i_Investigation.txt"),
        path_input_annotation,
    ]

    res = main(argv)
    assert not res

    # add files
    path_test = os.path.join(os.path.dirname(__file__), "data", "isa_tab_annotate_result1")
    #fs.add_real_directory(path_test)

    # tests
    # fs.add_real_directory(output_path)
    files = glob.glob(os.path.join(path_test, "*"))

    match, mismatch, errors = filecmp.cmpfiles(
        path_test,
        path_input,
        (os.path.basename(f) for f in files),
        shallow=False
        )
    print([match, mismatch, errors])
    assert len(mismatch) == 0
    assert len(errors) == 0


def test_run_isatab_annotate_case2(tmp_path):
    path_input = os.path.join(os.path.dirname(__file__), "data", "isa_tab_annotate")
    copytree(path_input, tmp_path)
    path_input = os.path.join(tmp_path, "isa_tab_annotate")
    #fs.add_real_directory(path_input, read_only=False)

    path_input_annotation = os.path.join(os.path.dirname(__file__), "data", "isa_tab_annotate",
                                         "isa_tab_annotation.csv")

    argv = [
        "isa-tab",
        "annotate",
        "--force-update",
        "--yes",
        os.path.join(path_input, "i_Investigation.txt"),
        path_input_annotation,
    ]

    res = main(argv)
    assert not res

    # add files
    path_test = os.path.join(os.path.dirname(__file__), "data", "isa_tab_annotate_result2")
    #fs.add_real_directory(path_test)

    # tests
    # fs.add_real_directory(output_path)
    files = glob.glob(os.path.join(path_test, "*"))

    match, mismatch, errors = filecmp.cmpfiles(
        path_test,
        path_input,
        (os.path.basename(f) for f in files),
        shallow=False
        )
    print([match, mismatch, errors])
    assert len(mismatch) == 0
    assert len(errors) == 0
