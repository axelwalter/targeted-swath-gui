import streamlit as st
import subprocess
from pathlib import Path
import pandas as pd


def run_openswath(mzML_files, rt_window, library, windows, out_dir):
    Path(out_dir).mkdir(exist_ok=True)
    for file in mzML_files:
        out_file = Path(out_dir, f"{Path(file).stem}_{Path(library).stem}_{rt_window}s_nones.tsv")
        # Set up command for OpenSwathWorkflow
        command = [
            "OpenSwathWorkflow",
            "-in",
            file,
            "-out_tsv",
            str(out_file),
            "-tr",
            library,
            "-swath_windows_file",
            windows,
            "-force",
            "-rt_extraction_window",
            rt_window,
            # "-Scoring:TransitionGroupPicker:min_peak_width",
            # str(30.0),
        ]

        print("Running command:", subprocess.list2cmdline(command))

        result = subprocess.run(command, capture_output=True, text=True)

        print(result.stdout)

        if out_file.exists():
            st.success("OpenSWATH run was successful.")
            df = pd.read_csv(out_file, sep="\t")
            if df.empty:
                st.warning("Results are empty, no metabolites detected. This run will not be shown in results.")
        else:
            st.error(
                "Something went wrong during OpenSWATH run, check your inputs and terminal output."
            )
