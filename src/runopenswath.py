import streamlit as st
import subprocess
from pathlib import Path
import shutil


@st.cache_resource
def run_openswath(mzML_files, rt_window, library, windows, out_dir):
    if Path(out_dir).exists():
        shutil.rmtree(out_dir)
    Path(out_dir).mkdir()
    for file in mzML_files:
        # Set up command for OpenSwathWorkflow
        command = [
            "OpenSwathWorkflow",
            "-in",
            file,
            "-out_tsv",
            str(Path(out_dir, Path(file).stem + ".tsv")),
            "-tr",
            library,
            "-swath_windows_file",
            windows,
            "-force",
            "-rt_extraction_window",
            rt_window,
        ]

        print("Running command:", subprocess.list2cmdline(command))

        result = subprocess.run(command, capture_output=True, text=True)

        print(result.stdout)

    if any(Path(out_dir).iterdir()):
        st.success("OpenSWATH run was successful.")
    else:
        st.error(
            "Something went wrong during OpenSWATH run, check your inputs and terminal output."
        )
