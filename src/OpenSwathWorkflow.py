import streamlit as st
import argparse
import subprocess
import shutil
from pathlib import Path
import pandas as pd

# Initialize parser
parser = argparse.ArgumentParser()

parser.add_argument(
    "-input",
    help="Input mzML files.",
    default=[],
)

parser.add_argument(
    "-input_directory",
    help="Process all mzML files in this directory.",
    default="",
)

parser.add_argument(
    "-output_directory",
    help="Output directory for tsv files with identifications.",
    default="results",
)

parser.add_argument(
    "-library",
    help="Transition library.",
    default="",
)

parser.add_argument(
    "-swath_windows_file",
    help="Swath windows in tsv format (columns: start mz, stop mz).",
    default="",
)

parser.add_argument(
    "-rt_extraction_window",
    help="Only extract RT around this value (-1 means extract over the whole range, a value of 600 means to extract around +/- 300 s of the expected elution).",
    default="10.0",
)

# Read arguments from command line
args = parser.parse_args()

# get mzML files as individual files
mzML_files = args.input

# OR read all mzML from a given directory
if args.input_directory:
    mzML_files = [str(path) for path in Path(args.input_directory).glob("*.mzML")]

# override results directory
out = Path(args.output_directory)
if out.is_dir():
    shutil.rmtree(out)
out.mkdir()

for file in mzML_files:
    out = str(Path(args.output_directory, Path(file).stem + ".tsv"))
    # Set up command for OpenSwathWorkflow
    command = [
        "OpenSwathWorkflow",
        "-in",
        file,
        "-out_tsv",
        out,
        "-tr",
        args.library,
        "-swath_windows_file",
        args.swath_windows_file,
        "-force",
        "-rt_extraction_window",
        args.rt_extraction_window,
        # "--use_ms1_traces"
        "-out_chrom",
        "results/chrom.mzML"
    ]

    print("Running command:", subprocess.list2cmdline(command))

    result = subprocess.run(command, capture_output=True, text=True)

    print(result.stdout)

    df = pd.read_csv(out, sep="\t")
    if df.empty:
        st.warning(f"Empty output for {file} generated!")
        Path(out).unlink()


# RT:
# rt_extraction_window: Only extract RT around this value (-1 means extract over the whole range, a value of 600 means to extract around +/- 300 s of the expected elution). (default: '600.0')
