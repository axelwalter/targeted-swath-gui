import argparse
import subprocess
import shutil
from pathlib import Path

# Initialize parser
parser = argparse.ArgumentParser()

parser.add_argument(
    "-input",
    help="Input mzML files.",
    default=[
        "/home/a/data/1_4_Spectra_perSec_mzML/20230307_AA_01.mzML",
        "/home/a/data/1_4_Spectra_perSec_mzML/20230307_AA_02.mzML",
    ],
)

parser.add_argument(
    "-input_directory",
    help="Process all mzML files in this directory.",
    default="/home/a/data/5spectra_persec_mzML",
)

parser.add_argument(
    "-output_directory",
    help="Output directory for tsv files with identifications.",
    default="results",
)

parser.add_argument(
    "-library",
    help="Transition library.",
    default="library.tsv",
)

parser.add_argument(
    "-swath_windows_file",
    help="Swath windows in tsv format (columns: start mz, stop mz).",
    default="SWATH-windows.tsv",
)

parser.add_argument(
    "-rt_extraction_window",
    help="Only extract RT around this value (-1 means extract over the whole range, a value of 600 means to extract around +/- 300 s of the expected elution).",
    default="600.0",
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
    # Set up command for OpenSwathWorkflow
    command = [
        "OpenSwathWorkflow",
        "-in",
        file,
        "-out_tsv",
        str(Path(args.output_directory, Path(file).stem + ".tsv")),
        "-tr",
        args.library,
        "-swath_windows_file",
        args.swath_windows_file,
        "-force",
        "-rt_extraction_window",
        args.rt_extraction_window,
    ]

    print("Running command:", subprocess.list2cmdline(command))

    result = subprocess.run(command, capture_output=True, text=True)

    print(result.stdout)

# RT:
# rt_extraction_window: Only extract RT around this value (-1 means extract over the whole range, a value of 600 means to extract around +/- 300 s of the expected elution). (default: '600.0')
