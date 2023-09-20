import streamlit as st
from pathlib import Path
import plotly.express as px
import pandas as pd
import subprocess
import shutil
from src.openswathresults import plot_openswath_results
from src.eic import get_extracted_ion_chromatogram
from src.common import show_fig, show_table

# st.set_page_config(layout="wide")

mzML_file = str(st.selectbox("mzML file", Path("mzML-files").glob("*.mzML")))
assay_library = str(st.selectbox("assay library", Path("assay-libraries").glob("*.tsv"), 4))
swath_window = str(st.selectbox("swath windows", Path("SWATH-windows").glob("*.tsv"), 2))
additional = st.text_input("additional commands", "-rt_extraction_window 10")


_, c2, _ = st.columns(3)
if c2.button("Run OpenSwathWorkflow", type="primary"):
    with st.status("Running...", expanded=True) as status:
        out_dir = Path("validator-results")
        if out_dir.exists():
            shutil.rmtree(out_dir)
            out_dir.mkdir()
        command = ["OpenSwathWorkflow", "-in", mzML_file, "-tr", assay_library,
                    "-out_tsv", str(Path(out_dir, Path(mzML_file).stem+".tsv")),
                    "-swath_windows_file", swath_window] + additional.split(" ") + ["-force"]
        st.write("Executing OpenSwathWorkflow command:")
        st.write(" ".join(command))
        subprocess.run(command)

        df = pd.read_csv(Path(out_dir, Path(mzML_file).stem+".tsv"), sep="\t")
        metabolites_found = list(set([m.split("_run")[0] for m in df["transition_group_id"]]))

        eic = get_extracted_ion_chromatogram(mzML_file, assay_library, 1000, 10, 10, metabolites_found)
        eic.to_pickle(Path("validator-results", "eic.pkl"))
        st.write("✅ Done exracting chromatograms.")
        eic_intys = eic[["area"]].rename(columns={"area": "EIC intensity"})
        eic_intys.index.name = "CompoundName"
        dfs = [eic_intys]

        for file in out_dir.glob("*.tsv"):
            if file.stem != "summary":
                df = pd.read_csv(file, sep="\t").loc[:, ["transition_group_id", "Intensity"]]
                if df.empty:
                    continue
                df["CompoundName"] = df.loc[:, "transition_group_id"].apply(lambda x: "_".join(x.split("_")[:-1]))
                df = df.drop(columns=["transition_group_id"])
                df = df.groupby("CompoundName").mean().sort_values(by="Intensity")
                df = df.rename(columns={"Intensity": Path(file).stem})
                dfs.append(df)

        df = pd.concat(dfs, axis=1)
        df = df.sort_values(by=df.columns[0])
        df.to_csv(Path("validator-results", "summary.tsv"), sep="\t")
        st.write("✅ Done combining XIC and OpenswathResults.")
        show_table(df, "combined-intensities")
        status.update(label="✅ Complete!", state="complete", expanded=False)

else:
    path = Path("validator-results", "summary.tsv")
    if path.exists():
        df = pd.read_csv(path, sep="\t", index_col="CompoundName")
    else:
        df = pd.DataFrame()
    path = Path("validator-results", "eic.pkl")
    if path.exists():
        eic = pd.read_pickle(path)
    else:
        eic = pd.DataFrame()

if not df.empty:
    fig = px.bar(df, barmode="group")
    show_fig(fig, "summary-fig")

# single metabolite EIC plots
if not eic.empty:
    metabolite = st.selectbox("metabolite", eic["area"].sort_values(ascending=False).index)
    fig = px.line(x=eic.loc[metabolite, "times"], y=eic.loc[metabolite, "intensities"])
    fig.update_layout(showlegend=False, xaxis_title="retention time (s)", yaxis_title="counts per second (cps)", title=metabolite)
    show_fig(fig, metabolite)

