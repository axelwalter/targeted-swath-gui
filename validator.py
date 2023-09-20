import streamlit as st
from pathlib import Path
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import subprocess
import shutil
from src.openswathresults import plot_openswath_results
from src.eic import get_extracted_ion_chromatogram
from src.common import show_fig, show_table

st.set_page_config(layout="wide")

mzML_files = st.multiselect("mzML files", [p.stem for p in Path("mzML-files").glob("*.mzML")])
assay_library = str(st.selectbox("assay library", Path("assay-libraries").glob("*.tsv"), 4))
swath_window = str(st.selectbox("swath windows", Path("SWATH-windows").glob("*.tsv"), 2))
create_eics = st.checkbox("additionally quantify with extracted ion chromatogram area", True)
additional = st.text_input("additional commands", "-rt_extraction_window 10 -Scoring:TransitionGroupPicker:min_peak_width 2 -Scoring:TransitionGroupPicker:use_precursors -mz_extraction_window_ms1 10 -mz_extraction_window 10")


_, c2, _ = st.columns(3)
if c2.button("Run OpenSwathWorkflow", type="primary"):
    with st.status("Running...", expanded=True) as status:
        out_dir = Path("validator-results")
        if out_dir.exists():
            shutil.rmtree(out_dir)
            out_dir.mkdir()
        
        dfs = []
        for mzML_file in mzML_files:
            st.markdown(f"**Processing file: {mzML_file}...**")
            mzML_file = str(Path("mzML-files", mzML_file+".mzML"))
            command = ["OpenSwathWorkflow", "-in", mzML_file, "-tr", assay_library,
                        "-out_tsv", str(Path(out_dir, Path(mzML_file).stem+".tsv")),
                        "-swath_windows_file", swath_window] + additional.split(" ") + ["-force"]
            st.write("Executing OpenSwathWorkflow command:")
            st.write(" ".join(command))
            subprocess.run(command)

            result_file_path = Path(out_dir, Path(mzML_file).stem+".tsv")
            if not result_file_path.exists():
                st.warning(f"Results empty for {mzML_file}")
                continue

            df = pd.read_csv(result_file_path, sep="\t")
            df["CompoundName"] = df.loc[:, "transition_group_id"].apply(lambda x: "_".join(x.split("_")[:-1]))
            df = df.drop(columns=["transition_group_id"])
            df = df.groupby("CompoundName")[["Intensity"]].mean().sort_values(by="Intensity")
            df = df.rename(columns={"Intensity": Path(mzML_file).stem})
            if not df.empty:
                dfs.append(df)
                if create_eics:
                    eic = get_extracted_ion_chromatogram(mzML_file, assay_library, 1000, 10, 10, df.index.tolist())
                    eic.to_pickle(Path("validator-results", f"{Path(mzML_file).stem}_eic.pkl"))
                    eic_intys = eic[["area"]].rename(columns={"area": f"{Path(mzML_file).stem} EIC"})
                    eic_intys.index.name = "CompoundName"
                    dfs.append(eic_intys)
                    st.write("✅ Done exracting chromatograms.")

        if dfs:
            df = pd.concat(dfs, axis=1)
            df = df.sort_index()
            df.to_csv(Path("validator-results", "summary.tsv"), sep="\t")
            st.write("✅ Done combining results.")
            show_table(df, "combined-intensities")
            status.update(label="✅ Complete!", state="complete", expanded=False)
        else:
            status.update(label="No results with selected settings.", state="error")

path = Path("validator-results", "summary.tsv")
if path.exists():
    df = pd.read_csv(path, sep="\t", index_col="CompoundName")
else:
    df = pd.DataFrame()

if not df.empty:
    fig = px.bar(df, barmode="group")
    show_fig(fig, "summary-fig")


dfs = [pd.read_pickle(f) for f in Path("validator-results").glob("*.pkl")]
if dfs:
    # Selectbox to choose a metabolite
    options = []
    for df in dfs:
        options += df.index.tolist()
    options = sorted(set(options))
    selected_metabolite = st.selectbox('Select Metabolite', options)
    filtered_dataframes = [df[df.index == selected_metabolite] for df in dfs]

    fig = go.Figure()
    for name, df in zip([f.stem[:-4] for f in Path("validator-results").glob("*.pkl")], filtered_dataframes):
        fig.add_trace(
            go.Scatter(
                x=df["times"][0],
                y=df["intensities"][0],
                name=name,
            )
        )
    fig.update_layout(
        title=selected_metabolite,
        xaxis_title=f"time (min)",
        yaxis_title="intensity (counts per second)",
        legend_title="sample",
        plot_bgcolor="rgb(255,255,255)",
    )
    show_fig(fig, "eics")
    # Create a list of DataFrames containing only the selected metabolite
    # filtered_dataframes = [df[df.index == selected_metabolite] for df in dfs]

    # print(filtered_dataframes)
# single metabolite EIC plots
# if not eic.empty:
#     metabolite = st.selectbox("metabolite", eic["area"].sort_values(ascending=False).index)
#     fig = px.line(x=eic.loc[metabolite, "times"], y=eic.loc[metabolite, "intensities"])
#     fig.update_layout(showlegend=False, xaxis_title="retention time (s)", yaxis_title="counts per second (cps)", title=metabolite)
#     show_fig(fig, metabolite)

