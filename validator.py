import streamlit as st
from pathlib import Path
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import numpy as np
import subprocess
import shutil
from src.eic import get_extracted_ion_chromatogram
from src.common import show_fig, show_table
import pyopenms as poms


st.markdown("""
    <style>
        .stMultiSelect [data-baseweb=select] span{
            max-width: 500px;
            font-size: 1rem;
        }
    </style>
    """, unsafe_allow_html=True)

# st.set_page_config(layout="wide")
mix2_001 = [f"{conc}uM_Mix2_Bioblank_pos_001" for conc in ("01", "05", "1", "5", "25")]

mzML_files = st.multiselect("mzML files", [p.stem for p in Path("mzML-files").glob("*.mzML")], mix2_001)
assay_library = str(st.selectbox("assay library", Path("assay-libraries").glob("*.tsv"), 4))
swath_window = str(st.selectbox("swath windows", Path("SWATH-windows").glob("*.tsv"), 2))
additional = st.text_area("additional commands", """-rt_extraction_window 60
-Scoring:TransitionGroupPicker:min_peak_width 2
-Scoring:TransitionGroupPicker:use_precursors
-mz_extraction_window_ms1 50
-mz_extraction_window 50
""", height=200)

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
                        "-out_chrom", str(Path(out_dir, Path(mzML_file).stem+"_chrom.mzML")),
                        "-ms1_isotopes", "0",
                        "-Scoring:TransitionGroupPicker:compute_peak_shape_metrics",
                        "-swath_windows_file", swath_window] + additional.strip().replace("\n", " ").split(" ") + ["-force"]
            subprocess.run(command)

            result_file_path = Path(out_dir, Path(mzML_file).stem+".tsv")
            if not result_file_path.exists():
                st.warning(f"Results empty for {mzML_file}")
                continue

            df = pd.read_csv(result_file_path, sep="\t")
            # df["CompoundName"] = df["peptide_group_label"] # df["transition_group_id"].apply(lambda x: x.split("_")[0])
            # df = df.drop(columns=["transition_group_id"])
            df = df.groupby("peptide_group_label")[["aggr_prec_Peak_Area"]].max() #.sort_values(by="aggr_prec_Peak_Area")
            df["CompoundName"] = [x.split("_")[0] for x in df.index.tolist()]
            df = df.groupby("CompoundName")[["aggr_prec_Peak_Area"]].mean()
            df = df.rename(columns={"aggr_prec_Peak_Area": Path(mzML_file).stem})
            if not df.empty:
                dfs.append(df)
                names = []
                rts = []
                intys = []
                exp = poms.MSExperiment()
                poms.MzMLFile().load(str(Path("validator-results", Path(mzML_file).stem +  "_chrom.mzML")), exp)
                for chrom in exp.getChromatograms():
                    if chrom.getChromatogramType() == 3: # 3 == BASEPEAK_CHROMATOGRAM, 5 = SELECTED_REACTION_MONITORING_CHROMATOGRAM
                        name = chrom.getPrecursor().getMetaValue("peptide_sequence")
                        if name not in names:
                            names.append(name) 
                            rt, inty = chrom.get_peaks()
                            rts.append(rt)
                            intys.append(inty)
                chroms = pd.DataFrame({"name": names, "times": rts, "intensities": intys})
                chroms = chroms.set_index("name")
                chroms = chroms[chroms.index.isin(df.index)]

                chroms.to_pickle(Path("validator-results", f"{Path(mzML_file).stem}_chrom.pkl"))
                
                eic = get_extracted_ion_chromatogram(mzML_file, assay_library, 100, 60, 25)
                eic.to_pickle(Path("validator-results", f"{Path(mzML_file).stem}_eic.pkl"))
                eic_intys = eic[["area"]].rename(columns={"area": f"{Path(mzML_file).stem} EIC"})
                eic_intys.index.name = "CompoundName"
                dfs.append(eic_intys)
            else:
                st.warning(f"No results for file {mzML_file}")

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
    df = df.fillna(0)
    show_eic_auc = st.checkbox("show EIC intensities", True)
    if not show_eic_auc:
        df = df[[col for col in df if "EIC" not in col]]
    # Remove rows where all values are 0 or None
    # df = df.dropna(how="all")
    df = df[(df != 0).any(axis=1)]
    fig = px.bar(df, barmode="group")
    show_fig(fig, "summary-fig")

if any(Path("validator-results").glob("*_chrom.pkl")):
    c1, c2 = st.columns(2)
    file = c1.selectbox("show chromatograms for file", [f.stem[:-6] for f in Path("validator-results").glob("*_chrom.mzML")])
    df_openswath = pd.read_pickle(Path("validator-results", file+"_chrom.pkl"))
    metabolite_options = df_openswath.index
    eic_path = Path("validator-results", file+"_eic.pkl")
    if eic_path.exists():
        df_eic = pd.read_pickle(eic_path)
        metabolite_options = df_eic.index
    metabolite = c2.selectbox("metabolite", sorted(metabolite_options))

    # Add OpenSWATH chromatogram
    df = df_openswath[df_openswath.index == metabolite]
    fig = go.Figure()
    if not df.empty:
        fig.add_trace(
            go.Scatter(
                x=df["times"][0],
                y=df["intensities"][0],
                name="OpenSWATH",
                line={"color": "#636efa", "width": 4}
            )
        )
    else:
        st.warning(f"No OpenSWATH result for {metabolite}")
    # Add EIC chromatogram
    df = df_eic[df_eic.index == metabolite]
    if not df.empty:
        fig.add_trace(
            go.Scatter(
                x=df["times"][0],
                y=df["intensities"][0],
                name="EIC",
                line={"color": "#ef553b"}
            )
        )
    fig.update_layout(
        title=file + ": " +metabolite,
        xaxis_title=f"time (min)",
        yaxis_title="intensity (counts per second)",
        legend_title="sample",
        plot_bgcolor="rgb(255,255,255)",
    )
    show_fig(fig, "eics")



# dfs = [pd.read_pickle(f) for f in Path("validator-results").glob("*_eic.pkl")]
# if dfs:
#     # Selectbox to choose a metabolite
#     options = []
#     for df in dfs:
#         options += df.index.tolist()
#     options = sorted(set(options))
#     selected_metabolite = st.selectbox('Select Metabolite', options)
#     filtered_dataframes = [df[df.index == selected_metabolite] for df in dfs]

#     fig = go.Figure()
#     for name, df in zip([f.stem[:-4] for f in Path("validator-results").glob("*_eic.pkl")], filtered_dataframes):
#         fig.add_trace(
#             go.Scatter(
#                 x=df["times"][0],
#                 y=df["intensities"][0],
#                 name=name,
#             )
#         )
#     fig.update_layout(
#         title=selected_metabolite,
#         xaxis_title=f"time (min)",
#         yaxis_title="intensity (counts per second)",
#         legend_title="sample",
#         plot_bgcolor="rgb(255,255,255)",
#     )
#     show_fig(fig, "eics")