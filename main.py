import streamlit as st
from pathlib import Path
from src.runopenswath import *
from src.openswathresults import *
from src.eic import *
from src.ms2 import *
from src.librarygeneration import *
from src.common import *

st.set_page_config(layout="wide")
st.session_state.mzML_options = [
    p.name for p in Path("mzML-files").glob("*.mzML")]
st.session_state.library_options = [
    p.name for p in Path("assay-libraries").glob("*.tsv")
]
st.session_state.window_options = [
    p.name for p in Path("SWATH-windows").glob("*.tsv")]
st.session_state.mass_list_options = [
    p.name for p in Path("precursor-lists").glob("*.tsv")
]

st.title("OpenSWATH Metabolomics")
t1, t2, t3, t4, t5 = st.tabs(
    [
        "OpenSWATHWorkflow",
        "OpenSWATH Results",
        "Extracted Ion Chromatogram",
        "View MS2 Spectra",
        "Library Generation",
    ]
)

with t1:
    st.multiselect(
        label="mzML files with SWATH data",
        options=st.session_state.mzML_options,
        default=[],
        key="openswath_mzML",
    )
    st.number_input(
        "RT extraction window in seconds",
        2.0,
        60.0,
        10.0,
        1.0,
        key="openswath_rt_window",
        help="Only extract RT around this value (-1 means extract over the whole range, a value of 6.0 means to extract around +/- 3.0 s of the expectedelution). (default: '10.0')",
    )
    st.selectbox(
        label="assay library file",
        options=st.session_state.library_options,
        key="openswath_library",
    )
    st.selectbox(
        label="SWATH window file",
        options=st.session_state.window_options,
        key="openswath_windows",
    )
    _, c1, _ = st.columns(3)
    if c1.button(label="Run OpenSWATH Workflow", type="primary"):
        run_openswath(
            [str(Path("mzML-files", f))
             for f in st.session_state.openswath_mzML],
            str(st.session_state.openswath_rt_window),
            str(Path("assay-libraries", st.session_state.openswath_library)),
            str(Path("SWATH-windows", st.session_state.openswath_windows)),
            "results",
        )

with t2:
    if any(Path("results").glob("*.tsv")):
        c1, c2 = st.columns(2)
        title = c1.text_input(label="custom plot title", value="")

        c1, c2 = st.columns(2)
        df = pd.DataFrame({"filename": [f.name for f in Path("results").glob("*.tsv")]})
        df["time changed"] = [Path("results", f).stat().st_mtime for f in df["filename"]]
        df = df.sort_values("time changed", ascending=False)
        file1 = c1.selectbox("select result file", df["filename"])
        remaining = df["filename"].tolist()
        remaining.remove(file1)
        remaining.insert(0, "none")
        file2 = c2.selectbox("select result file for visual comparison", remaining)
        files = [str(Path("results", file1))]
        if file2 != "none":
            files.append(str(Path("results", file2)))

        fig = plot_openswath_results(files, title)
        show_fig(fig, "openswath-results")
        st.markdown(file1)
        df = pd.read_csv(Path("results", file1), sep="\t")
        show_table(df, "openswath-results")
    else:
        st.warning("No results to show.")

with t3:
    c1, c2 = st.columns(2)
    file = c1.selectbox(label="mzML file", options=st.session_state.mzML_options)
    library = c2.selectbox("assay library", options=st.session_state.library_options)
    c1, c2, c3 = st.columns(3)
    eic_noise = c1.number_input("noise threshhold", 0, 10000, 1000, 1000)
    eic_ppm = c2.number_input("tolerance (ppm)", 1, 50, 10, 1)
    eic_rt_window = c3.number_input("RT window", 1, 300, 5)

    if "eic_df" not in st.session_state:
        st.session_state.eic_df = pd.DataFrame()
    
    _, c2, _ = st.columns(3)
    if c2.button("Extract Ion Chromatograms", type="primary"):
        st.session_state.eic_df = get_extracted_ion_chromatogram(str(Path("mzML-files", file)),
                                            str(Path("assay-libraries", library)),
                                            eic_noise,
                                            eic_rt_window,
                                            eic_ppm)

    if not st.session_state.eic_df.empty:
        fig = px.bar(st.session_state.eic_df["area"])
        fig.update_layout(showlegend=False, xaxis_title="", yaxis_title="intensity")
        show_fig(fig, "eic-area-plot")
        metabolite = st.selectbox("metabolite", st.session_state.eic_df["area"].sort_values(ascending=False).index)
        fig = px.line(x=st.session_state.eic_df.loc[metabolite, "times"], y=st.session_state.eic_df.loc[metabolite, "intensities"])
        fig.update_layout(showlegend=False, xaxis_title="retention time (s)", yaxis_title="counts per second (cps)", title=metabolite)
        show_fig(fig, metabolite)
        show_table(st.session_state.eic_df[["mz", "RT", "area"]], "eic-areas")

with t4:
    df = pd.DataFrame()
    st.selectbox(
        label="mzML file",
        options=st.session_state.mzML_options,
        key="ms2_file",
    )
    if st.session_state.ms2_file:
        df = get_ms2_df(str(Path("mzML-files", st.session_state.ms2_file)))
    if not df.empty:
        st.selectbox(
            "select MS2 spectrum",
            options=[
                f"{i} @{round(mz, 4)} m/z @{round(rt,2)} s"
                for i, mz, rt in zip(df.index, df["precursormz"], df["RT"])
            ],
            key="ms2_spec",
        )
        fig = get_ms2_spec_plot(df, st.session_state.ms2_spec)
        show_fig(fig, st.session_state.ms2_spec)

    else:
        st.warning("No MS2 spectra in data!")

with t5:
    mzML_file = st.selectbox(
        label="mzML file with DDA data for library generation",
        options=st.session_state.mzML_options
    )
    precursor_file = st.selectbox(
        "precursor mass list file",
        options=st.session_state.mass_list_options
    )
    c1, c2 = st.columns(2)
    tolerance_ppm = c1.number_input("mass error in ppm", 1, 50, 10)
    top_n = c2.number_input("take top n MS2 peaks as transitions", 1, 10, 4)
    collision_energy = c1.number_input("collision energy", 10, 50, 10, 5)
    v_space(1, c2)
    exclude_precursor_mass = c2.checkbox("exclude precursor mass", True)
    _, c, _ = st.columns(3)
    df = pd.DataFrame()
    if c.button("Generate Library", type="primary"):
        st.session_state.genlib_filename = f"{mzML_file[:-5]}_{precursor_file[:-4]}_top{top_n}_{exclude_precursor_mass}_ppm{tolerance_ppm}_ce{collision_energy}"
        df = pd.read_csv(str(Path("precursor-lists", precursor_file)), sep="\t", names=["name", "mz", "sum formula"])
        df = generate_library(
                str(Path("precursor-lists", precursor_file)),
                str(Path("mzML-files", mzML_file)),
                top_n,
                exclude_precursor_mass,
                tolerance_ppm,
                collision_energy
        )
        df.to_csv(Path("assay-libraries", st.session_state.genlib_filename+".tsv"), sep="\t")
    if "genlib_filename" in st.session_state:
        df = pd.read_csv(Path("assay-libraries", st.session_state.genlib_filename+".tsv"), sep="\t")
        if not df.empty:
            show_table(df, st.session_state.genlib_filename)