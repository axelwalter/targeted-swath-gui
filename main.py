import streamlit as st
from pathlib import Path
from src.runopenswath import *
from src.openswathresults import *
from src.eic import *
from src.ms2 import *
from src.librarygeneration import *

st.session_state.mzML_options = [p.name for p in Path("mzML-files").glob("*.mzML")]
st.session_state.library_options = [
    p.name for p in Path("assay-libraries").glob("*.tsv")
]
st.session_state.window_options = [p.name for p in Path("SWATH-windows").glob("*.tsv")]
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
    if c1.button(label="Run OpenSWATH Workflow"):
        run_openswath(
            [str(Path("mzML-files", f)) for f in st.session_state.openswath_mzML],
            str(st.session_state.openswath_rt_window),
            str(Path("assay-libraries", st.session_state.openswath_library)),
            str(Path("SWATH-windows", st.session_state.openswath_windows)),
            "results",
        )

with t2:
    if any(Path("results").glob("*.tsv")):
        st.text_input(label="custom plot title", value="", key="openswath_plot_title")
        st.markdown("#")
        title = f"library: {st.session_state.openswath_library} RT window: {st.session_state.openswath_rt_window}"
        if st.session_state.openswath_plot_title:
            title = st.session_state.openswath_plot_title
        fig = plot_openswath_results(
            [str(f) for f in Path("results").glob("*.tsv")], title
        )
        st.plotly_chart(
            fig,
            config={
                "displaylogo": False,
                "modeBarButtonsToRemove": [
                    "zoom",
                    "pan",
                    "select",
                    "lasso",
                    "zoomin",
                    "autoscale",
                    "zoomout",
                    "resetscale",
                ],
                "toImageButtonOptions": {"filename": title.replace(" ", "_")},
            },
        )
    else:
        st.warning("No results to show.")

with t3:
    st.selectbox(
        label="mzML file",
        options=st.session_state.mzML_options,
        key="eic_file",
    )
    c1, c2 = st.columns(2)
    c1.number_input(
        "m/z for extracted ion chromatogram", 50.0, 2000.0, 200.0, 1.0, key="eic_mz"
    )
    c2.number_input("tolerance (ppm)", 1.0, 50.0, 10.0, 1.0, key="eic_ppm")

    df = pd.DataFrame()

    _, c1, _ = st.columns(3)
    if c1.button("Extract Ion Chromatogram") and st.session_state.eic_file:
        df = get_extracted_ion_chromatogram(
            str(Path("mzML-files", st.session_state.eic_file)),
            st.session_state.eic_mz,
            st.session_state.eic_ppm,
        )

    if not df.empty:
        fig = get_eic_plot(df)
        st.plotly_chart(
            fig,
            config={
                "displaylogo": False,
                "modeBarButtonsToRemove": [
                    "zoom",
                    "pan",
                    "select",
                    "lasso",
                    "zoomin",
                    "autoscale",
                    "zoomout",
                    "resetscale",
                ],
                "toImageButtonOptions": {
                    "filename": f"{Path(st.session_state.eic_file).stem.replace(' ', '_')}_{st.session_state.eic_mz}"
                },
            },
        )

with t4:
    df = pd.DataFrame()
    st.selectbox(
        label="mzML file",
        options=st.session_state.mzML_options,
        key="ms2_file",
    )
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
        st.plotly_chart(
            fig,
            config={
                "displaylogo": False,
                "modeBarButtonsToRemove": [
                    "zoom",
                    "pan",
                    "select",
                    "lasso",
                    "zoomin",
                    "autoscale",
                    "zoomout",
                    "resetscale",
                ],
                "toImageButtonOptions": {"filename": st.session_state.ms2_spec},
            },
        )

    else:
        st.warning("No MS2 spectra in data!")

with t5:
    st.selectbox(
        label="mzML file with DDA data for library generation",
        options=st.session_state.mzML_options,
        key="generate_file",
    )
    st.selectbox(
        "mass list file",
        options=st.session_state.mass_list_options,
        key="generate_mass_list",
    )
    c1, c2 = st.columns(2)
    c1.number_input("mass error in ppm", 1, 50, 10, key="generate_ppm")
    c2.number_input(
        "take top n MS2 peaks as transitions", 1, 10, 4, key="generate_top_n"
    )
    _, c, _ = st.columns(3)
    if c.button("Generate Library"):
        generate_library(
            str(Path("mzML-files", st.session_state.generate_file)),
            str(Path("precursor-lists", st.session_state.generate_mass_list)),
            st.session_state.generate_top_n,
            st.session_state.generate_ppm,
        )
        st.experimental_rerun()
