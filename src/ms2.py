import streamlit as st
import pandas as pd
import plotly.express as px
from pyopenms import *
import numpy as np


@st.cache_data
def get_ms2_df(file):
    exp = MSExperiment()
    MzMLFile().load(file, exp)
    df = exp.get_df()
    df.insert(0, "mslevel", [spec.getMSLevel() for spec in exp])
    df.insert(
        0,
        "precursormz",
        [
            spec.getPrecursors()[0].getMZ() if spec.getPrecursors() else 0
            for spec in exp
        ],
    )
    df = df[df["mslevel"] == 2]
    df = df.reset_index()
    df = df.drop("index", axis=1)
    return df


@st.cache_resource
def get_ms2_spec_plot(df, spec):
    def create_spectra(x, y, zero=0):
        x = np.repeat(x, 3)
        y = np.repeat(y, 3)
        y[::3] = y[2::3] = zero
        return pd.DataFrame({"mz": x, "intensity": y})

    df = create_spectra(
        df.loc[int(spec.split(" ")[0]), "mzarray"],
        df.loc[int(spec.split(" ")[0]), "intarray"],
    )
    fig = px.line(df, x="mz", y="intensity")
    fig.update_layout(
        showlegend=False,
        title_text=spec,
        xaxis_title="m/z",
        yaxis_title="intensity",
        plot_bgcolor="rgb(255,255,255)",
    )
    fig.layout.template = "plotly_white"
    fig.update_yaxes(fixedrange=True)
    return fig
