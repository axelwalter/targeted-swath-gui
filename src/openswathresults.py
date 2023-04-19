import streamlit as st
from pathlib import Path
import pandas as pd
import plotly.express as px


@st.cache_data
def get_openswath_table(path):
    return pd.read_csv(path, sep="\t")


@st.cache_resource
def plot_openswath_results(files, title=""):
    dfs = []
    for file in files:
        df = pd.read_csv(file, sep="\t").loc[:, ["transition_group_id", "Intensity"]]
        df = df.rename(
            columns={"transition_group_id": "group_id", "Intensity": "intensity"}
        )
        df["name"] = df.loc[:, "group_id"].apply(lambda x: x.split("_")[0])
        df = df.groupby("name").mean().sort_values(by="intensity", ascending=False)
        df = df.rename(columns={"intensity": Path(file).stem})
        dfs.append(df)

    df = pd.concat(dfs, axis=1)
    fig = px.bar(df, barmode="group")
    fig.update_layout(title=title, xaxis_title="", yaxis_title="intensity")
    return fig
