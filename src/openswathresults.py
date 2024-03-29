import streamlit as st
from pathlib import Path
import pandas as pd
import plotly.express as px


@st.cache_data
def get_openswath_table(path):
    return pd.read_csv(path, sep="\t")


# @st.cache_resource
def plot_openswath_results(files, title):
    dfs = []
    for file in files:
        df = pd.read_csv(file, sep="\t").loc[:, ["transition_group_id", "Intensity"]]
        if df.empty:
            continue
        df["name"] = df.loc[:, "transition_group_id"].apply(lambda x: "_".join(x.split("_")[:-1]))
        df = df.drop(columns=["transition_group_id"])
        df = df.groupby("name").mean().sort_values(by="Intensity")
        df = df.rename(columns={"Intensity": Path(file).stem})
        dfs.append(df)

    df = pd.concat(dfs, axis=1)
    fig = px.bar(df, barmode="group")
    fig.update_layout(title=title, xaxis_title="", yaxis_title="intensity")
    return fig
