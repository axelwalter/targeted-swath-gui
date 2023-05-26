import streamlit as st
import pandas as pd
import plotly.express as px
from pyopenms import *


def get_extracted_ion_chromatogram(file, mz, tolerance_ppm):
    # create an empty dataframe to collect chromatogram data in
    df = pd.DataFrame()
    # load mzML file into exp
    exp = MSExperiment()
    MzMLFile().load(str(file), exp)

    # get BPC and time always for each file
    time = []
    intensity = []
    for spec in exp:
        if spec.getMSLevel() == 2:
            continue
        _, intensities = spec.get_peaks()
        rt = spec.getRT()
        time.append(rt)
        i = int(max(intensities))
        intensity.append(i)
    df["time"] = time
    intensity = []
    for spec in exp:
        if spec.getMSLevel() == 2:
            continue
        _, intensities = spec.get_peaks()
        delta = float((tolerance_ppm / 1000000) * mz)
        index_highest_peak_within_window = spec.findHighestInWindow(
            mz, delta, delta)
        # get peak intensity at highest peak index
        peak_intensity = 0
        if index_highest_peak_within_window > -1:
            peak_intensity = int(
                spec[index_highest_peak_within_window].getIntensity())
        intensity.append(peak_intensity)
    df["eic"] = intensity

    return df


def get_eic_plot(df):
    return px.line(df, x="time", y="eic")
