import streamlit as st
import pandas as pd
import numpy as np
from pyopenms import *

@st.cache_data
def get_extracted_ion_chromatogram(file, library, noise, rt_window, tolerance_ppm, openswath_metabolites=[]):
    # load compound names, mz and RT values from library
    lib = pd.read_csv(library, sep="\t").groupby("CompoundName").mean("PrecursorMz")[["PrecursorMz", "NormalizedRetentionTime"]]
    lib.index = pd.Index([x.replace(",", "") for x in lib.index])

    if openswath_metabolites:
        lib = lib[lib.index.isin(openswath_metabolites)]
    # load mzML file into exp
    exp = MSExperiment()
    MzMLFile().load(str(file), exp)

    def extract_intensities(row):
        intys = []
        for spec in exp:
            if spec.getMSLevel() == 2:
                continue
            if spec.getRT() < row["NormalizedRetentionTime"] - rt_window/2 or spec.getRT() > row["NormalizedRetentionTime"] + rt_window/2:
                peak_intensity = 0
            else:
                delta = float((tolerance_ppm / 1000000) * row["PrecursorMz"])
                index_highest_peak_within_window = spec.findHighestInWindow(
                    row["PrecursorMz"], delta, delta)
                # get peak intensity at highest peak index if above noise
                peak_intensity = 0
                if index_highest_peak_within_window > -1:
                    peak_intensity = int(
                        spec[index_highest_peak_within_window].getIntensity())
                    if peak_intensity < noise:
                        peak_intensity = 0
            intys.append(peak_intensity)
        return np.array(intys)
    
    times = []
    for spec in exp:
        if spec.getMSLevel() == 1:
            times.append(spec.getRT())
    
    lib["intensities"] = lib.apply(extract_intensities, axis=1)
    lib["times"] = [np.array(times) for _ in range(lib.shape[0])]
    lib["area"] = lib["intensities"].apply(lambda x: int(np.trapz(x)))

    lib = lib.rename(columns={"NormalizedRetentionTime": "RT", "PrecursorMz": "mz"})
    lib.index.name = "name"

    return lib.sort_values("area")

