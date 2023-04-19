import streamlit as st
from pyopenms import *
from massql import msql_engine
import pandas as pd
from pathlib import Path


def get_transitions(mzML_file, ms1, ms2, prec_mz, top_n, ppm_error):
    # get all MS2 spectra with precursor mass
    df = msql_engine.process_query(
        f"QUERY scaninfo(MS2DATA) WHERE MS2PREC={prec_mz}:TOLERANCEPPM={ppm_error}",
        mzML_file,
        ms1_df=ms1,
        ms2_df=ms2,
    )
    if df.empty:
        return []
    # get scan number of MS2 spectrum with highest TIC intensity (best quality spectrum)
    scan = df.loc[df["i"].idxmax(), "scan"]
    # now query that spectrum for the highest intensity peak
    df = ms2[ms2["scan"] == scan].sort_values("i_norm", ascending=False)
    # get the two highest intensity peaks with a buffer of 2 Da (to not peak peaks very close together)
    # only peaks equal or smaller precursor mass are considered
    df = df[df["mz"] < (prec_mz + 1)]
    intensities = []
    product_mzs = []
    # iteratively pick the top_n peaks and exclude them with the buffer
    for _ in range(top_n):
        i = df["i"].idxmax()
        mz = df.loc[i, "mz"]
        product_mzs.append(mz)
        intensities.append(df.loc[i, "i"])
        df = df[~((df["mz"] >= mz - 2) & (df["mz"] <= mz + 2))]
    # normalize intensities
    intensities = [x / max(intensities) for x in intensities]
    return list(zip(product_mzs, intensities))


def get_RT(ms1, prec_mz, ppm_error):
    # calculate mass error at given ppm value
    mass_error = (ppm_error / 1000000) * prec_mz
    # filter ms1 peaks for mz values within the error range and sort values by intensity in descending order
    df = ms1[
        (ms1["mz"] > (prec_mz - mass_error)) & (ms1["mz"] < (prec_mz + mass_error))
    ].sort_values("i", ascending=False)
    # calculate the RT from the mean of the top peaks in seconds
    rt = round(df["rt"].iloc[:10].mean() * 60)
    return rt


@st.cache_resource
def generate_library(mzML_file, mass_list_file, top_n, ppm_error):
    # load MSExperiment
    exp = MSExperiment()
    MzMLFile().load(mzML_file, exp)

    # get MS1 and MS2 dataframes
    ms1, ms2 = exp.get_massql_df()

    with open(mass_list_file, "r") as f:
        queries = []
        for line in f.readlines():
            l = line.strip().split("\t")
            name, mz = l[0], float(l[1])
            queries.append((name, mz))

    compound_name = []
    precursor_mz = []
    product_mz = []
    library_intensity = []
    normalized_retention_time = []
    transition_group_id = []
    collision_energy = []

    for query in queries:
        transitions = get_transitions(mzML_file, ms1, ms2, query[1], top_n, ppm_error)
        for t in transitions:
            compound_name.append(query[0])
            precursor_mz.append(query[1])
            product_mz.append(t[0])
            library_intensity.append(t[1])
            normalized_retention_time.append(get_RT(ms1, query[1], ppm_error))
            transition_group_id.append(query[0])
            collision_energy.append(10)

    df = pd.DataFrame(
        {
            "CompoundName": compound_name,
            "PrecursorMz": precursor_mz,
            "ProductMz": product_mz,
            "LibraryIntensity": library_intensity,
            "NormalizedRetentionTime": normalized_retention_time,
            "TransitionGroupId": transition_group_id,
            "CollisionEnergy": collision_energy,
        }
    )
    path = Path(
        "assay-libraries",
        f"{Path(mzML_file).stem}_{Path(mass_list_file).stem}_top{top_n}.tsv",
    )
    if not df.empty:
        df.to_csv(
            path,
            sep="\t",
        )
        st.success(f"Done! {str(path)}")
    else:
        st.warning(
            "Generated library was empty, please check your inputs. DDA data file required."
        )
