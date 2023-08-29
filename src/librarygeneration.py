import streamlit as st
from pyopenms import *
import pandas as pd
import numpy as np

@st.cache_data
def generate_library(precursor_file, mzML_file, top_n, exclude_precursor_mass, tolerance_ppm, collision_energy):
    """
    Generate a library of transitions for metabolites.

    Args:
        precursor_file (str): Path to the precursor file (CSV format).
        mzML_file (str): Path to the mzML file.
        top_n (int): Number of top intensity transitions to select.
        exclude_precursor_mass (bool): Whether to exclude precursor masses from transitions.
        tolerance_ppm (float): Mass tolerance in parts per million.
        collision_energy (int): Collision energy value.

    Returns:
        pd.DataFrame: DataFrame containing transition information.
    """
    # Load precursor information from CSV file
    df = pd.read_csv(precursor_file, sep="\t", names=["name", "mz", "sum formula"])

    # Load mzML file into exp
    exp = MSExperiment()
    MzMLFile().load(mzML_file, exp)

    def get_transitions(metabolite):
        """
        Generate transitions for a metabolite.

        Args:
            metabolite (pd.Series): Metabolite information.

        Returns:
            pd.Series: Transition information.
        """
        tic = 0
        mzs = np.array([])
        intys = np.array([])
        delta = (tolerance_ppm / 1000000) * metabolite["mz"]
        print(f"{metabolite['name']} mz: {metabolite['mz']}...")
        ms1_rt = 0
        ms1_rt_tmp = 0  # in case the first spec is MS2
        for spec in exp:
            if spec.getMSLevel() == 1:
                ms1_rt_tmp = spec.getRT()
            if spec.getMSLevel() == 2:
                prec = spec.getPrecursors()[0].getMZ()
                if metabolite["mz"] - delta < prec < metabolite["mz"] + delta:
                    mzs_tmp, intys_tmp = spec.get_peaks()
                    if intys_tmp.sum() > tic:
                        tic = intys_tmp.sum()
                        mzs = mzs_tmp
                        intys = intys_tmp
                        ms1_rt = ms1_rt_tmp
        if mzs.any():
            # Normalize intensity values
            intys = intys / intys.max()
            # Exclude masses, depends on keeping unfractionated precursor mass
            if exclude_precursor_mass:
                condition = mzs < metabolite["mz"] - 1
            else:
                condition = mzs < metabolite["mz"] + delta
            filtered_indices = np.where(condition)[0]
            mzs = mzs[filtered_indices]
            intys = intys[filtered_indices]
            # Get indices of top_n highest intensity peaks
            sorted_indeces = np.argsort(intys)[-top_n:]
            mzs = mzs[sorted_indeces]
            intys = intys[sorted_indeces]

        return pd.Series({"transition mzs": mzs, "transition intys": intys, "transition ms1 rt": ms1_rt})

    # Apply the get_transitions function to each row
    df[["transition mzs", "transition intys", "transition ms1 rt"]] = df.apply(get_transitions, axis=1)

    # Build transition table
    def build_transition_table(df):
        for _, row in df.iterrows():
            for mz, inty in zip(row["transition mzs"], row["transition intys"]):
                yield row["name"], row["mz"], mz, inty, row["transition ms1 rt"], row["name"], collision_energy

    # Create a DataFrame from the built transition table
    df = pd.DataFrame(np.fromiter(build_transition_table(df), dtype=[("CompoundName", "U100"), ("PrecursorMz", "f"),
                                                                     ("ProductMz", "f"), ("LibraryIntensity", "f"),
                                                                     ("NormalizedRetentionTime", "f"),
                                                                     ("TransitionGroupId", "U100"),
                                                                     ("CollisionEnergy", "i")]))
    
    return df