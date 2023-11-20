from pyopenms import *
import numpy as np
import pandas as pd


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
    df = pd.read_csv(precursor_file, sep="\t", names=[
                     "name", "mz", "sum formula"])

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
            # Exclude masses, depends on keeping unfractionated precursor mass
            if exclude_precursor_mass:
                condition = mzs < metabolite["mz"] - 1
            else:
                condition = mzs < metabolite["mz"] + delta
            mzs = mzs[condition]
            intys = intys[condition]
            # # Remove peaks with inty < 5% of top max
            # condition = intys > 1000
            # intys = intys[condition]
            # mzs = mzs[condition]
            if mzs.any():
                # Get indices of top_n highest intensity peaks
                sorted_indeces = np.argsort(intys)[-top_n:]
                mzs = mzs[sorted_indeces]
                intys = intys[sorted_indeces]
                # Normalize intensity values
                intys = intys / intys.max()

        return pd.Series({"transition mzs": mzs, "transition intys": intys, "transition ms1 rt": ms1_rt})

    # Apply the get_transitions function to each row
    df[["transition mzs", "transition intys", "transition ms1 rt"]
       ] = df.apply(get_transitions, axis=1)

    # Build transition table
    def build_transition_table(df):
        for _, row in df.iterrows():
            for mz, inty in zip(row["transition mzs"], row["transition intys"]):
                yield row["name"], row["mz"], mz, inty, row["transition ms1 rt"], row["name"], collision_energy

    # Create a DataFrame from the built transition table
    df = pd.DataFrame(np.fromiter(build_transition_table(df), dtype=[("CompoundName", "U100"), ("PrecursorMz", "f"),
                                                                     ("ProductMz", "f"), (
                                                                         "LibraryIntensity", "f"),
                                                                     ("NormalizedRetentionTime", "f"),
                                                                     ("TransitionGroupId",
                                                                      "U100"),
                                                                     ("CollisionEnergy", "i")]))

    # Merge rows where fragment masses are very similar

    return df


def generate_transitions_from_json_data(data, top_n = 0, exclude_precursor_mass=True):
    """data is list of dicts with filtered entries from MassBank"""
    for m in data:

        normalized_merged_intensities = [i/max(m["normalized intensity"]) for i in m["normalized intensity"]]

        # Combine mz and inty into a list of tuples
        combined_data = zip(m["m/z"], normalized_merged_intensities)

        if exclude_precursor_mass:
            combined_data = [d for d in combined_data if d[0] < m["precursor mz"] - 1]

        # Sort the list of tuples based on the values in inty in descending order
        sorted_data = sorted(
            combined_data, key=lambda x: x[1], reverse=True)

        if top_n:
            sorted_data = sorted_data[:top_n]

        for mz, inty in sorted_data:
            yield (m["name"],
                   m["precursor mz"],
                   mz,
                   inty,
                   60,  # Normalized Retention Time, todo
                   # f"{m['library name']}",
                   m["formula"],
                   m["SMILES"],
                   m["InChI"]
                   )


def generate_library_from_json_data(data, top_n, exclude_precursor_mass=False):
    df = pd.DataFrame(np.fromiter(generate_transitions_from_json_data(data, top_n, exclude_precursor_mass),
                                  dtype=[("CompoundName", "U100"),
                                         ("PrecursorMz", "f"),
                                         ("ProductMz", "f"),
                                         ("LibraryIntensity", "f"),
                                         ("NormalizedRetentionTime", "f"),
                                         # ("TransitionGroupId", "U100"),
                                         ("SumFormula", "U100"),
                                         ("SMILES", "U300"),
                                         ("InChI", "U300"),
                                         ]))

    # add unique TransitionGroupId
    df["TransitionGroupId"] = [f"{n}_{i}" for i,
                               n in enumerate(df["CompoundName"])]
    return df.sort_values(by="CompoundName")


def calculate_ppm_distance(value1, value2):
    ppm = abs((value1 - value2) / value1) * 1e6
    return ppm

def filter_duplicate_transitions(df, threshold_ppm=50):
    """Remove transitions where transition m/s's are not unique."""

    indeces_to_drop = set()
    for i in range(len(df)):
        for j in range(i + 1, len(df)):
            ppm_distance_ms1 = calculate_ppm_distance(df['PrecursorMz'][i], df['PrecursorMz'][j])
            ppm_distance_ms2 = calculate_ppm_distance(df['ProductMz'][i], df['ProductMz'][j])
            # Remove rows if ppm distance is within the threshold
            if ppm_distance_ms1 < threshold_ppm and ppm_distance_ms2 < threshold_ppm:
                indeces_to_drop.add(i)
    return df.drop(list(indeces_to_drop))