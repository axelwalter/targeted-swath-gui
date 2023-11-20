import json
import pandas as pd
import pyopenms as poms
import numpy as np
from src.librarygeneration import generate_library_from_json_data, filter_duplicate_transitions

mgf_file = "resources/Supplementary_Table3.mgf"

with open(mgf_file, "r") as f:
    content = f.read()

compounds = content.split("Name: ")[1:]

compound_list = []
for compound in compounds:
    new_compound = {}
    data = compound.split("\n")
    new_compound["name"] = data[0]

    mzs = []
    intys = []

    for line in data[1:]:
        if "Num Peaks" in line:
            peaks_start_index = data.index(line) + 1
            num_peaks = int(line.split(": ")[1])
        elif ": " in line:
            md_key, md_value = line.split(": ")
            try:
                md_value = float(md_value)
            except ValueError:
                pass
            new_compound[md_key] = md_value

    peak_data = data[peaks_start_index:peaks_start_index+num_peaks]

    for peak in peak_data:
        if not peak:
            continue
        mz, inty = peak.split(",")
        intys.append(int(inty) / 1000)
        mzs.append(float(mz))

    new_compound["m/z"] = mzs
    new_compound["normalized intensity"] = intys

    compound_list.append(new_compound)

key_mapping = {
    'Precursor_mz': 'precursor mz',
    'Formula': 'formula',
    'Precursor_type': 'adduct',
    'Spectrum_type': 'MS level',
    'Instrument_type': 'instrument type',
    'InChIKey': 'InChI',
    'Instrument': 'instrument',
    'Ion_mode': 'ion mode',
    'Collision_energy': 'CE',
    'Comments': 'comments'
}

compound_list = [{key_mapping.get(old_key, old_key): value for old_key,
                  value in original_dict.items()} for original_dict in compound_list]

data = [x for x in compound_list if "LC-ESI-QTOF" in x["instrument type"] and x["adduct"] == "[M+H]+"]

template = data[0]
data_merged = []

unique_metabolites = set([d["name"] for d in data])

for unique_metabolite in unique_metabolites:
    exp = poms.MSExperiment()
    for d in [entry for entry in data if entry["name"] == unique_metabolite]:
        spec = poms.MSSpectrum()
        spec.setMSLevel(2)
        spec.setRT(d["precursor mz"])
        p = poms.Precursor()
        p.setMZ(d["precursor mz"])
        spec.setPrecursors([p])
        spec.setMetaValue("name", d["name"])
        spec.setMetaValue("formula", d["formula"])
        spec.setMetaValue("SMILES", d["SMILES"])
        spec.setMetaValue("InChI", d["InChI"])

        for mz, i in zip(d["m/z"], d["normalized intensity"]):
            p = poms.Peak1D()
            p.setIntensity(i)
            p.setMZ(mz)
            spec.push_back(p)
        exp.addSpectrum(spec)
    exp.updateRanges()
    exp.sortSpectra()
    sm = poms.SpectraMerger()
    smp = sm.getParameters()
    smp.setValue("mz_binning_width", 50.0)  # 50 ppm binning width

    # smp.setValue("average_gaussian:ms_level", 2)
    # smp.setValue("mz_binning_width_unit", "Da")
    sm.setParameters(smp)
    try:
        sm.mergeSpectraPrecursors(exp)
    except RuntimeError:
        print(f"WARNING: can't merge spectra for {unique_metabolite}")
    for spec in exp:
        # window mower
        wm = poms.WindowMower()
        wm_params = wm.getDefaults()
        wm_params[b"peakcount"] = 1
        wm_params[b"windowsize"] = 25.0
        wm.setParameters(wm_params)
        wm.filterPeakSpectrumForTopNInJumpingWindow(spec)

        tmp = template.copy()
        tmp["name"] = spec.getMetaValue("name").split("_")[0]  # remove CE info
        tmp["InChI"] = spec.getMetaValue("InChI")
        tmp["SMILES"] = spec.getMetaValue("SMILES")
        tmp["formula"] = spec.getMetaValue("formula")
        tmp["precursor mz"] = spec.getPrecursors()[0].getMZ()
        m, i = spec.get_peaks()
        tmp["m/z"], tmp["normalized intensity"] = m.tolist(), i.tolist()
        data_merged.append(tmp)


# generate library
df = generate_library_from_json_data(data_merged, 5, True)

# filter duplicates
df = filter_duplicate_transitions(df, threshold_ppm=50)

df.to_csv("assay-libraries/Sonnenburg-ESI-QTOF-M+H-SpectraMerger-WindowMower.tsv",
          sep="\t", index=False)
