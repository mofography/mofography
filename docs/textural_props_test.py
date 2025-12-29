# %%

import numpy as np
import os
import mofography as mgr
import plotly.graph_objects as go
import dask
from scipy import stats
from dask import delayed
from dask.distributed import Client, progress
from pathlib import Path
import stackprinter
import pandas as pd
import pickle

stackprinter.set_excepthook(style="full")
# %%

# get dat from mofdb for iza
iza_df = await mgr.get_data_from_mofdb(
    params=dict(database="IZA", Vf_min=0, Vf_max=1.0)
)
iza_atoms = mgr.atoms_from_mofdb_df(iza_df)

dealyed_reults = []

# Dask cluster locally
# client.shutdown()
client = Client()

spacing = 0.25
mask = 0.5
tabs = 1.0


# workflow
delayed_results = []
for mat_atoms in iza_atoms:
    dgrid = delayed(mgr.dgrid_from_atoms_rust)(spacing=spacing, ase_atoms=mat_atoms)
    results = delayed(mgr.do_periodic_segmentation_rust)(
        dgrid=dgrid,
        mat_atoms=mat_atoms,
        spacing=spacing,
        mask_thickness=mask,
        threshold_abs=tabs,
        h=0.5,
        return_dgrid=True,
    )
    delayed_results.append(results)

# test
delayed_results[0].compute()

futures = dask.persist(*delayed_results)

progress(futures)

collected_results = client.compute(futures)
iza_results = [
    r.result() if r.status == "finished" else None for r in collected_results
]

print(f"these amany are none: {(np.array(iza_results) == None).sum()}")

import pickle

path_to_data = Path("/Users/arungopalan/Datasets/")
# write the iza feather file
iza_df.to_feather(path_to_data / "iza.feather")
# writ the atoms as pkl
print("Done with writing iza df to data folder")
pickle.dump(iza_atoms, open(path_to_data / "iza_atoms.pkl", "wb"))
print("Done with writing iza atoms to pickle folder")
# dump the iza results to pickle as well
pickle.dump(iza_results, open(path_to_data / "iza_results.pkl", "wb"))
print("Done with writing iza results to data folder")


# compute lcd pld dfsd  from dgrid, super fast
pore_dias = []
for i, results in enumerate(iza_results):
    if results is not None:
        mat_atoms = iza_atoms[i]
        dgrid = results["dgrid"]
        lcd, pld, dfsd, _ = mgr.get_pore_dia_mgr(mat_atoms, dgrid=dgrid)
        pore_dias.append((lcd, pld, dfsd))
    else:
        pore_dias.append((None, None, None))

pore_dias = np.vstack(pore_dias)
lcds = pore_dias[:, 0]
plds = pore_dias[:, 1]
dfsds = pore_dias[:, 2]


variable = lcds
variable = variable.astype(float)
variable = variable[np.isfinite(variable)]
hist, bin_edges = np.histogram(
    variable, bins=50, range=(0, variable.max()), density=True
)
binC = 0.5 * (bin_edges[1:] + bin_edges[:-1])
# hist trace
hist_trace = go.Bar(
    x=binC,
    y=hist,
    name="IZA LCD values (\u212b)",
    marker_color="skyblue",
    opacity=0.9,
)
# kde
kde = stats.gaussian_kde(variable)  # let's add a kde
kde_x = np.linspace(0, variable.max(), 200)
kde_y = kde(kde_x)
kde_trace = go.Scatter(
    x=kde_x,
    y=kde_y,
    mode="lines",
    name="KDE",
    line=dict(color="blue", width=2, dash="dash"),
)

fig = go.Figure()
fig.update_layout(template="simple_white")
fig.update_layout(
    xaxis=dict(
        zeroline=True,
        # range=[0,variable.max()],
        visible=True,
        showline=True,
        mirror=True,
        showgrid=False,
        title="plot title",
    ),
    yaxis=dict(
        zeroline=True,
        showgrid=False,
        visible=True,
        showline=True,
        range=[0, 1],
        mirror=True,
        title="Density",
    ),
    showlegend=True,
    legend=dict(x=0.7, y=0.75, groupclick="togglegroup"),
    width=1000,
    height=600,
)
fig.add_trace(hist_trace)
fig.add_trace(kde_trace)
fig.show(renderer="firefox")
# fig.write_image("file.svg")
#


# write iza cifs to files for zeo++
import pickle

iza_df = pd.read_feather("~/Datasets/iza.feather")
path_to_cifs = Path("~/Datasets/IZA_cifs").expanduser()
for i in range(0, len(iza_df)):
    name = iza_df["name"][i]
    cif = iza_df["cif"][i]
    with open(path_to_cifs / f"{name}.cif", "w") as f:
        f.write(cif)
