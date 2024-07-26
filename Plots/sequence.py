# Imports

from pathlib import Path

import pandas as pd

import matplotlib as mpl
from matplotlib import pyplot as plt

# Set fixed parameters

mpl.use("pdf")

mpl.rcParams["text.usetex"] = True
mpl.rcParams["font.family"] = "serif"

cm = 1 / 2.54
mpl.rcParams["figure.figsize"] = [14.00 * cm, 8.00 * cm]
mpl.rcParams["figure.constrained_layout.use"] = True

mpl.rcParams["legend.frameon"] = False

# Set auxiliary variables

bp_part = 33 * 200

chr_len = {
    "chr2L": 23_011_544,
    "chr2R": 21_146_708,
    "chr3L": 24_543_557,
    "chr3R": 27_905_053,
    "chr4": 1_351_857,
    "chrX": 22_422_827,
}

# Load data into a dataframe

simrootdir = Path("Bibliography")
datafilepath = simrootdir / "LADs.txt"

df_seq = pd.read_csv(datafilepath, sep="\\s+", header=None)

df_seq.columns = ["chr", "beg", "end", "_"]
df_seq = df_seq.drop(columns="_")

# Transform sequence data

chr_len_sum = 0
chr_lim = []

for chr in chr_len:
    df_seq.loc[df_seq["chr"] == chr, "beg"] += chr_len_sum
    df_seq.loc[df_seq["chr"] == chr, "end"] += chr_len_sum
    chr_len_sum += chr_len[chr]
    chr_lim.append(chr_len_sum)

print("chr_len_sum = {:_}".format(chr_len_sum))

df_cg_seq = df_seq.copy()
df_cg_seq["beg"] = df_cg_seq["beg"] // bp_part
df_cg_seq["end"] = df_cg_seq["end"] // bp_part

n_part = chr_len_sum // bp_part

print("n_part = {:_}".format(n_part))

df_seq["len"] = df_seq["end"] - df_seq["beg"]
df_seq["beg_len"] = df_seq[["beg", "len"]].apply(tuple, axis=1)

df_cg_seq["len"] = df_cg_seq["end"] - df_cg_seq["beg"]
df_cg_seq["beg_len"] = df_cg_seq[["beg", "len"]].apply(tuple, axis=1)

seq_ratio = df_seq["len"].sum() / chr_len_sum
print("seq_ratio = {:.4f}".format(seq_ratio))

cg_seq_ratio = df_cg_seq["len"].sum() / n_part
print("cg_seq_ratio = {:.4f}".format(cg_seq_ratio))

# Make sequence plot

plotsrootdir = Path("Plots")

fig, ax = plt.subplots(2)

ax[0].tick_params(left=False, labelleft=False)
ax[1].tick_params(left=False, labelleft=False)

ax[0].set_xlabel("par de bases")
ax[1].set_xlabel("índice de la partícula")

ax[0].set_ylabel("Drosophila LADs")
ax[1].set_ylabel("secuencia modelo")

ax[0].broken_barh(df_seq["beg_len"], (0, 1), linewidth=0.5, color="#d37b81")
ax[1].broken_barh(df_cg_seq["beg_len"], (0, 1), linewidth=0.5, color="#d37b81")

for lim in chr_lim:
    ax[0].axvline(lim, linestyle="--", linewidth=0.5, color="gray")
    ax[1].axvline(lim // bp_part, linestyle="--", linewidth=0.5, color="gray")

ax[0].set_xlim(0, chr_len_sum)
ax[0].set_ylim(0, 1)

ax[1].set_xlim(0, n_part)
ax[1].set_ylim(0, 1)

fig.savefig(plotsrootdir / "sequence.pdf")
