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
mpl.rcParams["figure.figsize"] = [12.00 * cm, 8.00 * cm]
mpl.rcParams["figure.constrained_layout.use"] = True

# Set auxiliary variables

bp_part = 33 * 200
chr_len = {
    "chrI": 15_072_434,
    "chrII": 15_279_421,
    "chrIII": 13_783_801,
    "chrIV": 17_493_829,
    "chrV": 20_924_180,
    "chrX": 17_718_942,
}

# Load data into a dataframe

datadir = Path("Bibliography")

hetfilepath = datadir / "hets.txt"
df_het_seq = pd.read_csv(hetfilepath, sep="\\s+", header=None)
df_het_seq.columns = ["chr", "beg", "end"]

LADfilepath = datadir / "LADs.txt"
df_LAD_seq = pd.read_csv(LADfilepath, sep="\\s+", header=None)
df_LAD_seq.columns = ["chr", "beg", "end"]

# Transform sequence data

chr_len_sum = 0
chr_lim = [0]
for chr in chr_len:
    df_het_seq.loc[df_het_seq["chr"] == chr, "beg"] += chr_len_sum
    df_het_seq.loc[df_het_seq["chr"] == chr, "end"] += chr_len_sum
    df_LAD_seq.loc[df_LAD_seq["chr"] == chr, "beg"] += chr_len_sum
    df_LAD_seq.loc[df_LAD_seq["chr"] == chr, "end"] += chr_len_sum
    chr_len_sum += chr_len[chr]
    chr_lim.append(chr_len_sum)
print("chr_len_sum = {:_}".format(chr_len_sum))
print("chr_lim :", chr_lim)

df_cg_seq = df_LAD_seq.copy()
df_cg_seq["beg"] = round(df_cg_seq["beg"] / bp_part)
df_cg_seq["end"] = round(df_cg_seq["end"] / bp_part)
n_part = round(chr_len_sum / bp_part)
cg_chr_lim = [round(lim / bp_part) for lim in chr_lim]
print("n_part = {:_}".format(n_part))
print("cg_chr_lim :", cg_chr_lim)

df_het_seq["len"] = df_het_seq["end"] - df_het_seq["beg"]
df_het_seq["beg_len"] = df_het_seq[["beg", "len"]].apply(tuple, axis=1)
het_seq_ratio = df_het_seq["len"].sum() / chr_len_sum
print("het_seq_ratio = {:.4f}".format(het_seq_ratio))

df_LAD_seq["len"] = df_LAD_seq["end"] - df_LAD_seq["beg"]
df_LAD_seq["beg_len"] = df_LAD_seq[["beg", "len"]].apply(tuple, axis=1)
LAD_seq_ratio = df_LAD_seq["len"].sum() / chr_len_sum
print("LAD_seq_ratio = {:.4f}".format(LAD_seq_ratio))

df_cg_seq["len"] = df_cg_seq["end"] - df_cg_seq["beg"]
df_cg_seq["beg_len"] = df_cg_seq[["beg", "len"]].apply(tuple, axis=1)
cg_seq_ratio = df_cg_seq["len"].sum() / n_part
print("cg_seq_ratio = {:.4f}".format(cg_seq_ratio))

avg_cg_dom_len = n_part / (2 * df_cg_seq.shape[0])
print("avg_cg_dom_len = {:.5f}".format(avg_cg_dom_len))

# Save coarse grained (diploid) sequence

outputdir = Path("Program/data")

n_part_d = 2 * n_part
cg_chr_lim_d = cg_chr_lim + [n_part + lim for lim in cg_chr_lim[1:]]

print("n_part_d = {:_}".format(n_part_d))
print("cg_chr_lim_d :", cg_chr_lim_d)

filename = outputdir / "sequence.txt"
file = open(filename, "w")
for i_p in range(n_part_d):
    after_beg = df_cg_seq["beg"] <= (i_p % n_part)
    before_end = (i_p % n_part) < df_cg_seq["end"]
    if (after_beg & before_end).sum() > 0:
        file.write("A")
    else:
        file.write("B")
file.close()

# Make sequence plot

plotsdir = Path("Plots")

fig, axl = plt.subplots()

axu = axl.twiny()

axl.tick_params(left=False, labelleft=False)

axu.set_xlabel("par de bases")
axl.set_xlabel("índice de la partícula")

axu.broken_barh(
    df_het_seq["beg_len"],
    (2, 1),
    color="#d81e2c",
    linewidth=0.0,
    label="C. elegans: heterocromatina",
)
axu.broken_barh(
    df_LAD_seq["beg_len"],
    (1, 1),
    color="#a31cc5",
    linewidth=0.0,
    label="C. elegans: LADs",
)
axl.broken_barh(
    df_cg_seq["beg_len"],
    (0, 1),
    color="#194bb2",
    linewidth=0.0,
    label="secuencia del modelo",
)

axu.hlines([1, 2], xmin=0, xmax=chr_len_sum, linewidth=0.5, color="black")

axu.vlines(chr_lim, ymin=1, ymax=3, linestyle="--", linewidth=0.5, color="black")
axl.vlines(cg_chr_lim, ymin=0, ymax=1, linestyle="--", linewidth=0.5, color="black")

axu.set_xlim(0, chr_len_sum)
axl.set_xlim(0, n_part)

axl.set_ylim(0, 3)

axu.legend(loc="upper left")
axl.legend(loc="lower left")

fig.savefig(plotsdir / "sequence.pdf")
