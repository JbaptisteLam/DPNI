# func utils

from functools import lru_cache
from os.path import join as osj
from pyfiglet import Figlet
from io import StringIO
import plotly.graph_objects as go
import plotly.express as px
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import subprocess
import statistics as sts
import sys


# git combo FF, dossier TEST TODO
def fancystdout(style, text):
    subprocess.call("pyfiglet -f " + style + " '" + text + "'", shell=True)


def scatter_vaf(tsv, output, name, dico):
    if not isinstance(tsv, pd.DataFrame):
        header = []
        with open(tsv, "r") as f:
            for lines in f:
                if lines.startswith("##"):
                    header.append(lines.strip())
                else:
                    break
        # if len(header) > 0:
        #    print("#[INFO] header " + tsv, header)
        df = pd.read_csv(tsv, sep="\t", skiprows=len(header), header=0)
    else:
        df = tsv

    dico = "<br>".join("{}={}".format(*i) for i in dico.items())
    fig = go.Figure()
    fig = px.scatter(
        df,
        x=df.index,
        y=df["varReadPercent"],
        trendline="lowess",
        trendline_options=dict(frac=0.1),
    )
    fig.update_layout(showlegend=False)
    fig.add_trace(
        go.Scatter(
            y=df["varReadPercent"],
            x=df.index,
            customdata=np.stack(
                (
                    df["variantID"],
                    df["cNomen"],
                    df["gene"],
                    df["varReadPercent"],
                    df["alleleFrequency"],
                ),
                axis=-1,
            ),
            mode="markers",
            marker=dict(
                color=df["varReadPercent"],
                colorbar=dict(title="VAF colorscale"),
                colorscale="Turbo",
                size=6,
                symbol="square",
            ),
            hovertemplate="<br>".join(
                [
                    "Chr: %{customdata[0]}",
                    "cNomen: %{customdata[1]}",
                    "Gene: %{customdata[2]}",
                    "AF: %{customdata[3]}",
                    "popFreq: %{customdata[4]}",
                ]
            ),
        )
    )
    fig.update_layout(
        title="Estimation of fetal fraction with allele frequency",
        legend_title_text="VAF",
    )
    fig.add_annotation(
        text="<b>Metrics</b> (all values *2)<br>" + dico,
        xref="paper",
        yref="paper",
        x=0.5,
        y=0.3,
        showarrow=False,
        font=dict(
            family="Courier New, monospace",
            size=20,
        ),
    )
    fig.update_xaxes(title_text="index")
    fig.update_yaxes(title_text="VAF")
    fig.write_html(osj(output, name + ".html"))
    return fig


def series_to_stats(series):
    # Remove multiplication by 2
    # values = series.apply(lambda x: x * 2).to_list()
    values = series.to_list()
    dico = {
        "Mean": sts.mean(values),
        "Median": sts.median(values),
        "Mode": sts.mode(values),
        "Stdev": sts.stdev(values),
        "Variance": sts.variance(values),
    }
    return dico


@lru_cache
def vcfTodataframe(file, rheader=False):

    """
    Take in input vcf file, or tsv and return a dataframe
    I"m gonna build my own vcf parser et puis c'est tout
    return 3 Dataframe, full, only sample, only info
    """
    name, extension = os.path.splitext(file)

    header = []
    variants_tmp = []
    variants = []
    with open(file) as f:
        for lines in f:
            if lines.startswith("#"):
                header.append(lines.strip())
            else:
                variants_tmp.append(lines)

    col = header[-1].strip().split("\t")
    for v in variants_tmp:
        variants.append(v.strip().split("\t"))

    dfVar = pd.DataFrame(columns=col)

    print("#[INFO] Whole VCF to Dataframe")
    for i, var in enumerate(variants):
        print_progress_bar(i, len(variants) - 1)
        rows = pd.Series(var, index=dfVar.columns)
        dfVar = dfVar.append(rows, ignore_index=True)

    print("\n")
    if rheader:
        return dfVar, header
    else:
        name = os.path.basename(file).split(".")[0]
        dfVar.to_pickle(osj(os.path.dirname(file), name + ".pickle"))
        return dfVar


def tsvTodataframe(file):
    df = pd.read_csv(file, sep="\t", skiprows=2, header=0)
    df.to_pickle(file + ".pickle")
    return df


def print_progress_bar(i, max):
    """
    i is your current iteration max should be equal to your maximum number of iterations - 1 (otherwise the bar stop at 99%)
    50 and 100 are arbutrary. 50 defines the size of the bar in the terminal. 100 makes this a percentage.
    The block using this function shoudl be followed by a print("\n") to skip to the next line.
    """
    j = int(i * 50 / max)
    sys.stdout.write("\r")
    sys.stdout.write("[%-50s]%d%%" % ("=" * j, 100 / max * i))
    sys.stdout.flush()


def systemcall(command):
    """
    *passing command to the shell*
    *return list containing stdout lines*
    command - shell command (string)
    """
    print("#[INFO] " + command)
    p = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
    return p.stdout.read().decode("utf8").strip().split("\n")


def average(iterable):
    return round(sum(iterable) / len(iterable), 4)


def getheader(file):
    with open(file, "r") as f:
        header = []
        for lines in f:
            line = lines.strip()
            header.append(line)
    return header


def plotvaf(df, output):
    df["varReadPercent"] = df["varReadPercent"].map(
        lambda varReadPercent: varReadPercent / 100
    )
    # print(df['varReadPercent'].head())
    print("#[INFO] Length paternal " + str(len(df.index)))
    df_plot = df.pivot(columns="varType", values="varReadPercent")
    fig, ax = plt.subplots(figsize=(8, 6))
    df_plot.plot.hist(
        bins=100,
        alpha=0.9,
        ax=ax,
        title="VAF of foetal variations",
        grid=True,
        xticks=np.linspace(0, 1, 11).tolist(),
    )
    ax.annotate(
        "VAF med: " + str(sts.median(df["varReadPercent"].tolist())),
        xy=(2, 1),
        xytext=(3, 1.5),
    )
    ax.annotate(
        "VAF mean: " + str(average(df["varReadPercent"].tolist())),
        xy=(2, 2),
        xytext=(3, 1.5),
    )
    ax.annotate(
        "VAF std: " + str(sts.pstdev(df["varReadPercent"].tolist())),
        xy=(2, 2),
        xytext=(3, 1.5),
    )
    fig.savefig(osj(output, "VAF_plot.jpeg"))
    return osj(output, "VAF_plot.jpeg")


def col_type_float(df, col):
    if df.dtypes[col] == "float64":
        return df
    elif df.dtypes[col] == "int64":
        df[col].astype("float64")
        return df
    else:
        df[col].str.replace(",", ".").astype("float")
        return df
