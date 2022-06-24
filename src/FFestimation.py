#!/usr/bin/python

"""
Python 3.8
Mother homozigous for reference allele and foetus got alternate allele from the father 
VAFp = 1/2 FF

Mother homozygous for a alternate allele and fetus got reference allele
VAFm = 1 - (0,5 * FF)

FF = VAFpaverage+(1-VAFm average)
https://doi.org/10.3390/biotech10030017

"""

__authors__ = "Jean-Baptiste Lamouche"
__contact__ = "jean-baptiste.lamouche@chru-strabsourg.fr"
__copyright__ = "MIT"
__date__ = "2022-02-01"
__version__ = "0.1"


from functools import lru_cache
from os.path import join as osj
from pyfiglet import Figlet
from io import StringIO
import argparse
import plotly.graph_objects as go
import plotly.express as px
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pysam
import time
import subprocess
import statistics as sts
import sys
import json

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


class Process:
    def __init__(self, mother, father, foetus, filetype, output, filterqual):
        m, d, f = self.preprocess(mother, father, foetus, filetype, output, filterqual)
        self.mother = m.loc[m["varType"] == "substitution"]
        self.father = d.loc[d["varType"] == "substitution"]
        self.foetus = f.loc[f["varType"] == "substitution"]
        self.filetype = filetype
        print("#[INFO] Length mother variants ", len(self.mother.index))
        print("#[INFO] Length father variants ", len(self.father.index))
        print("#[INFO] Length foetus variants ", len(self.foetus.index))
        # self.header = {}
        # fam = ["mother", "father", "foetus"]
        # for i, path in locals().items():
        #    if i in fam:
        #        name = os.path.basename(path).split(".")[0]
        #        self.header[name] = getheader(path)
        self.output = output
        self.filterqual = filterqual

    def preprocess(self, mother, father, foetus, filetype, output, filterqual):
        """
        input: tsv of variants from trio family
        output: pickle object for whole family
        """
        print("#[INFO] Locals init items ", locals())
        # make a pickle object in output folder if it not exists, even where the vcf are located
        for j, values in locals().items():
            if j in ["mother", "father", "foetus"]:
                path = os.path.basename(values)
                print("#[INFO] " + j + " " + path)
                filename = osj(output, path + ".pickle")
                if not os.path.exists(filename):
                    if filetype == "vcf":
                        vcfTodataframe(filename)
                    else:
                        tsvTodataframe(values)

        # FOR TSV only #TODO
        columns_list = ["alleleFrequency"]
        if filterqual:
            # print(
            #    "#[INFO] Filterquality remove varReadPercent < 30 varReadDepth < 30 and qualphred < 300"
            # )
            mother = self.filterquality(
                self.col_type(pd.read_pickle(mother + ".pickle"), columns_list), output
            )
            father = self.filterquality(
                self.col_type(pd.read_pickle(father + ".pickle"), columns_list), output
            )
        else:
            mother = self.col_type(pd.read_pickle(mother + ".pickle"), columns_list)
            father = self.col_type(pd.read_pickle(father + ".pickle"), columns_list)

        foetus = self.col_type(pd.read_pickle(foetus + ".pickle"), columns_list)

        # filter_foetus, filter_foetus_homo = self.processtsv(self.filterquality(mother), self.filterquality(father), foetus)
        # father.to_csv(
        #    osj(output, "paternal_test.tsv"), sep="\t", header=True, index=False
        # )
        # print(foetus[["alleleFrequency"]].head())
        return mother, father, foetus

    def col_type(self, df, columns_list):
        # print(columns_list)
        for col in columns_list:
            # print(col)
            # print(df.dtypes[col])
            if df.dtypes[col] != "float64":
                df.loc[:, col] = (
                    df[col].astype(str).str.replace(",", ".").astype("float")
                )
        return df

    def filterquality(self, df, output):
        # old
        # filter = df.loc[
        #    (df["varReadPercent"] > 30)
        #    & (df["varReadDepth"] > 30)
        #    & (df["QUALphred"] > 300)
        #    & (df["alleleFrequency"] < 0.02)
        # ]
        # print(df.dtypes["alleleFrequency"])
        # print(df.dtypes["totalReadDepth"])
        filter = df.loc[
            (df["varReadPercent"] > 20)
            & (df["totalReadDepth"] > 20)
            & (df["QUALphred"] > 300)
            # & (df["alleleFrequency"] < 0.02)
        ]
        # print(df["alleleFrequency"])
        config_filter = {}
        cols = ["varReadPercent", "totalReadDepth", "QUALphred", "alleleFrequency"]
        for col in cols:
            if col in filter.columns and filter[col].dtypes in ["int64", "float"]:
                config_filter[col] = [
                    str(min(filter[col].to_list())),
                    str(max(filter[col].to_list())),
                ]
        with open(osj(output, "filter.json"), "w+") as j:
            json.dump(config_filter, j)
        return filter


class Paternalidentification(Process):
    def __init__(self, mother, father, foetus, filetype, output, filterqual):
        super().__init__(mother, father, foetus, filetype, output, filterqual)

    def subtract_maternal(self):
        """
        input: dataframe of maternal blood (unknown FF) and mother as Index Case --> 100 %
        output: dataframe keep only variants from Father, foetus and potentially denovo
        """

        maternal_id_SNV = self.mother.loc[self.mother["varType"] == "substitution"][
            "variantID"
        ].to_list()
        # print(maternal_id_SNV)
        maternal_id_InDels = self.mother.loc[self.mother["varType"] != "substitution"][
            "cNomen"
        ].to_list()
        df_subtract_snp = self.foetus.loc[
            (~self.foetus["variantID"].isin(maternal_id_SNV))
            & (self.foetus["varType"] == "substitution")
        ]
        # print(df_subtract_snp.head())
        df_subtract_indels = self.foetus.loc[
            (~self.foetus["cNomen"].isin(maternal_id_InDels))
            & (self.foetus["varType"] != "substitution")
        ]
        df_subtract = pd.concat(
            [df_subtract_snp, df_subtract_indels], ignore_index=True
        )
        # print(df_subtract.head())
        # print(df_subtract.columns)
        df_subtract.sort_values("varReadPercent", inplace=True)

        # foetus_from_m = foetus.loc[(foetus['variantID'].isin(maternal_id_SNV)) | (foetus['cNomen'].isin(maternal_id_InDels))]
        print(
            "#[INFO] Length after removing mother variations ", len(df_subtract.index)
        )
        print("TSV after removal ", len(df_subtract.index))
        # df_subtract.to_csv('SUBTRACT.tsv', sep='\t', index=False, columns=df_subtract.columns)

        return df_subtract

    def identify_paternal(self, foetus_filter):
        """
        input: fetal dataframe where commom variant with mother have been discarded
        output: try to estimate denovo variant and FF
        """

        # TEST only var which are in 40 - 60 AF in father or supp to 90 to avoid dilution artefact
        # paternal_filter = self.father.loc[
        #    (self.father["varReadPercent"] > 40) & (self.father["varReadPercent"] < 60)
        #    | (self.father["varReadPercent"] > 90)
        # ]
        paternal_filter = self.father

        # select variants which are present in father
        paternal_id_SNV = paternal_filter.loc[
            paternal_filter["varType"] == "substitution"
        ]["variantID"].to_list()
        paternal_id_InDels = paternal_filter.loc[
            paternal_filter["varType"] != "substitution"
        ]["cNomen"].to_list()

        # select snp then indels and finally merge both filter, select only het in pool
        paternal_snv = foetus_filter.loc[
            (foetus_filter["variantID"].isin(paternal_id_SNV))
            & (foetus_filter["varType"] == "substitution")
            & (foetus_filter["zygosity"] == "het")
        ]
        paternal_indels = foetus_filter.loc[
            (foetus_filter["cNomen"].isin(paternal_id_InDels))
            & (foetus_filter["varType"] != "substitution")
            & (foetus_filter["zygosity"] == "het")
        ]

        paternal = pd.concat([paternal_snv, paternal_indels], ignore_index=True)

        print(
            "#[INFO] Variants common between denovo pool removing maternal and 100pct foetal "
            + str(len(paternal.index))
        )
        # Probably denovo
        # denovo = foetus_filter.loc[(foetus_filter["varReadPercent"] > 0.075)]

        # print("#[INFO] Variants comming from father between 4 and 15percent of FF ", len(paternal.index))
        # print("#[INFO] After all filter, probably denovo variants (normally high for now cuz using 100percent ES foetus ", len(denovo.index),)
        print(paternal.head())
        print(paternal["varReadPercent"].max())
        print(paternal["varReadPercent"])
        return paternal

    def getFF(self, paternal, foetus_filter):
        # equal to paternal if above func TODO
        # df = pd.concat([foetus_from_m, foetus_from_p], axis=1, ignore_index=True)
        # var = []
        # for i, variants in paternal.iterrows():
        # 	if variants['variantID'] in foetus_filter['variantID'].to_list():
        # 		if variants['varReadPercent'] < 0.075:
        # 			var.append(variants)
        #
        # print("#[INFO] LEN VAR ",len(var))
        ##print(var[0:5])
        ##df = foetus_from_p.loc[foetus_from_p['variantID']]
        # df = pd.DataFrame(var)
        # remove double line for homozygote variant
        # df_filter = df.drop_duplicates(subset=['variantID', 'cNomen'], keep='Last')
        ff = round(average(paternal["varReadPercent"]), 2) * 2
        print("#[INFO] FF estimation: " + str(ff))
        # paternal.to_csv('test_ff.tsv', sep='\t', columns=paternal.columns, index=False)
        return ff

    def main_paternal(self):
        sub = self.subtract_maternal()
        paternal_var = self.identify_paternal(sub)
        paternal_var.sort_values("varReadPercent", ascending=True, inplace=True)

        # In maternal blood keep only rare variant(probably patho) popfreq < 0.02 #TODO
        # paternal_var_rare = paternal_var.loc[paternal_var["alleleFrequency"] < 0.02]
        paternal_var_rare = paternal_var
        # paternal_var_rare.to_csv(
        #    osj(self.output, "paternal.tsv"), sep="\t", header=True, index=False
        # )
        ff = self.getFF(paternal_var_rare, sub)
        # paternal_var_rare["varReadPercent"] = paternal_var_rare["varReadPercent"].apply(
        #    lambda x: round(x / 100, 2)
        # )
        print(paternal_var_rare["varReadPercent"].head())
        plotvaf(paternal_var_rare, self.output)
        # Generate dico for annotation in plot

        dico = series_to_stats(paternal_var_rare["varReadPercent"])
        print(paternal_var_rare["varReadPercent"].head())
        scatter_vaf(paternal_var_rare, self.output, "vafFFestim", dico)
        print(paternal_var_rare["varReadPercent"].head())
        paternal_var_rare.to_csv(
            osj(self.output, "paternalcheck.tsv"), header=True, index=False, sep="\t"
        )


class Homozygotebased(Process):
    def __init__(
        self, mother, father, foetus, filetype, output, filterqual, bedtools, rmasker
    ):
        super().__init__(mother, father, foetus, filetype, output, filterqual)

        # self.mother = self.mother.loc[self.mother["varType"] == "substitution"]
        # self.father = self.father.loc[self.father["varType"] == "substitution"]
        # self.foetus = self.foetus.loc[self.foetus["varType"] == "substitution"]
        # self.dataframe_list = [mother, father, foetus]
        self.output = output
        self.bedtools = bedtools
        self.rmasker = rmasker
        print(self.mother.varType.unique())

    # PURE FETAL FRACTION ESTIMATION based on publciation below
    def globalfilter(self, df, rmasker, pattern):
        """
        Prenatal Testing. BioTech 2021, 10, 17.https://doi.org/10.3390/biotech10030017, parameters read depth and MAF SNP should be common enough to be detected
        Sims, D.; Sudbery, I.; Ilot, N.E.; Heger, A.; Ponting, C.P. Sequencing depth and coverage: Key considerations in genomic analyses.
        Nat. Rev. Genet. 2014, 15, 121–132.
        MAF > 5% and got dbSNP ID (could also use gnomAD stats)
        """

        # MAF dbsnp 5% so common variant in population, and got dbsnp variant btw
        df = col_type_float(df, "rsMAF")
        df["totalReadDepth"] = df["totalReadDepth"].astype("int")
        # tmp = df.loc[(df['rsMAF'] > 0.05) & (~df['rsId'].isnull()) & (df['totalReadDepth'] > 30)]
        # tmp.loc[:, 'chr'] = 'chr'+tmp.chr

        df.loc[:, "chr"] = "chr" + df.chr

        # name of dataframe of variants to bed
        bedname = osj(self.output, pattern)
        foetusfilter = osj(self.output, pattern + ".out")
        # repeatmasker
        bed = self.dataframetoregions(df, bedname, True)
        if not os.path.exists(foetusfilter):
            print("#[INFO] BEDTOOLS Processing ... ")
            systemcall(
                self.bedtools
                + " intersect -v -a "
                + bed
                + " -b "
                + rmasker
                + " -wa > "
                + foetusfilter
            )
        if os.path.exists(foetusfilter):
            foetus_df = pd.read_csv(foetusfilter, sep="\t")
            foetus_df = foetus_df.set_axis(["chr", "start", "stop"], axis="columns")
            filter = df.loc[df["start"].isin(foetus_df["start"].to_list())]
            return filter
        else:
            print(
                "ERROR "
                + foetusfilter
                + " does not exists check your bedtools install exit !"
            )
            exit()

    def getsensibility(self, dfoetal, dfparents):
        # TODO
        return

    def estimateFF(self, filter):
        if filter.dtypes["varReadPercent"] == "float64":
            VAF = filter["varReadPercent"].mean()
            return VAF

        elif filter.dtypes["varReadPercent"] == "int64":
            VAF = filter["varReadPercent"].astype("float").mean()
            return VAF
        else:
            print(filter.dtypes["varReadPercent"])
            VAF = filter["varReadPercent"].str.replace(",", ".").astype("float").mean()
            return VAF

    def processvcf(self):
        # for datafs in self.dataframe_list:
        #    datafs = datafs.loc[
        #        (datafs["REF"] == ".")
        #        | (datafs["ALT"] == ".")
        #        | (datafs["REF"].str.len() > 1)
        #        | (datafs["ALT"].str.len() > 1)
        #    ]
        ##foetus heterozygote avec alternate allele provenant du père
        filter_foetus = foetus.loc[(~foetus["POS"].isin(mother["POS"].to_list()))]
        print(
            "#[INFO] Length alternate variant provide by father (mother is homozygote for reference allele)",
            len(filter_foetus),
        )

        # foetus heterozygote avec alternate allele provenant de la mère sachant variant homozygote et père ayant donné un allèle de reference
        homotmp = mother.loc[mother["ASG2104747"].str.partition(":")[0] == "1/1"]
        homo = homotmp["POS"].to_list()
        print(
            "#[INFO] Length allele from homozygote alternate variant provide by mother (father gave ref allele)",
            len(homotmp),
        )

        filter_foetus_homo = foetus.loc[
            (foetus["POS"].isin(homo))
            & (
                (foetus["FCL2104751"].str.partition(":")[0] == "1/0")
                | (foetus["FCL2104751"].str.partition(":")[0] == "0/1")
            )
        ]
        return filter_foetus, filter_foetus_homo

    def processtsv(self, mother, father):

        ##Discard variants non present in father and mother but present in pool
        # filter_foetus = self.foetus.loc[(~self.foetus['start'].isin(self.mother['start'].to_list())) & (~self.foetus['start'].isin(self.father['start'].to_list()))]
        # print("INFO after no in both remove ", len(filter_foetus))

        ##foetus heterozygote avec alternate allele provenant du père ( second condition keep only heterozygote)
        # print("INFO list mother ", len(list_mother))
        # Check if there are really present in father and with "good AF"
        # convert values to float
        father = col_type_float(father, "varReadPercent")

        # father_tmp = father.loc[
        #    (father["varReadPercent"] > 0.45) & (father["varReadPercent"] < 55)
        #    | (father["varReadPercent"] > 95)
        # ]
        list_father = self.father["start"].to_list()

        list_mother = self.mother["start"].to_list()
        vafp = self.foetus.loc[
            (~self.foetus.start.isin(list_mother))
            & (self.foetus.start.isin(list_father))
        ]  # & (self.foetus['start'].isin(father['start'].to_list())) & (self.foetus['zygosity'] == 'het')]
        vafp = vafp.loc[vafp["zygosity"] == "het"]

        print(
            "#[INFO] Length alternate variant provide by father (mother is homozygote for reference allele)",
            len(vafp.index),
        )

        # foetus heterozygote avec alternate allele provenant de la mère sachant variant homozygote et père ayant donné un allèle de reference
        homotmp = mother.loc[mother["zygosity"] == "hom"]
        homo = homotmp["start"].to_list()

        vafm = self.foetus.loc[
            (self.foetus["start"].isin(homo)) & (self.foetus["zygosity"] == "het")
        ]
        print(
            "#[INFO] Length allele from homozygote alternate variant provide by mother (father gave ref allele)",
            len(vafm.index),
        )

        # print("MOTHER ", vafp.loc[(mother['start'] == 9783147)])
        return vafp, vafm

    def get_ff(self):
        VAFp_tmp = self.globalfilter(self.father, self.rmasker, "filter_father")
        VAFm_tmp = self.globalfilter(self.mother, self.rmasker, "filter_mother")

        VAFp, VAFm = self.processtsv(VAFm_tmp, VAFp_tmp)
        VAFp.to_csv(osj(self.output, "VAFp.tsv"), sep="\t", header=True, index=False)
        # print("MOTHER ", VAFp.loc[(VAFp['start'] == 9783147)])
        VAFm.to_csv(osj(self.output, "VAFm.tsv"), sep="\t", header=True, index=False)
        print("#[INFO] VAFp df " + str(len(VAFp.index)))
        print("#[INFO] VAFm df " + str(len(VAFm.index)))

        print("#[INFO] VAFp estimate " + str(self.estimateFF(VAFp) / 100))
        print("#[INFO] VAFm estimate " + str(self.estimateFF(VAFm) / 100))

        FF = self.estimateFF(VAFp) / 100 + (1 - self.estimateFF(VAFm) / 100)
        print("#[INFO] Estimation of Foetal fraction : ", FF)

    def get_ff_jiang(self):
        VAFp_tmp = self.globalfilter(self.father, self.rmasker, "filter_father")
        VAFm_tmp = self.globalfilter(self.mother, self.rmasker, "filter_mother")
        VAFp, VAFm = self.processtsv(VAFm_tmp, VAFp_tmp)
        res = {}
        for i, var in VAFp.iterrows():
            # print(var)
            res[var["variantID"]] = float(
                2
                * var["varReadDepth"]
                / (var["varReadDepth"] + (var["totalReadDepth"] - var["varReadDepth"]))
            )
        print(json.dumps(res, indent=4))
        print("FFestimation: ", average(res.values()))
        return average(res.values())

    def dataframetoregions(self, dataframe, bedname, save):
        bed = dataframe.loc[:, ["chr", "start", "end"]]
        if len(bed.index) == 0:
            print("ERROR col are missing from file exit")
            exit()
        # bed.set_axis(['#CHROM', 'START', 'END'], axis=1, inplace=True)
        if save:
            bed.to_csv(bedname, index=False, header=False, sep="\t")
        return bedname

    def getUMI(self, bamfile, position):
        """
        get number of different UMI carring alternate base for a given position #TODO

        """
        zscore = {}
        # bam = pysam.AlignmentFile(bamfile, 'rb')
        for i, variants in position.iterrows():  # "chr1:876499-876499"
            print(variants)
            zscore[variants["START"]] = []
            pos = (
                str(variants["#CHROM"])
                + ":"
                + str(variants["START"])
                + "-"
                + str(variants["END"])
            )
            read = pysam.view(bamfile, pos)
            for reads in read.split("\n"):
                fields = reads.split("\t")
                # print(type(fields[0]))
                for items in fields:
                    if (
                        items.startswith("RX")
                        and items not in zscore[variants["START"]]
                    ):
                        print(items)
                        zscore[variants["START"]].append(items.split(":")[-1])

            print(zscore)
            exit()
        return


class Seqff:
    def __init__(
        self, bamfoetus, output, samtools, seqff
    ):  # mount, image, k_fold, container_name
        self.bamfoetus = bamfoetus
        self.output = output
        self.samtools = samtools
        self.seqff = seqff
        # self.client = docker.APIClient(base_url='unix://var/run/docker.sock', timeout=10000)
        # self.config = self.client.create_host_config(binds=mount)
        # self.image = image
        # self.k_fold = k_fold
        # self.container_name = container_name

    def readprofile(self):
        bam = self.bamfoetus
        bamname = os.path.basename(bam).split(".")[0]
        stats = osj(self.output, "stats_samtools_" + bamname)
        if not os.path.exists(stats):
            print("#[INFO] Generate Reads profile " + stats)
            systemcall(self.samtools + " stats " + bam + " > " + stats)
        else:
            print("Warning " + stats + " already exists !")
        # col stats file | insert size, pairs total, inward oriented pairs, outward oriented pairs, other pairs
        dico = {}
        if os.path.exists(stats):
            with open(stats, "r") as s:
                for lines in s:
                    fields = lines.strip()
                    if fields.startswith("IS"):
                        col = fields.split("\t")
                        # taking pairs total column
                        dico[col[1]] = [col[2]]
        tmp = pd.DataFrame.from_dict(dico, dtype="float64")

        # For ML model we need 170 first features, read size from 50 to 220
        df_stats = tmp.iloc[:, 50:220].copy()
        # df_stats.loc[0, '220'] = 0.0910165153462862
        bamname = os.path.basename(bam).split(".")[0]
        df_stats.to_csv(
            osj(self.output, bamname + ".tsv"), sep="\t", header=False, index=False
        )
        return osj(self.output, bamname + ".tsv")

    def analysis(self):
        if self.bamfoetus.endswith(".sam"):
            self.launch_seqff(self.bamfoetus, self.output, "sam")
        elif self.bamfoetus.endswith(".tsv"):
            self.launch_seqff(self.bamfoetus, self.output, "counts")
        else:
            bf = self.readprofile()
            self.launch_seqff(bf, self.output, "counts")

    def launch_seqff(self, foetus, output, format):
        name = foetus.split(".", 1)[0]
        systemcall(
            "Rscript "
            + self.seqff
            + " --i="
            + os.path.dirname(foetus)
            + " --j="
            + os.path.basename(foetus)
            + " --o="
            + name
            + " --d="
            + output
            + " --t="
            + format
            + " --s="
            + os.path.dirname(self.seqff)
            + " &"
        )

        # SEQFF analysis
        # ml_folder = osj(output, 'data')
        # if not os.path.exists(ml_folder):
        # 	print("#[INFO] Create "+ml_folder)
        # 	os.makedirs(ml_folder)
        #
        #
        # 	pre = "python3 preproces.py /DATA/TEST/combo_ff/example/dataset.tsv -k "+self.k_fold+" -s 19 -O "+ml_folder
        # 	self.runcontainer(pre)
        #
        # 	#Training ML
        # 	train = self.generatetrain(self, ml_folder)
        # 	self.runcontainer(train)
        #
        # 	#Aggregate testing results
        # 	systemcall('cat '+osj(ml_folder, "result")+"/*.tsv > "+osj(ml_folder, "result_fl.tsv"))
        # 	systemcall('cat '+osj(ml_folder, "combo", "result")+"/*.tsv > "+osj(ml_folder, "result_combo.tsv"))
        #
        #
        # 	#Predict sample
        # 	predict = self.generatepredict(self, output, sample)
        # 	self.runcontainer(predict)
        #
        # return df_stats

    def runcontainer(self, cmd):
        self.container = self.client.create_container(
            self.image,
            command=cmd,
            user="root",
            detach=True,
            tty=True,
            name=self.container_name,
            entrypoint="/bin/bash",
            volumes=self.mount,
            host_config=self.config,
        )
        self.client.start(container=self.container)
        while (
            self.client.inspect_container(container=self.container)["State"]["Status"]
            == "running"
        ):
            time.sleep(30)
        self.client.remove_container(container=self.container, force=True)
        return

    def generatetrain(self, output):
        val = []
        for i in range(self.k_fold):
            val.append(
                "python3 train_fl.py "
                + osj(output, "train_dataset_" + i + ".tsv")
                + " "
                + osj(output, "test_dataset_" + i + ".tsv")
                + " -o "
                + osj(output, "fl", "model", "model_" + i)
                + " -c "
                + osj(output, "combo", "coeff", "coeffs_" + i + ".txt")
                + " -r "
                + osj(output, "combo", "result", "result_" + i + ".tsv")
                + " >> "
                + osj(output, "predict.tsv")
            )
        return " ".join(val)

    def generatepredict(self, output, sample):
        val = []
        for i in range(self.k_fold):
            val.append(
                "python3 /DATA/TEST/predict.py "
                + sample
                + " "
                + osj(output, "fl", "model", "model_" + i)
                + " -c "
                + osj(output, "combo", "model", "model_" + i)
                + " -m "
                + osj(output, "data", "train_dataset_mean_" + i + ".txt")
                + " -s "
                + osj(output, "data", "train_dataset_std_" + i + ".txt")
                + " -v"
            )
        return " ".join(val)


def parseargs():  # TODO continue subparser and add ML docker in script
    parser = argparse.ArgumentParser(
        description="FFestimation.py aims to analyze family variants especially in case of prenatal diagnosis. Developped in Strasbourg CHRU in the genetical diagnosis lab by bioinformatic team, it will be use as a routine tool integrated in STARK module"
    )
    subparsers = parser.add_subparsers()
    parser_a = subparsers.add_parser(
        "standard", help="Standard use of FFestimation.py TRIO"
    )
    parser_a.add_argument(
        "-d",
        "--dad",
        type=str,
        help="Absolute path of vcf variant from father",
        required=True,
    )
    parser_a.add_argument(
        "-f",
        "--foetus",
        type=str,
        help="Absolute path of vcf variant from cell free DNA, (maternal blood)",
        required=True,
    )
    parser_a.add_argument(
        "-m",
        "--mum",
        type=str,
        help="Absolute path of vcf variant from mother",
        required=True,
    )
    parser_a.add_argument(
        "-t", "--type", default="tsv", type=str, help="vcf or tsv, default tsv"
    )
    parser_a.add_argument(
        "-r",
        "--rmasker",
        default="/home1/BAS/lamouchj/scripts/DPNI/rmsk_norm.txt",
        type=str,
        help="repeatmaskerfile",
    )
    parser_a.set_defaults(func="standard")

    parser_b = subparsers.add_parser(
        "paternal", help="Estimation of FF based on variant comming from father TRIO"
    )
    parser_b.add_argument(
        "-d",
        "--dad",
        type=str,
        help="Absolute path of vcf variant from father",
        required=True,
    )
    parser_b.add_argument(
        "-f",
        "--foetus",
        type=str,
        help="Absolute path of vcf variant from cell free DNA, (maternal blood)",
        required=True,
    )
    parser_b.add_argument(
        "-m",
        "--mum",
        type=str,
        help="Absolute path of vcf variant from mother",
        required=True,
    )
    parser_b.add_argument(
        "-t", "--type", default="tsv", type=str, help="vcf or tsv, default tsv"
    )
    parser_b.add_argument(
        "-b",
        "--bedtools",
        default="/home1/TOOLS/tools/bedtools/current/bin/bedtools",
        type=str,
        help="Abs path of bedtools executable",
    )
    parser_b.add_argument(
        "-r",
        "--rmasker",
        default="/home1/BAS/lamouchj/scripts/DPNI/rmsk_norm.txt",
        type=str,
        help="repeatmaskerfile",
    )
    parser_b.set_defaults(func="paternal")

    parser_c = subparsers.add_parser(
        "seqff",
        help="ML model to estimate FF based on seqFF and regression model using read depth profile Maternal BLOOD only",
    )
    parser_c.add_argument(
        "-bf",
        "--bfoetus",
        type=str,
        help="Absolute path of bam variant from cell free DNA, (maternal blood)",
        required=True,
    )
    parser_c.add_argument(
        "-s",
        "--seqff",
        type=str,
        help="Absolute path of seqFF executable, default /app/seqff",
        default="/app/seqff/pd4615-sup-0002-seqff.r",
    )
    parser_c.add_argument(
        "-m",
        "--mount",
        type=dict,
        default={
            "/home1/BAS/lamouchj/scripts": {"bind": "/DATA"},
            "/home1/data/STARK/databases": {"bind": "/STARK/databases", "mode": "ro"},
        },
        help="Dict of values specify volume to mount",
    )
    parser_c.add_argument(
        "-i",
        "--image",
        type=str,
        default="r-base:4.0.2",
        help="Image of seqFF container",
    )
    parser_c.add_argument(
        "-k",
        "--kfold",
        type=int,
        default=5,
        help="Number of k_fold for cross validation",
    )
    parser_c.add_argument(
        "-c",
        "--container",
        type=str,
        default="seqff_DPNI",
        help="container_name for seqFF ML analysis",
    )
    parser_c.set_defaults(func="seqff")
    parser_d = subparsers.add_parser(
        "plotvaf",
        help="from dataframe with SNP variants, generate scatter plot in HTML format (Plotly library)",
    )
    parser_d.add_argument(
        "-d", "--dataframe", type=str, help="Abs path of TSV file containing variants"
    )
    parser_d.add_argument("-n", "--name", type=str, help="Name of output plot")
    parser_d.set_defaults(func="plotvaf")
    parser.add_argument(
        "-q",
        "--quality",
        action="store_true",
        help="Activate filtering on parents variants file, discarded variants with varReadDepth < 30, varReadPercent < 30 and qual Phred < 300, default True, set arg to remove filtering",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="/home1/data/STARK/data/DPNI/trio/TWIST",
        type=str,
        help="name of outputfolder",
        required=True,
    )
    parser.add_argument(
        "-s",
        "--samtools",
        type=str,
        default="/home1/TOOLS/tools/samtools/current/bin/samtools",
        help="Abs path of samtools executable",
    )
    # TODO create personnal class to personnalize help and so on
    # parser.add_argument('-h', '--help', action=, default=argparse.SUPPRESS,
    #                help='Show this help message and exit.')
    args = parser.parse_args()
    return args


def main():
    if (
        "-h" in sys.argv
        or "--help" in sys.argv
        or len(sys.argv) == 1
        or len(sys.argv) == 2
    ):
        fancystdout("big", "DPNI help")
    args = parseargs()

    fancystdout("speed", "DPNI module")
    # ffname = os.path.basename(args.foetus).split('.')[0]
    if not os.path.exists(args.output):
        print("#[INFO] Create output folder " + args.output)
        os.mkdir(args.output)
    print(args)
    if args.func == "standard":
        # 1) From CDC of DPNI study
        pi = Paternalidentification(
            args.mum, args.dad, args.foetus, args.type, args.output, args.quality
        )
        pi.main_paternal()

    # 2) From combo_FF seqFF modele
    elif args.func == "seqff":
        print("#[INFO] FF estimation by seqff !")
        pseq = Seqff(args.bfoetus, args.output, args.samtools, args.seqff)
        pseq.analysis()

    # 3) From publication with UMI standard deviation
    elif args.func == "paternal":
        p = Homozygotebased(
            args.mum,
            args.dad,
            args.foetus,
            args.type,
            args.output,
            args.quality,
            args.bedtools,
            args.rmasker,
        )
        print(p)
        p.get_ff()
        p.get_ff_jiang()

    elif args.func == "plotvaf":
        scatter_vaf(args.dataframe, args.output, args.name)

    # getUMI("/home1/data/STARK/data/DPNI/trio/TWIST/FCL2104751.bwamem.bam", pos)
    # FF = VAFp + (1 - VAFm)
    # print("#[INFO] Estimation of Foetal fraction : ", FF)


if __name__ == "__main__":
    main()
