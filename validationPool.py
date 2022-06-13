#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Aim: multisampleanalysis is a STARK module which allowed users to perform STARK analysis with parents or family member have been sequenced in pool sample.
This script generate metrics to assess quality of the analysis. Reference sample: foetus sample has beend diluted in mother sample at different conceentration. In that way variants present in foetus should be present in mother/foetus  sample
"""

import mimetypes
import time
import os
import re
from turtle import color, title
import pandas as pd
import numpy as np
import sys
import datetime
from collections import OrderedDict, defaultdict
from matplotlib.ticker import LinearLocator
import matplotlib.pyplot as plt
import matplotlib as mlp
import math
import argparse
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from os.path import join as osj


def plotSensibilityHisto(
    snvSensibList, indelSensibList, totalSNVList, totalINDELSList, maxFreq, output
):
    # TODO
    # ADD mean of the detection rate

    x = list(range(1, maxFreq + 1))
    y = snvSensibList
    y2 = indelSensibList
    meanSnv = round(sum(snvSensibList) / len(snvSensibList), 2)
    meanIndel = round(sum(indelSensibList) / len(indelSensibList), 2)
    oneAlleleSNV = snvSensibList[0]
    oneAlleleInDels = indelSensibList[0]

    fig, ax = plt.subplots()
    bar1 = ax.bar(x, totalSNVList, edgecolor="blue", color="None")
    bar2 = ax.bar(x, totalINDELSList, edgecolor="red", color="None")

    ax2 = ax.twinx()
    ax3 = ax2.twiny()
    # ax4 = ax2.twinx()
    # ax1 = ax.twinx()
    ax.set_title(
        "Precision of variant calling depending on allele and variant number",
        color="black",
    )
    # ax.grid()

    (line1,) = ax2.plot(x, y, color="blue")
    (line2,) = ax2.plot(x, y2, color="red")
    ax.set_yticks(np.linspace(0, 1.0, 11))
    ax2.get_xaxis().set_visible(False)

    ax.set_yticks(np.linspace(ax.get_ybound()[0], ax.get_ybound()[1], 15))
    ax3.set_yticks(np.linspace(0, 1, 15))
    # ax3.grid()

    ax.set_xticks(list(range(0, maxFreq + 1)))
    ax3.get_xaxis().set_visible(False)

    # Legend
    ax.set_xlabel("Allele Count")
    ax2.set_ylabel("Detection rate")
    ax.set_ylabel("Variant Count")
    plt.legend(
        (line1, line2, bar1[0], bar2[0]),
        ("SNV", "InDels", "SNV count", "InDels count"),
        loc="center right",
    )
    # plt.text(0.75, 0.75, "meanSnv = "+str(meanSnv)+"\nmeanIndels = "+str(meanIndel)+"\noneAlleleSNV = "+str(oneAlleleSNV)+"\noneAlleleInDels = "+str(oneAlleleInDels), ha = 'right', va='bottom')
    plt.text(
        0.75,
        0.75,
        "meanSnv = " + str(meanSnv) + "\nmeanIndels = " + str(meanIndel),
        ha="right",
        va="bottom",
    )

    # ax.set_zorder(3)
    # ax2.set_zorder(2)
    # ax.set_zorder(1)
    plt.style.use("seaborn")
    if not os.path.exists(output):
        os.mkdir(output)

    # if os.path.exists(os.path.join(output, 'plot.png')):
    # 	out = os.path.join(output, 'plot_'+datetime.datetime.today().strftime('%Y%m%d-%H%M%S'))
    # 	print("INFO "+output+" already exists, create file in same folder with date: "+out)
    # else:
    # 	out = os.path.join(output, 'plot.png')

    out = os.path.join(
        output, "plot_" + datetime.datetime.today().strftime("%Y%m%d-%H%M%S") + ".png"
    )
    plt.savefig(out)


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


def vcfTodataframe(file, kind):
    """
    new version of vcf to dataframe faster and handle better huge vcf
    """
    if kind == "vcf":
        data = {}
        skip = []
        with open(file, "r") as f:
            data["header"] = []
            for i, lines in enumerate(f):
                data["header"].append(lines.strip())
                if lines.split("\t")[0] == "#CHROM":
                    print(lines)
                    data["head"] = lines.strip().split("\t")
                    skip.append(i)
                    skiprows = skip[0]
                    break
    else:
        skiprows = 0
    # print(skip)
    chunk = pd.read_csv(
        file, sep="\t", skiprows=skiprows, header=0, chunksize=100000, low_memory=False
    )
    # dtype={'#CHROM': 'str', 'POS': 'int', 'ID':'str', 'REF':'str', 'ALT': 'str', 'QUAL':'str', 'FILTER': 'str', 'INFO': 'str', 'FORMAT': 'str'}
    df = pd.concat(chunk)
    return df


# def vcfTodataframe(file, rheader=False):
#    """
#    DEPRECIATED
#    Take in input vcf file, or tsv and return a dataframe
#    I"m gonna build my own vcf parser et puis c'est tout
#    return 3 Dataframe, full, only sample, only info
#    """
#    name, extension = os.path.splitext(file)
#    header = []
#    variants_tmp = []
#    variants = []
#
#    print("#[INFO] " + file)
#    if extension == ".vcf":
#        # print('[#INFO] VCF: '+file)
#        with open(file) as f:
#            for lines in f:
#                if lines.startswith("#"):
#                    header.append(lines.strip())
#                else:
#                    variants_tmp.append(lines)
#
#    print("#[INFO]", header[-1])
#    # print(variants[-1])
#
#    col = header[-1].strip().split("\t")
#    for v in variants_tmp:
#        variants.append(v.strip().split("\t"))
#
#    # headerCol = [res.replace('#', '') for res in colTemp]
#    dfVar = pd.DataFrame(columns=col)
#    # print(variants[0:3])
#
#    # Creating Dataframe from the whole VCF
#    print("#[INFO] Whole VCF to Dataframe")
#    for i, var in enumerate(variants):
#        print_progress_bar(i, len(variants))
#        rows = pd.Series(var, index=dfVar.columns)
#        # print(rows[0])
#        dfVar = dfVar.append(rows, ignore_index=True)
#
#    print("\n")
#    if rheader:
#        return dfVar, header
#    else:
#        return dfVar


def parseInfoField(dfVar):
    """
    input: take a dataframe (from vcf)

    output: return a dataframe of the vcf when the info field is parsed
    """

    ############
    # Parsing INFO field from dfVar dataframe containing all informations from vcf
    ############

    # print(dfVar.head())
    infoList = []
    dicoInfo = []
    headers = []

    print("#[INFO] Parsing INFO field")
    for i, elems in dfVar.iterrows():
        # print_progress_bar(i, len(dfVar.index)-1)
        infoList.append([x.split("=") for x in elems["INFO"].split(";")])

    print("\n")

    [headers.append(elems[0]) for ite in infoList for elems in ite]
    dfInfo = pd.DataFrame(columns=np.unique(np.array(headers)))
    print(np.unique(np.array(headers)))

    print(infoList[0:5])
    print("#[INFO] From INFO field to Dataframe")
    for j, elems in enumerate(infoList):
        # print_progress_bar(j, len(infoList)-1)
        add = {}
        for fields in elems:
            if len(fields) <= 1:
                f = {fields[0]: "TRUE"}
                add.update(f)
            else:
                f = dict([fields])
                add.update(f)

        dicoInfo.append(add)

    print("\n")
    # print(dicoInfo.keys())
    # print(dict(list(dicoInfo.items())[0:2]))

    df_final = pd.DataFrame(dicoInfo, columns=np.unique(np.array(headers)))

    dfInfo = dfVar.join(df_final, how="inner")

    return dfInfo


def parseField(vcf, df, field, sep):
    with open(vcf, "r") as f:
        targets = [line for line in f if "#CRHOM" in line]

    sampleList = []
    for keys, values in df.iterrows():
        var = {}
        if field == "sample":
            var.fromkeys(targets[-2].split(":"))
            for fi in values:
                print("ok")
    return


def parseSampleField(dfVar):
    #############
    ### Parsing Sample Field in VCF
    #############

    dico = []
    # dfTest = pd.Series(dfVar.TEM195660.values,index=dfVar.FORMAT).to_dict()
    bad_annotation = []

    # Parsing FORMAT field in VCF
    print("\n")
    print("[#INFO] Parsing FORMAT field")

    # index: line where the caller identify an event somethings
    for col in dfVar.columns[9:]:
        print("#[INFO]" + col + "\n")
        for i, row in dfVar.iterrows():
            # print_progress_bar(i, len(dfVar.index)-1)
            # print('\n')
            # print(row['FORMAT'].split(':'), row['bwamem.VarScan_HUSTUMSOL.howard'].split(':'))
            if len(row["FORMAT"].split(":")) != len(row[col].split(":")):
                bad_annotation.append(pd.Series(row[:], index=dfVar.columns))
                continue
            else:
                toadd = pd.Series(
                    row[col].split(":"), index=row["FORMAT"].split(":")
                ).to_dict()
                # toadd.update({"caller":col})
                dico.append(toadd)

    dfSample = pd.DataFrame(dico)
    print("\n")
    df_bad_anno = pd.DataFrame(bad_annotation)
    print("\n")
    df_final = dfVar.join(dfSample, how="inner")

    print(df_final.head())
    print(df_bad_anno.head())
    time.sleep(10)

    return df_final, df_bad_anno


def preprocessVcf(folder, samples):
    """ """

    dfs = []
    j = 1
    for vcf in os.listdir(folder):
        print(j)
        # print(os.path.join(os.getcwd(), "varank"))
        if vcf.endswith(".vcf") and vcf not in samples and "POOL" not in vcf:
            vcf = os.path.join(folder, vcf)
            print("#[INFO] vcf ", vcf)
            header = []

            df_tmp = vcfTodataframe(vcf, "vcf")
            df, bad = parseSampleField(df_tmp)
            # Sample ID add in column

            df["SampleID"] = os.path.basename(vcf)
            df["AlleleCountCohorte"] = np.nan
            for i, var in df.iterrows():
                # not in varank need to check on type or on sample field
                if var["GT"] == "1/1":
                    # print(var["zygosity"])
                    df.loc[i, "AlleleCountCohorte"] = 2
                elif var["GT"] == "1/0" or var["GT"] == "0/1":
                    df.loc[i, "AlleleCountCohorte"] = 1

            dfs.append(df)
            j += 1
        # number of allele it corresponds
    print(len(dfs))
    time.sleep(10)
    return dfs


def read_tsv(tsv):
    print("#[INFO] tsv ", tsv)
    header = []
    with open(tsv) as f:
        for lines in f:
            if lines.startswith("##"):
                header.append(lines)
            else:
                continue
    df = pd.read_csv(tsv, sep="\t", skiprows=len(header))
    return df, header


def preprocessTsv(files, samples, sample_name):
    """ """
    dfs = []
    tsv_list = []
    # if files is a folder containing many index case
    if os.path.isdir(files):
        [tsv_list.append(osj(files, file)) for file in os.listdir(files)]
        # only one TSV index case
    else:
        tsv_list.append(os.path.abspath(files))
    for tsv in tsv_list:
        # print(os.path.join(os.getcwd(), "varank"))
        if tsv.endswith(".tsv") and tsv not in samples and "POOL" not in tsv:
            # tsv = os.path.join(os.getcwd(), folder, tsv)
            df, header = read_tsv(tsv)
            # Sample ID add in column
            for items in header:
                if "Family" in items:
                    df["SampleID"] = items.strip().split(":")[-1]
            # If tsv does not come from varank
            if not "SampleID" in df.columns:
                if sample_name:
                    df["SampleID"] = sample_name
                else:
                    print(
                        "ERROR sample ID columns not found specify option --name exit"
                    )
                    exit()
            df["AlleleCountCohorte"] = np.nan
            for i, var in df.iterrows():
                if var["zygosity"] == "hom":
                    # print(var["zygosity"])
                    df.loc[i, "AlleleCountCohorte"] = 2
                elif var["zygosity"] == "het":
                    df.loc[i, "AlleleCountCohorte"] = 1
            dfs.append(df)
        else:
            continue
            # print(header)
    return dfs


def splitVcfAllele(final_df, length_df):
    """
    input: take a merge vcf from multiple sample, when there are some duplicata,
    variants are merged on the AlleleCountCohrote which represent the number of allele carring the variation inside a cohorte of patiens
    output: return a list of vcf without duplicata split by number of allele
    """
    df_split = []
    # for i in range(len(filterDataframe(final_df).AlleleCountCohorte.unique())):
    alleles = []

    df = final_df.groupby(final_df.AlleleCountCohorte)
    # test = df_test.get_group(24)
    # print(len(test.index))
    # for i in range(1,len([name for name in os.listdir(folder) if os.path.isfile(os.path.join(folder, name) and not 'POOL' in name)])*2+1):
    for i in range(1, length_df * 2 + 1):
        if i in final_df.AlleleCountCohorte.unique():
            frame = df.get_group(i)
        else:
            print("#[INFO] No variants to retrive for this number of allele: ", i)
            frame = pd.DataFrame(columns=final_df.columns)
        df_split.append(frame)

    for dataf in df_split:
        print("#[INFO] df ", dataf.AlleleCountCohorte.unique())
        print("#[INFO] number of variants ", len(dataf.index))
        if dataf.AlleleCountCohorte.unique():
            alleles.append(int(dataf.AlleleCountCohorte.unique()[0]))
        else:
            alleles.append("No data")
    return df_split, alleles


def splitVariants(df_list):
    """
    input: take a list of dataframe
    output: return list SNPs dataframe and a list of INDELs, MNPs dataframe DEPRECIATED
    output return SNPs df and indels df
    indels could be tricky to compare due to the different normalisation of callers (compare on cNOMEN)
    """
    dfs_snv_fin = []
    dfs_others_fin = []
    # Fore each dataframe
    for i, elems in enumerate(df_list):
        # print_progress_bar(i, len(df_list)-1)
        snv = []
        others = []
        # Filter 30 30 300# TODO filter 30 30 300
        # df_filtered = elems.loc[(elems['varReadPercent'] > 30) & (elems['varReadDepth'] > 30) & (elems['QUALphred'] > 300)]
        # print("#[INFO] DF_filtered "+str(len(df_filtered.index)))
        # for each variants
        for i, var in elems.iterrows():
            # print("INFO", var["ref"], var["alt"])
            headers = " ".join(elems.columns)
            # print("INFO headers "+headers)
            # print(re.search("ref", headers, re.IGNORECASE).group())
            if (
                var[re.search("ref", headers, re.IGNORECASE).group()].isalpha()
                and len(var[re.search("ref", headers, re.IGNORECASE).group()]) == 1
                and var[re.search("alt", headers, re.IGNORECASE).group()].isalpha()
                and len(var[re.search("alt", headers, re.IGNORECASE).group()]) == 1
            ):
                # print(var["variantID"])
                # print(type(var))
                snv.append(var)
            else:
                others.append(var)
        df_snv = pd.DataFrame(snv)
        dfs_snv_fin.append(df_snv)
        df_others = pd.DataFrame(others)
        dfs_others_fin.append(df_others)
    print("\n")
    return dfs_snv_fin[0], dfs_others_fin[0]


def getSensi(dfci, dfpool, col):
    """
    input: take a dataframe, containing SNV variants, a dataframe containing all cas index (pool)
    output: return tuple, containing percent of sensitivity and total number of variants to retrieve
    """
    ci = []
    pool = []
    miss = []
    count = 0
    total = len(dfci.index)
    dfci = dfci.reset_index()
    # print(dfci.head())
    for i, pos in dfci.iterrows():
        # print("INFO filter ", col)
        # print(len(col))
        if len(col) > 1:
            ci.append([pos[value] for value in col])
        else:
            ci.append(pos[col[0]])

    for j, variants in dfpool.iterrows():
        # snv_pool.append(["chr"+str(variants['chr']), variants['start'] ,variants['ref'], variants['alt']])
        if len(col) > 1:
            pool.append([variants[value] for value in col])
        else:
            pool.append(variants[col[0]])

    for va in ci:
        # if any(pool == va for pool in snv_pool):
        if va in pool:
            count += 1
        else:
            miss.append(va)

    print("#[INFO] Variants to retrieve " + str(col), len(ci))
    print("#[INFO] Variants miss", len(miss))
    if not ci:
        df_miss = pd.DataFrame(columns=dfci.columns)
        dfound = pd.DataFrame(columns=dfci.columns)
        mean = 0
    else:
        if "DP" in dfci.columns:
            col_depth = "DP"
        else:
            col_depth = "totalReadDepth"
        dfci.reset_index(inplace=True)

        df_miss = dfci.loc[dfci[col[0]].isin(miss)]
        dfci.loc[:, col_depth] = pd.to_numeric(dfci.loc[:, col_depth], downcast="float")

        if len(df_miss.index) != 0:
            # if some variants are missing
            # print("#[INFO] dfmiss", df_miss.head(), len(df_miss.index))
            dfound = dfci[~dfci[col[0]].isin(miss)]

        else:
            dfound = pd.DataFrame(dfci, columns=dfci.columns)
        if dfci.empty:
            mean = 0
        else:
            mean = dfci[col_depth].to_numpy().mean()
    try:
        percent = count / total
    except ZeroDivisionError:
        print("\n")
        print("Cannot divide by zero")
        percent = 0

    return (percent, total, dfound, df_miss, mean)


def dataframeMetrics(
    snvRate,
    indelRate,
    snvfound,
    snvList,
    indelsfound,
    indelList,
    meanDepthSNV,
    meanDepthIndels,
    missVariantsSNV,
    missVariantsInDels,
    fp_snv,
    fp_indels,
    totalvar,
    output,
):

    pd.set_option("colheader_justify", "center")  # FOR TABLE <th>

    html_string = """
	<html>
	  <head><title>HTML Pandas Dataframe with CSS</title></head>
	  <link rel="stylesheet" type="text/css" href="validationPoolStyle.css"/>
	  <body>
	    {table}
	  </body>
	</html>
	"""

    tail_html = """
	<style>.dataframe th{
	background: rgb(237,225,32);
	background: radial-gradient(circle, rgba(237,225,32,1) 0%, rgba(233,202,60,1) 100%);
	padding: 10px;
	font-family: arial black;
	font-size: 105%;
	color: #313131;
	border:2px solid white;
	text-align:left !important;
	}</style
	"""
    # print("#[INFO] LOCALS args")
    # print(locals())
    # TODO KEEP
    lengthmissSNV = len(missVariantsSNV.index)
    lengthmissIndels = len(missVariantsInDels.index)
    print("#[INFO] LengthmissSNV", lengthmissSNV)

    # PPV calcs
    PPV_snv = snvfound / (snvfound + fp_snv)
    PPV_indels = indelsfound / (indelsfound + fp_indels)

    dict_miss = {
        "DetectionRateSNV": snvRate,
        "SNV found": snvfound,
        "LengthSNV": math.trunc(snvList),
        "meanDepthSNV": meanDepthSNV,
        "missVariantsSNV": math.trunc(lengthmissSNV),
        "DetectionRateInDels": indelRate,
        "InDels found": indelsfound,
        "LengthInDels": math.trunc(indelList),
        "meanDepthInDels": meanDepthIndels,
        "missVariantsInDels": math.trunc(lengthmissIndels),
        "PPV SNV": PPV_snv,
        "PPV indels": PPV_indels,
    }
    print("#[INFO] Dict miss", dict_miss)
    df = pd.DataFrame(dict_miss, index=["sample"])
    df = df.astype(
        {
            "SNV found": "int64",
            "LengthSNV": "int64",
            "InDels found": "int64",
            "LengthInDels": "int64",
            "missVariantsSNV": "int64",
            "missVariantsInDels": "int64",
        }
    )
    # OUTPUT AN HTML FILE
    print(df.head())
    out = os.path.join(
        output,
        "metrics_" + datetime.datetime.today().strftime("%Y%m%d-%H%M%S") + ".html",
    )
    with open(out, "w+") as f:
        f.write(html_string.format(table=df.T.to_html(classes="mystyle")))
        f.write(tail_html)


# 09/05/2022   remove variant s which are not in this filter and
def filter_fp(df, name, output):
    df["alleleFrequency"] = df["alleleFrequency"].apply(lambda x: x.replace(",", "."))
    df["alleleFrequency"] = df["alleleFrequency"].astype(float)
    df_filter = df.loc[
        (df["varReadPercent"] > 20)
        & (df["varReadDepth"] > 20)
        & (df["alleleFrequency"] < 2)
    ]
    df_opposite = df.loc[
        ~(
            (df["varReadPercent"] > 20)
            & (df["varReadDepth"] > 20)
            & (df["alleleFrequency"] < 2)
        )
    ]
    print("#[INFO] variants total ", len(df.index))
    print("#[INFO] variants removed during filtering ", len(df_opposite.index))
    print("#[INFO] variants after filtering :", len(df_filter.index))
    writeTsv(df_opposite, name, output)
    return df_filter, df_opposite


def toPlot(
    dfci_snv, dfci_indels, dfpool, filt_snv, filt_indels, output, dfmother, covfile
):
    """
    input: a list of dataframe of SNV from 1 allele to the y patients x2 alleles ex: 12 patients from 1 to 24 dataframe,
    1 dataframe containing pool of all the sample in the cohorte
    """
    # percent total found miss mean
    print("#[INFO] Stats SNV")
    snv, total_snv, dfoundSNV, snvmiss, meanS = getSensi(dfci_snv, dfpool, filt_snv)
    print("#[INFO] Stats InDels")
    others, total_indels, dfoundInDels, indelsmiss, meanI = getSensi(
        dfci_indels, dfpool, filt_indels
    )

    print("\n")

    # Allelecountcohorte
    if output:
        # merge miss
        df_miss = pd.concat([snvmiss, indelsmiss], ignore_index=True)
        df_miss.sort_values(by="varankVarScore", inplace=True, ascending=False)

        writeTsv(dfoundSNV, "SNVfound", output)
        writeTsv(dfoundInDels, "InDelsfound", output)

        # Faux pos pool
        positif = pd.concat([dfoundSNV, dfoundInDels, df_miss], ignore_index=True)
        print("#[INFO] positif ", len(positif.index))
        positif.drop_duplicates(inplace=True)
        print("#[INFO] after removing duplicate positif ", len(positif.index))
        print("#[INFO] Length POOL ", len(dfpool.index))

        # TODO apply func on pool and mother to have good falsepos
        # total variant to retrieve
        total_ci = pd.concat([dfci_snv, dfci_indels], ignore_index=True)
        # filter mother 20 20 2
        dfmother_proper, dfmother_opposite = filter_fp(
            dfmother, "mother_filter_opposite", output
        )
        writeTsv(dfmother_opposite, "mother_filter_opposite", output)
        # print(dfmother_opposite.head())
        # print(*dfmother_proper.columns)

        # obtain variant of mother only
        variants_mother_only = merge_df(
            total_ci, dfmother_proper, "variantID", "cNomen", "outer", "right_only"
        )
        writeTsv(variants_mother_only, "mother_only", output)
        # remove variants of mother in pool
        dfpool_filter = merge_df(
            dfpool, variants_mother_only, "variantID", "cNomen", "outer", "left_only"
        )

        writeTsv(dfpool_filter, "pool_without_mother", output)
        false_pos = merge_df(
            dfpool_filter, positif, "variantID", "cNomen", "outer", "left_only"
        )
        print("#[INFO] Length falsepos ", len(false_pos.index))

        # Split snv indels
        false_pos_snv = false_pos.loc[false_pos["varType"] == "substitution"]
        false_pos_indels = false_pos.loc[false_pos["varType"] != "substitution"]
        writeTsv(false_pos_snv, "FalsePositive_SNV", output)
        writeTsv(false_pos_indels, "FalsePositive_InDels", output)

        # Metrics to return then present in html page

        # print(snv_list_percent,others_list_percent, total_snv_list, total_others_list, len(snv_list_percent))

        # generate informations on miss variants
        df_miss_process = add_cov(covfile, df_miss)
        df_miss_final = identify_miss(df_miss_process)
        viz_go(df_miss_final, osj(output, "res_miss.html"))

        writeTsv(positif, "all_positif_var", output)
        writeTsv(df_miss_final, "missVariants", output)

    # uncomment to plot TODO
    # plotSensibilityHisto(snv_list_percent,others_list_percent, total_snv_list, total_others_list, len(snv_list_percent), output)
    return (
        snv,
        others,
        len(dfoundSNV.index),
        total_snv,
        len(dfoundInDels.index),
        total_indels,
        meanS,
        meanI,
        snvmiss,
        indelsmiss,
        len(false_pos_snv.index),
        len(false_pos_indels.index),
    )


def add_cov(tsv, df):
    """
    from gatk3.8 depth of coverage
    """
    cov = tsvtodataframe(tsv, 0, False)
    # print(cov.head())
    # print(*cov.columns)
    # df = tsvtodataframe(variants, 0, False)
    list_var = []
    base_counts = {"A": [], "C": [], "G": [], "T": [], "N": [], "D": []}
    for i, var in df.iterrows():
        v = "chr" + var["chr"] + ":" + str(var["start"])
        list_var.append(v)
    df["merge_on"] = list_var
    true_var = cov.loc[cov["Locus"].isin(list_var)]
    merge = pd.merge(df, true_var, left_on="merge_on", right_on="Locus")
    for j, vari in merge.iterrows():
        explode = vari["FF5_base_counts"].split(" ")
        for items in explode:
            base_counts[items.split(":")[0]].append(items.split(":")[1])
    ex = pd.DataFrame.from_dict(base_counts)

    final = pd.concat([merge.reset_index(drop=True), ex.reset_index(drop=True)], axis=1)

    lookto = []
    for i, pos in final.iterrows():
        if len(pos["alt"]) == 1 and pos["alt"] != "-":
            lookto.append(pos[pos["alt"]])
        else:
            lookto.append(0)

    final["intend"] = lookto
    final["VAF_miss"] = round(
        final["intend"].astype(int) / final["Total_Depth"].astype(int), 2
    )
    return final


def writeTsv(df, name, output):
    out = os.path.join(
        output,
        name + "_" + datetime.datetime.today().strftime("%Y%m%d-%H%M%S") + ".tsv",
    )
    df.to_csv(out, sep="\t", index=False)
    print("#[INFO] Write TSV " + out)
    return out


def groupBy(df, togroup):
    """
    input: list of dataframe from vcf
    output: return dataframe groupby var 'togroup' which represent a column of the dataframe
    """
    # togroup = ["cNomen"] #and maybe by OMIM ID
    acc = (
        pd.concat([dataf for dataf in df])
        .groupby(togroup, as_index=True)["AlleleCountCohorte"]
        .sum()
    )
    si = (
        pd.concat([dataf for dataf in df])
        .groupby(togroup, as_index=True)["SampleID"]
        .agg(lambda x: ",".join(x))
    )
    final = pd.merge(acc, si, left_index=True, right_index=True)
    return final


def createFinaldf(df_whole, df_newcol, togroup):
    """
    input: list of dataframe from vcf, list of column to add (sample ID, AllelCountCohorte)
    output: return a dataframe without duplicate and containing the new columns
    """
    # merge list of dataframe in once
    all_df = pd.concat([dataf for dataf in df_whole])
    all_df = all_df.set_index(togroup)
    # drop columns to add, which will be replace by the columns for the whole sample
    all_df = all_df.drop(["SampleID", "AlleleCountCohorte"], axis=1)
    all_df = all_df[~all_df.index.duplicated(keep="first")]
    final_df = pd.merge(all_df, df_newcol, left_index=True, right_index=True)
    print(len(final_df))
    # final_df.head()
    return final_df


def tsvtodataframe(file, skiprows, columns):
    """
    from tsv file to pandas dataframe, skiprows number of row to reach header, columns: col which need change (from , to .) to allowed excel filter
    """
    if skiprows:
        df = pd.read_csv(
            file,
            skiprows=skiprows,
            sep="\t",
            header=0,
            chunksize=10000,
            low_memory=False,
        )
    else:
        df = pd.read_csv(file, sep="\t", header=0)
    if "index" in df:
        df = df.drop(columns="index")
    if columns:
        for col in columns:
            if col in df.columns:
                df[col] = df[col].apply(lambda x: x.replace(",", "."))
                df[col] = df[col].astype(float)
    return df


def parse(foetus, mother, colfilter):
    print("INFO foetus ", len(foetus.index))
    filter = foetus.loc[
        (foetus["varReadPercent"] > 30)
        & (foetus["varReadDepth"] > 30)
        & (foetus["QUALphred"] > 300)
    ]
    print("INFO foetus filter ", len(filter.index))
    tmp = pd.merge(filter, mother, how="left", on=colfilter, indicator=True)
    print("INFO testmerge ", len(tmp.index))
    keep = tmp.loc[tmp["_merge"] == "left_only"]
    keep_list = keep[colfilter].to_list()
    print(keep_list)
    print(len(keep_list))
    final = foetus.loc[foetus[colfilter].isin(keep_list)]
    print("INFO final ", len(final.index))
    return final


def splitypeofvar(df):
    # Split snv and indels/mnps
    snv = []
    indels = []
    for i, var in df.iterrows():
        # print("INFO", var["ref"], var["alt"])
        headers = " ".join(df.columns)
        # print("INFO headers "+headers)
        # print(re.search("ref", headers, re.IGNORECASE).group())
        if (
            var[re.search("ref", headers, re.IGNORECASE).group()].isalpha()
            and len(var[re.search("ref", headers, re.IGNORECASE).group()]) == 1
            and var[re.search("alt", headers, re.IGNORECASE).group()].isalpha()
            and len(var[re.search("alt", headers, re.IGNORECASE).group()]) == 1
        ):
            snv.append(var)
        else:
            indels.append(var)
    df_snv = pd.DataFrame(snv)
    df_indels = pd.DataFrame(indels)
    print(df_snv.head())
    print(len(df_snv.index))
    print(df_indels.head())
    print(len(df_indels.index))
    return df_snv, df_indels


def parseargs():
    parser = argparse.ArgumentParser(
        description="Filter tsv in DPNI and POOL context, basically tsv have to come from varank analysis "
    )
    subparsers = parser.add_subparsers(dest="command")

    # subparser A
    parser_stats = subparsers.add_parser(
        "stats", help="a help", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser_stats.add_argument(
        "-f",
        "--files",
        type=str,
        help="Absolute path of folder containing vcf cas index or only one vcf",
    )
    parser_stats.add_argument(
        "-p", "--pool", type=str, help="Absolute path of pool vcf or tsv"
    )
    parser_stats.add_argument("-o", "--output", type=str, help="Output Folder")
    parser_stats.add_argument(
        "-n",
        "--name",
        type=str,
        default=None,
        help="Sample name for index case",
        dest="sample_name",
    )
    parser_stats.add_argument(
        "-s",
        "--skip",
        type=str,
        help="vcf to skip in folder separate by ',', default take all vcf in folder",
    )
    parser_stats.add_argument(
        "-t", "--tsv", action="store_true", help="If file are tsv instead of vcf"
    )
    parser_stats.add_argument(
        "-m",
        "--mother",
        required=True,
        type=str,
        help="Abs path of of tsv containing variant comming from mother, they will be discarded from pool sample",
    )
    parser_stats.add_argument(
        "-c",
        "--covfile",
        type=str,
        help="Absolute path of covfile generated with gatk view depth split by bases for each variant position",
    )
    parser_stats.set_defaults(func=main_stats)
    # Subparser 2
    parser_filter = subparsers.add_parser(
        "filter", help="b help", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser_filter.add_argument(
        "-c", "--indexcase", type=str, help="Abs path of index case tsv"
    )
    parser_filter.add_argument(
        "-m", "--mother", type=str, help="Abs path of mother tsv"
    )
    parser_filter.add_argument(
        "-s",
        "--skiprows",
        type=int,
        default=2,
        help="row to skip in tsv, do not count columns name",
    )
    parser_filter.add_argument("-out", type=str, help="Output tsv")
    parser_filter.add_argument(
        "-indels",
        "--indelfilter",
        type=str,
        default="cNomen",
        help="col filter for indels comparison",
    )
    parser_filter.add_argument(
        "-snv",
        "--snvfilter",
        type=str,
        default="variantID",
        help="col filter for snv comparison",
    )
    parser_filter.set_defaults(func=main_filter)

    # subparser 3
    parser_pool = subparsers.add_parser(
        "pool", help="c help", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser_pool.add_argument(
        "-fd",
        "--foetdad",
        type=str,
        help="Abs Path of tsv containing fetal variants or dad variant at 100",
    )
    parser_pool.add_argument(
        "-p",
        "--pool",
        type=str,
        help="Abs path of tsv containing variants comming from maternal blood, or pooled sample",
    )
    parser_pool.add_argument(
        "-o", "--output", type=str, help="Folder for filtered output"
    )
    parser_pool.add_argument(
        "-c",
        "--columns",
        type=str,
        help="list of numerical column with , instead of . which need switch ex: --columns=varReadDepth,QUALphred",
        default="alleleFrequency",
    )
    parser_pool.set_defaults(func=main_pool)

    # ---------------------------
    # TEST args
    # ---------------------------
    parser_test = subparsers.add_parser(
        "test", help="c help", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser_test.add_argument(
        "-c",
        "--cov",
        type=str,
        help="gatk depth of coverage file",
        default="/home1/BAS/lamouchj/scripts/DPNI/test/FF5_depthofcov.cov",
    )
    parser_test.add_argument(
        "-m",
        "--miss",
        type=str,
        help="tsvmiss file",
        default="/home1/BAS/lamouchj/scripts/pool_script/test_5_2305_all_modif_split_snv_indels/missVariants_20220525-162907.tsv",
    )
    parser_test.add_argument(
        "-o",
        "--output",
        type=str,
        help="Folder for output",
        default="/home1/BAS/lamouchj/scripts/DPNI/test/FF5_depthofcov_miss_"
        + datetime.datetime.today().strftime("%Y%m%d-%H%M%S")
        + ".tsv",
    )
    parser_test.set_defaults(func=main_test)

    # return parser
    args = parser.parse_args()
    return args


def main_pool():
    args = parseargs()
    # replace , by . in column specified
    columns = args.columns.split(",")
    foetus = tsvtodataframe(args.foetdad, 2, columns)
    pool = tsvtodataframe(args.pool, 2, columns)
    foetus_filtered = foetus.loc[
        (foetus["varReadPercent"] > 20)
        & (foetus["varReadDepth"] > 20)
        & (foetus["alleleFrequency"] < 0.2)
    ]
    pool_filtered = pool.loc[
        (pool["varReadPercent"] > 20)
        & (pool["varReadDepth"] > 20)
        & (pool["alleleFrequency"] < 0.2)
    ]

    if not os.path.exists(args.output):
        os.mkdir(args.output)
    writeTsv(foetus_filtered, "foetus_filtered", args.output)
    writeTsv(pool_filtered, "pool_filtered", args.output)


def splitvar(df):
    if "varType" in df.columns:
        snv = df.loc[df["varType"] == "substitution"]
        indels = df.loc[df["varType"] != "substitution"]
    return snv, indels


def rework_cols(df, colid):
    # rework columns
    df.dropna(axis=1, how="all", inplace=True)
    new_col = []
    for values in df.columns:
        if values.endswith(colid):
            values = values.split("_")[:-1]
            new_col.append("_".join(values))
        else:
            new_col.append(values)
    # print(new_col)

    df.columns = new_col
    return df


def merge_query(left, right, on, how, keep):
    print("#[INFO] Params", {"on": on, "how": how, "keep": keep})
    return (
        left.merge(right, on=on, how=how, indicator="i")
        .query("i ==  '" + keep + "'")
        .drop("i", 1)
    )


def merge_df(left, right, on_snv, on_indels, how, keep):
    """
    left:dataframe to merge, right dataframe to merge, on column field, how to merge, keep row to keep default left_only
    so only values present in left or right are keept regarding args
    """
    if keep == "left_only":
        colid = "_x"
    elif keep == "right_only":
        colid = "_y"
    print("\n#[INFO] Merge step")
    left_snv, left_indels = splitvar(left)
    right_snv, right_indels = splitvar(right)
    snv = merge_query(left_snv, right_snv, on_snv, how, keep)

    finalsnv = rework_cols(snv, colid)

    indels = merge_query(left_indels, right_indels, on_indels, how, keep)
    finalindels = rework_cols(indels, colid)

    df = pd.concat([finalsnv, finalindels], ignore_index=True)
    # print(df.head())
    # print(*df.columns)
    # exit()

    df.dropna(axis=1, how="all", inplace=True)
    print("#[INFO] Left df ", len(left.index))
    print("#[INFO] Right df ", len(right.index))
    print("#[INFO] DF rows only in " + keep + " :", len(df.index))
    return df


def main_stats():
    args = parseargs()
    if not args.output:
        print("ERROR specify an output folder")
    elif not os.path.exists(args.output):
        os.mkdir(args.output)
        print("#[INFO) create output folder " + os.path.abspath(args.output))

    if args.skip:
        samples_skip = args.skip.split(",")
    else:
        samples_skip = []
    if not args.tsv:
        dfpool = vcfTodataframe(args.pool, "vcf")
        list_df = preprocessVcf(args.files, samples_skip)
        filter_snv = filter_indels = "POS"
    else:
        # Assuming tsv comes from varank
        pool = pd.read_csv(args.pool, sep="\t", skiprows=2)
        # 09/05 filter to have good false pos (popfreq and so on)
        # pool_tmp, dfpool_opposite = filter_fp(pool)
        list_df_tmp = preprocessTsv(args.files, samples_skip, args.sample_name)
        list_df = []
        for ci in list_df_tmp:
            ci_fp, ci_fp_opposite = filter_fp(ci, "cas_index_opposite", args.output)
            list_df.append(ci_fp)
        filter_snv = "variantID"
        filter_indels = "cNomen"

    # keep opposite df
    writeTsv(ci_fp_opposite, "cas_index_opposite", args.output)

    # remove variants in pool which are opposite of variant with filter in fasle_fp func 23/05, VAF < 20 AD < 20 AFbaseinterne > 2
    # cas index 100 opposite
    opposite = ci_fp_opposite
    dfpool = merge_df(pool, opposite, "variantID", "cNomen", "outer", "left_only")

    final_snv_tmp, final_indels_tmp = splitVariants(list_df)

    mother, head = read_tsv(args.mother)

    (
        snvRate,
        indelRate,
        snvfound,
        snvList,
        indelsfound,
        indelList,
        meanSNV,
        meanIndels,
        missSNV,
        missInDels,
        false_pos_snv,
        false_pos_indels,
    ) = toPlot(
        final_snv_tmp,
        final_indels_tmp,
        dfpool,
        (filter_snv,),
        (filter_indels,),
        args.output,
        mother,
        args.covfile,
    )
    dataframeMetrics(
        snvRate,
        indelRate,
        snvfound,
        snvList,
        indelsfound,
        indelList,
        meanSNV,
        meanIndels,
        missSNV,
        missInDels,
        false_pos_snv,
        false_pos_indels,
        len(dfpool.index),
        args.output,
    )


def main_filter():
    args = parseargs()
    if not args.out:
        print("ERROR no output specified")
        exit()
    else:
        output = args.out
    foetus = tsvtodataframe(args.indexcase, args.skiprows)
    mother = tsvtodataframe(args.mother, args.skiprows)
    snv_f, indels_f = splitypeofvar(foetus)
    # snv_m, indels_m = splitypeofvar(mother)
    snv = parse(snv_f, mother, args.indelfilter)
    indels = parse(indels_f, mother, args.snvfilter)
    last = pd.concat([snv, indels])
    last = last.sort_values(by="varankVarScore", ascending=False)
    last.to_csv(output, sep="\t", header=True, index=False)


def viz(df, output):
    df = df.sort_values(by="Total_Depth", ascending=False)
    print(df.head())
    print(df["Total_Depth"].head())
    fig = px.bar(
        df,
        x="variantID",
        y=["A", "C", "G", "T"],
        title="check repre in population",
        hover_name="variantID",
        hover_data=["FF5_base_counts", "ref", "alt", "Total_Depth"],
        barmode="stack",
    )


def viz_go(input, output):
    df = input.loc[input["varType"] == "substitution"]
    print(*df.columns)
    print(len(df.index))
    print(df.head())

    df.sort_values(by="Total_Depth", ascending=False, inplace=True)
    # print(df.head())

    fig = make_subplots(rows=3, cols=1, vertical_spacing=0.25)
    # test = [i for i in range(0, len(df.index) + 1, 1)]
    fig.update_layout(paper_bgcolor="#E3E3E3")
    # REMINDER CHANGE TYPE COL FROM OBJECT TO INT FOR PLOTLY
    # fig.add_trace(go.Bar(x=df["variantID"], y=df["A"].astype(int), name="test2"))

    for col in ["A", "C", "G", "T"]:
        fig.add_trace(
            go.Bar(x=df["variantID"], y=df[col].astype(int), name=col, legendgroup="1"),
            row=1,
            col=1,
        )
    # discard_base_interne = df.loc[df["alleleFrequency"] < 0.02]
    custom = df[
        ["varReadPercent", "varReadDepth", "varankVarScore", "alleleFrequency"]
    ].copy()
    custom["varReadPercent"] = custom["varReadPercent"] / 100
    print("#[INFO] HOVER scatter plot ", len(custom.index))
    # print(df[["varReadPercent", "varReadDepth", "varankVarScore", "alleleFrequency"]].head())
    # customdata=np.dstack((df["varReadPercent"] / 100, df["varReadDepth"], df["varankVarScore"], df["alleleFrequency"]))

    fig.add_trace(
        go.Scatter(
            x=(df["varReadPercent"] / 100).round(2),
            y=df["VAF_miss"],
            marker=dict(colorbar=dict(title="VarankScore", len=0.4), colorscale="Jet"),
            mode="markers",
            name="VAF foetus / maternal blood",
            text=df["variantID"],
            # customdata=np.dstack((df["varReadPercent"] / 100, df["varReadDepth"], df["varankVarScore"], df["alleleFrequency"])),
            customdata=custom,
            hovertemplate="<b>%{text}</b><br>varReadPercent: %{customdata[0]} <br>vRP maternalblood: %{y} <br>varReadDepth: %{customdata[1]} <br>varankVarScore: %{customdata[2]} <br>FrequencyCohorte: %{customdata[3]:.2f}",
            marker_color=df["varankVarScore"],
            legendgroup="2",
        ),
        row=2,
        col=1,
    )
    fig.add_trace(
        go.Histogram(
            histfunc="count", x=df["Explanations"], name="count", legendgroup="3"
        ),
        row=3,
        col=1,
    )
    fig.update_yaxes(title_text="Total depth split by bases", row=1, col=1)
    fig.update_xaxes(title_text="Variant ID", row=1, col=1)
    fig.update_yaxes(title_text="AlleleFrequency in maternal blood", row=2, col=1)
    fig.update_xaxes(
        title_text="AlleleFrequency in foetus invasive prelevement", row=2, col=1
    )
    fig.update_yaxes(title_text="Number of variants", row=3, col=1)
    fig.update_xaxes(title_text="Explanations", row=3, col=1)
    fig.update_layout(
        height=1000,
        width=1000,
        title="Missings variants exploration SNV",
        barmode="stack",
    )
    fig.write_html(output)
    # fig.show()


def identify_miss(df):
    dico = {"Explanations": []}
    # search why variants has not been called in maternalblood sample
    for i, variants in df.iterrows():
        range_var = []
        for j in range(1, 3, 1):
            range_var.append(variants["start"] - j)
        # case of variants nearby
        # if df["start"].isin(range_var).any():
        #    test = ",".join(df[df["start"].isin(range_var)]["variantID"].to_list())
        #    dico["Explanations"].append("Annotations issues " + str(test) + " nearby")
        # GC-rich or repeat zone arbitrary set at VAF allele muté 5%
        if variants["VAF_miss"] > 0.05:
            dico["Explanations"].append("GC-rich or repeat zone")
        # no explanation
        else:
            dico["Explanations"].append("Dilution effect")
    expla = pd.DataFrame.from_dict(dico)
    print(expla["Explanations"].unique())
    last = pd.concat([df, expla], axis=1)
    return last


def main_test():
    args = parseargs()
    df = add_cov(args.cov, args.miss)
    viz_go(df, osj(os.path.dirname(args.output), "res_plot.html"))


def main():
    args = parseargs()
    if args.command == "stats":
        main_stats()
    elif args.command == "filter":
        main_filter()
    elif args.command == "pool":
        main_pool()
    else:
        main_test()


if __name__ == "__main__":
    main()
