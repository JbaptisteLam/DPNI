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

    print("#[INFO] " + file)
    if extension == ".vcf":
        # print('[#INFO] VCF: '+file)
        with open(file) as f:
            for lines in f:
                if lines.startswith("#"):
                    header.append(lines.strip())
                else:
                    variants_tmp.append(lines)

    print("#[INFO]", header[-1])
    # print(variants[-1])

    col = header[-1].strip().split("\t")
    for v in variants_tmp:
        variants.append(v.strip().split("\t"))

    # headerCol = [res.replace('#', '') for res in colTemp]
    dfVar = pd.DataFrame(columns=col)
    # print(variants[0:3])

    # Creating Dataframe from the whole VCF
    print("#[INFO] Whole VCF to Dataframe")
    for i, var in enumerate(variants):
        print_progress_bar(i, len(variants))
        rows = pd.Series(var, index=dfVar.columns)
        # print(rows[0])
        dfVar = dfVar.append(rows, ignore_index=True)

    print("\n")
    if rheader:
        return dfVar, header
    else:
        return dfVar


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

            df_tmp = vcfTodataframe(vcf)
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

    # if len(df_miss) != 0:
    # 	#if some variants are missing
    # 	print("#[INFO] dfmiss", df_miss.head(), len(df_miss.index))
    # 	dfound = dfci[~dfci[col[0]].isin(miss)]
    # 	dfound.loc[:, "DP"] = pd.to_numeric(dfound.loc[:, "DP"], downcast="float")
    # else:
    # 	dfound = pd.DataFrame(dfci, columns=dfci.columns)
    #
    # if ci:
    # 	mean = dfound['DP'].to_numpy().mean()

    # meanDepth

    # print("INFO DP")
    ## print(dfound['DP'].head())
    # print("MEANDEPTH", mean)
    # print(df_miss.head())
    # print("LENGTH MISS", len(df_miss.index))
    # print("\n\n")
    # print(dfound.head())
    # print("LENGTH FOUND", len(dfound.index))

    # time.sleep(10)
    # print("#[INFO] length of the miss variants", len(df.index))
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


# ef toPlot(dfci_snv, dfci_indels, dfpool, filt_snv, filt_indels, output):
# 	'''
# 	input: a list of dataframe of SNV from 1 allele to the y patients x2 alleles ex: 12 patients from 1 to 24 dataframe,
# 	1 dataframe containing pool of all the sample in the cohorte
# 	'''
# 	snv_list_percent = []
# 	others_list_percent = []
# 	total_snv_list = []
# 	total_others_list = []
# 	snv_list_miss = []
# 	indels_list_miss = []
# 	found_snv = []
# 	found_indels = []
# 	meanSNV = []
# 	meanInDels = []
# 	print("#[INFO] Stats SNV ...")
# 	for i, variants in enumerate(dfci_snv):
# 		#print_progress_bar(i, len(dfci_snv))
# 		snv, total, dfoundSNV, snvmiss, meanS = getSensi(variants, dfpool, filt_snv)
# 		snv_list_percent.append(snv)
# 		total_snv_list.append(total)
# 		snv_list_miss.append(snvmiss)
# 		meanSNV.append(meanS)
# 		found_snv.append(dfoundSNV)
#
# 	print("\n")
# 	print("#[INFO] Stats InDels ...")
# 	for j, var in enumerate(dfci_indels):
# 		#print_progress_bar(j, len(dfci_indels) - 1)
# 		others, total, dfoundInDels, indelsmiss, meanI = getSensi(var, dfpool, filt_indels)
# 		others_list_percent.append(others)
# 		total_others_list.append(total)
# 		indels_list_miss.append(indelsmiss)
# 		meanInDels.append(meanI)
# 		found_indels.append(dfoundInDels)
#
# 	print("\n")
#
# 	#merge all variants miss
# 	dfsnv = []
# 	dfindels = []
# 	for df in snv_list_miss:
# 		df = df.reset_index()
# 		if 'index' in list(df.columns):
# 			df = df.drop(columns=['index'])
# 		dfsnv.append(df)
# 	for df in indels_list_miss:
# 		df = df.reset_index()
# 		if 'index' in list(df.columns):
# 			df = df.drop(columns=['index'])
# 		dfindels.append(df)
#
# 	#Allelecountcohorte
# 	if output:
# 		df_miss_snv = pd.concat(dfsnv)
# 		df_miss_indels = pd.concat(dfindels)
# 		df_miss = pd.concat([df_miss_snv,df_miss_indels])
# 		if 'level_0' in list(df_miss.columns):
# 			df_miss = df_miss.drop(columns=['level_0'])
#
# 		df_snv = pd.concat(found_snv)
# 		df_indels = pd.concat(found_indels)
# 		if '#CHROM' in found_snv[0].columns:
# 			col = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'AlleleCountCohorte', 'SampleID' ]
# 			df_miss = df_miss.loc[:, col]
# 			df_snv = df_snv.loc[:, col]
# 			df_indels = df_indels.loc[:, col]
# 		writeTsv(df_miss, 'missVariants', output)
# 		writeTsv(df_snv, 'SNVfound', output)
# 		writeTsv(df_indels, 'InDelsfound', output)
#
# 		#Faux pos pool
# 		positif = pd.concat([df_snv, df_indels, df_miss], ignore_index=True)
# 		print("#[INFO] positif ", len(positif.index))
# 		positif.drop_duplicates(inplace=True)
# 		print("#[INFO] after removing duplicate positif ", len(positif.index))
# 		print("#[INFO] Length POOL ", len(dfpool.index))
#
#
# 		false_pos = dfpool.merge(positif, indicator='i', on='variantID', how='outer').query('i == "left_only"').drop('i', 1)
# 		false_pos.dropna(axis=1, how='all', inplace=True)
# 		#remove _x
# 		new_col = []
# 		for values in df.columns:
# 			if values.endswith('_x'):
# 				values = values.split('_')[0]
# 				new_col.append(values)
# 			else:
# 				new_col.append(values)
# 		false_pos.columns = new_col
# 		print("#[INFO] Length falsepos ", len(false_pos.index))
#
# 		#Split snv indels
# 		false_pos_snv = false_pos.loc[false_pos['varType'] == 'substitution']
# 		false_pos_indels = false_pos.loc[false_pos['varType'] != 'substitution']
# 		writeTsv(false_pos_snv, 'FalsePositive_SNV', output)
# 		writeTsv(false_pos_indels, 'FalsePositive_InDels', output)
#
# 		#Metrics to return then present in html page
#
#
#
#
# 		#print(snv_list_percent,others_list_percent, total_snv_list, total_others_list, len(snv_list_percent))
#
# 	#uncomment to plot TODO
# 	#plotSensibilityHisto(snv_list_percent,others_list_percent, total_snv_list, total_others_list, len(snv_list_percent), output)
# 	return snv_list_percent, others_list_percent, total_snv_list, total_others_list, meanSNV, meanInDels, snv_list_miss, indels_list_miss

# 09/05/2022   remove variant s which are not in this filter and
def filter_fp(df):
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
    return df_filter, df_opposite


def toPlot(dfci_snv, dfci_indels, dfpool, filt_snv, filt_indels, output, dfmother):
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

        writeTsv(df_miss, "missVariants", output)
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
        dfmother_proper, dfmother_opposite = filter_fp(dfmother)
        # print(dfmother_opposite.head())
        # print(*dfmother_proper.columns)
        variants_mother_only = merge_df(
            total_ci, dfmother_proper, "variantID", "cNomen", "outer", "right_only"
        )
        #
        dfpool_filter = merge_df(
            dfpool, variants_mother_only, "variantID", "cNomen", "outer", "left_only"
        )

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
    if skiprows:
        df = pd.read_csv(file, skiprows=skiprows, sep="\t", header=0)
    else:
        df = pd.read_csv(file, sep="\t", header=0)
    if "index" in df:
        df = df.drop(columns="index")
    for col in columns:
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
    so only values present in left are keep
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
        dfpool = vcfTodataframe(args.pool)
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
            ci_fp, ci_fp_opposite = filter_fp(ci)
            list_df.append(ci_fp)
        filter_snv = "variantID"
        filter_indels = "cNomen"

    # remove variatns in pool which are opposite of variant with filter in fasle_fp func 23/05, VAF < 20 AD < 20 AFpop > 2
    # cas index 100 opposite
    opposite = ci_fp_opposite
    dfpool = merge_df(pool, opposite, "variantID", "cNomen", "outer", "left_only")
    # print(dfpool.head())
    # print(*dfpool.columns)
    # exit()

    final_snv_tmp, final_indels_tmp = splitVariants(list_df)

    # snv_add = groupBy(snv_tmp, filter_snv)
    # indels_add = groupBy(indels_tmp, filter_indels)
    #
    # final_snv_tmp = createFinaldf(snv_tmp, snv_add, filter_snv)
    # final_indels_tmp = createFinaldf(indels_tmp, indels_add, filter_indels)

    # final_snv, alleleSNV = splitVcfAllele(final_snv_tmp, len(list_df))
    # final_indels, alleleInDels = splitVcfAllele(final_indels_tmp, len(list_df))
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


def test():
    args = parseargs()
    for vcf in os.listdir(args.folder):
        df = vcfTodataframe(os.path.join(args.folder, vcf))
        sample, bad = parseSampleField(df)
        print(sample.head())
        # print(tabulate(sample.loc[:, ['#CHROM', 'POS', 'GT']].head(), headers='keys', tablefmt = 'psql'))
        # print(bad.head())
        time.sleep(30)
    return


def main():
    args = parseargs()
    if args.command == "stats":
        main_stats()
    elif args.command == "filter":
        main_filter()
    elif args.command == "pool":
        main_pool()


if __name__ == "__main__":
    main()
