import pandas as pd
from os.path import join as osj
import os
import json
import sys

sys.path.append("..")
from utils.utils import vcfTodataframe, tsvTodataframe


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
