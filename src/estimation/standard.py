import pandas as pd
from .preset import Process
from os.path import join as osj
import os
import sys


sys.path.append("..")
from utils.utils import average, plotvaf, series_to_stats, scatter_vaf


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
        # try jinja template test
        fig, js = scatter_vaf(paternal_var_rare, self.output, "vafFFestim", dico)

        print(paternal_var_rare["varReadPercent"].head())
        paternal_var_rare.to_csv(
            osj(self.output, "paternalcheck.tsv"), header=True, index=False, sep="\t"
        )
        return paternal_var_rare, dico, ff, js
