import pandas as pd
from .preset import Process
from os.path import join as osj
import os
import sys
import json
import pysam

sys.path.append("..")
from utils.utils import systemcall, average, col_type_float, systemcall


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
