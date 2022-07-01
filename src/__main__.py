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

from estimation.seqff import Seqff
from estimation.standard import Paternalidentification
from estimation.hombased import Homozygotebased
from utils.utils import fancystdout, fill_metrics, generate_report
from parseargs.parseargs import parseargs
from os.path import join as osj
import sys
import os
import json
import pandas as pd
import plotly.utils as pu
import plotly.io


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
        # 1) From CDC of DPNI studyff
        pi = Paternalidentification(
            args.mum, args.dad, args.foetus, args.type, args.output, args.quality
        )
        metrics_file = osj(args.output, "metrics.html")
        metrics_csv = osj(args.output, "metrics.tsv")
        paternal_var_rare, dico, ff, figjs = pi.main_paternal()
        dico_val = {
            "Method": [
                "mean of allele frequency of heterozygous variants comming from father"
            ],
            "Variants count": [len(paternal_var_rare.index)],
            "Foetal fraction": [ff],
        }
        df = fill_metrics(metrics_csv, dico_val)
        # df.to_pickle(osj(args.output, "metrics.pickle"))
        # generate_report(metrics_file, df, figjs)

    # 2) From combo_FF seqFF modele
    elif args.func == "seqff":
        print("#[INFO] FF estimation by seqff !")
        if not args.seqff:
            args.seqff = "/app/src/estimation/seqff/pd4615-sup-0002-seqff.r"
        pseq = Seqff(args.bfoetus, args.output, args.samtools, args.seqff)
        metrics_csv = osj(args.output, "metrics.tsv")
        dico_val = pseq.analysis()
        print(dico_val)
        fill_metrics(metrics_csv, dico_val)

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
        ff_zhang, vafp, vafm = p.get_ff()
        ff_jhiang, vars = p.get_ff_jiang()
        out = "FF Zhang et al. :" + str(ff_zhang)
        out2 = "FF Jiang et al. :" + str(ff_zhang)
        metrics_file = osj(args.output, "metrics.html")
        metrics_csv = osj(args.output, "metrics.tsv")
        dico_val = {
            "Method": [
                "FF = 1 - mean(VAFm) + mean(VAFp)",
                "FFi = 2p / p + q\n q as reference reads and p as alternate reads",
            ],
            "Variants count": [len(vafp.index) + len(vafm.index), len(vars.index)],
            "Foetal fraction": [ff_zhang, ff_jhiang],
        }
        fill_metrics(metrics_csv, dico_val)

    # elif args.func == "plotvaf":
    #    scatter_vaf(args.dataframe, args.output, args.name)

    # getUMI("/home1/data/STARK/data/DPNI/trio/TWIST/FCL2104751.bwamem.bam", pos)
    # FF = VAFp + (1 - VAFm)
    # print("#[INFO] Estimation of Foetal fraction : ", FF)
    elif args.func == "metrics":
        with open(args.js) as json_file:
            datatmp = json.load(json_file)
            data = json.dumps(datatmp, cls=pu.PlotlyJSONEncoder)
        df = pd.read_csv(args.df, sep="\t", header=0)
        print(df)
        print(type(plotly.io.from_json(data)))
        generate_report(args.output, df, data)


if __name__ == "__main__":
    main()
