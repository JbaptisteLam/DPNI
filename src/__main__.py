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
from utils.utils import fancystdout
from parseargs.parseargs import parseargs
import sys
import os


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
        if not "seqff" in args:
            args.seqff = "estimation/seqff/pd4615-sup-0002-seqff.r"
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

    # elif args.func == "plotvaf":
    #    scatter_vaf(args.dataframe, args.output, args.name)

    # getUMI("/home1/data/STARK/data/DPNI/trio/TWIST/FCL2104751.bwamem.bam", pos)
    # FF = VAFp + (1 - VAFm)
    # print("#[INFO] Estimation of Foetal fraction : ", FF)


if __name__ == "__main__":
    main()
