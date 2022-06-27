import argparse


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
        default="bedtools",
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
