import pandas as pd
import os
from os.path import join as osj
import sys

sys.path.append("..")
from utils.utils import systemcall
import numpy as np


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
            return self.launch_seqff(self.bamfoetus, self.output, "sam")
        elif self.bamfoetus.endswith(".tsv"):
            return self.launch_seqff(self.bamfoetus, self.output, "counts")
        else:
            bf = self.readprofile()
            return self.launch_seqff(bf, self.output, "counts")

    def launch_seqff(self, foetus, output, format):
        name = os.path.basename(foetus).split(".", 1)[0] + ".tsv"
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
        df = pd.read_csv(osj(output, name), sep=",", header=0)
        df.columns = ["Method", "Foetal fraction"]
        df["Variants count"] = np.nan
        return df[["Method", "Variants count", "Foetal fraction"]]

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
