# Report the status of all chains in the specified directory.
import os
import yaml
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("path")
args = parser.parse_args()

for model in os.listdir(args.path):
    if model.startswith("."): continue  # don't peek hidden files
    print(model)
    path1 = os.path.join(args.path, model)
    if model == "gyksrA": continue  # this is a symlink - will replace later
    for model_i in os.listdir(path1):
        print(f" {model_i}")
        path2 = os.path.join(path1, model_i)
        if not "cobaya.checkpoint" in os.listdir(path2):
            print("    no chains found")
            continue
        path3 = os.path.join(path2, "cobaya.checkpoint")
        with open(path3) as f:
            info = yaml.safe_load(f)
            print(f"    converged: {info['sampler']['mcmc']['converged']}")


"""
# List of unconverged chains.

gyksrA_0

gksA_0
gksA_1

gykrA_0
gykrA_1
gykrA_2
gykrA_3
gykrA_4
gykrA_5

gyksr_0
gyksr_1
gyksrAAA_0
gyksrAAA_1

gyksrA_SZ_0
gyksrA_T10_0

gyksrA_lmin40_1
gyksrA_kmax025_0
gyksrA_kmax025_1
gyksrA_kmax025_2
"""
