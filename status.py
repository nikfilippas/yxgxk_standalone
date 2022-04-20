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

gyksrA_0  # running (wrong dir: chains/gyksrA/gyksrA_0)
gyksrA_1  # running

gyksrA_bG073_0  # running

gksA_0  # running
gksA_1  # running

gykrA_0  # running

gyksr_0  # running
gyksr_1  # running

gyksrAAA_0  # running
gyksrAAA_1  # running
gyksrAAA_2  # running
gyksrAAA_3  # running

gyksrA_SZ_0  # running
gyksrA_SZ_1  # running

gyksrA_T10_0  # running
gyksrA_T10_1  # running

gyksrA_D16_0  # running
gyksrA_D16_1  # running

gyksrA_B16_0
gyksrA_B16_1

gyksrA_lmin40_0  # running
gyksrA_lmin40_1  # running
gyksrA_lmin40_2  # running

gyksrA_kmax025_0
gyksrA_kmax025_1
gyksrA_kmax025_2
gyksrA_kmax025_3
gyksrA_kmax025_4
gyksrA_kmax025_5

gyksrAw_0  # running
gyksrAw_1  # running
gyksrAw_2  # running
gyksrAw_3  # running
gyksrAw_4  # running
gyksrAw_5  # running

gyksrA_bU075
"""
