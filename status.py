# Report the status of all chains in the specified directory.
import os
import yaml
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("path")
args = parser.parse_args()

for model in os.listdir(args.path):
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
