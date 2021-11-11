#!/bin/bash

run_name=$1

for zbin in {0..5}
do
  addqueue -s -n 1x12 -q cmb -m 0.5 -c "$run_name"_"$zbin" /users/nikfilippas/.local/bin/cobaya-run chains/$run_name/"$run_name"_"$zbin"/params.yml
done
