#!/bin/bash

run_name=$1

for zbin in {0..5}
do
  addqueue -n 8 -q cmb -m 0.5 -c "$run_name"_"$zbin" -o "$run_name"_"$zbin".out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/$run_name/"$run_name"_"$zbin"/params.yml
done

# *** EXAMPLE ***
# addqueue -n 8 -q cmb -m 0.5 -c gyksrA_0 -o gyksrA_0.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gyksrA/gyksrA_0/params.yml
