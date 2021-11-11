#!/bin/bash

run_name=$1

for zbin in {0..5}
do
  addqueue -s -n 1x12 -q cmb -m 0.5 -c $run_name_$z_bin /users/nikfilippas/.local/bin/cobaya-run chains/$run_name/$run_name_$zbin/params.yml
done
