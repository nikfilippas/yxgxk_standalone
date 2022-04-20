#!/bin/bash

run_name=$1

for zbin in {0..5}
do
  addqueue -n 8 -q cmb -m 0.5 -c "$run_name"_"$zbin" -o "$run_name"_"$zbin".out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/$run_name/"$run_name"_"$zbin"/params.yml
done

# *** EXAMPLE ***
"""
addqueue -n 8 -q cmb -m 1.5 -c gyksrA_0 -o gyksrA_0.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gyksrA/gyksrA_0/params.yml
addqueue -n 8 -q cmb -m 1.5 -c gyksrA_1 -o gyksrA_1.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gyksrA/gyksrA_1/params.yml
addqueue -n 8 -q cmb -m 1.5 -c gyksrA_bG073_0 -o gyksrA_bG073_0.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gyksrA_bG073/gyksrA_bG073_0/params.yml
addqueue -n 8 -q cmb -m 1.5 -c gksA_0 -o gksA_0.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gksA/gksA_0/params.yml
addqueue -n 8 -q cmb -m 1.5 -c gksA_1 -o gksA_1.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gksA/gksA_1/params.yml
addqueue -n 8 -q cmb -m 1.5 -c gykrA_0 -o gykrA_0.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gykrA/gykrA_0/params.yml
addqueue -n 8 -q cmb -m 1.5 -c gyksr_0 -o gyksr_0.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gyksr/gyksr_0/params.yml
addqueue -n 8 -q cmb -m 1.5 -c gyksr_1 -o gyksr_1.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gyksr/gyksr_1/params.yml
addqueue -n 8 -q cmb -m 1.5 -c gyksrAAA_0 -o gyksrAAA_0.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gyksrAAA/gyksrAAA_0/params.yml
addqueue -n 8 -q cmb -m 1.5 -c gyksrAAA_1 -o gyksrAAA_1.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gyksrAAA/gyksrAAA_1/params.yml
addqueue -n 8 -q cmb -m 1.5 -c gyksrAAA_2 -o gyksrAAA_2.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gyksrAAA/gyksrAAA_2/params.yml
addqueue -n 8 -q cmb -m 1.5 -c gyksrAAA_3 -o gyksrAAA_3.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gyksrAAA/gyksrAAA_3/params.yml
addqueue -n 8 -q cmb -m 1.5 -c gyksrA_SZ_0 -o gyksrA_SZ_0.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gyksrA_SZ/gyksrA_SZ_0/params.yml
addqueue -n 8 -q cmb -m 1.5 -c gyksrA_SZ_1 -o gyksrA_SZ_1.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gyksrA_SZ/gyksrA_SZ_1/params.yml
addqueue -n 8 -q cmb -m 1.5 -c gyksrA_T10_0 -o gyksrA_T10_0.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gyksrA_T10/gyksrA_T10_0/params.yml
addqueue -n 8 -q cmb -m 1.5 -c gyksrA_T10_1 -o gyksrA_T10_1.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gyksrA_T10/gyksrA_T10_1/params.yml
addqueue -n 8 -q cmb -m 1.5 -c gyksrA_lmin40_0 -o gyksrA_lmin40_0.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gyksrA_lmin40/gyksrA_lmin40_0/params.yml
addqueue -n 8 -q cmb -m 1.5 -c gyksrA_lmin40_1 -o gyksrA_lmin40_1.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gyksrA_lmin40/gyksrA_lmin40_1/params.yml
addqueue -n 8 -q cmb -m 1.5 -c gyksrA_lmin40_2 -o gyksrA_lmin40_2.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gyksrA_lmin40/gyksrA_lmin40_2/params.yml
addqueue -n 8 -q berg -m 1.5 -c gyksrAw_0 -o gyksrAw_0.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gyksrAw/gyksrAw_0/params.yml
addqueue -n 8 -q berg -m 1.5 -c gyksrAw_1 -o gyksrAw_1.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gyksrAw/gyksrAw_1/params.yml
addqueue -n 8 -q berg -m 1.5 -c gyksrAw_2 -o gyksrAw_2.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gyksrAw/gyksrAw_2/params.yml
addqueue -n 8 -q berg -m 1.5 -c gyksrAw_3 -o gyksrAw_3.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gyksrAw/gyksrAw_3/params.yml
addqueue -n 8 -q berg -m 1.5 -c gyksrAw_4 -o gyksrAw_4.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gyksrAw/gyksrAw_4/params.yml
addqueue -n 8 -q berg -m 1.5 -c gyksrAw_5 -o gyksrAw_5.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gyksrAw/gyksrAw_5/params.yml
addqueue -n 8 -q cmb -m 1.5 -c gyksrA_D16_0 -o gyksrA_D16_0.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gyksrA_D16/gyksrA_D16_0/params.yml
addqueue -n 8 -q cmb -m 1.5 -c gyksrA_D16_1 -o gyksrA_D16_1.out /users/nikfilippas/.local/bin/cobaya-run chains/chains_new/gyksrA_D16/gyksrA_D16_1/params.yml
"""
