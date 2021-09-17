#!/bin/bash
declare -a arr=("2mpz" "wisc1" "wisc2" "wisc3" "wisc4" "wisc5")

for i in "${arr[@]}"
do
  addqueue -s -n 1x12 -q cmb -m 0.5 /users/nikfilippas/.local/bin/cobaya-run yxg_0/params.yml
  addqueue -s -n 1x12 -q cmb -m 0.5 /users/nikfilippas/.local/bin/cobaya-run yxg_1/params.yml
  addqueue -s -n 1x12 -q cmb -m 0.5 /users/nikfilippas/.local/bin/cobaya-run yxg_2/params.yml
  addqueue -s -n 1x12 -q cmb -m 0.5 /users/nikfilippas/.local/bin/cobaya-run yxg_3/params.yml
  addqueue -s -n 1x12 -q cmb -m 0.5 /users/nikfilippas/.local/bin/cobaya-run yxg_4/params.yml
  addqueue -s -n 1x12 -q cmb -m 0.5 /users/nikfilippas/.local/bin/cobaya-run yxg_5/params.yml
done
