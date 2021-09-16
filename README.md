A few instructions

- `master.yml` contains a yaml file with all possible free parameters of the model.
These are the usual cosmological parameters, the HOD parameters, mass_bias, rho_gy/gm,
Ahm_gy/gm/gg, alpha_gy/gm/gg, and the shift and width parameters for each redshift bins.
I'd suggest starting only with a minimal subset of these (HOD+rho_gy+mass_bias and sigma8
whenever you have lensing data).
- You'll see that there are a few other settings (mass function type, whether M0 tracks Mmin etc.).
I wondln't worry about these for now.
- Note that the different redshift bins are now called `"LOWZ__X"` where X goes from 0 to 5.
`"LOWZ__0"` is 2MPZ, and the rest are the WIxSC bins.
- The script `ymling.py` makes use of `master.yml` and a function defined in `make_yml.py`
to generate the yaml file needed to launch the cobaya chains for a given combination of
free parameter and set of power spectra. I've left only an example that shows how to create
 the yaml files for the minimal yxgxk runs in each redshift bin.
- After generating those yml files, this is the command I've been using to lauch chains on
glamdring (e.g. for the 4th redshift bin:
```
addqueue -c yxgxk_4 -s -n 1x12 -q cmb -m 0.5 /users/damonge/.local/bin/cobaya-run yxgxk_4/params.yml
```
If you need to run this on several nodes (although I didn't have to for these runs), change -s for -O,
and select the right number of nodes and cores per node after -n.
- The directory `yxgxk_2` contains an example of the result of running `ymling.py` and then running
the chain with cobaya-run.
- I've put the sacc file in `/mnt/extraspace/damonge/yxgxk/cls_cov.fits`.
- The likelihood code lives in `yxgxk_like`. This contains two python files. One is just a CCL wrapper.
 The real meat is in `yxgxk_like.py`. You'll need to install the likelihood (`python3 setup.py install --user`)
 before being able to run `cobaya-run`. Remember to do this again every time the likelihood code changes.
- `plot_triangles.py` shows how to plot the results of a cobaya chain.
