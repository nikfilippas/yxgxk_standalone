from make_yml import make_yml

run_name = 'gyksrA_T08_new'
fname_data = 'data/saccfiles/cls_cov_new.fits'

kmax_arr = [0.5, 1., 1., 1., 1., 1.]

for i in range(6):
    tname = f'LOWZ__{i}'
    make_yml(params_vary=[f'ygk_{tname}_lMmin_0',
                          f'ygk_{tname}_lM1_0',
                          'ygk_rhogy',
                          'ygk_mass_bias',
                          'sigma8',
                          'ygk_Ahmc'],
             corrs=[(tname, tname),
                    (tname, 'YMILCA'),
                    (tname, 'KAPPA')],
             bias_model='HaloModel',
             lmin=2,
             kmax=kmax_arr[i],
             mass_function="Tinker08",
             concentration="Ishiyama21",
             hm_correction="halofit",
             ns_independent=False,
             fname_data=fname_data,
             dirname_out=f'chains/{run_name}/{run_name}_{i}',
             sampler='mcmc', nsamples=30000, debug=False)
