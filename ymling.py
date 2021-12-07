from make_yml import make_yml

run_name = 'yxgxksig_kmax05'
fname_data = 'cls_cov.fits'

for i in range(6):
    tname = f'LOWZ__{i}'
    make_yml(params_vary=[f'ygk_{tname}_lMmin_0',
                          f'ygk_{tname}_lM1_0',
                          'ygk_rhogy',
                          'ygk_mass_bias',
                          'sigma8'],
             corrs=[(tname, tname),
                    (tname, 'YMILCA'),
                    (tname, 'KAPPA')],
             bias_model='HaloModel',
             kmax=0.5,
             mass_function="Tinker08",
             hm_correction="halofit",
             ns_independent=False,
             fname_data=fname_data,
             dirname_out=f'chains/{run_name}/{run_name}_{i}',
             sampler='mcmc', nsamples=30000, debug=False)
