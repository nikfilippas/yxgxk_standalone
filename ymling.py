from make_yml import make_yml

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
             kmax=1.0,
             fname_data=fname_data,
             dirname_out=f'yxgxk_b_uniform_{i}',
             sampler='mcmc', nsamples=30000, debug=False)
