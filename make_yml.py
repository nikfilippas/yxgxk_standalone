import yaml
import os

def make_yml(params_vary, corrs, bias_model, lmin, kmax,
             mass_function, concentration, hm_correction, ns_independent,
             fname_data, dirname_out, sampler,
             debug=True, nsamples=10000):
    samplers = {
        'minimize': {'ignore_prior': True, 'max_evals': nsamples},
        'mcmc': {'learn_proposal': True, 'burn_in': 10, 'max_samples': nsamples}}

    # Open master file
    with open('master.yml', "r") as f:
        dflt = yaml.load(f, Loader=yaml.FullLoader)

    # Update parameters
    pars = {k: d['ref']['loc'] for k, d in dflt['params'].items()}
    for p in params_vary:
        if p not in pars:
            raise ValueError(f"Unknown param {p}")
        pars[p] = dflt['params'][p]
    # Update two-point functions
    twopts = [{'bins': [t1, t2]} for t1, t2 in corrs]
    # Create dictionary
    dout = dflt.copy()
    dout['params'] = pars
    dout['likelihood']['ygk_like.ygkLike']['twopoints'] = twopts
    dout['likelihood']['ygk_like.ygkLike']['bz_model'] = bias_model
    dout['likelihood']['ygk_like.ygkLike']['mf_name'] = mass_function
    dout['likelihood']['ygk_like.ygkLike']['cm_name'] = concentration
    dout['likelihood']['ygk_like.ygkLike']['HM_correction'] = hm_correction
    dic = dout["likelihood"]["ygk_like.ygkLike"]["defaults"]
    for item in dic:
        if item.startswith("LOWZ"):
            dic[item]["lmin"] = lmin
    dout['likelihood']['ygk_like.ygkLike']['defaults']['lmin'] = lmin
    dout['likelihood']['ygk_like.ygkLike']['defaults']['kmax'] = kmax
    dout['likelihood']['ygk_like.ygkLike']['ns_independent'] = ns_independent
    dout['likelihood']['ygk_like.ygkLike']['input_file'] = fname_data
    dout['debug'] = debug
    dout['output'] = f'{dirname_out}/cobaya'
    dout['sampler'] = {sampler: samplers[sampler]}

    # Save output
    os.system(f'mkdir -p {dirname_out}')
    fname_out = f'{dirname_out}/params.yml'
    with open(fname_out, 'w') as outfile:
        yaml.dump(dout, outfile, default_flow_style=False)
