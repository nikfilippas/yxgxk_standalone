import copy
from tqdm import tqdm
import yaml
import numpy as np
from cobaya.model import get_model
import getdist.mcsamples as gmc

from utils import Container
import _names


class ChainCalculator(Container):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.chains = {}
        self.summary = {}
        self.best_fits = {}
        self._zarr = np.asarray(list(self.zmid.values()))

    def _load_chains(self, model, ibin):
        # TODO: This could be merged using `self.get_sacc_file`.
        fsetup = f"{self._chains_dir}/{model}/{model}_{ibin}/cobaya.input.yaml"
        with open(fsetup, "r") as f:
            config = yaml.safe_load(f)
        fchains = f"{self._chains_dir}/{model}/{model}_{ibin}/cobaya"
        samples = gmc.loadMCSamples(fchains, settings={'ignore_rows': 0.3})
        pars = samples.getParams()
        return config, samples, pars

    def _get_derived_chain(self, config, pars, par, tracer):
        mass_function = config["likelihood"]["ygk_like.ygkLike"]["mf_name"]
        rel = self.interps.dic[par][mass_function][tracer]
        if hasattr(pars, "sigma8"):  # free sigma8
            return (
                np.array([rel(s8, by) for s8, by
                          in zip(pars.sigma8, pars.ygk_mass_bias)]).squeeze())
        else:  # fixed sigma8
            return np.array([rel(self.cosmo["sigma8"], by)
                             for by in pars.ygk_mass_bias]).squeeze()

    def _get_summary(self, samples, par):
        """Get 0.68 watershed posterior boundaries."""
        dens = samples.get1DDensity(par)
        vmin, vmax = dens.getLimits(0.68)[:2]
        vbf = np.average(dens.getLimits(0.01)[:2])
        return np.asarray([vbf, vbf-vmin, vmax-vbf])  # med, min, max

    def _calculate_summary(self, model, par):
        if not model in self.summary:
            # create empty entry if it's not there
            self.summary[model] = {}

        bf_arr = np.zeros((len(self.tracers), 3))  # holding median, min, max
        for ibin, tracer in enumerate(tqdm(self.tracers)):
            try:
                config, samples, pars = self._load_chains(model, ibin)
            except FileNotFoundError:
                bf_arr[ibin] = np.nan
                continue
            if par in self.interps.parameters:
                der_chain = self._get_derived_chain(config, pars, par, tracer)
                samples.addDerived(
                    der_chain, name=par, label=_names.latex_names[par])
            bf_arr[ibin] = self._get_summary(samples, par)

        # normalize by the growth factor for sigma8
        if par == "sigma8":
            a_arr = np.asarray(list(1/(1 + self._zarr)))
            bf_arr *= self.cosmo.growth_factor(a_arr)[:, None]

        self.summary[model][par] = bf_arr

    def _get_chains(self, model):
        if model in self.chains:
            return

        out = {}
        for ibin, tracer in enumerate(self.tracers):
            fname = f"{self._chains_dir}/{model}/{model}_{ibin}/cobaya"
            s = gmc.loadMCSamples(fname, settings={'ignore_rows': 0.3})
            p = s.getParams()
            out[tracer] = {k: v
                           for k, v in p.__dict__.items()
                           if isinstance(v, np.ndarray)}
        self.chains[model] = out

    def _calculate_best_fit(self, model):
        if not model in self.best_fits:
            # create empty entry if it's not there
            self.best_fits[model] = {}
        self._get_chains(model)
        out = copy.deepcopy(self.chains[model])
        for tracer in self.tracers:
            chains = out[tracer]
            argbf = np.argmin(chains["chi2"])
            for par in chains:
                out[tracer][par] = out[tracer][par][argbf]
        self.best_fits[model] = out

    def get_chains(self, model):
        try:
            return self.chains[model]
        except KeyError:
            self._get_chains(model)
            return self.chains[model]

    def get_summary(self, model, par):
        try:
            return self.summary[model][par]
        except (KeyError, TypeError):
            self._calculate_summary(model, par)
            return self.summary[model][par]

    def get_best_fit(self, model):
        try:
            return self.best_fits[model]
        except (KeyError, TypeError):
            self._calculate_best_fit(model)
            return self.best_fits[model]

    def get_sacc_file(self, model, tracer):
        ibin = self.tracers.index(tracer)
        fname = f"{self._chains_dir}/{model}/{model}_{ibin}/cobaya.input.yaml"
        with open(fname, "r") as stream:
            info = yaml.safe_load(stream)
        mod = get_model(info)

        # get best fit
        bf = self.get_best_fit(model)[tracer]
        p_bf = {k: v for k, v in bf.items() if k != "chi2"}

        # get theory sacc object
        loglikes, derived = mod.loglikes(p_bf)
        l = mod.likelihood['ygk_like.ygkLike']
        params = l.current_state['params'].copy()
        params.update(p_bf)
        s_pred = l.get_sacc_file(**params)
        return s_pred
