from tqdm import tqdm
import yaml
import numpy as np
import getdist.mcsamples as gmc

from utils import Container


class ChainCalculator(Container):
    def __init__(self, *, base_model="gyksrA_T08", new_interps=False):
        super().__init__(base_model=base_model)
        self.best_fits = {}
        self._zarr = np.asarray(list(self.zmid.values()))

    def _load_chains(self, model, ibin):
        fsetup = f"chains/{model}/{model}_{ibin}/cobaya.input.yaml"
        with open(fsetup, "r") as f:
            config = yaml.safe_load(f)
        fchains = f"chains/{model}/{model}_{ibin}/cobaya"
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
        if not model in self.best_fits:
            # create empty entry if it's not there
            self.best_fits[model] = {}

        bf_arr = np.zeros((6, 3))
        for ibin, tracer in enumerate(tqdm(self.tracers)):
            config, samples, pars = self._load_chains(model, ibin)
            if par in self.interps.parameters:
                der_chain = self._get_derived_chain(config, pars, par, tracer)
                samples.addDerived(
                    der_chain, name=par, label=self._latex_names[par])
            bf_arr[ibin] = self._get_summary(samples, par)

        # normalize by the growth factor for sigma8
        if par == "sigma8":
            a_arr = np.asarray(list(1/(1 + self._zarr)))
            bf_arr *= self.cosmo.growth_factor(a_arr)[:, None]

        self.best_fits[model][par] = bf_arr

    def get_summary(self, model, par):
        try:
            return self.best_fits[model][par]
        except (KeyError, TypeError):
            self._calculate_summary(model, par)
            return self.best_fits[model][par]
