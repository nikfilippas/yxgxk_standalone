import os
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib import cm

fnames = {"mask_g": 'mask_v3.fits.gz',
          "mask_y": 'mask_planck60.fits.gz',
          "mask_k": 'COM_Lensing_4096_R3.00_mask.fits.gz'}

def plot_map(fname, nside=256, return_data=False):
    fname_use = os.path.join("../data/maps", fname)
    data = hp.ud_grade(hp.read_map(fname_use, dtype=float), nside_out=nside)

    fig = plt.figure(num="mollview", figsize=(12, 7))
    kw = {"map"   : data,
          "fig"   : "mollview",
          "title" : None,
          "xsize" : 1200,
          "cmap"  : cm.gray,
          "cbar"  : False}

    hp.mollview(**kw)
    hp.graticule(c="maroon", alpha=0.4)
    if not return_data:
        return fig
    else:
        return fig, data

nside = 256
data = np.ones(hp.nside2npix(nside))
for mask_name, fname in fnames.items():
    fig, new_data = plot_map(fname, nside=nside, return_data=True)
    fig.savefig(f"../figs/{mask_name}.pdf")
    plt.close()
    data *= new_data
    print(mask_name, np.mean(new_data))

print(np.mean(data))
