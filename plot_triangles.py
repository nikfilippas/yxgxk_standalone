import matplotlib.pyplot as plt
import getdist.plots as gplot
import getdist.mcsamples as gmc

burn = 0.3
for ibin in range(6):
    s = gmc.loadMCSamples(f'yxgxk_{ibin}/cobaya', settings={'ignore_rows': burn})

    gdplot = gplot.get_subplot_plotter()
    gdplot.triangle_plot([s], ["sigma8", "ygk_mass_bias", "ygk_rhogy", f"ygk_LOWZ__{ibin}_lMmin_0", f"ygk_LOWZ__{ibin}_lM1_0"], filled=True)
    plt.savefig("triang.pdf")
    plt.show()
