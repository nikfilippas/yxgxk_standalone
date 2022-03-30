import numpy as np

# Measurements are from:
#   - Pandey19         - Yan20
#   - Chiang20         - Koukoufilippas20
#   - van Waerbeke13   - Vikram17
data = np.genfromtxt("data/bpe_digitized.csv", delimiter=",", skip_header=2)


def auto_extractor(N):
    """Automatically extract data in a useful format
    from the digitized CSV file.

    Arguments:
        N (``int``):
            Increasing number indicating the data set.
    """
    # Columns contain (x, y, y_lo, y_hi).
    useful_cols = [0, 1, 3, 5]
    # There are 6 columns in each model, so pick up those,
    # and from those, only select the useful ones.
    d = data[:, 6*N: 6*(N+1)][:, useful_cols].T
    # If not all rows contain data, remove them.
    d = d[~np.isnan(d)].reshape((len(d), -1))
    return d


bpe_data = [auto_extractor(N) for N in range(6)]
