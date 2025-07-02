import numpy as np

# Divide x into chunks with gaps longer than x_gap
# Returns a list of index pairs (i, j) such that x[i:j] are the chunks


def divide_into_chunks(x, x_gap):

    i_breaks = np.where(np.diff(x) > x_gap)[0]
    i_lower = [0, *(i_breaks+1)]
    i_upper = [*(i_breaks+1), len(x)]

    return list(zip(i_lower, i_upper))
