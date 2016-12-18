#! /usr/local/bin/python
"""
mean_size.py
Created by Tim Stuart
"""
 
import numpy as np
 
 
def get_data(inp):
    lengths = []
    for line in inp:
        if line.startswith('@'):
            pass
        else:
            line = line.rsplit()
            length = int(line[8])
            if length > 0:
                lengths.append(length)
            else:
                pass
    return lengths
 
 
def reject_outliers(data, m=2.):
    """
    rejects outliers more than 2
    standard deviations from the median
    """
    median = np.median(data)
    std = np.std(data)
    for item in data:
        if abs(item - median) > m * std:
            data.remove(item)
        else:
            pass
 
 
def calc_size(data):
    mn = int(np.mean(data))
    std = int(np.std(data))
    minval = mn - std
    maxval = mn + std
    return minval, maxval, mn, std
 
 
if __name__ == "__main__":
    import sys
    lengths = get_data(sys.stdin)
    reject_outliers(lengths)
    #print calc_size(lengths) 
    minval, maxval, mn, std = calc_size(lengths)
    print minval, maxval, mn, std
