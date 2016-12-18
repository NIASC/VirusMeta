# NOTE: This script is not part of the standard 'khmer' distribution.
"""
    Plot one or more k-mer abundance distributions.
"""


import sys
import argparse

import operator

import numpy                    as np
import matplotlib		as mpl
mpl.use( "svg" )  # TODO? Allow user to choose output format.
import matplotlib.pyplot        as plt


count_and_num = operator.itemgetter( 0, 1 )


def main( ):
    
    parser = \
    argparse.ArgumentParser(
        description = "Plot k-mer abundance distribution(s) to a graphics file."
    )
    parser.add_argument(
        "-o", "--output-file", metavar = "OUTPUT-FILE",
        type = str, default = "kmer-abund-plot.svg",
        help = "Name of output file containing figure."
    )
    parser.add_argument(
        "-x", "--x-axis-max", metavar = "XMAX", type = int, default = -1,
	help = "Largest value on x-axis."
    )
    parser.add_argument(
        "-y", "--y-axis-max", metavar = "YMAX", type = int, default = -1,
	help = "Largest value on y-axis."
    )
    parser.add_argument(
        "abundance_dists", metavar = "INPUT-FILE", nargs = "+",
        type = argparse.FileType( "r" ),
        help = "Name of file containing a 'khmer' abundance distribution."
    )
    # TODO? Add option for plotting histograms under the curves.

    args = parser.parse_args( )
    ifile_names = args.abundance_dists
    
    for ifile_name in ifile_names:
        points = np.loadtxt( ifile_name, dtype = 'uint32', usecols = [ 0, 1 ] )
        plt.plot( points )
	if -1 < args.x_axis_max: plt.axis( xmax = args.x_axis_max )
	if -1 < args.y_axis_max: plt.axis( ymax = args.y_axis_max )

    plt.savefig( args.output_file )
    plt.close( )


if "__main__" == __name__:

    main( )

# vim: set ft=python ts=4 sts=4 sw=4 et tw=79:
