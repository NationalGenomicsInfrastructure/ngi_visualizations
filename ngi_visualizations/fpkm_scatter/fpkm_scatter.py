#!/usr/bin/python
"""
fpkm_scatter.py

Takes FPKM counts from two conditions and makes a scatter plot.
Also calculates a r-squared correlation score.

Either takes two Cufflinks FPKM summary files or a summary file
with multiple samples.
"""

from __future__ import print_function

import argparse
import coloredlogs
import logging
import numpy as np
import os
import re
import sys

# Import matplot lib but avoid default X environment
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def make_fpkm_scatter_plots (input_files, summary=False, output_fn='gene_counts', linear=False, heatmap_fn='heatmap_fn', force_plot=False):
    """
    Main function. Takes input files and makes a plot.
    """

    R_dict = {}
    logging.debug("Starting execution")

    # First off, assume that we have two FPKM
    if summary is False:

        #iterate over all uniqe pairs of input files
        for x in range(len(input_files)):
            for y in range(x+1, len(input_files)):
                # Get proper paths for input files
                input_1 = os.path.realpath(input_files[x])
                input_2 = os.path.realpath(input_files[y])

                # Parse the input files
                sample_1 = load_fpkm_counts(input_1)
                sample_2 = load_fpkm_counts(input_2)

                # What are the sample names?
                sample_1_name = os.path.splitext(os.path.basename(input_1))[0]
                sample_2_name = os.path.splitext(os.path.basename(input_2))[0]

                # Make the plot
                plot_filenames = plot_fpkm_scatter(sample_1, sample_2, sample_1_name, sample_2_name, output_fn=output_fn, linear=linear, force_plot=force_plot)

                if plot_filenames is not None:
                    # Extract the R-value
                    # Get the key from the sample names
                    fn_key = "{}-{}".format(sample_1_name, sample_2_name)
                    R_dict[fn_key] = plot_filenames[fn_key]

    # We have a summary file
    else:
       if len(input_files) !=2:
            logger.critical("'--summary' currently only works for two files")
            sys.exit()

       input_1 = os.path.realpath(input_files[0])
       input_2 = os.path.realpath(input_files[1])

       # Parse the input files
       # returns counts[sample_name][gene_id] = fpkm_count
       condition_1 = load_summary_fpkm_counts(input_1)
       condition_2 = load_summary_fpkm_counts(input_2)

       # Condition basenames
       cond_1_basename = os.path.splitext(os.path.basename(input_1))[0]
       cond_2_basename = os.path.splitext(os.path.basename(input_2))[0]

       # Find project names from sample ID (P1234_*)
       pname_re = re.compile('^P\d+_')
       cond_1_proj = pname_re.search(condition_1.keys()[0]).group(0)
       cond_2_proj = pname_re.search(condition_2.keys()[0]).group(0)

       # Go through each sample pair
       for cond_1_sample in condition_1.keys():
           try:
               cond_2_sample = cond_1_sample.replace(cond_1_proj, cond_2_proj)
               condition_2[cond_2_sample]
               outfile = "{}_{}-{}_{}".format(cond_1_basename, cond_1_sample, cond_2_basename, cond_2_sample)
               plot_filenames = plot_fpkm_scatter(condition_1[cond_1_sample], condition_2[cond_2_sample], cond_1_sample, cond_2_sample, output_fn=outfile, linear=linear)

               # Add to R_dict
               fn_key = "{}-{}".format(cond_1_sample,cond_2_sample )
               R_dict[fn_key]= plot_filenames[outfile]
           except KeyError:
               logger.error("Sample {} not found in {}".format(cond_2_sample, cond_2_basename))
               continue

    # Save R2 values in a matrix to file
    # Extract the sample name, i.e unique keys 'a & b' instead of 'a-b'
    keys = set()
    for i in R_dict.keys():
        keys.update(i.split('-'))
    keys = sorted(list(keys))
    with open('R2_values.tsv', 'w') as f:
        f.write("\t{}\n".format("\t".join(keys)))
        # r for rows and c for columns in the matrix
        for r in keys:
            f.write(r)
            for c in keys:
                # a=a is not a part of the dict (since always R2=1, need to add it to the file to get the rows/columns to match
                if c == r:
                    f.write("\t1")
                    continue
                try:
                    f.write("\t{}".format(R_dict['{}-{}'.format(r,c)]))
                except KeyError:
                    try:
                        f.write("\t{}".format(R_dict['{}-{}'.format(c,r)]))
                    except KeyError:
                        f.write("\t")
                        logger.critical("Couldn't write data for the combination {} and {}".format(c,r))
            f.write("\n")

    # Call the heatmap function
    if len(R_dict.keys()) > 2:
        make_heatmap(R_dict,heatmap_fn)

def load_fpkm_counts (file):
    """
    Reads FPKM file and returns a dict of counts.
    Input: Path to FPKM file
    Output: Dictionary. Keys are transcript names, values are counts.
    """

    counts = {}
    logger.info("Processing {}".format(file))
    try:
        with open(file, 'r') as fh:
            # Read the header
            headers = fh.readline().split("\t")
            try:
                FPKM_idx = headers.index('FPKM')
            except ValueError:
                FPKM_idx = 9
            # Iterate through the file
            try:
                for line in fh:
                    line = line.strip()
                    cols = line.split("\t")
                    gene_id = cols[0]
                    FPKM = cols[FPKM_idx]
                    counts[gene_id] = FPKM
            except IndexError:
                logger.critical ("One of the input files is not in a valid format, are you sure that it's a FPKM file?\nExiting")
                sys.exit()

    except IOError as e:
        logger.error("Error loading cufflinks FPKM file: {}\n{}".format(file, e))
        raise IOError(e)

    return counts


def load_summary_fpkm_counts (file):
    """
    Reads FPKM summary file and returns counts.
    Input: Path to summary FPKM file
    Output: Nested dict. First keys are sample names,
            second keys are gene names, values are counts.
            counts[sample_name][gene_id] = fpkm_count
    """
    # ENSEMBL_ID    Gene_ID    Sample_1    Sample_2    Sample_3    Sample_4
    # ENSG00000000003    TSPAN6    0.1234    0.1234    0.1234    0.1234

    # Read the counts file
    counts = dict()
    try:

        with open(file, 'r') as fh:
            header = fh.readline()
            sample_names = header.split("\t")[2:]
            for name in sample_names:
                counts[name] = {}
            for line in fh:
                sections = line.split("\t")
                gene_id = sections[0]
                for i in range(2, len(sections)-2):
                    counts[sample_names[i]][gene_id] = sections[i]
    except IOError as e:
        logger.error("Error - could not read FPKM summary file: {}\n{}".format(file, e))
        return None

    return counts



def plot_fpkm_scatter (sample_1, sample_2, x_lab, y_lab, output_fn=False, linear=False, force_plot=False):
    """
    Plots scatter plot of FPKM counts. Also calculates r-squared correlation
    and prints this to STDOUT as well as including it in the graph.
    Input: 2 x counts[gene_id] = fpkm_count and 2 x sample names
    Input: title=Plot title
    Returns a dict containing filenames of PNG and PDF graphs as a dict for
    each sample.
    """
    # Output filename
    if output_fn is None:
        output_fn = "{}-{}".format(x_lab, y_lab)

    # SET UP PLOT
    fig = plt.figure()
    axes = fig.add_subplot(111, aspect=1.0)
    x_vals = []
    y_vals = []

    # Collect the paired FPKM counts
    missing_genes = 0
    for gene in sample_1.keys():
        try:
            sample_2[gene]
            x_vals.append(float(sample_1[gene]))
            y_vals.append(float(sample_2[gene]))
        except KeyError:
            missing_genes += 1
            logger.debug("Could not find gene '{}' in sample 2".format(gene))

    missing_genes_pct = float(missing_genes) / float(len(sample_1.keys()))
    if missing_genes > 0:
        logger.error("{: >8} / {: <8} ({:4.1f}%) genes mentioned in '{}' not found in '{}'".format(missing_genes, len(sample_1.keys()), missing_genes_pct*100.0, x_lab, y_lab))

    # Check how many mismatched genes go the other way
    missing_s2_genes = 0
    for gene in sample_2.keys():
        try:
            sample_1[gene]
        except KeyError:
            missing_s2_genes += 1

    missing_s2_genes_pct = float(missing_s2_genes) / float(len(sample_2.keys()))
    if missing_s2_genes > 0:
        logger.error("{: >8} / {: <8} ({:4.1f}%) genes mentioned in '{}' not found in '{}'".format(missing_s2_genes, len(sample_2.keys()), missing_s2_genes_pct*100.0, y_lab, x_lab))


    if max(missing_genes_pct, missing_s2_genes_pct) > 0.3 and not force_plot:
        logger.critical("Percentage of missing genes too high (over 30%)! Aborting '{}' and '{}'.".format(x_lab, y_lab))
        return None

    # Calculate the r squared
    corr = np.corrcoef(x_vals, y_vals)[0,1]
    r_squared = corr ** 2
    logger.warn("R squared for {} = {}".format(output_fn, r_squared))

    # Make the plot
    axes.plot(x_vals, y_vals, 'o', markersize=1)
    plt.xlabel(x_lab)
    plt.ylabel(y_lab)
    plt.title("FPKM Counts for {} and {}".format(x_lab, y_lab))

    # Axes scales
    x1,x2,y1,y2 = axes.axis()
    max_xy = max(x2, y2)
    if linear is True:
        axes.set(xscale='log', xlim=[0,max_xy])
        axes.set(yscale='log', ylim=[0,max_xy])
    else:
        axes.set(xscale='log', xlim=[1,max_xy])
        axes.set(yscale='log', ylim=[1,max_xy])

    # Tidy axes
    axes.set_axisbelow(True)
    axes.tick_params(which='both', labelsize=8, direction='out', top=False, right=False)

    # Add a label about r squared
    plt.subplots_adjust(bottom=0.15)
    plt.text(1, -0.15, r'$r^2$ = {:2f}'.format(r_squared),
                horizontalalignment='right', fontsize=8, transform=axes.transAxes)

    # SAVE OUTPUT
    png_fn = "{}.png".format(output_fn)
    pdf_fn = "{}.pdf".format(output_fn)
    logger.debug("Saving to {} and {}".format(png_fn, pdf_fn))
    plt.savefig(png_fn)
    plt.savefig(pdf_fn)
    plt.close(fig)

    # Return the filenames and R2
    return {'png': png_fn, 'pdf': pdf_fn, output_fn:r_squared }


def make_heatmap(data, heatmap_fn):
    """
    Takes a dict of r2 values and generates a heatmap using matplotlib
    """

    # list the sample names
    names = set()
    for i in data.keys():
        names.update(i.split('-'))
    names = sorted(list(names))

    # If the sample looks like it's form a NGI project, do specific cleaning
    clean_names = []
    for i in names:
        i = i.split(".")[0]   #Remove everything after the first dot
        if i.endswith('Aligned'):
            i = i[:-7]
        clean_names.append(i)

    matrix = []
    # Generate the data matrix (list of lists)
    for idx, x in enumerate(names):
        matrix.append([])
        for idy, y in enumerate(names):
            v = 0
            try:
                v = data["{}-{}".format(x,y)]
            except KeyError:
                try:
                    v = data["{}-{}".format(y,x)]
                except KeyError:
                    if x == y:
                        v = 1
            matrix[idx].append(v)

    ### Turn data into numpy aray:
    data = np.array(matrix)
    fig, ax = plt.subplots()


    ## Dynamically change font size
    if len(clean_names) >= 20 :
        ax.xaxis.label.set_fontsize(4)
        ax.yaxis.label .set_fontsize(4)
    elif len(clean_names) >= 10:
        ax.xaxis.label.set_fontsize(8)
        ax.yaxis.label .set_fontsize(8)


    #set size
    heatmap = ax.pcolor(data, cmap='YlOrRd', vmin=0, vmax=1)
    fig = plt.gcf()
    ax.set_frame_on(False)

    #make the plot square
    plt.gca().set_aspect('equal', adjustable='box')

    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(data.shape[0]) +0.5 , minor=False)
    ax.set_xticks(np.arange(data.shape[1]) +0.5 , minor=False)

    # want a more natural, table-like display
    ax.invert_yaxis()

    # adjust plot position to make room for longer sample names.
    # This does creat a bit of whitespace
    plt.subplots_adjust(bottom=0.4, right=None, left=0.3, top = None)

    plt.colorbar(heatmap, ax=ax)

    # Set the labels
    labels = clean_names
    ax.set_xticklabels(labels, minor=False)
    ax.set_yticklabels(labels, minor=False)

    # rotate the plot
    plt.xticks(rotation=90)
    ax.grid(False)
    # Turn off all the ticks
    ax = plt.gca()

    plt.title("Heatmap of R2 values")
    #remove small ticks
    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    logger.info("Saving heatmap to {}.png and {}.pdf".format(heatmap_fn,heatmap_fn))
    pdf_fn = "{}.pdf".format(heatmap_fn)
    png_fn = "{}.png".format(heatmap_fn)
    plt.savefig(pdf_fn)
    plt.savefig(png_fn)
    plt.close(fig)


if __name__ == "__main__":
    # Command line arguments
    parser = argparse.ArgumentParser("Make a scatter plot of FPKM counts between conditions")
    parser.add_argument("-s", "--summary", dest="summary", action='store_true',
                        help="Input files are summary FPKM files.")
    parser.add_argument("-o", "--output", dest="output_fn", default=None,
                        help="Plot output filename base. Default: sample1_sample2.png / .pdf")
    parser.add_argument("-n", "--linear", dest="linear", action='store_true',
                        help="Plot using linear axes instead of log10")
    parser.add_argument("-u", "--log-output", dest="log_output", default='stdout',
                        help="Log output filename. Default: stdout")
    parser.add_argument("input_files", metavar='<input files>', nargs='+',
                        help="List of summary FPKM / cufflinks FPKM results files. See README.MD for more information.")
    parser.add_argument("-f", "--heatmap", dest="heatmap_fn", default='heatmap',
                        help="Name of heatmap output")
    parser.add_argument("-v", "--verbose", dest="verbose", action='store_true',
                        help="Use verbose logging")
    parser.add_argument("-a", "--force", dest="force_plot", action='store_true',
                        help="Force generation of plot, even if mismatched gene names")
    kwargs = vars(parser.parse_args())

    # Initialise logger
    loglevel = 'DEBUG' if kwargs['verbose'] else 'INFO'
    if kwargs['log_output'] != 'stdout':
        logging.basicConfig(filename=kwargs['log_output'])
        logger = logging.getLogger('fpkm_scatter')
        logger.setLevel(loglevel)
    else:
        logger = logging.getLogger('fpkm_scatter')
        coloredlogs.DEFAULT_FIELD_STYLES['levelname'] = {'color': 'blue' }
        coloredlogs.install(level=loglevel, fmt='[ %(levelname)8s ] %(message)s')

    # Remove logging parameters
    kwargs.pop('verbose', None)
    kwargs.pop('log_output', None)

    # Call plot_observed_genes()
    make_fpkm_scatter_plots(**kwargs)
