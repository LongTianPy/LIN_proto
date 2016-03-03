"""
This script takes the advantage of PyANI (https://github.com/widdowquinn/pyani) to calculate Average Nucleotide
Identity (ANI). And forward the result to AssignID.
Of the methods provided by PyANI, we are only using ANIb.
"""
# ANIb: FASTA sequences describing 1000nt fragments of each input sequence;
#       BLAST nucleotide databases - one for each set of fragments; and BLASTN
#       output files (tab-separated tabular format plain text) - one for each
#       pairwise comparison of input sequences. There are potentially a lot of
#       intermediate files.
#
# DEPENDENCIES
# ============
#
# o Biopython (http://www.biopython.org)
#
# o BLAST+ executable in the $PATH, or available on the command line (ANIb)
#       (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
#
# USAGE
# =====
# Options:
#   -h, --help            show this help message and exit
#   -o OUTDIRNAME, --outdir=OUTDIRNAME
#                         Output directory
#   -i INDIRNAME, --indir=INDIRNAME
#                         Input directory name
#   -f, --force           Force file overwriting
#   -s, --fragsize        Sequence fragment size for ANIb
#   --blast_exe=BLAST_EXE
#                         Path to BLASTN+ executable
#   --makeblastdb_exe=MAKEBLASTDB_EXE
#                         Path to BLAST+ makeblastdb executable





# IMPORT
import json
import logging
import logging.handlers
import os
import shutil
import sys
import traceback

from argparse import ArgumentParser

from pyani import average_nucleotide_identity
from pyani import anib
from pyani import run_multiprocessing as run_mp
from pyani.pyani_config import params_mpl, params_r

# FUNCTIONS

# Process command-line arguments
def parse_cmdline(args):
    parser = ArgumentParser(prog="average_nucleotide_identity.py")
    parser.add_argument("-o", "--outdir", dest="outdirname", action="store",
                        default="None", help="Output directory")
    parser.add_argument("-i", "--indir", dest="indirname", action="store",
                        default="None", help="Input directory")
    parser.add_argument("-f", "--force", dest="force", action="store_true",
                        default=False, help="Force file overwriting")
    parser.add_argument("-s", "--fragsize", const="fragsize",
                        action="store_const", default=pyani_config.FRAGSIZE,
                        help="Sequence fragment size for ANIb")
    parser.add_argument("--blastn_exe", dest="blastn_exe",
                        action="store", default=pyani_config.BLASTN_DEFAULT,
                        help="Path to BLASTN+ executable")
    parser.add_argument("--makeblastdb_exe", dest="makeblastdb_exe",
                        action="store",
                        default=pyani_config.MAKEBLASTDB_DEFAULT,
                        help="Path to BLAST+ makeblastdb executable")
    parser.add_argument("-v", "--verbose", dest="verbose",
                        action="store_true", default=False,
                        help="Give verbose output")
    parser.add_argument("-l", "--logfile", dest="logfile",
                        action="store", default=None,
                        help="Logfile location")
    parser.add_argument("-m", "--method", dest="method",
                        action="store", default="ANIb",
                        help="Do not change this option")

# Report the last exception as string
def last_exception():
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value, exc_traceback))

# Create output directory if it doesn't exist
def make_dir():
    """
    By default, this function abort the script if the output directory already exists,
    but we can force proceed by overwriting if the --force option is switched on.
    """
    if os.path.exists(args.outdirname):
        if not args.force:
            logger.error("Output directory %s would overwrite existing files (exiting)"%args.outdirname)
            sys.exit(1)
        else:
            logger.info("Removing directory %s and everything under it."%args.outdirname)
            shutil.rmtree(args.outdirname)
    logger.info("Creating directory %s"%args.outdirname)
    try:
        os.makedirs(args.outdirname)
    except:
        if args.force:
            logger.info("FORCE: not creating directory")
        else:
            logger.error(last_exception())
            sys.exit(1)

# Calculate ANIb for input
# This is originally from PyANI, with small modifications according to our case
def unified_anib(infiles, org_lengths):
    """Calculate ANIb for files in input directory.
    - infiles - paths to each input file
    - org_lengths - dictionary of input sequence lengths, keyed by sequence
    Calculates ANI by the ANIb method, as described in Goris et al. (2007)
    Int J Syst Evol Micr 57: 81-91. doi:10.1099/ijs.0.64483-0. There are
    some minor differences depending on whether BLAST+ or legacy BLAST
    (BLASTALL) methods are used.
    All FASTA format files (selected by suffix) in the input directory are
    used to construct BLAST databases, placed in the output directory.
    Each file's contents are also split into sequence fragments of length
    options.fragsize, and the multiple FASTA file that results written to
    the output directory. These are BLASTNed, pairwise, against the
    databases.
    The BLAST output is interrogated for all fragment matches that cover
    at least 70% of the query sequence, with at least 30% nucleotide
    identity over the full length of the query sequence. This is an odd
    choice and doesn't correspond to the twilight zone limit as implied by
    Goris et al. We persist with their definition, however.  Only these
    qualifying matches contribute to the total aligned length, and total
    aligned sequence identity used to calculate ANI.
    The results are processed to give matrices of aligned sequence length
    (aln_lengths.tab), similarity error counts (sim_errors.tab), ANIs
    (perc_ids.tab), and minimum aligned percentage (perc_aln.tab) of
    each genome, for each pairwise comparison. These are written to the
    output directory in plain text tab-separated format.
    """
    logger.info("Running ANIb")
    # Build BLAST databases and run pairwise BLASTN
    # Make sequence fragments
    logger.info("Fragmenting input files, and writing to %s" %
                args.outdirname)
    # Fraglengths does not get reused with BLASTN
    fragfiles, fraglengths = anib.fragment_FASTA_files(infiles,
                                                       args.outdirname,
                                                       args.fragsize)
    format_exe = args.makeblastdb_exe
    blast_exe = args.blastn_exe

    # Run BLAST database-building and executables from a jobgraph
    logger.info("Creating job dependency graph")
    jobgraph = anib.make_job_graph(infiles, fragfiles, args.outdirname,
                                   format_exe, blast_exe, "ANIb")

    logger.info("Running jobs with multiprocessing")
    logger.info("Running job dependency graph")
    cumval = run_mp.run_dependency_graph(jobgraph, verbose=args.verbose,
                                         logger=logger)
    if 0 < cumval:
        logger.warning("At least one BLAST run failed. " +
                       "ANIb may fail.")
    else:
        logger.info("All multiprocessing jobs complete.")


    # Process pairwise BLASTN output
    logger.info("Processing pairwise ANIb BLAST output.")
    try:
        data = anib.process_blast(args.outdirname, org_lengths,
                                  fraglengths=fraglengths, mode="ANIb")
    except ZeroDivisionError:
        logger.error("One or more BLAST output files has a problem.")
        if 0 < cumval:
            logger.error("This is possibly due to BLASTN run failure, " +
                         "please investigate")
        else:
            logger.error("This is possibly due to a BLASTN comparison " +
                         "being too distant for use.")
        logger.error(last_exception())
    return data

# Write ANIb output
def write(results, filestems):
    """Write ANIb/ANIm/TETRA results to output directory.
    - results - tuple of dataframes from analysis
    Each dataframe is written to an Excel-format file (if args.write_excel is
    True), and plain text tab-separated file in the output directory. The
    order of result output must be reflected in the order of filestems.
    """
    logger.info("Writing %s results to %s" % ("ANIb", args.outdirname))
    for df, filestem in zip(results, filestems):
        logger.info("\t%s" % filestem)
        df.to_csv(os.path.join(args.outdirname, filestem) + '.tab',
                  index=True, sep="\t")

# Run as script
if __name__ == '__main__':

    # Parse command-line
    args = parse_cmdline(sys.argv)

    # Set up logging
    logger = logging.getLogger('average_nucleotide_identity.py')
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)

    # Was a logfile specified? If so, use it
    if args.logfile is not None:
        try:
            logstream = open(args.logfile, 'w')
            err_handler_file = logging.StreamHandler(logstream)
            err_handler_file.setFormatter(err_formatter)
            err_handler_file.setLevel(logging.INFO)
            logger.addHandler(err_handler_file)
        except:
            logger.error("Could not open %s for logging" %
                         args.logfile)
            sys.exit(1)

    # Do we need verbosity?
    if args.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)

    # Report arguments, if verbose
    logger.info(args)

    # Have we got an input and output directory? If not, exit.
    if args.indirname is None:
        logger.error("No input directory name (exiting)")
        sys.exit(1)
    logger.info("Input directory: %s" % args.indirname)
    if args.outdirname is None:
        logger.error("No output directory name (exiting)")
        sys.exit(1)
    make_outdir()
    logger.info("Output directory: %s" % args.outdirname)

    # Check for the presence of space characters in any of the input filenames
    # or output directory. If we have any, abort here and now.
    filenames = [args.outdirname] + os.listdir(args.indirname)
    for fname in filenames:
        if ' ' in  os.path.abspath(fname):
            logger.error("File or directory '%s' contains whitespace" % fname)
            logger.error("This will cause issues with MUMmer and BLAST")
            logger.error("(exiting)")
            sys.exit(1)

    # Have we got a valid method choice?
    # Dictionary below defines analysis function, and output presentation
    # functions/settings, dependent on selected method.
    methods = {"ANIm": (calculate_anim, pyani_config.ANIM_FILESTEMS),
               "ANIb": (unified_anib, pyani_config.ANIB_FILESTEMS),
               "TETRA": (calculate_tetra, pyani_config.TETRA_FILESTEMS),
               "ANIblastall": (unified_anib,
                               pyani_config.ANIBLASTALL_FILESTEMS)}
    if args.method not in methods:
        logger.error("ANI method %s not recognised (exiting)" % args.method)
        logger.error("Valid methods are: %s" % methods.keys())
        sys.exit(1)
    logger.info("Using ANI method: %s" % args.method)

    # Get input files
    logger.info("Identifying FASTA files in %s" % args.indirname)
    infiles = pyani_files.get_fasta_files(args.indirname)
    logger.info("Input files:\n\t%s" % '\n\t'.join(infiles))

    # Get lengths of input sequences
    logger.info("Processing input sequence lengths")
    org_lengths = pyani_files.get_sequence_lengths(infiles)
    logger.info("Sequence lengths:\n" +
                os.linesep.join(["\t%s: %d" % (k, v) for
                                 k, v in org_lengths.items()]))

    # Run appropriate method on the contents of the input directory,
    # and write out corresponding results.
    logger.info("Carrying out %s analysis" % args.method)
    results = methods[args.method][0](infiles, org_lengths)
    write(results, methods[args.method][1])

    # Report that we've finished
    logger.info("Done.")