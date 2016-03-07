#!/usr/bin/python
"""
This script wraps the ANIb calculation part of PyANI into functions. (The original edition uses a lot of lines outside
of functions, which is fine for commandline usage.) I want this script be able to take input from the forms of front
end, as well as the k-mer result. So we are not using the module argparse to load options from input, and probably no
output directory needed.
"""

# IMPORT
from pyani import anib, pyani_config, pyani_files
from pyani import run_multiprocessing as run_mp
from pyani.pyani_config import params_mpl, params_r
import json
import logging
import logging.handlers
import os
import shutil
import sys
import traceback

def unified_anib(indirname):
    """

    :param indirnames: The directory where the new genome as well as the top n genomes are in
    :param org_lengths:
    :return:
    """
    logger = logging.getLogger('ANI_Wrapper_2.py')
    logger.setLevel(logging.DEBUG)
    infiles = pyani_files.get_fasta_files(indirname)
    org_lengths = pyani_files.get_sequence_lengths(infiles)
    fragsize = pyani_config.FRAGSIZE
    filestems = pyani_config.ANIB_FILESTEMS
    logger.info("Running ANIb")
    # Build BLAST databases and run pairwise BLASTN
    # Make sequence fragments
    # Fraglengths does not get reused with BLASTN
    filenames = os.listdir(indirname)
    for fname in filenames:
        if ' ' in  os.path.abspath(fname):
            logger.error("File or directory '%s' contains whitespace" % fname)
            logger.error("This will cause issues with MUMmer and BLAST")
            logger.error("(exiting)")
            sys.exit(1)
    fragfiles, fraglengths = anib.fragment_FASTA_files(infiles,
                                                       indirname,
                                                       fragsize)
    format_exe = pyani_config.MAKEBLASTDB_DEFAULT
    blast_exe = pyani_config.BLASTN_DEFAULT

    # Run BLAST database-building and executables from a jobgraph
    logger.info("Creating job dependency graph")
    jobgraph = anib.make_job_graph(infiles, fragfiles, indirname,
                                   format_exe, blast_exe, "ANIb")

    logger.info("Running jobs with multiprocessing")
    logger.info("Running job dependency graph")
    cumval = run_mp.run_dependency_graph(jobgraph, verbose=False,
                                         logger=logger)
    if 0 < cumval:
        logger.warning("At least one BLAST run failed. " +
                       "ANIb may fail.")
    else:
        logger.info("All multiprocessing jobs complete.")


    # Process pairwise BLASTN output
    logger.info("Processing pairwise ANIb BLAST output.")
    try:
        data = anib.process_blast(indirname, org_lengths,
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
    return data[1]

if __name__ == "__main__":
    indirname = sys.argv[1]
    print unified_anib(indirname)