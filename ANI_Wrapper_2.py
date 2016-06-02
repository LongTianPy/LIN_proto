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

# def unified_anib(indirname):
#     """
#
#     :param indirnames: The directory where the new genome as well as the top n genomes are in
#     :param org_lengths:
#     :return:
#     """
#     logging = logging.getLogger('ANI_Wrapper_2.py')
#     logging.setLevel(logging.DEBUG)
#     infiles = pyani_files.get_fasta_files(indirname)
#     org_lengths = pyani_files.get_sequence_lengths(infiles)
#     fragsize = pyani_config.FRAGSIZE
#     filestems = pyani_config.ANIB_FILESTEMS
#     logging.info("Running ANIb")
#     # Build BLAST databases and run pairwise BLASTN
#     # Make sequence fragments
#     # Fraglengths does not get reused with BLASTN
#     filenames = os.listdir(indirname)
#     for fname in filenames:
#         if ' ' in  os.path.abspath(fname):
#             logging.error("File or directory '%s' contains whitespace" % fname)
#             logging.error("This will cause issues with MUMmer and BLAST")
#             logging.error("(exiting)")
#             sys.exit(1)
#     fragfiles, fraglengths = anib.fragment_FASTA_files(infiles,
#                                                        indirname,
#                                                        fragsize)
#     format_exe = pyani_config.MAKEBLASTDB_DEFAULT
#     blast_exe = pyani_config.BLASTN_DEFAULT
#
#     # Run BLAST database-building and executables from a jobgraph
#     logging.info("Creating job dependency graph")
#     jobgraph = anib.make_job_graph(infiles, fragfiles, indirname,
#                                    format_exe, blast_exe, "ANIb")
#
#     logging.info("Running jobs with multiprocessing")
#     logging.info("Running job dependency graph")
#     cumval = run_mp.run_dependency_graph(jobgraph, verbose=False,
#                                          logging=logging)
#     if 0 < cumval:
#         logging.warning("At least one BLAST run failed. " +
#                        "ANIb may fail.")
#     else:
#         logging.info("All multiprocessing jobs complete.")
#
#
#     # Process pairwise BLASTN output
#     logging.info("Processing pairwise ANIb BLAST output.")
#     try:
#         data = anib.process_blast(indirname, org_lengths,
#                                   fraglengths=fraglengths, mode="ANIb")
#     except ZeroDivisionError:
#         logging.error("One or more BLAST output files has a problem.")
#         if 0 < cumval:
#             logging.error("This is possibly due to BLASTN run failure, " +
#                          "please investigate")
#         else:
#             logging.error("This is possibly due to ara BLASTN comparison " +
#                          "being too distant for use.")
#         logging.error(last_exception())
#     return data[1]

def unified_anib(indirname,User_ID):
    # Build BLAST databases and run pairwise BLASTN
    # Fraglengths does not get reused with BLASTN
    os.mkdir(indirname+'{0}_out/'.format(User_ID))
    os.system("chmod 777 {0}".format(indirname+'{0}_out'.format(User_ID)))
    logging.basicConfig(level=logging.DEBUG, filename="/home/linproject/Workspace/LIN_log/logfile_{0}".format(User_ID),
                        filemode="a+", format="%(asctime)-15s %(levelname)-8s %(message)s")
    infiles = pyani_files.get_fasta_files(indirname)
    org_lengths = pyani_files.get_sequence_lengths(infiles)
    fragsize = pyani_config.FRAGSIZE
    filestems = pyani_config.ANIB_FILESTEMS
    filenames = os.listdir(indirname)
    for fname in filenames:
        if ' ' in  os.path.abspath(fname):
            logging.error("File or directory '%s' contains whitespace" % fname)
            logging.error("This will cause issues with MUMmer and BLAST")
            logging.error("(exiting)")
            sys.exit(1)
    fragfiles, fraglengths = anib.fragment_FASTA_files(infiles, indirname+'{0}_out/'.format(User_ID), fragsize)
    # Export fragment lengths as JSON, in case we re-run BLASTALL with
    # --skip_blastn
    with open(os.path.join(indirname+'{0}_out/'.format(User_ID), 'fraglengths.json'), 'w') as outfile:
        json.dump(fraglengths, outfile)
    # Which executables are we using?
    format_exe = pyani_config.FORMATDB_DEFAULT
    blast_exe = pyani_config.BLASTALL_DEFAULT
    # Run BLAST database-building and executables from a jobgraph
    logging.info("Creating job dependency graph")
    jobgraph = anib.make_job_graph(infiles, fragfiles, indirname+'{0}_out/'.format(User_ID), format_exe, blast_exe, 'ANIblastall')

    logging.info("Running jobs with multiprocessing")
    logging.info("Running job dependency graph")
    cumval = run_mp.run_dependency_graph(jobgraph, verbose=False,
                                         logging=logging)
    if 0 < cumval:
        logging.warning("At least one BLAST run failed. " +
                       "%s may fail." % 'ANIblastall')
    else:
        logging.info("All multiprocessing jobs complete.")

    # Process pairwise BLASTN output
    logging.info("Processing pairwise %s BLAST output." % 'ANIblastall')
    try:
        data = anib.process_blast(indirname+'{0}_out/'.format(User_ID), org_lengths,
                                  fraglengths=fraglengths, mode='ANIblastall')
    except ZeroDivisionError:
        logging.error("One or more BLAST output files has a problem.")
        if 0 < cumval:
            logging.error("This is possibly due to BLASTN run failure, " +
                         "please investigate")
        else:
            logging.error("This is possibly due to ara BLASTN comparison " +
                         "being too distant for use.")
        logging.error(last_exception())
    return data[1]


if __name__ == "__main__":
    indirname = sys.argv[1]
    print unified_anib(indirname)