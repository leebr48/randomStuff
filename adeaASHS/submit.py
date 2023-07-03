#!/bin/python3

# This script was written by Brandon F. Lee (bflee@princeton.edu).

# Import necessary modules
import pandas as pd
import argparse as ap
import os
import subprocess as sp

# Set up command line arguments
desc = 'Automatically read patient IDs from the first column of a CSV file (with no header) and submit ASHS jobs for them to Slurm. ' + \
       'Please note that the "[]" symbols around the defaults listed below are visual only - ' + \
       'they should not be present in your inputs. ' + \
       'ALSO NOTE: to run this script on Oscar, you must call it using python3 after loading the python/3.11.0 and openssl/3.0.0 modules. ' + \
       'You may also need to pip install pandas.'
parser = ap.ArgumentParser(description=desc, formatter_class=ap.ArgumentDefaultsHelpFormatter)
parser.add_argument('csv', type=str, nargs=1, help='CSV file with subject IDs to be run in the first column (with no header). Include path if necessary.')
parser.add_argument('-nProcs', type=int, nargs=1, required=False, default=[1], help='Total number of processors allocated for each ASHS run.')
parser.add_argument('-mem', type=str, nargs=1, required=False, default=['4G'], help='Total amount of memory allocated for each ASHS run. This should be a number follwed by "G" (for gigabytes) with no space.')
parser.add_argument('-time', type=str, nargs=1, required=False, default=['04:00:00'], help='Wall clock time limit each ASHS run. Format is HH:MM:SS.')
parser.add_argument('-email', type=str, nargs=1, required=False, default=[None], help='If specified, this address will receive an email when something goes wrong with an ASHS calculation (that the computer automatically detects).')
parser.add_argument('-noSubmit', action='store_true', default=False, help='If this option is used, sbatch scripts will be written but not submitted to Slurm. This is useful for testing and troubleshooting.')
args = parser.parse_args()

# Validate command line arguments and assign internal variable names
stdComplain = ' Please format it like the default setting.'

patientsToRunFile = os.path.abspath(args.csv[0])

nProcs = args.nProcs[0]

if args.mem[0][-1] != 'G':
    raise IOError('Your memory input is invalid.' + stdComplain)
mem = args.mem[0]

stripTime = args.time[0].strip()
timeHourSplit = stripTime.split(':')
try:
    _ = [int(elem) for elem in timeHourSplit]
except ValueError:
    raise IOError('It appears that at least one of the "hours", "minutes", or "seconds" portions of the time input is incorrect.' + stdComplain)
time = stripTime

if (args.email[0] is not None) and ('@' not in args.email[0]):
    raise IOError('It appears that the input email address is incorrectly formatted.')
email = args.email[0]

noSubmit = args.noSubmit

# Read first column of CSV file and extract patient IDs
df = pd.read_csv(patientsToRunFile, usecols=[0], header=None)
patients = list(df.to_numpy().flatten())

# Write files for Slurm and submit them to the scheduler
s = '#SBATCH '
e = '\n'

for patient in patients:

    # Declare some variables
    filename = 'sbatch_%s.script'%patient
    ashsCmd = 'srun $ASHS_ROOT/bin/ashs_main.sh ' + \
              '-a $ASHS_ROOT/ashs_atlas_upennpmc_20170810 ' + \
              '-g /oscar/data/hoh23/oh-ctnlab/Projects/ADNI_ASHS/data/ADNI2_YEAR2/%s/MPRAGE_[0-9]*.nii.gz '%patient + \
              '-f /oscar/data/hoh23/oh-ctnlab/Projects/ADNI_ASHS/data/ADNI2_YEAR2/%s/HighResHippo_[0-9]*.nii.gz '%patient + \
              '-w /oscar/data/hoh23/oh-ctnlab/Projects/ADNI_ASHS/AR_output/%s '%patient + \
              '-I %s '%patient + \
              '-T'

    with open(filename, 'w') as f:
        
        f.write('#!/bin/bash' + e)
        f.write(e)
        f.write('#Slurm directives:' + e)
        f.write(e)
        f.write(s + '-n ' + str(nProcs) + e)
        f.write(s + '--mem=%s'%mem + e)
        f.write(s + '--time ' + time + e)
        f.write(e)
        f.write(s + '-J ' + 'ASHS_%s'%patient + e)
        f.write(s + '-o ' + 'ASHS_%s-%%j.out'%patient + e)
        f.write(s + '-e ' + 'ASHS_%s-%%j.err'%patient + e)
        f.write(e)

        if email is not None:
            f.write(s + '--mail-type=FAIL,REQUEUE' + e)
            f.write(s + '--mail-user=' + email + e)
            f.write(e)

        f.write('#Commands to run:' + e)
        f.write(e)
        f.write(ashsCmd + e)

    # Submit the job to the scheduler, if desired
    if not noSubmit:
        p = sp.Popen(['sbatch', filename])
