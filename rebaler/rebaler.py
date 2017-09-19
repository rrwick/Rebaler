# -*- coding: utf-8 -*-
"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Rebaler

This file is part of Rebaler. Rebaler is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Rebaler is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Rebaler. If
not, see <http://www.gnu.org/licenses/>.
"""

    
import argparse
import multiprocessing
import os
import sys
import subprocess
from .misc import MyHelpFormatter
from .alignment import Alignment


__version__ = "0.1.0"


def main():
    args = get_arguments()

    # Align reads to the reference with minimap2.
    alignments = []
    command = ['minimap2', '-c', '-x', 'map-ont', '-t', str(args.threads),
               args.reference, args.reads]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    while process.poll() is None:
        paf_line = process.stdout.readline().rstrip().decode()
        if paf_line:
            alignments.append(Alignment(paf_line))


def get_arguments():
    """
    Parse the command line arguments.
    """
    default_threads = min(multiprocessing.cpu_count(), 16)

    parser = argparse.ArgumentParser(description='Rebaler: reference-based long read assemblies '
                                                 'of bacterial genomes',
                                     formatter_class=MyHelpFormatter)
    parser.add_argument('-t', '--threads', type=int, default=default_threads,
                        help='Number of threads to use for alignment')
    parser.add_argument('reference', type=str,
                        help='FASTA file of reference assembly')
    parser.add_argument('reads', type=str,
                        help='FASTA/FASTQ file of long reads')

    if len(sys.argv) == 1:
        parser.print_help(file=sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    if not os.path.isfile(args.reference):
        sys.exit('Error: could not find ' + args.reference)
    if not os.path.isfile(args.reads):
        sys.exit('Error: could not find ' + args.reference)

    return args

