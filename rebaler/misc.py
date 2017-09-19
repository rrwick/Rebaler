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

import sys
import os
import gzip
import shutil
import argparse


class MyHelpFormatter(argparse.HelpFormatter):
    """
    This is a custom formatter class for argparse. It allows for some custom formatting,
    in particular for the help texts with multiple options (like bridging mode and verbosity level).
    http://stackoverflow.com/questions/3853722
    """
    def __init__(self, prog):
        terminal_width = shutil.get_terminal_size().columns
        os.environ['COLUMNS'] = str(terminal_width)
        max_help_position = min(max(24, terminal_width // 3), 40)
        super().__init__(prog, max_help_position=max_help_position)

    def _get_help_string(self, action):
        help_text = action.help
        if action.default != argparse.SUPPRESS and 'default' not in help_text.lower() and \
                action.default is not None:
            help_text += ' (default: ' + str(action.default) + ')'
        return help_text


def float_to_str(num, decimals, max_num=0):
    """
    Converts a number to a string. Will add left padding based on the max value to ensure numbers
    align well.
    """
    if decimals == 0:
        return int_to_str(int(round(num)), max_num=max_num)
    if num is None:
        num_str = 'n/a'
    else:
        num_str = '%.' + str(decimals) + 'f'
        num_str = num_str % num
        parts = num_str.split('.')
        before_decimal = parts[0]
        after_decimal = parts[1]
        num_str = int_to_str(int(before_decimal)) + '.' + after_decimal
    if max_num > 0:
        max_str = float_to_str(max_num, decimals)
        num_str = num_str.rjust(len(max_str))
    return num_str


def int_to_str(num, max_num=0):
    """
    Converts a number to a string. Will add left padding based on the max value to ensure numbers
    align well.
    """
    if num is None:
        num_str = 'n/a'
    else:
        num_str = '{:,}'.format(num)
    max_str = '{:,}'.format(int(max_num))
    return num_str.rjust(len(max_str))


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for filetype, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = filetype
    if compression_type == 'bz2':
        sys.exit('Error: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('Error: cannot use zip format - use gzip instead')
    return compression_type


def get_sequence_file_type(filename):
    """
    Determines whether a file is FASTA or FASTQ.
    """
    if not os.path.isfile(filename):
        sys.exit('Error: could not find ' + filename)
    if get_compression_type(filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open

    with open_func(filename, 'rt') as seq_file:
        try:
            first_char = seq_file.read(1)
        except UnicodeDecodeError:
            first_char = ''

    if first_char == '>':
        return 'FASTA'
    elif first_char == '@':
        return 'FASTQ'
    else:
        raise ValueError('File is neither FASTA or FASTQ')


def load_fasta_or_fastq(filename):
    """
    Returns a list of tuples (header, seq) for each record in the fasta/fastq file.
    """
    try:
        file_type = get_sequence_file_type(filename)
        if file_type == 'FASTA':
            return load_fasta(filename), 'FASTA'
        else:  # FASTQ
            return load_fastq(filename), 'FASTQ'
    except IndexError:
        sys.exit('\nError: ' + filename + ' could not be parsed - is it formatted correctly?')


def load_fasta(fasta_filename):
    """
    Returns a list of tuples (header, seq) for each record in the fasta file.
    """
    if get_compression_type(fasta_filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    fasta_seqs = []
    with open_func(fasta_filename, 'rt') as fasta_file:
        name = ''
        sequence = ''
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    fasta_seqs.append((name.split()[0], sequence, name))
                    sequence = ''
                name = line[1:]
            else:
                sequence += line
        if name:
            fasta_seqs.append((name.split()[0], sequence, name))
    return fasta_seqs


def load_fastq(fastq_filename):
    """
    Returns a list of tuples (header, seq) for each record in the fastq file.
    """
    if get_compression_type(fastq_filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    reads = []
    with open_func(fastq_filename, 'rt') as fastq:
        for line in fastq:
            full_name = line.strip()[1:]
            short_name = full_name.split()[0]
            sequence = next(fastq).strip()
            spacer = next(fastq).strip()
            qualities = next(fastq).strip()
            reads.append((short_name, sequence, spacer, qualities, full_name))
    return reads
