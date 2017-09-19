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
import re
from .misc import reverse_complement


class Alignment(object):

    def __init__(self, paf_line):
        line_parts = paf_line.strip().split('\t')

        self.read_name = line_parts[0]
        self.read_length = int(line_parts[1])
        self.read_start = int(line_parts[2])
        self.read_end = int(line_parts[3])
        self.read_strand = line_parts[4]
        self.read_seq = ''

        self.ref_name = line_parts[5]
        self.ref_length = int(line_parts[6])
        self.ref_start = int(line_parts[7])
        self.ref_end = int(line_parts[8])

        self.matching_bases = int(line_parts[9])
        self.num_bases = int(line_parts[10])
        self.percent_identity = 100.0 * self.matching_bases / self.num_bases

        self.cigar, self.alignment_score = None, None
        for part in line_parts:
            if part.startswith('cg:Z:'):
                self.cigar = part[5:]
            if part.startswith('AS:i:'):
                self.alignment_score = int(part[5:])
        if self.cigar is None:
            sys.exit('Error: no CIGAR string found')
        if self.alignment_score is None:
            sys.exit('Error: no alignment score')

        self.max_indel = 0
        cigar_parts = re.findall(r'\d+\w', self.cigar)
        for cigar_part in cigar_parts:
            num = int(cigar_part[:-1])
            letter = cigar_part[-1]
            if (letter == 'I' or letter == 'D') and num > self.max_indel:
                self.max_indel = num

        self.read_end_gap = self.read_length - self.read_end
        self.ref_end_gap = self.ref_length - self.ref_end

        # Quality offers a way to compare two alignments against each other.
        # length quality:   https://www.desmos.com/calculator/nman8btt8k
        # identity quality: https://www.desmos.com/calculator/5bcvlafzo7
        # indel quality:    https://www.desmos.com/calculator/2qjk1hinad
        min_identity = 70.0
        half_length_score = 5000.0
        half_indel_score = 25.0
        ref_length = abs(self.ref_end - self.ref_start)
        if self.percent_identity <= min_identity:
            self.quality = 0.0
        else:
            length_quality = 100.0 * (1.0 + (-half_length_score / (ref_length + half_length_score)))
            identity_quality = (self.percent_identity - min_identity) / (100.0 - min_identity)
            indel_quality = half_indel_score / (self.max_indel + half_indel_score)
            self.quality = length_quality * identity_quality * indel_quality

    def __repr__(self):
        return self.read_name + ':' + str(self.read_start) + '-' + str(self.read_end) + \
               '(' + self.read_strand + '),' + \
               self.ref_name + ':' + str(self.ref_start) + '-' + str(self.ref_end) + \
               '(' + ('%.3f' % self.percent_identity) + '%)'

    def add_read_sequence(self, sequence):
        """
        This function takes the whole read sequence and stores the relevant part/strand in the
        object.
        """
        sequence = sequence[self.read_start:self.read_end]
        if self.read_strand == '-':
            sequence = reverse_complement(sequence)
        self.read_seq = sequence

    def get_read_seq_by_ref_coords(self, ref_start, ref_end, ref_seq):
        """
        This function returns the piece of the read sequence which corresponds to the given
        reference coordinates. It steps through the CIGAR string to find the appropriate locations.
        """
        read_seq_start, read_seq_end = None, None
        read_pos, ref_pos = 0, self.ref_start

        # gapped_ref = []  # TEMP
        # gapped_read = []  # TEMP

        match_count = 0  # TEMP
        for i in self.get_expanded_cigar():
            if ref_pos == ref_start:
                read_seq_start = read_pos
            if ref_pos == ref_end:
                read_seq_end = read_pos

            if i == 'M':
                # gapped_ref.append(ref_seq[ref_pos])  # TEMP
                # gapped_read.append(self.read_seq[read_pos])  # TEMP
                if ref_seq[ref_pos] == self.read_seq[read_pos]:  # TEMP
                    match_count += 1  # TEMP
                read_pos += 1
                ref_pos += 1
            elif i == 'I':
                # gapped_ref.append('-')  # TEMP
                # gapped_read.append(self.read_seq[read_pos])  # TEMP
                read_pos += 1
            elif i == 'D':
                # gapped_ref.append(ref_seq[ref_pos])  # TEMP
                # gapped_read.append('-')  # TEMP
                ref_pos += 1

        # print(''.join(gapped_read))  # TEMP
        # print(''.join(gapped_ref))  # TEMP

        if ref_pos == ref_start:
            read_seq_start = read_pos
        if ref_pos == ref_end:
            read_seq_end = read_pos

        assert match_count == self.matching_bases  # TEMP
        assert read_pos + self.read_start == self.read_end
        assert ref_pos == self.ref_end
        assert read_seq_start is not None
        assert read_seq_end is not None
        seq = self.read_seq[read_seq_start:read_seq_end]
        assert len(seq) == read_seq_end - read_seq_start
        return seq

    def get_expanded_cigar(self):
        expanded_cigar = []
        cigar_parts = re.findall(r'\d+\w', self.cigar)
        for cigar_part in cigar_parts:
            num = int(cigar_part[:-1])
            letter = cigar_part[-1]
            if letter == 'M' or letter == 'I' or letter == 'D':
                expanded_cigar.append(letter * num)
        return ''.join(expanded_cigar)
