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

import itertools
from collections import defaultdict
from .misc import add_line_breaks_to_sequence, load_fasta
from Bio import pairwise2


class UnitigGraph(object):

    def __init__(self, names, seqs, circularity):
        self.segments = {}                      # unsigned seg name -> GraphSegment
        self.forward_links = defaultdict(list)  # signed seg name -> list of signed segment name
        self.reverse_links = defaultdict(list)  # signed seg name <- list of signed segment name
        self.links = {}                         # tuple (start, end) -> StringGraphLink

        for name in names:
            sequence = seqs[name]
            self.segments[name] = GraphSegment(name, sequence)
            if circularity[name]:
                pos_name = name + '+'
                self.add_link(pos_name, pos_name, 0, 0)

    def save_to_gfa(self, filename, include_depth=True):
        """
        Saves whole graph to a GFA file.
        """
        with open(filename, 'w') as gfa:
            for segment in sorted(self.segments.values(), key=lambda x: x.full_name):
                gfa.write(segment.gfa_segment_line(include_depth))
            for link in sorted(self.links.keys()):
                gfa.write(self.links[link].gfa_link_line())

    def save_to_fasta(self, filename, min_length=1):
        with open(filename, 'w') as fasta:
            for segment in sorted(self.segments.values(), reverse=True,
                                  key=lambda x: x.get_length()):
                if segment.get_length() >= min_length:
                    fasta.write(segment.fasta_record())

    def print_fasta_to_stdout(self, names):
        for name in names:
            segment = self.segments[name]
            header = '>' + name + \
                     ' length=' + str(segment.get_length()) + \
                     ' depth=' + '%.2f' % segment.depth + 'x'
            pos_name = name + '+'
            if pos_name in self.links and self.links[pos_name] == [pos_name]:
                header += ' circular=true'
            print(header)
            print(add_line_breaks_to_sequence(segment.forward_sequence, 70), end='')

    def get_preceding_segments(self, seg_name):
        if seg_name not in self.reverse_links:
            return []
        return self.reverse_links[seg_name]

    def get_following_segments(self, seg_name):
        if seg_name not in self.forward_links:
            return []
        return self.forward_links[seg_name]

    def add_link(self, start, end, overlap_1, overlap_2):
        """
        Adds a link to the graph in all necessary ways: forward and reverse, and for reverse
        complements too.
        """
        rev_start = flip_segment_name(start)
        rev_end = flip_segment_name(end)

        if start not in self.forward_links:
            self.forward_links[start] = []
        if end not in self.forward_links[start]:
            self.forward_links[start].append(end)

        if end not in self.reverse_links:
            self.reverse_links[end] = []
        if start not in self.reverse_links[end]:
            self.reverse_links[end].append(start)

        if rev_start not in self.reverse_links:
            self.reverse_links[rev_start] = []
        if rev_end not in self.reverse_links[rev_start]:
            self.reverse_links[rev_start].append(rev_end)

        if rev_end not in self.forward_links:
            self.forward_links[rev_end] = []
        if rev_start not in self.forward_links[rev_end]:
            self.forward_links[rev_end].append(rev_start)

        link_tuple = (start, end)
        self.links[link_tuple] = StringGraphLink(start, end)
        self.links[link_tuple].seg_1_overlap = overlap_1
        self.links[link_tuple].seg_2_overlap = overlap_2

        rev_link_tuple = (rev_end, rev_start)
        self.links[rev_link_tuple] = StringGraphLink(rev_end, rev_start)
        self.links[rev_link_tuple].seg_1_overlap = overlap_2
        self.links[rev_link_tuple].seg_2_overlap = overlap_1

    def segment_is_circular(self, seg_name):
        """
        Returns whether or not the segment has a circularising link.
        """
        pos_seg_name = seg_name + '+'
        preceding_segments = self.get_preceding_segments(pos_seg_name)
        following_segments = self.get_following_segments(pos_seg_name)
        if len(preceding_segments) != 1 or len(following_segments) != 1:
            return False
        preceding_seg_name = preceding_segments[0]
        following_seg_name = following_segments[0]
        return preceding_seg_name == pos_seg_name and following_seg_name == pos_seg_name

    def get_total_segment_length(self):
        return sum(s.get_length() for s in self.segments.values())

    def replace_with_polished_sequences(self, polished_fasta):
        """
        Swaps out the current sequences with polished versions from Racon.
        """
        polished_seqs = load_fasta(polished_fasta)
        for seg_name, segment in self.segments.items():
            try:
                polished_seq = [x[1] for x in polished_seqs if 'Consensus_' + seg_name == x[0]][0]

                # Racon sometimes drops the start or end of sequences, so we do some semi-global
                # alignments to see if bases have been lost. If so, we put them back!
                gap = 500
                unpolished_seq_start = segment.forward_sequence[:gap]
                unpolished_seq_end = segment.forward_sequence[-gap:]
                polished_seq_start = polished_seq[:gap]
                polished_seq_end = polished_seq[-gap:]

                missing_start_seq, missing_end_seq = '', ''

                start_align = pairwise2.align.globalms(unpolished_seq_start, polished_seq_start,
                                                       1, -1, -1, -1, penalize_end_gaps=False,
                                                       one_alignment_only=True)[0][1]
                missing_start_count = sum(1 for _ in itertools.takewhile(lambda c: c == '-',
                                                                         start_align))
                if missing_start_count > 0:
                    missing_start_seq = unpolished_seq_start[:missing_start_count]

                end_align = pairwise2.align.globalms(unpolished_seq_end, polished_seq_end,
                                                     1, -1, -1, -1, penalize_end_gaps=False,
                                                     one_alignment_only=True)[0][1]
                missing_end_count = sum(1 for _ in itertools.takewhile(lambda c: c == '-',
                                                                       reversed(end_align)))
                if missing_end_count > 0:
                    missing_end_seq = unpolished_seq_end[-missing_end_count:]

                if missing_start_seq or missing_end_seq:
                    polished_seq = missing_start_seq + polished_seq + missing_end_seq

                segment.forward_sequence = polished_seq
            except IndexError:
                pass

    def rotate_circular_sequences(self, shift_fraction=0.70710678118655):
        """
        Rotates the sequence to a new starting point. It shifts by a non-rational (well, almost)
        fraction of the sequence length so repeated executions of this function don't result in
        repeated starting positions.
        """
        for seg_name, segment in self.segments.items():
            if self.segment_is_circular(seg_name):
                seq = segment.forward_sequence
                shift = int(len(seq) * shift_fraction)
                seq = seq[shift:] + seq[:shift]
                segment.forward_sequence = seq

    def get_median_read_depth(self):
        """
        Returns the assembly graph's median read depth (by base).
        """
        sorted_segments = sorted(self.segments.values(), key=lambda x: x.depth)
        total_length = 0
        for segment in sorted_segments:
            total_length += segment.get_length()
        halfway_length = total_length // 2
        length_so_far = 0
        for segment in sorted_segments:
            length_so_far += segment.get_length()
            if length_so_far >= halfway_length:
                return segment.depth
        return 0.0

    def normalise_read_depths(self):
        median_depth = self.get_median_read_depth()
        if median_depth == 0.0:
            return
        for segment in self.segments.values():
            segment.depth /= median_depth


class GraphSegment(object):

    def __init__(self, full_name, sequence):
        self.full_name = full_name
        self.forward_sequence = sequence
        self.depth = 1.0

    def __repr__(self):
        if len(self.forward_sequence) > 6:
            seq_string = (self.forward_sequence[:3] + '...' + self.forward_sequence[-3:] + ', ' +
                          str(len(self.forward_sequence)) + ' bp')
        else:
            seq_string = self.forward_sequence
        return self.full_name + ' (' + seq_string + ')'

    def get_length(self):
        return len(self.forward_sequence)

    def gfa_segment_line(self, include_depth=True):
        s_line_parts = ['S', self.full_name, self.forward_sequence,
                        'LN:i:' + str(self.get_length())]
        if include_depth:
            s_line_parts += ['dp:f:' + str(self.depth)]
        return '\t'.join(s_line_parts) + '\n'

    def fasta_record(self):
        return ''.join(['>', self.full_name, '\n',
                        add_line_breaks_to_sequence(self.forward_sequence, 70)])


class StringGraphLink(object):

    def __init__(self, seg_1_signed_name, seg_2_signed_name):
        self.seg_1_signed_name = seg_1_signed_name
        self.seg_2_signed_name = seg_2_signed_name
        self.seg_1_overlap = None
        self.seg_2_overlap = None

    def __repr__(self):
        return (self.seg_1_signed_name + ' -> ' + self.seg_2_signed_name + ' (' +
                str(self.seg_1_overlap) + ', ' + str(self.seg_2_overlap) + ')')

    def gfa_link_line(self):
        seg_1_name = get_unsigned_seg_name(self.seg_1_signed_name)
        seg_1_sign = self.seg_1_signed_name[-1]
        seg_2_name = get_unsigned_seg_name(self.seg_2_signed_name)
        seg_2_sign = self.seg_2_signed_name[-1]
        overlap = str(self.seg_1_overlap) + 'M'
        return '\t'.join(['L', seg_1_name, seg_1_sign, seg_2_name, seg_2_sign, overlap]) + '\n'


def flip_segment_name(seg_name):
    assert(seg_name.endswith('+') or seg_name.endswith('-'))
    if seg_name.endswith('+'):
        return get_unsigned_seg_name(seg_name) + '-'
    else:
        return get_unsigned_seg_name(seg_name) + '+'


def get_unsigned_seg_name(seg_name):
    assert(seg_name.endswith('+') or seg_name.endswith('-'))
    return seg_name[:-1]
