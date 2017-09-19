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
from collections import deque, defaultdict
from .misc import reverse_complement, add_line_breaks_to_sequence, bold, load_fasta
from . import log
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

    def save_to_gfa(self, filename, verbosity=1, newline=False, include_depth=True):
        """
        Saves whole graph to a GFA file.
        """
        log.log(('\n' if newline else '') + 'Saving ' + filename, verbosity)
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

    def remove_segment(self, seg_name_to_remove):
        """
        Removes a segment from the graph and all of its related links.
        """
        def remove_signed_segment(graph, seg_name):
            for preceding_seg_name in graph.get_preceding_segments(seg_name):
                del graph.links[(preceding_seg_name, seg_name)]
                graph.forward_links[preceding_seg_name].remove(seg_name)
            for following_seg_name in graph.get_following_segments(seg_name):
                del graph.links[(seg_name, following_seg_name)]
                graph.reverse_links[following_seg_name].remove(seg_name)
            graph.forward_links.pop(seg_name, None)
            graph.reverse_links.pop(seg_name, None)

        remove_signed_segment(self, seg_name_to_remove + '+')
        remove_signed_segment(self, seg_name_to_remove + '-')
        self.segments.pop(seg_name_to_remove, None)

    def remove_branching_paths(self):
        log.log_explanation('Unicycler removes any links from the string graph which create '
                            'branches. I.e. if any segment has two or more links connected to one '
                            'end, those links are removed. This will result in a graph with only '
                            'simple linear paths that is suitable for creating unambiguous '
                            'bridges.', verbosity=2)
        # Put together a set of all links to be deleted.
        links_to_delete = set()
        for seg_name, segment in self.segments.items():
            pos_seg_name = seg_name + '+'
            neg_seg_name = seg_name + '-'
            following_segments = self.get_following_segments(pos_seg_name)
            preceding_segments = self.get_preceding_segments(pos_seg_name)
            if len(following_segments) > 1:
                for f in following_segments:
                    links_to_delete.add((pos_seg_name, f))
                    links_to_delete.add((flip_segment_name(f), neg_seg_name))
            if len(preceding_segments) > 1:
                for p in preceding_segments:
                    links_to_delete.add((p, pos_seg_name))
                    links_to_delete.add((neg_seg_name, flip_segment_name(p)))

        # Delete all links in the set in each possible way.
        deleted_links = []
        for link in sorted(links_to_delete):
            if link in self.links:
                deleted_links.append(link)
                seg_1, seg_2 = link
                rev_seg_1 = flip_segment_name(seg_1)
                rev_seg_2 = flip_segment_name(seg_2)
                del self.links[(seg_1, seg_2)]
                self.forward_links[seg_1].remove(seg_2)
                self.reverse_links[seg_2].remove(seg_1)
                del self.links[(rev_seg_2, rev_seg_1)]
                self.forward_links[rev_seg_2].remove(rev_seg_1)
                self.reverse_links[rev_seg_1].remove(rev_seg_2)

        if deleted_links:
            log.log('Removed links:', verbosity=2)
            for seg_1, seg_2 in deleted_links:
                log.log('  ' + seg_1 + ' ' + get_right_arrow() + ' ' + seg_2, verbosity=2)
            log.log('', verbosity=2)
        else:
            log.log('No links needed removal', verbosity=2)


    def segment_leads_directly_to_contig_in_both_directions(self, seg_name):
        if self.segments[seg_name].contig:
            return True
        return (self.segment_leads_directly_to_contig(seg_name + '+') and
                self.segment_leads_directly_to_contig(seg_name + '-'))

    def segment_leads_directly_to_contig(self, signed_seg_name):
        """
        Tests whether a given segment leads to a contig via a simple unbranching path. Only tests
        in a single direction.
        """
        starting_seg_name = signed_seg_name
        current_seg_name = signed_seg_name
        while True:
            following_segments = self.get_following_segments(current_seg_name)
            preceding_segments = self.get_preceding_segments(current_seg_name)
            if len(following_segments) != 1 or len(preceding_segments) != 1:
                return False
            if self.segments[get_unsigned_seg_name(current_seg_name)].contig:
                return True
            current_seg_name = following_segments[0]
            if current_seg_name == starting_seg_name:  # Check if we've looped back to the start!
                return False

    def get_bridging_paths(self):
        """
        Returns a list of all bridging paths. The contigs being bridged are included at the start
        and end of each path.
        """
        paths = []
        used_segments = set()
        for seg_name in sorted(self.segments.keys()):
            segment = self.segments[seg_name]
            if not segment.contig and seg_name not in used_segments and \
                    self.segment_leads_directly_to_contig_in_both_directions(seg_name):
                starting_seg = seg_name + '+'
                current_seg = starting_seg
                path = [current_seg]
                while True:
                    current_seg = self.get_following_segments(current_seg)[0]
                    path.append(current_seg)
                    if self.segments[get_unsigned_seg_name(current_seg)].contig:
                        break
                current_seg = starting_seg
                while True:
                    current_seg = self.get_preceding_segments(current_seg)[0]
                    path.insert(0, current_seg)
                    if self.segments[get_unsigned_seg_name(current_seg)].contig:
                        break
                for seg in path:
                    used_segments.add(get_unsigned_seg_name(seg))
                paths.append(path)
        return paths

    def seq_from_signed_seg_name(self, signed_name):
        assert(signed_name.endswith('+') or signed_name.endswith('-'))
        unsigned_seg_name = get_unsigned_seg_name(signed_name)
        if signed_name.endswith('+'):
            return self.segments[unsigned_seg_name].forward_sequence
        else:
            return self.segments[unsigned_seg_name].reverse_sequence

    def save_non_contigs_to_file(self, filename, min_length):
        """
        Saves all graph segments which are not short read contigs to a FASTA file.
        """
        log.log('Saving ' + filename, 1)
        with open(filename, 'w') as fasta:
            for segment in sorted(self.segments.values(), reverse=True,
                                  key=lambda x: x.get_length()):
                if segment.contig or segment.get_length() < min_length:
                    continue
                fasta.write(segment.fasta_record())

    def check_graph_has_no_overlaps(self):
        """
        Asserts that the graph has no branching structures and no overlaps.
        """
        for seg_name in self.segments.keys():
            pos_seg_name = seg_name + '+'
            neg_seg_name = seg_name + '-'
            preceding_segments = self.get_preceding_segments(pos_seg_name)
            following_segments = self.get_following_segments(pos_seg_name)
            assert len(preceding_segments) < 2
            assert len(following_segments) < 2
            if len(preceding_segments) == 1:
                preceding_seg_name = preceding_segments[0]
                start_link = self.links[(preceding_seg_name, pos_seg_name)]
                rev_start_link = self.links[(neg_seg_name, flip_segment_name(preceding_seg_name))]
                assert start_link.seg_1_overlap == 0
                assert start_link.seg_2_overlap == 0
                assert rev_start_link.seg_1_overlap == 0
                assert rev_start_link.seg_2_overlap == 0
            if len(following_segments) == 1:
                following_seg_name = following_segments[0]
                end_link = self.links[(pos_seg_name, following_seg_name)]
                rev_end_link = self.links[(flip_segment_name(following_seg_name), neg_seg_name)]
                assert end_link.seg_1_overlap == 0
                assert end_link.seg_2_overlap == 0
                assert rev_end_link.seg_1_overlap == 0
                assert rev_end_link.seg_2_overlap == 0

    def check_segment_names_and_ranges(self, read_dict, assembly_graph):
        """
        This function looks at the string graph segment names and makes sure that their ranges
        match up with the original sequences. It is to ensure that we haven't screwed anything up
        in the various string graph manipulations we've done.
        """
        for seg_name, seg in self.segments.items():
            assert seg_name == seg.full_name

            # If the string graph segment name contains a range, make sure it matches up with the
            # range in the object.
            try:
                range_in_name = [int(x) for x in seg_name.rsplit(':', 1)[1].split('-')]
            except (IndexError, ValueError):
                range_in_name = [None, None]
            if range_in_name[0] is None:
                assert seg.short_name == seg.full_name
            else:
                assert range_in_name[0] == seg.start_pos
                assert range_in_name[1] == seg.end_pos
                assert seg.short_name == seg_name.rsplit(':', 1)[0]

            if seg_name.startswith('CONTIG_'):
                seg_num = int(seg_name[7:].split(':')[0])
                full_seq = assembly_graph.seq_from_signed_seg_num(seg_num)
            elif seg.short_name in read_dict:  # the segment is a long read
                full_seq = read_dict[seg.short_name].sequence
            else:  # We can't check split or merged reads, so they are skipped.
                full_seq = None

            if full_seq is not None:
                # Miniasm uses 1-based inclusive ranges
                assert seg.forward_sequence == full_seq[seg.start_pos-1:seg.end_pos]

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


    def completed_circular_replicons(self):
        completed_components = []
        single_segment_components = [x for x in self.get_connected_components() if len(x) == 1]
        for component in single_segment_components:
            seg = component[0]
            if self.segment_is_circular(seg):
                completed_components.append(seg)
        return completed_components


    def get_connected_components(self):
        """
        Returns a list of lists, where each inner list is the segment names of one connected
        component of the graph.
        E.g. [[1, 2], [3, 4, 5]] would mean that segments 1 and 2 are in a connected component
        and segments 3, 4 and 5 are in another connected component.
        """
        visited = set()
        components = []
        for v in self.segments:
            if v not in visited:
                component = []
                q = deque()
                q.append(v)
                visited.add(v)
                while q:
                    w = q.popleft()
                    component.append(w)
                    connected_segments = self.get_connected_segments(w)
                    for k in connected_segments:
                        if k not in visited:
                            visited.add(k)
                            q.append(k)
                components.append(sorted(component))

        # Sort (just for consistency from one run to the next)
        return sorted(components)

    def get_connected_segments(self, seg_name):
        """
        Given a segment number, this function returns a list of all other names for segments that
        are directly connected.
        It only returns unsigned segment names (i.e. is not strand-specific).
        """
        connected_segments = set()
        pos_seg_name = seg_name + '+'
        if pos_seg_name in self.forward_links:
            downstream_segments = self.forward_links[pos_seg_name]
            for segment in downstream_segments:
                connected_segments.add(get_unsigned_seg_name(segment))
        if pos_seg_name in self.reverse_links:
            upstream_segments = self.reverse_links[pos_seg_name]
            for segment in upstream_segments:
                connected_segments.add(get_unsigned_seg_name(segment))
        return list(connected_segments)

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
                segment.reverse_sequence = reverse_complement(polished_seq)
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
                segment.reverse_sequence = reverse_complement(seq)

    def get_total_segment_length(self):
        return sum(s.get_length() for s in self.segments.values())

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

    def get_circular_segment_count(self):
        circular_count = 0
        for seg_name in self.segments.keys():
            circular_count += (1 if self.segment_is_circular(seg_name) else 0)
        return circular_count

    def get_linear_segment_count(self):
        linear_count = 0
        for seg_name in self.segments.keys():
            linear_count += (0 if self.segment_is_circular(seg_name) else 1)
        return linear_count


class GraphSegment(object):

    def __init__(self, full_name, sequence, qual=None):
        self.full_name = full_name
        self.forward_sequence = sequence
        self.reverse_sequence = reverse_complement(sequence)
        self.depth = 1.0

        # Miniasm trims reads and puts the start/end positions in the name...
        try:
            name_parts = full_name.rsplit(':', 1)
            self.short_name = name_parts[0]
            self.start_pos, self.end_pos = (int(x) for x in name_parts[1].split('-'))

        # ..but ranges don't apply to some graph segments, like 'merged_reads' segments.
        except (IndexError, ValueError):
            self.short_name = self.full_name
            self.start_pos, self.end_pos = 1, len(self.forward_sequence)

        if self.short_name.startswith('CONTIG_'):
            self.contig = True
            self.qual = settings.CONTIG_READ_QSCORE
        else:
            self.contig = False
            self.qual = None

        # If the constructor gave an explicit quality, we'll use that.
        if qual is not None:
            self.qual = qual

    def __repr__(self):
        if len(self.forward_sequence) > 6:
            seq_string = (self.forward_sequence[:3] + '...' + self.forward_sequence[-3:] + ', ' +
                          str(len(self.forward_sequence)) + ' bp')
        else:
            seq_string = self.forward_sequence
        return self.full_name + ' (' + seq_string + '), mean score = ' + str(self.qual)

    def get_length(self):
        return len(self.forward_sequence)

    def gfa_segment_line(self, include_depth=True):
        s_line_parts = ['S', self.full_name, self.forward_sequence, 'LN:i:' + str(self.get_length())]
        if include_depth:
            s_line_parts += ['dp:f:' + str(self.depth)]
        return '\t'.join(s_line_parts) + '\n'

    def fasta_record(self):
        return ''.join(['>', self.full_name, '\n',
                        add_line_breaks_to_sequence(self.forward_sequence, 70)])

    def rotate_sequence(self, start_pos, flip):
        """
        Rotates the sequence so it begins at start_pos. If flip is True, it also switches the
        forward and reverse strands. This function assumes that the segment is a circular
        completed replicon with no overlap.
        """
        unrotated_seq = self.forward_sequence
        rotated_seq = unrotated_seq[start_pos:] + unrotated_seq[:start_pos]
        rev_comp_rotated_seq = reverse_complement(rotated_seq)

        if flip:
            self.forward_sequence = rev_comp_rotated_seq
            self.reverse_sequence = rotated_seq
        else:
            self.forward_sequence = rotated_seq
            self.reverse_sequence = rev_comp_rotated_seq


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


def get_adjusted_contig_name_and_seq(contig_name, full_seq, start_pos, end_pos):
    """
    This function adjusts the start/end positions in the contig name. It's used when a partial
    contig (from miniasm) has a partial alignment (from minimap).
    """
    sign = contig_name[-1]
    name_parts = contig_name[:-1].rsplit(':', 1)
    base_name = name_parts[0]
    old_start, _ = (int(x) for x in name_parts[1].split('-'))
    new_start = old_start + start_pos
    new_end = (end_pos - start_pos) + new_start - 1

    if sign == '-':
        new_start_pos = len(full_seq) - end_pos
        new_end_pos = len(full_seq) - start_pos
        start_pos, end_pos = new_start_pos, new_end_pos
    adjusted_name = base_name + ':' + str(new_start) + '-' + str(new_end) + sign
    adjusted_seq = full_seq[start_pos:end_pos]

    return adjusted_name, adjusted_seq

def merge_string_graph_segments_into_unitig_graph(string_graph, read_nicknames):
    """
    Creates a unitig graph from a string graph. In essence, reimplements make_unitig_graph function
    in miniasm. Assumes that branching paths have already been removed from the string graph.
    """
    log.log('', verbosity=2)
    unitig_sequences = []
    for component in string_graph.get_connected_components():
        segments_with_dead_ends = []
        for seg_name in component:
            pos_seg_name = seg_name + '+'
            neg_seg_name = seg_name + '-'
            if not string_graph.get_preceding_segments(pos_seg_name):
                segments_with_dead_ends.append(pos_seg_name)
            if not string_graph.get_following_segments(pos_seg_name):
                segments_with_dead_ends.append(neg_seg_name)

        # We should have found either two dead ends (for a linear unitig) or zero dead ends (for a
        # circular unitig).
        assert len(segments_with_dead_ends) == 2 or len(segments_with_dead_ends) == 0
        circular = len(segments_with_dead_ends) == 0

        # If the unitig is circular, then we could start anywhere, so we'll choose the biggest
        # segment (positive strand).
        if circular:
            start_seg_name = sorted(component,
                                    key=lambda x: string_graph.segments[x].get_length())[0] + '+'

        # If the unitig is linear, then we have two possible starting locations. For consistency,
        # we'll take the larger of the two.
        else:
            option_1 = string_graph.segments[get_unsigned_seg_name(segments_with_dead_ends[0])]
            option_2 = string_graph.segments[get_unsigned_seg_name(segments_with_dead_ends[1])]
            if option_1.get_length() >= option_2.get_length():
                start_seg_name = segments_with_dead_ends[0]
            else:
                start_seg_name = segments_with_dead_ends[1]

        # Now we can build the unitig sequence by following the graph outward from the starting
        # segment, always trimming overlaps from the end of segments.
        unitig_seq = ''
        current_seg = start_seg_name
        name_list = []
        while True:
            name_list.append(get_string_graph_segment_nickname(current_seg, read_nicknames))
            current_seq = string_graph.seq_from_signed_seg_name(current_seg)
            next_seg = string_graph.get_following_segments(current_seg)
            if circular:
                assert len(next_seg) == 1  # no dead ends in circular unitigs
            if len(next_seg) == 0:  # no next segment means we've hit the end of a linear unitig
                unitig_seq += current_seq
                break
            else:
                assert len(next_seg) == 1
                overlap = string_graph.links[(current_seg, next_seg[0])].seg_1_overlap
                if overlap == 0:  # I don't think this will happen...
                    unitig_seq += current_seq
                else:
                    unitig_seq += current_seq[:-overlap]
            if circular and next_seg[0] == start_seg_name:  # don't loop endlessly in a circle
                break
            current_seg = next_seg[0]

        arrow = ' ' + get_right_arrow() + ' '
        if circular:
            unitig_sequences.append((unitig_seq, 'circular'))
            log.log(bold('Circular unitig: ') + arrow.join(name_list), verbosity=2)
        else:
            unitig_sequences.append((unitig_seq, 'linear'))
            log.log(bold('Linear unitig: ') + arrow.join(name_list), verbosity=2)
        log.log('', verbosity=2)

    # Build and return the unitig graph using the sequences we just made.
    unitig_sequences = sorted(unitig_sequences, key=lambda x: len(x[0]), reverse=True)
    unitig_graph = StringGraph(None)
    for i, unitig_sequence in enumerate(unitig_sequences):
        unitig_seq, circular_or_linear = unitig_sequence
        unitig_name = str(i+1)
        pos_unitig_name = unitig_name + '+'
        unitig_graph.segments[unitig_name] = GraphSegment(unitig_name, unitig_seq)
        if circular_or_linear == 'circular':
            unitig_graph.add_link(pos_unitig_name, pos_unitig_name, 0, 0)
    return unitig_graph


def get_string_graph_segment_nickname(seg_name, read_nicknames):
    """
    String graph segments are often made from long reads, which may have unwieldy names. This
    function shortens them with the read nicknames, if possible.
    """
    # Can't shorten contig names.
    if seg_name.startswith('CONTIG_'):
        return seg_name
    name_parts = seg_name.rsplit(':', 1)
    if name_parts[0] in read_nicknames:
        return read_nicknames[name_parts[0]] + ':' + name_parts[1]
    else:
        return seg_name
