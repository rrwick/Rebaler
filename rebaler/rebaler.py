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
import itertools
import collections

from .misc import MyHelpFormatter, load_fasta, load_fasta_or_fastq, int_to_str, float_to_str,\
    print_table, red, reverse_complement
from .alignment import Alignment
from . import log
from .unitig_graph import UnitigGraph


__version__ = "0.1.0"


def main():
    args = get_arguments()
    log.logger = log.Log(None, 1)

    log.log_section_header('Loading reference')
    log.log_explanation('This reference sequence will be used as a template for the Rebaler '
                        'assembly.')
    reference = load_fasta(args.reference)
    ref_names = [x[0] for x in reference]
    circularity = {x[0]: 'circular=true' in x[2].lower() for x in reference}
    ref_seqs = {x[0]: x[1] for x in reference}
    print_ref_info(ref_names, ref_seqs, circularity)

    log.log_section_header('Building unpolished assembly')
    log.log_explanation('Rebaler first aligns long reads to the reference using minimap2. It then '
                        'selects high quality alignments and replaces the reference sequence with '
                        'the corresponding read sequence. This creates an unpolished assembly '
                        'made directly from read fragments, similar to what would be produced by '
                        'miniasm.')
    print('Loading reads...                             ', flush=True, end='')
    reads, _ = load_fasta_or_fastq(args.reads)
    print(int_to_str(len(reads)) + ' reads', flush=True)

    print('Aligning reads to reference with minimap2... ', flush=True, end='')
    alignments = get_initial_alignments(args)
    print(int_to_str(len(alignments)) + ' initial alignments', flush=True)

    print('Culling alignments to a non-redundant set... ', flush=True, end='')
    alignments, depths = cull_alignments(alignments, reference)
    print(int_to_str(len(alignments)) + ' alignments remain', flush=True)

    print('Constructing unpolished sequences...         ', flush=True, end='')
    store_read_seqs_in_alignments(alignments, reads)
    partitions = partition_reference(reference, alignments)
    unpolished_sequences = get_unpolished_sequences(partitions, ref_seqs)
    print('done', flush=True)

    log.log_section_header('Polishing assembly')
    log.log_explanation('Rebaler now runs Racon to polish the miniasm assembly. It does '
                        'multiple rounds of polishing to get the best consensus. Circular unitigs '
                        'are rotated between rounds such that all parts (including the ends) are '
                        'polished well. Assembly quality is measured by the sum of all read '
                        'alignment scores.')
    polish_assembly_with_racon(ref_names, unpolished_sequences, circularity, args.reads,
                               args.threads)


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


def print_ref_info(names, seqs, circularity):
    table = [['Reference contig', 'Circular', 'Length']]
    for name in names:
        table.append([name, 'yes' if circularity[name] else 'no', int_to_str(len(seqs[name]))])
    print_table(table, left_align_header=False, alignments='LLR', indent=0)


def get_initial_alignments(args):
    alignments = []
    command = ['minimap2', '-c', '-x', 'map-ont', '-t', str(args.threads),
               args.reference, args.reads]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    while process.poll() is None:
        paf_line = process.stdout.readline().rstrip().decode()
        if paf_line:
            alignments.append(Alignment(paf_line))
    return alignments


def cull_alignments(alignments, reference):
    """
    Takes the alignments and removes redundant ones, preferentially keeping higher quality reads.
    """
    depths = get_reference_depths(alignments, reference)

    alignments = sorted(alignments, reverse=True, key=lambda x: x.quality)
    alignment_count = len(alignments)

    for i in range(alignment_count-1, 0, -1):
        a = alignments[i]
        name = a.ref_name
        can_delete = all(x > 1 for x in depths[name][a.ref_start:a.ref_end])
        if can_delete:
            alignments.pop(i)
            for j in range(a.ref_start, a.ref_end):
                depths[name][j] -= 1

    return alignments, depths


def get_reference_depths(alignments, reference):
    """
    Returns per-base depth of coverage for the reference contigs.
    """
    depths = {}
    for name, seq, _ in reference:
        depths[name] = [0] * len(seq)
    for a in alignments:
        name = a.ref_name
        for i in range(a.ref_start, a.ref_end):
            depths[name][i] += 1
    return depths


def store_read_seqs_in_alignments(alignments, reads):
    read_seqs = {}
    for r in reads:
        read_seqs[r[0]] = r[1]
    for a in alignments:
        a.add_read_sequence(read_seqs[a.read_name])


def partition_reference(reference, alignments):
    """
    This function splits the reference up into chunks based on the alignments. Each chunk has the
    start and end positions (in reference coordinates) and the associated alignment.
    """
    partitions = {}
    for name, seq, _ in reference:
        partitions[name] = []
        ref_len = len(seq)
        ref_alignments = [x for x in alignments if x.ref_name == name]
        ref_alignments = sorted(ref_alignments, key=lambda x: x.ref_start)

        for i in range(len(ref_alignments)):
            a = ref_alignments[i]
            if i > 0:
                a_prev = ref_alignments[i-1]
                prev_overlap = a_prev.ref_end > a.ref_start
            else:
                a_prev = None
                prev_overlap = False
            try:
                a_next = ref_alignments[i+1]
                next_overlap = a.ref_end > a_next.ref_start
            except IndexError:
                a_next = None
                next_overlap = False

            # Make the starting gap, if appropriate.
            if a_prev is None and a.ref_start > 0:
                partitions[name].append((0, a.ref_start, None))
            if a_prev is not None and a_prev.ref_end < a.ref_start:
                partitions[name].append((a_prev.ref_end, a.ref_start, None))

            # Make the depth=1 partition.
            if prev_overlap:
                start = a_prev.ref_end
            else:
                start = a.ref_start
            if next_overlap:
                end = a_next.ref_start
            else:
                end = a.ref_end
            partitions[name].append((start, end, a))

            # Make the depth=2 partition.
            if next_overlap:
                best_a = a if a.percent_identity > a_next.percent_identity else a_next
                partitions[name].append((a_next.ref_start, a.ref_end, best_a))

            # Make the ending gap, if appropriate.
            if a_next is None and a.ref_end < ref_len:
                partitions[name].append((a.ref_end, ref_len, None))

    return partitions


def get_unpolished_sequences(partitions, ref_seqs):
    """
    This function goes through the partitions and returns
    """
    unpolished_sequences = {}
    for name, ref_partitions in partitions.items():
        seq_parts = []
        # print(name)  # TEMP
        ref_seq = ref_seqs[name]
        for start, end, alignment in ref_partitions:
            # print('  ', str(start) + '-' + str(end))  # TEMP
            # print(alignment)  # TEMP
            # print(alignment.cigar)  # TEMP
            # print(ref_seqs[name][start:end])  # TEMP
            # print(alignment.get_read_seq_by_ref_coords(start, end, ref_seq))  # TEMP

            # If there is no alignment, then the reference sequence is used for this part.
            if alignment is None:
                seq_parts.append(ref_seq[start:end])

            # If there is an alignment, then the sequence is taken from the read.
            else:
                seq_parts.append(alignment.get_read_seq_by_ref_coords(start, end, ref_seq))
        # print('\n')  # TEMP

        seq = ''.join(seq_parts)
        unpolished_sequences[name] = seq
    return unpolished_sequences


def polish_assembly_with_racon(names, unpolished_sequences, circularity, polish_reads, threads):
    polish_dir = 'temp_polish_' + str(os.getpid())
    if not os.path.isdir(polish_dir):
        os.makedirs(polish_dir)

    unitig_graph = UnitigGraph(names, unpolished_sequences, circularity)

    col_widths = [6, 12, 14]
    racon_table_header = ['Polish round', 'Assembly size', 'Mapping quality']
    print_table([racon_table_header], fixed_col_widths=col_widths, left_align_header=False,
                alignments='LRR', indent=0)

    best_fasta = None
    best_unitig_sequences = {}
    best_mapping_quality = 0
    times_quality_failed_to_beat_best = 0

    counter = itertools.count(start=1)
    current_fasta = os.path.join(polish_dir, ('%03d' % next(counter)) + '_unpolished_unitigs.fasta')
    unitig_graph.save_to_fasta(current_fasta)
    current_sequences = unpolished_sequences

    racon_loop_count = 10
    for polish_round_count in range(racon_loop_count):
        mappings_filename = os.path.join(polish_dir, ('%03d' % next(counter)) + '_alignments.paf')
        racon_log = os.path.join(polish_dir, ('%03d' % next(counter)) + '_racon.log')
        polished_fasta = os.path.join(polish_dir, ('%03d' % next(counter)) + '_polished.fasta')
        fixed_fasta = os.path.join(polish_dir, ('%03d' % next(counter)) + '_fixed.fasta')
        rotated_fasta = os.path.join(polish_dir, ('%03d' % next(counter)) + '_rotated.fasta')
        rotated_gfa = os.path.join(polish_dir, ('%03d' % next(counter)) + '_rotated.gfa')

        mapping_quality, unitig_depths = \
            make_racon_polish_alignments(current_fasta, mappings_filename, polish_reads, threads)
        for unitig_name, unitig_seg in unitig_graph.segments.items():
            if unitig_name in unitig_depths:
                unitig_seg.depth = unitig_depths[unitig_name]

        racon_table_row = ['begin' if polish_round_count == 0 else str(polish_round_count),
                           int_to_str(get_total_sequence_length(current_sequences)),
                           float_to_str(mapping_quality, 2)]
        print_table([racon_table_row], fixed_col_widths=col_widths, left_align_header=False,
                    alignments='LRR', indent=0, header_format='normal', bottom_align_header=False)

        # Do we have a new best?
        if mapping_quality > best_mapping_quality:
            best_mapping_quality = mapping_quality
            best_fasta = current_fasta
            best_unitig_sequences = {name: seq for name, seq in current_sequences.items()}
            times_quality_failed_to_beat_best = 0
        else:
            times_quality_failed_to_beat_best += 1

        # If we've failed to improve on our best quality for a few rounds, then we're done!
        if times_quality_failed_to_beat_best > 2:
            break

        # Run Racon. It crashes sometimes, so repeat until its return code is 0.
        command = ['racon', '--verbose', '9', '-t', str(threads), '--bq', '-1',
                   polish_reads, mappings_filename, current_fasta, polished_fasta]
        return_code = 1
        for _ in range(100):  # Only try a fixed number of times, to prevent an infinite loop.
            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = process.communicate()
            with open(racon_log, 'wb') as log_file:
                log_file.write(out)
                log_file.write(err)
            return_code = process.returncode
            if return_code == 0 and os.path.isfile(polished_fasta):
                break
            if os.path.isfile(polished_fasta):
                os.remove(polished_fasta)
            if os.path.isfile(racon_log):
                os.remove(racon_log)

        # If even after all those tries Racon still didn't succeed, then we give up!
        if return_code != 0 or not os.path.isfile(polished_fasta):
            break

        unitig_graph.replace_with_polished_sequences(polished_fasta)
        unitig_graph.save_to_fasta(fixed_fasta)
        unitig_graph.rotate_circular_sequences()
        unitig_graph.save_to_fasta(rotated_fasta)
        unitig_graph.save_to_gfa(rotated_gfa)
        current_fasta = rotated_fasta

    log.log('')
    if best_fasta:
        log.log('Best polish: ' + best_fasta)
        for unitig_name, unitig_seq in best_unitig_sequences.items():
            segment = unitig_graph.segments[unitig_name]
            segment.forward_sequence = unitig_seq
        unitig_graph.normalise_read_depths()
    else:
        log.log(red('Polishing failed'))


def make_racon_polish_alignments(current_fasta, mappings_filename, polish_reads, threads):
    mapping_quality = 0.0
    unitig_depths = collections.defaultdict(float)

    with open(mappings_filename, 'wt') as mappings:
        command = ['minimap2', '-c', '-x', 'map-ont', '-t', str(threads),
                   current_fasta, polish_reads]
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        while process.poll() is None:
            paf_line = process.stdout.readline().rstrip().decode()
            if paf_line:
                a = Alignment(paf_line)
                mapping_quality += a.percent_identity
                mappings.write(paf_line.split('cg:Z:')[0].rstrip())
                mappings.write('\n')
                unitig_depths[a.ref_name] += a.fraction_ref_aligned()
    return mapping_quality, unitig_depths


def get_total_sequence_length(sequences):
    return sum(len(x) for x in sequences.values())
