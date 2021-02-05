# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
import os
import shutil
import subprocess
import tempfile

import qiime2
from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt,
    FastqGzFormat,
)
from q2_types.multiplexed_sequences import (
    MultiplexedSingleEndBarcodeInSequenceDirFmt,
    MultiplexedPairedEndBarcodeInSequenceDirFmt,
)

import pandas as pd
import numpy as np


def run_command(cmd, verbose=True):
    print('Running external command line application. This may print '
          'messages to stdout and/or stderr.')
    print('The command being run is below. This command cannot '
          'be manually re-run as it will depend on temporary files that '
          'no longer exist.')
    print('\nCommand:', end=' ')
    print(' '.join(cmd), end='\n\n')
    subprocess.run(cmd, check=True)


def _build_demux_command(seqs_dir_fmt, barcode_fhs, per_sample_dir_fmt,
                         untrimmed_dir_fmt, error_rate, minimum_length):
    cmd = ['cutadapt',
           '--front', 'file:%s' % barcode_fhs['fwd'].name,
           '--error-rate', str(error_rate),
           '--minimum-length', str(minimum_length),
           # {name} is a cutadapt convention for interpolating the sample id
           # into the filename.
           '-o', os.path.join(str(per_sample_dir_fmt), '{name}.1.fastq.gz'),
           '--untrimmed-output',
           os.path.join(str(untrimmed_dir_fmt), 'forward.fastq.gz'),
           ]
    if isinstance(seqs_dir_fmt, MultiplexedPairedEndBarcodeInSequenceDirFmt):
        # PAIRED-END
        if barcode_fhs['rev'] is not None:
            # Dual indices
            cmd += [
                '--pair-adapters',
                '-G', 'file:%s' % barcode_fhs['rev'].name,
            ]
        cmd += [
            '-p', os.path.join(str(per_sample_dir_fmt), '{name}.2.fastq.gz'),
            '--untrimmed-paired-output',
            os.path.join(str(untrimmed_dir_fmt), 'reverse.fastq.gz'),
            str(seqs_dir_fmt.forward_sequences.view(FastqGzFormat)),
            str(seqs_dir_fmt.reverse_sequences.view(FastqGzFormat)),
            ]
    else:
        # SINGLE-END
        cmd += [str(seqs_dir_fmt.file.view(FastqGzFormat))]
    return cmd


def _rename_files(seqs_dir_fmt, per_sample_dir_fmt, barcode_series):
    read_directions = [1]
    if isinstance(seqs_dir_fmt, MultiplexedPairedEndBarcodeInSequenceDirFmt):
        # PAIRED-END
        read_directions.append(2)

    for (sample_id, barcode_id) in barcode_series.iteritems():
        for read_direction in read_directions:
            out_fp = per_sample_dir_fmt.sequences.path_maker(
                sample_id=sample_id, barcode_id=barcode_id,
                lane_number=1, read_number=read_direction)
            src = os.path.join(str(per_sample_dir_fmt),
                               '%s.%d.fastq.gz' % (sample_id,
                                                   read_direction))

            # TODO: remove this outer guard when we upgrade to cutadapt 3
            if os.path.isfile(src):
                if out_fp.exists():
                    _merge_files(src, str(out_fp))
                    os.remove(src)
                else:
                    os.rename(src, str(out_fp))


def _merge_files(src, dst):
    with gzip.open(src, mode='rt', encoding='ascii') as src_fh, \
            gzip.open(dst, mode='at', encoding='ascii') as dst_fh:
        shutil.copyfileobj(src_fh, dst_fh)


def _write_barcode_fasta(barcode_series, barcode_fasta):
    with open(barcode_fasta.name, 'w') as fh:
        for (sample_id, barcode) in barcode_series.iteritems():
            fh.write('>%s\n%s\n' % (sample_id, barcode))


def _write_empty_fastq_to_mux_barcode_in_seq_fmt(seqs_dir_fmt):
    fastq = FastqGzFormat()
    with gzip.open(str(fastq), 'w') as fh:
        fh.write(b'')
    # PAIRED-END
    if isinstance(seqs_dir_fmt, MultiplexedPairedEndBarcodeInSequenceDirFmt):
        seqs_dir_fmt.forward_sequences.write_data(fastq, FastqGzFormat)
        seqs_dir_fmt.reverse_sequences.write_data(fastq, FastqGzFormat)
    # SINGLE-END
    else:
        seqs_dir_fmt.file.write_data(fastq, FastqGzFormat)


def _demux(seqs, per_sample_sequences, forward_barcodes, reverse_barcodes,
           error_tolerance, mux_fmt, batch_size, minimum_length):
    fwd_barcode_name = forward_barcodes.name
    forward_barcodes = forward_barcodes.drop_missing_values()
    barcodes = forward_barcodes.to_series().to_frame()
    if reverse_barcodes is not None:
        barcode_pairs = set()
        samples_w_missing_barcodes = set()
        samples_w_dup_barcode_pairs = set()
        rev_barcode_name = reverse_barcodes.name
        rev_barcodes = reverse_barcodes.to_series()
        # 'sort = false' below prevents a warning about future behavior changes
        # by selecting the future behavior explicitly
        barcodes = pd.concat([barcodes, rev_barcodes], axis=1, sort=False)

        for sample_id, f_barcode, r_barcode in barcodes.itertuples():
            if pd.isnull(f_barcode) or pd.isnull(r_barcode):
                samples_w_missing_barcodes.add(sample_id)
            if (f_barcode, r_barcode) in barcode_pairs:
                samples_w_dup_barcode_pairs.add(sample_id)
            barcode_pairs.add((f_barcode, r_barcode))

        if samples_w_missing_barcodes:
            raise ValueError('The following samples do not have both '
                             'forward and reverse barcodes (note: if your '
                             'reads are in single index mixed orientation, '
                             'try again with all of your barcodes in a single '
                             'metadata column): %s'
                             % ', '.join(sorted(samples_w_missing_barcodes)))
        if samples_w_dup_barcode_pairs:
            raise ValueError('The following samples have duplicate barcode'
                             ' pairs: %s' %
                             ', '.join(sorted(samples_w_dup_barcode_pairs)))

    n_samples = len(barcodes)
    if batch_size > n_samples:
        raise ValueError('The batch_size (%d) cannot be greater than the '
                         'number of samples (%d).' % (
                             batch_size, n_samples))
    batch_size = n_samples if batch_size == 0 else batch_size
    batches = np.arange(n_samples) // batch_size
    previous_untrimmed = seqs
    for _, barcode_batch in barcodes.groupby(batches):
        current_untrimmed = mux_fmt()
        _write_empty_fastq_to_mux_barcode_in_seq_fmt(current_untrimmed)
        open_fhs = {'fwd': tempfile.NamedTemporaryFile(), 'rev': None}
        _write_barcode_fasta(barcode_batch[fwd_barcode_name], open_fhs['fwd'])
        if reverse_barcodes is not None:
            open_fhs['rev'] = tempfile.NamedTemporaryFile()
            _write_barcode_fasta(barcode_batch[rev_barcode_name],
                                 open_fhs['rev'])
        cmd = _build_demux_command(previous_untrimmed, open_fhs,
                                   per_sample_sequences,
                                   current_untrimmed, error_tolerance,
                                   minimum_length)
        run_command(cmd)
        open_fhs['fwd'].close()
        if reverse_barcodes is not None:
            open_fhs['rev'].close()
        previous_untrimmed = current_untrimmed

    # Only use the forward barcode in the renamed files
    _rename_files(seqs, per_sample_sequences, barcodes[fwd_barcode_name])
    muxed = len(list(per_sample_sequences.sequences.iter_views(FastqGzFormat)))
    if muxed == 0:
        raise ValueError('No samples were demultiplexed.')
    return previous_untrimmed


def demux_single(seqs: MultiplexedSingleEndBarcodeInSequenceDirFmt,
                 barcodes: qiime2.CategoricalMetadataColumn,
                 error_rate: float = 0.1,
                 batch_size: int = 0,
                 minimum_length: int = 1) -> \
                    (CasavaOneEightSingleLanePerSampleDirFmt,
                     MultiplexedSingleEndBarcodeInSequenceDirFmt):
    per_sample_sequences = CasavaOneEightSingleLanePerSampleDirFmt()
    mux_fmt = MultiplexedSingleEndBarcodeInSequenceDirFmt

    untrimmed = _demux(
        seqs, per_sample_sequences, barcodes, None, error_rate, mux_fmt,
        batch_size, minimum_length)

    return per_sample_sequences, untrimmed


def demux_paired(seqs: MultiplexedPairedEndBarcodeInSequenceDirFmt,
                 forward_barcodes: qiime2.CategoricalMetadataColumn,
                 reverse_barcodes: qiime2.CategoricalMetadataColumn = None,
                 error_rate: float = 0.1,
                 batch_size: int = 0,
                 minimum_length: int = 1,
                 mixed_orientation: bool = False) -> \
                    (CasavaOneEightSingleLanePerSampleDirFmt,
                     MultiplexedPairedEndBarcodeInSequenceDirFmt):
    if mixed_orientation and reverse_barcodes is not None:
        raise ValueError('Dual-indexed barcodes for mixed orientation '
                         'reads are not supported.')

    per_sample_sequences = CasavaOneEightSingleLanePerSampleDirFmt()
    mux_fmt = MultiplexedPairedEndBarcodeInSequenceDirFmt

    untrimmed = _demux(
        seqs, per_sample_sequences, forward_barcodes, reverse_barcodes,
        error_rate, mux_fmt, batch_size, minimum_length)

    if mixed_orientation:
        fwd = untrimmed.forward_sequences.view(FastqGzFormat)
        rev = untrimmed.reverse_sequences.view(FastqGzFormat)

        remaining_seqs = MultiplexedPairedEndBarcodeInSequenceDirFmt()
        # fwd -> rev && rev -> fwd
        remaining_seqs.forward_sequences.write_data(rev, FastqGzFormat)
        remaining_seqs.reverse_sequences.write_data(fwd, FastqGzFormat)

        untrimmed = _demux(
            remaining_seqs, per_sample_sequences, forward_barcodes,
            reverse_barcodes, error_rate, mux_fmt, batch_size,
            minimum_length)

    return per_sample_sequences, untrimmed
