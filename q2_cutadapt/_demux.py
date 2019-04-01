# ----------------------------------------------------------------------------
# Copyright (c) 2017-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
import os
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


def _build_demux_command(seqs_dir_fmt, barcode_fasta, per_sample_dir_fmt,
                         untrimmed_dir_fmt, error_rate):
    cmd = ['cutadapt',
           '--front', 'file:%s' % barcode_fasta.name,
           '--error-rate', str(error_rate),
           # {name} is a cutadapt convention for interpolating the sample id
           # into the filename.
           '-o', os.path.join(str(per_sample_dir_fmt), '{name}.1.fastq.gz'),
           '--untrimmed-output',
           os.path.join(str(untrimmed_dir_fmt), 'forward.fastq.gz'),
           ]
    if isinstance(seqs_dir_fmt, MultiplexedPairedEndBarcodeInSequenceDirFmt):
        # PAIRED-END
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
            if os.path.isfile(src):
                os.rename(src, str(out_fp))


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


def _demux(seqs, mux_fmt, barcodes, error_tolerance, batch_size):
    barcodes = barcodes.to_series()
    per_sample_sequences = CasavaOneEightSingleLanePerSampleDirFmt()
    n_samples = len(barcodes)
    batch_size = n_samples if batch_size == 0 else batch_size
    batches = np.arange(n_samples) // batch_size
    previous_untrimmed = seqs
    for _, barcode_batch in barcodes.groupby(batches):
        current_untrimmed = mux_fmt()
        _write_empty_fastq_to_mux_barcode_in_seq_fmt(current_untrimmed)
        with tempfile.NamedTemporaryFile() as barcode_fasta:
            _write_barcode_fasta(barcode_batch, barcode_fasta)
            cmd = _build_demux_command(previous_untrimmed, barcode_fasta,
                                       per_sample_sequences,
                                       current_untrimmed, error_tolerance)
            run_command(cmd)
        previous_untrimmed = current_untrimmed
    _rename_files(seqs, per_sample_sequences, barcodes)
    muxed = len(list(per_sample_sequences.sequences.iter_views(FastqGzFormat)))
    if muxed == 0:
        raise ValueError('No samples were demultiplexed.')
    return per_sample_sequences, previous_untrimmed


def demux_single(seqs: MultiplexedSingleEndBarcodeInSequenceDirFmt,
                 barcodes: qiime2.CategoricalMetadataColumn,
                 error_rate: float = 0.1,
                 batch_size: int = 0) -> \
                    (CasavaOneEightSingleLanePerSampleDirFmt,
                     MultiplexedSingleEndBarcodeInSequenceDirFmt):
    mux_fmt = MultiplexedSingleEndBarcodeInSequenceDirFmt
    return _demux(seqs, mux_fmt, barcodes, error_rate, batch_size)


def demux_paired(seqs: MultiplexedPairedEndBarcodeInSequenceDirFmt,
                 forward_barcodes: qiime2.CategoricalMetadataColumn,
                 error_rate: float = 0.1,
                 batch_size: int = 0) -> \
                    (CasavaOneEightSingleLanePerSampleDirFmt,
                     MultiplexedPairedEndBarcodeInSequenceDirFmt):
    mux_fmt = MultiplexedPairedEndBarcodeInSequenceDirFmt
    return _demux(seqs, mux_fmt, forward_barcodes, error_rate, batch_size)
