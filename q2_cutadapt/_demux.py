# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
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
                         untrimmed_dir_fmt, error_tolerance):
    cmd = ['cutadapt',
           '-g', 'file:%s' % barcode_fasta.name,
           '-e', str(error_tolerance),
           ]
    # PAIRED-END
    if isinstance(seqs_dir_fmt, MultiplexedPairedEndBarcodeInSequenceDirFmt):
        cmd = cmd + [
            # {name} is a cutadapt convention for interpolating the sample id
            # into the filename.
            '-o', '%s/{name}.1.fastq.gz' % str(per_sample_dir_fmt),
            '-p', '%s/{name}.2.fastq.gz' % str(per_sample_dir_fmt),
            '--untrimmed-output',
            '%s/forward.fastq.gz' % str(untrimmed_dir_fmt),
            '--untrimmed-paired-output',
            '%s/reverse.fastq.gz' % str(untrimmed_dir_fmt),
            str(seqs_dir_fmt.forward_sequences.view(FastqGzFormat)),
            str(seqs_dir_fmt.reverse_sequences.view(FastqGzFormat)),
            ]
    # SINGLE-END
    else:
        cmd = cmd + [
            # {name} is a cutadapt convention for interpolating the sample id
            # into the filename.
            '-o', '%s/{name}.fastq.gz' % str(per_sample_dir_fmt),
            '--untrimmed-output',
            '%s/forward.fastq.gz' % str(untrimmed_dir_fmt),
            str(seqs_dir_fmt.file.view(FastqGzFormat)),
            ]
    return cmd


def _rename_files(seqs_dir_fmt, per_sample_dir_fmt, barcode_series):
    # PAIRED-END
    if isinstance(seqs_dir_fmt, MultiplexedPairedEndBarcodeInSequenceDirFmt):
        for (sample_id, barcode_id) in barcode_series.iteritems():
            for read_direction in [1, 2]:
                out_fp = per_sample_dir_fmt.sequences.path_maker(
                    sample_id=sample_id, barcode_id=barcode_id,
                    lane_number=1, read_number=read_direction)
                src = os.path.join(str(per_sample_dir_fmt),
                                   '%s.%d.fastq.gz' % (sample_id,
                                                       read_direction))
                if os.path.isfile(src):
                    os.rename(src, str(out_fp))
    # SINGLE-END
    else:
        for (sample_id, barcode_id) in barcode_series.iteritems():
            out_fp = per_sample_dir_fmt.sequences.path_maker(
                sample_id=sample_id, barcode_id=barcode_id, lane_number=1,
                read_number=1)
            src = os.path.join(str(per_sample_dir_fmt),
                               '%s.fastq.gz' % sample_id)
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


def _demux(seqs, barcodes, error_tolerance, untrimmed):
    barcodes = barcodes.to_series()
    per_sample_sequences = CasavaOneEightSingleLanePerSampleDirFmt()

    _write_empty_fastq_to_mux_barcode_in_seq_fmt(untrimmed)

    with tempfile.NamedTemporaryFile() as barcode_fasta:
        _write_barcode_fasta(barcodes, barcode_fasta)
        cmd = _build_demux_command(seqs, barcode_fasta, per_sample_sequences,
                                   untrimmed, error_tolerance)
        run_command(cmd)

    _rename_files(seqs, per_sample_sequences, barcodes)
    muxed = len(list(per_sample_sequences.sequences.iter_views(FastqGzFormat)))
    if muxed == 0:
        raise ValueError('No samples were demultiplexed.')

    return per_sample_sequences, untrimmed


def demux_single(seqs: MultiplexedSingleEndBarcodeInSequenceDirFmt,
                 barcodes: qiime2.MetadataCategory,
                 error_tolerance: float=0.1) -> \
                    (CasavaOneEightSingleLanePerSampleDirFmt,
                     MultiplexedSingleEndBarcodeInSequenceDirFmt):

    untrimmed = MultiplexedSingleEndBarcodeInSequenceDirFmt()
    return _demux(seqs, barcodes, error_tolerance, untrimmed)


def demux_paired(seqs: MultiplexedPairedEndBarcodeInSequenceDirFmt,
                 forward_barcodes: qiime2.MetadataCategory,
                 error_tolerance: float=0.1) -> \
                    (CasavaOneEightSingleLanePerSampleDirFmt,
                     MultiplexedPairedEndBarcodeInSequenceDirFmt):

    untrimmed = MultiplexedPairedEndBarcodeInSequenceDirFmt()
    return _demux(seqs, forward_barcodes, error_tolerance, untrimmed)
