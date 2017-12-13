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
                         untrimmed_dir_fmt):
    cmd = ['cutadapt',
           '-g', 'file:%s' % barcode_fasta.name,
           # {name} is a cutadapt convention for interpolating the sample id
           # into the filename.
           '-o', '%s/{name}.fastq.gz' % str(per_sample_dir_fmt),
           '--untrimmed-output',
           '%s/forward.fastq.gz' % str(untrimmed_dir_fmt),
           str(seqs_dir_fmt.file.view(FastqGzFormat)),
           ]
    return cmd


def _rename_files(per_sample_dir_fmt, barcode_series):
    for (sample_id, barcode_id) in barcode_series.iteritems():
        out_fp = per_sample_dir_fmt.sequences.path_maker(sample_id=sample_id,
                                                         barcode_id=barcode_id,
                                                         lane_number=1,
                                                         read_number=1)
        src = os.path.join(str(per_sample_dir_fmt), '%s.fastq.gz' % sample_id)
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
    seqs_dir_fmt.file.write_data(fastq, FastqGzFormat)


def demux_single(seqs: MultiplexedSingleEndBarcodeInSequenceDirFmt,
                 barcodes: qiime2.MetadataCategory) -> \
                    (CasavaOneEightSingleLanePerSampleDirFmt,
                     MultiplexedSingleEndBarcodeInSequenceDirFmt):

    barcodes = barcodes.to_series()
    per_sample_sequences = CasavaOneEightSingleLanePerSampleDirFmt()
    untrimmed = MultiplexedSingleEndBarcodeInSequenceDirFmt()

    _write_empty_fastq_to_mux_barcode_in_seq_fmt(untrimmed)

    with tempfile.NamedTemporaryFile() as barcode_fasta:
        _write_barcode_fasta(barcodes, barcode_fasta)
        cmd = _build_demux_command(seqs, barcode_fasta, per_sample_sequences,
                                   untrimmed)
        run_command(cmd)

    _rename_files(per_sample_sequences, barcodes)
    muxed = len(list(per_sample_sequences.sequences.iter_views(FastqGzFormat)))
    if muxed == 0:
        raise ValueError('No samples were demultiplexed.')

    return per_sample_sequences, untrimmed
