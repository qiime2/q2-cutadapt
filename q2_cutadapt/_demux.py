# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import subprocess

import yaml
from qiime2 import Artifact, Metadata
from q2_types.per_sample_sequences import (
    FastqManifestFormat,
    SingleLanePerSampleSingleEndFastqDirFmt,
    YamlFormat,
    FastqGzFormat,
)
from q2_types.multiplexed_sequences import (
    MultiplexedSingleEndBarcodeInSequenceDirFmt,
)

from ._format import CutadaptStatsDirFmt


def run_command(cmd, verbose=True):
    print('Running external command line application. This may print '
          'messages to stdout and/or stderr.')
    print('The command being run is below. This command cannot '
          'be manually re-run as it will depend on temporary files that '
          'no longer exist.')
    print('\nCommand:', end=' ')
    print(' '.join(cmd), end='\n\n')
    subprocess.run(cmd, check=True)


def _build_demux_command(seqs_dir_fmt, barcode_series, per_sample_dir_fmt,
                         stats_dir_fmt, untrimmed_dir_fmt):
    cmd = ['cutadapt']
    for (sample_id, barcode) in barcode_series.iteritems():
        cmd = cmd + ['-g', '%s=%s' % (sample_id, barcode)]
    cmd = cmd + [
        '-o', '%s/{name}.fastq.gz' % str(per_sample_dir_fmt),
        '--info-file', '%s/stats.tsv' % str(stats_dir_fmt),
        '--untrimmed-output', '%s/forward.fastq.gz' % str(untrimmed_dir_fmt),
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


def _write_metadata_yaml_in_results(per_sample_dir_fmt):
    metadata = YamlFormat()
    metadata.path.write_text(yaml.dump({'phred-offset': 33}))
    per_sample_dir_fmt.metadata.write_data(metadata, YamlFormat)


def _write_manifest_in_results(per_sample_dir_fmt):
    manifest = FastqManifestFormat()
    with manifest.open() as fh:
        fh.write('sample-id,filename,direction\n')
        filenames = per_sample_dir_fmt.sequences.iter_views(FastqGzFormat)
        for filename, _ in filenames:
            filename = str(filename)
            sample_id, _, _, _, _ = filename.rsplit('_', maxsplit=4)
            fh.write('%s,%s,forward\n' % (sample_id, filename))
    per_sample_dir_fmt.manifest.write_data(manifest, FastqManifestFormat)


def demux_single(ctx, seqs, barcodes):
    barcode_series = barcodes.to_series()
    seqs = seqs.view(MultiplexedSingleEndBarcodeInSequenceDirFmt)
    per_sample_sequences = SingleLanePerSampleSingleEndFastqDirFmt()
    untrimmed = MultiplexedSingleEndBarcodeInSequenceDirFmt()
    stats = CutadaptStatsDirFmt()

    cmd = _build_demux_command(seqs, barcode_series, per_sample_sequences,
                               stats, untrimmed)
    run_command(cmd)

    _rename_files(per_sample_sequences, barcode_series)
    _write_manifest_in_results(per_sample_sequences)
    _write_metadata_yaml_in_results(per_sample_sequences)

    per_sample_sequences = Artifact.import_data(
        'SampleData[SequencesWithQuality]', per_sample_sequences)
    untrimmed = Artifact.import_data(
        'MultiplexedSingleEndBarcodeInSequence', untrimmed)
    stats = Artifact.import_data(
        'CutadaptStats', stats)

    metadata_tabulate = ctx.get_action('metadata', 'tabulate')
    stats = stats.view(Metadata)
    stats_viz = metadata_tabulate(input=stats)

    return per_sample_sequences, untrimmed, stats_viz.visualization
