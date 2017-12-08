# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
import os
import pathlib
import shutil
import tempfile
import unittest

import pandas as pd

from q2_cutadapt._demux import (_build_demux_command, _rename_files,
                                _write_metadata_yaml_in_results,
                                _write_manifest_in_results)
from q2_cutadapt import CutadaptStatsDirFmt
from q2_types.multiplexed_sequences import (
    MultiplexedSingleEndBarcodeInSequenceDirFmt)
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt, FastqGzFormat)
from qiime2 import Artifact, MetadataCategory
from qiime2.util import redirected_stdio
from qiime2.plugin.testing import TestPluginBase


class TestDemuxSingle(TestPluginBase):
    package = 'q2_cutadapt.tests'

    def setUp(self):
        super().setUp()
        self.demux_single_fn = self.plugin.pipelines['demux_single']

        muxed_sequences_fp = self.get_data_path('forward.fastq.gz')
        self.muxed_sequences = Artifact.import_data(
            'MultiplexedSingleEndBarcodeInSequence', muxed_sequences_fp)

        self.metadata = MetadataCategory(pd.Series(['AAAA', 'CCCC'],
                                         index=['sample_a', 'sample_b'],
                                         name='Barcode'))

    def test_positive(self):
        with redirected_stdio(stderr=os.devnull):
            obs_demuxed, obs_untrimmed, obs_stats = self.demux_single_fn(
                self.muxed_sequences, self.metadata)

        obs_demuxed = obs_demuxed.view(SingleLanePerSampleSingleEndFastqDirFmt)
        exp_samples_and_barcodes = self.metadata.to_series().iteritems()
        obs_demuxed_seqs = obs_demuxed.sequences.iter_views(FastqGzFormat)
        zipped = zip(exp_samples_and_barcodes, obs_demuxed_seqs)
        for (sample_id, barcode), (filename, _) in zipped:
            filename = str(filename)
            self.assertTrue(sample_id in filename)
            self.assertTrue(barcode in filename)

        obs_untrimmed = obs_untrimmed.view(
            MultiplexedSingleEndBarcodeInSequenceDirFmt)
        obs_untrimmed = obs_untrimmed.file.view(FastqGzFormat)
        obs_untrimmed = gzip.decompress(obs_untrimmed.path.read_bytes())
        self.assertEqual(b'@id6\nGGGGACGTACGT\n+\nzzzzzzzzzzzz\n',
                         obs_untrimmed)

        with tempfile.TemporaryDirectory() as viz_dir:
            obs_stats.export_data(viz_dir)
            index_fp = os.path.join(viz_dir, 'index.html')
            index = open(index_fp, 'r').read()
            for id_ in ['id1', 'id2', 'id3', 'id4', 'id5', 'id6']:
                self.assertTrue(id_ in index)


class TestDemuxUtils(TestPluginBase):
    package = 'q2_cutadapt.tests'

    def setUp(self):
        super().setUp()

        self.fastq_fp = \
            self.get_data_path('forward.fastq.gz')
        self.fastq = FastqGzFormat(self.fastq_fp, mode='r')
        self.seqs_dir_fmt = MultiplexedSingleEndBarcodeInSequenceDirFmt()
        self.seqs_dir_fmt.file.write_data(self.fastq, FastqGzFormat)
        self.barcode_series = pd.Series(['A', 'G'],
                                        index=['sample_a', 'sample_b'])
        self.per_sample_dir_fmt = SingleLanePerSampleSingleEndFastqDirFmt()
        self.untrimmed_dir_fmt = MultiplexedSingleEndBarcodeInSequenceDirFmt()
        self.stats_dir_fmt = CutadaptStatsDirFmt()

    def test_build_demux_command(self):
        obs = _build_demux_command(self.seqs_dir_fmt, self.barcode_series,
                                   self.per_sample_dir_fmt, self.stats_dir_fmt,
                                   self.untrimmed_dir_fmt)
        exp_map = ['cutadapt', '-g', 'sample_a=A', '-g', 'sample_b=G']
        obs_map = obs[:5]

        self.assertEqual(exp_map, obs_map)

        self.assertTrue(str(self.per_sample_dir_fmt) in obs[6])
        self.assertTrue(str(self.stats_dir_fmt) in obs[8])
        self.assertTrue(str(self.untrimmed_dir_fmt) in obs[10])
        self.assertEqual(str(self.seqs_dir_fmt.file.view(FastqGzFormat)),
                         obs[11])

    def test_rename_files(self):
        for fn in ['sample_a.fastq.gz', 'sample_b.fastq.gz']:
            shutil.copy(self.fastq_fp,
                        str(self.per_sample_dir_fmt.path / pathlib.Path(fn)))

        _rename_files(self.per_sample_dir_fmt, self.barcode_series)

        seqs = self.per_sample_dir_fmt.sequences.iter_views(FastqGzFormat)
        for fn, (sample_id, barcode) in zip(seqs,
                                            self.barcode_series.iteritems()):
            self.assertTrue(sample_id in str(fn))
            self.assertTrue(barcode in str(fn))

    def test_rename_files_extra_samples_in_barcode_map(self):
        barcode_series = pd.Series(['A', 'G', 'C'],
                                   index=['sample_a', 'sample_b', 'sample_c'])

        for fn in ['sample_a.fastq.gz', 'sample_b.fastq.gz']:
            shutil.copy(self.fastq_fp,
                        str(self.per_sample_dir_fmt.path / pathlib.Path(fn)))

        _rename_files(self.per_sample_dir_fmt, barcode_series)

        seqs = self.per_sample_dir_fmt.sequences.iter_views(FastqGzFormat)
        for fn, (sample_id, barcode) in zip(seqs, barcode_series.iteritems()):
            self.assertTrue(sample_id in str(fn))
            self.assertTrue(barcode in str(fn))

    def test_write_metadata_yaml_in_results(self):
        _write_metadata_yaml_in_results(self.per_sample_dir_fmt)
        yml = os.path.join(str(self.per_sample_dir_fmt), 'metadata.yml')
        self.assertTrue(os.path.isfile(yml))
        yml_content = open(yml, 'r').read()
        self.assertEqual('{phred-offset: 33}\n', yml_content)

    def test_write_manifest_in_results(self):
        for fn in ['sample_a_A_L001_R1_001.fastq.gz',
                   'sample_b_G_L001_R1_001.fastq.gz']:
            shutil.copy(self.fastq_fp,
                        str(self.per_sample_dir_fmt.path / pathlib.Path(fn)))

        _write_manifest_in_results(self.per_sample_dir_fmt)

        manifest = os.path.join(str(self.per_sample_dir_fmt), 'MANIFEST')
        self.assertTrue(os.path.isfile(manifest))
        manifest_content = open(manifest, 'r').read()
        self.assertEqual('sample-id,filename,direction\n'
                         'sample_a,sample_a_A_L001_R1_001.fastq.gz,forward\n'
                         'sample_b,sample_b_G_L001_R1_001.fastq.gz,forward\n',
                         manifest_content)


if __name__ == '__main__':
    unittest.main()
