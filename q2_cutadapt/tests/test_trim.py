# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
import itertools
import os
import unittest

import pandas as pd

from q2_cutadapt._trim import _build_trim_command
from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    FastqGzFormat,
)
from qiime2 import Artifact
from qiime2.util import redirected_stdio
from qiime2.plugin.testing import TestPluginBase


class TestTrimSingle(TestPluginBase):
    package = 'q2_cutadapt.tests'

    # This test is really just to make sure that the command runs - the
    # detailed tests in the Util Tests below ensure the commands are crafted
    # appropriately.
    def test_typical(self):
        demuxed_art = Artifact.import_data('SampleData[SequencesWithQuality]',
                                           self.get_data_path('single-end'))
        adapter = ['TACGGAGGATCC']
        with redirected_stdio(stdout=os.devnull):
            obs_art, = self.plugin.methods['trim_single'](demuxed_art,
                                                          front=adapter)
        demuxed = demuxed_art.view(SingleLanePerSampleSingleEndFastqDirFmt)
        demuxed_seqs = demuxed.sequences.iter_views(FastqGzFormat)
        obs = obs_art.view(SingleLanePerSampleSingleEndFastqDirFmt)
        obs_seqs = obs.sequences.iter_views(FastqGzFormat)
        # Iterate over each sample, side-by-side
        for (_, exp_fp), (_, obs_fp) in zip(demuxed_seqs, obs_seqs):
            exp_fh = gzip.open(str(exp_fp), 'rt')
            obs_fh = gzip.open(str(obs_fp), 'rt')
            # Iterate over expected and observed reads, side-by-side
            for records in itertools.zip_longest(*[exp_fh] * 4, *[obs_fh] * 4):
                (exp_seq_h, exp_seq, _, exp_qual,
                 obs_seq_h, obs_seq, _, obs_qual) = records
                # Make sure cutadapt hasn't shuffled the read order
                self.assertEqual(exp_seq_h, obs_seq_h)
                self.assertTrue(obs_seq in exp_seq)
                # The adapter should not be present in the trimmed seqs
                self.assertTrue('TACGGAGGATCC' not in obs_seq)
                self.assertTrue(obs_qual in exp_qual)
                # Make sure cutadapt trimmed the quality scores, too
                self.assertEqual(len(obs_seq), len(obs_qual))
            exp_fh.close(), obs_fh.close()

    def test_min_length(self):
        demuxed_art = Artifact.import_data('SampleData[SequencesWithQuality]',
                                           self.get_data_path('single-end'))
        # The following "adapter" has been picked specifically to remove
        # the entire sequence with the ID @HWI-EAS440_0386:1:28:6491:1375#0/1.
        adapter = ['GGGGGGATCGGGGGCG']
        empty_seq_id = '@HWI-EAS440_0386:1:28:6491:1375#0/1'

        with redirected_stdio(stdout=os.devnull):
            obs_art, = self.plugin.methods['trim_single'](demuxed_art,
                                                          adapter=adapter)
        obs = obs_art.view(SingleLanePerSampleSingleEndFastqDirFmt)
        for _, obs_fp in obs.sequences.iter_views(FastqGzFormat):
            with gzip.open(str(obs_fp), 'rt') as obs_fh:
                for record in itertools.zip_longest(*[obs_fh] * 4):
                    self.assertTrue(record[0] != empty_seq_id)


class TestTrimPaired(TestPluginBase):
    package = 'q2_cutadapt.tests'

    # This test is really just to make sure that the command runs - the
    # detailed tests in the Util Tests below ensure the commands are crafted
    # appropriately.
    def test_typical(self):
        demuxed_art = Artifact.import_data(
            'SampleData[PairedEndSequencesWithQuality]',
            self.get_data_path('paired-end'))
        adapter = ['TACGGAGGATCC']
        with redirected_stdio(stdout=os.devnull):
            # The forward and reverse reads are identical in these data
            obs_art, = self.plugin.methods['trim_paired'](demuxed_art,
                                                          front_f=adapter,
                                                          front_r=adapter)
        demuxed = demuxed_art.view(SingleLanePerSampleSingleEndFastqDirFmt)
        demuxed_seqs = demuxed.sequences.iter_views(FastqGzFormat)
        obs = obs_art.view(SingleLanePerSampleSingleEndFastqDirFmt)
        obs_seqs = obs.sequences.iter_views(FastqGzFormat)
        # Iterate over each sample, side-by-side
        for (_, exp_fp), (_, obs_fp) in zip(demuxed_seqs, obs_seqs):
            exp_fh = gzip.open(str(exp_fp), 'rt')
            obs_fh = gzip.open(str(obs_fp), 'rt')
            # Iterate over expected and observed reads, side-by-side
            for records in itertools.zip_longest(*[exp_fh] * 4, *[obs_fh] * 4):
                (exp_seq_h, exp_seq, _, exp_qual,
                 obs_seq_h, obs_seq, _, obs_qual) = records
                # Make sure cutadapt hasn't shuffled the read order
                self.assertEqual(exp_seq_h, obs_seq_h)
                self.assertTrue(obs_seq in exp_seq)
                # The adapter should not be present in the trimmed seqs
                self.assertTrue('TACGGAGGATCC' not in obs_seq)
                self.assertTrue(obs_qual in exp_qual)
                # Make sure cutadapt trimmed the quality scores, too
                self.assertEqual(len(obs_seq), len(obs_qual))
            exp_fh.close(), obs_fh.close()

    def test_unordered(self):
        demuxed_art = Artifact.import_data(
            'SampleData[PairedEndSequencesWithQuality]',
            self.get_data_path('paired-end-unordered'))
        with redirected_stdio(stdout=os.devnull):
            # The forward and reverse reads are identical in these data
            obs_art, = self.plugin.methods['trim_paired'](demuxed_art,
                                                          front_f=['TTTT'],
                                                          front_r=['AAAA'])
        demuxed = demuxed_art.view(SingleLanePerSampleSingleEndFastqDirFmt)
        demuxed_seqs = demuxed.sequences.iter_views(FastqGzFormat)
        obs = obs_art.view(SingleLanePerSampleSingleEndFastqDirFmt)
        obs_seqs = obs.sequences.iter_views(FastqGzFormat)
        # Iterate over each sample, side-by-side
        for (_, exp_fp), (_, obs_fp) in zip(demuxed_seqs, obs_seqs):
            exp_fh = gzip.open(str(exp_fp), 'rt')
            obs_fh = gzip.open(str(obs_fp), 'rt')
            # Iterate over expected and observed reads, side-by-side
            for records in itertools.zip_longest(*[exp_fh] * 4, *[obs_fh] * 4):
                (exp_seq_h, exp_seq, _, exp_qual,
                 obs_seq_h, obs_seq, _, obs_qual) = records
                # The adapter should not be present in the trimmed seqs
                if 'R1_001.fastq' in str(obs_fp):
                    self.assertNotIn('TTTT', obs_seq)
                else:
                    self.assertNotIn('AAAA', obs_seq)
                self.assertTrue(obs_qual in exp_qual)
                # Make sure cutadapt trimmed the quality scores, too
                self.assertEqual(len(obs_seq), len(obs_qual))
            exp_fh.close(), obs_fh.close()


class TestTrimUtilsSingle(TestPluginBase):
    package = 'q2_cutadapt.tests'

    def setUp(self):
        super().setUp()

        self.demux_seqs = SingleLanePerSampleSingleEndFastqDirFmt(
            self.get_data_path('single-end'), mode='r')
        self.trimmed_seqs = CasavaOneEightSingleLanePerSampleDirFmt()

    def test_build_trim_command_typical(self):
        df = self.demux_seqs.manifest.view(pd.DataFrame)
        for _, fwd in df.itertuples():
            obs = _build_trim_command(fwd, None,
                                      self.trimmed_seqs,
                                      cores=0,
                                      adapter_f=['AAAA'],
                                      front_f=['GGGG'],
                                      anywhere_f=['CCCC'],
                                      error_rate=2,
                                      indels=False,
                                      times=3,
                                      overlap=4,
                                      match_read_wildcards=True,
                                      match_adapter_wildcards=False,
                                      minimum_length=2,
                                      discard_untrimmed=True)
            obs = ' '.join(obs)

            self.assertTrue('-o %s' % str(self.trimmed_seqs.path / fwd[0])
                            in obs)
            self.assertTrue('--cores 0' in obs)
            self.assertTrue('--adapter AAAA' in obs)
            self.assertTrue('--front GGGG' in obs)
            self.assertTrue('--anywhere CCCC' in obs)
            self.assertTrue('--error-rate 2' in obs)
            self.assertTrue('--times 3' in obs)
            self.assertTrue('--overlap 4' in obs)
            self.assertTrue('--no-indels' in obs)
            self.assertTrue('--match-read-wildcards' in obs)
            self.assertTrue('--no-match-adapter-wildcards' in obs)
            self.assertTrue('--minimum-length 2' in obs)
            self.assertTrue('--discard-untrimmed' in obs)

            self.assertTrue(str(self.demux_seqs) in obs)

    def test_build_trim_command_multiple_adapters(self):
        df = self.demux_seqs.manifest.view(pd.DataFrame)
        for _, fwd in df.itertuples():
            obs = _build_trim_command(fwd, None,
                                      self.trimmed_seqs,
                                      adapter_f=['AAAA', 'GGGG', 'CCCC'])
            obs = ' '.join(obs)

            self.assertTrue('--adapter AAAA' in obs)
            self.assertTrue('--adapter GGGG' in obs)
            self.assertTrue('--adapter CCCC' in obs)
            self.assertTrue('--front' not in obs)
            self.assertTrue('--anywhere' not in obs)

    def test_build_trim_command_no_adapters_or_flags(self):
        df = self.demux_seqs.manifest.view(pd.DataFrame)
        for _, fwd in df.itertuples():
            obs = _build_trim_command(fwd, None,
                                      self.trimmed_seqs)
            obs = ' '.join(obs)

            self.assertTrue('--adapter' not in obs)
            self.assertTrue('--front' not in obs)
            self.assertTrue('--anywhere' not in obs)
            self.assertTrue('--no-indels' not in obs)
            self.assertTrue('--match-read-wildcards' not in obs)
            self.assertTrue('--no-match-adapter-wildcards' not in obs)
            self.assertTrue('--minimum-length 1' in obs)
            self.assertTrue('--discard-untrimmed' not in obs)


class TestTrimUtilsPaired(TestPluginBase):
    package = 'q2_cutadapt.tests'

    def setUp(self):
        super().setUp()

        self.demux_seqs = SingleLanePerSamplePairedEndFastqDirFmt(
            self.get_data_path('paired-end'), mode='r')
        self.trimmed_seqs = CasavaOneEightSingleLanePerSampleDirFmt()

    def test_build_trim_command_typical(self):
        df = self.demux_seqs.manifest.view(pd.DataFrame)
        for _, fwd, rev in df.itertuples():
            obs = _build_trim_command(fwd, rev,
                                      self.trimmed_seqs,
                                      cores=0,
                                      adapter_f=['AAAA'],
                                      front_f=['GGGG'],
                                      anywhere_f=['CCCC'],
                                      adapter_r=['TTTT'],
                                      front_r=['CCCC'],
                                      anywhere_r=['GGGG'],
                                      error_rate=2,
                                      indels=False,
                                      times=3,
                                      overlap=4,
                                      match_read_wildcards=True,
                                      match_adapter_wildcards=False,
                                      minimum_length=2,
                                      discard_untrimmed=True)
            obs = ' '.join(obs)

            self.assertTrue('-o %s' % str(self.trimmed_seqs.path / fwd[0])
                            in obs)
            self.assertTrue('-p %s' % str(self.trimmed_seqs.path / rev[0])
                            in obs)
            self.assertTrue('--cores 0' in obs)
            self.assertTrue('--adapter AAAA' in obs)
            self.assertTrue('--front GGGG' in obs)
            self.assertTrue('--anywhere CCCC' in obs)
            self.assertTrue('-A TTTT' in obs)
            self.assertTrue('-G CCCC' in obs)
            self.assertTrue('-B GGGG' in obs)
            self.assertTrue('--error-rate 2' in obs)
            self.assertTrue('--times 3' in obs)
            self.assertTrue('--overlap 4' in obs)
            self.assertTrue('--no-indels' in obs)
            self.assertTrue('--match-read-wildcards' in obs)
            self.assertTrue('--no-match-adapter-wildcards' in obs)
            self.assertTrue('--minimum-length 2' in obs)
            self.assertTrue('--discard-untrimmed' in obs)

            self.assertTrue(str(self.demux_seqs) in obs)

    def test_build_trim_command_multiple_adapters(self):
        df = self.demux_seqs.manifest.view(pd.DataFrame)
        for _, fwd, rev in df.itertuples():
            obs = _build_trim_command(fwd, rev, self.trimmed_seqs,
                                      adapter_f=['AAAA', 'GGGG', 'CCCC'],
                                      adapter_r=['TTTT', 'CCCC', 'GGGG'])
            obs = ' '.join(obs)

            self.assertTrue('--adapter AAAA' in obs)
            self.assertTrue('--adapter GGGG' in obs)
            self.assertTrue('--adapter CCCC' in obs)
            self.assertTrue('-A TTTT' in obs)
            self.assertTrue('-A CCCC' in obs)
            self.assertTrue('-A GGGG' in obs)

            self.assertTrue('--front' not in obs)
            self.assertTrue('--anywhere' not in obs)
            self.assertTrue('-G' not in obs)
            self.assertTrue('-B' not in obs)
            self.assertTrue('--discard-trimmed' not in obs)


if __name__ == '__main__':
    unittest.main()
