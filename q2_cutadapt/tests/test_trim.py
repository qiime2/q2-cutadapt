# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
import itertools
import os
import unittest

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


class TestTrimUtilsSingle(TestPluginBase):
    package = 'q2_cutadapt.tests'

    def setUp(self):
        super().setUp()

        self.demux_seqs = SingleLanePerSampleSingleEndFastqDirFmt(
            self.get_data_path('single-end'), mode='r')
        self.trimmed_seqs = CasavaOneEightSingleLanePerSampleDirFmt()

    def test_build_trim_command_typical(self):
        for fwd in self.demux_seqs.sequences.iter_views(FastqGzFormat):
            obs = _build_trim_command(self.demux_seqs, fwd[0], None,
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
                                      match_adapter_wildcards=False)
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

            self.assertTrue(str(self.demux_seqs) in obs)

    def test_build_trim_command_multiple_adapters(self):
        for fwd in self.demux_seqs.sequences.iter_views(FastqGzFormat):
            obs = _build_trim_command(self.demux_seqs, fwd[0], None,
                                      self.trimmed_seqs,
                                      adapter_f=['AAAA', 'GGGG', 'CCCC'])
            obs = ' '.join(obs)

            self.assertTrue('--adapter AAAA' in obs)
            self.assertTrue('--adapter GGGG' in obs)
            self.assertTrue('--adapter CCCC' in obs)
            self.assertTrue('--front' not in obs)
            self.assertTrue('--anywhere' not in obs)

    def test_build_trim_command_no_adapters_or_flags(self):
        for fwd in self.demux_seqs.sequences.iter_views(FastqGzFormat):
            obs = _build_trim_command(self.demux_seqs, fwd[0], None,
                                      self.trimmed_seqs)
            obs = ' '.join(obs)

            self.assertTrue('--adapter' not in obs)
            self.assertTrue('--front' not in obs)
            self.assertTrue('--anywhere' not in obs)
            self.assertTrue('--no-indels' not in obs)
            self.assertTrue('--match-read-wildcards' not in obs)
            self.assertTrue('--no-match-adapter-wildcards' not in obs)


class TestTrimUtilsPaired(TestPluginBase):
    package = 'q2_cutadapt.tests'

    def setUp(self):
        super().setUp()

        self.demux_seqs = SingleLanePerSamplePairedEndFastqDirFmt(
            self.get_data_path('paired-end'), mode='r')
        self.trimmed_seqs = CasavaOneEightSingleLanePerSampleDirFmt()

    def test_build_trim_command_typical(self):
        itr = self.demux_seqs.sequences.iter_views(FastqGzFormat)
        for fwd in itr:
            rev = next(itr)
            obs = _build_trim_command(self.demux_seqs, fwd[0], rev[0],
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
                                      match_adapter_wildcards=False)
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

            self.assertTrue(str(self.demux_seqs) in obs)

    def test_build_trim_command_multiple_adapters(self):
        itr = self.demux_seqs.sequences.iter_views(FastqGzFormat)
        for fwd in itr:
            rev = next(itr)
            obs = _build_trim_command(self.demux_seqs, fwd[0], rev[0],
                                      self.trimmed_seqs,
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


if __name__ == '__main__':
    unittest.main()
