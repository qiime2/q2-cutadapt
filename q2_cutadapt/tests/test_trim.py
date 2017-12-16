# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import unittest

from q2_cutadapt._trim import _build_trim_command
from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
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
        demuxed = Artifact.import_data('SampleData[SequencesWithQuality]',
                                       self.get_data_path('single-end'))
        with redirected_stdio(stdout=os.devnull):
            self.plugin.methods['trim_single'](demuxed)


class TestTrimUtilsSingle(TestPluginBase):
    package = 'q2_cutadapt.tests'

    def setUp(self):
        super().setUp()

        self.demux_seqs = SingleLanePerSampleSingleEndFastqDirFmt(
            self.get_data_path('single-end'), mode='r')
        self.trimmed_sequences = CasavaOneEightSingleLanePerSampleDirFmt()

    def test_build_trim_command_typical(self):
        for fwd in self.demux_seqs.sequences.iter_views(FastqGzFormat):
            obs = _build_trim_command(self.demux_seqs, fwd[0], 1, ['AAAA'],
                                      ['GGGG'], ['CCCC'], 2, True, 3,
                                      4, True, True, self.trimmed_sequences)
            obs = ' '.join(obs)

            self.assertTrue('-o %s' % str(self.trimmed_sequences.path / fwd[0])
                            in obs)
            self.assertTrue('--cores 1' in obs)
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
            obs = _build_trim_command(self.demux_seqs, fwd[0], 1,
                                      ['AAAA', 'GGGG', 'CCCC'], [], [], 2,
                                      True, 3, 4, True, True,
                                      self.trimmed_sequences)
            obs = ' '.join(obs)

            self.assertTrue('--adapter AAAA' in obs)
            self.assertTrue('--adapter GGGG' in obs)
            self.assertTrue('--adapter CCCC' in obs)
            self.assertTrue('--front' not in obs)
            self.assertTrue('--anywhere' not in obs)

    def test_build_trim_command_no_bool_flags(self):
        for fwd in self.demux_seqs.sequences.iter_views(FastqGzFormat):
            obs = _build_trim_command(self.demux_seqs, fwd[0], 1,
                                      ['AAAA', 'GGGG', 'CCCC'], [], [], 2,
                                      False, 3, 4, False, False,
                                      self.trimmed_sequences)
            obs = ' '.join(obs)

            self.assertTrue('--no-indels' not in obs)
            self.assertTrue('--match-read-wildcards' not in obs)
            self.assertTrue('--no-match-adapter-wildcards' not in obs)


if __name__ == '__main__':
    unittest.main()
