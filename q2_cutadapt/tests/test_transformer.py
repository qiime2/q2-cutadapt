# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import pandas as pd
import numpy as np

from q2_cutadapt import CutadaptStatsFmt
from q2_cutadapt._transformer import STATS_COLS
from qiime2 import Metadata
from qiime2.plugin.testing import TestPluginBase


class TestTransformers(TestPluginBase):
    package = 'q2_cutadapt.tests'

    def test_cutadapt_stats_fmt_to_metadata(self):
        _, obs = self.transform_format(CutadaptStatsFmt, Metadata, 'stats.tsv')

        idx = pd.Index(['id1', 'id2', 'id3'], name=STATS_COLS[0])

        exp = Metadata(pd.DataFrame([[0, '0', '4', np.NaN, 'AAAA', 'ACGTACGT',
                                      'sample_a', np.NaN, 'zzzz', 'zzzzzzzz'],
                                     [0, '0', '4', np.NaN, 'CCCC', 'ACGTACGT',
                                      'sample_b', np.NaN, 'zzzz', 'zzzzzzzz'],
                                     [-1, 'GGGGACGTACGT', 'zzzzzzzzzzzz',
                                      np.NaN, np.NaN, np.NaN, np.NaN, np.NaN,
                                      np.NaN, np.NaN]],
                                    index=idx, columns=STATS_COLS[1:]))
        self.assertEqual(exp, obs)


if __name__ == '__main__':
    unittest.main()
