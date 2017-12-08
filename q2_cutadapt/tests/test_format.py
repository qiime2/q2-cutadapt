# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from q2_cutadapt import CutadaptStatsFmt
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugin import ValidationError


class TestFormats(TestPluginBase):
    package = 'q2_cutadapt.tests'

    def test_cutadapts_stats_format_validate_positive(self):
        filepath = self.get_data_path('stats.tsv')
        format_ = CutadaptStatsFmt(filepath, mode='r')

        format_.validate()

    def test_cutadapts_stats_format_validate_negative(self):
        filepath = self.get_data_path('stats-bad.tsv')
        format_ = CutadaptStatsFmt(filepath, mode='r')

        with self.assertRaisesRegex(ValidationError, 'CutadaptStatsFmt'):
            format_.validate()


if __name__ == '__main__':
    unittest.main()
