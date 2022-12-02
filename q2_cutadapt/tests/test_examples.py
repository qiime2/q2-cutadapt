# ----------------------------------------------------------------------------
# Copyright (c) 2017-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import unittest
from qiime2.plugin.testing import TestPluginBase


class TestExamples(TestPluginBase):
    package = 'q2_cutadapt.tests'

    def test_examples(self):
        self.execute_examples()


if __name__ == '__main__':
    unittest.main()
