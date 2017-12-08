# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from qiime2.plugin.testing import TestPluginBase


class TestDemuxSingle(TestPluginBase):
    package = 'q2_cutadapt.tests'

    def test_nothing(self):
        # This is just so that the tests won't error when run
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
