# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2.plugin.model as model


class CutadaptStatsFmt(model.TextFileFormat):
    def sniff(self):
        with self.open() as fh:
            count = 0
            while count < 10:
                line = fh.readline()

                if line == '':
                    # EOF
                    break
                else:
                    cells = line.split('\t')
                    if len(cells) < 4:
                        return False
                    count += 1
            return False if count == 0 else True


CutadaptStatsDirFmt = model.SingleFileDirectoryFormat(
    'CutadaptStatsDirFmt', 'stats.tsv', CutadaptStatsFmt)
