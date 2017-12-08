# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import qiime2

from .plugin_setup import plugin
from ._format import CutadaptStatsFmt


STATS_COLS = [
    'Read name',
    'Number of errors',
    'Start coordinate for adapter match or Read sequence',
    'End coordinate for adapter match or Read quality',
    'Sequence to left of adapter match',
    'Sequence of adapter match',
    'Sequence to right of adapter match',
    'Adapter name (sample ID)',
    'Quality scores to left of adapter match',
    'Quality scores of adapter match',
    'Quality scores to right of adapter match',
]


@plugin.register_transformer
def _0(ff: CutadaptStatsFmt) -> qiime2.Metadata:
    df = pd.read_csv(str(ff), index_col=0, header=None, names=STATS_COLS,
                     delimiter='\t')
    return qiime2.Metadata(df)
