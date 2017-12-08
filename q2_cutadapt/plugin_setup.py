# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib

from qiime2.plugin import Plugin, MetadataCategory, Visualization
from q2_types.multiplexed_sequences import (
    MultiplexedSingleEndBarcodeInSequence)
from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import SequencesWithQuality

import q2_cutadapt
import q2_cutadapt._demux
from ._format import CutadaptStatsFmt, CutadaptStatsDirFmt
from ._type import CutadaptStats


plugin = Plugin(
    name='cutadapt',
    version=q2_cutadapt.__version__,
    website='https://github.com/qiime2/q2-cutadapt',
    package='q2_cutadapt',
    description='This QIIME 2 plugin uses cutadapt to work with '
                'non-biological signal in sequences (e.g. barcodes, primers).',
    short_description='Plugin for handling non-biological signal in sequence '
                      'data.',
)

plugin.register_formats(CutadaptStatsFmt, CutadaptStatsDirFmt)
plugin.register_semantic_types(CutadaptStats)
plugin.register_semantic_type_to_format(
    CutadaptStats,
    artifact_format=CutadaptStatsDirFmt)
importlib.import_module('q2_cutadapt._transformer')

plugin.pipelines.register_function(
    function=q2_cutadapt._demux.demux_single,
    inputs={
        'seqs': MultiplexedSingleEndBarcodeInSequence,
    },
    parameters={
        'barcodes': MetadataCategory,
    },
    outputs=[
        ('per_sample_sequences', SampleData[SequencesWithQuality]),
        ('untrimmed_sequences', MultiplexedSingleEndBarcodeInSequence),
        ('demux_stats', Visualization)
    ],
    input_descriptions={
        'seqs': 'The single-end sequences to be demultiplexed.',
    },
    parameter_descriptions={
        'barcodes': 'The sample metadata category listing the per-sample '
                    'barcodes.',
    },
    output_descriptions={
        'per_sample_sequences': 'The resulting demultiplexed sequences.',
        'untrimmed_sequences': 'The sequences that were unmatched to '
                               'barcodes.',
        'demux_stats': 'The cutdapat-generated report for read-matching.',
    },
    name='Demultiplex single-end sequence data with barcodes in-sequence.',
    description='Demultiplex sequence data (i.e., map barcode reads to '
                'sample ids). Barcodes are expected to be located within the '
                'sequence data (versus the header, or a separate barcode '
                'file).',
)
