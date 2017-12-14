# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import Plugin, MetadataCategory, Float, Range
from q2_types.multiplexed_sequences import (
    MultiplexedSingleEndBarcodeInSequence)
from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import SequencesWithQuality

import q2_cutadapt
import q2_cutadapt._demux


plugin = Plugin(
    name='cutadapt',
    version=q2_cutadapt.__version__,
    website='https://github.com/qiime2/q2-cutadapt',
    package='q2_cutadapt',
    description='This QIIME 2 plugin uses cutadapt to work with '
                'adapters (e.g. barcodes, primers) in sequence data.',
    short_description='Plugin for removing adapter sequences, primers, and '
                      'other unwanted sequence from sequence data.',
    citation_text='Martin, M. (2011). Cutadapt removes adapter sequences from '
                  'high-throughput sequencing reads. EMBnet.Journal, 17(1), '
                  'pp. 10-12.\ndoi:http://dx.doi.org/10.14806/ej.17.1.200',
)

plugin.methods.register_function(
    function=q2_cutadapt._demux.demux_single,
    inputs={
        'seqs': MultiplexedSingleEndBarcodeInSequence,
    },
    parameters={
        'barcodes': MetadataCategory,
        'error_tolerance': Float % Range(0, 1, inclusive_start=True,
                                         inclusive_end=True),
    },
    outputs=[
        ('per_sample_sequences', SampleData[SequencesWithQuality]),
        ('untrimmed_sequences', MultiplexedSingleEndBarcodeInSequence),
    ],
    input_descriptions={
        'seqs': 'The single-end sequences to be demultiplexed.',
    },
    parameter_descriptions={
        'barcodes': 'The sample metadata category listing the per-sample '
                    'barcodes.',
        'error_tolerance': 'The level of error tolerance, specified as the '
                           'maximum allowable error rate.',
    },
    output_descriptions={
        'per_sample_sequences': 'The resulting demultiplexed sequences.',
        'untrimmed_sequences': 'The sequences that were unmatched to '
                               'barcodes.',
    },
    name='Demultiplex single-end sequence data with barcodes in-sequence.',
    description='Demultiplex sequence data (i.e., map barcode reads to '
                'sample ids). Barcodes are expected to be located within the '
                'sequence data (versus the header, or a separate barcode '
                'file).',
)
