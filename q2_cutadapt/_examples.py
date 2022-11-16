# ----------------------------------------------------------------------------
# Copyright (c) 2022-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


# the examples are registered in plugin setup!
import os
import pkg_resources

import qiime2


#def _retrieve_data():
#    pass


def cutadapt_demux_single(use):
    seqs_fp = 
    sequence_url = 'https://data.qiime2.org/2022.8/tutorials/moving-pictures/emp-single-end-sequences/sequences.fastq.gz'
    barcode_url = 'https://data.qiime2.org/2022.8/tutorials/moving-pictures/emp-single-end-sequences/barcodes.fastq.gz'

    seqs = use.init_artifact_from_url(
            'seqs',
            sequence_url
            )

    barcodes = use.init_metadata_from_url(
            'barcodes',
            barcode_url
            )

    per_sample_sequences, untrimmed_sequences = use.action(

            use.UsageAction(plugin_id='cutadapt', action_id='demux_single'),
            use.UsageInputs(seqs=seqs, barcodes=barcodes), # barcodes_column='barcode-sequence'
            use.UsageOutputNames(per_sample_sequences='per_sample_sequences',
                                 untrimmed_sequences='untrimmed_sequences'))
    per_sample_sequences.assert_output_type('SampleData[SequencesWithQuality')
    untrimmed_sequences.assert_output_type('MultiplexedSingleEndBarcodeInSequence')


#   example of how to get the artifact from the url getter
#       def feature_table_filter_samples(use):
#           feature_table = use.init_artifact_from_url(
#               'feature_table',
#               'https://docs.qiime2.org/2022.8/data/tutorials/moving-pictures/table.qza'
#               )
#           sample_metadata = use.init_metadata_from_url(
#               'sample_metadata',
#               'https://data.qiime2.org/2022.8/tutorials/moving-pictures/sample_metadata.tsv'
#               )
#
#           filtered_table, = use.action(
#               use.UsageAction(plugin_id='feature_table', action_id='filter_samples'),
#               use.UsageInputs(table=feature_table, metadata=sample_metadata,
#                               where="[body-site] IN ('left palm', 'right palm')"),
#               use.UsageOutputNames(filtered_table='filtered_table')
#           )


# def _get_data_from_tests(path):
#    return pkg_resources.resource_filename('q2_cutadapt.tests',
#                                            os.path.join('data', path))

# -->
# def alpha_md_factory():
#    return qiime2.Metadata.load(
#        _get_data_from_tests('sample_metadata_alpha_div.tsv'))
#
#
# def beta_md_factory():
#    return qiime2.Metadata.load(
#        _get_data_from_tests('sample_metadata_donors.tsv'))
#
#
# def load_single_end_multiplexed_artifact_factory():
#     return qiime2.Artifact.load(_get_data_from_tests(os.path.join('single-end',
#                                                      'emp-single-end-sequences.qza')))

# def barcodes_artifact_factory():
#         return qiime2.Artifact.load(_get_data_from_tests(os.path.join(os.path.join('single-end', 'barcodes.fastq.gz'))))

# def beta_div_factory():
#    return qiime2.Artifact.import_data(
#        'DistanceMatrix', _get_data_from_tests('dist_matrix_donors.tsv'))
#
#
# def faithpd_md_factory():
#    return qiime2.Metadata.load(
#        _get_data_from_tests('metadata-faithpd.tsv')
#    )
#
#
# def faithpd_div_factory():
#    return qiime2.Artifact.import_data(
#        'SampleData[AlphaDiversity]', _get_data_from_tests('faithpd.tsv')
#    )
#
#
# def group_timepoints_alpha_independent(use):
#    alpha = use.init_artifact('alpha', load_single_end_multiplexed_data_factory)
#    metadata = use.init_metadata('metadata', alpha_md_factory)
#
#    timepoints, references = use.action(
#        use.UsageAction('fmt', 'group_timepoints'),
#        use.UsageInputs(
#            diversity_measure=alpha,
#            metadata=metadata,
#            time_column='days_post_transplant',
#            reference_column='relevant_donor',
#            subject_column=False
#        ),
#        use.UsageOutputNames(
#            timepoint_dists='timepoint_dists',
#            reference_dists='reference_dists'
#        )
#    )
#
#    timepoints.assert_output_type('GroupDist[Ordered, Independent]')
#    references.assert_output_type('GroupDist[Unordered, Independent]')
#
#
# def group_timepoints_beta(use):
#     beta = use.init_artifact('beta', beta_div_factory)
#     metadata = use.init_metadata('metadata', beta_md_factory)
#
#     timepoints, references = use.action(
#         use.UsageAction('fmt', 'group_timepoints'),
#         use.UsageInputs(
#             diversity_measure=beta,
#             metadata=metadata,
#             time_column='days_post_transplant',
#             reference_column='relevant_donor',
#             subject_column='subject',
#         ),
#         use.UsageOutputNames(
#             timepoint_dists='timepoint_dists',
#             reference_dists='reference_dists'
#         )
#     )
#
#     timepoints.assert_output_type('GroupDist[Ordered, Matched]')
#     references.assert_output_type('GroupDist[Unordered, Independent]')
#
#
# # Engraftment example using faith PD, baseline0 comparison
# def engraftment_baseline(use):
#     md = use.init_metadata('md', faithpd_md_factory)
#     div_measure = use.init_artifact('div_measure', faithpd_div_factory)
#
#     stats_table, raincloud = use.action(
#         use.UsageAction('fmt', 'engraftment'),
#         use.UsageInputs(
#             diversity_measure=div_measure,
#             metadata=md,
#             compare='baseline',
#             time_column='week',
#             reference_column='InitialDonorSampleID',
#             subject_column='SubjectID',
#             where='SampleType="stool"',
#             filter_missing_references=True,
#             against_group='0',
#             p_val_approx='asymptotic',
#         ),
#         use.UsageOutputNames(
#             stats='stats',
#             raincloud_plot='raincloud_plot'
#         )
#     )
#
#     stats_table.assert_output_type('StatsTable[Pairwise]')
#     raincloud.assert_output_type('Visualization')
