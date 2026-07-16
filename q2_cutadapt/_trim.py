# ----------------------------------------------------------------------------
# Copyright (c) 2017-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pandas as pd
from pathlib import Path
import tempfile
import warnings

import qiime2
from qiime2.plugin.util import run_commands
from qiime2.core.exceptions import RachisWarning

from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
)
from rachis.metadata import CategoricalMetadataColumn

from q2_cutadapt._stats import _summarize_cutadapt_json_reports


_trim_defaults = {
    'cores': 1,
    'adapter_f': None,
    'adapter_r': None,
    'front_f': None,
    'front_r': None,
    'anywhere_f': None,
    'anywhere_r': None,
    'forward_cut': 0,
    'reverse_cut': 0,
    'error_rate': 0.1,
    'indels': True,
    'times': 1,
    'overlap': 3,
    'match_read_wildcards': False,
    'match_adapter_wildcards': True,
    'minimum_length': 1,
    'discard_untrimmed': False,
    'discard_trimmed': False,
    'max_expected_errors': None,
    'max_n': None,
    'quality_cutoff_5end': 0,
    'quality_cutoff_3end': 0,
    'quality_base': 33,
    'pair_filter': None,
    'nextseq_trim': 0,
}


def _build_trim_command(
    f_read,
    r_read,
    trimmed_seqs,
    json_report_path,
    cores=_trim_defaults['cores'],
    adapter_f=_trim_defaults['adapter_f'],
    front_f=_trim_defaults['front_f'],
    anywhere_f=_trim_defaults['anywhere_f'],
    adapter_r=_trim_defaults['adapter_r'],
    front_r=_trim_defaults['front_r'],
    anywhere_r=_trim_defaults['anywhere_r'],
    forward_cut=_trim_defaults['forward_cut'],
    reverse_cut=_trim_defaults['reverse_cut'],
    error_rate=_trim_defaults['error_rate'],
    indels=_trim_defaults['indels'],
    times=_trim_defaults['times'],
    overlap=_trim_defaults['overlap'],
    match_read_wildcards=_trim_defaults['match_read_wildcards'],
    match_adapter_wildcards=_trim_defaults['match_adapter_wildcards'],
    minimum_length=_trim_defaults['minimum_length'],
    discard_untrimmed=_trim_defaults['discard_untrimmed'],
    discard_trimmed=_trim_defaults['discard_trimmed'],
    max_expected_errors=_trim_defaults['max_expected_errors'],
    max_n=_trim_defaults['max_n'],
    quality_cutoff_5end=_trim_defaults['quality_cutoff_5end'],
    quality_cutoff_3end=_trim_defaults['quality_cutoff_3end'],
    quality_base=_trim_defaults['quality_base'],
    nextseq_trim=_trim_defaults['nextseq_trim'],
    pair_filter=_trim_defaults['pair_filter'],
):
    if quality_cutoff_3end and nextseq_trim:
        warnings.warn(
            'Quality cutoff for 3 prime and NextSeq trimming is not supported '
            'ignoring quality cutoff for 3 prime.',
            RachisWarning
        )

    cmd = [
        'cutadapt',
        '-u', str(forward_cut),
        '--error-rate', str(error_rate),
        '--times', str(times),
        '--overlap', str(overlap),
        '--minimum-length', str(minimum_length)
    ]

    if nextseq_trim:
        cmd += [
            f'--nextseq-trim={nextseq_trim}',
            '-q', ','.join([str(quality_cutoff_5end), str(0)])
        ]
    else:
        cmd += [
            '-q', ','.join(
                [str(quality_cutoff_5end), str(quality_cutoff_3end)]
            )]

    cmd += [
        '--quality-base', str(quality_base),
        '--cores', str(cores),
        '-o', str(trimmed_seqs.path / os.path.basename(f_read)),
    ]

    if r_read is not None:
        cmd += ['-p', str(trimmed_seqs.path / os.path.basename(r_read))]

    if adapter_f is not None:
        for adapter in adapter_f:
            cmd += ['--adapter', adapter]
    if front_f is not None:
        for adapter in front_f:
            cmd += ['--front', adapter]
    if anywhere_f is not None:
        for adapter in anywhere_f:
            cmd += ['--anywhere', adapter]

    if adapter_r is not None:
        for adapter in adapter_r:
            cmd += ['-A', adapter]  # cutadapt doesn't have a long-form flag
    if front_r is not None:
        for adapter in front_r:
            cmd += ['-G', adapter]  # cutadapt doesn't have a long-form flag
    if anywhere_r is not None:
        for adapter in anywhere_r:
            cmd += ['-B', adapter]  # cutadapt doesn't have a long-form flag

    if reverse_cut is not None:
        cmd += ['-U', str(reverse_cut)]

    if not indels:
        cmd += ['--no-indels']
    if match_read_wildcards:
        cmd += ['--match-read-wildcards']
    if not match_adapter_wildcards:
        cmd += ['--no-match-adapter-wildcards']
    if discard_untrimmed:
        cmd += ['--discard-untrimmed']
    if discard_trimmed:
        cmd += ['--discard-trimmed']

    cmd += ['--json', str(json_report_path)]

    cmd += [f_read]

    if r_read is not None:
        cmd += [r_read]
    if max_expected_errors is not None:
        cmd += ['--max-expected-errors', str(max_expected_errors)]
    if max_n is not None:
        cmd += ['--max-n', str(max_n)]

    if pair_filter is not None:
        cmd += ['--pair-filter', pair_filter]

    return cmd


def _normalize_adapter(adapter):
    if isinstance(adapter, CategoricalMetadataColumn):
        adapter = adapter.to_dataframe().values.tolist()
        adapter = [adapt for sublist in adapter for adapt in sublist]

    return adapter


def trim_single(
    demultiplexed_sequences: SingleLanePerSampleSingleEndFastqDirFmt,
    adapter: str | CategoricalMetadataColumn = _trim_defaults['adapter_f'],
    front: str | CategoricalMetadataColumn = _trim_defaults['front_f'],
    anywhere: str | CategoricalMetadataColumn = _trim_defaults['anywhere_f'],
    cut: int = _trim_defaults['forward_cut'],
    error_rate: float = _trim_defaults['error_rate'],
    indels: bool = _trim_defaults['indels'],
    times: int = _trim_defaults['times'],
    overlap: int = _trim_defaults['overlap'],
    match_read_wildcards: bool = _trim_defaults['match_read_wildcards'],
    match_adapter_wildcards: bool = _trim_defaults['match_adapter_wildcards'],
    minimum_length: int = _trim_defaults['minimum_length'],
    discard_untrimmed: bool = _trim_defaults['discard_untrimmed'],
    discard_trimmed: bool = _trim_defaults['discard_trimmed'],
    max_expected_errors: float = _trim_defaults['max_expected_errors'],
    max_n: float = _trim_defaults['max_n'],
    quality_cutoff_5end: int = _trim_defaults['quality_cutoff_5end'],
    quality_cutoff_3end: int = _trim_defaults['quality_cutoff_3end'],
    quality_base: int = _trim_defaults['quality_base'],
    cores: int = _trim_defaults['cores'],
    nextseq_trim: int = _trim_defaults['nextseq_trim'],
) -> (CasavaOneEightSingleLanePerSampleDirFmt, qiime2.Metadata):

    adapter = _normalize_adapter(adapter)
    front = _normalize_adapter(front)
    anywhere = _normalize_adapter(anywhere)

    trimmed_sequences = CasavaOneEightSingleLanePerSampleDirFmt()
    cmds = []
    json_reports = {}
    df = demultiplexed_sequences.manifest.view(pd.DataFrame)
    with tempfile.TemporaryDirectory(prefix='q2-cutadapt-') as report_dir:
        for sample_id, fwd in df.itertuples():
            json_report_fp = Path(report_dir) / f'{len(json_reports)}.json'
            cmd = _build_trim_command(
                f_read=fwd,
                r_read=None,
                trimmed_seqs=trimmed_sequences,
                json_report_path=json_report_fp,
                adapter_f=adapter,
                front_f=front,
                anywhere_f=anywhere,
                adapter_r=None,
                front_r=None,
                anywhere_r=None,
                forward_cut=cut,
                reverse_cut=None,
                error_rate=error_rate,
                indels=indels,
                times=times,
                overlap=overlap,
                match_read_wildcards=match_read_wildcards,
                match_adapter_wildcards=match_adapter_wildcards,
                minimum_length=minimum_length,
                discard_untrimmed=discard_untrimmed,
                discard_trimmed=discard_trimmed,
                max_expected_errors=max_expected_errors,
                max_n=max_n,
                quality_cutoff_5end=quality_cutoff_5end,
                quality_cutoff_3end=quality_cutoff_3end,
                quality_base=quality_base,
                cores=cores,
                nextseq_trim=nextseq_trim,
            )
            cmds.append(cmd)
            json_reports[sample_id] = json_report_fp

        run_commands(cmds)
        stats = _summarize_cutadapt_json_reports(json_reports)

    return trimmed_sequences, stats


def trim_paired(
    demultiplexed_sequences: SingleLanePerSamplePairedEndFastqDirFmt,
    adapter_f: str | CategoricalMetadataColumn = _trim_defaults['adapter_f'],
    front_f: str | CategoricalMetadataColumn = _trim_defaults['front_f'],
    anywhere_f: str = _trim_defaults['anywhere_f'],
    adapter_r: str | CategoricalMetadataColumn = _trim_defaults['adapter_r'],
    front_r: str | CategoricalMetadataColumn = _trim_defaults['front_r'],
    anywhere_r: str = _trim_defaults['anywhere_r'],
    forward_cut: int = _trim_defaults['forward_cut'],
    reverse_cut: int = _trim_defaults['reverse_cut'],
    error_rate: float = _trim_defaults['error_rate'],
    indels: bool = _trim_defaults['indels'],
    times: int = _trim_defaults['times'],
    overlap: int = _trim_defaults['overlap'],
    match_read_wildcards: bool = _trim_defaults['match_read_wildcards'],
    match_adapter_wildcards: bool = _trim_defaults['match_adapter_wildcards'],
    minimum_length: int = _trim_defaults['minimum_length'],
    discard_untrimmed: bool = _trim_defaults['discard_untrimmed'],
    discard_trimmed: bool = _trim_defaults['discard_trimmed'],
    max_expected_errors: float = _trim_defaults['max_expected_errors'],
    max_n: float = _trim_defaults['max_n'],
    quality_cutoff_5end: int = _trim_defaults['quality_cutoff_5end'],
    quality_cutoff_3end: int = _trim_defaults['quality_cutoff_3end'],
    quality_base: int = _trim_defaults['quality_base'],
    cores: int = _trim_defaults['cores'],
    nextseq_trim: int = _trim_defaults['nextseq_trim'],
    pair_filter: str = 'any',
) -> (CasavaOneEightSingleLanePerSampleDirFmt, qiime2.Metadata):

    front_f = _normalize_adapter(front_f)
    front_r = _normalize_adapter(front_r)
    anywhere_f = _normalize_adapter(anywhere_f)
    anywhere_r = _normalize_adapter(anywhere_r)
    adapter_f = _normalize_adapter(adapter_f)
    adapter_r = _normalize_adapter(adapter_r)

    trimmed_sequences = CasavaOneEightSingleLanePerSampleDirFmt()
    cmds = []
    json_reports = {}
    df = demultiplexed_sequences.manifest.view(pd.DataFrame)
    with tempfile.TemporaryDirectory(prefix='q2-cutadapt-') as report_dir:
        for sample_id, fwd, rev in df.itertuples():
            json_report_fp = Path(report_dir) / f'{len(json_reports)}.json'
            cmd = _build_trim_command(
                f_read=fwd,
                r_read=rev,
                trimmed_seqs=trimmed_sequences,
                json_report_path=json_report_fp,
                adapter_f=adapter_f,
                front_f=front_f,
                anywhere_f=anywhere_f,
                adapter_r=adapter_r,
                front_r=front_r,
                anywhere_r=anywhere_r,
                forward_cut=forward_cut,
                reverse_cut=reverse_cut,
                error_rate=error_rate,
                indels=indels,
                times=times,
                overlap=overlap,
                match_read_wildcards=match_read_wildcards,
                match_adapter_wildcards=match_adapter_wildcards,
                minimum_length=minimum_length,
                discard_untrimmed=discard_untrimmed,
                discard_trimmed=discard_trimmed,
                max_expected_errors=max_expected_errors,
                max_n=max_n,
                quality_cutoff_5end=quality_cutoff_5end,
                quality_cutoff_3end=quality_cutoff_3end,
                quality_base=quality_base,
                cores=cores,
                nextseq_trim=nextseq_trim,
                pair_filter=pair_filter
            )
            cmds.append(cmd)
            json_reports[sample_id] = json_report_fp

        run_commands(cmds)
        stats = _summarize_cutadapt_json_reports(json_reports)

    return trimmed_sequences, stats
