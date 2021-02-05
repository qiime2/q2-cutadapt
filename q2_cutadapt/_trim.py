# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os

import subprocess
import pandas as pd

from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
)


def run_commands(cmds, verbose=True):
    print('Running external command line application. This may print '
          'messages to stdout and/or stderr.')
    print('The commands to be run are below. These commands cannot '
          'be manually re-run as they will depend on temporary files that '
          'no longer exist.')
    for cmd in cmds:
        print('\nCommand:', end=' ')
        print(' '.join(cmd), end='\n\n')
        subprocess.run(cmd, check=True)


_trim_defaults = {
    'cores': 1,
    'adapter_f': None,
    'adapter_r': None,
    'front_f': None,
    'front_r': None,
    'anywhere_f': None,
    'anywhere_r': None,
    'error_rate': 0.1,
    'indels': True,
    'times': 1,
    'overlap': 3,
    'match_read_wildcards': False,
    'match_adapter_wildcards': True,
    'minimum_length': 1,
    'discard_untrimmed': False,
}


def _build_trim_command(f_read, r_read, trimmed_seqs,
                        cores=_trim_defaults['cores'],
                        adapter_f=_trim_defaults['adapter_f'],
                        front_f=_trim_defaults['front_f'],
                        anywhere_f=_trim_defaults['anywhere_f'],
                        adapter_r=_trim_defaults['adapter_r'],
                        front_r=_trim_defaults['front_r'],
                        anywhere_r=_trim_defaults['anywhere_r'],
                        error_rate=_trim_defaults['error_rate'],
                        indels=_trim_defaults['indels'],
                        times=_trim_defaults['times'],
                        overlap=_trim_defaults['overlap'],
                        match_read_wildcards=_trim_defaults[
                            'match_read_wildcards'],
                        match_adapter_wildcards=_trim_defaults[
                            'match_adapter_wildcards'],
                        minimum_length=_trim_defaults['minimum_length'],
                        discard_untrimmed=_trim_defaults['discard_untrimmed'],
                        ):
    cmd = [
        'cutadapt',
        '--cores', str(cores),
        '--error-rate', str(error_rate),
        '--times', str(times),
        '--overlap', str(overlap),
        '--minimum-length', str(minimum_length),
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

    if not indels:
        cmd += ['--no-indels']
    if match_read_wildcards:
        cmd += ['--match-read-wildcards']
    if not match_adapter_wildcards:
        cmd += ['--no-match-adapter-wildcards']
    if discard_untrimmed:
        cmd += ['--discard-untrimmed']

    cmd += [f_read]

    if r_read is not None:
        cmd += [r_read]

    return cmd


def trim_single(demultiplexed_sequences:
                SingleLanePerSampleSingleEndFastqDirFmt,
                cores: int = _trim_defaults['cores'],
                adapter: str = _trim_defaults['adapter_f'],
                front: str = _trim_defaults['front_f'],
                anywhere: str = _trim_defaults['anywhere_f'],
                error_rate: float = _trim_defaults['error_rate'],
                indels: bool = _trim_defaults['indels'],
                times: int = _trim_defaults['times'],
                overlap: int = _trim_defaults['overlap'],
                match_read_wildcards:
                bool = _trim_defaults['match_read_wildcards'],
                match_adapter_wildcards:
                bool = _trim_defaults['match_adapter_wildcards'],
                minimum_length: int = _trim_defaults['minimum_length'],
                discard_untrimmed:
                bool = _trim_defaults['discard_untrimmed']) -> \
                    CasavaOneEightSingleLanePerSampleDirFmt:
    trimmed_sequences = CasavaOneEightSingleLanePerSampleDirFmt()
    cmds = []
    df = demultiplexed_sequences.manifest.view(pd.DataFrame)
    for _, fwd in df.itertuples():
        cmd = _build_trim_command(fwd, None,
                                  trimmed_sequences, cores, adapter, front,
                                  anywhere, None, None, None, error_rate,
                                  indels, times, overlap, match_read_wildcards,
                                  match_adapter_wildcards, minimum_length,
                                  discard_untrimmed)
        cmds.append(cmd)

    run_commands(cmds)
    return trimmed_sequences


def trim_paired(demultiplexed_sequences:
                SingleLanePerSamplePairedEndFastqDirFmt,
                cores: int = _trim_defaults['cores'],
                adapter_f: str = _trim_defaults['adapter_f'],
                front_f: str = _trim_defaults['front_f'],
                anywhere_f: str = _trim_defaults['anywhere_f'],
                adapter_r: str = _trim_defaults['adapter_r'],
                front_r: str = _trim_defaults['front_r'],
                anywhere_r: str = _trim_defaults['anywhere_r'],
                error_rate: float = _trim_defaults['error_rate'],
                indels: bool = _trim_defaults['indels'],
                times: int = _trim_defaults['times'],
                overlap: int = _trim_defaults['overlap'],
                match_read_wildcards:
                bool = _trim_defaults['match_read_wildcards'],
                match_adapter_wildcards:
                bool = _trim_defaults['match_adapter_wildcards'],
                minimum_length: int = _trim_defaults['minimum_length'],
                discard_untrimmed:
                bool = _trim_defaults['discard_untrimmed']) -> \
                    CasavaOneEightSingleLanePerSampleDirFmt:
    trimmed_sequences = CasavaOneEightSingleLanePerSampleDirFmt()
    cmds = []
    df = demultiplexed_sequences.manifest.view(pd.DataFrame)
    for _, fwd, rev in df.itertuples():
        cmd = _build_trim_command(fwd, rev, trimmed_sequences, cores,
                                  adapter_f, front_f,
                                  anywhere_f, adapter_r, front_r, anywhere_r,
                                  error_rate, indels, times, overlap,
                                  match_read_wildcards,
                                  match_adapter_wildcards, minimum_length,
                                  discard_untrimmed)
        cmds.append(cmd)

    run_commands(cmds)
    return trimmed_sequences
