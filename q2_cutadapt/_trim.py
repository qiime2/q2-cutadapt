# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import subprocess

from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    FastqGzFormat,
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


def _build_trim_command(demux_seqs, fwd_read, rev_read, cores, adapter_fwd,
                        front_fwd, anywhere_fwd, adapter_rev, front_rev,
                        anywhere_rev, error_rate, indels, times, overlap,
                        match_read_wildcards, match_adapter_wildcards,
                        trimmed_sequences):
    cmd = [
        'cutadapt',
        '--cores', str(cores),
        '--error-rate', str(error_rate),
        '--times', str(times),
        '--overlap', str(overlap),
        '-o', str(trimmed_sequences.path / fwd_read),
    ]

    if rev_read is not None:
        cmd += ['-p', str(trimmed_sequences.path / rev_read)]

    if adapter_fwd:
        for a in adapter_fwd:
            cmd += ['--adapter', a]
    if front_fwd:
        for f in front_fwd:
            cmd += ['--front', f]
    if anywhere_fwd:
        for a in anywhere_fwd:
            cmd += ['--anywhere', a]

    if adapter_rev:
        for a in adapter_rev:
            cmd += ['-A', a]  # cutadapt doesn't have a long-form flag
    if front_rev:
        for f in front_rev:
            cmd += ['-G', f]  # cutadapt doesn't have a long-form flag
    if anywhere_rev:
        for a in anywhere_rev:
            cmd += ['-B', a]  # cutadapt doesn't have a long-form flag

    if not indels:
        cmd += ['--no-indels']
    if match_read_wildcards:
        cmd += ['--match-read-wildcards']
    if not match_adapter_wildcards:
        cmd += ['--no-match-adapter-wildcards']

    cmd += [str(demux_seqs.path / fwd_read)]

    if rev_read is not None:
        cmd += [str(demux_seqs.path / rev_read)]

    return cmd


def trim_single(demultiplexed_sequences:
                SingleLanePerSampleSingleEndFastqDirFmt, cores: int=1,
                adapter: str=None, front: str=None, anywhere: str=None,
                error_rate: float=0.1, indels: bool=True, times: int=1,
                overlap: int=3, match_read_wildcards: bool=False,
                match_adapter_wildcards: bool=True) -> \
                    CasavaOneEightSingleLanePerSampleDirFmt:
    trimmed_sequences = CasavaOneEightSingleLanePerSampleDirFmt()
    cmds = []
    for fwd in demultiplexed_sequences.sequences.iter_views(FastqGzFormat):
        cmd = _build_trim_command(demultiplexed_sequences, fwd[0], None, cores,
                                  adapter, front, anywhere, None, None, None,
                                  error_rate, indels, times, overlap,
                                  match_read_wildcards,
                                  match_adapter_wildcards, trimmed_sequences)
        cmds.append(cmd)

    run_commands(cmds)
    return trimmed_sequences


def trim_paired(demultiplexed_sequences:
                SingleLanePerSamplePairedEndFastqDirFmt, cores: int=1,
                adapter_f: str=None, front_f: str=None, anywhere_f: str=None,
                adapter_r: str=None, front_r: str=None, anywhere_r: str=None,
                error_rate: float=0.1, indels: bool=True, times: int=1,
                overlap: int=3, match_read_wildcards: bool=False,
                match_adapter_wildcards: bool=True) -> \
                    CasavaOneEightSingleLanePerSampleDirFmt:
    trimmed_sequences = CasavaOneEightSingleLanePerSampleDirFmt()
    cmds = []
    seqs = demultiplexed_sequences.sequences.iter_views(FastqGzFormat)
    for fwd in seqs:
        rev = next(seqs)
        cmd = _build_trim_command(demultiplexed_sequences, fwd[0], rev[0],
                                  cores, adapter_f, front_f, anywhere_f,
                                  adapter_r, front_r, anywhere_r, error_rate,
                                  indels, times, overlap, match_read_wildcards,
                                  match_adapter_wildcards, trimmed_sequences)
        cmds.append(cmd)

    run_commands(cmds)
    return trimmed_sequences
