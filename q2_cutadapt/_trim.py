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


def _build_trim_command(demux_seqs, fwd_read, cores, adapter, front, anywhere,
                        error_rate, no_indels, times, overlap,
                        match_read_wildcards, no_match_adapter_wildcards,
                        trimmed_sequences):
    cmd = [
        'cutadapt',
        '--cores', str(cores),
        '--error-rate', str(error_rate),
        '--times', str(times),
        '--overlap', str(overlap),
        '-o', str(trimmed_sequences.path / fwd_read),
    ]

    if adapter:
        for a in adapter:
            cmd += ['--adapter', a]
    if front:
        for f in front:
            cmd += ['--front', f]
    if anywhere:
        for a in anywhere:
            cmd += ['--anywhere', a]

    if no_indels:
        cmd += ['--no-indels']
    if match_read_wildcards:
        cmd += ['--match-read-wildcards']
    if no_match_adapter_wildcards:
        cmd += ['--no-match-adapter-wildcards']

    cmd += [str(demux_seqs.path / fwd_read)]

    return cmd


def trim_single(demultiplexed_sequences:
                SingleLanePerSampleSingleEndFastqDirFmt, cores: int=1,
                adapter: str=None, front: str=None, anywhere: str=None,
                error_rate: float=0.1, no_indels: bool=False, times: int=1,
                overlap: int=3, match_read_wildcards: bool=False,
                no_match_adapter_wildcards: bool=False) -> \
                    CasavaOneEightSingleLanePerSampleDirFmt:
    trimmed_sequences = CasavaOneEightSingleLanePerSampleDirFmt()
    cmds = []
    for fwd in demultiplexed_sequences.sequences.iter_views(FastqGzFormat):
        cmd = _build_trim_command(demultiplexed_sequences, fwd[0], cores,
                                  adapter, front, anywhere, error_rate,
                                  no_indels, times, overlap,
                                  match_read_wildcards,
                                  no_match_adapter_wildcards,
                                  trimmed_sequences)
        cmds.append(cmd)

    run_commands(cmds)
    return trimmed_sequences
