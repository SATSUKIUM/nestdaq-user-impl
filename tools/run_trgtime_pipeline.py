#!/usr/bin/env python3
"""Replay the Tracking topology over isolated IPC sockets, without Redis/tmux.

Defaults reproduce the trigger settings in nestdaq-run/LogicFilter/triggerlogic.sh.
Each invocation records exact commands in its output directory.
"""
import argparse
import json
from pathlib import Path
import re
import signal
import struct
import subprocess
import tempfile
import time

from compare_trgtime import records

SIGNALS = '(0 0xc0a802a9 8 0) (1 0xc0a802a9 10 0) ' + ' '.join(
    f'(2-{i} 0xc0a802aa {a} -12) (2-{i} 0xc0a802aa {b} -12)'
    for i, (a, b) in enumerate([(16, 17), (18, 19), (20, 21), (22, 23), (24, 25), (27, 28)]))


def count_input_frames(path):
    """Count v1 TFs by validated record boundaries without reading their payloads."""
    count = 0
    size = path.stat().st_size
    with path.open('rb') as stream:
        while stream.tell() < size:
            offset = stream.tell()
            header = stream.read(16)
            if len(header) != 16:
                raise ValueError(f'Truncated record header at {offset}')
            magic, length, hlength, _ = struct.unpack('<QIHH', header)
            if magic == 0x004d5246454d4954:  # TimeFrame::v1
                minimum = 24
                count += 1
            elif magic in (0x004b4e53454c4946, 0x004c5254454c4946):
                minimum = 304  # FileSink header/trailer
                if length != minimum or hlength != minimum:
                    raise ValueError(f'Unsupported FileSink record at {offset}')
            elif magic == 0x0049474f4c544c46:  # TFBFilePlayer skips the Filter header only
                minimum = length = 56
            else:
                raise ValueError(f'Unsupported record magic {magic:#x} at {offset}')
            if not minimum <= hlength <= length or offset + length > size:
                raise ValueError(f'Invalid/truncated record at {offset}')
            stream.seek(offset + length)
    if not count:
        raise ValueError('Input contains no TFs')
    return count


class TraceProgress:
    """Read complete CSV lines incrementally while the receiver is still running."""
    def __init__(self):
        self.pending = ''
        self.messages = 0
        self.frames = self.expected_frames = 0
        self.metas = self.expected_metas = 0
        self.triggers = self.expected_triggers = 0

    def feed(self, chunk):
        lines = (self.pending + chunk).split('\n')
        self.pending = lines.pop()
        for line in lines:
            row = line.split(',')
            if row[0] == 'MESSAGE':
                if self.messages and not self.message_complete():
                    raise RuntimeError('Previous trace message is incomplete')
                if int(row[1]) != self.messages:
                    raise RuntimeError('Nonconsecutive trace message')
                self.messages += 1
                self.expected_frames, self.expected_metas = map(int, row[2:4])
                self.frames = self.metas = self.triggers = self.expected_triggers = 0
            elif row[0] == 'META':
                self.frames += 1
                self.metas += 1
                self.expected_triggers += int(row[5])
            elif row[0] == 'SLICE':
                self.frames += 1
            elif row[0] == 'TRG':
                self.triggers += 1

    def message_complete(self):
        return (self.frames == self.expected_frames and self.metas == self.expected_metas
                and self.triggers == self.expected_triggers)

    def complete(self, expected):
        return self.messages == expected and self.message_complete() and not self.pending


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build-dir', required=True, type=Path)
    parser.add_argument('--run-dir', required=True, type=Path)
    parser.add_argument('--input', type=Path, help='Default: RUN_DIR/data/run000400.dat')
    parser.add_argument('--output-dir', type=Path, help='Must not exist; default: new /tmp directory')
    parser.add_argument('--iterations', type=int, help='Optional TF limit; default: replay the entire input file')
    parser.add_argument('--wait', type=int, default=1, help='Player pacing in ms (default: 1; use 0 for unpaced replay)')
    parser.add_argument('--timeout', type=float, default=180)
    parser.add_argument('--gpu-exec-cuda-events', action='store_true', help='Also insert CUDA events (implies --gpu-exec-details, adds overhead)')
    parser.add_argument('--gpu-exec-details', action='store_true', help='Profile host stages within GPU Exec')
    parser.add_argument('--gpu-hbf-batch-size', type=int, choices=range(1, 5), help='GPU HBF batch size (default: 4; use 1 for sequential reference)')
    parser.add_argument('--performance', action='store_true', help='Save LogicFilter timing/data TTree as logicfilter.root')
    args = parser.parse_args()
    if args.gpu_exec_cuda_events:
        args.gpu_exec_details = True
    if (args.iterations is not None and args.iterations < 1) or args.wait < 0 or args.timeout <= 0:
        parser.error('Require iterations>0, wait>=0, timeout>0')
    build = args.build_dir.resolve()
    # Consult the cache rather than accidentally selecting an old binary.
    cache = (build / 'CMakeCache.txt').read_text()
    mode = re.search(r'^CMAKE_LOGICFILTER_CUDA_ENABLE:[^=]+=(.*)$', cache, re.M)
    if not mode:
        parser.error('Build cache has no CUDA option')
    variant = 'GPU' if mode.group(1).upper() in {'1', 'ON', 'TRUE', 'YES'} else 'CPU'
    if args.gpu_exec_details and (not args.performance or variant != 'GPU'):
        parser.error('--gpu-exec-details requires --performance and a GPU build')
    if args.gpu_hbf_batch_size is not None and variant != 'GPU':
        parser.error('--gpu-hbf-batch-size requires a GPU build')
    logic = build / f'LogicFilter/{variant}/LogicFilter'
    track = build / 'FilterTimeFrameSliceByTrack/FilterTimeFrameSliceByTrack'
    for binary in [logic, track, *(build / name for name in ['FileSink', 'TFBFilePlayer', 'TimeFrameSlicerByLogicTiming'])]:
        if not binary.is_file():
            parser.error(f'Missing executable: {binary}')
    run = args.run_dir.resolve()
    input_file = (args.input or run / 'data/run000400.dat').resolve()
    if not input_file.is_file():
        parser.error(f'Missing input: {input_file}')
    expected = args.iterations if args.iterations is not None else count_input_frames(input_file)
    print(f'INPUT_TFS={expected}; replay={"limited" if args.iterations else "entire file"}', flush=True)
    if args.output_dir:
        work = args.output_dir.resolve()
        work.mkdir(parents=True, exist_ok=False)
    else:
        work = Path(tempfile.mkdtemp(prefix='ndq-trgtime-'))
    (work / 'fileout/tracking').mkdir(parents=True)
    print(f'RESULT_DIR={work}', flush=True)
    processes = []
    files = []

    def channel(name, kind, method, index):
        return f'name={name},type={kind},method={method},address=ipc://{work}/q{index}'

    def launch(name, binary, options, channels):
        log = (work / f'{name}.log').open('w')
        files.append(log)
        command = [str(binary), '--id', name + '-0', '--control', 'static',
                   '--severity', 'info', '--transport', 'zeromq', *options,
                   '--channel-config', *channels]
        (work / f'{name}.command.json').write_text(json.dumps(command, indent=2))
        process = subprocess.Popen(command, cwd=work, stdout=log, stderr=subprocess.STDOUT)
        processes.append((name, process))
        return process

    def check_receivers():
        for name, process in processes:
            if name == 'FileSink' and process.poll() == 0:
                continue  # The finite sink exits only after writing all requested messages.
            if name != 'TFBFilePlayer' and process.poll() is not None:
                raise RuntimeError(f'{name} exited early: {process.returncode}')

    try:
        sink = launch('FileSink', build / 'FileSink', ['--path', str(work / 'sink.dat'), '--n-threads', '1',
                      '--max-iteration', str(expected)],
               [channel('in', 'pull', 'bind', 3)])
        launch('FilterTimeFrameSliceByTrack', track, [
            '--chmap-data-file', str(run / 'mapdata_T103_20260909.csv'),
            '--geometry-file', str(run / 'detector_configurations/param_DetGeom001'),
            '--dctdc-calib-file', str(run / 'detector_configurations/param_DCTdcCalib_run0602'),
            '--dc-drift-param-file', str(run / 'detector_configurations/param_DCDrift_run0602'),
            '--trgtime-output', str(work / 'trgtime.csv')],
            [channel('in', 'pull', 'bind', 2), channel('out', 'push', 'connect', 3)])
        launch('TimeFrameSlicerByLogicTiming', build / 'TimeFrameSlicerByLogicTiming',
               ['--time-offset-begin', '-250', '--time-offset-end', '250'],
               [channel('in', 'pull', 'bind', 1), channel('out', 'push', 'connect', 2)])
        performance_options = ['--performance-output', str(work / 'logicfilter.root')] if args.performance else []
        if variant == 'GPU':
            performance_options += ['--hbf-batch-size', str(args.gpu_hbf_batch_size or 4)]
        if args.gpu_exec_details:
            performance_options += ['--exec-details', 'true']
        if args.gpu_exec_cuda_events:
            performance_options += ['--exec-cuda-events', 'true']
        launch('LogicFilter', logic, [*performance_options, '--trigger-signals', SIGNALS,
               '--trigger-expression', '(0 & 1) & 2', '--trigger-width', '30'],
               [channel('in', 'pull', 'bind', 0), channel('out', 'push', 'connect', 1)])
        time.sleep(1)
        check_receivers()
        ack = work / 'downstream-complete.ack'
        player_limit = ['--max-iterations', str(args.iterations)] if args.iterations is not None else []
        player = launch('TFBFilePlayer', build / 'TFBFilePlayer', ['--in-file', str(input_file),
                        *player_limit, '--wait', str(args.wait),
                        '--completion-ack-file', str(ack)],
                        [channel('out', 'push', 'connect', 0)])
        deadline = time.monotonic() + args.timeout
        progress = TraceProgress()
        trace_stream = None
        try:
            while time.monotonic() < deadline:
                check_receivers()
                trace = work / 'trgtime.csv'
                if trace_stream is None and trace.exists():
                    trace_stream = trace.open()
                if trace_stream is not None:
                    progress.feed(trace_stream.read())
                if progress.complete(expected) and not ack.exists():
                    ack.write_text('All requested TFs processed by ByTrack\n')
                if player.poll() is not None:
                    if player.returncode or not ack.exists():
                        raise RuntimeError(f'TFBFilePlayer exited before acknowledged completion: {player.returncode}')
                    if sink.poll() == 0:
                        break
                time.sleep(.05)
            else:
                raise RuntimeError(f'Pipeline timeout: received {progress.messages}/{expected} messages; '
                                   f'last message complete={progress.message_complete()}')
        finally:
            if trace_stream is not None:
                trace_stream.close()
        print(f'{variant}: all {expected} messages processed and written by FileSink', flush=True)
    finally:
        for name, process in reversed(processes):
            if process.poll() is None:
                process.send_signal(signal.SIGINT)
                try:
                    process.wait(timeout=8)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
            print(f'{name} exit={process.returncode}', flush=True)
        for stream in files:
            stream.close()

    # END alone does not establish that upstream queues were drained.
    for record in records(work / 'trgtime.csv'):
        last = record
    if last[1] != expected:
        raise RuntimeError(f'Received {last[1]} messages, expected {expected}; check logs/--timeout')
    for name, process in processes:
        if process.returncode != 0:
            raise RuntimeError(f'{name} failed: {process.returncode}')
    print(f'COMPLETE {variant}: messages={last[1]}, metas={last[2]}, TrgTime={last[3]}, slices={last[4]}')


if __name__ == '__main__':
    main()
