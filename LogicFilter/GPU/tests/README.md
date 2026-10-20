# GPU regression test

```sh
nvcc -std=c++17 -arch=sm_86 -rdc=true -I. -ILogicFilter/GPU -I$NESTDAQ/include \
  LogicFilter/GPU/tests/TriggerGpuTest.cu "$NESTDAQ/lib/libcudamap.a" \
  -o /tmp/trigger-gpu-test
/tmp/trigger-gpu-test
```

Use the CUDA architecture supported by the installed cudamap library. Exit 77
means that no CUDA device is available. The test checks 1/2/32-bit subgroup
initialization, subgroup-to-main bit mapping, concurrent overlapping marks,
four TDC formats, and reset between frames.

The mapping passed to SetfiSubTimeRegionMainTimeRegion contains bit masks
(`1u << group_id`), not bit indices. Call it after SetSubTimeRegion.
Each subgroup supports 1..32 signals; decoder subgroup indices are 0..254.
The decoder stores one configuration per FEM/channel; duplicate registrations
are rejected. Window width must be 1..255; marking uses inclusive +/- width/2
windows as in the CPU implementation. The LUT is copied to device byte storage. Each bin detects a rising edge,
CUB DeviceScan::ExclusiveSum assigns output offsets, and a compaction kernel
writes timestamps in ascending order. Bin 0 is excluded, matching CPU behavior.
The device output has floor(N/2) slots. Only a uint32_t count and exactly that
many timestamps are copied to the host; GetTimeRegion explicitly downloads the
full array for diagnostics. Scratch buffers are allocated when the region is
sized and reused across frames. The compaction test compares against CPU LUT
lookup for empty/full/alternating/random patterns, out-of-range LUT inputs,
negation, odd/even sizes, block boundaries, resizing and repeated scans. Mark reuses its device input buffer but currently copies and synchronizes
once per HBF segment; end-to-end speedup requires measurements on a CUDA host.
