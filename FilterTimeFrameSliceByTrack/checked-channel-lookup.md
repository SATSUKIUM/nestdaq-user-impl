# Initialization-only readiness checks (2026-10-30)

Requires channel-map d0e88a8f971ad713e6698c90b76efe6697973c27 or later, from main.
Rebuild both the library and this consumer against matching headers.

InitTask validates isInitialized(), hasInverseMap(), and the KLDC dictionary
entry once. Lookup functions no longer check readiness per call. Do not call
findDET before successful initialization, or findFE without an inverse map.
Detector configuration registration retains its readiness/input checks.
ProcessSlice no longer repeats the KLDC presence check. Hit address range and
hole checks remain. Tracking is still disabled by the existing #if 0; the
active ProcessSlice returns false. This change does not enable tracking.

Library regression tests and optimized assembly inspection passed. Full
consumer compilation/DAQ performance measurement remain unavailable due to
missing dependencies in this environment.

Earlier migration history follows.

---

# Map-owned configuration migration (2026-10-30)

Requires channel-map `7b2762346b4d3d41679e17f96e7af6cce7547a1a` (branch `feat/map-owned-config`).
Rebuild/install channel-map, then rebuild this consumer with matching headers
and library. Previous raw-pointer/template registration APIs have been removed.

This branch combines a8e6946 (checked lookups) and origin/LogicFilterMore32GPU
at e14424d (CUDA visibility), which had one unique commit each. main at 22e01d0
and LogicFilterMore32 at a4ef1c8 were already ancestors of the checked-lookup
branch. Neither existing branch alone contained both updates. No existing branch
was rebased or force-pushed; the latest GPU change was merged normally.

Configuration registration now calls setGeometry, setTdcCalibration and
setDriftCalibration. The map owns every accepted object until its destruction;
borrowed const pointers are returned by getConfiguration and its typed getters.
Duplicate assignments return false without replacing existing objects. Failed
arguments are consumed/destroyed. Tracking registration checks failures and
counts only successful calibration registrations. DCHit uses the const accessors.

Validation: channel-map build/install, skeleton, ownership and regression tests
passed. Consumer CMake configuration remains blocked by unavailable FairMQ,
hiredis, redis++ and UHBOOK dependencies. Full tracking/CUDA build and runtime
validation must be performed in the deployment environment. Do not mix this
source with an older installed channel-map.

The record below describes the previous lookup-only migration, not the new API.

---

# Checked channel lookup migration

Base: LogicFilterMore32GPU, 17bc9438967feda5ff0e659021d6576c0e53d582.
Requires channel-map commit 981bd850c9e8723a77c74e896c2461f9a9146382
(branch fix/checked-channel-lookup).

Rebuild and install that channel-map version, then rebuild this consumer against
its headers and library. An older installed ChannelMap does not provide findDET
or findFE. Rebuild both together; detector name indices have shifted by one.

FilterTimeFrameSliceByTrack now uses findDET/findFE for the hit loop,
diagnostics, and detector geometry/calibration registration. Missing entries
return nullptr, including holes inside allocated ranges. KLDC's dictionary ID
is resolved once at initialization instead of hard-coding 0x06. Hit channels
above 255 are rejected before conversion to the map's uint8_t input. Raw hits
retain const pointers to map items. Drift configuration lookup passes integer
segment/channel values without premature narrowing.

Detector name ID 0 and the all-zero FE address tuple are reserved for unmapped
slots. CSV input containing the latter is rejected by channel-map. Individual
zero fields, including CH0, remain valid.

Validation on macOS (2026-10-25): channel-map's three functional cases passed.
The consumer's CMake configuration could not complete because FairMQ, hiredis,
redis++ and UHBOOK were unavailable. The installed system ChannelMap was also
an older version. Full FilterTimeFrameSliceByTrack compilation and runtime
validation must therefore be performed in the deployment environment.
Direct syntax checks of tracking helper files were also blocked by missing
fairmq/Device.h and pre-existing acos-based constexpr definitions rejected by
AppleClang. Those unrelated portability issues were left unchanged.
