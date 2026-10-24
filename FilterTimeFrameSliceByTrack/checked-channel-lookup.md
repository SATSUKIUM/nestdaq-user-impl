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
