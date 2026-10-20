#ifndef NESTDAQ_LOGICFILTER_EXEC_PROFILE_DATA_H
#define NESTDAQ_LOGICFILTER_EXEC_PROFILE_DATA_H
#include <cstdint>

namespace nestdaq {
// Host stages are disjoint. CUDA-event intervals are a separate view and must
// never be added to the host stages (they describe overlapping execution).
struct ExecProfileData {
#define LF_EXEC_HOST_FIELDS(X) \
    X(buffer_host_us) X(h2d_host_us) X(mark_launch_host_us) X(mark_sync_host_us) \
    X(subtr_launch_host_us) X(lut_launch_host_us) X(scan_host_us) X(compact_launch_host_us) \
    X(count_d2h_host_us) X(hits_resize_host_us) X(hits_d2h_host_us) \
    X(instrumentation_host_us) X(other_host_us) X(total_host_us)
#define LF_EXEC_DEVICE_FIELDS(X) \
    X(h2d_device_us) X(mark_device_us) X(subtr_device_us) X(lut_device_us) \
    X(scan_device_us) X(compact_device_us) X(count_d2h_device_us) X(hits_d2h_device_us)
#define LF_EXEC_COUNT_FIELDS(X) \
    X(mark_calls) X(h2d_bytes) X(count_d2h_bytes) X(hits_d2h_bytes) X(buffer_growths) X(mark_launches) X(metadata_h2d_bytes)
#define LF_EXEC_DOUBLE(name) double name{};
    LF_EXEC_HOST_FIELDS(LF_EXEC_DOUBLE)
    LF_EXEC_DEVICE_FIELDS(LF_EXEC_DOUBLE)
#undef LF_EXEC_DOUBLE
#define LF_EXEC_COUNT(name) std::uint64_t name{};
    LF_EXEC_COUNT_FIELDS(LF_EXEC_COUNT)
#undef LF_EXEC_COUNT
    void Add(const ExecProfileData& other) {
#define LF_EXEC_ADD(name) name += other.name;
        LF_EXEC_HOST_FIELDS(LF_EXEC_ADD)
        LF_EXEC_DEVICE_FIELDS(LF_EXEC_ADD)
        LF_EXEC_COUNT_FIELDS(LF_EXEC_ADD)
#undef LF_EXEC_ADD
    }
};
}
#endif
