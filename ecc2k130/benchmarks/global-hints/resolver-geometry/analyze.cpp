#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>

namespace {
constexpr uint64_t kWorkers = 385024;
constexpr uint64_t kBatch = 16;
constexpr uint64_t kQueueEntries = kWorkers * kBatch;
constexpr uint64_t kQueueBytes = kQueueEntries * sizeof(uint32_t);
constexpr uint64_t kCounterBytes = sizeof(uint32_t);
constexpr uint64_t kSharedBytes = 57052;
constexpr uint64_t kSharedWords = kSharedBytes / sizeof(uint32_t);
constexpr uint64_t kSms = 188;
constexpr uint64_t kResolverBlocks = kSms;
constexpr uint64_t kRegistersPerThread = 124;
constexpr uint64_t kRegistersPerSm = 65536;

struct Geometry { unsigned threads; bool admitted; };
constexpr std::array<Geometry,4> kGeometries{{{64,false},{128,true},{256,true},{512,true}}};

uint64_t ceilDiv(uint64_t value, uint64_t divisor) {
    return (value + divisor - 1) / divisor;
}

void need(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    try {
        need(kSharedBytes % sizeof(uint32_t) == 0, "shared table is not word framed");
        need(kQueueEntries == 6160384, "queue entry ledger drift");
        need(kQueueBytes == 24641536, "queue byte ledger drift");
        need(kSharedWords == 14263, "shared table word ledger drift");
        need(kSharedBytes * kResolverBlocks == 10725776, "grid table-copy byte drift");

        std::cout << "threads,admitted,grid_lanes,warps_per_sm,registers_per_block,"
                     "register_headroom,table_words_per_lane,table_bytes_per_grid,"
                     "queue_entries,queue_bytes,counter_bytes\n";
        for (const Geometry geometry : kGeometries) {
            const uint64_t registers = kRegistersPerThread * geometry.threads;
            need(registers <= kRegistersPerSm, "geometry exceeds register file");
            std::cout << geometry.threads << ',' << (geometry.admitted ? 1 : 0) << ','
                      << kResolverBlocks * geometry.threads << ',' << geometry.threads / 32 << ','
                      << registers << ',' << kRegistersPerSm - registers << ','
                      << ceilDiv(kSharedWords, geometry.threads) << ','
                      << kSharedBytes * kResolverBlocks << ',' << kQueueEntries << ','
                      << kQueueBytes << ',' << kCounterBytes << '\n';
        }
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
