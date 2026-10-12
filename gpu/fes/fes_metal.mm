// gpu/fes/fes_metal.mm -- Metal FES worker (Apple silicon). Reads a system in
// the shared contract (fes_io.hpp), runs the fes.metal compute kernel, and
// prints the solutions in the same format as fes_solve / fes_cuda, so the Rust
// `icx` binary drives it through one code path. Built by `make metal-worker`
// (macOS + xcrun + clang++); compile-checked in CI. A real run needs an Apple
// GPU; the Rust side re-verifies every returned solution regardless.
//
//   fes_metal --in <system-file>
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <libgen.h>
#include <string>
#include <vector>

#include "fes_io.hpp"

// Locate fes.metallib next to the executable, then in the CWD.
static NSString *metallib_path(const char *argv0) {
  std::string dir(argv0);
  char *d = dirname(&dir[0]);
  std::string beside = std::string(d) + "/fes.metallib";
  if ([[NSFileManager defaultManager]
          fileExistsAtPath:[NSString stringWithUTF8String:beside.c_str()]]) {
    return [NSString stringWithUTF8String:beside.c_str()];
  }
  return @"fes.metallib";
}

int main(int argc, char **argv) {
  @autoreleasepool {
    const char *in_path = nullptr;
    for (int i = 1; i < argc; ++i)
      if (!strcmp(argv[i], "--in") && i + 1 < argc) in_path = argv[++i];
    if (!in_path) {
      fprintf(stderr, "usage: %s --in <system-file>\n", argv[0]);
      return 2;
    }
    FesSystem sys;
    if (!fes_read_system(in_path, sys)) {
      fprintf(stderr, "fes_metal: could not parse %s\n", in_path);
      return 2;
    }

    id<MTLDevice> device = MTLCreateSystemDefaultDevice();
    if (!device) {
      fprintf(stderr, "fes_metal: no Metal device\n");
      return 3;
    }
    NSError *err = nil;
    id<MTLLibrary> lib =
        [device newLibraryWithURL:[NSURL fileURLWithPath:metallib_path(argv[0])]
                            error:&err];
    if (!lib) {
      fprintf(stderr, "fes_metal: cannot load fes.metallib: %s\n",
              err ? err.localizedDescription.UTF8String : "?");
      return 3;
    }
    id<MTLFunction> fn = [lib newFunctionWithName:@"fes_kernel"];
    id<MTLComputePipelineState> pso =
        [device newComputePipelineStateWithFunction:fn error:&err];
    id<MTLCommandQueue> queue = [device newCommandQueue];

    const int qlen = sys.n * (sys.n + 1) / 2;
    const int max_out = 1 << 20;
    // Shard the cube: up to 2^k threads over the top k variables.
    const int k = sys.n < 8 ? sys.n : 8;
    int low_bits = sys.n - k;
    if (low_bits < 0) low_bits = 0;
    const uint64_t nprefix = 1ULL << (sys.n - low_bits);

    auto shared = MTLResourceStorageModeShared;
    id<MTLBuffer> bLin = [device newBufferWithBytes:sys.lin.data()
                                             length:sizeof(uint64_t) * (sys.n ? sys.n : 1)
                                            options:shared];
    id<MTLBuffer> bQuad = [device newBufferWithBytes:sys.quad_tri.data()
                                              length:sizeof(uint64_t) * (qlen ? qlen : 1)
                                             options:shared];
    uint64_t cst = sys.cst;
    int n = sys.n, lb = low_bits, mo = max_out;
    id<MTLBuffer> bCst = [device newBufferWithBytes:&cst length:sizeof(uint64_t) options:shared];
    id<MTLBuffer> bN = [device newBufferWithBytes:&n length:sizeof(int) options:shared];
    id<MTLBuffer> bLow = [device newBufferWithBytes:&lb length:sizeof(int) options:shared];
    id<MTLBuffer> bOut = [device newBufferWithLength:sizeof(uint64_t) * max_out options:shared];
    uint32_t zero = 0;
    id<MTLBuffer> bCount = [device newBufferWithBytes:&zero length:sizeof(uint32_t) options:shared];
    id<MTLBuffer> bMaxOut = [device newBufferWithBytes:&mo length:sizeof(int) options:shared];

    id<MTLCommandBuffer> cb = [queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
    [enc setComputePipelineState:pso];
    [enc setBuffer:bLin offset:0 atIndex:0];
    [enc setBuffer:bQuad offset:0 atIndex:1];
    [enc setBuffer:bCst offset:0 atIndex:2];
    [enc setBuffer:bN offset:0 atIndex:3];
    [enc setBuffer:bLow offset:0 atIndex:4];
    [enc setBuffer:bOut offset:0 atIndex:5];
    [enc setBuffer:bCount offset:0 atIndex:6];
    [enc setBuffer:bMaxOut offset:0 atIndex:7];
    NSUInteger tgw = pso.maxTotalThreadsPerThreadgroup;
    if (tgw > nprefix) tgw = nprefix ? nprefix : 1;
    [enc dispatchThreads:MTLSizeMake(nprefix, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(tgw, 1, 1)];
    [enc endEncoding];
    [cb commit];
    [cb waitUntilCompleted];

    uint32_t count = *(uint32_t *)bCount.contents;
    uint32_t listed = count < (uint32_t)max_out ? count : (uint32_t)max_out;
    const uint64_t *o = (const uint64_t *)bOut.contents;
    std::vector<uint64_t> sols(o, o + listed);
    fes_write_solutions(sols);
    return 0;
  }
}
