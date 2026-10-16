"""IC advice bytes (PROTOCOL.md formula) for compact-orbit K, n."""
import sys
def ic_bytes(K, n):
    cand = K * K * n
    slots = 1 << (4 * cand - 1).bit_length()
    return 16 * slots + 24 * cand + 32 * 2 * K * n + 8 * K * n
if __name__ == "__main__":
    print(ic_bytes(int(sys.argv[1]), int(sys.argv[2])))
