#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <x86intrin.h>
#include "kem.h"
#include "params.h"
#define N 20000
int main(void) {
  uint8_t pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES], ct[CRYPTO_CIPHERTEXTBYTES], ss[32], coins[64];
  memset(coins, 7, 64);
  crypto_kem_keypair_derand(pk, sk, coins);
  crypto_kem_enc_derand(ct, ss, pk, coins);
  for (int w = 0; w < 1000; w++) crypto_kem_enc_derand(ct, ss, pk, coins);
  uint64_t t = __rdtsc();
  for (int i = 0; i < N; i++) { coins[0] = i; crypto_kem_keypair_derand(pk, sk, coins); }
  uint64_t kg = (__rdtsc() - t) / N;
  t = __rdtsc();
  for (int i = 0; i < N; i++) { coins[0] = i; crypto_kem_enc_derand(ct, ss, pk, coins); }
  uint64_t en = (__rdtsc() - t) / N;
  t = __rdtsc();
  for (int i = 0; i < N; i++) { crypto_kem_dec(ss, ct, sk); }
  uint64_t de = (__rdtsc() - t) / N;
  printf("keygen\t%lu\nencaps\t%lu\ndecaps\t%lu\n", kg, en, de);
  return ss[0] == 0xff;
}
