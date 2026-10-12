#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "kem.h"
#include "params.h"
int main(int argc, char **argv) {
  const char *op = argv[1]; int n = atoi(argv[2]);
  uint8_t pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES], ct[CRYPTO_CIPHERTEXTBYTES], ss[32], coins[64];
  memset(coins, 7, 64);
  crypto_kem_keypair_derand(pk, sk, coins);
  crypto_kem_enc_derand(ct, ss, pk, coins);
  for (int i = 0; i < n; i++) {
    coins[0] = (uint8_t)i;
    if (!strcmp(op, "keygen")) crypto_kem_keypair_derand(pk, sk, coins);
    else if (!strcmp(op, "encaps")) crypto_kem_enc_derand(ct, ss, pk, coins);
    else crypto_kem_dec(ss, ct, sk);
  }
  printf("%02x\n", ss[0]);
  return 0;
}
