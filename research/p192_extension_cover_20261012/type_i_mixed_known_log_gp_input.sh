#!/bin/sh
set -eu

input=research/p192_extension_cover_20261012/type_i_mixed_known_log_support.json
output=research/p192_extension_cover_20261012/type_i_mixed_known_log_gp_queries.gp

# The receipt stores all scalars and point coordinates as decimal strings.
# Transcribe them without passing 192-bit values through jq's number type.
jq -r '"steered=[" + ([.steered_known_only_queries[] |
  "[" + .scalar + "," + (.point_key|gsub(":";",")) + "]"]|join(",")) +
  "];\nfresh=[" + ([.independent_accepted_queries[].scalar]|join(",")) + "];"' \
  "$input" > "$output"
