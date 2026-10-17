#!/bin/sh
set -eu

input=research/p192_extension_cover_20261012/type_i_guided_negative_fixed_origin_candidates.gp
output=research/p192_extension_cover_20261012/type_i_guided_negative_fixed_origin_reduction.json

# GP's integer vectors have JSON array syntax. Quote every decimal before jq
# parses it, so 192-bit field coefficients survive byte-for-byte as strings.
sed -e 's/^negative_cases=//' -e 's/;$//' "$input" |
  perl -pe 's/([0-9]+)/"$1"/g' |
  jq '{
    status: "independently_replayed_fixed_origin_guided_negative_classes",
    field_prime: "6277101735386680763835789423207666416083908700390324961279",
    class_identity: "[D-T]=-[B-A]",
    origin_base_indices: [0,1,2],
    classes_reduced: length,
    representatives: map({
      input_index: (.[0] | tonumber),
      original_representative_index: ((.[0] | tonumber) - 1),
      denominator_cubic_coefficients: .[1],
      intermediate_x_polynomial: .[2],
      intermediate_y_polynomial: .[3],
      numerator_cubic_coefficients: .[4],
      representative_x_polynomial: .[5],
      representative_y_polynomial: .[6],
      denominator_trial: (.[7] | map(tonumber))
    })
  }' > "$output"
