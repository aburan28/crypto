#!/bin/bash
# phase-medians.sh RUNS ARM...: per size and arm, the median over the arm's
# clean processes of cold time, build and collect, in ms, with the count.
RUNS=$1; shift
slug() {
  case $1 in
    k0n53) echo icv1-f2m53-tm56619371-dac20a85 ;;
    k1n59) echo icv1-f2m59-tm943548413-98844ecc ;;
    k0n61) echo icv1-f2m61-t158598901-ab42b6c5 ;;
  esac
}
printf "curve\tarm\tprocesses\tcold_ms\tbuild_ms\tcollect_ms\n"
for s in k0n53 k1n59 k0n61; do
  for arm in "$@"; do
    files=()
    for f in $RUNS/$arm/$s/*/r*.price.json; do
      [ -s "$f" ] || continue
      grep -q '"contended": false' ${f%.price.json}.isolation.jsonl 2>/dev/null || continue
      files+=("$f")
    done
    [ ${#files[@]} -gt 0 ] || continue
    jq -s -r --arg a $arm --arg c "$(slug $s)" --arg n ${#files[@]} '
      def med: sort | .[length/2|floor];
      [.[] | ([.repetitions[] | (.setup_ns + .ic_online.wall_ns)] | med)] as $cold
      | [.[] | ([.repetitions[] | .setup_phases_ns.build] | med)] as $build
      | [.[] | ([.repetitions[] | .setup_phases_ns.collect] | med)] as $coll
      | "\($c)\t\($a)\t\($n)\t\($cold | med / 1e6 * 10 | round / 10)\t\($build | med / 1e6 * 10 | round / 10)\t\($coll | med / 1e6 * 10 | round / 10)"' "${files[@]}"
  done
done
