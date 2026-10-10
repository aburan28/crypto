#!/usr/bin/env bash
# Sparse checkout for building and running the harness (ecbench, cargo test)
# without the 4 GB of research archives it never reads. docs/sparse-checkout.md.
#
#   scripts/sparse-checkout.sh status
#   scripts/sparse-checkout.sh apply [--path P]...     # the harness profile
#   scripts/sparse-checkout.sh add --path P...         # widen the current profile
#   scripts/sparse-checkout.sh patterns [--path P]...  # print the rules, change nothing
#   scripts/sparse-checkout.sh disable                 # full checkout again
#   scripts/sparse-checkout.sh clone URL DIR [--branch B] [--path P]...
#
# The profile is everything except the large research/ directories below,
# plus every research/ path the Rust code needs: include_str!/include_bytes!
# targets anywhere (compile time) and research/ string literals in src/ and
# tests/ (read by cargo test). Those are recomputed from the sources on every
# apply. New research/<topic>_<date>/ directories, where ecbench sessions land,
# are included by default, so `git add` accepts them.
#
# Excluded files stay in HEAD and in the index (skip-worktree): `git ls-files`
# still lists them. Off disk is not absent; `add --path` materializes a path.
# Thin git orchestration only, per AGENTS.md ("Implementation language").
set -euo pipefail

# research/ directories over ~20 MB on 2026-10-09, largest first (4.2 GB).
HEAVY=(
  research/ic_candidate_tournament_20260915
  research/notes
  research/sat_factor_base_review_20260908
  research/boolean_halfword_20260923
  research/boolean_linear_tail_20260923
  research/boolean_byte_sieve_20260923
  research/ic_autolab_evidence_20260915
  research/boolean_linear_fibers_20260923
  research/boolean_projected_fibers_20260923
  research/boolean_affine_substitution_20260923
  research/boolean_initial_restriction_20260923
  research/boolean_packed_construction_20260922
  research/ic_tool_program
  research/boolean_construction_reduction_20260922
  research/boolean_envelope_indexed_20261002
  research/boolean_graded_tail_reuse_20261003
  research/rho_parity_20260915
  research/boolean_basis_transport_20260923
)
STATE="# sparse-checkout.sh"
SEG='[A-Za-z0-9_.@+=,-]'

die() { echo "sparse-checkout: $*" >&2; exit 2; }

# Collapse "." and ".." in a relative path.
normpath() {
  local IFS=/ part out=()
  for part in $1; do
    case "$part" in
      ''|.) ;;
      ..) [ ${#out[@]} -gt 0 ] && unset 'out[${#out[@]}-1]' ;;
      *) out+=("$part") ;;
    esac
  done
  echo "${out[*]}"
}

# research/ paths the Rust code compiles in or reads, one per line.
candidates() {
  local roots=() r file target
  for r in src tests benches examples build.rs; do [ -e "$r" ] && roots+=("$r"); done
  [ ${#roots[@]} -gt 0 ] || return 0
  grep -rnoE --include='*.rs' 'include_(str|bytes)!\([[:space:]]*"[^"]+"' "${roots[@]}" 2>/dev/null |
    while IFS= read -r line; do
      file=${line%%:*}
      target=${line#*\"}; target=${target%\"}
      normpath "$(dirname "$file")/$target"
    done
  for r in src tests; do
    [ -d "$r" ] && grep -rhoE --include='*.rs' "(\.\./)*research/$SEG+(/$SEG+)*" "$r" 2>/dev/null |
      sed -E 's#^(\.\./)+##'
  done
  return 0
}

# Rule lines for the tracked paths among the candidates: a file, or dir/.
referenced() {
  local tracked=$1 cand
  candidates | sed -E 's/[.,:;]+$//' | grep '^research/' | sort -u |
    while IFS= read -r cand; do
      if grep -qxF -- "$cand" "$tracked"; then
        echo "/$cand"
      elif awk -v p="${cand%/}/" 'index($0, p) == 1 { found = 1; exit } END { exit !found }' "$tracked"; then
        echo "/${cand%/}/"
      fi
    done
}

# Rules for a named path, re-including excluded directories beneath it
# (a parent rule does not override a more specific exclusion).
include() {
  local path=${1#/} h
  echo "/$path"
  case "$path" in */) ;; *) return 0 ;; esac
  for h in "${HEAVY[@]}"; do
    case "$h/" in "$path"?*) echo "/$h/" ;; esac
  done
}

render() {  # render [paths...]
  local tracked h p
  tracked=$(mktemp)
  git ls-files > "$tracked"
  [ -s "$tracked" ] || git ls-tree -r --name-only HEAD > "$tracked"
  echo "$STATE profile=harness"
  for p in "$@"; do echo "$STATE path=$p"; done
  echo "/*"
  for h in "${HEAVY[@]}"; do echo "!/$h/"; done
  echo "# research paths the Rust code compiles in or reads"
  referenced "$tracked"
  if [ $# -gt 0 ]; then
    echo "# named paths"
    for p in "$@"; do include "$p"; done
  fi
  rm -f "$tracked"
}

rules_file() { git rev-parse --git-path info/sparse-checkout; }

is_sparse() { [ "$(git config --bool core.sparseCheckout 2>/dev/null)" = true ]; }

saved_paths() {
  is_sparse || return 0
  sed -n "s/^$STATE path=//p" "$(rules_file)" 2>/dev/null || true
}

apply() {
  render "$@" | git sparse-checkout set --no-cone --stdin
  status
}

status() {
  local total skipped
  total=$(git ls-files | wc -l | tr -d ' ')
  if is_sparse; then
    skipped=$(git ls-files -t | grep -c '^S ' || true)
    if grep -q "^$STATE profile=" "$(rules_file)" 2>/dev/null; then
      echo "sparse checkout (profile harness): $((total - skipped)) of $total tracked files on disk, $skipped skipped"
    else
      echo "sparse checkout (rules not set by this script): $((total - skipped)) of $total tracked files on disk"
    fi
    saved_paths | sed 's/^/  + path /'
  else
    echo "full checkout: $total tracked files on disk"
  fi
  echo "history: $( [ "$(git rev-parse --is-shallow-repository)" = true ] && echo shallow || echo complete ); partial clone filter: $(git config remote.origin.partialclonefilter || echo none)"
}

parse_paths() {
  PATHS=()
  while [ $# -gt 0 ]; do
    case "$1" in
      --path) [ $# -ge 2 ] || die "--path needs a value"; PATHS+=("$2"); shift 2 ;;
      *) die "unexpected argument: $1" ;;
    esac
  done
}

cmd=${1:-status}
[ $# -gt 0 ] && shift

if [ "$cmd" = clone ]; then
  [ $# -ge 2 ] || die "usage: clone URL DIR [--branch B] [--path P]..."
  url=$1 dir=$2; shift 2
  branch=()
  if [ "${1:-}" = --branch ]; then
    [ $# -ge 2 ] || die "--branch needs a value"; branch=(--branch "$2"); shift 2
  fi
  parse_paths "$@"
  git clone --filter=blob:none --no-checkout ${branch[@]+"${branch[@]}"} "$url" "$dir"
  cd "$dir"
  # The base rules first, so the sources exist to compute the references from.
  { echo "/*"; for h in "${HEAVY[@]}"; do echo "!/$h/"; done; } |
    git sparse-checkout set --no-cone --stdin
  git checkout -f -q HEAD
  apply ${PATHS[@]+"${PATHS[@]}"}
  exit 0
fi

top=$(git rev-parse --show-toplevel 2>/dev/null) || die "not inside a git worktree"
cd "$top"
case "$cmd" in
  status) status ;;
  apply) parse_paths "$@"; apply ${PATHS[@]+"${PATHS[@]}"} ;;
  patterns) parse_paths "$@"; render ${PATHS[@]+"${PATHS[@]}"} ;;
  add)
    { is_sparse && grep -q "^$STATE profile=" "$(rules_file)"; } || die "no profile applied here; run 'apply' first"
    parse_paths "$@"
    merged=()
    while IFS= read -r p; do
      [ -n "$p" ] && merged+=("$p")
    done < <({ saved_paths; printf '%s\n' ${PATHS[@]+"${PATHS[@]}"}; } | sort -u)
    apply ${merged[@]+"${merged[@]}"}
    ;;
  disable) git sparse-checkout disable; status ;;
  -h|--help|help) sed -n '2,21p' "$0" | sed 's/^# \{0,1\}//' ;;
  *) die "unknown command '$cmd' (status, apply, add, patterns, disable, clone)" ;;
esac
