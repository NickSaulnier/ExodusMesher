#!/usr/bin/env bash
# Run OpenFOAM's checkMesh on each OpenFOAM case directory under the given
# root (as the synthetic suite writes them with SUITE_FOAM_DIR=<root>) and
# print one line per case: its result, cells, worst non-orthogonality and
# skewness, and any check that failed or warned. Exits 1 if any case fails a
# check (or checkMesh itself fails), 0 otherwise.
#
# checkMesh must be on PATH (OpenFOAM's environment sourced). OpenFOAM is
# GPL: it runs here as a separate process and is never a dependency. Locally,
# in OpenCFD's image:
#
#   docker run --rm -v "$PWD/target/tmp/foam:/cases" -v "$PWD/tests/foam:/check" \
#       opencfd/openfoam-run:2512 bash /check/check.sh /cases
set -u
root="${1:?usage: check.sh <directory of cases>}"
failed=0
total=0
for case in "$root"/*/; do
    case="${case%/}"
    [ -f "$case/constant/polyMesh/owner" ] || continue
    total=$((total + 1))
    log="$case/log.checkMesh"
    checkMesh -case "$case" > "$log" 2>&1
    status=$?
    cells=$(awk '/^ *cells:/ {print $2; exit}' "$log")
    nonorth=$(awk '/Mesh non-orthogonality Max:/ {print $4; exit}' "$log")
    skew=$(awk '/Max skewness =/ {print $4; exit}' "$log")
    # Lines checkMesh marks: `***` a failed check, a leading ` *` a warning.
    marked=$(grep -E '^\s*\*' "$log" | sed -E 's/^\s+//' | cut -c1-80 | paste -sd'|' -)
    if [ $status -ne 0 ] || ! grep -q '^Mesh OK\.' "$log"; then
        result=FAILED
        failed=$((failed + 1))
    else
        result=ok
    fi
    printf '%s\t%s\tcells %s\tnon-orthogonality %s\tskewness %s\t%s\n' \
        "$result" "$(basename "$case")" "${cells:-?}" "${nonorth:-?}" "${skew:-?}" "$marked"
done
echo "checkMesh: $((total - failed)) of $total cases OK"
[ "$total" -gt 0 ] && [ "$failed" -eq 0 ]
