#!/usr/bin/env bash
# Run OpenFOAM's checkMesh on each OpenFOAM case directory under the given
# root (as the synthetic suite writes them with SUITE_FOAM_DIR=<root>) and
# print one line per case: its result, cells, worst non-orthogonality and
# skewness, and any check that failed or warned.
#
# Where a case also holds `mesh_fluent.msh` (our Fluent file of the same
# mesh), the Fluent file goes round through OpenFOAM: `fluent3DMeshToFoam`
# converts it into the case's `fluent/` subcase, checkMesh checks that, and
# its counts (points, faces, internal faces, cells, patches, cell zones, and
# each patch's name and faces) must equal the case's own.
#
# A case holding `allowed_failures` lists checks it may fail (one a line;
# only `skewness` is known: the polyhedral duals of Milestone 21, which the
# user decided are written though checkMesh's skewness limit fails them at
# the wall). Such a case passes when checkMesh ran to its end and every
# check it failed is one of those; it is counted apart.
#
# Exits 1 if any case fails a check or a round trip (or a tool fails), 0
# otherwise.
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
round_trips=0
allowed_cases=0

# Whether a checkMesh log passes: `Mesh OK`, or (the case allowing it)
# failed checks that are all skewness.
passes() {
    local log="$1" case="$2"
    grep -q '^Mesh OK\.' "$log" && return 0
    [ -f "$case/allowed_failures" ] && grep -qx 'skewness' "$case/allowed_failures" || return 1
    grep -qE '^Failed [0-9]+ mesh checks\.' "$log" || return 1
    ! grep -E '^\s*\*\*\*' "$log" | grep -qv 'skew'
}

# The counts to compare from a checkMesh log: mesh stats and the patch table
# (names lowered, as Fluent zone names are lowercase).
counts() {
    grep -E '^ +(points|faces|internal faces|cells|boundary patches|cell zones):' "$1" | tr -s ' '
    awk '/^ +Patch +Faces +Points/ {on = 1; next} on && NF == 0 {exit} on && $1 != "\".*\"" {print tolower($1), $2}' "$1"
}

for case in "$root"/*/; do
    case="${case%/}"
    [ -f "$case/constant/polyMesh/owner" ] || continue
    total=$((total + 1))
    log="$case/log.checkMesh"
    checkMesh -case "$case" > "$log" 2>&1
    cells=$(awk '/^ *cells:/ {print $2; exit}' "$log")
    nonorth=$(awk '/Mesh non-orthogonality Max:/ {print $4; exit}' "$log")
    skew=$(awk '/Max skewness =/ {print $4; exit}' "$log")
    # Lines checkMesh marks: `***` a failed check, a leading ` *` a warning.
    marked=$(grep -E '^\s*\*' "$log" | sed -E 's/^\s+//' | cut -c1-80 | paste -sd'|' -)
    result=ok
    if grep -q '^Mesh OK\.' "$log"; then
        :
    elif passes "$log" "$case"; then
        result="ok (skewness allowed)"
        allowed_cases=$((allowed_cases + 1))
    else
        result=FAILED
    fi
    fluent=-
    if [ -f "$case/mesh_fluent.msh" ]; then
        round_trips=$((round_trips + 1))
        sub="$case/fluent"
        rm -rf "$sub"
        mkdir -p "$sub"
        cp -r "$case/system" "$sub/"
        if ! (cd "$case" && fluent3DMeshToFoam -case fluent mesh_fluent.msh > fluent/log.fluent3DMeshToFoam 2>&1); then
            fluent="conversion FAILED"
        elif ! { checkMesh -case "$sub" > "$sub/log.checkMesh" 2>&1; passes "$sub/log.checkMesh" "$case"; }; then
            fluent="checkMesh FAILED"
        elif ! diff <(counts "$log") <(counts "$sub/log.checkMesh") > "$sub/counts.diff"; then
            fluent="counts differ ($(paste -sd' ' "$sub/counts.diff" | cut -c1-120))"
        else
            fluent="fluent ok"
        fi
        case "$fluent" in *FAILED*|*differ*) result=FAILED ;; esac
    fi
    case "$result" in ok*) ;; *) failed=$((failed + 1)) ;; esac
    printf '%s\t%s\tcells %s\tnon-orthogonality %s\tskewness %s\t%s\t%s\n' \
        "$result" "$(basename "$case")" "${cells:-?}" "${nonorth:-?}" "${skew:-?}" "$fluent" "$marked"
done
echo "checkMesh: $((total - failed)) of $total cases OK ($allowed_cases of them failing only an allowed check; $round_trips with a Fluent round trip)"
[ "$total" -gt 0 ] && [ "$failed" -eq 0 ]
