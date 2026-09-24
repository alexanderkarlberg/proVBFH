#!/usr/bin/env bash
#
# Run a command, showing its output as usual; if it fails, also report
# the most relevant lines of the output (lines mentioning errors, or
# else the last lines) as GitHub Actions error annotations, so that the
# cause of a failure is visible in the run summary.
#
# Usage: ci/annotate.sh <command> [args...]
#
set -uo pipefail

log=$(mktemp)
"$@" 2>&1 | tee "$log"
status=${PIPESTATUS[0]}
if [ "$status" -ne 0 ]; then
    lines=$(grep -iE '(^|[^a-z])(error|fatal|undefined|not found|failed)' "$log" | head -8)
    [ -n "$lines" ] || lines=$(tail -8 "$log")
    printf '%s\n' "$lines" | while IFS= read -r l; do
        # annotations are one line each; % and newlines must be escaped
        l=${l//'%'/'%25'}
        echo "::error title=$1 failed::$l"
    done
    echo "::error title=$1 failed::exit code $status (see the step log)"
fi
rm -f "$log"
exit "$status"
