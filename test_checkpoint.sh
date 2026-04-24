#!/bin/bash
# Test script for FastTree checkpoint/restart functionality.
# Verifies bit-identical output for single-threaded runs and
# successful completion for multi-threaded runs.
#
# Usage: ./test_checkpoint.sh [path_to_alignment]
# Default alignment: BigCOGs/COG1011.500.p

set -euo pipefail

ALN="${1:-BigCOGs/COG1011.500.p}"
FASTTREE="./FastTree"
TMPDIR=$(mktemp -d)
PASS=0
FAIL=0

cleanup() {
    rm -rf "$TMPDIR"
}
trap cleanup EXIT

fail_test() {
    echo "FAIL: $1"
    FAIL=$((FAIL + 1))
}

pass_test() {
    echo "PASS: $1"
    PASS=$((PASS + 1))
}

echo "=== Compiling FastTree ==="
gcc -DOPENMP -O3 -fopenmp -fopenmp-simd -funsafe-math-optimizations -march=native \
    -o "$FASTTREE" FastTree.c -lm
echo "Compiled OK"

echo ""
echo "=== Single-threaded tests (OMP_NUM_THREADS=1) ==="
echo "Using alignment: $ALN"
echo ""

# 1. Reference run (single-threaded, no checkpoint)
echo "--- 1. Reference run ---"
OMP_NUM_THREADS=1 "$FASTTREE" "$ALN" > "$TMPDIR/ref.tree" 2>"$TMPDIR/ref.stderr"
REF_TIME=$(grep "Total time:" "$TMPDIR/ref.stderr" | sed 's/Total time: \([0-9.]*\).*/\1/')
echo "Reference run complete: ${REF_TIME}s ($(wc -c < "$TMPDIR/ref.tree") bytes)"

# 2. Checkpoint run (single-threaded, full run) - should match reference
echo "--- 2. Checkpoint run (full, no restart) ---"
OMP_NUM_THREADS=1 "$FASTTREE" -checkpoint "$TMPDIR/ckpt_full.bin" "$ALN" \
    > "$TMPDIR/ckpt_full.tree" 2>"$TMPDIR/ckpt_full.stderr"
CKPT_TIME=$(grep "Total time:" "$TMPDIR/ckpt_full.stderr" | sed 's/Total time: \([0-9.]*\).*/\1/')
echo "Checkpoint run complete: ${CKPT_TIME}s"
if diff -q "$TMPDIR/ref.tree" "$TMPDIR/ckpt_full.tree" > /dev/null 2>&1; then
    pass_test "Checkpoint run matches reference (no restart)"
else
    fail_test "Checkpoint run differs from reference (no restart)"
fi

# 3. Restart from early checkpoint (after NJ + a few ME-NNI rounds)
echo "--- 3. Restart from early phase ---"
# Kill after 15% of reference time to get an early checkpoint
EARLY_TIMEOUT=$(echo "$REF_TIME" | awk '{t=int($1*0.15); if(t<15) t=15; print t}')
echo "  Timeout: ${EARLY_TIMEOUT}s"
timeout "$EARLY_TIMEOUT" bash -c \
    "OMP_NUM_THREADS=1 $FASTTREE -checkpoint $TMPDIR/ckpt_early.bin $ALN > /dev/null 2>$TMPDIR/ckpt_early.stderr" || true
if [ -f "$TMPDIR/ckpt_early.bin" ]; then
    PHASE=$(grep "Checkpoint saved" "$TMPDIR/ckpt_early.stderr" | tail -1 | sed 's/.*phase \([0-9]*\).*/\1/')
    ROUND=$(grep "Checkpoint saved" "$TMPDIR/ckpt_early.stderr" | tail -1 | sed 's/.*round \([0-9]*\).*/\1/')
    echo "  Checkpointed at phase $PHASE round $ROUND"
    OMP_NUM_THREADS=1 "$FASTTREE" -restart "$TMPDIR/ckpt_early.bin" "$ALN" \
        > "$TMPDIR/restart_early.tree" 2>"$TMPDIR/restart_early.stderr"
    if diff -q "$TMPDIR/ref.tree" "$TMPDIR/restart_early.tree" > /dev/null 2>&1; then
        pass_test "Restart from phase $PHASE round $ROUND matches reference"
    else
        fail_test "Restart from phase $PHASE round $ROUND differs from reference"
    fi
else
    fail_test "Early checkpoint was not created (NJ took longer than ${EARLY_TIMEOUT}s)"
fi

# 4. Restart from mid checkpoint (after ME phases, into ML)
echo "--- 4. Restart from mid phase ---"
# Kill after 50% of reference time to get a mid-run checkpoint
MID_TIMEOUT=$(echo "$REF_TIME" | awk '{t=int($1*0.50); if(t<60) t=60; print t}')
echo "  Timeout: ${MID_TIMEOUT}s"
timeout "$MID_TIMEOUT" bash -c \
    "OMP_NUM_THREADS=1 $FASTTREE -checkpoint $TMPDIR/ckpt_mid.bin $ALN > /dev/null 2>$TMPDIR/ckpt_mid.stderr" || true
if [ -f "$TMPDIR/ckpt_mid.bin" ]; then
    PHASE=$(grep "Checkpoint saved" "$TMPDIR/ckpt_mid.stderr" | tail -1 | sed 's/.*phase \([0-9]*\).*/\1/')
    ROUND=$(grep "Checkpoint saved" "$TMPDIR/ckpt_mid.stderr" | tail -1 | sed 's/.*round \([0-9]*\).*/\1/')
    echo "  Checkpointed at phase $PHASE round $ROUND"
    OMP_NUM_THREADS=1 "$FASTTREE" -restart "$TMPDIR/ckpt_mid.bin" "$ALN" \
        > "$TMPDIR/restart_mid.tree" 2>"$TMPDIR/restart_mid.stderr"
    if diff -q "$TMPDIR/ref.tree" "$TMPDIR/restart_mid.tree" > /dev/null 2>&1; then
        pass_test "Restart from phase $PHASE round $ROUND matches reference"
    else
        fail_test "Restart from phase $PHASE round $ROUND differs from reference"
    fi
else
    fail_test "Mid checkpoint was not created (timeout ${MID_TIMEOUT}s too short)"
fi

# 5. Restart from late checkpoint (last checkpoint of full run)
echo "--- 5. Restart from late phase ---"
PHASE=$(grep "Checkpoint saved" "$TMPDIR/ckpt_full.stderr" | tail -1 | sed 's/.*phase \([0-9]*\).*/\1/')
ROUND=$(grep "Checkpoint saved" "$TMPDIR/ckpt_full.stderr" | tail -1 | sed 's/.*round \([0-9]*\).*/\1/')
echo "  Last checkpoint at phase $PHASE round $ROUND"
OMP_NUM_THREADS=1 "$FASTTREE" -restart "$TMPDIR/ckpt_full.bin" "$ALN" \
    > "$TMPDIR/restart_late.tree" 2>"$TMPDIR/restart_late.stderr"
if diff -q "$TMPDIR/ref.tree" "$TMPDIR/restart_late.tree" > /dev/null 2>&1; then
    pass_test "Restart from phase $PHASE round $ROUND (late ML) matches reference"
else
    fail_test "Restart from phase $PHASE round $ROUND (late ML) differs from reference"
fi

echo ""
echo "=== Multi-threaded test ==="
echo ""

# 6. Multi-threaded checkpoint + restart (verify it completes without error)
echo "--- 6. Multi-threaded checkpoint run ---"
"$FASTTREE" -checkpoint "$TMPDIR/ckpt_mt.bin" "$ALN" \
    > "$TMPDIR/mt_full.tree" 2>"$TMPDIR/mt_full.stderr"
MT_PHASE=$(grep "Checkpoint saved" "$TMPDIR/mt_full.stderr" | tail -1 | sed 's/.*phase \([0-9]*\).*/\1/')
MT_TIME=$(grep "Total time:" "$TMPDIR/mt_full.stderr" | sed 's/Total time: \([0-9.]*\).*/\1/')
echo "  Completed in ${MT_TIME}s, last checkpoint at phase $MT_PHASE"
echo "  Tree size: $(wc -c < "$TMPDIR/mt_full.tree") bytes"

echo "--- 7. Multi-threaded restart ---"
"$FASTTREE" -restart "$TMPDIR/ckpt_mt.bin" "$ALN" \
    > "$TMPDIR/restart_mt.tree" 2>"$TMPDIR/restart_mt.stderr"
MT_RESTART_SIZE=$(wc -c < "$TMPDIR/restart_mt.tree")
if [ "$MT_RESTART_SIZE" -gt 0 ] && grep -q "Total time:" "$TMPDIR/restart_mt.stderr"; then
    pass_test "Multi-threaded restart completed successfully ($MT_RESTART_SIZE bytes)"
else
    fail_test "Multi-threaded restart failed or produced empty output"
fi
# Note: multi-threaded trees are not compared because OpenMP scheduling
# causes non-determinism between separate runs.

echo ""
echo "=== Results ==="
echo "Passed: $PASS"
echo "Failed: $FAIL"
echo ""

if [ "$FAIL" -gt 0 ]; then
    echo "SOME TESTS FAILED"
    exit 1
else
    echo "ALL TESTS PASSED"
    exit 0
fi
