#!/bin/bash

INFILE=$1
FLAGS=$2
if [[ -z $INFILE ]]; then
    echo "Usage: $0 <input_alignment> <flags>"
    exit 1
fi

INFILE=$(realpath "$INFILE")
OUTPREFIX=$(basename "$INFILE" | sed 's/\.[^.]*$//')


CFLAGS_BASE=(-O3 -fopenmp-simd -funsafe-math-optimizations)
CFLAGS_OPENMP=('-fopenmp' '-DOPENMP')
compile() {
    local out=$1
    shift
    if [[ $(stat --printf="%Y" "$out") -gt $(stat --printf="%Y" ../FastTree.c) ]]; then
        echo "Skipping $out, already compiled."
        return 2
    fi
    gcc -march=native "${CFLAGS_BASE[@]}" "$@" -o "$out" -lm ../FastTree.c
}

EXE=''
if [[ $(uname) == *_NT* ]]; then
    EXE='.exe'
fi

mkdir -p test
pushd test

TO_RUN=()
to_run() {
    TO_RUN+=("ft-$1$EXE $FLAGS > $OUTPREFIX.$1 2> err.$1 || rm $OUTPREFIX.$1")
}

compile_and_to_run() {
    local variant=$1
    shift
    compile "ft-head-$variant$EXE" "$@"
    local status=$?
    if [[ $status -eq 0 ]]; then
        to_run "head-$variant"
    elif [[ $status -eq 2 ]]; then
        if [[ -e $OUTPREFIX.head-$variant ]]; then
            echo "Skipping ft-head-$variant, not compiled and already run on this input."
        else
            to_run "head-$variant"
        fi
    else
        echo "Compilation failed for $variant"
    fi
}

cp ../FastTree$EXE ft-old$EXE &&
[[ -e $OUTPREFIX.old ]] || to_run old
hyperfine --input "$INFILE" --export-markdown "old.md" --runs 1 "${TO_RUN[@]}"
TO_RUN=()

compile_and_to_run 'base'
compile_and_to_run 'unroll' -funroll-loops
compile_and_to_run 'ipa' -fipa-reorder-for-locality -fipa-pta -fwhole-program
compile_and_to_run 'double' -DUSE_DOUBLE
compile_and_to_run 'sse3' -DUSE_SSE3
compile_and_to_run 'openmp' "${CFLAGS_OPENMP[@]}"

if (( ${#TO_RUN[@]})); then
    hyperfine --input "$INFILE" --export-markdown "head.md" --runs 3 "${TO_RUN[@]}"
fi
popd

