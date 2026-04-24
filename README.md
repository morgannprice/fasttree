FastTree infers approximately-maximum-likelihood phylogenetic trees from alignments of nucleotide or protein sequences. FastTree can handle alignments with up to a million of sequences in a reasonable amount of time and memory.

See [documentation](https://morgannprice.github.io/fasttree/)
and
[compilation instructions](https://morgannprice.github.io/fasttree#install),
or download executables for
Linux ([double-precision](https://morgannprice.github.io/fasttree/FastTreeDbl)
or [single-precision](https://morgannprice.github.io/fasttree/FastTree))
or Windows
([single-precision command-line executable](https://morgannprice.github.io/fasttree/FastTree.exe))

## Checkpointing

For long-running analyses, FastTree can save checkpoints so that interrupted
runs can be resumed without starting over.

### Saving checkpoints

Use `-checkpoint` to periodically save state to a binary file:

    FastTree -checkpoint myrun.ckpt alignment.fasta > tree

Checkpoints are saved at natural boundaries between algorithm phases (after
neighbor-joining, after each NNI round, after each SPR round, etc.). The
checkpoint file is overwritten at each boundary, so it always contains the
most recent state. Writes are atomic (write to a temp file, then rename) so
the checkpoint is never corrupted by a crash.

Checkpointing has negligible overhead — it only writes a few MB at phase
boundaries, which take seconds compared to the minutes or hours of
computation between them.

### Resuming from a checkpoint

Use `-restart` to resume from a saved checkpoint:

    FastTree -restart myrun.ckpt alignment.fasta > tree

The same alignment file must be provided. FastTree will skip all completed
phases and resume where it left off. With single-threaded runs
(`OMP_NUM_THREADS=1`), the output is bit-identical to an uninterrupted run.
With multi-threaded runs, the tree will be fully optimized but may differ
slightly from an uninterrupted run due to OpenMP scheduling non-determinism.

### Testing

Run the checkpoint test suite:

    ./test_checkpoint.sh alignment.fasta
