/*
 * FastTree -- inferring minimum-evolution trees from multiple sequence alignments
 * by using profiles.
 *
 * Morgan N. Price, 2008-2009
 * http://www.microbesonline.org/fasttree/
 *
 *  Copyright (C) 2008 The Regents of the University of California
 *  All rights reserved.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  or visit http://www.gnu.org/copyleft/gpl.html
 *
 *  Disclaimer
 *
 *  NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY,
 *  NOR ANY OF THEIR EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
 *  OR ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY,
 *  COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT,
 *  OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
 *  PRIVATELY OWNED RIGHTS.
 */

/*
 * To compile FastTree, do:
 * cc -O2 -lm -o FastTree FastTree.c
 * Use -DTRACK_MEMORY if you want it to report its memory usage,
 * but results are not correct above 4GB because mallinfo stores int values
 *
 * To get usage guidance, do:
 * FastTree -help
 *
 * FastTree uses profiles instead of a distance matrix, and computes
 * support values for each split from the profiles of the 4 nodes
 * around the split. It stores a profile for each node and a average
 * profile over all active nodes (the "out-profile" for computing the
 * total sum of distance to other nodes).  The neighbor joining phase
 * requires O(N*L*a) space, where N is the number of sequences, L is
 * the alignment width, and a is the alphabet size. The top-hits
 * heuristic requires an additional O(N sqrt(N)) memory. After
 * neighbor-joining, FastTree improves the topology with
 * nearest-neighbor interchanges (NNIs), which does not have a
 * significant additional memory requirement.
 *
 * Although FastTree uses an after-the-fact log correction during the
 * NNIs, the operations on profiles require "additive" distances --
 * either %different or by using an amino acid similarity matrix
 * (-matrix option) If we are using %different as our distance matrix
 * then
 *
 * Profile_distance(A,B) = 1 - sum over characters of freq(A)*freq(B)
 *
 * and we can average this value over positions. Positions with gaps
 * are weighted by %ungapped(A) * %ungapped(B).
 *
 * If we are using an amino acid dissimilarity matrix D(i,j) then at
 * each position
 *
 * Profile_distance(A,B) = sum(i,j) freq(A==i) * freq(B==j) * D(i,j)
 * = sum(k) Ak * Bk * Lambda(k)
 *
 * where k iterates over 20 eigenvectors, Lambda(k) is the eigenvalue,
 * and if A==i, then Ak is the kth column of the inverse of the
 * eigenvector matrix.
 *
 * The exhaustive approach (-slow) takes O(N**3*L*a) time, but
 * this can be reduced to as little as O(N**(3/2)*log(N)*L*a) time
 * by using heuristics.
 *
 * It uses a combination of three heuristics: a visible set similar to
 * that of FastTree (Elias & Lagergren 2005), a local hill-climbing
 * search for a better join (as in relaxed neighbor-joining, Evans et
 * al. 2006), and a top-hit list to reduce the search space (see
 * below).
 *
 * The "visible" set stores, for each node, the best join for that
 * node, as identified at some point in the past
 *
 * If top-hits are not being used, then the method can be summarized
 * as:
 *
 * Compute the out-profile by averaging the leaves
 * Compute the out-distance of each leaf quickly, using the out-profile
 * Compute the visible set (or approximate it using top-hits, see below)
 * Until we're down to 3 active nodes:
 *   Find the best join in the visible set
 *	(This involves recomputing the neighbor-joining criterion,
 *      as out-distances and #active nodes may have changed)
 *   Follow a chain of best hits (again recomputing the criterion)
 *  	until we find a locally best join, as in relaxed neighbor joining
 *   Create a profile of the parent node, either using simple averages
 *	or using weighted joining as in BIONJ
 *   Update the out-profile and the out-distances
 *   Update the visible set:
 *      find the best join for the new joined node
 *      replace hits to the joined children with hits to the parent
 *      if we stumble across a join for the new node that is better
 *          than the corresponding entry in the visible set, "reset"
 *          that entry.
 *
 * For each iteration, this method does
 * O(N) work to find the best hit in the visible set
 * O(L*N*a*log(N)) work to do the local search, where log(N)
 *	is a pessimistic estimate of the number of iterations. In
 *      practice, we average <1 iteration for 2,000 sequences.
 * O(N*a) work to compute the joined profile and update the out-profile
 * O(L*N*a) work to update the out-distances
 * O(L*N*a) work to compare the joined profile to the other nodes
 *      (to find the new entry in the visible set)
 *
 * and there are N-3 iterations, so it takes O(N**2 * L * log(N) * a) time.
 *
 * The profile distances give exactly the same result as matrix
 * distances in neighbor-joining or BIONJ would if there are no gaps
 * in the alignment. If there are gaps, then it is an
 * approximation. To get the same result we also store a "diameter"
 * for each node (diameter is 0 for leaves).
 *
 * In the simpler case (NJ rather than BIONJ), when we join A and B to
 * give a new node AB,
 *
 * Profile(AB) = (A+B)/2
 * Profile_distance(AB,C) = (Profile_distance(A,C)+Profile_distance(B,C))/2
 * because the formulas above are linear
 *
 * And according to the neighor-joining rule,
 * d(AB,C) = (d(A,C)+d(B,C)-d(A,B))/2
 *
 * and we can achieve the same value by writing
 * diameter(AB) = pd(A,B)/2
 * diameter(leaf) = 0
 * d(A,B) = pd(A,B) - diameter(A) - diameter(B)
 *
 * because
 * d(AB,C) = (d(A,C)+d(B,C)-d(A,B))/2
 * = (pd(A,C)-diam(A)-diam(C)+pd(B,C)-diam(B)-diam(C)-d(A,B)+diam(A)+diam(B))/2
 * = (pd(A,C)+pd(B,C))/2 - diam(C) - pd(A,B)
 * = pd(AB,C) - diam(AB) - diam(C)
 *
 * If we are using BIONJ, with weight lambda for the join:
 * Profile(AB) = lambda*A + (1-lambda)*B
 * then a similar argument gives
 * diam(AB) = lambda*diam(A) + (1-lambda)*diam(B) + lambda*d(A,AB) + (1-lambda)*d(B,AB),
 *
 * where, as in neighbor joining,
 * d(A,AB) = d(A,B) + (total out_distance(A) - total out_distance(B))/(n-2)
 *
 * A similar recursion formula works for the "variance" matrix of BIONJ,
 * var(AB,C) = lambda*var(A,C) + (1-lambda)*var(B,C) - lambda*(1-lambda)*var(A,B)
 * is equivalent to
 * var(A,B) = pv(A,B) - vd(A) - vd(B), where
 * pv(A,B) = pd(A,B)
 * vd(A) = 0 for leaves
 * vd(AB) = lambda*vd(A) + (1-lambda)*vd(B) + lambda*(1-lambda)*var(A,B)
 *
 * The top-hist heuristic to reduce the work below O(N**2*L) stores a top-hit
 * list of size m=sqrt(N) for each active node.
 *
 * The list can be initialized for all the leaves in sub (N**2 * L) time as follows:
 * Pick a "seed" sequence and compare it to all others
 * Store the top m hits of the seed as its top-hit list
 * Take "close" hits of the seed(within the top m, and see the "close" parameter),
 *    and assume that their top m hits lie within the top 2*m hits of the seed.
 *    So, compare them to the seed's neighors (if they do not already
 *    have a top hit list) and set their top hits.
 *
 * This method does O(N*L) work for each seed, or O(N**(3/2)*L) work total.
 *
 * To avoid doing O(N*L) work at each iteration, we need to avoid
 * updating the visible set and the out-distances. So, we use "stale"
 * out-distances, and when searching the visible set for the best hit,
 * we only inspect the top m=sqrt(N) entries. We then update those
 * out-distances (up to 2*m*L*a work) and then find the best hit.
 *
 * To avoid searching the entire visible set, FastTree keeps
 * and updates a list of the top sqrt(N) entries in the visible set.
 * This costs O(sqrt(N)) time per join to find the best entry and to
 * update, or (N sqrt(N)) time overall.
 *
 * Similarly, when doing the local hill-climbing, we avoid O(N*L) work
 * by only considering the top-hits for the current node. So this adds
 * O(m*a*log(N)) work per iteration.
 *
 * When we join two nodes, we compute profiles and update the
 * out-profile as before. We need to compute the best hits of the node
 * -- we merge the lists for the children and select the best up-to-m
 * hits. If the top hit list contains a stale node we replace it with
 * its parent. If we still have <m/2 entries, we do a "refresh".
 *
 * In a "refresh", similar to the fast top-hit computation above, we
 * compare the "seed", in this case the new joined node, to all other
 * nodes. We compare its close neighbors (the top m hits) to all
 * neighbors (the top 2*m hits) and update the top-hit lists of all
 * neighbors (by merging to give a list of 3*m entries and then
 * selecting the best m entries).
 *
 * Finally, during these processes we update the visible sets for
 * other nodes with better hits if we find them, and we set the
 * visible entry for the new joined node to the best entry in its
 * top-hit list. (And whenever we update a visible entry, we
 * do O(sqrt(N)) work to update the top-visible list.)
 * These udpates are not common so they do not alter the
 * O(N sqrt(N) log(N) L a) total running time for the joining phase.
 */

#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h>
#ifdef TRACK_MEMORY
/* malloc.h apparently doesn't exist on MacOS */
#include <malloc.h>
#endif

char *usage =
  "Usage for FastTree version 1.1.0:\n"
  "  FastTree protein_alignment > tree\n"
  "  FastTree -nt nucleotide_alignment > tree\n"
  "  FastTree -nt < nucleotide_alignment > tree\n"
  "FastTree accepts alignments in fasta or phylip interleaved formats\n"
  "\n"
  "Common options (must be before the alignment file):\n"
  "  -quiet to suppress most reporting information\n"
  "  -pseudo to use pseudocounts (useful for highly gapped sequences)\n"
  "  -constraints constraintAlignment to constrain the topology search\n"
  "       constraintAlignment should have 1s or 0s to indicates splits\n"
  "  -noboot to not compute the 'local bootstrap'\n"
  "  -n <number> if the input has multiple alignments (phylip format only)\n"
  "  -intree newick_file to set the starting tree(s)\n"
  "  -intree1 newick_file to use this starting tree for all the alignments\n"
  "Use FastTree -expert to see detailed documentation and more options\n";

char *expertUsage =
  "FastTree [ -nt] [-n 100] [-pseudo | -pseudo 1.0]  [-boot 1000] [-quiet]\n"
  "           [-intree starting_trees_file | -intree1 starting_tree_file]\n"
  "           [-nni 10] [-spr 2] [-slow | -fastest] [-seed 1253] \n"
  "           [-top | -notop] [-topm 1.0 [-close 0.75] [-refresh 0.8]]\n"
  "           [-matrix Matrix | -nomatrix] [-nj | -bionj]\n"
  "           [-nt] \n"
  "           [ -constraints constraintAlignment [ -constraintWeight 10.0 ] ]\n"
  "         [ alignment_file ]\n"
  "        > newick_tree\n"
  "\n"
  "or\n"
  "\n"
  "FastTree [-nt] [-matrix Matrix | -nomatrix] [-rawdist] -makematrix [alignment]\n"
  "    [-n 100] > phylip_distance_matrix\n"
  "\n"
  "  FastTree supports fasta or phylip interleaved alignments\n"
  "  By default FastTree expects protein alignments,  use -nt for nucleotides\n"
  "  FastTree reads standard input if no alignment file is given\n"
  "\n"
  "  Use -n if you want to read multiple alignments in. This only\n"
  "  works with phylip interleaved format. For example, you can\n"
  "  use it with the output from phylip's seqboot. If you use -n, FastTree\n"
  "  will write 1 tree per line to standard output. You might also\n"
  "  want to use -quiet to eliminate status messages to standard error.\n"
  "  If you use -n together with -intree starting_tree_file,\n"
  "  then FastTree will also read that many trees from the file\n"
  "  (Use -intree1 if you want to use the same starting tree each time)\n"
  "\n"
  "Distances:\n"
  "  Default: For protein sequences, log-corrected distances and an\n"
  "     amino acid dissimilarity matrix derived from BLOSUM45\n"
  "  or for nucleotide sequences, Jukes-Cantor distances\n"
  "  To specify a different matrix, use -matrix FilePrefix or -nomatrix\n"
  "  Use -rawdist to turn the log-correction off\n"
  "  or to use %different instead of Jukes-Cantor\n"
  "\n"
  "  -pseudo [weight] -- Use pseudocounts to estimate distances between\n"
  "      sequences with little or no overlap. (Off by default.) Recommended\n"
  "      if analyzing the alignment has sequences with little or no overlap.\n"
  "      If the weight is not specified, it is 1.0\n"
  "\n"
  "Topology search:\n"
  "  By default, FastTree tries to improve the tree by doing 2*log2(N)\n"
  "  rounds of nearest-neighbor interchanges (NNI), where N is the number of\n"
  "  unique sequences, and 2 rounds of subtree-prune-regraft (SPR) moves\n"
  "  Use -nni to set the number of rounds of NNI and -spr to set the rounds of SPRs.\n"
  "  -nni 0 will turn off both NNIs and SPRs\n"
  "  If FastTree reports 'Bad splits: 0' then additional rounds of NNIs will probably\n"
  "  have no effect.\n"
  "  -sprlength set the maximum length of a SPR move (default 10)\n"
  "\n"
  "Support value options:\n"
  "  by default, FastTree computes a local bootstrap with 1,000 resamples.\n"
  "  The support values are proportions ranging from 0 to 1\n"
  "  The local bootstrap does not recompute topologies, so it is very fast\n"
  "\n"
  "  Use -noboot to turn off bootstrap or -boot 100 to use 100 resamples\n"
  "  Use -seed to initialize the random number generator\n"
  "\n"
  "Searching for the best join:\n"
  "  By default, FastTree combines the 'visible set' of fast neighbor-joining with\n"
  "      local hill-climbing as in relaxed neighbor-joining\n"
  "  -slow -- exhaustive search (like NJ or BIONJ, but different gap handling)\n"
  "      -slow takes half an hour instead of 8 seconds for 1,250 proteins\n"
  "  -fastest -- search the visible set (the top hit for each node) only\n"
  "      Unlike the original fast neighbor-joining, -fastest updates visible(C)\n"
  "      after joining A and B if join(AB,C) is better than join(C,visible(C))\n"
  "      -fastest also updates out-distances in a very lazy way\n"
  "      -fastest also sets the top-hits heuristic (-close) to be more aggressive\n"
  "\n"
  "Top-hit heuristics:\n"
  "  by default, FastTree uses a top-hit list to speed up search\n"
  "  use -notop (or -slow) to turn this feature off\n"
  "         and compare all leaves to each other,\n"
  "         and all new joined nodes to each other\n"
  "  -topm 1.0 -- set the top-hit list size to parameter*sqrt(N)\n"
  "         FastTree estimates the top m hits of a leaf from the\n"
  "         top 2*m hits of a 'close' neighbor, where close is\n"
  "         defined as d(seed,close) < 0.75 * d(seed, hit of rank 2*m),\n"
  "         and updates the top-hits as joins proceed\n"
  "  -close 0.75 -- modify the close heuristic, lower is more conservative\n"
  "  -refresh 0.8 -- compare a joined node to all other nodes if its\n"
  "         top-hit list is less than 80% of the desired length,\n"
  "         or if the age of the top-hit list is log2(m) or greater\n"
  "\n"
  "Join options:\n"
  "  -nj: regular (unweighted) neighbor-joining (default)\n"
  "  -bionj: weighted joins as in BIONJ\n"
  "          FastTree will also weight joins during NNIs\n"
  "\n"
  "Constrained topology search options:\n"
  "  -constraints alignmentfile -- an alignment with values of 0, 1, and -\n"
  "       Not all sequences need be present. A column of 0s and 1s defines a\n"
  "       constrained split. Some constraints may be violated\n"
  "       (see 'violating constraints:' in standard error).\n"
  "  -constraintWeight -- how strongly to weight the constraints. A value of 1\n"
  "       means a penalty of 1 in tree length for violating a constraint\n"
  "       Default: 10.0\n"
;


#define MAXCODES 20
#define NOCODE 127
/* Note -- sequence lines longer than BUFFER_SIZE are
   allowed, but FASTA header lines must be within this limit */
#define BUFFER_SIZE 5000

typedef struct {
  int nPos;
  int nSeq;
  char **names;
  char **seqs;
  int nSaved; /* actual allocated size of names and seqs */
} alignment_t;

/* For each position in a profile, we have a weight (% non-gapped) and a
   frequency vector. (If using a matrix, the frequency vector is in eigenspace).
   We also store codes for simple profile positions (all gaps or only 1 value)
   If weight[pos] > 0 && codes[pos] == NOCODE then we store the vector
   vectors itself is sets of nCodes long, so the vector for the ith nonconstant position
   starts at &vectors[nCodes*i]
   
   To speed up comparison of outprofile to a sequence or other simple profile, we also
   (for outprofiles) store codeDist[iPos*nCodes+k] = dist(k,profile[iPos])

   For constraints, we store a vector of nOn and nOff
   If not using constraints, those will be NULL
*/
typedef struct {
  /* alignment profile */
  float *weights;
  unsigned char *codes;
  float *vectors;		/* NULL if no non-constant positions, e.g. for leaves */
  int nVectors;
  float *codeDist;		/* Optional -- distance to each code at each position */

  /* constraint profile */
  int *nOn;
  int *nOff;
} profile_t;

/* A visible node is a pair of nodes i, j such that j is the best hit of i,
   using the neighbor-joining criterion, at the time the comparison was made,
   or approximately so since then.

   Note that variance = dist because in BIONJ, constant factors of variance do not matter,
   and because we weight ungapped sequences higher naturally when averaging profiles,
   so we do not take this into account in the computation of "lambda" for BIONJ.

   For the top-hit list heuristic, if the top hit list becomes "too short",
   we store invalid entries with i=j=-1 and dist/criterion very high.
*/
typedef struct {
  int i, j;
  float weight;			/* Total product of weights (maximum value is nPos) */
  float dist;			/* The uncorrected distance (includes diameter correction) */
  float criterion;		/* changes when we update the out-profile or change nActive */
} besthit_t;

typedef struct {
  int nChild;
  int child[3];
} children_t;

typedef struct {
  /* Distances between amino acids */
  float distances[MAXCODES][MAXCODES];

  /* Inverse of the eigenvalue matrix, for rotating a frequency vector
     into eigenspace so that profile similarity computations are
     O(alphabet) not O(alphabet*alphabet) time.
  */
  float eigeninv[MAXCODES][MAXCODES];
  float eigenval[MAXCODES];	/* eigenvalues */


  /* eigentot=eigeninv times the all-1s frequency vector
     useful for normalizing rotated frequency vectors
  */
  float eigentot[MAXCODES];	

  /* codeFreq is the transpose of the eigeninv matrix is
     the rotated frequency vector for each code */
  float codeFreq[MAXCODES][MAXCODES];
} distance_matrix_t;

typedef struct {
  /* The input */
  int nSeq;
  int nPos;
  char **seqs;			/* the aligment sequences array (not reallocated) */
  distance_matrix_t *distance_matrix; /* a pointer (not reallocated), or NULL if using %identity distance */
  /* Topological constraints are represented for each sequence as binary characters
     with values of '0', '1', or '-' (for missing data)
     Sequences that have no constraint may have a NULL string
  */
  int nConstraints;
  char **constraintSeqs;

  /* The profile data structures */
  int maxnode;			/* The next index to allocate */
  int maxnodes;			/* Space allocated in data structures below */
  profile_t **profiles;         /* Profiles of leaves and intermediate nodes */
  float *diameter;		/* To correct for distance "up" from children (if any) */
  float *varDiameter;		/* To correct variances for distance "up" */
  float *selfdist;		/* Saved for use in some formulas */
  float *selfweight;		/* Saved for use in some formulas */

  /* Average profile of all active nodes, the "outprofile"
   * If all inputs are ungapped, this has weight 1 (not nSequences) at each position
   * The frequencies all sum to one (or that is implied by the eigen-representation)
   */
  profile_t *outprofile;
  double totdiam;

  /* We sometimes use stale out-distances, so we remember what nActive was  */
  float *outDistances;		/* Sum of distances to other active (parent==-1) nodes */
  int *nOutDistActive;		/* What nActive was when this outDistance was computed */

  /* the inferred tree */
  int root;			/* index of the root. Unlike other internal nodes, it has 3 children */
  int *parent;			/* -1 or index of parent */
  children_t *child;
  float *branchlength;		/* Distance to parent */
  float *support;		/* 1 for high-confidence nodes */
} NJ_t;

/* Uniquify sequences in an alignment -- map from indices
   in the alignment to unique indicies in a NJ_t
*/
typedef struct {
  int nSeq;
  int nUnique;
  int *uniqueFirst;		/* iUnique -> iAln */
  int *alnNext;			/* iAln -> next, or -1  */
  int *alnToUniq;		/* iAln -> iUnique, or -1 if another was the exemplar */
  char **uniqueSeq;		/* indexed by iUniq -- points to strings allocated elsewhere */
} uniquify_t;

/* Describes which switch to do */
typedef enum {ABvsCD,ACvsBD,ADvsBC} nni_t;

/* A list of these describes a chain of NNI moves in a rooted tree,
   making up, in total, an SPR move
*/
typedef struct {
  int nodes[2];
  double deltaLength;		/* change in tree length for this step (lower is better) */
} spr_step_t;

/* Global variables */
/* Options */
int verbose = 1;
int slow = 0;
int fastest = 0;
int bionj = 0;
double tophitsMult = 1.0;	/* 0 means compare nodes to all other nodes */
double tophitsClose = -1.0;	/* Parameter for how close is close; also used as a coverage req. */
double tophitsRefresh = 0.8;	/* Refresh if fraction of top-hit-length drops to this */
double staleOutLimit = 0.01;	/* nActive changes by at most this amount before we recompute 
				   an out-distance. (Only applies if using the top-hits heuristic) */
double fResetOutProfile = 0.02;	/* Recompute out profile from scratch if nActive has changed
				   by more than this proportion, and */
int nResetOutProfile = 200;	/* nActive has also changed more than this amount */
int nBootstrap = 1000;		/* If set, number of replicates of local bootstrap to do */
int nCodes=20;			/* 20 if protein, 4 if nucleotide */
bool useMatrix=true;		/* If false, use %different as the uncorrected distance */
bool logdist = true;		/* If true, do a log-correction (scoredist-like or Jukes-Cantor)
				   but only during NNIs and support values, not during neighbor-joining */
double pseudoWeight = 0.0;      /* The weight of pseudocounts to avoid artificial long branches when
				   nearby sequences in the tree have little or no overlap
				   (off by default). The prior distance is based on
				   all overlapping positions among the quartet or triplet under
				   consideration. The log correction takes place after the
				   pseudocount is used. */
double constraintWeight = 10.0; /* Cost of violation of a topological constraint in evolutionary distance */

/* Performance and memory usage */
long profileOps = 0;		/* Full profile-based distance operations */
long outprofileOps = 0;		/* How many of profileOps are comparisons to outprofile */
long seqOps = 0;		/* Faster leaf-based distance operations */
long profileAvgOps = 0;		/* Number of profile-average steps */
long nBetter = 0;		/* Number of hill-climbing steps */
long nCloseUsed = 0;		/* Number of "close" neighbors we avoid full search for */
long nRefreshTopHits = 0;	/* Number of full-blown searches (interior nodes) */
long nVisibleReset = 0;		/* Number of resets of the visible set */
long nNNI = 0;			/* Number of NNI changes performed */
long nSPR = 0;			/* Number of SPR changes performed */
long nSuboptimalSplits = 0;	/* # of splits that are rejected given final tree (during bootstrap) */
long nSuboptimalConstrained = 0; /* Bad splits that are due to constraints */
long nConstraintViolations = 0;	/* Number of constraint violations */
long nProfileFreqAlloc = 0;
long nProfileFreqAvoid = 0;
long szAllAlloc = 0;
long mymallocUsed = 0;		/* useful allocations by mymalloc */
long maxmallocHeap = 0;		/* Maximum of mi.arena+mi.hblkhd from mallinfo (actual mem usage) */

/* Protein character set */
unsigned char *codesStringAA = (unsigned char*) "ARNDCQEGHILKMFPSTWYV";
unsigned char *codesStringNT = (unsigned char*) "ACGT";
unsigned char *codesString = NULL;

distance_matrix_t *ReadDistanceMatrix(char *prefix);
void SetupDistanceMatrix(/*IN/OUT*/distance_matrix_t *); /* set eigentot, codeFreq */
void ReadMatrix(char *filename, /*OUT*/float codes[MAXCODES][MAXCODES], bool check_codes);
void ReadVector(char *filename, /*OUT*/float codes[MAXCODES]);
alignment_t *ReadAlignment(/*READ*/FILE *fp); /* Returns a list of strings (exits on failure) */
alignment_t *FreeAlignment(alignment_t *); /* returns NULL */

NJ_t *InitNJ(char **sequences, int nSeqs, int nPos,
	     /*IN OPTIONAL*/char **constraintSeqs, int nConstraints,
	     /*IN OPTIONAL*/distance_matrix_t *); /* Allocates memory, initializes */
NJ_t *FreeNJ(NJ_t *NJ); /* returns NULL */
void FastNJ(/*IN/OUT*/NJ_t *NJ); /* Does the joins */
void ReliabilityNJ(/*IN/OUT*/NJ_t *NJ);	  /* Estimates the reliability of the joins */

/* One round of nearest-neighbor interchanges */
void NNI(/*IN/OUT*/NJ_t *NJ, int iRound, int nRounds); 

/* One round of subtree-prune-regraft moves */
void SPR(/*IN/OUT*/NJ_t *NJ, int maxSPRLength, int iRound, int nRounds);

void UpdateBranchLengths(/*IN/OUT*/NJ_t *NJ); /* Recomputes all branch lengths */
double TreeLength(/*IN/OUT*/NJ_t *NJ, bool recomputeProfiles); /*Also branch lengths*/

typedef struct {
  int nBadSplits;
  int nConstraintViolations;
  int nBadBoth;
  int nSplits;
} SplitCount_t;

void TestSplits(NJ_t *NJ, /*OUT*/SplitCount_t *splitcount);

/* Use out-profile and NJ->totdiam to recompute out-distance for node iNode
   Only does this computation if the out-distance is "stale" (nOutDistActive[iNode] != nActive)
 */
void SetOutDistance(/*IN/UPDATE*/NJ_t *NJ, int iNode, int nActive);

/* Always sets join->criterion; may update NJ->outDistance and NJ->nOutDistActive,
   assumes join's weight and distance are already set,
   and that the constraint penalty (if any) is included in the distance
*/
void SetCriterion(/*IN/UPDATE*/NJ_t *NJ, int nActive, /*IN/OUT*/besthit_t *join);

/* Compute the constraint penalty for a join. This is added to the "distance"
   by SetCriterion */
int JoinConstraintPenalty(/*IN*/NJ_t *NJ, int node1, int node2);
int JoinConstraintPenaltyPiece(NJ_t *NJ, int node1, int node2, int iConstraint);

/* Helper function for computing the number of constraints violated by
   a split, represented as counts of on and off on each side */
int SplitConstraintPenalty(int nOn1, int nOff1, int nOn2, int nOff2);

/* Computes weight and distance (which includes the constraint penalty)
   and then sets the criterion (maybe update out-distances)
*/
void SetDistCriterion(/*IN/UPDATE*/NJ_t *NJ, int nActive, /*IN/OUT*/besthit_t *join);

/* Reports the support for the (1,2) vs. (3,4) split
   col[iBoot*nPos+j] is column j for bootstrap iBoot
*/
double SplitSupport(profile_t *p1, profile_t *p2, profile_t *p3, profile_t *p4,
		    /*OPTIONAL*/distance_matrix_t *dmat,
		    int nPos,
		    int *col);

profile_t *SeqToProfile(/*IN/OUT*/NJ_t *NJ,
			char *seq, int nPos,
			/*OPTIONAL*/char *constraintSeqs, int nConstraints,
			int iNode,
			unsigned long counts[256]);

/* ProfileDist and SeqDist only set the dist and weight fields
   If using an outprofile, use the second argument of ProfileDist
   for better performance.

   These produce uncorrected distances.
*/
void ProfileDist(profile_t *profile1, profile_t *profile2, int nPos,
		 /*OPTIONAL*/distance_matrix_t *distance_matrix,
		 /*OUT*/besthit_t *hit);
void SeqDist(unsigned char *codes1, unsigned char *codes2, int nPos,
	     /*OPTIONAL*/distance_matrix_t *distance_matrix,
	     /*OUT*/besthit_t *hit);

/* Computes all pairs of profile distances, applies pseudocounts
   if pseudoWeight > 0, and applies log-correction if logdist is true.
   The lower index is compared to the higher index, e.g. for profiles
   A,B,C,D the comparison will be as in quartet_pair_t
*/
typedef enum {qAB,qAC,qAD,qBC,qBD,qCD} quartet_pair_t;
void CorrectedPairDistances(profile_t **profiles, int nProfiles,
			    /*OPTIONAL*/distance_matrix_t *distance_matrix,
			    int nPos,
			    /*OUT*/double *distances);

/* output is indexed by nni_t
   To ensure good behavior while evaluating a subtree-prune-regraft move as a series
   of nearest-neighbor interchanges, this uses a distance-ish model of constraints,
   as given by PairConstraintDistance(), rather than
   counting the number of violated splits (which is what FastTree does
   during neighbor-joining).
   Thus, penalty values may well be >0 even if no constraints are violated, but the
   relative scores for the three NNIs will be correct.
 */
void QuartetConstraintPenalties(profile_t *profiles[4], int nConstraints, /*OUT*/double d[3]);

double PairConstraintDistance(int nOn1, int nOff1, int nOn2, int nOff2);

/* the split is consistent with the constraint if any of the profiles have no data
   or if three of the profiles have the same uniform value (all on or all off)
   or if AB|CD = 00|11 or 11|00 (all uniform)
 */
bool SplitViolatesConstraint(profile_t *profiles[4], int iConstraint);

/* If false, no values were set because this constraint was not relevant.
   output is for the 3 splits
*/
bool QuartetConstraintPenaltiesPiece(profile_t *profiles[4], int iConstraint, /*OUT*/double penalty[3]);

/* Apply Jukes-Cantor or scoredist-like log(1-d) transform
   to correct the distance for multiple substitutions.
*/
double LogCorrect(double distance);

/* AverageProfile is used to do a weighted combination of nodes
   when doing a join. If weight is negative, then the value is ignored and the profiles
   are averaged. The weight is *not* adjusted for the gap content of the nodes.
   Also, the weight does not affect the representation of the constraints
*/
profile_t *AverageProfile(profile_t *profile1, profile_t *profile2,
			  int nPos, int nConstraints,
			  distance_matrix_t *distance_matrix,
			  double weight1);

/* Set a node's profile from its children.
   Deletes the previous profile if it exists
   Use -1.0 for a balanced join
   Fails unless the node has two children (e.g., no leaves or root)
*/
void SetProfile(/*IN/OUT*/NJ_t *NJ, int node, double weight1);

/* OutProfile does an unweighted combination of nodes to create the
   out-profile. It always sets code to NOCODE so that UpdateOutProfile
   can work.
*/
profile_t *OutProfile(profile_t **profiles, int nProfiles,
		      int nPos, int nConstraints,
		      distance_matrix_t *distance_matrix);

void UpdateOutProfile(/*UPDATE*/profile_t *out, profile_t *old1, profile_t *old2,
		      profile_t *new, int nActiveOld,
		      int nPos, int nConstraints,
		      distance_matrix_t *distance_matrix);

profile_t *NewProfile(int nPos, int nConstraints); /* returned has no vectors */
profile_t *FreeProfile(profile_t *profile, int nPos, int nConstraints); /* returns NULL */


/* f1 can be NULL if code1 != NOCODE, and similarly for f2
   Or, if (say) weight1 was 0, then can have code1==NOCODE *and* f1==NULL
   In that case, returns an arbitrary large number.
*/
double ProfileDistPiece(unsigned int code1, unsigned int code2,
			float *f1, float *f2, 
			/*OPTIONAL*/distance_matrix_t *dmat,
			/*OPTIONAL*/float *codeDist2);

/* Adds (or subtracts, if weight is negative) fIn/codeIn from fOut
   fOut is assumed to exist (as from an outprofile)
   do not call unless weight of input profile > 0
 */
void AddToFreq(/*IN/OUT*/float *fOut, double weight,
	       unsigned int codeIn, /*OPTIONAL*/float *fIn,
	       /*OPTIONAL*/distance_matrix_t *dmat);

/* Divide the vector (of length nCodes) by a constant
   so that the total (unrotated) frequency is 1.0 */
void NormalizeFreq(/*IN/OUT*/float *freq, distance_matrix_t *distance_matrix);

/* Allocate, if necessary, and recompute the codeDist*/
void SetCodeDist(/*IN/OUT*/profile_t *profile, int nPos, distance_matrix_t *dmat);

/* The allhits list contains the distances of the node to all other active nodes
   This is useful for the "reset" improvement to the visible set
   Note that the following routines do not handle the tophits heuristic
   and assume that out-distances are up to date.
*/
void SetBestHit(int node, NJ_t *NJ, int nActive,
		/*OUT*/besthit_t *bestjoin,
		/*OUT OPTIONAL*/besthit_t *allhits);
void ExhaustiveNJSearch(NJ_t *NJ, int nActive, /*OUT*/besthit_t *bestjoin);

/* Searches the visible set */
void FastNJSearch(NJ_t *NJ, int nActive, /*UPDATE*/besthit_t *visible, /*OUT*/besthit_t *bestjoin);

/* Subroutines for handling the tophits heuristic
   NJ may be modified because of updating of out-distances
*/

/* Before we do any joins -- sets tophits */
void SetAllLeafTopHits(NJ_t *NJ, int m, /*OUT*/besthit_t **tophits);

/* Find the best join to do.  topvisible is a list of active nodes (or
   -1) that stores the m known joins. If it shrinks
   to m/2 active entries then TopHitNJSearch will recompute it.
*/
void TopHitNJSearch(/*IN/UPDATE*/NJ_t *NJ,
		    int nActive,
		    int m,
		    /*IN/UPDATE*/besthit_t *visible,
		    /*IN/OUT*/int *topvisible, /* a list */
		    /*IN/OUT*/int *topVisibleAge, /* an integer */
		    /*IN/UPDATE*/besthit_t **tophits,
		    /*OUT*/besthit_t *bestjoin);

/* Returns the index of the best hit within the tophits.
   NJ may be modified because it updates out-distances if they are too stale
   tophits may be modified because it walks "up" from hits to joined nodes
   Does *not* update visible set
*/
int GetBestFromTopHits(int iNode, /*IN/UPDATE*/NJ_t *NJ, int nActive,
			/*IN/UPDATE*/besthit_t *tophits, /* of iNode */
		       int nTopHits);

/* visible set is modifiable so that we can reset it more globally when we do
   a "refresh", but we also set the visible set for newnode and do any
   "reset" updates too. And, we update many outdistances.
 */
void TopHitJoin(/*IN/UPDATE*/NJ_t *NJ, int newnode, int nActive, int m,
		/*IN/OUT*/besthit_t **tophits,
		/*IN/OUT*/int *tophitAge,
		/*IN/OUT*/besthit_t *visible,
		/*IN/OUT*/int *topvisible);

/* Sort by criterion and save the best nOut hits as a new array,
   which is returned.
   Does not update criterion or out-distances
   Ignores (silently removes) hit to self
   Pads the list with invalid entries so that it is always of length nOut
*/
besthit_t *SortSaveBestHits(/*IN/UPDATE*/besthit_t *besthits, int iNode, int nIn, int nOut);

/* Given candidate hits from one node, "transfer" them to another node:
   Stores them in a new place in the same order
   searches up to active nodes if hits involve non-active nodes
   If update flag is set, it also recomputes distance and criterion
   (and ensures that out-distances are updated)

 */
void TransferBestHits(/*IN/UPDATE*/NJ_t *NJ, int nActive,
		      int iNode,
		      /*IN*/besthit_t *oldhits,
		      int nOldHits,
		      /*OUT*/besthit_t *newhits,
		      bool updateDistance);

/* Given a top-hit list, look for improvements to the visible set of (j).
   Updates out-distances as it goes.
   If visible[j] is stale, then it set the current node to visible
   (visible[j] is usually stale because of the join that created this node...)
*/
void ResetVisible(/*IN/UPDATE*/NJ_t *NJ, int nActive,
		  /*IN*/besthit_t *tophitsNode,
		  int nTopHits,
		  /*IN/UPDATE*/besthit_t *visible,
		  /*IN/UPDATE*/int *topvisible);

/* Update the top-visible list to perhaps include visible[iNode] */
void UpdateTopVisible(/*IN*/NJ_t * NJ,
		      /*IN*/besthit_t *visible,
		      /*IN/UPDATE*/int *topvisible,
		      int m, 	/* length of topvisible */
		      int iNode);

/* Recompute the topvisible list */
void ResetTopVisible(/*IN*/NJ_t *NJ,
		     int nActive,
		     /*IN*/besthit_t *visible,
		     /*OUT*/int *topvisible,
		     int m);	/* length of topvisible */

/* Make a shorter list with only unique entries.
   Ignores self-hits to iNode and "dead" hits to
   nodes that have parents.
*/
besthit_t *UniqueBestHits(NJ_t *NJ, int iNode, besthit_t *combined, int nCombined, /*OUT*/int *nUniqueOut);

nni_t ChooseNNI(profile_t *profiles[4],
		/*OPTIONAL*/distance_matrix_t *dmat,
		int nPos, int nConstraints,
		/*OUT*/double criteria[3]);  /* The three potential internal branch lengths */

/* Returns the number of steps considered, with the actual steps in steps[]
   Modifies the tree by this chain of NNIs
*/
int FindSPRSteps(/*IN/OUT*/NJ_t *NJ, 
		 int node,
		 int parent,	/* sibling or parent of node to NNI to start the chain */
		 /*IN/OUT*/profile_t **upProfiles,
		 /*OUT*/spr_step_t *steps,
		 int maxSteps,
		 bool bFirstAC);

/* Undo a single NNI */
void UnwindSPRStep(/*IN/OUT*/NJ_t *NJ,
	       /*IN*/spr_step_t *step,
	       /*IN/OUT*/profile_t **upProfiles);


/* Update the profile of node and its ancestor, and delete nearby out-profiles */
void UpdateForNNI(/*IN/OUT*/NJ_t *NJ, int node, /*IN/OUT*/profile_t **upProfiles);

/* Sets NJ->parent[newchild] and replaces oldchild with newchild
   in the list of children of parent
*/
void ReplaceChild(/*IN/OUT*/NJ_t *NJ, int parent, int oldchild, int newchild);

int CompareHitsByCriterion(const void *c1, const void *c2);
int CompareHitsByJ(const void *c1, const void *c2);

int NGaps(NJ_t *NJ, int node);	/* only handles leaf sequences */

/* node is the parent of AB, sibling of C, nephew of D
   node cannot be root or a leaf
*/
void SetupABCD(NJ_t *NJ, int node,
	       /* the 4 profiles for ABCD; the last one is an outprofile */
	       /*OUT*/profile_t *profiles[4], 
	       /*IN/OUT*/profile_t **upProfiles,
	       /*OUT*/int nodeABC[3]);

int Sibling(NJ_t *NJ, int node); /* At root, no unique sibling so returns -1 */
void RootSiblings(NJ_t *NJ, int node, /*OUT*/int sibs[2]);

/* Print a progress report if more than 0.1 second has gone by since the progress report */
/* Format should include 0-4 %d references and no newlines */
void ProgressReport(char *format, int iArg1, int iArg2, int iArg3, int iArg4);

void *mymalloc(size_t sz);       /* Prints "Out of memory" and exits on failure */
void *myfree(void *, size_t sz); /* Always returns NULL */

void ran_start(long seed);
double knuth_rand();		/* Random number between 0 and 1 */

/* Like mymalloc; duplicates the input (returns NULL if given NULL) */
void *mymemdup(void *data, size_t sz);
void *myrealloc(void *data, size_t szOld, size_t szNew);

double pnorm(double z);		/* Probability(value <=z)  */

/* Hashtable functions */
typedef struct
{
  char *string;
  int nCount;			/* number of times this entry was seen */
  int first;			/* index of first entry with this value */
} hashbucket_t;

typedef struct {
  int nBuckets;
  /* hashvalue -> bucket. Or look in bucket + 1, +2, etc., till you hit a NULL string */
  hashbucket_t *buckets;
} hashstrings_t;
typedef int hashiterator_t;

hashstrings_t *MakeHashtable(char **strings, int nStrings);
hashstrings_t *FreeHashtable(hashstrings_t* hash); /*returns NULL*/
hashiterator_t FindMatch(hashstrings_t *hash, char *string);

/* Return NULL if we have run out of values */
char *GetHashString(hashstrings_t *hash, hashiterator_t hi);
int HashCount(hashstrings_t *hash, hashiterator_t hi);
int HashFirst(hashstrings_t *hash, hashiterator_t hi);

void PrintNJ(/*WRITE*/FILE *, NJ_t *NJ, char **names, uniquify_t *unique);
void PrintNJInternal(/*WRITE*/FILE *, NJ_t *NJ); /* Print topology using node indices as node names */

uniquify_t *UniquifyAln(/*IN*/alignment_t *aln);
uniquify_t *FreeUniquify(uniquify_t *);	/* returns NULL */

/* Convert a constraint alignment to a list of sequences. The returned array is indexed
   by iUnique and points to values in the input alignment
*/
char **AlnToConstraints(alignment_t *constraints, uniquify_t *unique, hashstrings_t *hashnames);

/* ReadTree ignores non-unique leaves after the first instance.
   At the end, it prunes the tree to ignore empty children and it
   unroots the tree if necessary.
*/
void ReadTree(/*IN/OUT*/NJ_t *NJ,
	      /*IN*/uniquify_t *unique,
	      /*IN*/hashstrings_t *hashnames,
	      /*READ*/FILE *fpInTree);
char *ReadTreeToken(/*READ*/FILE *fp); /* returns a static array, or NULL on EOF */
void ReadTreeAddChild(int parent, int child, /*IN/OUT*/int *parents, /*IN/OUT*/children_t *children);
/* Do not add the leaf if we already set this unique-set to another parent */
void ReadTreeMaybeAddLeaf(int parent, char *name,
			  hashstrings_t *hashnames, uniquify_t *unique,
			  /*IN/OUT*/int *parents, /*IN/OUT*/children_t *children);
void ReadTreeRemove(/*IN/OUT*/int *parents, /*IN/OUT*/children_t *children, int node);

/* Routines to support tree traversal and prevent visiting a node >1 time
   (esp. if topology changes).
*/
typedef bool *traversal_t;
traversal_t InitTraversal(NJ_t*);
traversal_t FreeTraversal(traversal_t, NJ_t*); /*returns NULL*/

/* returns new node, or -1 if nothing left to do. Use root for the first call.
   Will return every node and then root.
   Uses postorder tree traversal (depth-first search going down to leaves first)
*/
int TraversePostorder(int lastnode, NJ_t *NJ, /*IN/OUT*/traversal_t);

/* Routines to support storing up-profiles during tree traversal
   Eventually these should be smart enough to do weighted joins and
   to minimize memory usage
*/
profile_t **UpProfiles(NJ_t *NJ);
profile_t *GetUpProfile(/*IN/OUT*/profile_t **upProfiles, NJ_t *NJ, int node);
profile_t *DeleteUpProfile(/*IN/OUT*/profile_t **upProfiles, NJ_t *NJ, int node); /* returns NULL */
profile_t **FreeUpProfiles(profile_t **upProfiles, NJ_t *NJ); /* returns NULL */

/* Recomputes the profile for a node, presumably to reflect topology changes
   If bionj is set, does a weighted join -- which requires using upProfiles
 */
void RecomputeProfile(/*IN/OUT*/NJ_t *NJ, /*IN/OUT*/profile_t **upProfiles, int node);

/* If bionj is set, computes the weight to be given to A when computing the
   profile for the ancestor of A and B. C and D are the other profiles in the quartet
   If bionj is not set, returns -1 (which means unweighted in AverageProfile).
   (A and B are the first two profiles in the array)
*/
double QuartetWeight(profile_t *profiles[4], distance_matrix_t *dmat, int nPos);

/* Returns a list of nodes, starting with node and ending with root */
int *PathToRoot(NJ_t *NJ, int node, /*OUT*/int *depth);
int *FreePath(int *path, NJ_t *NJ); /* returns NULL */

/* The default amino acid distance matrix, derived from the BLOSUM45 similarity matrix */
distance_matrix_t matrixBLOSUM45;

int main(int argc, char **argv) {
  int nAlign = 1; /* number of alignments to read */
  int iArg;
  char *matrixPrefix = NULL;
  distance_matrix_t *distance_matrix = NULL;
  bool make_matrix = false;
  char *constraintsFile = NULL;
  char *intreeFile = NULL;
  bool intree1 = false;		/* the same starting tree each round */
  int nni = -1;			/* number of rounds of NNI, defaults to log2(n)+1 */
  int spr = 2;			/* number of rounds of SPR */
  int maxSPRLength = 10;	/* maximum distance to move a node */

  if (isatty(STDIN_FILENO) && argc == 1) {
      fprintf(stderr,"%s",usage);
      exit(0);
  }    
  for (iArg = 1; iArg < argc; iArg++) {
    if (strcmp(argv[iArg],"-makematrix") == 0) {
      make_matrix = true;
    } else if (strcmp(argv[iArg],"-logdist") == 0) {
      fprintf(stderr, "Warning: logdist is now on by default and obsolete\n");
    } else if (strcmp(argv[iArg],"-rawdist") == 0) {
      logdist = false;
    } else if (strcmp(argv[iArg],"-verbose") == 0 && iArg < argc-1) {
      verbose = atoi(argv[++iArg]);
    } else if (strcmp(argv[iArg],"-quiet") == 0) {
      verbose = 0;
    } else if (strcmp(argv[iArg],"-slow") == 0) {
      slow = 1;
    } else if (strcmp(argv[iArg],"-fastest") == 0) {
      fastest = 1;
    } else if (strcmp(argv[iArg], "-matrix") == 0 && iArg < argc-1) {
      iArg++;
      matrixPrefix = argv[iArg];
    } else if (strcmp(argv[iArg], "-nomatrix") == 0) {
      useMatrix = false;
    } else if (strcmp(argv[iArg], "-n") == 0 && iArg < argc-1) {
      iArg++;
      nAlign = atoi(argv[iArg]);
      if (nAlign < 1) {
	fprintf(stderr, "-n argument for #input alignments must be > 0 not %s\n", argv[iArg]);
	exit(1);
      }
    } else if (strcmp(argv[iArg], "-nt") == 0) {
      nCodes = 4;
    } else if (strcmp(argv[iArg], "-intree") == 0 && iArg < argc-1) {
      iArg++;
      intreeFile = argv[iArg];
    } else if (strcmp(argv[iArg], "-intree1") == 0 && iArg < argc-1) {
      iArg++;
      intreeFile = argv[iArg];
      intree1 = true;
    } else if (strcmp(argv[iArg], "-nj") == 0) {
      bionj = 0;
    } else if (strcmp(argv[iArg], "-bionj") == 0) {
      bionj = 1;
    } else if (strcmp(argv[iArg], "-boot") == 0 && iArg < argc-1) {
      iArg++;
      nBootstrap = atoi(argv[iArg]);
    } else if (strcmp(argv[iArg], "-noboot") == 0) {
      nBootstrap = 0;
    } else if (strcmp(argv[iArg], "-seed") == 0 && iArg < argc-1) {
      iArg++;
      long seed = atol(argv[iArg]);
      ran_start(seed);
    } else if (strcmp(argv[iArg],"-top") == 0) {
      if(tophitsMult < 0.01)
	tophitsMult = 1.0;
    } else if (strcmp(argv[iArg],"-notop") == 0) {
      tophitsMult = 0.0;
    } else if (strcmp(argv[iArg], "-topm") == 0 && iArg < argc-1) {
      iArg++;
      tophitsMult = atof(argv[iArg]);
    } else if (strcmp(argv[iArg], "-close") == 0 && iArg < argc-1) {
      iArg++;
      tophitsClose = atof(argv[iArg]);
      if (tophitsMult <= 0) {
	fprintf(stderr, "Cannot use -close unless -top is set above 0\n");
	exit(1);
      }
      if (tophitsClose <= 0 || tophitsClose >= 1) {
	fprintf(stderr, "-close argument must be between 0 and 1\n");
	exit(1);
      }
    } else if (strcmp(argv[iArg], "-refresh") == 0 && iArg < argc-1) {
      iArg++;
      tophitsRefresh = atof(argv[iArg]);
      if (tophitsMult <= 0) {
	fprintf(stderr, "Cannot use -refresh unless -top is set above 0\n");
	exit(1);
      }
      if (tophitsRefresh <= 0 || tophitsRefresh >= 1) {
	fprintf(stderr, "-refresh argument must be between 0 and 1\n");
	exit(1);
      }
    } else if (strcmp(argv[iArg],"-nni") == 0 && iArg < argc-1) {
      iArg++;
      nni = atoi(argv[iArg]);
      if (nni == 0)
	spr = 0;
    } else if (strcmp(argv[iArg],"-spr") == 0 && iArg < argc-1) {
      iArg++;
      spr = atoi(argv[iArg]);
    } else if (strcmp(argv[iArg],"-sprlength") == 0 && iArg < argc-1) {
      iArg++;
      maxSPRLength = atoi(argv[iArg]);
    } else if (strcmp(argv[iArg],"-help") == 0) {
      fprintf(stderr,"%s",usage);
      exit(0);
    } else if (strcmp(argv[iArg],"-expert") == 0) {
      fprintf(stderr, "%s", expertUsage);
      exit(0);
    } else if (strcmp(argv[iArg],"-pseudo") == 0) {
      if (iArg < argc-1 && isdigit(argv[iArg+1][0])) {
	iArg++;
	pseudoWeight = atof(argv[iArg]);
	if (pseudoWeight < 0.0) {
	  fprintf(stderr,"Illegal argument to -pseudo: %s\n", argv[iArg]);
	  exit(1);
	}
      } else {
	pseudoWeight = 1.0;
      }
    } else if (strcmp(argv[iArg],"-constraints") == 0 && iArg < argc-1) {
      iArg++;
      constraintsFile = argv[iArg];
    } else if (strcmp(argv[iArg],"-constraintWeight") == 0 && iArg < argc-1) {
      iArg++;
      constraintWeight = atof(argv[iArg]);
      if (constraintWeight <= 0.0) {
	fprintf(stderr, "Illegal argument to -constraintWeight (must be greater than zero): %s\n", argv[iArg]);
	exit(1);
      }
    } else if (argv[iArg][0] == '-') {
      fprintf(stderr, "Unknown or incorrect use of option %s\n%s", argv[iArg], usage);
      exit(1);
    } else
      break;
  }
  if(iArg < argc-1) {
    fprintf(stderr, usage);
    exit(1);
  }

  codesString = nCodes == 20 ? codesStringAA : codesStringNT;
  if (nCodes == 4 && matrixPrefix == NULL)
    useMatrix = false; 		/* no default nucleotide matrix */

  char *fileName = iArg == (argc-1) ?  argv[argc-1] : NULL;

  if (slow && fastest) {
    fprintf(stderr,"Cannot be both slow and fastest\n");
    exit(1);
  }
  if (slow && tophitsMult > 0) {
    tophitsMult = 0.0;
  }

  if (verbose && !make_matrix) {		/* Report settings */
    char tophitString[100] = "no";
    char tophitsCloseStr[100] = "default";
    if(tophitsClose > 0) sprintf(tophitsCloseStr,"%.2f",tophitsClose);
    if(tophitsMult>0) sprintf(tophitString,"%.2f*sqrtN close=%s refresh=%.2f",
			      tophitsMult, tophitsCloseStr, tophitsRefresh);
    char supportString[100] = "none";
    if (nBootstrap>0) sprintf(supportString,"Local boot %d",nBootstrap);
    char nniString[100] = "(no NNI)";
    if (nni > 0)
      sprintf(nniString, "+NNI (%d rounds)", nni);
    if (nni == -1)
      strcpy(nniString, "+NNI");
    char sprString[100] = "(no SPR)";
    if (spr > 0)
      sprintf(sprString, "+SPR (%d rounds, max-distance %d)", spr, maxSPRLength);

    fprintf(stderr,"Alignment: %s", fileName != NULL ? fileName : "standard input");
    if (nAlign>1)
      fprintf(stderr, " (%d alignments)", nAlign);
    fprintf(stderr,"\n%s distances: %s Joins: %s Support: %s\n",
	    nCodes == 20 ? "Amino acid" : "Nucleotide",
	    matrixPrefix ? matrixPrefix : (useMatrix? "BLOSUM45"
					   : (nCodes==4 && logdist ? "Jukes-Cantor" : "%different")),
	    bionj ? "weighted" : "balanced" ,
	    supportString);
    if (intreeFile == NULL)
      fprintf(stderr, "Search: %s %s %s\nTopHits: %s\n",
	      slow?"Exhaustive (slow)" : (fastest ? "Fastest" : "Normal (fast/relax)"),
	      nniString, sprString,
	      tophitString);
    else
      fprintf(stderr, "Start at tree from %s %s %s\n", intreeFile, nniString, sprString);

    if (constraintsFile != NULL)
      fprintf(stderr, "Constraints: %s Weight: %.3f\n", constraintsFile, constraintWeight);
    if (pseudoWeight > 0)
      fprintf(stderr, "Pseudocount weight for comparing sequences with little overlap: %.3lf\n",pseudoWeight);
  }

  if (matrixPrefix != NULL) {
    if (!useMatrix) {
      fprintf(stderr,"Cannot use both -matrix and -nomatrix arguments!");
      exit(1);
    }
    distance_matrix = ReadDistanceMatrix(matrixPrefix);
  } else if (useMatrix) { 	/* use default matrix */
    assert(nCodes==20);
    distance_matrix = &matrixBLOSUM45;
    SetupDistanceMatrix(distance_matrix);
  } else {
    distance_matrix = NULL;
  }

  int iAln;
  FILE *fpIn = fileName != NULL ? fopen(fileName, "r") : stdin;
  if (fpIn == NULL) {
    fprintf(stderr, "Cannot read %s\n", fileName);
    exit(1);
  }
  FILE *fpConstraints = NULL;
  if (constraintsFile != NULL) {
    fpConstraints = fopen(constraintsFile, "r");
    if (fpConstraints == NULL) {
      fprintf(stderr, "Cannot read %s\n", constraintsFile);
      exit(1);
    }
  }

  FILE *fpInTree = NULL;
  if (intreeFile != NULL) {
    fpInTree = fopen(intreeFile,"r");
    if (fpInTree == NULL) {
      fprintf(stderr, "Cannot read %s\n", intreeFile);
      exit(1);
    }
  }

  for(iAln = 0; iAln < nAlign; iAln++) {
    alignment_t *aln = ReadAlignment(fpIn);
    if (aln->nSeq < 1) {
      fprintf(stderr, "No alignment sequences\n");
      exit(1);
    }

    clock_t clock_start = clock();
    ProgressReport("Read alignment",0,0,0,0);

    /* Check that all names in alignment are unique */
    hashstrings_t *hashnames = MakeHashtable(aln->names, aln->nSeq);
    int i;
    for (i=0; i<aln->nSeq; i++) {
      hashiterator_t hi = FindMatch(hashnames,aln->names[i]);
      if (HashCount(hashnames,hi) != 1) {
	fprintf(stderr,"Non-unique name %s in the alignment\n",aln->names[i]);
	exit(1);
      }
    }

    /* Make a list of unique sequences -- note some lists are bigger than required */
    ProgressReport("Hashed the names",0,0,0,0);
    if (make_matrix) {
      NJ_t *NJ = InitNJ(aln->seqs, aln->nSeq, aln->nPos,
			/*constraintSeqs*/NULL, /*nConstraints*/0,
			distance_matrix);
      printf("   %d\n",aln->nSeq);
      int i,j;
      for(i = 0; i < NJ->nSeq; i++) {
	printf("%s",aln->names[i]);
	for (j = 0; j < NJ->nSeq; j++) {
	  besthit_t hit;
	  SeqDist(NJ->profiles[i]->codes,NJ->profiles[j]->codes,NJ->nPos,NJ->distance_matrix,/*OUT*/&hit);
	  if (logdist)
	    hit.dist = LogCorrect(hit.dist);
	  printf(" %f", hit.dist);
	}
	printf("\n");
      }
    } else {
      /* reset counters*/
      profileOps = 0;
      outprofileOps = 0;
      seqOps = 0;
      profileAvgOps = 0;
      nBetter = 0;
      nCloseUsed = 0;
      nRefreshTopHits = 0;
      nVisibleReset = 0;
      nNNI = 0;
      nProfileFreqAlloc = 0;
      nProfileFreqAvoid = 0;
      szAllAlloc = 0;
      mymallocUsed = 0;
      maxmallocHeap = 0;

      uniquify_t *unique = UniquifyAln(aln);
      ProgressReport("Identified unique sequences",0,0,0,0);

      /* read constraints */
      alignment_t *constraints = NULL;
      char **uniqConstraints = NULL;
      if (constraintsFile != NULL) {
	constraints = ReadAlignment(fpConstraints);
	if (constraints->nSeq < 4) {
	  fprintf(stderr, "Warning: constraints file with less than 4 sequences ignored:\nalignment #%d in %s\n",
		  iAln+1, constraintsFile);
	  constraints = FreeAlignment(constraints);
	} else {
	  uniqConstraints = AlnToConstraints(constraints, unique, hashnames);
	  ProgressReport("Read the constraints",0,0,0,0);
	}
      }	/* end load constraints */

      NJ_t *NJ = InitNJ(unique->uniqueSeq, unique->nUnique, aln->nPos,
			uniqConstraints,
			uniqConstraints != NULL ? constraints->nPos : 0, /* nConstraints */
			distance_matrix);
      if (verbose>1) fprintf(stderr, "read %s seqs %d (%d unique) positions %d nameLast %s seqLast %s\n",
			     fileName ? fileName : "standard input",
			     aln->nSeq, unique->nUnique, aln->nPos, aln->names[aln->nSeq-1], aln->seqs[aln->nSeq-1]);
      if (fpInTree != NULL) {
	if (intree1)
	  fseek(fpInTree, 0L, SEEK_SET);
	ReadTree(/*IN/OUT*/NJ, /*IN*/unique, /*IN*/hashnames, /*READ*/fpInTree);
	if (verbose > 1)
	  fprintf(stderr, "Read tree from %s\n", intreeFile);
	if (verbose >= 2)
	  PrintNJ(stderr, NJ, aln->names, unique);
      } else {
	FastNJ(NJ);
      }

      /* profile-frequencies for the "up-profiles" in ReliabilityNJ take only diameter(Tree)*L*a
	 space not N*L*a space, because we can free them as we go.
	 And up-profile by their nature tend to be complicated.
	 So save the profile-frequency memory allocation counters now to exclude later results.
      */
#ifdef TRACK_MEMORY
      long svProfileFreqAlloc = nProfileFreqAlloc;
      long svProfileFreqAvoid = nProfileFreqAvoid;
#endif
      int nniToDo = nni == -1 ? 1 + (int)(0.5 + 2.0 * log(NJ->nSeq)/log(2)) : nni;
      int sprRemaining = spr;
      if(verbose>0) {
	if (isatty(STDERR_FILENO)) fprintf(stderr,"\n");
	if (fpInTree == NULL)
	    fprintf(stderr, "Initial topology in %.2f seconds\n", (clock()-clock_start)/(double)CLOCKS_PER_SEC);
	if (spr > 0 || nniToDo > 0)
	  fprintf(stderr,"Refining topology: %d rounds NNIs, %d rounds SPRs\n", nniToDo, spr);
	}

      if (nniToDo>0) {
	int i;
	for (i=0; i < nniToDo; i++) {
	  if(verbose>1) {
	    fprintf(stderr, "Topology before NNI round %d\n",i);
	    fflush(stderr);
	    PrintNJ(stderr, NJ, aln->names, unique);
	  }
	  NNI(/*IN/OUT*/NJ, i, nniToDo);
	  /* Interleave SPRs with NNIs (typically 1/3rd NNI, SPR, 1/3rd NNI, SPR, 1/3rd NNI */
	  if (sprRemaining > 0 && nniToDo/(spr+1) > 0 && ((i+1) % (nniToDo/(spr+1))) == 0) {
	    SPR(/*IN/OUT*/NJ, maxSPRLength, spr-sprRemaining, spr);
	    sprRemaining--;
	  }
	}
      }
      while(sprRemaining > 0) {	/* do any remaining SPR rounds */
	SPR(/*IN/OUT*/NJ, maxSPRLength, spr-sprRemaining, spr);
	sprRemaining--;
      }

      /* Always update branch lengths (even if nni=0) so that they are log-corrected,
	 do not include penalties from constraints,
	 and avoid errors due to approximation of out-distances.
      */
      UpdateBranchLengths(/*IN/OUT*/NJ);
      SplitCount_t splitcount;
      TestSplits(NJ, /*OUT*/&splitcount);
      if (nBootstrap > 0) {
	if(verbose>0) {
	  double total_len = 0;
	  int iNode;
	  for (iNode = 0; iNode < NJ->maxnode; iNode++)
	    total_len += fabs(NJ->branchlength[iNode]);
	  if (isatty(STDERR_FILENO)) fprintf(stderr,"\n");
	  fprintf(stderr, "Total branch-length %.3f after %.2f sec -- computing support values\n",
		  total_len,
		  (clock()-clock_start)/(double)CLOCKS_PER_SEC);
	}
	fflush(stderr);
	ReliabilityNJ(NJ);
      }
      if(verbose) {
	if (isatty(STDERR_FILENO)) fprintf(stderr,"\n");
	fprintf(stderr, "Unique: %d/%d", NJ->nSeq, aln->nSeq);
	if(!slow) fprintf(stderr, " Hill-climb: %ld Update-best: %ld", nBetter, nVisibleReset);
	if (nniToDo > 0 || spr > 0)
	  fprintf(stderr, " NNI: %ld SPR: %ld", nNNI, nSPR);
	fprintf(stderr,"\n");
	if (nCloseUsed>0 || nRefreshTopHits>0)
	  fprintf(stderr, "Top hits: close neighbors %ld/%d refreshes %ld\n",
		  nCloseUsed, NJ->nSeq, nRefreshTopHits);
	if (NJ->nSeq > 3) {
	  if (NJ->nConstraints == 0) {
	    fprintf(stderr, "Bad splits: %d of %d\n", splitcount.nBadSplits, splitcount.nSplits);
	  } else {
	    fprintf(stderr, "Bad splits: %d of %d violating constraints: %d both bad: %d\n",
		    splitcount.nBadSplits,
		    splitcount.nSplits,
		    splitcount.nConstraintViolations,
		    splitcount.nBadBoth
		    );
	  }
	}
	double dN2 = NJ->nSeq*(double)NJ->nSeq;
	fprintf(stderr, "Time %.2f Dist/N**2: by-profile %.3f (out %.3f) by-leaf %.3f Avg %.3f\n",
		(clock()-clock_start)/(double)CLOCKS_PER_SEC,
		profileOps/dN2, outprofileOps/dN2, seqOps/dN2, profileAvgOps/dN2);
#ifdef TRACK_MEMORY
	fprintf(stderr, "Memory: %.2f MB (%.1f byte/pos) ",
		maxmallocHeap/1.0e6, maxmallocHeap/(double)(aln->nSeq*(double)aln->nPos));
	/* Only report numbers from before we do reliability estimates */
	fprintf(stderr, "profile-freq-alloc %ld avoided %.2f%%\n", 
		svProfileFreqAlloc,
		svProfileFreqAvoid > 0 ?
		100.0*svProfileFreqAvoid/(double)(svProfileFreqAlloc+svProfileFreqAvoid)
		: 0);
#endif
      }
      fflush(stderr);
      PrintNJ(stdout, NJ, aln->names, unique);
      fflush(stdout);
      FreeNJ(NJ);
      if (uniqConstraints != NULL)
	uniqConstraints = myfree(uniqConstraints, sizeof(char*) * unique->nUnique);
      constraints = FreeAlignment(constraints);
      unique = FreeUniquify(unique);
    } /* end build tree */
    hashnames = FreeHashtable(hashnames);
    aln = FreeAlignment(aln);
  } /* end loop over alignments */
  exit(0);
}

void ProgressReport(char *format, int i1, int i2, int i3, int i4) {
  static bool clock_set = false;
  static clock_t clock_last;
  static clock_t clock_begin;

  if (!verbose)
    return;

  clock_t clock_now = clock();
  if (!clock_set) {
    clock_begin = clock_last = clock_now;
    clock_set = true;
  }
  if ((clock_now-clock_last)/(double)CLOCKS_PER_SEC > 0.1 || verbose > 1) {
    fprintf(stderr, "%9.2f seconds: ", (clock_now-clock_begin)/(double)CLOCKS_PER_SEC);
    fprintf(stderr, format, i1, i2, i3, i4);
    if (verbose > 1 || !isatty(STDERR_FILENO)) {
      fprintf(stderr, "\n");
    } else {
      fprintf(stderr, "   \r");
    }
    fflush(stderr);
    clock_last = clock_now;
  }
}

NJ_t *InitNJ(char **sequences, int nSeq, int nPos,
	     /*OPTIONAL*/char **constraintSeqs, int nConstraints,
	     /*OPTIONAL*/distance_matrix_t *distance_matrix) {
  int iNode;

  NJ_t *NJ = (NJ_t*)mymalloc(sizeof(NJ_t));
  NJ->root = -1; 		/* set at end of FastNJ() */
  NJ->maxnode = NJ->nSeq = nSeq;
  NJ->nPos = nPos;
  NJ->maxnodes = 2*nSeq;
  NJ->seqs = sequences;
  NJ->distance_matrix = distance_matrix;
  NJ->nConstraints = nConstraints;
  NJ->constraintSeqs = constraintSeqs;

  NJ->profiles = (profile_t **)mymalloc(sizeof(profile_t*) * NJ->maxnodes);

  unsigned long counts[256];
  int i;
  for (i = 0; i < 256; i++)
    counts[i] = 0;
  for (iNode = 0; iNode < NJ->nSeq; iNode++) {
    NJ->profiles[iNode] = SeqToProfile(NJ, NJ->seqs[iNode], nPos,
				       constraintSeqs != NULL ? constraintSeqs[iNode] : NULL,
				       nConstraints,
				       iNode,
				       /*IN/OUT*/counts);
  }
  unsigned long totCount = 0;
  for (i = 0; i < 256; i++)
    totCount += counts[i];

  /* warnings about unknown characters */
  for (i = 0; i < 256; i++) {
    if (counts[i] == 0 || i == '.' || i == '-')
      continue;
    unsigned char *codesP;
    bool bMatched = false;
    for (codesP = codesString; *codesP != '\0'; codesP++) {
      if (*codesP == i || tolower(*codesP) == i) {
	bMatched = true;
	break;
      }
    }
    if (!bMatched)
      fprintf(stderr, "Ignored unknown character %c (seen %lu times)\n", i, counts[i]);
  }
    

  /* warnings about the counts */
  double fACGTUN = (counts['A'] + counts['C'] + counts['G'] + counts['T'] + counts['U'] + counts['N']
		    + counts['a'] + counts['c'] + counts['g'] + counts['t'] + counts['u'] + counts['n'])
    / (double)(totCount - counts['-'] - counts['.']);
  if (nCodes == 4 && fACGTUN < 0.9)
    fprintf(stderr, "WARNING! ONLY %.1f%% NUCLEOTIDE CHARACTERS -- IS THIS REALLY A NUCLEOTIDE ALIGNMENT?\n",
	    100.0 * fACGTUN);
  else if (nCodes == 20 && fACGTUN >= 0.9)
    fprintf(stderr, "WARNING! %.1f%% NUCLEOTIDE CHARACTERS -- IS THIS REALLY A PROTEIN ALIGNMENT?\n",
	    100.0 * fACGTUN);

  if(verbose>10) fprintf(stderr,"Made sequence profiles\n");
  for (iNode = NJ->nSeq; iNode < NJ->maxnodes; iNode++) 
    NJ->profiles[iNode] = NULL; /* not yet exists */

  NJ->outprofile = OutProfile(NJ->profiles, NJ->nSeq,
			      NJ->nPos, NJ->nConstraints,
			      NJ->distance_matrix);
  if(verbose>10) fprintf(stderr,"Made out-profile\n");

  NJ->totdiam = 0.0;

  NJ->diameter = (float *)mymalloc(sizeof(float)*NJ->maxnodes);
  for (iNode = 0; iNode < NJ->maxnodes; iNode++) NJ->diameter[iNode] = 0;

  NJ->varDiameter = (float *)mymalloc(sizeof(float)*NJ->maxnodes);
  for (iNode = 0; iNode < NJ->maxnodes; iNode++) NJ->varDiameter[iNode] = 0;

  NJ->selfdist = (float *)mymalloc(sizeof(float)*NJ->maxnodes);
  for (iNode = 0; iNode < NJ->maxnodes; iNode++) NJ->selfdist[iNode] = 0;

  NJ->selfweight = (float *)mymalloc(sizeof(float)*NJ->maxnodes);
  for (iNode = 0; iNode < NJ->nSeq; iNode++)
    NJ->selfweight[iNode] = NJ->nPos - NGaps(NJ,iNode);

  NJ->outDistances = (float *)mymalloc(sizeof(float)*NJ->maxnodes);
  NJ->nOutDistActive = (int *)mymalloc(sizeof(int)*NJ->maxnodes);
  for (iNode = 0; iNode < NJ->maxnodes; iNode++)
    NJ->nOutDistActive[iNode] = NJ->nSeq * 10; /* unreasonably high value */
  NJ->parent = NULL;		/* so SetOutDistance ignores it */
  for (iNode = 0; iNode < NJ->nSeq; iNode++)
    SetOutDistance(/*IN/UPDATE*/NJ, iNode, /*nActive*/NJ->nSeq);

  if (verbose>2) {
    for (iNode = 0; iNode < 4 && iNode < NJ->nSeq; iNode++)
      fprintf(stderr, "Node %d outdist %f\n", iNode, NJ->outDistances[iNode]);
  }

  NJ->parent = (int *)mymalloc(sizeof(int)*NJ->maxnodes);
  for (iNode = 0; iNode < NJ->maxnodes; iNode++) NJ->parent[iNode] = -1;

  NJ->branchlength = (float *)mymalloc(sizeof(float)*NJ->maxnodes); /* distance to parent */
  for (iNode = 0; iNode < NJ->maxnodes; iNode++) NJ->branchlength[iNode] = 0;

  NJ->support = (float *)mymalloc(sizeof(float)*NJ->maxnodes);
  for (iNode = 0; iNode < NJ->maxnodes; iNode++) NJ->support[iNode] = -1.0;

  NJ->child = (children_t*)mymalloc(sizeof(children_t)*NJ->maxnodes);
  for (iNode= 0; iNode < NJ->maxnode; iNode++) NJ->child[iNode].nChild = 0;

  return(NJ);
}

NJ_t *FreeNJ(NJ_t *NJ) {
  if (NJ==NULL)
    return(NJ);

  int i;
  for (i=0; i < NJ->maxnode; i++)
    NJ->profiles[i] = FreeProfile(NJ->profiles[i], NJ->nPos, NJ->nConstraints);
  NJ->profiles = myfree(NJ->profiles, sizeof(profile_t*) * NJ->maxnodes);
  NJ->outprofile = FreeProfile(NJ->outprofile, NJ->nPos, NJ->nConstraints);
  NJ->diameter = myfree(NJ->diameter, sizeof(float)*NJ->maxnodes);
  NJ->varDiameter = myfree(NJ->varDiameter, sizeof(float)*NJ->maxnodes);
  NJ->selfdist = myfree(NJ->selfdist, sizeof(float)*NJ->maxnodes);
  NJ->selfweight = myfree(NJ->selfweight, sizeof(float)*NJ->maxnodes);
  NJ->outDistances = myfree(NJ->outDistances, sizeof(float)*NJ->maxnodes);
  NJ->nOutDistActive = myfree(NJ->nOutDistActive, sizeof(int)*NJ->maxnodes);
  NJ->parent = myfree(NJ->parent, sizeof(int)*NJ->maxnodes);
  NJ->branchlength = myfree(NJ->branchlength, sizeof(float)*NJ->maxnodes);
  NJ->support = myfree(NJ->support, sizeof(float)*NJ->maxnodes);
  NJ->child = myfree(NJ->child, sizeof(children_t)*NJ->maxnodes);
  return(myfree(NJ, sizeof(NJ_t)));
}

void FastNJ(NJ_t *NJ) {
  int iNode;

  assert(NJ->nSeq >= 1);
  if (NJ->nSeq < 3) {
    NJ->root = NJ->maxnode++;
    NJ->child[NJ->root].nChild = NJ->nSeq;
    for (iNode = 0; iNode < NJ->nSeq; iNode++) {
      NJ->parent[iNode] = NJ->root;
      NJ->child[NJ->root].child[iNode] = iNode;
    }
    if (NJ->nSeq == 1) {
      NJ->branchlength[0] = 0;
    } else {
      assert (NJ->nSeq == 2);
      besthit_t hit;
      SeqDist(NJ->profiles[0]->codes,NJ->profiles[1]->codes,NJ->nPos,NJ->distance_matrix,/*OUT*/&hit);
      NJ->branchlength[0] = hit.dist/2.0;
      NJ->branchlength[1] = hit.dist/2.0;
    }
    return;
  }

  /* else 3 or more sequences */

  /* The visible set stores the best hit of each node */
  besthit_t *visible = NULL;
  besthit_t *besthitNew = NULL;	/* All hits of new node -- not used if doing top-hits */
  int *topvisible = NULL;	/* The top m candidates in the visible set */

  /* The top-hits lists, with the key parameter m = length of each top-hit list */
  besthit_t **tophits = NULL;	/* Up to top m hits for each node; i and j are -1 if past end of list */
  int *tophitAge = NULL;	/* #Joins since list was refreshed, 1 value per node */
  int m = 0;			/* length of each list */
  if (tophitsMult > 0) {
    m = (int)(0.5 + tophitsMult*sqrt(NJ->nSeq));
    if(m<4 || 2*m >= NJ->nSeq) {
      m=0;
      if(verbose>1) fprintf(stderr,"Too few leaves, turning off top-hits\n");
    } else {
      if(verbose>2) fprintf(stderr,"Top-hit-list size = %d of %d\n", m, NJ->nSeq);
    }
  }


  assert(!(slow && m>0));

  if (m>0) {
      tophits = (besthit_t**)mymalloc(sizeof(besthit_t*) * NJ->maxnodes);
      for(iNode=0; iNode < NJ->maxnodes; iNode++) tophits[iNode] = NULL;
      SetAllLeafTopHits(NJ, m, /*OUT*/tophits);
      tophitAge = (int*)mymalloc(sizeof(int) * NJ->maxnodes);
      for(iNode=0; iNode < NJ->maxnodes; iNode++) tophitAge[iNode] = 0;
#ifdef TRACK_MEMORY
      if (verbose>1) {
	struct mallinfo mi = mallinfo();
	fprintf(stderr, "Memory after SetAllLeafTopHits(): %.2f MB (%.1f byte/pos) useful %.2f\n",
		(mi.arena+mi.hblkhd)/1.0e6, (mi.arena+mi.hblkhd)/(double)(NJ->nSeq*(double)NJ->nPos),
		mi.uordblks/1.0e6);
      }
#endif
      topvisible = (int*)mymalloc(sizeof(int)*m);
      int i;
      for (i = 0; i < m; i++)
	topvisible[i] = -1;
  }

  /* Initialize visible set */
  if (!slow) {
    visible = (besthit_t*)mymalloc(sizeof(besthit_t)*NJ->maxnodes);
    if(m==0) besthitNew = (besthit_t*)mymalloc(sizeof(besthit_t)*NJ->maxnodes);
    for (iNode = 0; iNode < NJ->nSeq; iNode++) {
      if (m>0)
	visible[iNode] = tophits[iNode][GetBestFromTopHits(iNode, NJ, /*nActive*/NJ->nSeq, tophits[iNode], /*nTop*/m)];
      else
	SetBestHit(iNode, NJ, /*nActive*/NJ->nSeq, /*OUT*/&visible[iNode], /*OUT IGNORED*/NULL);
    }
  }

  /* Iterate over joins */
  int nActiveOutProfileReset = NJ->nSeq;
  int nActive;
  int topVisibleAge = 0;
  for (nActive = NJ->nSeq; nActive > 3; nActive--) {
    int nJoinsDone = NJ->nSeq - nActive;
    if (nJoinsDone > 0 && (nJoinsDone % 100) == 0)
      ProgressReport("Joined %6d of %6d", nJoinsDone, NJ->nSeq-3, 0, 0);
    
    besthit_t join; 		/* the join to do */
    if (slow) {
      ExhaustiveNJSearch(NJ,nActive,/*OUT*/&join);
    } else if (m>0) {
      TopHitNJSearch(/*IN/OUT*/NJ, nActive, m,
		     /*IN/OUT*/visible, /*IN/OUT*/topvisible, /*IN/OUT*/&topVisibleAge,
		     /*IN/OUT*/tophits, /*OUT*/&join);
    } else {
      FastNJSearch(NJ, nActive, /*IN/OUT*/visible, /*OUT*/&join);
    }

    if (verbose>2) {
      double penalty = constraintWeight
	* (double)JoinConstraintPenalty(NJ, join.i, join.j);
      if (penalty > 0.001) {
	fprintf(stderr, "Constraint violation during neighbor-joining %d %d into %d penalty %.3f\n",
		join.i, join.j, NJ->maxnode, penalty);
	int iC;
	for (iC = 0; iC < NJ->nConstraints; iC++) {
	  int local = JoinConstraintPenaltyPiece(NJ, join.i, join.j, iC);
	  if (local > 0)
	    fprintf(stderr, "Constraint %d piece %d %d/%d %d/%d %d/%d\n", iC, local,
		    NJ->profiles[join.i]->nOn[iC],
		    NJ->profiles[join.i]->nOff[iC],
		    NJ->profiles[join.j]->nOn[iC],
		    NJ->profiles[join.j]->nOff[iC],
		    NJ->outprofile->nOn[iC] - NJ->profiles[join.i]->nOn[iC] - NJ->profiles[join.j]->nOn[iC],
		    NJ->outprofile->nOff[iC] - NJ->profiles[join.i]->nOff[iC] - NJ->profiles[join.j]->nOff[iC]);
	}
      }
    }

    /* because of the stale out-distance heuristic, make sure that these are up-to-date */
    SetOutDistance(NJ, join.i, nActive);
    SetOutDistance(NJ, join.j, nActive);
    assert(NJ->nOutDistActive[join.i] == nActive);
    assert(NJ->nOutDistActive[join.j] == nActive);

    int newnode = NJ->maxnode++;
    NJ->parent[join.i] = newnode;
    NJ->parent[join.j] = newnode;
    NJ->child[newnode].nChild = 2;
    NJ->child[newnode].child[0] = join.i < join.j ? join.i : join.j;
    NJ->child[newnode].child[1] = join.i > join.j ? join.i : join.j;

    double rawIJ = join.dist + NJ->diameter[join.i] + NJ->diameter[join.j];
    double distIJ = join.dist;

    double deltaDist = (NJ->outDistances[join.i]-NJ->outDistances[join.j])/(double)(nActive-2);
    NJ->branchlength[join.i] = (distIJ + deltaDist)/2;
    NJ->branchlength[join.j] = (distIJ - deltaDist)/2;

    double bionjWeight = 0.5;	/* IJ = bionjWeight*I + (1-bionjWeight)*J */
    double varIJ = rawIJ - NJ->varDiameter[join.i] - NJ->varDiameter[join.j];

    if (bionj && join.weight > 0.01 && varIJ > 0.001) {
      /* Set bionjWeight according to the BIONJ formula, where
	 the variance matrix is approximated by

	 Vij = ProfileVar(i,j) - varDiameter(i) - varDiameter(j)
	 ProfileVar(i,j) = distance(i,j) = top(i,j)/weight(i,j)

	 (The node's distance diameter does not affect the variances.)

	 The BIONJ formula is equation 9 from Gascuel 1997:

	 bionjWeight = 1/2 + sum(k!=i,j) (Vjk - Vik) / ((nActive-2)*Vij)
	 sum(k!=i,j) (Vjk - Vik) = sum(k!=i,j) Vik - varDiameter(j) + varDiameter(i)
	 = sum(k!=i,j) ProfileVar(j,k) - sum(k!=i,j) ProfileVar(i,k) + (nActive-2)*(varDiameter(i)-varDiameter(j))

	 sum(k!=i,j) ProfileVar(i,k)
	 ~= (sum(k!=i,j) distance(i,k) * weight(i,k))/(mean(k!=i,j) weight(i,k))
	 ~= (N-2) * top(i, Out-i-j) / weight(i, Out-i-j)

	 weight(i, Out-i-j) = N*weight(i,Out) - weight(i,i) - weight(i,j)
	 top(i, Out-i-j) = N*top(i,Out) - top(i,i) - top(i,j)
      */
      besthit_t outI;
      besthit_t outJ;
      ProfileDist(NJ->profiles[join.i],NJ->outprofile,NJ->nPos,NJ->distance_matrix,/*OUT*/&outI);
      ProfileDist(NJ->profiles[join.j],NJ->outprofile,NJ->nPos,NJ->distance_matrix,/*OUT*/&outJ);
      outprofileOps += 2;

      double varIWeight = (nActive * outI.weight - NJ->selfweight[join.i] - join.weight);
      double varJWeight = (nActive * outJ.weight - NJ->selfweight[join.j] - join.weight);

      double varITop = outI.dist * outI.weight * nActive
	- NJ->selfdist[join.i] * NJ->selfweight[join.i] - rawIJ * join.weight;
      double varJTop = outJ.dist * outJ.weight * nActive
	- NJ->selfdist[join.j] * NJ->selfweight[join.j] - rawIJ * join.weight;

      double deltaProfileVarOut = (nActive-2) * (varJTop/varJWeight - varITop/varIWeight);
      double deltaVarDiam = (nActive-2)*(NJ->varDiameter[join.i] - NJ->varDiameter[join.j]);
      if (varJWeight > 0.01 && varIWeight > 0.01)
	bionjWeight = 0.5 + (deltaProfileVarOut+deltaVarDiam)/(2*(nActive-2)*varIJ);
      if(bionjWeight<0) bionjWeight=0;
      if(bionjWeight>1) bionjWeight=1;
      if (verbose>2) fprintf(stderr,"dVarO %f dVarDiam %f varIJ %f from dist %f weight %f (pos %d) bionjWeight %f %f\n",
			     deltaProfileVarOut, deltaVarDiam,
			     varIJ, join.dist, join.weight, NJ->nPos,
			     bionjWeight, 1-bionjWeight);
      if (verbose>3 && (newnode%5) == 0) {
	/* Compare weight estimated from outprofiles from weight made by summing over other nodes */
	double deltaProfileVarTot = 0;
	for (iNode = 0; iNode < newnode; iNode++) {
	  if (NJ->parent[iNode] < 0) { /* excludes join.i, join.j */
	    besthit_t di, dj;
	    ProfileDist(NJ->profiles[join.i],NJ->profiles[iNode],NJ->nPos,NJ->distance_matrix,/*OUT*/&di);
	    ProfileDist(NJ->profiles[join.j],NJ->profiles[iNode],NJ->nPos,NJ->distance_matrix,/*OUT*/&dj);
	    deltaProfileVarTot += dj.dist - di.dist;
	  }
	}
	double lambdaTot = 0.5 + (deltaProfileVarTot+deltaVarDiam)/(2*(nActive-2)*varIJ);
	if (lambdaTot < 0) lambdaTot = 0;
	if (lambdaTot > 1) lambdaTot = 1;
	if (fabs(bionjWeight-lambdaTot) > 0.01 || verbose > 4)
	  fprintf(stderr, "deltaProfileVar actual %.6f estimated %.6f lambda actual %.3f estimated %.3f\n",
		  deltaProfileVarTot,deltaProfileVarOut,lambdaTot,bionjWeight);
      }
    }
    if (verbose>1) fprintf(stderr, "Join\t%d\t%d\t%.6f\tlambda\t%.6f\tselfw\t%.3f\t%.3f\tnew\t%d\n",
			   join.i < join.j ? join.i : join.j,
			   join.i < join.j ? join.j : join.i,
			   join.criterion, bionjWeight,
			   NJ->selfweight[join.i < join.j ? join.i : join.j],
			   NJ->selfweight[join.i < join.j ? join.j : join.i],
			   newnode);
    
    NJ->diameter[newnode] = bionjWeight * (NJ->branchlength[join.i] + NJ->diameter[join.i])
      + (1-bionjWeight) * (NJ->branchlength[join.j] + NJ->diameter[join.j]);
    NJ->varDiameter[newnode] = bionjWeight * NJ->varDiameter[join.i]
      + (1-bionjWeight) * NJ->varDiameter[join.j]
      + bionjWeight * (1-bionjWeight) * varIJ;

    NJ->profiles[newnode] = AverageProfile(NJ->profiles[join.i],NJ->profiles[join.j],
					   NJ->nPos, NJ->nConstraints,
					   NJ->distance_matrix,
					   bionj ? bionjWeight : /*noweight*/-1.0);

    /* Update out-distances and total diameters */
    int changedActiveOutProfile = nActiveOutProfileReset - (nActive-1);
    if (changedActiveOutProfile >= nResetOutProfile
	&& changedActiveOutProfile >= fResetOutProfile * nActiveOutProfileReset) {
      /* Recompute the outprofile from scratch to avoid roundoff error */
      profile_t **activeProfiles = (profile_t**)mymalloc(sizeof(profile_t*)*(nActive-1));
      int nSaved = 0;
      NJ->totdiam = 0;
      for (iNode=0;iNode<NJ->maxnode;iNode++) {
	if (NJ->parent[iNode]<0) {
	  assert(nSaved < nActive-1);
	  activeProfiles[nSaved++] = NJ->profiles[iNode];
	  NJ->totdiam += NJ->diameter[iNode];
	}
      }
      assert(nSaved==nActive-1);
      FreeProfile(NJ->outprofile, NJ->nPos, NJ->nConstraints);
      if(verbose>1) fprintf(stderr,"Recomputing outprofile %d %d\n",nActiveOutProfileReset,nActive-1);
      NJ->outprofile = OutProfile(activeProfiles, nSaved,
				  NJ->nPos, NJ->nConstraints,
				  NJ->distance_matrix);
      activeProfiles = myfree(activeProfiles, sizeof(profile_t*)*(nActive-1));
      nActiveOutProfileReset = nActive-1;
    } else {
      UpdateOutProfile(/*OUT*/NJ->outprofile,
		       NJ->profiles[join.i], NJ->profiles[join.j], NJ->profiles[newnode],
		       nActive,
		       NJ->nPos, NJ->nConstraints,
		       NJ->distance_matrix);
      NJ->totdiam += NJ->diameter[newnode] - NJ->diameter[join.i] - NJ->diameter[join.j];
    }

    /* Store self-dist for use in other computations */
    besthit_t selfdist;
    ProfileDist(NJ->profiles[newnode],NJ->profiles[newnode],NJ->nPos,NJ->distance_matrix,/*OUT*/&selfdist);
    NJ->selfdist[newnode] = selfdist.dist;
    NJ->selfweight[newnode] = selfdist.weight;

    /* Find the best hit of the joined node IJ */
    if (m>0) {
      TopHitJoin(/*IN/OUT*/NJ, newnode, nActive-1, m,
		 /*IN/OUT*/tophits, /*IN/OUT*/tophitAge,
		 /*IN/OUT*/visible, /*IN/OUT*/topvisible);
    } else {
      /* Not using top-hits, so we update all out-distances */
      for (iNode = 0; iNode < NJ->maxnode; iNode++) {
	if (NJ->parent[iNode] < 0) {
	  /* True nActive is now nActive-1 */
	  SetOutDistance(/*IN/UPDATE*/NJ, iNode, nActive-1);
	}
      }
    
      if(visible != NULL) {
	SetBestHit(newnode, NJ, nActive-1, /*OUT*/&visible[newnode], /*OUT OPTIONAL*/besthitNew);
	if (verbose>2)
	  fprintf(stderr,"Visible %d %d %f %f\n",
		  visible[newnode].i, visible[newnode].j,
		  visible[newnode].dist, visible[newnode].criterion);
	if (besthitNew != NULL) {
	  /* Use distances to new node to update visible set entries that are non-optimal */
	  for (iNode = 0; iNode < NJ->maxnode; iNode++) {
	    if (NJ->parent[iNode] >= 0 || iNode == newnode)
	      continue;
	    int iOldVisible = visible[iNode].j;
	    assert(iOldVisible>=0);
	    assert(visible[iNode].i == iNode);
	      
	    /* Update the criterion; use nActive-1 because haven't decremented nActive yet */
	    if (NJ->parent[iOldVisible] < 0)
	      SetCriterion(/*IN/OUT*/NJ, nActive-1, &visible[iNode]);
	    
	    if (NJ->parent[iOldVisible] >= 0
		|| besthitNew[iNode].criterion < visible[iNode].criterion) {
	      if(verbose>3) fprintf(stderr,"Visible %d reset from %d to %d (%f vs. %f)\n",
				     iNode, iOldVisible, 
				     newnode, visible[iNode].criterion, besthitNew[iNode].criterion);
	      if(NJ->parent[iOldVisible] < 0) nVisibleReset++;
	      visible[iNode].j = newnode;
	      visible[iNode].dist = besthitNew[iNode].dist;
	      visible[iNode].criterion = besthitNew[iNode].criterion;
	    }
	  } /* end loop over all nodes */
	} /* end if recording all hits of new node */
      } /* end if keeping a visible set */
    } /* end else (m==0) */
  } /* end loop over nActive */

#ifdef TRACK_MEMORY
  if (verbose>1) {
    struct mallinfo mi = mallinfo();
    fprintf(stderr, "Memory @ end of FastNJ(): %.2f MB (%.1f byte/pos) useful %.2f expected %.2f\n",
	    (mi.arena+mi.hblkhd)/1.0e6, (mi.arena+mi.hblkhd)/(double)(NJ->nSeq*(double)NJ->nPos),
	    mi.uordblks/1.0e6, mymallocUsed/1e6);
  }
#endif

  /* We no longer need the tophits, visible set, etc. */
  if (visible != NULL) visible = myfree(visible,sizeof(besthit_t)*NJ->maxnodes);
  if (besthitNew != NULL) besthitNew = myfree(besthitNew,sizeof(besthit_t)*NJ->maxnodes);
  if (tophits != NULL) {
    for (iNode = 0; iNode < NJ->maxnode; iNode++) {
      if (tophits[iNode] != NULL)
	tophits[iNode] = myfree(tophits[iNode],sizeof(besthit_t)*m);
    }
    tophits = myfree(tophits, sizeof(besthit_t*)*NJ->maxnodes);
    tophitAge = myfree(tophitAge,sizeof(int)*NJ->maxnodes);
  }
  topvisible = myfree(topvisible, sizeof(int)*m);

  /* Add a root for the 3 remaining nodes */
  int top[3];
  int nTop = 0;
  for (iNode = 0; iNode < NJ->maxnode; iNode++) {
    if (NJ->parent[iNode] < 0) {
      assert(nTop <= 2);
      top[nTop++] = iNode;
    }
  }
  assert(nTop==3);
  
  NJ->root = NJ->maxnode++;
  NJ->child[NJ->root].nChild = 3;
  for (nTop = 0; nTop < 3; nTop++) {
    NJ->parent[top[nTop]] = NJ->root;
    NJ->child[NJ->root].child[nTop] = top[nTop];
  }

  besthit_t dist01, dist02, dist12;
  ProfileDist(NJ->profiles[top[0]], NJ->profiles[top[1]], NJ->nPos, NJ->distance_matrix, /*OUT*/&dist01);
  ProfileDist(NJ->profiles[top[0]], NJ->profiles[top[2]], NJ->nPos, NJ->distance_matrix, /*OUT*/&dist02);
  ProfileDist(NJ->profiles[top[1]], NJ->profiles[top[2]], NJ->nPos, NJ->distance_matrix, /*OUT*/&dist12);

  double d01 = dist01.dist - NJ->diameter[top[0]] - NJ->diameter[top[1]];
  double d02 = dist02.dist - NJ->diameter[top[0]] - NJ->diameter[top[2]];
  double d12 = dist12.dist - NJ->diameter[top[1]] - NJ->diameter[top[2]];
  NJ->branchlength[top[0]] = (d01 + d02 - d12)/2;
  NJ->branchlength[top[1]] = (d01 + d12 - d02)/2;
  NJ->branchlength[top[2]] = (d02 + d12 - d01)/2;

  /* Check how accurate the outprofile is */
  if (verbose>1) {
    profile_t *p[3] = {NJ->profiles[top[0]], NJ->profiles[top[1]], NJ->profiles[top[2]]};
    profile_t *out = OutProfile(p, 3, NJ->nPos, NJ->nConstraints, NJ->distance_matrix);
    int i;
    double freqerror = 0;
    double weighterror = 0;
    for (i=0;i<NJ->nPos;i++) {
      weighterror += fabs(out->weights[i] - NJ->outprofile->weights[i]);
      int k;
      for(k=0;k<nCodes;k++)
	freqerror += fabs(out->vectors[nCodes*i+k] - NJ->outprofile->vectors[nCodes*i+k]);
    }
    fprintf(stderr,"Roundoff error in outprofile@end: WeightError %f FreqError %f\n", weighterror, freqerror);
    FreeProfile(out, NJ->nPos, NJ->nConstraints);
  }
  return;
}

void ExhaustiveNJSearch(NJ_t *NJ, int nActive, /*OUT*/besthit_t *join) {
  join->i = -1;
  join->j = -1;
  join->weight = 0;
  join->dist = 1e20;
  join->criterion = 1e20;
  double bestCriterion = 1e20;

  int i, j;
  for (i = 0; i < NJ->maxnode-1; i++) {
    if (NJ->parent[i] < 0) {
      for (j = i+1; j < NJ->maxnode; j++) {
	if (NJ->parent[j] < 0) {
	  besthit_t hit;
	  hit.i = i;
	  hit.j = j;
	  SetDistCriterion(NJ, nActive, /*IN/OUT*/&hit);
	  if (hit.criterion < bestCriterion) {
	    *join = hit;
	    bestCriterion = hit.criterion;
	  }
	}
      }
    }
  }
  assert (join->i >= 0 && join->j >= 0);
}

void FastNJSearch(NJ_t *NJ, int nActive, /*IN/OUT*/besthit_t *besthits, /*OUT*/besthit_t *join) {
  join->i = -1;
  join->j = -1;
  join->dist = 1e20;
  join->weight = 0;
  join->criterion = 1e20;
  int iNode;
  for (iNode = 0; iNode < NJ->maxnode; iNode++) {
    int jNode = besthits[iNode].j;
    if (NJ->parent[iNode] < 0 && NJ->parent[jNode] < 0) { /* both i and j still active */
      /* recompute criterion to reflect the current out-distances */
      SetCriterion(NJ, nActive, /*IN/OUT*/&besthits[iNode]);
      if (besthits[iNode].criterion < join->criterion)
	*join = besthits[iNode];      
    }
  }

  if(!fastest) {
    int changed;
    do {
      changed = 0;
      assert(join->i >= 0 && join->j >= 0);
      SetBestHit(join->i, NJ, nActive, /*OUT*/&besthits[join->i], /*OUT IGNORED*/NULL);
      if (besthits[join->i].j != join->j) {
	changed = 1;
	if (verbose>2)
	  fprintf(stderr,"BetterI\t%d\t%d\t%d\t%d\t%f\t%f\n",
		  join->i,join->j,besthits[join->i].i,besthits[join->i].j,
		  join->criterion,besthits[join->i].criterion);
      }
      
      // Save the best hit either way, because the out-distance has probably changed
      // since we started the computation.
      join->j = besthits[join->i].j;
      join->weight = besthits[join->i].weight;
      join->dist = besthits[join->i].dist;
      join->criterion = besthits[join->i].criterion;
      
      SetBestHit(join->j, NJ, nActive, /*OUT*/&besthits[join->j], /*OUT IGNORE*/NULL);
      if (besthits[join->j].j != join->i) {
	changed = 1;
	if (verbose>2)
	  fprintf(stderr,"BetterJ\t%d\t%d\t%d\t%d\t%f\t%f\n",
		  join->i,join->j,besthits[join->j].i,besthits[join->j].j,
		  join->criterion,besthits[join->j].criterion);
	join->i = besthits[join->j].j;
	join->weight = besthits[join->j].weight;
	join->dist = besthits[join->j].dist;
	join->criterion = besthits[join->j].criterion;
      }
      if(changed) nBetter++;
    } while(changed);
  }
}

/* A token is one of ():;, or an alphanumeric string without whitespace
   Any whitespace between tokens is ignored */
char *ReadTreeToken(FILE *fp) {
  static char buf[BUFFER_SIZE];
  int len = 0;
  int c;
  for (c = fgetc(fp); c != EOF; c = fgetc(fp)) {
    if (c == '(' || c == ')' || c == ':' || c == ';' || c == ',') {
      /* standalone token */
      if (len == 0) {
	buf[len++] = c;
	buf[len] = '\0';
	return(buf);
      } else {
	ungetc(c, fp);
	buf[len] = '\0';
	return(buf);
      }
    } else if (isspace(c)) {
      if (len > 0) {
	buf[len] = '\0';
	return(buf);
      }
      /* else ignore whitespace at beginning of token */
    } else {
      /* not whitespace or standalone token */
      buf[len++] = c;
      if (len >= BUFFER_SIZE) {
	buf[BUFFER_SIZE-1] = '\0';
	fprintf(stderr, "Token too long in tree file, token begins with\n%s\n", buf);
	exit(1);
      }
    }
  }
  if (len > 0) {
    /* return the token we have so far */
    buf[len] = '\0';
    return(buf);
  }
  /* else */
  return(NULL);
}

void ReadTreeError(char *err, char *token) {
  fprintf(stderr, "Tree parse error: unexpected token '%s' -- %s\n",
	  token == NULL ? "(End of file)" : token,
	  err);
  exit(1);
}

void ReadTreeAddChild(int parent, int child, /*IN/OUT*/int *parents, /*IN/OUT*/children_t *children) {
  assert(parent >= 0);
  assert(child >= 0);
  assert(parents[child] < 0);
  assert(children[parent].nChild < 3);
  parents[child] = parent;
  children[parent].child[children[parent].nChild++] = child;
}

void ReadTreeMaybeAddLeaf(int parent, char *name,
			  hashstrings_t *hashnames, uniquify_t *unique,
			  /*IN/OUT*/int *parents, /*IN/OUT*/children_t *children) {
  hashiterator_t hi = FindMatch(hashnames,name);
  if (HashCount(hashnames,hi) != 1)
    ReadTreeError("not recognized as a sequence name", name);

  int iSeqNonunique = HashFirst(hashnames,hi);
  assert(iSeqNonunique >= 0 && iSeqNonunique < unique->nSeq);
  int iSeqUnique = unique->alnToUniq[iSeqNonunique];
  assert(iSeqUnique >= 0 && iSeqUnique < unique->nUnique);
  /* Either record this leaves' parent (if it is -1) or ignore this leaf (if already seen) */
  if (parents[iSeqUnique] < 0) {
    ReadTreeAddChild(parent, iSeqUnique, /*IN/OUT*/parents, /*IN/OUT*/children);
    if(verbose > 5)
      fprintf(stderr, "Found leaf uniq%d name %s child of %d\n", iSeqUnique, name, parent);
  } else {
    if (verbose > 5)
      fprintf(stderr, "Skipped redundant leaf uniq%d name %s\n", iSeqUnique, name);
  }
}

void ReadTreeRemove(/*IN/OUT*/int *parents, /*IN/OUT*/children_t *children, int node) {
  if(verbose > 5)
    fprintf(stderr,"Removing node %d parent %d\n", node, parents[node]);
  assert(parents[node] >= 0);
  int parent = parents[node];
  parents[node] = -1;
  children_t *pc = &children[parent];
  int oldn;
  for (oldn = 0; oldn < pc->nChild; oldn++) {
    if (pc->child[oldn] == node)
      break;
  }
  assert(oldn < pc->nChild);

  /* move successor nodes back in child list and shorten list */
  int i;
  for (i = oldn; i < pc->nChild-1; i++)
    pc->child[i] = pc->child[i+1];
  pc->nChild--;

  /* add its children to parent's child list */
  children_t *nc = &children[node];
  if (nc->nChild > 0) {
    assert(nc->nChild<=2);
    assert(pc->nChild < 3);
    assert(pc->nChild + nc->nChild <= 3);
    int j;
    for (j = 0; j < nc->nChild; j++) {
      if(verbose > 5)
	fprintf(stderr,"Repointing parent %d to child %d\n", parent, nc->child[j]);
      pc->child[pc->nChild++] = nc->child[j];
      parents[nc->child[j]] = parent;
    }
    nc->nChild = 0;
  }
}  

void ReadTree(/*IN/OUT*/NJ_t *NJ,
	      /*IN*/uniquify_t *unique,
	      /*IN*/hashstrings_t *hashnames,
	      /*READ*/FILE *fpInTree) {
  assert(NJ->nSeq == unique->nUnique);
  /* First, do a preliminary parse of the tree to with non-unique leaves ignored
     We need to store this separately from NJ because it may have too many internal nodes
     (matching sequences show up once in the NJ but could be in multiple places in the tree)
     Will use iUnique as the index of nodes, as in the NJ structure
  */
  int maxnodes = unique->nSeq*2;
  int maxnode = unique->nSeq;
  int *parent = (int*)mymalloc(sizeof(int)*maxnodes);
  children_t *children = (children_t *)mymalloc(sizeof(children_t)*maxnodes);
  int root = maxnode++;
  int i;
  for (i = 0; i < maxnodes; i++) {
    parent[i] = -1;
    children[i].nChild = 0;
  }

  /* The stack is the current path to the root, with the root at the first (top) position */
  int stack_size = 1;
  int *stack = (int*)mymalloc(sizeof(int)*maxnodes);
  stack[0] = root;
  int nDown = 0;
  int nUp = 0;

  char *token;
  token = ReadTreeToken(fpInTree);
  if (token == NULL || *token != '(')
    ReadTreeError("No '(' at start", token);
  /* nDown is still 0 because we have created the root */

  while ((token = ReadTreeToken(fpInTree)) != NULL) {
    if (nDown > 0) {		/* In a stream of parentheses */
      if (*token == '(')
	nDown++;
      else if (*token == ',' || *token == ';' || *token == ':' || *token == ')')
	ReadTreeError("while reading parentheses", token);
      else {
	/* Add intermediate nodes if nDown was > 1 (for nDown=1, the only new node is the leaf) */
	while (nDown-- > 0) {
	  int new = maxnode++;
	  assert(new < maxnodes);
	  ReadTreeAddChild(stack[stack_size-1], new, /*IN/OUT*/parent, /*IN/OUT*/children);
	  if(verbose > 5)
	    fprintf(stderr, "Added internal child %d of %d, stack size increase to %d\n",
		    new, stack[stack_size-1],stack_size+1);
	  stack[stack_size++] = new;
	  assert(stack_size < maxnodes);
	}
	ReadTreeMaybeAddLeaf(stack[stack_size-1], token,
			     hashnames, unique,
			     /*IN/OUT*/parent, /*IN/OUT*/children);
      }
    } else if (nUp > 0) {
      if (*token == ';') {	/* end the tree? */
	if (nUp != stack_size)
	  ReadTreeError("unbalanced parentheses", token);
	else
	  break;
      } else if (*token == ')')
	nUp++;
      else if (*token == '(')
	ReadTreeError("unexpected '(' after ')'", token);
      else if (*token == ':') {
	token = ReadTreeToken(fpInTree);
	/* Read the branch length and ignore it */
	if (token == NULL || (*token != '-' && !isdigit(*token)))
	  ReadTreeError("not recognized as a branch length", token);
      } else if (*token == ',') {
	/* Go back up the stack the correct #times */
	while (nUp-- > 0) {
	  stack_size--;
	  if(verbose > 5)
	    fprintf(stderr, "Up to nUp=%d stack size %d at %d\n",
		    nUp, stack_size, stack[stack_size-1]);
	  if (stack_size <= 0)
	    ReadTreeError("too many ')'", token);
	}
	nUp = 0;
      } else if (*token == '-' || isdigit(*token))
	; 			/* ignore bootstrap value */
      else
	fprintf(stderr, "Warning while parsing tree: non-numeric label %s for internal node\n",
		token);
    } else if (*token == '(') {
      nDown = 1;
    } else if (*token == ')') {
      nUp = 1;
    } else if (*token == ':') {
      token = ReadTreeToken(fpInTree);
      if (token == NULL || (*token != '-' && !isdigit(*token)))
	ReadTreeError("not recognized as a branch length", token);
    } else if (*token == ',') {
      ;				/* do nothing */
    } else if (*token == ';')
      ReadTreeError("unexpected token", token);
    else
      ReadTreeMaybeAddLeaf(stack[stack_size-1], token,
			   hashnames, unique,
			   /*IN/OUT*/parent, /*IN/OUT*/children);
  }

  /* Verify that all sequences were seen */
  for (i = 0; i < unique->nUnique; i++) {
    if (parent[i] < 0) {
      fprintf(stderr, "Alignment sequence %d (unique %d) absent from input tree\n", unique->uniqueFirst[i], i);
      exit(1);
    }
  }

  /* Simplify the tree -- remove all internal nodes with < 2 children
     Keep trying until no nodes get removed
  */
  int nRemoved;
  do {
    nRemoved = 0;
    /* Here stack is the list of nodes we haven't visited yet while doing
       a tree traversal */
    stack_size = 1;
    stack[0] = root;
    while (stack_size > 0) {
      int node = stack[--stack_size];
      if (node >= unique->nUnique) { /* internal node */
	if (children[node].nChild <= 1) {
	  if (node != root) {
	    ReadTreeRemove(/*IN/OUT*/parent,/*IN/OUT*/children,node);
	    nRemoved++;
	  } else if (node == root && children[node].nChild == 1) {
	    int newroot = children[node].child[0];
	    parent[newroot] = -1;
	    children[root].nChild = 0;
	    nRemoved++;
	    if(verbose > 5)
	      fprintf(stderr,"Changed root from %d to %d\n",root,newroot);
	    root = newroot;
	    stack[stack_size++] = newroot;
	  }
	} else {
	  int j;
	  for (j = 0; j < children[node].nChild; j++) {
	    assert(stack_size < maxnodes);
	    stack[stack_size++] = children[node].child[j];
	    if(verbose > 5)
	      fprintf(stderr,"Added %d to stack\n", stack[stack_size-1]);
	  }
	}
      }
    }
  } while (nRemoved > 0);

  /* Simplify the root node to 3 children if it has 2 */
  if (children[root].nChild == 2) {
    for (i = 0; i < 2; i++) {
      int child = children[root].child[i];
      assert(child >= 0 && child < maxnodes);
      if (children[child].nChild == 2) {
	ReadTreeRemove(parent,children,child); /* replace root -> child -> A,B with root->A,B */
	break;
      }
    }
  }

  for (i = 0; i < maxnodes; i++)
    if(verbose > 5)
      fprintf(stderr,"Simplfied node %d has parent %d nchild %d\n",
	      i, parent[i], children[i].nChild);

  /* Map the remaining internal nodes to NJ nodes */
  int *map = (int*)mymalloc(sizeof(int)*maxnodes);
  for (i = 0; i < unique->nUnique; i++)
    map[i] = i;
  for (i = unique->nUnique; i < maxnodes; i++)
    map[i] = -1;
  stack_size = 1;
  stack[0] = root;
  while (stack_size > 0) {
    int node = stack[--stack_size];
    if (node >= unique->nUnique) { /* internal node */
      assert(node == root || children[node].nChild > 1);
      map[node] =  NJ->maxnode++;
      for (i = 0; i < children[node].nChild; i++) {
	assert(stack_size < maxnodes);
	stack[stack_size++] = children[node].child[i];
      }
    }
  }
  for (i = 0; i < maxnodes; i++)
    if(verbose > 5)
      fprintf(stderr,"Map %d to %d (parent %d nchild %d)\n",
	      i, map[i], parent[i], children[i].nChild);

  /* Set NJ->parent, NJ->children, NJ->root */
  NJ->root = map[root];
  int node;
  for (node = 0; node < maxnodes; node++) {
    int njnode = map[node];
    if (njnode >= 0) {
      NJ->child[njnode].nChild = children[node].nChild;
      for (i = 0; i < children[node].nChild; i++) {
	assert(children[node].child[i] >= 0 && children[node].child[i] < maxnodes);
	NJ->child[njnode].child[i] = map[children[node].child[i]];
      }
      if (parent[node] >= 0)
	NJ->parent[njnode] = map[parent[node]];
    }
  }

  /* Make sure that parent/child relationships match */
  for (i = 0; i < NJ->maxnode; i++) {
    children_t *c = &NJ->child[i];
    int j;
    for (j = 0; j < c->nChild;j++)
      assert(c->child[j] >= 0 && c->child[j] < NJ->maxnode && NJ->parent[c->child[j]] == i);
  }
  assert(NJ->parent[NJ->root] < 0);

  map = myfree(map,sizeof(int)*maxnodes);
  stack = myfree(stack,sizeof(int)*maxnodes);
  children = myfree(children,sizeof(children_t)*maxnodes);
  parent = myfree(parent,sizeof(int)*maxnodes);

  /* Compute profiles as balanced -- the NNI stage will recompute these
     profiles anyway
  */
  traversal_t traversal = InitTraversal(NJ);
  node = NJ->root;
  while((node = TraversePostorder(node, NJ, /*IN/OUT*/traversal)) >= 0) {
    if (node >= NJ->nSeq && node != NJ->root)
      SetProfile(/*IN/OUT*/NJ, node, /*noweight*/-1.0);
  }
  traversal = FreeTraversal(traversal,NJ);
}

/* Print topology using node indices as node names */
void PrintNJInternal(FILE *fp, NJ_t *NJ) {
  if (NJ->nSeq < 4) {
    return;
  }
  typedef struct { int node; int end; } stack_t;
  stack_t *stack = (stack_t *)mymalloc(sizeof(stack_t)*NJ->maxnodes);
  int stackSize = 1;
  stack[0].node = NJ->root;
  stack[0].end = 0;

  while(stackSize>0) {
    stack_t *last = &stack[stackSize-1];
    stackSize--;
    /* Save last, as we are about to overwrite it */
    int node = last->node;
    int end = last->end;

    if (node < NJ->nSeq) {
      if (NJ->child[NJ->parent[node]].child[0] != node) fputs(",",fp);
      fprintf(fp, "%d", node);
    } else if (end) {
      fprintf(fp, ")%d", node);
    } else {
            if (node != NJ->root && NJ->child[NJ->parent[node]].child[0] != node) fprintf(fp, ",");
      fprintf(fp, "(");
      stackSize++;
      stack[stackSize-1].node = node;
      stack[stackSize-1].end = 1;
      children_t *c = &NJ->child[node];
      // put children on in reverse order because we use the last one first
      int i;
      for (i = c->nChild-1; i >=0; i--) {
	stackSize++;
	stack[stackSize-1].node = c->child[i];
	stack[stackSize-1].end = 0;
      }
    }
  }
  fprintf(fp, ";\n");
  stack = myfree(stack, sizeof(stack_t)*NJ->maxnodes);
}

void PrintNJ(FILE *fp, NJ_t *NJ, char **names, uniquify_t *unique) {
  /* And print the tree: depth first search
   * The stack contains
   * list of remaining children with their depth
   * parent node, with a flag of -1 so I know to print right-paren
   */
  if (NJ->nSeq==1 && unique->alnNext[unique->uniqueFirst[0]] >= 0) {
    /* Special case -- otherwise we end up with double parens */
    int first = unique->uniqueFirst[0];
    assert(first >= 0 && first < unique->nSeq);
    fprintf(fp,"(%s:0.0",names[first]);
    int iName = unique->alnNext[first];
    while (iName >= 0) {
      assert(iName < unique->nSeq);
      fprintf(fp,",%s:0.0",names[iName]);
      iName = unique->alnNext[iName];
    }
    fprintf(fp,");\n");
    return;
  }

  typedef struct { int node; int end; } stack_t;
  stack_t *stack = (stack_t *)mymalloc(sizeof(stack_t)*NJ->maxnodes);
  int stackSize = 1;
  stack[0].node = NJ->root;
  stack[0].end = 0;

  while(stackSize>0) {
    stack_t *last = &stack[stackSize-1];
    stackSize--;
    /* Save last, as we are about to overwrite it */
    int node = last->node;
    int end = last->end;

    if (node < NJ->nSeq) {
      if (NJ->child[NJ->parent[node]].child[0] != node) fputs(",",fp);
      int first = unique->uniqueFirst[node];
      assert(first >= 0 && first < unique->nSeq);
      /* Print the name, or the subtree of duplicate names */
      if (unique->alnNext[first] == -1) {
	fprintf(fp, names[first]);
      } else {
	fprintf(fp,"(%s:0.0",names[first]);
	int iName = unique->alnNext[first];
	while (iName >= 0) {
	  assert(iName < unique->nSeq);
	  fprintf(fp,",%s:0.0",names[iName]);
	  iName = unique->alnNext[iName];
	}
	fprintf(fp,")");
      }
      /* Print the branch length */
      fprintf(fp, ":%.5f", NJ->branchlength[node]);
    } else if (end) {
      if (node == NJ->root)
	fprintf(fp, ")");
      else if (nBootstrap > 0)
	fprintf(fp, ")%.3f:%.5f", NJ->support[node], NJ->branchlength[node]);
      else
	fprintf(fp, "):%.5f", NJ->branchlength[node]);
    } else {
      if (node != NJ->root && NJ->child[NJ->parent[node]].child[0] != node) fprintf(fp, ",");
      fprintf(fp, "(");
      stackSize++;
      stack[stackSize-1].node = node;
      stack[stackSize-1].end = 1;
      children_t *c = &NJ->child[node];
      // put children on in reverse order because we use the last one first
      int i;
      for (i = c->nChild-1; i >=0; i--) {
	stackSize++;
	stack[stackSize-1].node = c->child[i];
	stack[stackSize-1].end = 0;
      }
    }
  }
  fprintf(fp, ";\n");
  stack = myfree(stack, sizeof(stack_t)*NJ->maxnodes);
}

alignment_t *ReadAlignment(/*IN*/FILE *fp) {
  int nSeq = 0;
  int nPos = 0;
  char **names = NULL;
  char **seqs = NULL;
  char buf[BUFFER_SIZE] = "";
  if (fgets(buf,sizeof(buf),fp) == NULL) {
    fprintf(stderr, "Error reading header line\n");
    exit(1);
  }
  int nSaved = 100;
  if (buf[0] == '>') {
    /* FASTA, truncate names at any of these */
    char *nameStop = "(),: \t\r\n";
    char *seqSkip = " \t\r\n";
    seqs = (char**)mymalloc(sizeof(char*) * nSaved);
    names = (char**)mymalloc(sizeof(char*) * nSaved);

    do {
      /* loop over lines */
      if (buf[0] == '>') {
	/* truncate the name */
	char *p, *q;
	for (p = buf+1; *p != '\0'; p++) {
	  for (q = nameStop; *q != '\0'; q++) {
	    if (*p == *q) {
	      *p = '\0';
	      break;
	    }
	  }
	  if (*p == '\0') break;
	}

	/* allocate space for another sequence */
	nSeq++;
	if (nSeq > nSaved) {
	  int nNewSaved = nSaved*2;
	  seqs = myrealloc(seqs,sizeof(char*)*nSaved,sizeof(char*)*nNewSaved);
	  names = myrealloc(names,sizeof(char*)*nSaved,sizeof(char*)*nNewSaved);
	  nSaved = nNewSaved;
	}
	names[nSeq-1] = (char*)mymemdup(buf+1,strlen(buf));
	seqs[nSeq-1] = NULL;
      } else {
	/* count non-space characters and append to sequence */
	int nKeep = 0;
	char *p, *q;
	for (p=buf; *p != '\0'; p++) {
	  for (q=seqSkip; *q != '\0'; q++) {
	    if (*p == *q)
	      break;
	  }
	  if (*p != *q)
	    nKeep++;
	}
	int nOld = (seqs[nSeq-1] == NULL) ? 0 : strlen(seqs[nSeq-1]);
	seqs[nSeq-1] = (char*)myrealloc(seqs[nSeq-1], nOld, nOld+nKeep+1);
	if (nOld+nKeep > nPos)
	  nPos = nOld + nKeep;
	char *out = seqs[nSeq-1] + nOld;
	for (p=buf; *p != '\0'; p++) {
	  for (q=seqSkip; *q != '\0'; q++) {
	    if (*p == *q)
	      break;
	  }
	  if (*p != *q) {
	    *out = *p;
	    out++;
	  }
	}
	assert(out-seqs[nSeq-1] == nKeep + nOld);
	*out = '\0';
      }
    } while(fgets(buf,sizeof(buf),fp) != NULL);

    if (seqs[nSeq-1] == NULL) {
      fprintf(stderr, "No sequence data for last entry %s\n",names[nSeq-1]);
      exit(1);
    }
    names = myrealloc(names,sizeof(char*)*nSaved,sizeof(char*)*nSeq);
    seqs = myrealloc(seqs,sizeof(char*)*nSaved,sizeof(char*)*nSeq);
  } else {
    /* PHYLIP interleaved-like format
       Allow arbitrary length names, require spaces between names and sequences
       Allow multiple alignments, either separated by a single empty line (e.g. seqboot output)
       or not.
     */
    if (buf[0] == '\n' || buf[0] == '\r') {
      if (fgets(buf,sizeof(buf),fp) == NULL) {
	fprintf(stderr, "Empty header line followed by EOF\n");
	exit(1);
      }
    }
    if (sscanf(buf, "%d%d", &nSeq, &nPos) != 2
      || nSeq < 1 || nPos < 1) {
      fprintf(stderr, "Error parsing header line:%s\n", buf);
      exit(1);
    }
    names = (char **)mymalloc(sizeof(char*) * nSeq);
    seqs = (char **)mymalloc(sizeof(char*) * nSeq);
    nSaved = nSeq;

    int i;
    for (i = 0; i < nSeq; i++) {
      names[i] = NULL;
      seqs[i] = (char *)mymalloc(nPos+1);	/* null-terminate */
      seqs[i][0] = '\0';
    }
    int iSeq = 0;
    
    while(fgets(buf,sizeof(buf),fp)) {
      if ((buf[0] == '\n' || buf[0] == '\r') && (iSeq == nSeq || iSeq == 0)) {
	iSeq = 0;
      } else {
	int j = 0; /* character just past end of name */
	if (buf[0] == ' ') {
	  if (names[iSeq] == NULL) {
	    fprintf(stderr, "No name in phylip line %s", buf);
	    exit(1);
	  }
	} else {
	  while (buf[j] != '\n' && buf[j] != '\0' && buf[j] != ' ')
	    j++;
	  if (buf[j] != ' ' || j == 0) {
	    fprintf(stderr, "No sequence in phylip line %s", buf);
	    exit(1);
	  }
	  if (iSeq >= nSeq) {
	    fprintf(stderr, "No empty line between sequence blocks (is the sequence count wrong?)\n");
	    exit(1);
	  }
	  if (names[iSeq] == NULL) {
	    /* save the name */
	    names[iSeq] = (char *)mymalloc(j+1);
	    int k;
	    for (k = 0; k < j; k++) names[iSeq][k] = buf[k];
	    names[iSeq][j] = '\0';
	  } else {
	    /* check the name */
	    int k;
	    int match = 1;
	    for (k = 0; k < j; k++) {
	      if (names[iSeq][k] != buf[k]) {
		match = 0;
		break;
	      }
	    }
	    if (!match || names[iSeq][j] != '\0') {
	      fprintf(stderr, "Wrong name in phylip line %s\nExpected %s\n", buf, names[iSeq]);
	      exit(1);
	    }
	  }
	}
	int seqlen = strlen(seqs[iSeq]);
	for (; buf[j] != '\n' && buf[j] != '\0'; j++) {
	  if (buf[j] != ' ') {
	    if (seqlen >= nPos) {
	      fprintf(stderr, "Too many characters (expected %d) for sequence named %s\nSo far have:\n%s\n",
		      nPos, names[iSeq], seqs[iSeq]);
	      exit(1);
	    }
	    seqs[iSeq][seqlen++] = toupper(buf[j]);
	  }
	}
	seqs[iSeq][seqlen] = '\0'; /* null-terminate */
	if(verbose>10) fprintf(stderr,"Read iSeq %d name %s seqsofar %s\n", iSeq, names[iSeq], seqs[iSeq]);
	iSeq++;
	if (iSeq == nSeq && strlen(seqs[0]) == nPos)
	  break; /* finished alignment */
      } /* end else non-empty phylip line */
    }
    if (iSeq != nSeq && iSeq != 0) {
      fprintf(stderr, "Wrong number of sequences: expected %d\n", nSeq);
      exit(1);
    }
  }
  /* Check lengths of sequences */
  int i;
  for (i = 0; i < nSeq; i++) {
    int seqlen = strlen(seqs[i]);
    if (seqlen != nPos) {
      fprintf(stderr, "Wrong number of characters for %s: expected %d have %d\n", names[i], nPos, seqlen);
      exit(1);
    }
  }
  /* Replace "." with "-" and warn if we find any */
  /* If nucleotide sequences, replace U with T and N with X */
  bool findDot = false;
  for (i = 0; i < nSeq; i++) {
    char *p;
    for (p = seqs[i]; *p != '\0'; p++) {
      if (*p == '.') {
	findDot = true;
	*p = '-';
      }
      if (nCodes == 4 && *p == 'U')
	*p = 'T';
      if (nCodes == 4 && *p == 'N')
	*p = 'X';
    }
  }
  if (findDot)
    fprintf(stderr, "Warning! Found \".\" character(s). These are treated as gaps\n");

  if (ferror(fp)) {
    fprintf(stderr, "Error reading input file\n");
    exit(1);
  }

  alignment_t *align = (alignment_t*)mymalloc(sizeof(alignment_t));
  align->nSeq = nSeq;
  align->nPos = nPos;
  align->names = names;
  align->seqs = seqs;
  align->nSaved = nSaved;
  return(align);
}

alignment_t *FreeAlignment(alignment_t *aln) {
  if(aln==NULL)
    return(NULL);
  int i;
  for (i = 0; i < aln->nSeq; i++) {
    aln->names[i] = myfree(aln->names[i],strlen(aln->names[i])+1);
    aln->seqs[i] = myfree(aln->seqs[i], aln->nPos+1);
  }
  aln->names = myfree(aln->names, sizeof(char*)*aln->nSaved);
  aln->seqs = myfree(aln->seqs, sizeof(char*)*aln->nSaved);
  myfree(aln, sizeof(alignment_t));
  return(NULL);
}

char **AlnToConstraints(alignment_t *constraints, uniquify_t *unique, hashstrings_t *hashnames) {
  /* look up constraints as names and map to unique-space */
  char **  uniqConstraints = (char**)mymalloc(sizeof(char*) * unique->nUnique);	
  int i;
  for (i = 0; i < unique->nUnique; i++)
    uniqConstraints[i] = NULL;
  for (i = 0; i < constraints->nSeq; i++) {
    char *name = constraints->names[i];
    char *constraintSeq = constraints->seqs[i];
    hashiterator_t hi = FindMatch(hashnames,name);
    if (HashCount(hashnames,hi) != 1) {
      fprintf(stderr, "Sequence %s from constraints file is not in the alignment\n", name);
      exit(1);
    }
    int iSeqNonunique = HashFirst(hashnames,hi);
    assert(iSeqNonunique >= 0 && iSeqNonunique < unique->nSeq);
    int iSeqUnique = unique->alnToUniq[iSeqNonunique];
    assert(iSeqUnique >= 0 && iSeqUnique < unique->nUnique);
    if (uniqConstraints[iSeqUnique] != NULL) {
      /* Already set a constraint for this group of sequences!
	 Warn that we are ignoring this one unless the constraints match */
      if (strcmp(uniqConstraints[iSeqUnique],constraintSeq) != 0) {
	fprintf(stderr,
		"Warning: ignoring constraints for %s:\n%s\n"
		"Another sequence has the same sequence but different constraints\n",
		name, constraintSeq);
      }
    } else {
      uniqConstraints[iSeqUnique] = constraintSeq;
    }
  }
  return(uniqConstraints);
}


profile_t *SeqToProfile(/*IN/OUT*/NJ_t *NJ,
			char *seq, int nPos,
			/*OPTIONAL*/char *constraintSeq, int nConstraints,
			int iNode,
			unsigned long counts[256]) {
  static unsigned char charToCode[256];
  static int codeSet = 0;
  int c, i;

  if (!codeSet) {
    for (c = 0; c < 256; c++) {
      charToCode[c] = nCodes;
    }
    for (i = 0; codesString[i]; i++) {
      charToCode[codesString[i]] = i;
      charToCode[tolower(codesString[i])] = i;
    }
    charToCode['-'] = NOCODE;
    codeSet=1;
  }

  assert(strlen(seq) == nPos);
  profile_t *profile = NewProfile(nPos,nConstraints);

  for (i = 0; i < nPos; i++) {
    unsigned int character = (unsigned int) seq[i];
    counts[character]++;
    c = charToCode[character];
    if(verbose>10 && i < 2) fprintf(stderr,"pos %d char %c code %d\n", i, seq[i], c);
    /* treat unknowns as gaps */
    if (c == nCodes || c == NOCODE) {
      profile->codes[i] = NOCODE;
      profile->weights[i] = 0.0;
    } else {
      profile->codes[i] = c;
      profile->weights[i] = 1.0;
    }
  }
  if (nConstraints > 0) {
    for (i = 0; i < nConstraints; i++) {
      profile->nOn[i] = 0;
      profile->nOff[i] = 0;
    }
    bool bWarn = false;
    if (constraintSeq != NULL) {
      assert(strlen(constraintSeq) == nConstraints);
      for (i = 0; i < nConstraints; i++) {
	if (constraintSeq[i] == '1') {
	  profile->nOn[i] = 1;
	} else if (constraintSeq[i] == '0') {
	  profile->nOff[i] = 1;
	} else if (constraintSeq[i] != '-') {
	  if (!bWarn) {
	    fprintf(stderr, "Constraint characters in unique sequence %d replaced with gap:", iNode+1);
	    bWarn = true;
	  }
	  fprintf(stderr, " %c%d", constraintSeq[i], i+1);
	  /* For the benefit of ConstraintSequencePenalty -- this is a bit of a hack, as
	     this modifies the value read from the alignment
	  */
	  constraintSeq[i] = '-';
	}
      }
      if (bWarn)
	fprintf(stderr, "\n");
    }
  }
  return profile;
}

void SeqDist(unsigned char *codes1, unsigned char *codes2, int nPos,
	     distance_matrix_t *dmat, 
	     /*OUT*/besthit_t *hit) {
  double top = 0;		/* summed over positions */
  int nUse = 0;
  int i;
  if (dmat==NULL) {
    int nDiff = 0;
    for (i = 0; i < nPos; i++) {
      if (codes1[i] != NOCODE && codes2[i] != NOCODE) {
	nUse++;
	if (codes1[i] != codes2[i]) nDiff++;
      }
    }
    top = (double)nDiff;
  } else {
    for (i = 0; i < nPos; i++) {
      if (codes1[i] != NOCODE && codes2[i] != NOCODE) {
	nUse++;
	top += dmat->distances[(unsigned int)codes1[i]][(unsigned int)codes2[i]];
      }
    }
  }
  hit->weight = (double)nUse;
  hit->dist = nUse > 0 ? top/(double)nUse : 1.0;
  seqOps++;
}

void CorrectedPairDistances(profile_t **profiles, int nProfiles,
			    /*OPTIONAL*/distance_matrix_t *distance_matrix,
			    int nPos,
			    /*OUT*/double *distances) {
  assert(distances != NULL);
  assert(profiles != NULL);
  assert(nProfiles>1 && nProfiles <= 4);
  besthit_t hit[6];
  int iHit,i,j;

  for (iHit=0, i=0; i < nProfiles; i++) {
    for (j=i+1; j < nProfiles; j++, iHit++) {
      ProfileDist(profiles[i],profiles[j],nPos,distance_matrix,/*OUT*/&hit[iHit]);
      distances[iHit] = hit[iHit].dist;
    }
  }
  if (pseudoWeight > 0) {
    /* Estimate the prior distance */
    double dTop = 0;
    double dBottom = 0;
    for (iHit=0; iHit < (nProfiles*(nProfiles-1))/2; iHit++) {
      dTop += hit[iHit].dist * hit[iHit].weight;
      dBottom += hit[iHit].weight;
    }
    double prior = (dBottom > 0.01) ? dTop/dBottom : 3.0;
    for (iHit=0; iHit < (nProfiles*(nProfiles-1))/2; iHit++)
      distances[iHit] = (distances[iHit] * hit[iHit].weight + prior * pseudoWeight)
	/ (hit[iHit].weight + pseudoWeight);
  }
  if (logdist) {
    for (iHit=0; iHit < (nProfiles*(nProfiles-1))/2; iHit++)
      distances[iHit] = LogCorrect(distances[iHit]);
  }
}

/* During the neighbor-joining phase, a join only violates our constraints if
   node1, node2, and other are all represented in the constraint
   and if one of the 3 is split and the other two do not agree
 */
int JoinConstraintPenalty(/*IN*/NJ_t *NJ, int node1, int node2) {
  if (NJ->nConstraints == 0)
    return(0.0);
  int penalty = 0;
  int iC;
  for (iC = 0; iC < NJ->nConstraints; iC++)
    penalty += JoinConstraintPenaltyPiece(NJ, node1, node2, iC);
  return(penalty);
}

int JoinConstraintPenaltyPiece(NJ_t *NJ, int node1, int node2, int iC) {
  profile_t *pOut = NJ->outprofile;
  profile_t *p1 = NJ->profiles[node1];
  profile_t *p2 = NJ->profiles[node2];
  int nOn1 = p1->nOn[iC];
  int nOff1 = p1->nOff[iC];
  int nOn2 = p2->nOn[iC];
  int nOff2 = p2->nOff[iC];
  int nOnOut = pOut->nOn[iC] - nOn1 - nOn2;
  int nOffOut = pOut->nOff[iC] - nOff1 - nOff2;

  if ((nOn1+nOff1) > 0 && (nOn2+nOff2) > 0 && (nOnOut+nOffOut) > 0) {
    /* code is -1 for split, 0 for off, 1 for on */
    int code1 = (nOn1 > 0 && nOff1 > 0) ? -1 : (nOn1 > 0 ? 1 : 0);
    int code2 = (nOn2 > 0 && nOff2 > 0) ? -1 : (nOn2 > 0 ? 1 : 0);
    int code3 = (nOnOut > 0 && nOffOut) > 0 ? -1 : (nOnOut > 0 ? 1 : 0);
    int nSplit = (code1 == -1 ? 1 : 0) + (code2 == -1 ? 1 : 0) + (code3 == -1 ? 1 : 0);
    int nOn = (code1 == 1 ? 1 : 0) + (code2 == 1 ? 1 : 0) + (code3 == 1 ? 1 : 0);
    if (nSplit == 1 && nOn == 1)
      return(SplitConstraintPenalty(nOn1+nOn2, nOff1+nOff2, nOnOut, nOffOut));
  }
  /* else */
  return(0);
}

void QuartetConstraintPenalties(profile_t *profiles[4], int nConstraints, /*OUT*/double penalty[3]) {
  int i;
  for (i=0; i < 3; i++)
    penalty[i] = 0.0;
  if(nConstraints == 0)
    return;
  int iC;
  for (iC = 0; iC < nConstraints; iC++) {
    double part[3];
    if (QuartetConstraintPenaltiesPiece(profiles, iC, /*OUT*/part)) {
      for (i=0;i<3;i++)
	penalty[i] += part[i];

      if (verbose>2
	  && (fabs(part[ABvsCD]-part[ACvsBD]) > 0.001 || fabs(part[ABvsCD]-part[ADvsBC]) > 0.001))
	fprintf(stderr, "Constraint Penalties at %d: ABvsCD %.3f ACvsBD %.3f ADvsBC %.3f %d/%d %d/%d %d/%d %d/%d\n",
		iC, part[ABvsCD], part[ACvsBD], part[ADvsBC],
		profiles[0]->nOn[iC], profiles[0]->nOff[iC],
		profiles[1]->nOn[iC], profiles[1]->nOff[iC],
		profiles[2]->nOn[iC], profiles[2]->nOff[iC],
		profiles[3]->nOn[iC], profiles[3]->nOff[iC]);
    }
  }
  if (verbose>2)
    fprintf(stderr, "Total Constraint Penalties: ABvsCD %.3f ACvsBD %.3f ADvsBC %.3f\n",
	    penalty[ABvsCD], penalty[ACvsBD], penalty[ADvsBC]);
}

double PairConstraintDistance(int nOn1, int nOff1, int nOn2, int nOff2) {
  double f1 = nOn1/(double)(nOn1+nOff1);
  double f2 = nOn2/(double)(nOn2+nOff2);
  /* 1 - f1 * f2 - (1-f1)*(1-f2) = 1 - f1 * f2 - 1 + f1 + f2 - f1 * f2 */
  return(f1 + f2 - 2.0 * f1 * f2);
}

bool QuartetConstraintPenaltiesPiece(profile_t *profiles[4], int iC, /*OUT*/double piece[3]) {
  int nOn[4];
  int nOff[4];
  int i;
  int nSplit = 0;
  int nPlus = 0;
  int nMinus = 0;
  
  for (i=0; i < 4; i++) {
    nOn[i] = profiles[i]->nOn[iC];
    nOff[i] = profiles[i]->nOff[iC];
    if (nOn[i] + nOff[i] == 0)
      return(false);		/* ignore */
    else if (nOn[i] > 0 && nOff[i] > 0)
      nSplit++;
    else if (nOn[i] > 0)
      nPlus++;
    else
      nMinus++;
  }
  /* If just one of them is split or on the other side and the others all agree, also ignore */
  if (nPlus >= 3 || nMinus >= 3)
    return(false);
  piece[ABvsCD] = constraintWeight
    * (PairConstraintDistance(nOn[0],nOff[0],nOn[1],nOff[1])
       + PairConstraintDistance(nOn[2],nOff[2],nOn[3],nOff[3]));
  piece[ACvsBD] = constraintWeight
    * (PairConstraintDistance(nOn[0],nOff[0],nOn[2],nOff[2])
       + PairConstraintDistance(nOn[1],nOff[1],nOn[3],nOff[3]));
  piece[ADvsBC] = constraintWeight
    * (PairConstraintDistance(nOn[0],nOff[0],nOn[3],nOff[3])
       + PairConstraintDistance(nOn[2],nOff[2],nOn[1],nOff[1]));
  return(true);
}

/* Minimum number of constrained leaves that need to be moved
   to satisfy the constraint (or 0 if constraint is satisfied)
   Defining it this way should ensure that SPR moves that break
   constraints get a penalty
*/
int SplitConstraintPenalty(int nOn1, int nOff1, int nOn2, int nOff2) {
  return(nOn1 + nOff2 < nOn2 + nOff1 ?
	 (nOn1 < nOff2 ? nOn1 : nOff2)
	 : (nOn2 < nOff1 ? nOn2 : nOff1));
}

bool SplitViolatesConstraint(profile_t *profiles[4], int iConstraint) {
  int i;
  int codes[4]; /* 0 for off, 1 for on, -1 for split (quit if not constrained at all) */
  for (i = 0; i < 4; i++) {
    if (profiles[i]->nOn[iConstraint] + profiles[i]->nOff[iConstraint] == 0)
      return(false);
    else if (profiles[i]->nOn[iConstraint] > 0 && profiles[i]->nOff[iConstraint] == 0)
      codes[i] = 1;
    else if (profiles[i]->nOn[iConstraint] == 0 && profiles[i]->nOff[iConstraint] > 0)
      codes[i] = 0;
    else
      codes[i] = -1;
  }
  int n0 = 0;
  int n1 = 0;
  for (i = 0; i < 4; i++) {
    if (codes[i] == 0)
      n0++;
    else if (codes[i] == 1)
      n1++;
  }
  /* 3 on one side means no violation, even if other is code -1
     otherwise must have code != -1 and agreement on the split
   */
  if (n0 >= 3 || n1 >= 3)
    return(false);
  if (n0==2 && n1==2 && codes[0] == codes[1] && codes[2] == codes[3])
    return(false);
  return(true);
}

double LogCorrect(double dist) {
  const double maxscore = 3.0;
  if (nCodes == 4 && !useMatrix) { /* Jukes-Cantor */
    dist = dist < 0.74 ? -0.75*log(1.0 - dist * 4.0/3.0) : maxscore;
  } else {			/* scoredist-like */
    dist = dist < 0.99 ? -1.3*log(1.0 - dist) : maxscore;
  }
  return (dist < maxscore ? dist : maxscore);
}

/* A helper function -- f1 and f2 can be NULL if the corresponding code != NOCODE
*/
double ProfileDistPiece(unsigned int code1, unsigned int code2,
			float *f1, float *f2, 
			/*OPTIONAL*/distance_matrix_t *dmat,
			/*OPTIONAL*/float *codeDist2) {
  if (dmat) {
    if (code1 != NOCODE && code2 != NOCODE) { /* code1 vs code2 */
      return(dmat->distances[code1][code2]);
    } else if (codeDist2 != NULL && code1 != NOCODE) { /* code1 vs. codeDist2 */
      return(codeDist2[code1]);
    } else { /* f1 vs f2 */
      int k;
      if (f1 == NULL) {
	if(code1 == NOCODE) return(10.0);
	f1 = &dmat->codeFreq[code1][0];
      }
      if (f2 == NULL) {
	if(code2 ==NOCODE) return(10.0);
	f2 = &dmat->codeFreq[code2][0];
      }
      double piece = 0;
      for (k = 0; k < nCodes; k++)
	piece += f1[k] * f2[k] * dmat->eigenval[k];
      return(piece);
    }
  } else {
    /* no matrix */
    if (code1 != NOCODE) {
      if (code2 != NOCODE) {
	return(code1 == code2 ? 0.0 : 1.0); /* code1 vs code2 */
      } else {
	if(f2 == NULL) return(10.0);
	return(1.0 - f2[code1]); /* code1 vs. f2 */
      }
    } else {
      if (code2 != NOCODE) {
	if(f1 == NULL) return(10.0);
	return(1.0 - f1[code2]); /* f1 vs code2 */
      } else { /* f1 vs. f2 */
	if (f1 == NULL || f2 == NULL) return(10.0);
	double piece = 1.0;
	int k;
	for (k = 0; k < nCodes; k++) {
	  piece -= f1[k] * f2[k];
	}
	return(piece);
      }
    }
  }
  assert(0);
}

/* E.g. GET_FREQ(profile,iPos,iVector)
   Gets the next element of the vectors (and updates iVector), or
   returns NULL if we didn't store a vector
*/
#define GET_FREQ(P,I,IVECTOR) \
(P->weights[I] > 0 && P->codes[I] == NOCODE ? &P->vectors[nCodes*(IVECTOR++)] : NULL)

void ProfileDist(profile_t *profile1, profile_t *profile2, int nPos,
		 /*OPTIONAL*/distance_matrix_t *dmat,
		 /*OUT*/besthit_t *hit) {
  double top = 0;
  double denom = 0;
  int iFreq1 = 0;
  int iFreq2 = 0;
  int i = 0;
  for (i = 0; i < nPos; i++) {
      float *f1 = GET_FREQ(profile1,i,/*IN/OUT*/iFreq1);
      float *f2 = GET_FREQ(profile2,i,/*IN/OUT*/iFreq2);
      if (profile1->weights[i] > 0 && profile2->weights[i] > 0) {
	double weight = profile1->weights[i] * profile2->weights[i];
	denom += weight;
	double piece = ProfileDistPiece(profile1->codes[i],profile2->codes[i],f1,f2,dmat,
					profile2->codeDist ? &profile2->codeDist[i*nCodes] : NULL);
	top += weight * piece;
      }
  }
  assert(iFreq1 == profile1->nVectors);
  assert(iFreq2 == profile2->nVectors);
  hit->weight = denom > 0 ? denom : 0.01; /* 0.01 is an arbitrarily low value of weight (normally >>1) */
  hit->dist = denom > 0 ? top/denom : 1;
  profileOps++;
}

/* This should not be called if the update weight is 0, as
   in that case code==NOCODE and in=NULL is possible, and then
   it will fail.
*/
void AddToFreq(/*IN/OUT*/float *fOut,
	       double weight,
	       unsigned int codeIn, /*OPTIONAL*/float *fIn,
	       /*OPTIONAL*/distance_matrix_t *dmat) {
  int k;
  assert(fOut != NULL);
  if (fIn != NULL) {
    for (k = 0; k < nCodes; k++)
      fOut[k] += fIn[k] * weight;
  } else if (dmat) {
    assert(codeIn != NOCODE);
    for (k = 0; k < nCodes; k++)
      fOut[k] += dmat->codeFreq[codeIn][k] * weight;
  } else {
    assert(codeIn != NOCODE);
    fOut[codeIn] += weight;
  }
}

void SetProfile(/*IN/OUT*/NJ_t *NJ, int node, double weight1) {
    children_t *c = &NJ->child[node];
    assert(c->nChild == 2);
    assert(NJ->profiles[c->child[0]] != NULL);
    assert(NJ->profiles[c->child[1]] != NULL);
    if (NJ->profiles[node] != NULL)
      FreeProfile(NJ->profiles[node], NJ->nPos, NJ->nConstraints);
    NJ->profiles[node] = AverageProfile(NJ->profiles[c->child[0]],
					NJ->profiles[c->child[1]],
					NJ->nPos, NJ->nConstraints,
					NJ->distance_matrix,
					weight1);
}

/* bionjWeight is the weight of the first sequence (between 0 and 1),
   or -1 to do the average.
   */
profile_t *AverageProfile(profile_t *profile1, profile_t *profile2,
			  int nPos, int nConstraints,
			  distance_matrix_t *dmat,
			  double bionjWeight) {
  int i;
  if (bionjWeight < 0) {
    bionjWeight = 0.5;
  }

  /* First, set codes and weights and see how big vectors will be */
  profile_t *out = NewProfile(nPos, nConstraints);

  for (i = 0; i < nPos; i++) {
    out->weights[i] = bionjWeight * profile1->weights[i]
      + (1-bionjWeight) * profile2->weights[i];
    out->codes[i] = NOCODE;
    if (out->weights[i] > 0) {
      if (profile1->weights[i] > 0 && profile1->codes[i] != NOCODE
	  && (profile2->weights[i] <= 0 || profile1->codes[i] == profile2->codes[i])) {
	out->codes[i] = profile1->codes[i];
      } else if (profile1->weights[i] <= 0
		 && profile2->weights[i] > 0
		 && profile2->codes[i] != NOCODE) {
	out->codes[i] = profile2->codes[i];
      }
      if (out->codes[i] == NOCODE) out->nVectors++;
    }
  }

  /* Allocate and set the vectors */
  out->vectors = (float*)mymalloc(sizeof(float)*nCodes*out->nVectors);
  for (i = 0; i < nCodes * out->nVectors; i++) out->vectors[i] = 0;
  nProfileFreqAlloc += out->nVectors;
  nProfileFreqAvoid += nPos - out->nVectors;
  int iFreqOut = 0;
  int iFreq1 = 0;
  int iFreq2 = 0;
  for (i=0; i < nPos; i++) {
    float *f = GET_FREQ(out,i,/*IN/OUT*/iFreqOut);
    float *f1 = GET_FREQ(profile1,i,/*IN/OUT*/iFreq1);
    float *f2 = GET_FREQ(profile2,i,/*IN/OUT*/iFreq2);
    if (f != NULL) {
      if (profile1->weights[i] > 0)
	AddToFreq(/*IN/OUT*/f, profile1->weights[i] * bionjWeight,
		  profile1->codes[i], f1, dmat);
      if (profile2->weights[i] > 0)
	AddToFreq(/*IN/OUT*/f, profile2->weights[i] * (1.0-bionjWeight),
		  profile2->codes[i], f2, dmat);
      NormalizeFreq(/*IN/OUT*/f, dmat);
    } /* end if computing f */
    if (verbose > 10 && i < 5) {
      fprintf(stderr,"Average profiles: pos %d in-w1 %f in-w2 %f bionjWeight %f to weight %f code %d\n",
	      i, profile1->weights[i], profile2->weights[i], bionjWeight,
	      out->weights[i], out->codes[i]);
      if (f!= NULL) {
	int k;
	for (k = 0; k < nCodes; k++)
	  fprintf(stderr, "\t%c:%f", codesString[k], f ? f[k] : -1.0);
	fprintf(stderr,"\n");
      }
    }
  } /* end loop over positions */
  assert(iFreq1 == profile1->nVectors);
  assert(iFreq2 == profile2->nVectors);
  assert(iFreqOut == out->nVectors);

  /* compute total constraints */
  for (i = 0; i < nConstraints; i++) {
    out->nOn[i] = profile1->nOn[i] + profile2->nOn[i];
    out->nOff[i] = profile1->nOff[i] + profile2->nOff[i];
  }
  profileAvgOps++;
  return(out);
}

/* Make the (unrotated) frequencies sum to 1
   Simply dividing by total_weight is not ideal because of roundoff error
   So compute total_freq instead
*/
void NormalizeFreq(/*IN/OUT*/float *freq, distance_matrix_t *dmat) {
  double total_freq = 0;
  int k;
  if (dmat != NULL) {
    /* The total frequency is dot_product(true_frequencies, 1)
       So we rotate the 1 vector by eigeninv (stored in eigentot)
    */
    for (k = 0; k < nCodes; k++) {
      total_freq += freq[k] * dmat->eigentot[k];
    }
  } else {
    for (k = 0; k < nCodes; k++)
      total_freq += freq[k];
  }
  if (total_freq > 1e-10) {
    double inverse_weight = 1.0/total_freq;
    for (k = 0; k < nCodes; k++)
      freq[k] *= inverse_weight;
  }
}

/* OutProfile() computes the out-profile */
profile_t *OutProfile(profile_t **profiles, int nProfiles,
		      int nPos, int nConstraints,
		      distance_matrix_t *dmat) {
  int i;			/* position */
  int in;			/* profile */
  profile_t *out = NewProfile(nPos, nConstraints);

  double inweight = 1.0/(double)nProfiles;   /* The maximal output weight is 1.0 */

  /* First, set weights -- code is always NOCODE, prevent weight=0 */
  for (i = 0; i < nPos; i++) {
    out->weights[i] = 0;
    for (in = 0; in < nProfiles; in++)
      out->weights[i] += profiles[in]->weights[i] * inweight;
    if (out->weights[i] <= 0) out->weights[i] = 1e-20; /* always store a vector */
    out->nVectors++;
    out->codes[i] = NOCODE;		/* outprofile is normally complicated */
  }

  /* Initialize the frequencies to 0 */
  out->vectors = (float*)mymalloc(sizeof(float)*nCodes*out->nVectors);
  for (i = 0; i < nCodes*out->nVectors; i++)
    out->vectors[i] = 0;

  /* Add up the weights, going through each sequence in turn */
  for (in = 0; in < nProfiles; in++) {
    int iFreqOut = 0;
    int iFreqIn = 0;
    for (i = 0; i < nPos; i++) {
      float *fIn = GET_FREQ(profiles[in],i,/*IN/OUT*/iFreqIn);
      float *fOut = GET_FREQ(out,i,/*IN/OUT*/iFreqOut);
      if (profiles[in]->weights[i] > 0)
	AddToFreq(/*IN/OUT*/fOut, profiles[in]->weights[i],
		  profiles[in]->codes[i], fIn, dmat);
    }
    assert(iFreqOut == out->nVectors);
    assert(iFreqIn == profiles[in]->nVectors);
  }

  /* And normalize the frequencies to sum to 1 */
  int iFreqOut = 0;
  for (i = 0; i < nPos; i++) {
    float *fOut = GET_FREQ(out,i,/*IN/OUT*/iFreqOut);
    if (fOut)
      NormalizeFreq(/*IN/OUT*/fOut, dmat);
  }
  assert(iFreqOut == out->nVectors);
  if (verbose > 10) fprintf(stderr,"Average %d profiles\n", nProfiles);
  if(dmat)
    SetCodeDist(/*IN/OUT*/out, nPos, dmat);

  /* Compute constraints */
  for (i = 0; i < nConstraints; i++) {
    out->nOn[i] = 0;
    out->nOff[i] = 0;
    for (in = 0; in < nProfiles; in++) {
      out->nOn[i] += profiles[in]->nOn[i];
      out->nOff[i] += profiles[in]->nOff[i];
    }
  }
  return(out);
}

void UpdateOutProfile(/*IN/OUT*/profile_t *out, profile_t *old1, profile_t *old2,
		      profile_t *new, int nActiveOld,
		      int nPos, int nConstraints,
		      distance_matrix_t *dmat) {
  int i, k;
  int iFreqOut = 0;
  int iFreq1 = 0;
  int iFreq2 = 0;
  int iFreqNew = 0;
  assert(nActiveOld > 0);

  for (i = 0; i < nPos; i++) {
    float *fOut = GET_FREQ(out,i,/*IN/OUT*/iFreqOut);
    float *fOld1 = GET_FREQ(old1,i,/*IN/OUT*/iFreq1);
    float *fOld2 = GET_FREQ(old2,i,/*IN/OUT*/iFreq2);
    float *fNew = GET_FREQ(new,i,/*IN/OUT*/iFreqNew);

    assert(out->codes[i] == NOCODE && fOut != NULL); /* No no-vector optimization for outprofiles */
    if (verbose > 2 && i < 3) {
      fprintf(stderr,"Updating out-profile position %d weight %f (mult %f)\n",
	      i, out->weights[i], out->weights[i]*nActiveOld);
    }
    double originalMult = out->weights[i]*nActiveOld;
    double newMult = originalMult + new->weights[i] - old1->weights[i] - old2->weights[i];
    out->weights[i] = newMult/(nActiveOld-1);
    if (out->weights[i] <= 0) out->weights[i] = 1e-20; /* always use the vector */

    for (k = 0; k < nCodes; k++) fOut[k] *= originalMult;
    
    if (old1->weights[i] > 0)
      AddToFreq(/*IN/OUT*/fOut, -old1->weights[i], old1->codes[i], fOld1, dmat);
    if (old2->weights[i] > 0)
      AddToFreq(/*IN/OUT*/fOut, -old2->weights[i], old2->codes[i], fOld2, dmat);
    if (new->weights[i] > 0)
      AddToFreq(/*IN/OUT*/fOut, new->weights[i], new->codes[i], fNew, dmat);

    /* And renormalize */
    NormalizeFreq(/*IN/OUT*/fOut, dmat);

    if (verbose > 2 && i < 3) {
      fprintf(stderr,"Updated out-profile position %d weight %f (mult %f)",
	      i, out->weights[i], out->weights[i]*nActiveOld);
      if(out->weights[i] > 0)
	for (k=0;k<nCodes;k++)
	  fprintf(stderr, " %c:%f", dmat?'?':codesString[k], fOut[k]);
      fprintf(stderr,"\n");
    }
  }
  assert(iFreqOut == out->nVectors);
  assert(iFreq1 == old1->nVectors);
  assert(iFreq2 == old2->nVectors);
  assert(iFreqNew == new->nVectors);
  if(dmat)
    SetCodeDist(/*IN/OUT*/out,nPos,dmat);

  /* update constraints -- note in practice this should be a no-op */
  for (i = 0; i < nConstraints; i++) {
    out->nOn[i] += new->nOn[i] - old1->nOn[i] - old2->nOn[i];
    out->nOff[i] += new->nOff[i] - old1->nOff[i] - old2->nOff[i];
  }
}

void SetCodeDist(/*IN/OUT*/profile_t *profile, int nPos,
			   distance_matrix_t *dmat) {
  if (profile->codeDist == NULL)
    profile->codeDist = (float*)mymalloc(sizeof(float)*nPos*nCodes);
  int i;
  int iFreq = 0;
  for (i = 0; i < nPos; i++) {
    float *f = GET_FREQ(profile,i,/*IN/OUT*/iFreq);

    int k;
    for (k = 0; k < nCodes; k++)
      profile->codeDist[i*nCodes+k] = ProfileDistPiece(/*code1*/profile->codes[i], /*code2*/k,
						       /*f1*/f, /*f2*/NULL,
						       dmat, NULL);
  }
  assert(iFreq==profile->nVectors);
}


void SetBestHit(int node, NJ_t *NJ, int nActive,
		/*OUT*/besthit_t *bestjoin, /*OUT OPTIONAL*/besthit_t *allhits) {
  assert(NJ->parent[node] <  0);

  bestjoin->i = node;
  bestjoin->j = -1;
  bestjoin->dist = 1e20;
  bestjoin->criterion = 1e20;

  int j;
  besthit_t tmp;

  for (j = 0; j < NJ->maxnode; j++) {
    besthit_t *sv = allhits != NULL ? &allhits[j] : &tmp;
    sv->i = node;
    sv->j = j;
    if (NJ->parent[j] >= 0) {
      sv->weight = 0.0;
      sv->criterion = sv->dist = 1e20;
      continue;
    }
    /* Note that we compute self-distances (allow j==node) because the top-hit heuristic
       expects self to be within its top hits, but we exclude those from the bestjoin
       that we return...
    */
    SetDistCriterion(NJ, nActive, /*IN/OUT*/sv);
    if (sv->criterion < bestjoin->criterion && node != j)
      *bestjoin = *sv;
  }
  if (verbose>5) {
    fprintf(stderr, "SetBestHit %d %d %f %f\n", bestjoin->i, bestjoin->j, bestjoin->dist, bestjoin->criterion);
  }
}

void ReadMatrix(char *filename, /*OUT*/float codes[MAXCODES][MAXCODES], bool checkCodes) {
  char buf[BUFFER_SIZE] = "";
  FILE *fp = fopen(filename, "r");
  if (fp == NULL) {
    fprintf(stderr, "Cannot read %s\n",filename);
    exit(1);
  }
  if (fgets(buf,sizeof(buf),fp) == NULL) {
    fprintf(stderr, "Error reading header line for %s:\n%s\n", filename, buf);
    exit(1);
  }
  if (checkCodes) {
    int i;
    int iBufPos;
    for (iBufPos=0,i=0;i<nCodes;i++,iBufPos++) {
      if(buf[iBufPos] != codesString[i]) {
	fprintf(stderr,"Header line\n%s\nin file %s does not have expected code %c # %d in %s\n",
		buf, filename, codesString[i], i, codesString);
	exit(1);
      }
      iBufPos++;
      if(buf[iBufPos] != '\n' && buf[iBufPos] != '\r' && buf[iBufPos] != '\0' && buf[iBufPos] != '\t') {
	fprintf(stderr, "Header line in %s should be tab-delimited\n", filename);
	exit(1);
      }
      if (buf[iBufPos] == '\0' && i < nCodes-1) {
	fprintf(stderr, "Header line in %s ends prematurely\n",filename);
	exit(1);
      }
    } /* end loop over codes */
    /* Should be at end, but allow \n because of potential DOS \r\n */
    if(buf[iBufPos] != '\0' && buf[iBufPos] != '\n' && buf[iBufPos] != '\r') {
      fprintf(stderr, "Header line in %s has too many entries\n", filename);
      exit(1);
    }
  }
  int iLine;
  for (iLine = 0; iLine < nCodes; iLine++) {
    buf[0] = '\0';
    if (fgets(buf,sizeof(buf),fp) == NULL) {
      fprintf(stderr, "Cannot read line %d from file %s\n", iLine+2, filename);
      exit(1);
    }
    char *field = strtok(buf,"\t\r\n");
    field = strtok(NULL, "\t");	/* ignore first column */
    int iColumn;
    for (iColumn = 0; iColumn < nCodes && field != NULL; iColumn++, field = strtok(NULL,"\t")) {
      if(sscanf(field,"%f",&codes[iLine][iColumn]) != 1) {
	fprintf(stderr,"Cannot parse field %s in file %s\n", field, filename);
	exit(1);
      }
    }
  }
}

void ReadVector(char *filename, /*OUT*/float codes[MAXCODES]) {
  FILE *fp = fopen(filename,"r");
  if (fp == NULL) {
    fprintf(stderr, "Cannot read %s\n",filename);
    exit(1);
  }
  int i;
  for (i = 0; i < nCodes; i++) {
    if (fscanf(fp,"%f",&codes[i]) != 1) {
      fprintf(stderr,"Cannot read %d entry of %s\n",i+1,filename);
      exit(1);
    }
  }
  if (fclose(fp) != 0) {
    fprintf(stderr, "Error reading %s\n",filename);
    exit(1);
  }
}

distance_matrix_t *ReadDistanceMatrix(char *prefix) {
  char buffer[BUFFER_SIZE];
  distance_matrix_t *dmat = (distance_matrix_t*)mymalloc(sizeof(distance_matrix_t));

  if(strlen(prefix) > BUFFER_SIZE-20) {
    fprintf(stderr,"Filename %s too long\n", prefix);
    exit(1);
  }

  strcpy(buffer, prefix);
  strcat(buffer, ".distances");
  ReadMatrix(buffer, /*OUT*/dmat->distances, /*checkCodes*/true);

  strcpy(buffer, prefix);
  strcat(buffer, ".inverses");
  ReadMatrix(buffer, /*OUT*/dmat->eigeninv, /*checkCodes*/false);

  strcpy(buffer, prefix);
  strcat(buffer, ".eigenvalues");
  ReadVector(buffer, /*OUT*/dmat->eigenval);

  if(verbose>1) fprintf(stderr, "Read distance matrix from %s\n",prefix);
  SetupDistanceMatrix(/*IN/OUT*/dmat);
  return(dmat);
}

void SetupDistanceMatrix(/*IN/OUT*/distance_matrix_t *dmat) {
  /* Check that the eigenvalues and eigen-inverse are consistent with the
     distance matrix and that the matrix is symmetric */
  int i,j,k;
  for (i = 0; i < nCodes; i++) {
    for (j = 0; j < nCodes; j++) {
      if(fabs(dmat->distances[i][j]-dmat->distances[j][i]) > 1e-6) {
	fprintf(stderr,"Distance matrix not symmetric for %d,%d: %f vs %f\n",
		i+1,j+1,
		dmat->distances[i][j],
		dmat->distances[j][i]);
	exit(1);
      }
      double total = 0.0;
      for (k = 0; k < nCodes; k++)
	total += dmat->eigenval[k] * dmat->eigeninv[k][i] * dmat->eigeninv[k][j];
      if(fabs(total - dmat->distances[i][j]) > 1e-6) {
	fprintf(stderr,"Distance matrix entry %d,%d should be %f but eigen-representation gives %f\n",
		i+1,j+1,dmat->distances[i][j],total);
	exit(1);
      }
    }
  }
  
  /* And compute eigentot */
  for (k = 0; k < nCodes; k++) {
    dmat->eigentot[k] = 0.;
    int j;
    for (j = 0; j < nCodes; j++)
      dmat->eigentot[k] += dmat->eigeninv[k][j];
  }
  
  /* And compute codeFreq */
  int code;
  for(code = 0; code < nCodes; code++) {
    for (k = 0; k < nCodes; k++)
      dmat->codeFreq[code][k] = dmat->eigeninv[k][code];
  }
  if(verbose>10) fprintf(stderr, "Made codeFreq\n");
}


nni_t ChooseNNI(profile_t *profiles[4],
		/*OPTIONAL*/distance_matrix_t *dmat,
		int nPos, int nConstraints,
		/*OUT*/double criteria[3]) {
  double d[6];
  CorrectedPairDistances(profiles, 4, dmat, nPos, /*OUT*/d);
  double penalty[3]; 		/* indexed as nni_t */
  QuartetConstraintPenalties(profiles, nConstraints, /*OUT*/penalty);
  criteria[ABvsCD] = d[qAB] + d[qCD] + penalty[ABvsCD];
  criteria[ACvsBD] = d[qAC] + d[qBD] + penalty[ACvsBD];
  criteria[ADvsBC] = d[qAD] + d[qBC] + penalty[ADvsBC];

  nni_t choice = ABvsCD;
  if (criteria[ACvsBD] < criteria[ABvsCD] && criteria[ACvsBD] <= criteria[ADvsBC]) {
    choice = ACvsBD;
  } else if (criteria[ADvsBC] < criteria[ABvsCD] && criteria[ADvsBC] <= criteria[ACvsBD]) {
    choice = ADvsBC;
  }
  if (verbose > 1 && penalty[choice] > penalty[ABvsCD] + 1e-6) {
    fprintf(stderr, "Worsen constraint: from %.3f to %.3f distance %.3f to %.3f: ",
	    penalty[ABvsCD], penalty[choice],
	    criteria[ABvsCD], choice == ACvsBD ? criteria[ACvsBD] : criteria[ADvsBC]);
    int iC;
    for (iC = 0; iC < nConstraints; iC++) {
      double ppart[3];
      if (QuartetConstraintPenaltiesPiece(profiles, iC, /*OUT*/ppart)) {
	double old_penalty = ppart[ABvsCD];
	double new_penalty = ppart[choice];
	if (new_penalty > old_penalty + 1e-6)
	  fprintf(stderr, " %d (%d/%d %d/%d %d/%d %d/%d)", iC,
		  profiles[0]->nOn[iC], profiles[0]->nOff[iC],
		  profiles[1]->nOn[iC], profiles[1]->nOff[iC],
		  profiles[2]->nOn[iC], profiles[2]->nOff[iC],
		  profiles[3]->nOn[iC], profiles[3]->nOff[iC]);
      }
    }
    fprintf(stderr,"\n");
  }
  if (verbose > 2)
    fprintf(stderr, "NNI scores ABvsCD %.5f ACvsBD %.5f ADvsBC %.5f choice %s\n",
	    criteria[ABvsCD], criteria[ACvsBD], criteria[ADvsBC],
	    choice == ABvsCD ? "AB|CD" : (choice == ACvsBD ? "AC|BD" : "AD|BC"));
  return(choice);
}

double TreeLength(/*IN/OUT*/NJ_t *NJ, bool recomputeProfiles) {
  if (recomputeProfiles) {
    traversal_t traversal2 = InitTraversal(NJ);
    int j = NJ->root;
    while((j = TraversePostorder(j, NJ, /*IN/OUT*/traversal2)) >= 0) {
      /* nothing to do for leaves or root */
      if (j >= NJ->nSeq && j != NJ->root)
	SetProfile(/*IN/OUT*/NJ, j, /*noweight*/-1.0);
    }
    traversal2 = FreeTraversal(traversal2,NJ);
  }
  UpdateBranchLengths(/*IN/OUT*/NJ);
  double total_len = 0;
  int iNode;
  for (iNode = 0; iNode < NJ->maxnode; iNode++)
    total_len += NJ->branchlength[iNode];
  return(total_len);
}

void NNI(/*IN/OUT*/NJ_t *NJ, int iRound, int nRounds) {
  /* For each non-root node N, with children A,B, sibling C, and uncle D,
     we compare the current topology AB|CD to the alternate topologies
     AC|BD and AD|BC, by using log-corrected distances on profiles.
     (If logdist is false, then the log correction is not done.)
     It then selects the best topology.

     Regardless of whether it changes the topology, it recomputes the
     profile for the node, using the pairwise distances and BIONJ-like
     weightings. The parent's profile has changed, but recomputing it
     is not necessary because we will visit it before we need it (we
     use postorder, so we may visit the sibling and its children
     before we visit the parent, but we never consider an ancestor's
     profile, so that is OK). When we change the parent's profile,
     this alters the uncle's up-profile, so we remove that.  Finally,
     if the topology has changed, we remove the up-profiles of the
     nodes.

     NNI() does NOT fix up branch lengths.
  */
  int i;

  if (NJ->nSeq <= 3)
    return;			/* nothing to do */

  /* For each node the upProfile or NULL */
  profile_t **upProfiles = UpProfiles(NJ);
  
  traversal_t traversal = InitTraversal(NJ);
  int node = NJ->root;
  int iDone = 0;
  while((node = TraversePostorder(node, NJ, /*IN/OUT*/traversal)) >= 0) {
    if (node < NJ->nSeq || node == NJ->root)
      continue; /* nothing to do for leaves or root */
    if ((iDone % 100) == 0)
      ProgressReport("NNI round %3d of %3d, %d of %d splits",
		     iRound+1, nRounds, iDone+1, NJ->maxnode - NJ->nSeq);
    iDone++;

    profile_t *profiles[4];
    int nodeABC[3];
    SetupABCD(NJ, node, /*OUT*/profiles, /*IN/OUT*/upProfiles, /*OUT*/nodeABC);

    /* Given our 4 profiles, consider doing a swap */
    int nodeA = nodeABC[0];
    int nodeB = nodeABC[1];
    int nodeC = nodeABC[2];

    if (verbose > 2)
      fprintf(stderr,"Considering NNI around %d: Swap A=%d B=%d C=%d D=out(C) (node parent %d)\n",
	      node, nodeA, nodeB, nodeC, NJ->parent[node]);

    double criteria[3];
    nni_t choice = ChooseNNI(profiles, NJ->distance_matrix, NJ->nPos, NJ->nConstraints,
			     /*OUT*/criteria);
    if (choice != ABvsCD)
      nNNI++;
    if (verbose>1 && choice != ABvsCD)
      fprintf(stderr,"NNI around %d: Swap A=%d B=%d C=%d D=out(C) -- choose %s deltaLen %.4f\n",
	      node, nodeA, nodeB, nodeC, choice == ACvsBD ? "AC|BD" : "AD|BC",
	      criteria[choice] - criteria[ABvsCD]);

    if (choice == ACvsBD) {
      /* swap B and C */
      ReplaceChild(/*IN/OUT*/NJ, node, nodeB, nodeC);
      ReplaceChild(/*IN/OUT*/NJ, NJ->parent[node], nodeC, nodeB);
    } else if (choice == ADvsBC) {
      /* swap A and C */
      ReplaceChild(/*IN/OUT*/NJ, node, nodeA, nodeC);
      ReplaceChild(/*IN/OUT*/NJ, NJ->parent[node], nodeC, nodeA);
    }

    /* update profiles */
    if (choice == ABvsCD) {
      /* No longer needed */
      DeleteUpProfile(upProfiles, NJ, nodeA);
      DeleteUpProfile(upProfiles, NJ, nodeB);
      RecomputeProfile(/*IN/OUT*/NJ, /*IN/OUT*/upProfiles, node);
    } else {
      UpdateForNNI(NJ, node, upProfiles);
    }
  } /* end postorder traversal */
  traversal = FreeTraversal(traversal,NJ);
  if (verbose>1) {
    int nUp = 0;
    for (i = 0; i < NJ->maxnodes; i++)
      if (upProfiles[i] != NULL)
	nUp++;
    fprintf(stderr, "N up profiles at end of NNI:  %d\n", nUp);
  }
  upProfiles = FreeUpProfiles(upProfiles,NJ);
}

int FindSPRSteps(/*IN/OUT*/NJ_t *NJ, 
		 int nodeMove,	 /* the node to move multiple times */
		 int nodeAround, /* sibling or parent of node to NNI to start the chain */
		 /*IN/OUT*/profile_t **upProfiles,
		 /*OUT*/spr_step_t *steps,
		 int maxSteps,
		 bool bFirstAC) {
  int iStep;
  for (iStep = 0; iStep < maxSteps; iStep++) {
    if (NJ->child[nodeAround].nChild != 2)
      break;			/* no further to go */

    /* Consider the NNIs around nodeAround */
    profile_t *profiles[4];
    int nodeABC[3];
    SetupABCD(NJ, nodeAround, /*OUT*/profiles, /*IN/OUT*/upProfiles, /*OUT*/nodeABC);
    double criteria[3];
    (void) ChooseNNI(profiles, NJ->distance_matrix, NJ->nPos, NJ->nConstraints,
		     /*OUT*/criteria);

    /* Do & save the swap */
    spr_step_t *step = &steps[iStep];
    if (iStep == 0 ? bFirstAC : criteria[ACvsBD] < criteria[ADvsBC]) {
      /* swap B & C to put AC together */
      step->deltaLength = criteria[ACvsBD] - criteria[ABvsCD];
      step->nodes[0] = nodeABC[1];
      step->nodes[1] = nodeABC[2];
    } else {
      /* swap AC to put AD together */
      step->deltaLength = criteria[ADvsBC] - criteria[ABvsCD];
      step->nodes[0] = nodeABC[0];
      step->nodes[1] = nodeABC[2];
    }

    if (verbose>1) {
      fprintf(stderr, "SPR chain step %d for %d around %d swap %d %d deltaLen %.5f\n",
	      iStep+1, nodeAround, nodeMove, step->nodes[0], step->nodes[1], step->deltaLength);
      if (verbose>3)
	PrintNJInternal(stderr, NJ);
    }
    ReplaceChild(/*IN/OUT*/NJ, nodeAround, step->nodes[0], step->nodes[1]);
    ReplaceChild(/*IN/OUT*/NJ, NJ->parent[nodeAround], step->nodes[1], step->nodes[0]);
    UpdateForNNI(/*IN/OUT*/NJ, nodeAround, /*IN/OUT*/upProfiles);

    /* set the new nodeAround -- either parent(nodeMove) or sibling(nodeMove) --
       so that it different from current nodeAround
     */
    int newAround[2] = { NJ->parent[nodeMove], Sibling(NJ, nodeMove) };
    if (NJ->parent[nodeMove] == NJ->root)
      RootSiblings(NJ, nodeMove, /*OUT*/newAround);
    assert(newAround[0] == nodeAround || newAround[1] == nodeAround);
    assert(newAround[0] != newAround[1]);
    nodeAround = newAround[newAround[0] == nodeAround ? 1 : 0];
  }
  return(iStep);
}

void UnwindSPRStep(/*IN/OUT*/NJ_t *NJ,
		   /*IN*/spr_step_t *step,
		   /*IN/OUT*/profile_t **upProfiles) {
  int parents[2];
  int i;
  for (i = 0; i < 2; i++) {
    assert(step->nodes[i] >= 0 && step->nodes[i] < NJ->maxnodes);
    parents[i] = NJ->parent[step->nodes[i]];
    assert(parents[i] >= 0);
  }
  assert(parents[0] != parents[1]);
  ReplaceChild(/*IN/OUT*/NJ, parents[0], step->nodes[0], step->nodes[1]);
  ReplaceChild(/*IN/OUT*/NJ, parents[1], step->nodes[1], step->nodes[0]);
  int iYounger = 0;
  if (NJ->parent[parents[0]] == parents[1]) {
    iYounger = 0;
  } else {
    assert(NJ->parent[parents[1]] == parents[0]);
    iYounger = 1;
  }
  UpdateForNNI(/*IN/OUT*/NJ, parents[iYounger], /*IN/OUT*/upProfiles);
}

/* Update the profile of node and its ancestor, and delete nearby out-profiles */
void UpdateForNNI(/*IN/OUT*/NJ_t *NJ, int node, /*IN/OUT*/profile_t **upProfiles) {
  int i;
  if (slow) {
    /* exhaustive update */
    for (i = 0; i < NJ->maxnodes; i++)
      DeleteUpProfile(upProfiles, NJ, i);

    /* update profiles back to root */
    int ancestor;
    for (ancestor = node; ancestor >= 0; ancestor = NJ->parent[ancestor])
      RecomputeProfile(/*IN/OUT*/NJ, upProfiles, ancestor);

    /* remove any up-profiles made while doing that*/
    for (i = 0; i < NJ->maxnodes; i++)
      DeleteUpProfile(upProfiles, NJ, i);
  } else {
    /* if fast, only update around self */
    DeleteUpProfile(upProfiles, NJ, node);
    for (i = 0; i < NJ->child[node].nChild; i++)
      DeleteUpProfile(upProfiles, NJ, NJ->child[node].child[i]);
    assert(node != NJ->root);
    int parent = NJ->parent[node];
    int neighbors[2] = { parent, Sibling(NJ, node) };
    if (parent == NJ->root)
      RootSiblings(NJ, node, /*OUT*/neighbors);
    DeleteUpProfile(upProfiles, NJ, neighbors[0]);
    DeleteUpProfile(upProfiles, NJ, neighbors[1]);
    int uncle = Sibling(NJ, parent);
    if (uncle >= 0)
      DeleteUpProfile(upProfiles, NJ, uncle);
    RecomputeProfile(/*IN/OUT*/NJ, upProfiles, node);
    RecomputeProfile(/*IN/OUT*/NJ, upProfiles, parent);
  }
}

void SPR(/*IN/OUT*/NJ_t *NJ, int maxSPRLength, int iRound, int nRounds) {
  /* Given a non-root node N with children A,B, sibling C, and uncle D,
     we can try to move A by doing three types of moves (4 choices):
     "down" -- swap A with a child of B (if B is not a leaf) [2 choices]
     "over" -- swap B with C
     "up" -- swap A with D
     We follow down moves with down moves, over moves with down moves, and
     up moves with either up or over moves. (Other choices are just backing
     up and hence useless.)

     As with NNIs, we keep track of up-profiles as we go. However, some of the regular
     profiles may also become "stale" so it is a bit trickier.

     We store the traversal before we do SPRs to avoid any possible infinite loop
  */
  double last_tot_len = 0.0;
  if (NJ->nSeq <= 3 || maxSPRLength < 1)
    return;
  if (slow)
    last_tot_len = TreeLength(NJ, /*recomputeLengths*/true);
  int *nodeList = mymalloc(sizeof(int) * NJ->maxnodes);
  int nodeListLen = 0;
  traversal_t traversal = InitTraversal(NJ);
  int node = NJ->root;
  while((node = TraversePostorder(node, NJ, /*IN/OUT*/traversal)) >= 0) {
    nodeList[nodeListLen++] = node;
  }
  assert(nodeListLen == NJ->maxnode);
  traversal = FreeTraversal(traversal,NJ);

  profile_t **upProfiles = UpProfiles(NJ);
  spr_step_t *steps = mymalloc(sizeof(spr_step_t) * maxSPRLength); /* current chain of SPRs */

  int i;
  for (i = 0; i < nodeListLen; i++) {
    node = nodeList[i];
    if (node == NJ->root)
      continue; /* nothing to do for root */
    if ((i % 100) == 0)
      ProgressReport("SPR round %3d of %3d, %d of %d splits",
		     iRound+1, nRounds, i+1, nodeListLen);
    /* The nodes to NNI around */
    int nodeAround[2] = { NJ->parent[node], Sibling(NJ, node) };
    if (NJ->parent[node] == NJ->root) {
      /* NNI around both siblings instead */
      RootSiblings(NJ, node, /*OUT*/nodeAround);
    }
    bool bChanged = false;
    int iAround;
    for (iAround = 0; iAround < 2 && bChanged == false; iAround++) {
      int ACFirst;
      for (ACFirst = 0; ACFirst < 2 && bChanged == false; ACFirst++) {
	if(verbose > 3)
	  PrintNJInternal(stderr, NJ);
	int chainLength = FindSPRSteps(/*IN/OUT*/NJ, node, nodeAround[iAround],
				       upProfiles, /*OUT*/steps, maxSPRLength, (bool)ACFirst);
	double dMinDelta = 0.0;
	int iCBest = -1;
	double dTotDelta = 0.0;
	int iC;
	for (iC = 0; iC < chainLength; iC++) {
	  dTotDelta += steps[iC].deltaLength;
	  if (dTotDelta < dMinDelta) {
	    dMinDelta = dTotDelta;
	    iCBest = iC;
	  }
	}
      
	if (verbose>1) {
	  fprintf(stderr, "SPR %s %d around %d chainLength %d of %d deltaLength %.5f swaps:",
		  iCBest >= 0 ? "move" : "abandoned",
		  node,nodeAround[iAround],iCBest+1,chainLength,dMinDelta);
	  for (iC = 0; iC < chainLength; iC++)
	    fprintf(stderr, " (%d,%d)%.4f", steps[iC].nodes[0], steps[iC].nodes[1], steps[iC].deltaLength);
	  fprintf(stderr,"\n");
	}
	for (iC = chainLength - 1; iC > iCBest; iC--)
	  UnwindSPRStep(/*IN/OUT*/NJ, /*IN*/&steps[iC], /*IN/OUT*/upProfiles);
	if(verbose > 3)
	  PrintNJInternal(stderr, NJ);
	while (slow && iCBest >= 0) {
	  double expected_tot_len = last_tot_len + dMinDelta;
	  double new_tot_len = TreeLength(NJ, /*recompute*/true);
	  if (verbose > 1)
	    fprintf(stderr, "Total branch-length is now %.4f was %.4f expected %.4f\n",
		    new_tot_len, last_tot_len, expected_tot_len);
	  if (new_tot_len < last_tot_len) {
	    last_tot_len = new_tot_len;
	    break;		/* no rewinding necessary */
	  }
	  if (verbose > 1)
	    fprintf(stderr, "Rewinding SPR to %d\n",iCBest);
	  UnwindSPRStep(/*IN/OUT*/NJ, /*IN*/&steps[iCBest], /*IN/OUT*/upProfiles);
	  dMinDelta -= steps[iCBest].deltaLength;
	  iCBest--;
	}
	if (iCBest >= 0)
	  bChanged = true;
      }	/* loop over which step to take at 1st NNI */
    } /* loop over which node to pivot around */

    if (bChanged) {
      nSPR++;		/* the SPR move is OK */
      /* make sure all the profiles are OK */
      int j;
      for (j = 0; j < NJ->maxnodes; j++)
	DeleteUpProfile(upProfiles, NJ, j);
      int ancestor;
      for (ancestor = NJ->parent[node]; ancestor >= 0; ancestor = NJ->parent[ancestor])
	RecomputeProfile(/*IN/OUT*/NJ, upProfiles, ancestor);
    }
  } /* end loop over subtrees to prune & regraft */
  steps = myfree(steps, sizeof(spr_step_t) * maxSPRLength);
  upProfiles = FreeUpProfiles(upProfiles,NJ);
  nodeList = myfree(nodeList, sizeof(int) * NJ->maxnodes);
}

void RecomputeProfile(/*IN/OUT*/NJ_t *NJ, /*IN/OUT*/profile_t **upProfiles, int node) {
  if (node < NJ->nSeq || node == NJ->root)
    return;			/* no profile to compute */
  assert(NJ->child[node].nChild==2);

  profile_t *profiles[4];
  double weight = 0.5;
  if (!bionj) {
    profiles[0] = NJ->profiles[NJ->child[node].child[0]];
    profiles[1] = NJ->profiles[NJ->child[node].child[1]];
  } else {
    int nodeABC[3];
    SetupABCD(NJ, node, /*OUT*/profiles, /*IN/OUT*/upProfiles, /*OUT*/nodeABC);
    weight = QuartetWeight(profiles, NJ->distance_matrix, NJ->nPos);
  }
  if (verbose>2)
    fprintf(stderr, "Recompute %d from %d %d weight %.3f\n",
	    node, NJ->child[node].child[0], NJ->child[node].child[1], weight);
  NJ->profiles[node] = FreeProfile(NJ->profiles[node], NJ->nPos, NJ->nConstraints);
  NJ->profiles[node] = AverageProfile(profiles[0], profiles[1],
				      NJ->nPos, NJ->nConstraints,
				      NJ->distance_matrix, weight);
}

/* The BIONJ-like formula for the weight of A when building a profile for AB is
     1/2 + (avgD(B,CD) - avgD(A,CD))/(2*d(A,B))
*/
double QuartetWeight(profile_t *profiles[4], distance_matrix_t *dmat, int nPos) {
  if (!bionj)
    return(-1.0); /* even weighting */
  double d[6];
  CorrectedPairDistances(profiles, 4, dmat, nPos, /*OUT*/d);
  if (d[qAB] < 0.01)
    return -1.0;
  double weight = 0.5 + ((d[qBC]+d[qBD])-(d[qAC]+d[qAD]))/(4*d[qAB]);
  if (weight < 0)
    weight = 0;
  if (weight > 1)
    weight = 1;
  return (weight);
}

/* Resets the children entry of parent and also the parent entry of newchild */
void ReplaceChild(/*IN/OUT*/NJ_t *NJ, int parent, int oldchild, int newchild) {
  NJ->parent[newchild] = parent;

  int iChild;
  for (iChild = 0; iChild < NJ->child[parent].nChild; iChild++) {
    if (NJ->child[parent].child[iChild] == oldchild) {
      NJ->child[parent].child[iChild] = newchild;
      return;
    }
  }
  assert(0);
}

/* Recomputes all branch lengths

   For internal branches such as (A,B) vs. (C,D), uses the formula 

   length(AB|CD) = (d(A,C)+d(A,D)+d(B,C)+d(B,D))/4 - d(A,B)/2 - d(C,D)/2

   (where all distances are profile distances - diameters).

   For external branches (e.g. to leaves) A vs. (B,C), use the formula

   length(A|BC) = (d(A,B)+d(A,C)-d(B,C))/2
*/
void UpdateBranchLengths(/*IN/OUT*/NJ_t *NJ) {
  if (NJ->nSeq < 2)
    return;
  else if (NJ->nSeq == 2) {
    int root = NJ->root;
    int nodeA = NJ->child[root].child[0];
    int nodeB = NJ->child[root].child[1];
    besthit_t h;
    ProfileDist(NJ->profiles[nodeA],NJ->profiles[nodeB],
		NJ->nPos, NJ->distance_matrix, /*OUT*/&h);
    if (logdist)
      h.dist = LogCorrect(h.dist);
    NJ->branchlength[nodeA] = h.dist/2.0;
    NJ->branchlength[nodeB] = h.dist/2.0;
    return;
  }

  profile_t **upProfiles = UpProfiles(NJ);
  traversal_t traversal = InitTraversal(NJ);
  int node = NJ->root;

  while((node = TraversePostorder(node, NJ, /*IN/OUT*/traversal)) >= 0) {
    /* reset branch length of node (distance to its parent) */
    if (node == NJ->root)
      continue; /* no branch length to set */
    if (node < NJ->nSeq) { /* a leaf */
      profile_t *profileA = NJ->profiles[node];
      profile_t *profileB = NULL;
      profile_t *profileC = NULL;

      int sib = Sibling(NJ,node);
      if (sib == -1) { /* at root, have 2 siblings */
	int sibs[2];
	RootSiblings(NJ, node, /*OUT*/sibs);
	profileB = NJ->profiles[sibs[0]];
	profileC = NJ->profiles[sibs[1]];
      } else {
	profileB = NJ->profiles[sib];
	profileC = GetUpProfile(/*IN/OUT*/upProfiles, NJ, NJ->parent[node]);
      }
      profile_t *profiles[3] = {profileA,profileB,profileC};
      double d[3]; /*AB,AC,BC*/
      CorrectedPairDistances(profiles, 3, NJ->distance_matrix, NJ->nPos, /*OUT*/d);
      /* d(A,BC) = (dAB+dAC-dBC)/2 */
      NJ->branchlength[node] = (d[0]+d[1]-d[2])/2.0;
    } else {
      profile_t *profiles[4];
      int nodeABC[3];
      SetupABCD(NJ, node, /*OUT*/profiles, /*IN/OUT*/upProfiles, /*OUT*/nodeABC);
      double d[6];
      CorrectedPairDistances(profiles, 4, NJ->distance_matrix, NJ->nPos, /*OUT*/d);
      NJ->branchlength[node] = (d[qAC]+d[qAD]+d[qBC]+d[qBD])/4.0 - (d[qAB]+d[qCD])/2.0;
      
      /* no longer needed */
      DeleteUpProfile(upProfiles, NJ, nodeABC[0]);
      DeleteUpProfile(upProfiles, NJ, nodeABC[1]);
    }
  }
  traversal = FreeTraversal(traversal,NJ);
  upProfiles = FreeUpProfiles(upProfiles,NJ);
}

void ReliabilityNJ(/*IN/OUT*/NJ_t *NJ) {
  /* For each non-root node N, with children A,B, parent P, sibling C, and grandparent G,
     we test the reliability of the split (A,B) versus rest by comparing the profiles
     of A, B, C, and the "up-profile" of P.

     Each node's upProfile is the average of its sibling's (down)-profile + its parent's up-profile
     (If node's parent is the root, then there are two siblings and we don't need an up-profile)

     To save memory, we do depth-first-search down from the root, and we only keep
     up-profiles for nodes in the active path.
  */
  int i;

  if (NJ->nSeq <= 3 || nBootstrap <= 0)
    return;			/* nothing to do */

  /* Pick columns for resampling, stored as col[iBoot*nPos + j] */
  int *col = (int*)mymalloc(sizeof(int)*NJ->nPos*nBootstrap);
  for (i = 0; i < nBootstrap; i++) {
    int j;
    for (j = 0; j < NJ->nPos; j++) {
      int pos   = (int)(knuth_rand() * NJ->nPos);
      if (pos<0)
	pos = 0;
      else if (pos == NJ->nPos)
	pos = NJ->nPos-1;
      col[i*NJ->nPos + j] = pos;
    }
  }
  if (verbose > 5) {
    for (i=0; i < 3 && i < nBootstrap; i++) {
      fprintf(stderr,"Boot%d",i);
      int j;
      for (j = 0; j < NJ->nPos; j++) {
	fprintf(stderr,"\t%d",col[i*NJ->nPos+j]);
      }
      fprintf(stderr,"\n");
    }
  }

  profile_t **upProfiles = UpProfiles(NJ);
  traversal_t traversal = InitTraversal(NJ);
  int node = NJ->root;
  int iNodesDone = 0;
  while((node = TraversePostorder(node, NJ, /*IN/OUT*/traversal)) >= 0) {
    if (node < NJ->nSeq || node == NJ->root)
      continue; /* nothing to do for leaves or root */

    if(iNodesDone > 0 && (iNodesDone % 100) == 0)
      ProgressReport("Local bootstrap for %6d of %6d internal splits", iNodesDone, NJ->nSeq-3, 0, 0);
    iNodesDone++;

    profile_t *profiles[4];
    int nodeABC[3];
    SetupABCD(NJ, node, /*OUT*/profiles, /*IN/OUT*/upProfiles, /*OUT*/nodeABC);

    NJ->support[node] = SplitSupport(profiles[0], profiles[1], profiles[2], profiles[3],
				     NJ->distance_matrix,
				     NJ->nPos,
				     col);

    /* no longer needed */
    DeleteUpProfile(upProfiles, NJ, nodeABC[0]);
    DeleteUpProfile(upProfiles, NJ, nodeABC[1]);
  }
  traversal = FreeTraversal(traversal,NJ);
  upProfiles = FreeUpProfiles(upProfiles,NJ);
  col = myfree(col, sizeof(int)*NJ->nPos+nBootstrap);
}

profile_t *NewProfile(int nPos, int nConstraints) {
  profile_t *profile = (profile_t *)mymalloc(sizeof(profile_t));
  profile->weights = mymalloc(sizeof(float)*nPos);
  profile->codes = mymalloc(sizeof(unsigned char)*nPos);
  profile->vectors = NULL;
  profile->nVectors = 0;
  profile->codeDist = NULL;
  if (nConstraints == 0) {
    profile->nOn = NULL;
    profile->nOff = NULL;
  } else {
    profile->nOn = mymalloc(sizeof(int)*nConstraints);
    profile->nOff = mymalloc(sizeof(int)*nConstraints);
  }
  return(profile);
}

profile_t *FreeProfile(profile_t *profile, int nPos, int nConstraints) {
    if(profile==NULL) return(NULL);
    myfree(profile->codes, nPos);
    myfree(profile->weights, nPos);
    myfree(profile->vectors, sizeof(float)*nCodes*profile->nVectors);
    myfree(profile->codeDist, sizeof(float)*nCodes*nPos);
    if (nConstraints > 0) {
      myfree(profile->nOn, sizeof(int)*nConstraints);
      myfree(profile->nOff,  sizeof(int)*nConstraints);
    }
    return(myfree(profile, sizeof(profile_t)));
}

void SetupABCD(NJ_t *NJ, int node,
	       /* the 4 profiles; the last one is an outprofile */
	       /*OUT*/profile_t *profiles[4], 
	       /*IN/OUT*/profile_t **upProfiles,
	       /*OUT*/int nodeABC[3]) {
  int parent = NJ->parent[node];
  assert(parent >= 0);
  assert(NJ->child[node].nChild == 2);
  nodeABC[0] = NJ->child[node].child[0]; /*A*/
  nodeABC[1] = NJ->child[node].child[1]; /*B*/

  profile_t *profile4 = NULL;
  if (parent == NJ->root) {
    int sibs[2];
    RootSiblings(NJ, node, /*OUT*/sibs);
    nodeABC[2] = sibs[0];
    profile4 = NJ->profiles[sibs[1]];
  } else {
    nodeABC[2] = Sibling(NJ,node);
    assert(nodeABC[2] >= 0);
    profile4 = GetUpProfile(upProfiles,NJ,parent);
  }
  int i;
  for (i = 0; i < 3; i++)
    profiles[i] = NJ->profiles[nodeABC[i]];
  profiles[3] = profile4;
}


int Sibling(NJ_t *NJ, int node) {
  int parent = NJ->parent[node];
  if (parent < 0 || parent == NJ->root)
    return(-1);
  int iChild;
  for(iChild=0;iChild<NJ->child[parent].nChild;iChild++) {
    if(NJ->child[parent].child[iChild] != node)
      return (NJ->child[parent].child[iChild]);
  }
  assert(0);
  return(-1);
}

void RootSiblings(NJ_t *NJ, int node, /*OUT*/int sibs[2]) {
  assert(NJ->parent[node] == NJ->root);
  assert(NJ->child[NJ->root].nChild == 3);

  int nSibs = 0;
  int iChild;
  for(iChild=0; iChild < NJ->child[NJ->root].nChild; iChild++) {
    int child = NJ->child[NJ->root].child[iChild];
    if (child != node) sibs[nSibs++] = child;
  }
  assert(nSibs==2);
}

void TestSplits(NJ_t *NJ, /*OUT*/SplitCount_t *splitcount) {
  const double tolerance = 1e-6;
  splitcount->nBadSplits = 0;
  splitcount->nConstraintViolations = 0;
  splitcount->nBadBoth = 0;
  splitcount->nSplits = 0;

  profile_t **upProfiles = UpProfiles(NJ);
  traversal_t traversal = InitTraversal(NJ);
  int node = NJ->root;

  while((node = TraversePostorder(node, NJ, /*IN/OUT*/traversal)) >= 0) {
    if (node < NJ->nSeq || node == NJ->root)
      continue; /* nothing to do for leaves or root */

    profile_t *profiles[4];
    int nodeABC[3];
    SetupABCD(NJ, node, /*OUT*/profiles, /*IN/OUT*/upProfiles, /*OUT*/nodeABC);

    if (verbose>2)
      fprintf(stderr,"Testing Split around %d: A=%d B=%d C=%d D=out(C) (node parent %d)\n",
	      node, nodeABC[0], nodeABC[1], nodeABC[2], NJ->parent[node]);

    double d[6];		/* distances, perhaps log-corrected distances, no constraint penalties */
    CorrectedPairDistances(profiles, 4, NJ->distance_matrix, NJ->nPos, /*OUT*/d);

    /* alignment-based scores for each split (lower is better) */
    double sABvsCD = d[qAB] + d[qCD];
    double sACvsBD = d[qAC] + d[qBD];
    double sADvsBC = d[qAD] + d[qBC];

    /* constraint penalties, indexed by nni_t */
    double p[3];
    QuartetConstraintPenalties(profiles, NJ->nConstraints, /*OUT*/p);

    int nConstraintsViolated = 0;
    int iC;
    for (iC=0; iC < NJ->nConstraints; iC++) {
      if (SplitViolatesConstraint(profiles, iC)) {
	nConstraintsViolated++;
	if (verbose > 1) {
	  double penalty[3] = {0.0,0.0,0.0};
	  (void)QuartetConstraintPenaltiesPiece(profiles, iC, /*OUT*/penalty);
	  fprintf(stderr, "Violate constraint %d at %d (children %d %d) penalties %.3f %.3f %.3f %d/%d %d/%d %d/%d %d/%d\n",
		  iC, node, NJ->child[node].child[0], NJ->child[node].child[1],
		  penalty[ABvsCD], penalty[ACvsBD], penalty[ADvsBC],
		  profiles[0]->nOn[iC], profiles[0]->nOff[iC],
		  profiles[1]->nOn[iC], profiles[1]->nOff[iC],
		  profiles[2]->nOn[iC], profiles[2]->nOff[iC],
		  profiles[3]->nOn[iC], profiles[3]->nOff[iC]);
	}
      }
    }

    bool bBadDist = sABvsCD > sACvsBD + tolerance || sABvsCD > sADvsBC + tolerance;
    bool bBadConstr = p[ABvsCD] > p[ACvsBD] + tolerance || p[ABvsCD] > p[ADvsBC] + tolerance;

    splitcount->nSplits++;
    if (nConstraintsViolated > 0)
      splitcount->nConstraintViolations++; /* count splits with any violations, not #constraints in a splits */
    if (bBadDist)
      splitcount->nBadSplits++;
    if (bBadDist && bBadConstr)
      splitcount->nBadBoth++;
    if (bBadConstr && verbose > 2) {
      /* Which NNI would be better */
      double dist_advantage = 0;
      double constraint_penalty = 0;
      if (p[ACvsBD] < p[ADvsBC]) {
	dist_advantage = sACvsBD - sABvsCD;
	constraint_penalty = p[ABvsCD] - p[ACvsBD];
      } else {
	dist_advantage = sADvsBC - sABvsCD;
	constraint_penalty = p[ABvsCD] - p[ADvsBC];
      }
      fprintf(stderr, "Violate constraints %d distance_advantage %.3f constraint_penalty %.3f (children %d %d):",
	      node, dist_advantage, constraint_penalty,
	      NJ->child[node].child[0], NJ->child[node].child[1]);
      /* list the constraints with a penalty, meaning that ABCD all have non-zero
         values and that AB|CD worse than others */
      for (iC = 0; iC < NJ->nConstraints; iC++) {
	double ppart[6];
	if (QuartetConstraintPenaltiesPiece(profiles, iC, /*OUT*/ppart)) {
	  if (ppart[qAB] + ppart[qCD] > ppart[qAD] + ppart[qBC] + tolerance
	      || ppart[qAB] + ppart[qCD] > ppart[qAC] + ppart[qBD] + tolerance)
	    fprintf(stderr, " %d (%d/%d %d/%d %d/%d %d/%d)", iC,
		    profiles[0]->nOn[iC], profiles[0]->nOff[iC],
		    profiles[1]->nOn[iC], profiles[1]->nOff[iC],
		    profiles[2]->nOn[iC], profiles[2]->nOff[iC],
		    profiles[3]->nOn[iC], profiles[3]->nOff[iC]);
	}
      }
      fprintf(stderr, "\n");
    }
    
    /* no longer needed */
    DeleteUpProfile(upProfiles, NJ, nodeABC[0]);
    DeleteUpProfile(upProfiles, NJ, nodeABC[1]);
  }
  traversal = FreeTraversal(traversal,NJ);
  upProfiles = FreeUpProfiles(upProfiles,NJ);
}

/* Computes support for (A,B),(C,D) compared to that for (A,C),(B,D) and (A,D),(B,C) */
double SplitSupport(profile_t *pA, profile_t *pB, profile_t *pC, profile_t *pD,
		    /*OPTIONAL*/distance_matrix_t *dmat,
		    int nPos,
		    int *col) {
  int i,j;

  /* Note distpieces are weighted */
  double *distpieces[6];
  double *weights[6];
  for (j = 0; j < 6; j++) {
    distpieces[j] = (double*)mymalloc(sizeof(double)*nPos);
    weights[j] = (double*)mymalloc(sizeof(double)*nPos);
  }

  int iFreqA = 0;
  int iFreqB = 0;
  int iFreqC = 0;
  int iFreqD = 0;
  for (i = 0; i < nPos; i++) {
    float *fA = GET_FREQ(pA, i, /*IN/OUT*/iFreqA);
    float *fB = GET_FREQ(pB, i, /*IN/OUT*/iFreqB);
    float *fC = GET_FREQ(pC, i, /*IN/OUT*/iFreqC);
    float *fD = GET_FREQ(pD, i, /*IN/OUT*/iFreqD);

    weights[qAB][i] = pA->weights[i] * pB->weights[i];
    weights[qAC][i] = pA->weights[i] * pC->weights[i];
    weights[qAD][i] = pA->weights[i] * pD->weights[i];
    weights[qBC][i] = pB->weights[i] * pC->weights[i];
    weights[qBD][i] = pB->weights[i] * pD->weights[i];
    weights[qCD][i] = pC->weights[i] * pD->weights[i];

    distpieces[qAB][i] = weights[qAB][i] * ProfileDistPiece(pA->codes[i], pB->codes[i], fA, fB, dmat, NULL);
    distpieces[qAC][i] = weights[qAC][i] * ProfileDistPiece(pA->codes[i], pC->codes[i], fA, fC, dmat, NULL);
    distpieces[qAD][i] = weights[qAD][i] * ProfileDistPiece(pA->codes[i], pD->codes[i], fA, fD, dmat, NULL);
    distpieces[qBC][i] = weights[qBC][i] * ProfileDistPiece(pB->codes[i], pC->codes[i], fB, fC, dmat, NULL);
    distpieces[qBD][i] = weights[qBD][i] * ProfileDistPiece(pB->codes[i], pD->codes[i], fB, fD, dmat, NULL);
    distpieces[qCD][i] = weights[qCD][i] * ProfileDistPiece(pC->codes[i], pD->codes[i], fC, fD, dmat, NULL);
  }
  assert(iFreqA == pA->nVectors);
  assert(iFreqB == pB->nVectors);
  assert(iFreqC == pC->nVectors);
  assert(iFreqD == pD->nVectors);

  double totpieces[6];
  double totweights[6];
  double dists[6];
  for (j = 0; j < 6; j++) {
    totpieces[j] = 0.0;
    totweights[j] = 0.0;
    for (i = 0; i < nPos; i++) {
      totpieces[j] += distpieces[j][i];
      totweights[j] += weights[j][i];
    }
    dists[j] = totweights[j] > 0.01 ? totpieces[j]/totweights[j] : 3.0;
    if (logdist)
      dists[j] = LogCorrect(dists[j]);
  }

  /* Support1 = Support(AB|CD over AC|BD) = d(A,C)+d(B,D)-d(A,B)-d(C,D)
     Support2 = Support(AB|CD over AD|BC) = d(A,D)+d(B,C)-d(A,B)-d(C,D)
  */
  double support1 = dists[qAC] + dists[qBD] - dists[qAB] - dists[qCD];
  double support2 = dists[qAD] + dists[qBC] - dists[qAB] - dists[qCD];

  if (support1 < 0 || support2 < 0) {
    nSuboptimalSplits++;	/* Another split seems superior */
  }

  assert(nBootstrap > 0);
  int nSupport = 0;

  int iBoot;
  for (iBoot=0;iBoot<nBootstrap;iBoot++) {
    int *colw = &col[nPos*iBoot];

    for (j = 0; j < 6; j++) {
      double totp = 0;
      double totw = 0;
      double *d = distpieces[j];
      double *w = weights[j];
      for (i=0; i<nPos; i++) {
	int c = colw[i];
	totp += d[c];
	totw += w[c];
      }
      dists[j] = totw > 0.01 ? totp/totw : 3.0;
      if (logdist)
	dists[j] = LogCorrect(dists[j]);
    }
    support1 = dists[qAC] + dists[qBD] - dists[qAB] - dists[qCD];
    support2 = dists[qAD] + dists[qBC] - dists[qAB] - dists[qCD];
    if (support1 > 0 && support2 > 0)
      nSupport++;
  } /* end loop over bootstrap replicates */

  for (j = 0; j < 6; j++) {
    distpieces[j] = myfree(distpieces[j], sizeof(double)*nPos);
    weights[j] = myfree(weights[j], sizeof(double)*nPos);
  }
  return( nSupport/(double)nBootstrap );
}

void SetDistCriterion(/*IN/OUT*/NJ_t *NJ, int nActive, /*IN/OUT*/besthit_t *hit) {
  if (hit->i < NJ->nSeq && hit->j < NJ->nSeq) {
    SeqDist(NJ->profiles[hit->i]->codes,
	    NJ->profiles[hit->j]->codes,
	    NJ->nPos, NJ->distance_matrix, /*OUT*/hit);
  } else {
    ProfileDist(NJ->profiles[hit->i],
		NJ->profiles[hit->j],
		NJ->nPos, NJ->distance_matrix, /*OUT*/hit);
    hit->dist -= (NJ->diameter[hit->i] + NJ->diameter[hit->j]);
  }
  hit->dist += constraintWeight
    * (double)JoinConstraintPenalty(NJ, hit->i, hit->j);
  SetCriterion(NJ,nActive,/*IN/OUT*/hit);
}

void SetCriterion(/*IN/OUT*/NJ_t *NJ, int nActive, /*IN/OUT*/besthit_t *join) {
  if(join->i < 0
     || join->j < 0
     || NJ->parent[join->i] >= 0
     || NJ->parent[join->j] >= 0)
    return;
  assert(NJ->nOutDistActive[join->i] >= nActive);
  assert(NJ->nOutDistActive[join->j] >= nActive);

  int nDiffAllow = tophitsMult > 0 ? (int)(nActive*staleOutLimit) : 0;
  if (NJ->nOutDistActive[join->i] - nActive > nDiffAllow)
    SetOutDistance(NJ, join->i, nActive);
  if (NJ->nOutDistActive[join->j] - nActive > nDiffAllow)
    SetOutDistance(NJ, join->j, nActive);
  double outI = NJ->outDistances[join->i];
  if (NJ->nOutDistActive[join->i] != nActive)
    outI *= (nActive-1)/(double)(NJ->nOutDistActive[join->i]-1);
  double outJ = NJ->outDistances[join->j];
  if (NJ->nOutDistActive[join->j] != nActive)
    outJ *= (nActive-1)/(double)(NJ->nOutDistActive[join->j]-1);
  join->criterion = join->dist - (outI+outJ)/(double)(nActive-2);
  if (verbose > 2 && nActive <= 5) {
    fprintf(stderr, "Set Criterion to join %d %d with nActive=%d dist+penalty %.3f criterion %.3f\n",
	    join->i, join->j, nActive, join->dist, join->criterion);
  }
}

void SetOutDistance(NJ_t *NJ, int iNode, int nActive) {
  if (NJ->nOutDistActive[iNode] == nActive)
    return;

  /* May be called by InitNJ before we have parents */
  assert(iNode>=0 && (NJ->parent == NULL || NJ->parent[iNode]<0));
  besthit_t dist;
  ProfileDist(NJ->profiles[iNode], NJ->outprofile, NJ->nPos, NJ->distance_matrix, &dist);
  outprofileOps++;

  /* out(A) = sum(X!=A) d(A,X)
     = sum(X!=A) (profiledist(A,X) - diam(A) - diam(X))
     = sum(X!=A) profiledist(A,X) - (N-1)*diam(A) - (totdiam - diam(A))

     in the absence of gaps:
     profiledist(A,out) = mean profiledist(A, all active nodes)
     sum(X!=A) profiledist(A,X) = N * profiledist(A,out) - profiledist(A,A)

     With gaps, we need to take the weights of the comparisons into account, where
     w(Ai) is the weight of position i in profile A:
     w(A,B) = sum_i w(Ai) * w(Bi)
     d(A,B) = sum_i w(Ai) * w(Bi) * d(Ai,Bi) / w(A,B)

     sum(X!=A) profiledist(A,X) ~= (N-1) * profiledist(A, Out w/o A)
     profiledist(A, Out w/o A) = sum_X!=A sum_i d(Ai,Xi) * w(Ai) * w(Bi) / ( sum_X!=A sum_i w(Ai) * w(Bi) )
     d(A, Out) = sum_A sum_i d(Ai,Xi) * w(Ai) * w(Bi) / ( sum_X sum_i w(Ai) * w(Bi) )

     and so we get
     profiledist(A,out w/o A) = (top of d(A,Out) - top of d(A,A)) / (weight of d(A,Out) - weight of d(A,A))
     top = dist * weight
     with another correction of nActive because the weight of the out-profile is the average
     weight not the total weight.
  */
  double top = (nActive-1)
    * (dist.dist * dist.weight * nActive - NJ->selfweight[iNode] * NJ->selfdist[iNode]);
  double bottom = (dist.weight * nActive - NJ->selfweight[iNode]);
  double pdistOutWithoutA = top/bottom;
  NJ->outDistances[iNode] =  bottom > 0.01 ? 
    pdistOutWithoutA - NJ->diameter[iNode] * (nActive-1) - (NJ->totdiam - NJ->diameter[iNode])
    : 3.0;
  NJ->nOutDistActive[iNode] = nActive;

  if(verbose>3 && iNode < 5)
    fprintf(stderr,"NewOutDist for %d %f from dist %f selfd %f diam %f totdiam %f newActive %d\n",
	    iNode, NJ->outDistances[iNode], dist.dist, NJ->selfdist[iNode], NJ->diameter[iNode],
	    NJ->totdiam, nActive);
  if (verbose>6 && (iNode % 10) == 0) {
    /* Compute the actual out-distance and compare */
    double total = 0.0;
    double total_pd = 0.0;
    int j;
    for (j=0;j<NJ->maxnode;j++) {
      if (j!=iNode && (NJ->parent==NULL || NJ->parent[j]<0)) {
	besthit_t bh;
	ProfileDist(NJ->profiles[iNode], NJ->profiles[j], NJ->nPos, NJ->distance_matrix, /*OUT*/&bh);
	total_pd += bh.dist;
	total += bh.dist - (NJ->diameter[iNode] + NJ->diameter[j]);
      }
    }
    fprintf(stderr,"OutDist for Node %d %f truth %f profiled %f truth %f pd_err %f\n",
	    iNode, NJ->outDistances[iNode], total, pdistOutWithoutA, total_pd,fabs(pdistOutWithoutA-total_pd));
  }

}

/* Helper function for sorting in SetAllLeafTopHits,
   and the global variables it needs
*/
NJ_t *CompareSeedNJ = NULL;
int *CompareSeedGaps = NULL;
int CompareSeeds(const void *c1, const void *c2) {
  int seed1 = *(int *)c1;
  int seed2 = *(int *)c2;
  int gapdiff = CompareSeedGaps[seed1] - CompareSeedGaps[seed2];
  if (gapdiff != 0) return(gapdiff);	/* fewer gaps is better */
  double outdiff = CompareSeedNJ->outDistances[seed1] - CompareSeedNJ->outDistances[seed2];
  if(outdiff < 0) return(-1);	/* closer to more nodes is better */
  if(outdiff > 0) return(1);
  return(0);
}

/* Using the seed heuristic and the close global variable */
void SetAllLeafTopHits(NJ_t *NJ, int m, /*OUT*/besthit_t **tophits) {
  double close = tophitsClose;
  if (close < 0) {
    double logN = log((double)NJ->nSeq)/log(2.0);
    close = logN/(logN+2.0);
  }
  /* Sort the potential seeds, by a combination of nGaps and NJ->outDistances
     We don't store nGaps so we need to compute that
  */
  int *nGaps = (int*)mymalloc(sizeof(int)*NJ->nSeq);
  int iNode;
  for(iNode=0; iNode<NJ->nSeq; iNode++) {
    nGaps[iNode] = (int)(0.5 + NJ->nPos - NJ->selfweight[iNode]);
  }
  int *seeds = (int*)mymalloc(sizeof(int)*NJ->nSeq);
  for (iNode=0; iNode<NJ->nSeq; iNode++) seeds[iNode] = iNode;
  CompareSeedNJ = NJ;
  CompareSeedGaps = nGaps;
  qsort(/*IN/OUT*/seeds, NJ->nSeq, sizeof(int), CompareSeeds);
  CompareSeedNJ = NULL;
  CompareSeedGaps = NULL;

  besthit_t *besthitsSeed = (besthit_t*)mymalloc(sizeof(besthit_t)*NJ->nSeq);
  besthit_t *besthitsNeighbor = (besthit_t*)mymalloc(sizeof(besthit_t)*2*m);
  besthit_t bestjoin;

  /* For each seed, save its top 2*m hits and then look for close neighbors */
  assert(2*m <= NJ->nSeq);
  int iSeed;
  for(iSeed=0; iSeed < NJ->nSeq; iSeed++) {
    int seed = seeds[iSeed];
    if (iSeed > 0 && (iSeed % 100) == 0)
      ProgressReport("Top hits for %6d of %6d unique sequences", iSeed, NJ->nSeq, 0, 0);
    if (tophits[seed] != NULL) {
      if(verbose>2) fprintf(stderr, "Skipping seed %d\n", seed);
      continue;
    }
    if(verbose>2) fprintf(stderr,"Trying seed %d\n", seed);
    SetBestHit(seed, NJ, /*nActive*/NJ->nSeq, /*OUT*/&bestjoin, /*OUT*/besthitsSeed);

    /* sort & save top hits of self. besthitsSeed is now sorted. */
    tophits[seed] = SortSaveBestHits(besthitsSeed, seed, /*IN-SIZE*/NJ->nSeq, /*OUT-SIZE*/m);

    /* find "close" neighbors and compute their top hits */
    double neardist = besthitsSeed[2*m-1].dist * close;
    /* must have at least average weight, rem higher is better
       and allow a bit more than average, e.g. if we are looking for within 30% away,
       20% more gaps than usual seems OK
       Alternatively, have a coverage requirement in case neighbor is short
    */
    double nearweight = 0;
    int iClose;
    for (iClose = 0; iClose < 2*m; iClose++)
      nearweight += besthitsSeed[iClose].weight;
    nearweight = nearweight/(2.0*m); /* average */
    nearweight *= (1.0-2.0*neardist/3.0);
    double nearcover = 1.0 - neardist/2.0;

    if(verbose>2) fprintf(stderr,"Distance limit for close neighbors %f weight %f ungapped %d\n",
			  neardist, nearweight, NJ->nPos-nGaps[seed]);
    for (iClose = 0; iClose < m; iClose++) {
      besthit_t *closehit = &tophits[seed][iClose];
      int closeNode = closehit->j;
      /* If within close-distance, or identical, use as close neighbor */
      bool close = closehit->dist <= neardist
	&& (closehit->weight >= nearweight
	    || closehit->weight >= (NJ->nPos-nGaps[closeNode])*nearcover);
      bool identical = closehit->dist == 0
	&& fabs(closehit->weight - (NJ->nPos - nGaps[seed])) < 1e-5
	&& fabs(closehit->weight - (NJ->nPos - nGaps[closeNode])) < 1e-5;
      if (tophits[closeNode] == NULL && (close || identical)) {
	nCloseUsed++;
	if(verbose>2) fprintf(stderr, "Near neighbor %d (rank %d weight %f ungapped %d %d)\n",
			      closeNode, iClose, tophits[seed][iClose].weight,
			      NJ->nPos-nGaps[seed],
			      NJ->nPos-nGaps[closeNode]);

	/* compute top 2*m hits */
	TransferBestHits(NJ, /*nActive*/NJ->nSeq,
			 closeNode,
			 /*IN*/besthitsSeed, /*SIZE*/2*m,
			 /*OUT*/besthitsNeighbor,
			 /*updateDistance*/true);
	tophits[closeNode] = SortSaveBestHits(besthitsNeighbor, closeNode, /*IN-SIZE*/2*m, /*OUT-SIZE*/m);
	if (verbose>3 && (closeNode%10)==0) {
	  /* Double-check the top-hit list */
	  besthit_t best;
	  SetBestHit(closeNode, NJ, /*nActive*/NJ->nSeq, &best, /*OPTIONAL-ALL*/NULL);
	  int iBest;
	  int found = 0;
	  for (iBest=0; iBest<2*m; iBest++) {
	    if (tophits[closeNode][iBest].j == best.j) {
	      found = 1;
	      break;
	    }
	  }
	  if (found==0) fprintf(stderr,"Missed from %d to %d %f %f gaps %d %d seedgap %d\n",
				best.i,best.j,best.dist,best.criterion,
				nGaps[best.i],nGaps[best.j],nGaps[seed]);
	} /* end double-checking test of closeNode */
      }	/* end test if should transfer hits */
    } /* end loop over close candidates */
  } /* end loop over seeds */

  for (iNode=0;iNode<NJ->nSeq;iNode++) {
    assert(tophits[iNode] != NULL);
    assert(tophits[iNode][0].i == iNode);
    assert(tophits[iNode][0].j >= 0);
    assert(tophits[iNode][0].j < NJ->nSeq);
    assert(tophits[iNode][0].j != iNode);
  }
  if (verbose>1) fprintf(stderr, "#Close neighbors among leaves: %ld seeds %ld\n", nCloseUsed, NJ->nSeq-nCloseUsed);
  nGaps = myfree(nGaps, sizeof(int)*NJ->nSeq);
  seeds = myfree(seeds, sizeof(int)*NJ->nSeq);
  besthitsSeed = myfree(besthitsSeed, sizeof(besthit_t)*NJ->nSeq);
  besthitsNeighbor = myfree(besthitsNeighbor, sizeof(besthit_t)*2*m);
}

/* Updates out-distances but does not reset or update visible set */
int GetBestFromTopHits(int iNode,
			/*IN/OUT*/NJ_t *NJ,
			int nActive,
			/*IN/UPDATE*/besthit_t *tophits,
			int nTopHits) {
  assert(NJ->parent[iNode] < 0);
  int bestIndex = -1;
  if(!fastest)
    SetOutDistance(NJ, iNode, nActive); /* ensure out-distances are not stale */

  int iBest;
  for(iBest=0; iBest<nTopHits; iBest++) {
    besthit_t *hit = &tophits[iBest];
    if(hit->j < 0) continue;	/* empty slot */
    assert(hit->i == iNode);

    /* Walk up to active node and compute new distance value if necessary */
    int j = hit->j;
    while(NJ->parent[j] >= 0) j = NJ->parent[j];
    if (iNode == j)
      continue;
    if(!fastest)
      SetOutDistance(NJ, j, nActive); /* ensure out-distances are not stale */
    if (j != hit->j) {
      hit->j = j;
      SetDistCriterion(NJ, nActive, /*IN/OUT*/hit);
    } else {
      /* Update out distances if needed, and compute criterion */
      SetCriterion(/*IN/OUT*/NJ, nActive, /*IN/OUT*/hit);
    }
    if (bestIndex < 0)
      bestIndex = iBest;
    else if (hit->criterion < tophits[bestIndex].criterion)
      bestIndex = iBest;
  }
  assert(bestIndex >= 0);	/* a hit was found */
  assert(tophits[bestIndex].i == iNode);
  if (verbose > 5) fprintf(stderr, "BestHit %d %d %f %f\n",
			   tophits[bestIndex].i, tophits[bestIndex].j,
			   tophits[bestIndex].dist, tophits[bestIndex].criterion);
  return(bestIndex);
}

/* Make a shorter list with only unique entries
   Also removes "stale" hits to nodes that have parents
*/
besthit_t *UniqueBestHits(NJ_t *NJ, int iNode,
			  besthit_t *combined, int nCombined,
			  /*OUT*/int *nUniqueOut) {
  qsort(/*IN/OUT*/combined, nCombined, sizeof(besthit_t), CompareHitsByJ);

  besthit_t *uniqueList = (besthit_t*)mymalloc(sizeof(besthit_t)*nCombined);
  int nUnique = 0;
  int iHit = 0;
  for (iHit = 0; iHit < nCombined; iHit++) {
    besthit_t *hit = &combined[iHit];
    if(hit->j < 0 || hit->j == iNode || NJ->parent[hit->j] >= 0) continue;
    assert(hit->i == iNode);
    if (nUnique > 0 && hit->j == uniqueList[nUnique-1].j) continue;
    assert(nUnique < nCombined);
    uniqueList[nUnique++] = *hit;
  }
  *nUniqueOut = nUnique;
  return(uniqueList);
}


/*
  Create a top hit list for the new node, either
  from children (if there are enough best hits left) or by a "refresh"
  Also set visible set for newnode
  Also update visible set for other nodes if we stumble across a "better" hit
*/
 
void TopHitJoin(/*IN/OUT*/NJ_t *NJ,
		int newnode,
		int nActive,
		int m,
		/*IN/OUT*/besthit_t **tophits,
		/*IN/OUT*/int *tophitAge,
		/*IN/OUT*/besthit_t *visible,
		/*IN/OUT*/int *topvisible) {
  besthit_t *combinedList = (besthit_t*)mymalloc(sizeof(besthit_t)*2*m);
  assert(NJ->child[newnode].nChild == 2);
  assert(tophits[newnode] == NULL);

  /* Copy the hits */
  TransferBestHits(NJ, nActive, newnode, tophits[NJ->child[newnode].child[0]], m,
		   /*OUT*/combinedList,
		   /*updateDistance*/false);
  TransferBestHits(NJ, nActive, newnode, tophits[NJ->child[newnode].child[1]], m,
		   /*OUT*/combinedList+m,
		   /*updateDistance*/false);
  int nUnique;
  besthit_t *uniqueList = UniqueBestHits(NJ, newnode, combinedList, 2*m, /*OUT*/&nUnique);
  combinedList = myfree(combinedList, sizeof(besthit_t)*2*m);

  /* Forget the top-hit lists of the joined nodes */
  int c;
  for(c = 0; c < NJ->child[newnode].nChild; c++) {
    int child = NJ->child[newnode].child[c];
    tophits[child] = myfree(tophits[child], m*sizeof(besthit_t));
  }

  tophitAge[newnode] = tophitAge[NJ->child[newnode].child[0]];
  if (tophitAge[newnode] < tophitAge[NJ->child[newnode].child[1]])
    tophitAge[newnode] = tophitAge[NJ->child[newnode].child[1]];
  tophitAge[newnode]++;

  /* If top hit ages always match, then log2(m) would mean a refresh after
     m joins, which is about what we want.
  */
  int tophitAgeLimit = (int)(0.5 + log((double)m)/log(2.0));
  if (tophitAgeLimit < 1) tophitAgeLimit = 1;
  tophitAgeLimit++;		/* make it a bit conservative, we have the tophitsRefresh threshold too */

  /* Either use the merged list as candidate top hits or do a refresh
     UniqueBestHits eliminates hits to self, so if nUnique==nActive-1,
     we've already done the exhaustive search.

     Either way, we set tophits, visible(newnode), update visible of its top hits,
     and modify topvisible: if we do a refresh, then we reset it, otherwise we update
  */
  if (nUnique==nActive-1
      || (nUnique >= (int)(0.5+m*tophitsRefresh)
	  && tophitAge[newnode] <= tophitAgeLimit)) {
    if(verbose>2) fprintf(stderr,"Top hits for %d from combined %d nActive=%d tophitsage %d\n",
			  newnode,nUnique,nActive,tophitAge[newnode]);
    /* Update distances */
    int iHit;
    for (iHit = 0; iHit < nUnique; iHit++)
      SetDistCriterion(NJ, nActive, /*IN/OUT*/&uniqueList[iHit]);
    tophits[newnode] = SortSaveBestHits(uniqueList, newnode, /*nIn*/nUnique, /*nOut*/m);
    visible[newnode] = tophits[newnode][0];
    ResetVisible(NJ, nActive, tophits[newnode], m, /*IN/OUT*/visible, topvisible);
    UpdateTopVisible(NJ, visible, /*IN/UPDATE*/topvisible, m, newnode);
  } else {
    /* need to refresh: set top hits for node and for its top hits */
    if(verbose>1) fprintf(stderr,"Top hits for %d by refresh (%d unique age %d) nActive=%d\n",
			  newnode,nUnique,tophitAge[newnode],nActive);
    nRefreshTopHits++;
    tophitAge[newnode] = 0;

    int iNode;
    /* update all out-distances */
    for (iNode = 0; iNode < NJ->maxnode; iNode++) {
      if (NJ->parent[iNode] < 0)
	SetOutDistance(/*IN/OUT*/NJ, iNode, nActive);
    }

    /* exhaustively get the best 2*m hits for newnode, set visible, and save the top m */
    besthit_t *allhits = (besthit_t*)mymalloc(sizeof(besthit_t)*NJ->maxnode);
    assert(2*m <= NJ->maxnode);
    SetBestHit(newnode, NJ, nActive, /*OUT*/&visible[newnode], /*OUT*/allhits);
    qsort(/*IN/OUT*/allhits, NJ->maxnode, sizeof(besthit_t), CompareHitsByCriterion);
    tophits[newnode] = SortSaveBestHits(allhits, newnode, /*nIn*/NJ->maxnode, /*nOut*/m);

    /* Do not need to call ResetVisible becaues we will reset visible from its top hits, below */

    /* And use the top 2*m entries to expand other best-hit lists, but only for top m */
    besthit_t *bothList = (besthit_t*)mymalloc(sizeof(besthit_t)*3*m);
    int iHit;
    for (iHit=0; iHit < m; iHit++) {
      if (allhits[iHit].i < 0) continue;
      int iNode = allhits[iHit].j;
      assert(iNode>=0);
      if (NJ->parent[iNode] >= 0) continue;
      tophitAge[iNode] = 0;

      /* Merge */
      int i;
      for (i=0;i<m;i++) {
	bothList[i] = tophits[iNode][i];
	SetCriterion(/*IN/OUT*/NJ, nActive, /*IN/OUT*/&bothList[i]);
      }
      TransferBestHits(NJ, nActive, iNode, /*IN*/allhits, /*nOldHits*/2*m,
		       /*OUT*/&bothList[m],
		       /*updateDist*/true);
      int nUnique2;
      besthit_t *uniqueList2 = UniqueBestHits(NJ, iNode, bothList, 3*m, /*OUT*/&nUnique2);
      tophits[iNode] = myfree(tophits[iNode], m*sizeof(besthit_t));
      tophits[iNode] = SortSaveBestHits(uniqueList2, iNode, /*nIn*/nUnique2, /*nOut*/m);
      visible[iNode] = tophits[iNode][0]; /* will update topvisible below */
      uniqueList2 = myfree(uniqueList2, 3*m*sizeof(besthit_t));
    }
    ResetTopVisible(NJ, nActive, visible, /*OUT*/topvisible, m);
    bothList = myfree(bothList,3*m*sizeof(besthit_t));
    allhits = myfree(allhits,sizeof(besthit_t)*NJ->maxnode);
  }
  uniqueList = myfree(uniqueList, 2*m*sizeof(besthit_t));
}

void ResetVisible(NJ_t *NJ, int nActive,
		   /*IN*/besthit_t *tophits,
		   int nTopHits,
		  /*IN/UPDATE*/besthit_t *visible,
		  /*OPTIONAL IN/UPDATE*/int *topvisible) {
  int iHit;

  /* reset visible set for all top hits of node */
  for(iHit = 0; iHit < nTopHits; iHit++) {
    besthit_t *hit = &tophits[iHit];
    if (hit->i < 0) continue;
    assert(hit->j >= 0 && NJ->parent[hit->j] < 0);
    if (NJ->parent[visible[hit->j].j] >= 0) {
      /* Visible no longer active, so use this ("reset") */
      visible[hit->j] = *hit;
      visible[hit->j].j = visible[hit->j].i;
      visible[hit->j].i = hit->j;
      if (topvisible != NULL)
	UpdateTopVisible(NJ, visible, /*IN/UPDATE*/topvisible, nTopHits, hit->j);
      if(verbose>5) fprintf(stderr,"NewVisible %d %d %f %f\n",
			    hit->j,visible[hit->j].j,visible[hit->j].dist,
			    visible[hit->j].criterion);
    } else {
      /* see if this is a better hit -- if it is, "reset" */
      SetCriterion(/*IN/OUT*/NJ, nActive, /*IN/OUT*/&visible[hit->j]); /* may change out-dists */
      if (hit->criterion < visible[hit->j].criterion) {
	visible[hit->j] = *hit;
	visible[hit->j].j = visible[hit->j].i;
	visible[hit->j].i = hit->j;
	if (topvisible != NULL)
	  UpdateTopVisible(NJ, visible, /*IN/UPDATE*/topvisible, nTopHits, hit->j);
	if(verbose>5) fprintf(stderr,"ResetVisible %d %d %f %f\n",
			      hit->j,visible[hit->j].j,visible[hit->j].dist,
			      visible[hit->j].criterion);
	nVisibleReset++;
      }
    }
  } /* end loop over hits */
}

/* Update the top-visible list to perhaps include visible[iNode] */
void UpdateTopVisible(/*IN*/NJ_t * NJ,
		      /*IN*/besthit_t *visible,
		      /*IN/UPDATE*/int *topvisible,
		      int m, 	/* length of topvisible */
		      int iNewNode) {
  assert(topvisible != NULL);

  /* First, if the list is not full, put it in somewhere;
     otherwise, remember the worst hit
  */
  int iPosWorst = -1;
  double dCriterionWorst = -1e20;
  int i;
  for (i = 0; i < m; i++) {
    int iNode = topvisible[i];
    if (iNode < 0 || NJ->parent[iNode] >= 0) {
      topvisible[i] = iNewNode;
      return;
    } else if (visible[iNode].criterion >= dCriterionWorst) {
      iPosWorst = i;
    }
  }
  /* Otherwise, may replace worst element in the list */
  if (visible[iNewNode].criterion < dCriterionWorst) {
    assert(iPosWorst >= 0);
    topvisible[iPosWorst] = iNewNode;
  }
}

/* Recompute the topvisible list */
void ResetTopVisible(/*IN*/NJ_t *NJ,
		     int nActive,
		     /*IN*/besthit_t *visible,
		     /*OUT*/int *topvisible,
		     int m) {	/* length of topvisible */
  besthit_t *visibleSorted = mymalloc(sizeof(besthit_t)*nActive);
  int nVisible = 0;		/* #entries in visibleSorted */
  int iNode;
  for (iNode = 0; iNode < NJ->maxnode; iNode++) {
    /* skip joins involving stale nodes or joins we've already saved */
    if (NJ->parent[iNode] >= 0) continue;
    assert(visible[iNode].i == iNode);
    int j = visible[iNode].j;
    assert(j >= 0);
    while (NJ->parent[j] >= 0) {
      j = NJ->parent[j];
    }
    if (j < iNode && visible[j].j == iNode) continue;	/* keep only 1 direction of the hit */
    if (j != visible[iNode].j) {
      /* we moved up teh list */
      visible[iNode].j = j;
      SetDistCriterion(NJ, nActive, /*IN/OUT*/&visible[iNode]);
    } else {
      /* just make sure out-distances are up to date */
      SetCriterion(/*IN/OUT*/NJ, nActive, &visible[iNode]);
    }
    assert(nVisible < nActive);
    visibleSorted[nVisible++] = visible[iNode];
  }
  assert(nVisible > 0);
    
  qsort(/*IN/OUT*/visibleSorted,nVisible,sizeof(besthit_t),CompareHitsByCriterion);
    
  /* Only keep the top m items, and make sure their out-distances are up to date */
  if (verbose > 2)
    fprintf(stderr, "top-hit search: nActive %d nVisible %d considering up to %d items\n",
	    nActive, nVisible, m);
  if(nVisible > m) nVisible = m;
    
  int iBest;
  for (iBest = 0; iBest < nVisible; iBest++)
    SetCriterion(/*IN/OUT*/NJ, nActive, /*IN/OUT*/&visibleSorted[iBest]);
    
  qsort(/*IN/OUT*/visibleSorted,nVisible,sizeof(besthit_t),CompareHitsByCriterion);

  /* save the sorted indices in topvisible */
  int i;
  for (i = 0; i < nVisible; i++)
    topvisible[i] = visibleSorted[i].i;
  while(i < m)
    topvisible[i++] = -1;
  myfree(visibleSorted, sizeof(besthit_t)*nActive);
}

/*
  Find best hit to do in O(N*log(N) + m*L*log(N)) time, by
  copying and sorting the visible list
  updating out-distances for the top (up to m) candidates
  selecting the best hit
  if !fastest then
  	local hill-climbing for a better join,
	using best-hit lists only, and updating
	all out-distances in every best-hit list
*/
void TopHitNJSearch(/*IN/OUT*/NJ_t *NJ, int nActive, int m,
		    /*IN/OUT*/besthit_t *visible,
		    /*IN/OUT*/int *topvisible,
		    /*IN/OUT*/int *topvisibleAge,
		    /*IN/OUT*/besthit_t **tophits,
		    /*OUT*/besthit_t *join) {
  /* first, do we have at least m/2 candidates in topvisible?
     And remember the best one */
  int nCandidate = 0;
  int iNodeBestCandidate = -1;
  double dBestCriterion = 1e20;
  int i;
  for (i = 0; i < m; i++) {
    if (topvisible[i] >= 0 && NJ->parent[topvisible[i]] < 0) {
      nCandidate++;
      int iNode = topvisible[i];
      assert(visible[iNode].i == iNode);
      int j = visible[iNode].j;
      assert(j >= 0);
      while (NJ->parent[j] >= 0)
	j = NJ->parent[j];
      if (j != visible[iNode].j) {
	visible[iNode].j = j;
	SetDistCriterion(NJ, nActive, /*IN/OUT*/&visible[iNode]);
      } else {
	/* just make sure out-distances are up to date */
	SetCriterion(/*IN/OUT*/NJ, nActive, &visible[iNode]);
      }
      if (iNodeBestCandidate < 0 || visible[iNode].criterion < dBestCriterion) {
	iNodeBestCandidate = iNode;
	dBestCriterion = visible[iNode].criterion;
      }
    }
  }
  (*topvisibleAge)++;
  if (2 * (*topvisibleAge) > m ||  (2*nCandidate < m && 2*nCandidate < nActive)) {
    /* recompute top visible */
    if (verbose > 1)
      fprintf(stderr, "Resetting the top-visible list at nActive=%d\n",nActive);
    ResetTopVisible(NJ, nActive, visible, /*OUT*/topvisible, m);
    iNodeBestCandidate = topvisible[0];
    *topvisibleAge = 0;
  } else {
    if (verbose > 1)
      fprintf(stderr, "Top-visible list size %d (nActive %d m %d)\n",
	      nCandidate, nActive, m);
  }
  assert(iNodeBestCandidate >= 0 && NJ->parent[iNodeBestCandidate] < 0);
  *join = visible[iNodeBestCandidate];
  assert(join->j >= 0 && NJ->parent[join->j] < 0);

  if(fastest)
    return;

  int changed;
  do {
    changed = 0;

    besthit_t *bestI = &tophits[join->i][GetBestFromTopHits(join->i, NJ, nActive, tophits[join->i], m)];
    assert(bestI->i == join->i);
    if (bestI->j != join->j && bestI->criterion < join->criterion) {
      changed = 1;
      if (verbose>2)
	fprintf(stderr,"BetterI\t%d\t%d\t%d\t%d\t%f\t%f\n",
		join->i,join->j,bestI->i,bestI->j,
		join->criterion,bestI->criterion);
      *join = *bestI;
    }

    besthit_t *bestJ = &tophits[join->j][GetBestFromTopHits(join->j, NJ, nActive, tophits[join->j], m)];
    assert(bestJ->i == join->j);
    if (bestJ->j != join->i && bestJ->criterion < join->criterion) {
      changed = 1;
      if (verbose>2)
	fprintf(stderr,"BetterJ\t%d\t%d\t%d\t%d\t%f\t%f\n",
		join->i,join->j,bestJ->i,bestJ->j,
		join->criterion,bestJ->criterion);
      *join = *bestJ;
    }
    if(changed) nBetter++;
  } while(changed);
}

int NGaps(/*IN*/NJ_t *NJ, int iNode) {
  assert(iNode < NJ->nSeq);
  int nGaps = 0;
  int p;
  for(p=0; p<NJ->nPos; p++) {
    if (NJ->profiles[iNode]->codes[p] == NOCODE)
      nGaps++;
  }
  return(nGaps);
}

int CompareHitsByCriterion(const void *c1, const void *c2) {
  const besthit_t *hit1 = (besthit_t*)c1;
  const besthit_t *hit2 = (besthit_t*)c2;
  if (hit1->criterion < hit2->criterion) return(-1);
  if (hit1->criterion > hit2->criterion) return(1);
  return(0);
}

int CompareHitsByJ(const void *c1, const void *c2) {
  const besthit_t *hit1 = (besthit_t*)c1;
  const besthit_t *hit2 = (besthit_t*)c2;
  return hit1->j - hit2->j;
}

besthit_t *SortSaveBestHits(besthit_t *besthits, int iNode, int insize, int outsize) {
  qsort(/*IN/OUT*/besthits,insize,sizeof(besthit_t),CompareHitsByCriterion);

  besthit_t *saved = (besthit_t*)mymalloc(sizeof(besthit_t)*outsize);
  int nSaved = 0;
  int iBest;
  for (iBest = 0; iBest < insize && nSaved < outsize; iBest++) {
    assert(besthits[iBest].i == iNode);
    if (besthits[iBest].j != iNode)
      saved[nSaved++] = besthits[iBest];
  }
  /* pad saved list with invalid entries if necessary */
  for(; nSaved < outsize; nSaved++) {
    saved[nSaved].i = -1;
    saved[nSaved].j = -1;
    saved[nSaved].weight = 0;
    saved[nSaved].dist = 1e20;
    saved[nSaved].criterion = 1e20;
  }
  return(saved);
}

void TransferBestHits(/*IN/OUT*/NJ_t *NJ,
		       int nActive,
		      int iNode,
		      /*IN*/besthit_t *oldhits,
		      int nOldHits,
		      /*OUT*/besthit_t *newhits,
		      bool updateDistances) {
  assert(NJ->parent[iNode] < 0);

  int iBest;
  for(iBest = 0; iBest < nOldHits; iBest++) {
    int j = oldhits[iBest].j;
    besthit_t *new = &newhits[iBest];
    if(j<0) {			/* empty (invalid) entry */
      new->i = iNode;
      new->j = -1;
      new->weight = 0;
      new->dist = 1e20;
      new->criterion = 1e20;
    } else {
      /* Move up to an active node */
      while(NJ->parent[j] >= 0)
	j = NJ->parent[j];
      
      new->i = iNode;
      new->j = j;
      if (updateDistances)
	SetDistCriterion(NJ, nActive, /*IN/OUT*/new);
    }
  }
}

/* Algorithm 26.2.17 from Abromowitz and Stegun, Handbook of Mathematical Functions
   Absolute accuracy of only about 1e-7, which is enough for us
*/
double pnorm(double x)
{
  double b1 =  0.319381530;
  double b2 = -0.356563782;
  double b3 =  1.781477937;
  double b4 = -1.821255978;
  double b5 =  1.330274429;
  double p  =  0.2316419;
  double c  =  0.39894228;

  if(x >= 0.0) {
    double t = 1.0 / ( 1.0 + p * x );
    return (1.0 - c * exp( -x * x / 2.0 ) * t *
	    ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
  }
  /*else*/
  double t = 1.0 / ( 1.0 - p * x );
  return ( c * exp( -x * x / 2.0 ) * t *
	   ( t *( t * ( t * ( t * b5 + b4 ) + b3 ) + b2 ) + b1 ));
}

void *mymalloc(size_t sz) {
  if (sz == 0) return(NULL);
  void *new = malloc(sz);
  if (new == NULL) {
    fprintf(stderr, "Out of memory\n");
    exit(1);
  }
  szAllAlloc += sz;
  mymallocUsed += sz;
#ifdef TRACK_MEMORY
  struct mallinfo mi = mallinfo();
  if (mi.arena+mi.hblkhd > maxmallocHeap)
    maxmallocHeap = mi.arena+mi.hblkhd;
#endif
  return (new);
}

void *mymemdup(void *data, size_t sz) {
  if(data==NULL) return(NULL);
  void *new = mymalloc(sz);
  memcpy(/*to*/new, /*from*/data, sz);
  return(new);
}

void *myrealloc(void *data, size_t szOld, size_t szNew) {
  if (data == NULL && szOld == 0)
    return(mymalloc(szNew));
  if (data == NULL || szOld == 0 || szNew == 0) {
    fprintf(stderr,"Empty myrealloc\n");
    exit(1);
  }
  void *new = realloc(data,szNew);
  if (new == NULL) {
    fprintf(stderr, "Out of memory\n");
    exit(1);
  }
  szAllAlloc += (szNew-szOld);
  mymallocUsed += (szNew-szOld);
#ifdef TRACK_MEMORY
    struct mallinfo mi = mallinfo();
    if (mi.arena+mi.hblkhd > maxmallocHeap)
      maxmallocHeap = mi.arena+mi.hblkhd;
#endif
  return(new);
}

void *myfree(void *p, size_t sz) {
  if(p==NULL) return(NULL);
  free(p);
  mymallocUsed -= sz;
  return(NULL);
}

/* The random number generator is taken from D E Knuth 
   http://www-cs-faculty.stanford.edu/~knuth/taocp.html
*/

/*    This program by D E Knuth is in the public domain and freely copyable.
 *    It is explained in Seminumerical Algorithms, 3rd edition, Section 3.6
 *    (or in the errata to the 2nd edition --- see
 *        http://www-cs-faculty.stanford.edu/~knuth/taocp.html
 *    in the changes to Volume 2 on pages 171 and following).              */

/*    N.B. The MODIFICATIONS introduced in the 9th printing (2002) are
      included here; there's no backwards compatibility with the original. */

/*    This version also adopts Brendan McKay's suggestion to
      accommodate naive users who forget to call ran_start(seed).          */

/*    If you find any bugs, please report them immediately to
 *                 taocp@cs.stanford.edu
 *    (and you will be rewarded if the bug is genuine). Thanks!            */

/************ see the book for explanations and caveats! *******************/
/************ in particular, you need two's complement arithmetic **********/

#define KK 100                     /* the long lag */
#define LL  37                     /* the short lag */
#define MM (1L<<30)                 /* the modulus */
#define mod_diff(x,y) (((x)-(y))&(MM-1)) /* subtraction mod MM */

long ran_x[KK];                    /* the generator state */

#ifdef __STDC__
void ran_array(long aa[],int n)
#else
     void ran_array(aa,n)    /* put n new random numbers in aa */
     long *aa;   /* destination */
     int n;      /* array length (must be at least KK) */
#endif
{
  register int i,j;
  for (j=0;j<KK;j++) aa[j]=ran_x[j];
  for (;j<n;j++) aa[j]=mod_diff(aa[j-KK],aa[j-LL]);
  for (i=0;i<LL;i++,j++) ran_x[i]=mod_diff(aa[j-KK],aa[j-LL]);
  for (;i<KK;i++,j++) ran_x[i]=mod_diff(aa[j-KK],ran_x[i-LL]);
}

/* the following routines are from exercise 3.6--15 */
/* after calling ran_start, get new randoms by, e.g., "x=ran_arr_next()" */

#define QUALITY 1009 /* recommended quality level for high-res use */
long ran_arr_buf[QUALITY];
long ran_arr_dummy=-1, ran_arr_started=-1;
long *ran_arr_ptr=&ran_arr_dummy; /* the next random number, or -1 */

#define TT  70   /* guaranteed separation between streams */
#define is_odd(x)  ((x)&1)          /* units bit of x */

#ifdef __STDC__
void ran_start(long seed)
#else
     void ran_start(seed)    /* do this before using ran_array */
     long seed;            /* selector for different streams */
#endif
{
  register int t,j;
  long x[KK+KK-1];              /* the preparation buffer */
  register long ss=(seed+2)&(MM-2);
  for (j=0;j<KK;j++) {
    x[j]=ss;                      /* bootstrap the buffer */
    ss<<=1; if (ss>=MM) ss-=MM-2; /* cyclic shift 29 bits */
  }
  x[1]++;              /* make x[1] (and only x[1]) odd */
  for (ss=seed&(MM-1),t=TT-1; t; ) {       
    for (j=KK-1;j>0;j--) x[j+j]=x[j], x[j+j-1]=0; /* "square" */
    for (j=KK+KK-2;j>=KK;j--)
      x[j-(KK-LL)]=mod_diff(x[j-(KK-LL)],x[j]),
	x[j-KK]=mod_diff(x[j-KK],x[j]);
    if (is_odd(ss)) {              /* "multiply by z" */
      for (j=KK;j>0;j--)  x[j]=x[j-1];
      x[0]=x[KK];            /* shift the buffer cyclically */
      x[LL]=mod_diff(x[LL],x[KK]);
    }
    if (ss) ss>>=1; else t--;
  }
  for (j=0;j<LL;j++) ran_x[j+KK-LL]=x[j];
  for (;j<KK;j++) ran_x[j-LL]=x[j];
  for (j=0;j<10;j++) ran_array(x,KK+KK-1); /* warm things up */
  ran_arr_ptr=&ran_arr_started;
}

#define ran_arr_next() (*ran_arr_ptr>=0? *ran_arr_ptr++: ran_arr_cycle())
long ran_arr_cycle()
{
  if (ran_arr_ptr==&ran_arr_dummy)
    ran_start(314159L); /* the user forgot to initialize */
  ran_array(ran_arr_buf,QUALITY);
  ran_arr_buf[KK]=-1;
  ran_arr_ptr=ran_arr_buf+1;
  return ran_arr_buf[0];
}

/* end of code from Knuth */

double knuth_rand() {
  return(9.31322574615479e-10 * ran_arr_next()); /* multiply by 2**-30 */
}

hashstrings_t *MakeHashtable(char **strings, int nStrings) {
  hashstrings_t *hash = (hashstrings_t*)mymalloc(sizeof(hashstrings_t));
  hash->nBuckets = 8*nStrings;
  hash->buckets = (hashbucket_t*)mymalloc(sizeof(hashbucket_t) * hash->nBuckets);
  int i;
  for (i=0; i < hash->nBuckets; i++) {
    hash->buckets[i].string = NULL;
    hash->buckets[i].nCount = 0;
    hash->buckets[i].first = -1;
  }
  for (i=0; i < nStrings; i++) {
    hashiterator_t hi = FindMatch(hash, strings[i]);
    if (hash->buckets[hi].string == NULL) {
      /* save a unique entry */
      assert(hash->buckets[hi].nCount == 0);
      hash->buckets[hi].string = strings[i];
      hash->buckets[hi].nCount = 1;
      hash->buckets[hi].first = i;
    } else {
      /* record a duplicate entry */
      assert(hash->buckets[hi].string != NULL);
      assert(strcmp(hash->buckets[hi].string, strings[i]) == 0);
      assert(hash->buckets[hi].first >= 0);
      hash->buckets[hi].nCount++;
    }
  }
  return(hash);
}

hashstrings_t *FreeHashtable(hashstrings_t* hash) {
  if (hash != NULL) {
    myfree(hash->buckets, sizeof(hashbucket_t) * hash->nBuckets);
    myfree(hash, sizeof(hashstrings_t));
  }
  return(NULL);
}

#define MAXADLER 65521
hashiterator_t FindMatch(hashstrings_t *hash, char *string) {
  /* Adler-32 checksum */
  unsigned int hashA = 1;
  unsigned int hashB = 0;
  char *p;
  for (p = string; *p != '\0'; p++) {
    hashA = ((unsigned int)*p + hashA);
    hashB = hashA+hashB;
  }
  hashA %= MAXADLER;
  hashB %= MAXADLER;
  hashiterator_t hi = (hashB*65536+hashA) % hash->nBuckets;
  while(hash->buckets[hi].string != NULL
	&& strcmp(hash->buckets[hi].string, string) != 0) {
    hi++;
    if (hi >= hash->nBuckets)
      hi = 0;
  }
  return(hi);
}

char *GetHashString(hashstrings_t *hash, hashiterator_t hi) {
  return(hash->buckets[hi].string);
}

int HashCount(hashstrings_t *hash, hashiterator_t hi) {
  return(hash->buckets[hi].nCount);
}

int HashFirst(hashstrings_t *hash, hashiterator_t hi) {
  return(hash->buckets[hi].first);
}

uniquify_t *UniquifyAln(alignment_t *aln) {
    int nUniqueSeq = 0;
    char **uniqueSeq = (char**)mymalloc(aln->nSeq * sizeof(char*)); /* iUnique -> seq */
    int *uniqueFirst = (int*)mymalloc(aln->nSeq * sizeof(int)); /* iUnique -> iFirst in aln */
    int *alnNext = (int*)mymalloc(aln->nSeq * sizeof(int)); /* i in aln -> next, or -1 */
    int *alnToUniq = (int*)mymalloc(aln->nSeq * sizeof(int)); /* i in aln -> iUnique; many -> -1 */

    int i;
    for (i = 0; i < aln->nSeq; i++) {
      uniqueSeq[i] = NULL;
      uniqueFirst[i] = -1;
      alnNext[i] = -1;
      alnToUniq[i] = -1;
    }
    hashstrings_t *hashseqs = MakeHashtable(aln->seqs, aln->nSeq);
    for (i=0; i<aln->nSeq; i++) {
      hashiterator_t hi = FindMatch(hashseqs,aln->seqs[i]);
      int first = HashFirst(hashseqs,hi);
      if (first == i) {
	uniqueSeq[nUniqueSeq] = aln->seqs[i];
	uniqueFirst[nUniqueSeq] = i;
	alnToUniq[i] = nUniqueSeq;
	nUniqueSeq++;
      } else {
	int last = first;
	while (alnNext[last] != -1)
	  last = alnNext[last];
	assert(last>=0);
	alnNext[last] = i;
	assert(alnToUniq[last] >= 0 && alnToUniq[last] < nUniqueSeq);
	alnToUniq[i] = alnToUniq[last];
      }
    }
    assert(nUniqueSeq>0);
    hashseqs = FreeHashtable(hashseqs);

    uniquify_t *uniquify = (uniquify_t*)mymalloc(sizeof(uniquify_t));
    uniquify->nSeq = aln->nSeq;
    uniquify->nUnique = nUniqueSeq;
    uniquify->uniqueFirst = uniqueFirst;
    uniquify->alnNext = alnNext;
    uniquify->alnToUniq = alnToUniq;
    uniquify->uniqueSeq = uniqueSeq;
    return(uniquify);
}

uniquify_t *FreeUniquify(uniquify_t *unique) {
  if (unique != NULL) {
    myfree(unique->uniqueFirst, sizeof(int)*unique->nSeq);
    myfree(unique->alnNext, sizeof(int)*unique->nSeq);
    myfree(unique->alnToUniq, sizeof(int)*unique->nSeq);
    myfree(unique->uniqueSeq, sizeof(char*)*unique->nSeq);
    myfree(unique,sizeof(uniquify_t));
    unique = NULL;
  }
  return(unique);
}

traversal_t InitTraversal(NJ_t *NJ) {
  traversal_t worked = (bool*)mymalloc(sizeof(bool)*NJ->maxnodes);
  int i;
  for (i=0; i<NJ->maxnodes; i++)
    worked[i] = false;
  return(worked);
}

int TraversePostorder(int node, NJ_t *NJ, /*IN/OUT*/traversal_t traversal) {
  while(1) {
    assert(node >= 0);

    /* move to a child if possible */
    bool found = false;
    int iChild;
    for (iChild=0; iChild < NJ->child[node].nChild; iChild++) {
      int child = NJ->child[node].child[iChild];
      if (!traversal[child]) {
	node = child;
	found = true;
	break;
      }
    }
    if (found)
      continue; /* keep moving down */
    if (!traversal[node]) {
      traversal[node] = true;
      return(node);
    }
    /* If we've already done this node, need to move up */
    if (node == NJ->root)
      return(-1); /* nowhere to go -- done traversing */
    node = NJ->parent[node];
  }
}

traversal_t FreeTraversal(traversal_t traversal, NJ_t *NJ) {
  myfree(traversal, sizeof(bool)*NJ->maxnodes);
  return(NULL);
}

profile_t **UpProfiles(NJ_t *NJ) {
  profile_t **upProfiles = (profile_t**)mymalloc(sizeof(profile_t*)*NJ->maxnodes);
  int i;
  for (i=0; i<NJ->maxnodes; i++) upProfiles[i] = NULL;
  return(upProfiles);
}

profile_t *GetUpProfile(/*IN/OUT*/profile_t **upProfiles, NJ_t *NJ, int outnode) {
  assert(outnode != NJ->root && outnode >= NJ->nSeq); /* not for root or leaves */
  if (upProfiles[outnode] != NULL)
    return(upProfiles[outnode]);

  int depth;
  int *pathToRoot = PathToRoot(NJ, outnode, /*OUT*/&depth);
  int i;
  /* depth-1 is root */
  for (i = depth-2; i>=0; i--) {
    int node = pathToRoot[i];

    if (upProfiles[node] == NULL) {
      /* Note -- SetupABCD may call GetUpProfile, but it should do it farther
	 up in the path to the root
      */
      profile_t *profiles[4];
      int nodeABC[3];
      SetupABCD(NJ, node, /*OUT*/profiles, /*IN/OUT*/upProfiles, /*OUT*/nodeABC);
      profile_t *profilesCDAB[4] = { profiles[2], profiles[3], profiles[0], profiles[1] };
      double weight = QuartetWeight(profilesCDAB, NJ->distance_matrix, NJ->nPos);
      if (verbose>3)
	fprintf(stderr, "Compute upprofile of %d from %d and parents (vs. children %d %d) with weight %.3f\n",
		node, nodeABC[2], nodeABC[0], nodeABC[1], weight);
      upProfiles[node] = AverageProfile(profiles[2], profiles[3],
					NJ->nPos, NJ->nConstraints,
					NJ->distance_matrix,
					weight);
    }
  }
  FreePath(pathToRoot,NJ);
  assert(upProfiles[outnode] != NULL);
  return(upProfiles[outnode]);
}

profile_t *DeleteUpProfile(/*IN/OUT*/profile_t **upProfiles, NJ_t *NJ, int node) {
  assert(node>=0 && node < NJ->maxnodes);
  if (upProfiles[node] != NULL)
    upProfiles[node] = FreeProfile(upProfiles[node], NJ->nPos, NJ->nConstraints); /* returns NULL */
  return(NULL);
}

profile_t **FreeUpProfiles(profile_t **upProfiles, NJ_t *NJ) {
  int i;
  for (i=0; i < NJ->maxnodes; i++)
    DeleteUpProfile(upProfiles, NJ, i);
  myfree(upProfiles, sizeof(profile_t*)*NJ->maxnodes);
  return(NULL);
}

int *PathToRoot(NJ_t *NJ, int node, /*OUT*/int *outDepth) {
  int *pathToRoot = (int*)mymalloc(sizeof(int)*NJ->maxnodes);
  int depth = 0;
  int ancestor = node;
  while(ancestor >= 0) {
    pathToRoot[depth] = ancestor;
    ancestor = NJ->parent[ancestor];
    depth++;
  }
  *outDepth = depth;
  return(pathToRoot);
}

int *FreePath(int *path, NJ_t *NJ) {
  myfree(path, sizeof(int)*NJ->maxnodes);
  return(NULL);
}


distance_matrix_t matrixBLOSUM45 =
  {
    /*distances*/
    { 
      {0, 1.31097856157468, 1.06573001937323, 1.2682782988532, 0.90471293383305, 1.05855446876905, 1.05232790675508, 0.769574440593014, 1.27579668305679, 0.964604099952603, 0.987178199640556, 1.05007594438157, 1.05464162250736, 1.1985987403937, 0.967404475245526, 0.700490199584332, 0.880060189098976, 1.09748548316685, 1.28141710375267, 0.800038509951648},
      {1.31097856157468, 0, 0.8010890222701, 0.953340718498495, 1.36011107208122, 0.631543775840481, 0.791014908659279, 1.15694899265629, 0.761152570032029, 1.45014917711188, 1.17792001455227, 0.394661075648738, 0.998807558909651, 1.135143404599, 1.15432562628921, 1.05309036790541, 1.05010474413616, 1.03938321130789, 0.963216908696184, 1.20274751778601},
      {1.06573001937323, 0.8010890222701, 0, 0.488217214273568, 1.10567116937273, 0.814970207038261, 0.810176440932339, 0.746487413974582, 0.61876156253224, 1.17886558630004, 1.52003670190022, 0.808442678243754, 1.2889025816028, 1.16264109995678, 1.18228799147301, 0.679475681649858, 0.853658619686283, 1.68988558988005, 1.24297493464833, 1.55207513886163},
      {1.2682782988532, 0.953340718498495, 0.488217214273568, 0, 1.31581050011876, 0.769778474953791, 0.482077627352988, 0.888361752320536, 0.736360849050364, 1.76756333403346, 1.43574761894039, 0.763612910719347, 1.53386612356483, 1.74323672079854, 0.886347403928663, 0.808614044804528, 1.01590147813779, 1.59617804551619, 1.1740494822217, 1.46600946033173},
      {0.90471293383305, 1.36011107208122, 1.10567116937273, 1.31581050011876, 0, 1.3836789310481, 1.37553994252576, 1.26740695314856, 1.32361065635259, 1.26087264215993, 1.02417540515351, 1.37259631233791, 1.09416720447891, 0.986982088723923, 1.59321190226694, 0.915638787768407, 0.913042853922533, 1.80744143643002, 1.3294417177004, 0.830022143283238},
      {1.05855446876905, 0.631543775840481, 0.814970207038261, 0.769778474953791, 1.3836789310481, 0, 0.506942797642807, 1.17699648087288, 0.614595446514896, 1.17092829494457, 1.19833088638994, 0.637341078675405, 0.806490842729072, 1.83315144709714, 0.932064479113502, 0.850321696813199, 1.06830084665916, 1.05739353225849, 0.979907428113788, 1.5416250309563},
      {1.05232790675508, 0.791014908659279, 0.810176440932339, 0.482077627352988, 1.37553994252576, 0.506942797642807, 0, 1.17007322676118, 0.769786956320484, 1.46659942462342, 1.19128214039009, 0.633592151371708, 1.27269395724349, 1.44641491621774, 0.735428579892476, 0.845319988414402, 1.06201695511881, 1.324395996498, 1.22734387448031, 1.53255698189437},
      {0.769574440593014, 1.15694899265629, 0.746487413974582, 0.888361752320536, 1.26740695314856, 1.17699648087288, 1.17007322676118, 0, 1.1259007054424, 1.7025415585924, 1.38293205218175, 1.16756929156758, 1.17264582493965, 1.33271035269688, 1.07564768421292, 0.778868281341681, 1.23287107008366, 0.968539655354582, 1.42479529031801, 1.41208067821187},
      {1.27579668305679, 0.761152570032029, 0.61876156253224, 0.736360849050364, 1.32361065635259, 0.614595446514896, 0.769786956320484, 1.1259007054424, 0, 1.4112324673522, 1.14630894167097, 0.967795284542623, 0.771479459384692, 1.10468029976148, 1.12334774065132, 1.02482926701639, 1.28754326478771, 1.27439749294131, 0.468683841672724, 1.47469999960758},
      {0.964604099952603, 1.45014917711188, 1.17886558630004, 1.76756333403346, 1.26087264215993, 1.17092829494457, 1.46659942462342, 1.7025415585924, 1.4112324673522, 0, 0.433350517223017, 1.463460928818, 0.462965544381851, 0.66291968000662, 1.07010201755441, 1.23000200130049, 0.973485453109068, 0.963546200571036, 0.708724769805536, 0.351200119909572},
      {0.987178199640556, 1.17792001455227, 1.52003670190022, 1.43574761894039, 1.02417540515351, 1.19833088638994, 1.19128214039009, 1.38293205218175, 1.14630894167097, 0.433350517223017, 0, 1.49770950074319, 0.473800072611076, 0.538473125003292, 1.37979627224964, 1.5859723170438, 0.996267398224516, 0.986095542821092, 0.725310666139274, 0.570542199221932},
      {1.05007594438157, 0.394661075648738, 0.808442678243754, 0.763612910719347, 1.37259631233791, 0.637341078675405, 0.633592151371708, 1.16756929156758, 0.967795284542623, 1.463460928818, 1.49770950074319, 0, 1.0079761868248, 1.44331961488922, 0.924599080166146, 1.06275728888356, 1.05974425835993, 1.04892430642749, 0.972058829603409, 1.21378822764856},
      {1.05464162250736, 0.998807558909651, 1.2889025816028, 1.53386612356483, 1.09416720447891, 0.806490842729072, 1.27269395724349, 1.17264582493965, 0.771479459384692, 0.462965544381851, 0.473800072611076, 1.0079761868248, 0, 0.72479754849538, 1.1699868662153, 1.34481214251794, 1.06435197383538, 1.05348497728858, 0.774878150710318, 0.609532859331199},
      {1.1985987403937, 1.135143404599, 1.16264109995678, 1.74323672079854, 0.986982088723923, 1.83315144709714, 1.44641491621774, 1.33271035269688, 1.10468029976148, 0.66291968000662, 0.538473125003292, 1.44331961488922, 0.72479754849538, 0, 1.32968844979665, 1.21307373491949, 0.960087571600877, 0.475142555482979, 0.349485367759138, 0.692733248746636},
      {0.967404475245526, 1.15432562628921, 1.18228799147301, 0.886347403928663, 1.59321190226694, 0.932064479113502, 0.735428579892476, 1.07564768421292, 1.12334774065132, 1.07010201755441, 1.37979627224964, 0.924599080166146, 1.1699868662153, 1.32968844979665, 0, 0.979087429691819, 0.97631161216338, 1.21751652292503, 1.42156458605332, 1.40887880416009},
      {0.700490199584332, 1.05309036790541, 0.679475681649858, 0.808614044804528, 0.915638787768407, 0.850321696813199, 0.845319988414402, 0.778868281341681, 1.02482926701639, 1.23000200130049, 1.5859723170438, 1.06275728888356, 1.34481214251794, 1.21307373491949, 0.979087429691819, 0, 0.56109848274013, 1.76318885009194, 1.29689226231656, 1.02015839286433},
      {0.880060189098976, 1.05010474413616, 0.853658619686283, 1.01590147813779, 0.913042853922533, 1.06830084665916, 1.06201695511881, 1.23287107008366, 1.28754326478771, 0.973485453109068, 0.996267398224516, 1.05974425835993, 1.06435197383538, 0.960087571600877, 0.97631161216338, 0.56109848274013, 0, 1.39547634461879, 1.02642577026706, 0.807404666228614},
      {1.09748548316685, 1.03938321130789, 1.68988558988005, 1.59617804551619, 1.80744143643002, 1.05739353225849, 1.324395996498, 0.968539655354582, 1.27439749294131, 0.963546200571036, 0.986095542821092, 1.04892430642749, 1.05348497728858, 0.475142555482979, 1.21751652292503, 1.76318885009194, 1.39547634461879, 0, 0.320002937404137, 1.268589159299},
      {1.28141710375267, 0.963216908696184, 1.24297493464833, 1.1740494822217, 1.3294417177004, 0.979907428113788, 1.22734387448031, 1.42479529031801, 0.468683841672724, 0.708724769805536, 0.725310666139274, 0.972058829603409, 0.774878150710318, 0.349485367759138, 1.42156458605332, 1.29689226231656, 1.02642577026706, 0.320002937404137, 0, 0.933095433689795},
      {0.800038509951648, 1.20274751778601, 1.55207513886163, 1.46600946033173, 0.830022143283238, 1.5416250309563, 1.53255698189437, 1.41208067821187, 1.47469999960758, 0.351200119909572, 0.570542199221932, 1.21378822764856, 0.609532859331199, 0.692733248746636, 1.40887880416009, 1.02015839286433, 0.807404666228614, 1.268589159299, 0.933095433689795, 0}
    },
    /*eigeninv*/
    {
      {-0.216311217101265, -0.215171653035930, -0.217000020881064, -0.232890860601250, -0.25403526530177, -0.211569372858927, -0.218073620637049, -0.240585637190076, -0.214507049619293, -0.228476323330312, -0.223235445346107, -0.216116483840334, -0.206903836810903, -0.223553828183343, -0.236937609127783, -0.217652789023588, -0.211982652566286, -0.245995223308316, -0.206187718714279, -0.227670670439422},
      {-0.0843931919568687, -0.0342164464991033, 0.393702284928246, -0.166018266253027, 0.0500896782860136, -0.262731388032538, 0.030139964190519, -0.253997503551094, -0.0932603349591988, -0.32884667697173, 0.199966846276877, -0.117543453869516, 0.196248237055757, -0.456448703853250, 0.139286961076387, 0.241166801918811, -0.0783508285295053, 0.377438091416498, 0.109499076984234, 0.128581669647144},
      {-0.0690428674271772, 0.0133858672878363, -0.208289917312908, 0.161232925220819, 0.0735806288007248, -0.316269599838174, -0.0640708424745702, -0.117078801507436, 0.360805085405857, 0.336899760384943, 0.0332447078185156, 0.132954055834276, 0.00595209121998118, -0.157755611190327, -0.199839273133436, 0.193688928807663, 0.0970290928040946, 0.374683975138541, -0.478110944870958, -0.243290196936098},
      {0.117284581850481, 0.310399467781876, -0.143513477698805, 0.088808130300351, 0.105747812943691, -0.373871701179853, 0.189069306295134, 0.133258225034741, -0.213043549687694, 0.301303731259140, -0.182085224761849, -0.161971915020789, 0.229301173581378, -0.293586313243755, -0.0260480060747498, -0.0217953684540699, 0.0202675755458796, -0.160134624443657, 0.431950096999465, -0.329885160320501},
      {0.256496969244703, 0.0907408349583135, 0.0135731083898029, 0.477557831930769, -0.0727379669280703, 0.101732675207959, -0.147293025369251, -0.348325291603251, -0.255678082078362, -0.187092643740172, -0.177164064346593, -0.225921480146133, 0.422318841046522, 0.319959853469398, -0.0623652546300045, 0.0824203908606883, -0.102057926881110, 0.120728407576411, -0.156845807891241, -0.123528163091204},
      {-0.00906668858975576, -0.0814722888231236, -0.0762715085459023, 0.055819989938286, -0.0540516675257271, -0.0070589302769034, -0.315813159989213, -0.0103527463419808, -0.194634331372293, -0.0185860407566822, 0.50134169352609, 0.384531812730061, -0.0405008616742061, 0.0781033650669525, 0.069334900096687, 0.396455180448549, -0.204065801866462, -0.215272089630713, 0.171046818996465, -0.396393364716348},
      {0.201971098571663, 0.489747667606921, 0.00226258734592836, 0.0969514005747054, 0.0853921636903791, 0.0862068740282345, -0.465412154271164, -0.130516676347786, 0.165513616974634, 0.0712238027886633, 0.140746943067963, -0.325919272273406, -0.421213488261598, -0.163508199065965, 0.269695802810568, -0.110296405171437, -0.106834099902202, 0.00509414588152415, 0.00909215239544615, 0.0500401865589727},
      {0.515854176692456, -0.087468413428258, 0.102796468891449, -0.06046105990993, -0.212014383772414, -0.259853648383794, -0.0997372883043333, -0.109934574535736, 0.284891018406112, -0.250578342940183, 0.142174204994568, 0.210384918947619, 0.118803190788946, -0.0268434355996836, 0.0103721198836548, -0.355555176478458, 0.428042332431476, -0.150610175411631, 0.0464090887952940, -0.140238796382057},
      {-0.239392215229762, -0.315483492656425, 0.100205194952396, 0.197830195325302, 0.40178804665223, 0.195809461460298, -0.407817115321684, 0.0226836686147386, -0.169780276210306, 0.0818161585952184, -0.172886230584939, 0.174982644851064, 0.0868786992159535, -0.198450519980824, 0.168581078329968, -0.361514336004068, 0.238668430084722, 0.165494019791904, 0.110437707249228, -0.169592003035203},
      {-0.313151735678025, 0.10757884850664, -0.49249098807229, 0.0993472335619114, -0.148695715250836, 0.0573801136941699, -0.190040373500722, 0.254848437434773, 0.134147888304352, -0.352719341442756, 0.0839609323513986, -0.207904182300122, 0.253940523323376, -0.109832138553288, 0.0980084518687944, 0.209026594443723, 0.406236051871548, -0.0521120230935943, 0.0554108014592302, 0.134681046631955},
      {-0.102905214421384, 0.235803606800009, 0.213414976431981, -0.253606415825635, 0.00945656859370683, 0.259551282655855, 0.159527348902192, 0.083218761193016, -0.286815935191867, 0.0135069477264877, 0.336758103107357, -0.271707359524149, -0.0400009875851839, 0.0871186292716414, -0.171506310409388, -0.0954276577211755, 0.393467571460712, 0.111732846649458, -0.239886066474217, -0.426474828195231},
      {-0.0130795552324104, 0.0758967690968058, -0.165099404017689, -0.46035152559912, 0.409888158016031, -0.0235053940299396, 0.0699393201709723, -0.161320910316996, 0.226111732196825, -0.177811841258496, -0.219073917645916, -0.00703219376737286, 0.162831878334912, 0.271670554900684, 0.451033612762052, 0.0820942662443393, -0.0904983490498446, -0.0587000279313978, -0.0938852980928252, -0.306078621571843},
      {0.345092040577428, -0.257721588971295, -0.301689123771848, -0.0875212184538126, 0.161012613069275, 0.385104899829821, 0.118355290985046, -0.241723794416731, 0.083201920119646, -0.0809095291508749, -0.0820275390511991, -0.115569770103317, -0.250105681098033, -0.164197583037664, -0.299481453795592, 0.255906951902366, 0.129042051416371, 0.203761730442746, 0.347550071284268, -0.109264854744020},
      {0.056345924962239, 0.072536751679082, 0.303127492633681, -0.368877185781648, -0.343024497082421, 0.206879529669083, -0.413012709639426, 0.078538816203612, 0.103382383425097, 0.288319996147499, -0.392663258459423, 0.0319588502083897, 0.220316797792669, -0.0563686494606947, -0.0869286063283735, 0.323677017794391, 0.0984875197088935, -0.0303289828821742, 0.0450197853450979, -0.0261771221270139},
      {-0.253701638374729, -0.148922815783583, 0.111794052194159, 0.157313977830326, -0.269846001260543, -0.222989872703583, 0.115441028189268, -0.350456582262355, -0.0409581422905941, 0.174078744248002, -0.130673397086811, -0.123963802708056, -0.351609207081548, 0.281548012920868, 0.340382662112428, 0.180262131025562, 0.3895263830793, 0.0121546812430960, 0.214830943227063, -0.0617782909660214},
      {-0.025854479416026, 0.480654788977767, -0.138024550829229, -0.130191670810919, 0.107816875829919, -0.111243997319276, -0.0679814460571245, -0.183167991080677, -0.363355166018786, -0.183934891092050, -0.216097125080962, 0.520240628803255, -0.179616013606479, 0.0664131536100941, -0.178350708111064, 0.0352047611606709, 0.223857228692892, 0.128363679623513, -0.000403433628490731, 0.224972110977704},
      {0.159207394033448, -0.0371517305736114, -0.294302634912281, -0.0866954375908417, -0.259998567870054, 0.284966673982689, 0.205356416771391, -0.257613708650298, -0.264820519037270, 0.293359248624603, 0.0997476397434102, 0.151390539497369, 0.165571346773648, -0.347569523551258, 0.43792310820533, -0.0723248163210163, 0.0379214984816955, -0.0542758730251438, -0.258020301801603, 0.128680501102363},
      {0.316853842351797, -0.153950010941153, -0.13387065213508, -0.0702971390607613, -0.202558481846057, -0.172941438694837, -0.068882524588574, 0.524738203063889, -0.271670479920716, -0.112864756695310, -0.146831636946145, -0.0352336188578041, -0.211108490884767, 0.097857111349555, 0.276459740956662, 0.0231297536754823, -0.0773173324868396, 0.487208384389438, -0.0734191389266824, -0.113198765573319},
      {-0.274285525741087, 0.227334266052039, -0.0973746625709059, -0.00965256583655389, -0.402438444750043, 0.198586229519026, 0.0958135064575833, -0.108934376958686, 0.253641732094319, -0.0551918478254021, 0.0243640218331436, 0.181936272247179, 0.090952738347629, 0.0603352483029044, -0.0043821671755761, -0.347720824658591, -0.267879988539971, 0.403804652116592, 0.337654323971186, -0.241509293972297},
      {-0.0197089518344238, 0.139681034626696, 0.251980475788267, 0.341846624362846, -0.075141195125153, 0.2184951591319, 0.268870823491343, 0.150392399018138, 0.134592404015057, -0.337050200539163, -0.313109373497998, 0.201993318439135, -0.217140733851970, -0.337622749083808, 0.135253284365068, 0.181729249828045, -0.00627813335422765, -0.197218833324039, -0.194060005031698, -0.303055888528004}
    },
    /*eigenval*/
    {
      20.29131, 0.5045685, 0.2769945, 0.1551147, 0.03235484, -0.04127639, -0.3516426, -0.469973, -0.5835191, -0.6913107, -0.7207972, -0.7907875, -0.9524307, -1.095310, -1.402153, -1.424179, -1.936704, -2.037965, -3.273561, -5.488734 
    },
    /*eigentot and codeFreq left out, these are initialized elsewhere*/
  };
 
