/*
 * FastTree -- neighbor joining for multiple sequence alignments using profiles
 * Morgan N. Price, January-April 2008
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
 * Use -DTRACK_MEMORY if you want it to report its memory usage
 * (but results are not correct above 4GB because mallinfo stores int values)
 *
 * To get usage guidance, do:
 * FastTree -help
 *
 * FastTree uses profiles instead of a distance matrix, and computes
 * support values for each split from the profiles of the 4 nodes
 * around the split. It stores a profile for each node and a average
 * profile over all active nodes (the "out-profile" for computing the
 * total sum of distance to other nodes).  Requires O(N*L*a) space,
 * where N is the number of sequences, L is the alignment width, and a
 * is the alphabet size
 *
 * FastTree operates only on "additive" distances -- either %different
 * or by using an amino acid similarity matrix (-matrix option) If we
 * are using %different as our distance matrix then
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
 * this can be reduced to as little as O(N**2*log(N) + N**(3/2)*L*a) time
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
 * FastTree does sort the entire visible set, however, which is
 * O(N*log(N)) work per iteration, or O(N**2*log(N)) work total.
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
 * top-hit list.
 *
 * Morgan N. Price, March 2008
 */

#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <malloc.h>

char *usage = "Usage for FastTree 0.9.1:\n"
              "FastTree [-quiet] [-boot 1000] [-slow | -fastest]\n"
              "          [-top | -notop]\n"
              "          [-topm 1.0 [-close 0.75] [-refresh 0.8]]\n"
              "          [-matrix Matrix | -nomatrix] [-nj | -bionj]\n"
              "          [-nt] [alignment] > newick_tree\n"
              "\n"
              "or\n"
              "\n"
              "FastTree [-nt] [-matrix Matrix | -nomatrix] [-logdist] -makematrix [alignment]\n"
              "    > phylip_distance_matrix\n"
              "\n"
              "  FastTree supports fasta or phylip interleaved alignments\n"
              "  By default FastTree expects protein alignments,  use -nt for nucleotides\n"
              "  FastTree reads standard input if no alignment file is given\n"
              "\n"
              "Distances:\n"
              "  By default, FastTree uses the BLOSUM45 matrix for protein sequences\n"
              "  and fraction-different as a distance for nucleotides\n"
              "  To specify a different matrix, use -matrix FilePrefix or -nomatrix\n"
              "\n"
              "Searching for the best join:\n"
              "  by default, FastTree combines the 'visible set' of fast neighbor-joining with\n"
              "      local hill-climbing as in relaxed neighbor-joining\n"
              "  -slow -- exhaustive search (standard NJ or BIONJ, but different gap handling)\n"
              "  -fastest -- search the visible set (the top hit for each node) only\n"
              "      Unlike the original fast neighbor-joining, -fastest updates visible(C)\n"
              "      after joining A and B if join(AB,C) is better than join(C,visible(C))\n"
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
              "  -bionj: weighted joins as in BIONJ (default)\n"
              "  -nj: regular (unweighted) neighbor-joining\n"
              "\n"
              "Support value options:\n"
              "  by default, FastTree computes support by a maximum-of-normals test\n"
              "  This test approximates a local form of bootstrap.\n"
              "\n"
              "   -boot 1000 -- Compute local bootstrap, where the argument is the\n"
              "          number of replicates. Reported support values still range from\n"
              "          0 to 1. Local bootstrap does not recompute the topology for each\n" 
              "          replicate, so it is very fast.\n"
              "\n";

typedef struct {
  int nPos;
  int nSeq;
  char **names;
  char **seqs;
} alignment_t;

#define MAXCODES 20
#define NOCODE 127
#define BUFFER_SIZE 1000

/* For each position in a profile, we have a weight (% non-gapped) and a
   frequency vector. (If using a matrix, the frequency vector is in eigenspace).
   We also store codes for simple profile positions (all gaps or only 1 value)
   If weight[pos] > 0 && codes[pos] == NOCODE then we store the vector
   vectors itself is sets of nCodes long, so the vector for the ith nonconstant position
   starts at &vectors[nCodes*i]
   To speed up comparison of outprofile to a sequence or other simple profile, we also
   (for outprofiles) store codeDist[iPos*nCodes+k] = dist(k,profile[iPos])
*/
typedef struct {
  float *weights;
  unsigned char *codes;
  float *vectors;		/* NULL if no non-constant positions, e.g. for leaves */
  int nVectors;
  float *codeDist;		/* Optional -- distance to each code at each position */
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

typedef struct { int nChild; int child[3]; } children_t;

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
  char **seqs;			/* a copy of the aligment sequences array */
  distance_matrix_t *distance_matrix; /* or NULL if using %identity distance */

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

/* Global variables */
/* Options */
int verbose = 1;
int slow = 0;
int fastest = 0;
int bionj = 1;
double tophitsMult = 1.0;	/* 0 means compare nodes to all other nodes */
double tophitsClose = 0.75;	/* Parameter for how close is close; also used as a coverage req. */
double tophitsRefresh = 0.8;	/* Refresh if fraction of top-hit-length drops to this */
double staleOutLimit = 0.02;	/* nActive changes by at most this amount before we recompute */
int nRecomputeOutProfile = 200;	/* Recompute out-profile every this many joins (avoid roundoff) */
int nBootstrap = 0;		/* If set, number of replicates of local bootstrap to do */
int nCodes=20;			/* 20 if protein, 4 if nucleotide */
bool useMatrix=true;

/* Performance and memory usage */
long profileOps = 0;		/* Full profile-based distance operations */
long outprofileOps = 0;		/* How many of profileOps are comparisons to outprofile */
long seqOps = 0;		/* Faster leaf-based distance operations */
long nBetter = 0;		/* Number of hill-climbing steps */
long nCloseUsed = 0;		/* Number of "close" neighbors we avoid full search for */
long nRefreshTopHits = 0;	/* Number of full-blown searches (interior nodes) */
long nVisibleReset = 0;		/* Number of resets of the visible set */
long nSuboptimalSplits = 0;	/* Number of splits that are rejected given full topology */
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
alignment_t *ReadAlignment(/*IN*/char *filename); /* Returns a list of strings (exits on failure) */
NJ_t *InitNJ(char **sequences, int nSeqs, int nPos,
	     /*IN OPTIONAL*/distance_matrix_t *); /* Allocates memory, initializes */
void FastNJ(/*IN/OUT*/NJ_t *NJ); /* Does the joins */
void ReliabilityNJ(/*IN/OUT*/NJ_t *NJ);	  /* Estimates the reliability of the joins */
void PrintNJ(FILE *, NJ_t *NJ, char **names, int *uniqueFirst, int *nameNext);

/* Use out-profile and NJ->totdiam to recompute out-distance for node iNode
   Only does this computation if the out-distance is "stale" (nOutDistActive[iNode] != nActive)
 */
void SetOutDistance(/*IN/UPDATE*/NJ_t *NJ, int iNode, int nActive);

/* Always sets join->criterion; may update NJ->outDistance and NJ->nOutDistActive,
   assumes join's weight and distance are already set
*/
void SetCriterion(/*IN/UPDATE*/NJ_t *NJ, int nActive, /*IN/OUT*/besthit_t *join);

/* Computes weight and distance and then sets the criterion (maybe update out-distances) */
void SetDistCriterion(/*IN/UPDATE*/NJ_t *NJ, int nActive, /*IN/OUT*/besthit_t *join);

/* Reports the support for the (1,2) vs. (3,4) split */
double SplitSupport(profile_t *p1, profile_t *p2, profile_t *p3, profile_t *p4,
		    /*OPTIONAL*/distance_matrix_t *dmat,
		    int nPos);

profile_t *SeqToProfile(NJ_t *NJ, char *seq, int nPos, int iNode);

/* ProfileDist and SeqDist only set the dist and weight fields
   If using an outprofile, use the second argument of ProfileDist
   for better performance.
 */
void ProfileDist(profile_t *profile1, profile_t *profile2, int nPos,
		 /*OPTIONAL*/distance_matrix_t *distance_matrix,
		 /*OUT*/besthit_t *hit);
void SeqDist(unsigned char *codes1, unsigned char *codes2, int nPos,
	     /*OPTIONAL*/distance_matrix_t *distance_matrix,
	     /*OUT*/besthit_t *hit);

/* AverageProfile is used to do a weighted combination of nodes
   when doing a join. If weight is negative, then the value is ignored and the profiles
   are averaged. The weight is *not* adjusted for the gap content of the nodes...
*/
profile_t *AverageProfile(profile_t *profile1, profile_t *profile2, int nPos,
			  distance_matrix_t *distance_matrix,
			  double weight1);

/* OutProfile does an unweighted combination of nodes to create the
   out-profile. It always sets code to NOCODE so that UpdateOutProfile
   can work.
*/
profile_t *OutProfile(profile_t **profiles, int nProfiles, int nPos,
		      distance_matrix_t *distance_matrix);

void UpdateOutProfile(/*UPDATE*/profile_t *out, profile_t *old1, profile_t *old2,
		      profile_t *new, int nActiveOld, int nPos,
		      distance_matrix_t *distance_matrix);

profile_t *NewProfile(int nPos); /* returned has no vectors */
profile_t *FreeProfile(profile_t *profile, int nPos); /* returns NULL */


/* f1 can be NULL if code1 != NOCODE, and similarly for f2
   Do not call if either weight was 0
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

/* Find the best join to do */
void TopHitNJSearch(/*IN/UPDATE*/NJ_t *NJ,
		    int nActive,
		    int m,
		    /*IN/UPDATE*/besthit_t *visible,
		    /*IN/UPDATE*/besthit_t **tophits,
		    /*OUT*/besthit_t *bestjoin);

/* Returns the index of the best hit within the tophits.
   NJ may be modified because it updates out-distances if necessary
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
		/*IN/OUT*/besthit_t *visible);

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
		  /*IN/UPDATE*/besthit_t *visible);

/* Make a shorter list with only unique entries.
   Ignores self-hits to iNode and "dead" hits to
   nodes that have parents.
*/
besthit_t *UniqueBestHits(NJ_t *NJ, int iNode, besthit_t *combined, int nCombined, /*OUT*/int *nUniqueOut);


int CompareHitsByCriterion(const void *c1, const void *c2);
int CompareHitsByJ(const void *c1, const void *c2);

int NGaps(NJ_t *NJ, int node);	/* only handles leaf sequences */

int Sibling(NJ_t *NJ, int node); /* At root, no unique sibling so returns -1 */

void *mymalloc(size_t sz);       /* Prints "Out of memory" and exits on failure */
void *myfree(void *, size_t sz); /* Always returns NULL */

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
hashstrings_t *DeleteHashtable(hashstrings_t* hash); /*returns NULL*/
hashiterator_t FindMatch(hashstrings_t *hash, char *string);

/* Return NULL if we have run out of values */
char *GetHashString(hashstrings_t *hash, hashiterator_t hi);
int HashCount(hashstrings_t *hash, hashiterator_t hi);
int HashFirst(hashstrings_t *hash, hashiterator_t hi);

/* The default amino acid distance matrix, derived from the BLOSUM45 similarity matrix */
distance_matrix_t matrixBLOSUM45;

int main(int argc, char **argv) {
  int iArg;
  char *matrixPrefix = NULL;
  distance_matrix_t *distance_matrix = NULL;
  bool make_matrix = false;
  bool logdist = false;

  for (iArg = 1; iArg < argc; iArg++) {
    if (strcmp(argv[iArg],"-makematrix") == 0) {
      make_matrix = true;
    } else if (strcmp(argv[iArg],"-logdist") == 0) {
      logdist = true;
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
    } else if (strcmp(argv[iArg], "-nt") == 0) {
      nCodes = 4;
    } else if (strcmp(argv[iArg], "-nj") == 0) {
      bionj = 0;
    } else if (strcmp(argv[iArg], "-bionj") == 0) {
      bionj = 1;
    } else if (strcmp(argv[iArg], "-boot") == 0 && iArg < argc-1) {
      iArg++;
      nBootstrap = atoi(argv[iArg]);
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
    } else if (strcmp(argv[iArg],"-help") == 0) {
      fprintf(stderr,"%s",usage);
      exit(0);
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
  if (logdist && !make_matrix) {
    fprintf(stderr,"Log-corrected distances are only supported for -makematrix not neighbor joining\n");
    exit(1);
  }

  if (verbose && !make_matrix) {		/* Report settings */
    char tophitString[100] = "no";
    if(tophitsMult>0) sprintf(tophitString,"%.2f*sqrtN close=%.2f refresh=%.2f",
			      tophitsMult, tophitsClose, tophitsRefresh);
    char supportString[100] = "max-norm";
    if (nBootstrap>0) sprintf(supportString,"Boot %d",nBootstrap);
    fprintf(stderr,"Alignment: %s\n%s distances: %s Method: %s Support: %s\nSearch: %s TopHits: %s\n",
	    fileName ? fileName : "(stdin)",
	    nCodes == 20 ? "Amino acid" : "Nucleotide",
	    matrixPrefix ? matrixPrefix : (useMatrix? "BLOSUM45 (default)" : "%different"),
	    bionj ? "BIONJ" : "neighbor-joining" ,
	    supportString,
	    slow?"Exhaustive (slow)" : (fastest ? "Fastest" : "Normal (fast/relax)"),
	    tophitString);
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

  alignment_t *aln = ReadAlignment(fileName);
  if (aln->nSeq < 1) {
    fprintf(stderr, "No alignment sequences\n");
    exit(1);
  }
  /* Check that all names are unique */
  hashstrings_t *hashnames = MakeHashtable(aln->names, aln->nSeq);
  int i;
  for (i=0; i<aln->nSeq; i++) {
    hashiterator_t hi = FindMatch(hashnames,aln->names[i]);
    if (HashCount(hashnames,hi) != 1) {
      fprintf(stderr,"Non-unique name %s in the alignment\n",aln->names[i]);
      exit(1);
    }
  }
  hashnames = DeleteHashtable(hashnames);

  /* Make a list of unique sequences -- note some lists are bigger than required */
  int nUniqueSeq = 0;
  char **uniqueSeq = (char**)mymalloc(aln->nSeq * sizeof(char*)); /* iUnique -> seq */
  int *uniqueFirst = (int*)mymalloc(aln->nSeq * sizeof(int)); /* iUnique -> iFirst in aln */
  int *nameNext = (int*)mymalloc(aln->nSeq * sizeof(int)); /* i in aln -> next, or -1 */

  for (i = 0; i < aln->nSeq; i++) {
    uniqueSeq[i] = NULL;
    uniqueFirst[i] = -1;
    nameNext[i] = -1;
  }
  hashstrings_t *hashseqs = MakeHashtable(aln->seqs, aln->nSeq);
  for (i=0; i<aln->nSeq; i++) {
    hashiterator_t hi = FindMatch(hashseqs,aln->seqs[i]);
    int first = HashFirst(hashseqs,hi);
    if (first == i) {
      uniqueSeq[nUniqueSeq] = aln->seqs[i];
      uniqueFirst[nUniqueSeq] = i;
      nUniqueSeq++;
    } else {
      int last = first;
      while (nameNext[last] != -1)
	last = nameNext[last];
      assert(last>=0);
      nameNext[last] = i;
    }
  }
  assert(nUniqueSeq>0);
  hashseqs = DeleteHashtable(hashseqs);
  
  if (verbose>1) fprintf(stderr, "read %s seqs %d (%d unique) positions %d nameLast %s seqLast %s\n",
			 fileName ? fileName : "standard input",
			 aln->nSeq, nUniqueSeq, aln->nPos, aln->names[aln->nSeq-1], aln->seqs[aln->nSeq-1]);

  clock_t clock_start = clock();
  if (make_matrix) {
    NJ_t *NJ = InitNJ(aln->seqs, aln->nSeq, aln->nPos, distance_matrix);
    printf("   %d\n",aln->nSeq);
    int i,j;
    for(i = 0; i < NJ->nSeq; i++) {
      printf("%s",aln->names[i]);
      for (j = 0; j < NJ->nSeq; j++) {
	besthit_t hit;
	SeqDist(NJ->profiles[i]->codes,NJ->profiles[j]->codes,NJ->nPos,NJ->distance_matrix,/*OUT*/&hit);
	if(logdist) {
	  hit.dist = hit.dist < 0.99 ? -1.3*log(1.0-hit.dist) : 3.0;
	  if(hit.dist > 3.0) hit.dist = 3.0;
	}
	printf(" %f", hit.dist);
      }
      printf("\n");
    }
    exit(0);
  }
  /* else */
  NJ_t *NJ = InitNJ(uniqueSeq, nUniqueSeq, aln->nPos, distance_matrix);
  FastNJ(NJ);

  /* profile-frequencies for the "up-profiles" in ReliabilityNJ take only diameter(Tree)*L*a
     space not N*L*a space, because we can free them as we go.
     And up-profile by their nature tend to be complicated.
     So save the profile-frequency counters now to exclude later results.
  */
#ifdef TRACK_MEMORY
  long svProfileFreqAlloc = nProfileFreqAlloc;
  long svProfileFreqAvoid = nProfileFreqAvoid;
#endif

  if(verbose>1) fprintf(stderr, "Topology done after %.2f sec -- computing support values\n",
			(clock()-clock_start)/(double)CLOCKS_PER_SEC);
  ReliabilityNJ(NJ);
  fflush(stderr);
  PrintNJ(stdout, NJ, aln->names, uniqueFirst, nameNext);
  fflush(stdout);
  if(verbose) {
    fprintf(stderr, "Unique sequences: %d/%d Bad splits: %ld/%d",
	    NJ->nSeq, aln->nSeq,
	    nSuboptimalSplits, NJ->nSeq >= 3 ? NJ->nSeq-3 : 0);
    if(!slow) fprintf(stderr, " Hill-climb: %ld Update-best: %ld",
		      nBetter, nVisibleReset);
    fprintf(stderr,"\n");
    if (nCloseUsed>0 || nRefreshTopHits>0)
      fprintf(stderr, "Top hits: close neighbors %ld/%d refreshes %ld\n",
	      nCloseUsed, NJ->nSeq, nRefreshTopHits);
    double dN2 = NJ->nSeq*(double)NJ->nSeq;
    fprintf(stderr, "Time %.2f Distances per N*N: by-profile %.3f (out %.3f) by-leaf %.3f\n",
	    (clock()-clock_start)/(double)CLOCKS_PER_SEC,
	    profileOps/dN2, outprofileOps/dN2, seqOps/dN2);
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
  exit(0);
}

NJ_t *InitNJ(char **sequences, int nSeq, int nPos,
	     /*OPTIONAL*/distance_matrix_t *distance_matrix) {
  int iNode;

  NJ_t *NJ = (NJ_t*)mymalloc(sizeof(NJ_t));
  NJ->root = -1; 		/* set at end of FastNJ() */
  NJ->maxnode = NJ->nSeq = nSeq;
  NJ->nPos = nPos;
  NJ->maxnodes = 2*nSeq;
  NJ->seqs = sequences;
  NJ->distance_matrix = distance_matrix;

  NJ->profiles = (profile_t **)mymalloc(sizeof(profile_t*) * NJ->maxnodes);

  for (iNode = 0; iNode < NJ->nSeq; iNode++) {
    NJ->profiles[iNode] = SeqToProfile(NJ, NJ->seqs[iNode], nPos, iNode);
  }
  if(verbose>10) fprintf(stderr,"Made sequence profiles\n");
  for (iNode = NJ->nSeq; iNode < NJ->maxnodes; iNode++) 
    NJ->profiles[iNode] = NULL; /* not yet exists */

  NJ->outprofile = OutProfile(NJ->profiles, NJ->nSeq, NJ->nPos, NJ->distance_matrix);
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
    NJ->nOutDistActive[iNode] = NJ->nSeq * 2; /* unreasonably high value */
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

  NJ->child = (children_t*)malloc(sizeof(children_t)*NJ->maxnodes);
  for (iNode= 0; iNode < NJ->maxnode; iNode++) NJ->child[iNode].nChild = 0;

  return(NJ);
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
  int nActive;
  for (nActive = NJ->nSeq; nActive > 3; nActive--) {
    besthit_t join; 		/* the join to do */
    if (slow) {
      ExhaustiveNJSearch(NJ,nActive,/*OUT*/&join);
    } else if (m>0) {
      TopHitNJSearch(/*IN/OUT*/NJ, nActive, m, /*IN/OUT*/visible, /*IN/OUT*/tophits, /*OUT*/&join);
    } else {
      FastNJSearch(NJ, nActive, /*IN/OUT*/visible, /*OUT*/&join);
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

    if (bionj && join.weight > 0 && varIJ > 0) {
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

    NJ->profiles[newnode] = AverageProfile(NJ->profiles[join.i],NJ->profiles[join.j],NJ->nPos,
					   NJ->distance_matrix,
					   bionj ? bionjWeight : /*noweight*/-1.0);

    /* Update out-distances and total diameters */
    if (((NJ->nSeq - nActive + 1) % nRecomputeOutProfile) == 0) {
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
      FreeProfile(NJ->outprofile,NJ->nPos);
      if(verbose>1) fprintf(stderr,"Recomputing outprofile\n");
      NJ->outprofile = OutProfile(activeProfiles, nSaved, NJ->nPos, NJ->distance_matrix);
      activeProfiles = myfree(activeProfiles, sizeof(profile_t*)*(nActive-1));
    } else {
      UpdateOutProfile(/*OUT*/NJ->outprofile,
		       NJ->profiles[join.i], NJ->profiles[join.j], NJ->profiles[newnode],
		       nActive, NJ->nPos,
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
		 /*IN/OUT*/visible);
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
    profile_t *out = OutProfile(p, 3, NJ->nPos, NJ->distance_matrix);
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
    FreeProfile(out,NJ->nPos);
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

void PrintNJ(FILE *fp, NJ_t *NJ, char **names,
	     int *uniqueFirst,  /* index in NJ to first index in names */
	     int *nameNext	/* index in names to next index in names or -1 */
	     ) {
  /* And print the tree: depth first search
   * The stack contains
   * list of remaining children with their depth
   * parent node, with a flag of -1 so I know to print right-paren
   */
  if (NJ->nSeq==1 && nameNext[uniqueFirst[0]] >= 0) {
    /* Special case -- otherwise we end up with double parens */
    int first = uniqueFirst[0];
    fprintf(fp,"(%s:0.0",names[first]);
    int iName = nameNext[first];
    while (iName >= 0) {
      fprintf(fp,",%s:0.0",names[iName]);
      iName = nameNext[iName];
    }
    fprintf(fp,");\n");
    return;
  }

  typedef struct { int node; int end; } stack_t;
  stack_t *stack = (stack_t *)malloc(sizeof(stack_t)*NJ->maxnodes);
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
      int first = uniqueFirst[node];
      assert(first >= 0);
      /* Print the name, or the subtree of duplicate names */
      if (nameNext[first] == -1) {
	fprintf(fp, names[uniqueFirst[node]]);
      } else {
	fprintf(fp,"(%s:0.0",names[first]);
	int iName = nameNext[first];
	while (iName >= 0) {
	  fprintf(fp,",%s:0.0",names[iName]);
	  iName = nameNext[iName];
	}
	fprintf(fp,")");
      }
      /* Print the branch length */
      fprintf(fp, ":%.5f", NJ->branchlength[node]);
    } else if (end) {
      if (node == NJ->root) fprintf(fp, ")");
      else fprintf(fp, ")%.3f:%.5f", NJ->support[node], NJ->branchlength[node]);
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

alignment_t *ReadAlignment(/*IN*/char *filename) {
  FILE *fp = filename != NULL ? fopen(filename,"r") : stdin;
  if (fp == NULL) {
    fprintf(stderr, "Cannot read %s\n", filename);
    exit(1);
  }
  int nSeq = 0;
  int nPos = 0;
  char **names = NULL;
  char **seqs = NULL;
  char buf[BUFFER_SIZE] = "";
  if (fgets(buf,sizeof(buf),fp) == NULL) {
    fprintf(stderr, "Error reading header line from %s\n",
	    filename ? filename : "standard input");
    exit(1);
  }
  if (buf[0] == '>') {
    /* FASTA, truncate names at any of these */
    char *nameStop = "(),: \t\r\n";
    char *seqSkip = " \t\r\n";
    int nSaved = 100;
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
     */
    if (sscanf(buf, "%d%d", &nSeq, &nPos) != 2
      || nSeq < 1 || nPos < 1) {
      fprintf(stderr, "Error parsing header line from %s:\n%s\n",
	      filename ? filename : "standard input",
	      buf);
      exit(1);
    }
    names = (char **)mymalloc(sizeof(char*) * nSeq);
    seqs = (char **)mymalloc(sizeof(char*) * nSeq);

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

  if (filename != NULL) {
    if (fclose(fp) != 0) {
      fprintf(stderr, "Error reading %s\n",filename);
      exit(1);
    }
  }

  alignment_t *align = (alignment_t*)mymalloc(sizeof(alignment_t));
  align->nSeq = nSeq;
  align->nPos = nPos;
  align->names = names;
  align->seqs = seqs;
  return(align);
}


profile_t *SeqToProfile(NJ_t *NJ, char *seq, int nPos, int iNode) {
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

  int seqlen = strlen(seq);
  profile_t *profile = NewProfile(nPos);

  bool bWarn = false;

  for (i = 0; i < seqlen; i++) {
    c = charToCode[(unsigned int)seq[i]];
    if(verbose>10 && i < 2) fprintf(stderr,"pos %d char %c code %d\n", i, seq[i], c);
    /* treat unknowns as gaps, but warn if verbose and unknown isn't X */
    if (c == nCodes && verbose && seq[i] != 'X') {
      if (!bWarn) {
	fprintf(stderr, "Characters in unique sequence %d replaced with gap:", iNode+1);
	bWarn = true;
      }
      fprintf(stderr, " %c%d", seq[i], i+1);
    }
    if (c == nCodes || c == NOCODE) {
      profile->codes[i] = NOCODE;
      profile->weights[i] = 0.0;
    } else {
      profile->codes[i] = c;
      profile->weights[i] = 1.0;
    }
  }
  if (bWarn)
    fprintf(stderr, "\n");
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
  seqOps++;
  hit->weight = (double)nUse;
  hit->dist = nUse > 0 ? top/(double)nUse : 1.0;
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
	assert(code1 != NOCODE);
	f1 = &dmat->codeFreq[code1][0];
      }
      if (f2 == NULL) {
	assert(code2 != NOCODE);
	f2 = &dmat->codeFreq[code2][0];
      }
      double piece = 0;
      for (k = 0; k < nCodes; k++)
	piece += f1[k] * f2[k] * dmat->eigenval[k];
      return(piece);
    }
  } else {
    if (code1 != NOCODE) {
      if (code2 != NOCODE) {
	return(code1 == code2 ? 0.0 : 1.0); /* code1 vs code2 */
      } else {
	assert(f2 != NULL);
	return(1.0 - f2[code1]); /* code1 vs. f2 */
      }
    } else {
      if (code2 != NOCODE) {
	assert(f1 != NULL);
	return(1.0 - f1[code2]); /* f1 vs code2 */
      } else {
	assert(f1 != NULL && f2 != NULL);	/* f1 vs. f2 */
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
(P->weights[i] > 0 && P->codes[i] == NOCODE ? &P->vectors[nCodes*(IVECTOR++)] : NULL);

void ProfileDist(profile_t *profile1, profile_t *profile2, int nPos,
		 /*OPTIONAL*/distance_matrix_t *dmat,
		 /*OUT*/besthit_t *hit) {
  double top = 0;
  double denom = 0;
  int iFreq1 = 0;
  int iFreq2 = 0;
  int i = 0;
  for (i = 0; i < nPos; i++) {
      double weight = profile1->weights[i] * profile2->weights[i];

      float *f1 = GET_FREQ(profile1,i,/*IN/OUT*/iFreq1);
      float *f2 = GET_FREQ(profile2,i,/*IN/OUT*/iFreq2);
      if (profile1->weights[i] > 0 && profile2->weights[i] > 0) {
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

/* bionjWeight is the weight of the first sequence (between 0 and 1),
   or -1 to do the average.
   */
profile_t *AverageProfile(profile_t *profile1, profile_t *profile2, int nPos,
			  distance_matrix_t *dmat,
			  double bionjWeight) {
  int i;
  if (bionjWeight < 0) {
    bionjWeight = 0.5;
  }

  /* First, set codes and weights and see how big vectors will be */
  profile_t *out = NewProfile(nPos);

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
profile_t *OutProfile(profile_t **profiles, int nProfiles, int nPos,
		      distance_matrix_t *dmat) {
  int i;			/* position */
  int in;			/* profile */
  profile_t *out = NewProfile(nPos);

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
  return(out);
}

void UpdateOutProfile(/*IN/OUT*/profile_t *out, profile_t *old1, profile_t *old2,
		      profile_t *new, int nActiveOld, int nPos,
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
  if (fclose(fp) != 0) {
    fprintf(stderr, "Error reading %s\n",filename);
    exit(1);
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

void ReliabilityNJ(/*IN/OUT*/NJ_t *NJ) {
  /* For each non-root node N, with children A,B, parent P, sibling C, and grandparent G,
     we test the reliability of the split (A,B) versus rest by comparing the profiles
     of A, B, C, and the "up-profile" of P.

     Each node's upProfile is the average of its sibling's (down)-profile + its parent's up-profile
     (If node's parent is the root, then there are two siblings and we don't need an up-profile)

     To save memory, we do depth-first-search down from the root, and we only keep
     up-profiles for nodes in the active path.
  */
  int iChild;
  int i;

  if (NJ->nSeq <= 3)
    return;			/* nothing to do */

  int *path = (int*)mymalloc(sizeof(int)*NJ->maxnodes);
  int path_size = 0;
  int *work = (int*)mymalloc(sizeof(int)*NJ->maxnodes);
  int work_size = 0;

  /* For each node the upProfile, but we keep these as NULL for non-active nodes */
  profile_t **upProfiles = (profile_t**)mymalloc(sizeof(profile_t*)*NJ->maxnodes);
  for (i=0; i<NJ->maxnodes; i++) upProfiles[i] = NULL;

  /* Start with children of root, and with path as just the root */
  for (iChild = 0; iChild < NJ->child[NJ->root].nChild; iChild++) {
    work[work_size++] = NJ->child[NJ->root].child[iChild];
    if(verbose>10) fprintf(stderr,"Add to work %d\n",work[work_size-1]);
  }
  path[path_size++] = NJ->root;
  if(verbose>10) fprintf(stderr,"Add to path %d\n",path[path_size-1]);

  while(work_size>0) {
    int node = work[--work_size];
    int parent = NJ->parent[node];
    int child1,child2;

    if (node < NJ->nSeq) {
      if(verbose>10) fprintf(stderr,"Skipping leaf %d\n", node);
      continue;			/* for leaves, nothing to do */
    }

    /* last position in path should be parent */
    while(path[path_size-1] != NJ->parent[node]) {
      /* Forget profile of parent */
      int old = path[--path_size];
      upProfiles[old] = FreeProfile(upProfiles[old], NJ->nPos);
      if(verbose>10) fprintf(stderr,"Free upprofile %d\n",old);
      assert(path_size >= 1);
    }

    /* record child1,child2 -- these give 2 of our 4 profiles */
    assert(NJ->child[node].nChild==2);
    child1 = NJ->child[node].child[0];
    child2 = NJ->child[node].child[1];

    assert(upProfiles[node] == NULL);

    /* get the other 2 profiles */
    profile_t *outprofiles[2];
    if (parent == NJ->root) {
      int nSibs = 0;
      int sibs[3];
      for(iChild=0;iChild<NJ->child[NJ->root].nChild;iChild++) {
	int child = NJ->child[NJ->root].child[iChild];
	if (child != node) sibs[nSibs++] = child;
      }
      assert(nSibs==2);
      outprofiles[0] = NJ->profiles[sibs[0]];
      outprofiles[1] = NJ->profiles[sibs[1]];
      if(verbose>10) fprintf(stderr,"outprofiles for node %d -- %d %d (2 sibs)\n",node,sibs[0],sibs[1]);
    } else {
      int sibling = Sibling(NJ,node);
      assert(sibling>=0);
      assert(upProfiles[parent] != NULL);
      outprofiles[0] = upProfiles[parent];
      outprofiles[1] = NJ->profiles[sibling];
      if(verbose>10) fprintf(stderr,"outprofiles for node %d -- %d %d (parent & sibling)\n",node,parent,sibling);
    }

    /* Given our 4 profiles, compute support values and save
       average of the outprofiles in upProfiles
    */
    if(verbose>2) fprintf(stderr,"SplitSupport for %d\n",node);
    NJ->support[node] = SplitSupport(NJ->profiles[child1],
				     NJ->profiles[child2],
				     outprofiles[0],
				     outprofiles[1],
				     NJ->distance_matrix,
				     NJ->nPos);
    upProfiles[node] = AverageProfile(outprofiles[0], outprofiles[1], NJ->nPos,
				      NJ->distance_matrix,
				      /*no-weight*/-1.0);

    if(verbose>10) fprintf(stderr, "Set upprofile %d\n", node);

    /* Add children to work and add self to path */
    for(iChild=0;iChild<NJ->child[node].nChild;iChild++) {
      work[work_size++] = NJ->child[node].child[iChild];
      if(verbose>10) fprintf(stderr, "Add to work %d\n", work[work_size-1]);
    }
    path[path_size++] = node;
    if(verbose>10) fprintf(stderr, "Add to path %d\n", path[path_size-1]);
  }

  while(path_size>0) {
    int node = path[--path_size];
    upProfiles[node] = FreeProfile(upProfiles[node],NJ->nPos);
    if(verbose>10) fprintf(stderr, "Free upprofile %d (cleaning up path)\n", node);
  }
  for(i=0;i<NJ->maxnodes;i++) assert(upProfiles[i]==NULL);
  upProfiles = myfree(upProfiles, sizeof(profile_t*)*NJ->maxnodes);
  path = myfree(path, sizeof(int)*NJ->maxnodes);
  work = myfree(work, sizeof(int)*NJ->maxnodes);
}

profile_t *NewProfile(int nPos) {
  profile_t *profile = (profile_t *)mymalloc(sizeof(profile_t));
  profile->weights = mymalloc(sizeof(float)*nPos);
  profile->codes = mymalloc(sizeof(unsigned char)*nPos);
  profile->vectors = NULL;
  profile->nVectors = 0;
  profile->codeDist = NULL;
  return(profile);
}

profile_t *FreeProfile(profile_t *profile, int nPos) {
    if(profile==NULL) return(NULL);
    myfree(profile->codes, nPos);
    myfree(profile->weights, nPos);
    myfree(profile->vectors, sizeof(float)*nCodes*profile->nVectors);
    myfree(profile->codeDist, sizeof(float)*nCodes*nPos);
    return(myfree(profile, sizeof(profile_t)));
}

int Sibling(NJ_t *NJ, int node) {
  int parent = NJ->parent[node];
  int iChild;
  if(parent==NJ->root) return(-1);
  for(iChild=0;iChild<NJ->child[parent].nChild;iChild++) {
    if(NJ->child[parent].child[iChild] != node)
      return (NJ->child[parent].child[iChild]);
  }
  assert(0);
  return(-1);
}

/* Computes support for (A,B),(C,D) compared to that for (A,C),(B,D) and (A,D),(B,C) */
double SplitSupport(profile_t *pA, profile_t *pB, profile_t *pC, profile_t *pD,
		    /*OPTIONAL*/distance_matrix_t *dmat,
		    int nPos) {
  /* Support1 = Support(AB|CD over AC|BD) = d(A,C)+d(B,D)-d(A,B)-d(C,D)
     We keep track of sum, sumsq, and total weight
     Note sums are weighted by weight

     Support2 = Support(AB|CD over AD|BC) = d(A,D)+d(B,C)-d(A,B)-d(C,D)

     m1 is mean(Support1) and s1 is standard deviation of Support1
     cor12 is the correlation between Support1 and Support2
  */
  int i;
  double support1sum, support1sumsq, m1, s1;
  double support2sum, support2sumsq, m2, s2;
  double support12sum, cor12;
  double meanMin, sdMin, actualMin, actualZ, support;
  /* If we are doing local bootstrap */
  double *support1 = NULL;
  double *support2 = NULL;
  if (nBootstrap > 0) {
    support1 = (double*)mymalloc(sizeof(double)*nPos);
    support2 = (double*)mymalloc(sizeof(double)*nPos);
  }

  double total_weight = 0;
  double total_weight_sq = 0;
  support1sum = 0;
  support2sum = 0;
  support1sumsq = 0;
  support2sumsq = 0;
  support12sum = 0;

  int iFreqA = 0;
  int iFreqB = 0;
  int iFreqC = 0;
  int iFreqD = 0;
  for (i = 0; i < nPos; i++) {
    float *fA = GET_FREQ(pA, i, /*IN/OUT*/iFreqA);
    float *fB = GET_FREQ(pB, i, /*IN/OUT*/iFreqB);
    float *fC = GET_FREQ(pC, i, /*IN/OUT*/iFreqC);
    float *fD = GET_FREQ(pD, i, /*IN/OUT*/iFreqD);
    double support1piece = 0;
    double support2piece = 0;

    double weight = pA->weights[i];
    if(pB->weights[i] < weight) weight = pB->weights[i];
    if(pC->weights[i] < weight) weight = pC->weights[i];
    if(pD->weights[i] < weight) weight = pD->weights[i];

    if (weight > 0) {
      /* Use the minimum of all the weights. (Should I be using the product instead?) */
      total_weight += weight;
      total_weight_sq += weight*weight;

      double pieceAB = ProfileDistPiece(pA->codes[i], pB->codes[i], fA, fB, dmat, NULL);
      double pieceAC = ProfileDistPiece(pA->codes[i], pC->codes[i], fA, fC, dmat, NULL);
      double pieceAD = ProfileDistPiece(pA->codes[i], pD->codes[i], fA, fD, dmat, NULL);
      double pieceBC = ProfileDistPiece(pB->codes[i], pC->codes[i], fB, fC, dmat, NULL);
      double pieceBD = ProfileDistPiece(pB->codes[i], pD->codes[i], fB, fD, dmat, NULL);
      double pieceCD = ProfileDistPiece(pC->codes[i], pD->codes[i], fC, fD, dmat, NULL);

      support1piece = pieceAC+pieceBD-pieceAB-pieceCD;
      support1sum += weight*support1piece;
      support1sumsq += weight*support1piece*support1piece;

      support2piece = pieceAD+pieceBC-pieceAB-pieceCD;
      support2sum += weight*support2piece;
      support2sumsq += weight*support2piece*support2piece;

      support12sum += weight*support1piece*support2piece;
    }
    if (nBootstrap>0) {
      support1[i] = support1piece;
      support2[i] = support2piece;
    }
    if(verbose>11) fprintf(stderr, "Pos%d weight %f support1piece %f support2piece %f\n",
			   i, weight, support1piece, support2piece);
  }
  assert(iFreqA == pA->nVectors);
  assert(iFreqB == pB->nVectors);
  assert(iFreqC == pC->nVectors);
  assert(iFreqD == pD->nVectors);

  if (support1sum < 0 || support2sum < 0)
    nSuboptimalSplits++;	/* Another split seems superior */
  if (total_weight < 2.0) return(0);

  if (nBootstrap > 0) {
    int nSupport = 0;
    int iBoot;
    for (iBoot=0;iBoot<nBootstrap;iBoot++) {
      double s1 = 0;
      double s2 = 0;
      for (i=0;i<nPos;i++) {
	int pos = (int)(knuth_rand() * nPos);
	if(pos<0) pos=0;
	if(pos==nPos) pos=nPos-1;
	s1 += support1[pos];
	s2 += support2[pos];
      }
      if (s1>0 && s2>0) nSupport++;
    }
    support1 = myfree(support1,nPos*sizeof(double));
    support2 = myfree(support2,nPos*sizeof(double));
    return( nSupport/(double)nBootstrap );
  } else {
    /* mean((X-mean(X))**2) = mean(X*X) - mean(X)**2 
       var(X) = sqrt(N/(N-1)) * sqrt(mean((X-meanX)**2))
       covar(X,Y) = mean(XY) - mean(X)*mean(Y)
    */
    m1 = support1sum/total_weight;
    s1 = sqrt(support1sumsq/total_weight - m1*m1);
    if(s1 < 0.001) s1=0.001;
    
    m2 = support2sum/total_weight;
    s2 = sqrt(support2sumsq/total_weight - m2*m2);
    if(s2 < 0.001) s2=0.001;
    
    cor12 = (support12sum/total_weight - m1*m2)/(s1*s2);
    /* cor12 might be out of bounds if s1 or s2 was very small */
    if (cor12 > 0.999) cor12 = 0.999;
    if (cor12 < -0.999) cor12 = -0.999;
    
    /* An unbiased estimator of the variance of n variables is normally mean squared deviance / (sample_size-1)
       That can be derived from the following:
       Given independent standard normal variables xi, then
       SumSquaredDeviance(xi) = sum(square(x-mean(x))) = sum(square(x)) - square(sum(x))/n
       E(s.s.d.) = n*E(square(x)) - E(square(sum(x)))/n = n -  Var(sum(x))/n = n - n/n = n-1
       
       Similarly, if we have weighted observations with weights wi so that sum(wi) = 1, then
       MeanSquaredDeviance(xi) = sum(wi*square(xi-sum(wj*xj))) = sum(wi*square(xi)) - square(sum(wi*xi))
       E(m.s.d.) = sum(wi*E(square(xi))) - Var(sum(wi*xi)) = sum(wi) - sum(square(wi)) = 1 - sum(square(wi)) <= (n-1)/n
       
       If variables have a variance other than one then we can write
       E(m.s.d) = Var(xi) * (1 - sum(square(wi)))
       
       But we really want the variance of the mean, so we need another factor of total_weight_sq, as
       Var(mean) = Var(sum(wi*xi)) = sum(wi**2) * Var(xi) = sum(wi**2) * m.s.d. / (1 - sum(square(wi)))
    */
    total_weight_sq /= total_weight * total_weight; /* for summing to 1 */
    s1 *= sqrt(total_weight_sq/(1-total_weight_sq)); /* factor = 1/(n-1) if all weights are equal */
    s2 *= sqrt(total_weight_sq/(1-total_weight_sq));
    
    if(verbose>2) fprintf(stderr,"Mean1 %f SD1 %f Mean2 %f SD2 %f cor %f\n", m1, s1, m2, s2, cor12);
    
    /* The mean and variance of the minimum of two normal variables is given by
       Offer Kella (1986) Commun. Statistc.-Theory Meth. 15:3265-3276
    */
    meanMin = -sqrt((s1*s1 + s2*s2 - 2*cor12*s1*s2)/(2*M_PI));
    sdMin = sqrt((s1*s1 + s2*s2)*(1-1/M_PI)*0.5 + s1*s2*cor12/M_PI);
    actualMin = m1 < m2 ? m1 : m2;
    actualZ = sdMin > 0 ? (actualMin-meanMin)/sdMin : 0; /* higher is better */
    
    /* P(z > observed | selected this join)  
       = P(z > observed | P > 0)
       = 3 * P(z > observed) by symmetry
    */
    support = 1 - 3 * pnorm(-actualZ); /* Use -actualZ to get P(z >= actualZ) */
    if(support<0) support=0;
    
    if(verbose>2) fprintf(stderr,"MeanMin %f SDMin %f actualMin %f actualZ %f support %f\n",
			  meanMin,sdMin,actualMin,actualZ,support);
    return(support);
  }
  assert(0);
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
  if (NJ->nOutDistActive[join->j] > nDiffAllow)
    SetOutDistance(NJ, join->j, nActive);
  join->criterion = join->dist - (NJ->outDistances[join->i]+NJ->outDistances[join->j])/(double)(nActive-2);
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
  double pdistOutWithoutA = (nActive-1)
    * (dist.dist * dist.weight * nActive - NJ->selfweight[iNode] * NJ->selfdist[iNode])
    / (dist.weight * nActive - NJ->selfweight[iNode]);
  NJ->outDistances[iNode] =  pdistOutWithoutA
    - NJ->diameter[iNode] * (nActive-1) - (NJ->totdiam - NJ->diameter[iNode]);
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
    if (tophits[seed] != NULL) {
      if(verbose>2) fprintf(stderr, "Skipping seed %d\n", seed);
      continue;
    }
    if(verbose>2) fprintf(stderr,"Trying seed %d\n", seed);
    SetBestHit(seed, NJ, /*nActive*/NJ->nSeq, /*OUT*/&bestjoin, /*OUT*/besthitsSeed);

    /* sort & save top hits of self. besthitsSeed is now sorted. */
    tophits[seed] = SortSaveBestHits(besthitsSeed, seed, /*IN-SIZE*/NJ->nSeq, /*OUT-SIZE*/m);

    /* find "close" neighbors and compute their top hits */
    double neardist = besthitsSeed[2*m-1].dist * tophitsClose;
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

  int iBest;
  for(iBest=0; iBest<nTopHits; iBest++) {
    besthit_t *hit = &tophits[iBest];
    if(hit->j < 0) continue;	/* empty slot */
    assert(hit->i == iNode);

    /* Walk up to active node and compute new distance value if necessary */
    int j = hit->j;
    while(NJ->parent[j] >= 0) j = NJ->parent[j];
    if (iNode == j) continue;
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
		/*IN/OUT*/besthit_t *visible) {
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

  tophitAge[newnode] = tophitAge[NJ->child[newnode].child[0]];
  if (tophitAge[newnode] < tophitAge[NJ->child[newnode].child[1]])
    tophitAge[newnode] = tophitAge[NJ->child[newnode].child[1]];
  tophitAge[newnode]++;

  /* If top hit ages always match, then log2(m) would mean a refresh after
     m joins, which is about what we want.
  */
  int tophitAgeLimit = (int)(0.5 + log((double)m)/log(2.0));
  if (tophitAgeLimit < 1) tophitAgeLimit = 1;
  tophitAgeLimit++;		/* make it a bit more conservative, we have tophitsRefresh threshold too */

  /* UniqueBestHits eliminates hits to self, so if nUnique==nActive-1,
     we've already done the exhaustive search.
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
  } else {
    /* need to refresh: set top hits for node and for its top hits */
    if(verbose>1) fprintf(stderr,"Top hits for %d by refresh (%d unique age %d) nActive=%d\n",
			  newnode,nUnique,tophitAge[newnode],nActive);
    nRefreshTopHits++;
    tophitAge[newnode] = 0;

    /* update all out-distances */
    int iNode;
    for (iNode = 0; iNode < NJ->maxnode; iNode++) {
      if (NJ->parent[iNode] < 0)
	SetOutDistance(/*IN/OUT*/NJ, iNode, nActive);
    }

    /* exhaustively get the best 2*m hits for newnode */
    besthit_t *allhits = (besthit_t*)mymalloc(sizeof(besthit_t)*NJ->maxnode);
    assert(2*m <= NJ->maxnode);
    SetBestHit(newnode, NJ, nActive, /*OUT*/&visible[newnode], /*OUT*/allhits);
    qsort(/*IN/OUT*/allhits, NJ->maxnode, sizeof(besthit_t), CompareHitsByCriterion);

    /* set its top hit list  */
    tophits[newnode] = SortSaveBestHits(allhits, newnode, /*nIn*/NJ->maxnode, /*nOut*/m);

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
      uniqueList2 = myfree(uniqueList2, 3*m*sizeof(besthit_t));

      visible[iNode] = tophits[iNode][GetBestFromTopHits(iNode,NJ,nActive,tophits[iNode],m)];
      ResetVisible(NJ, nActive, tophits[iNode], m, /*IN/OUT*/visible);
    }
    bothList = myfree(bothList,3*m*sizeof(besthit_t));
    allhits = myfree(allhits,sizeof(besthit_t)*NJ->maxnode);
  }
  /* Still need to set visible[newnode] and reset */
  visible[newnode] = tophits[newnode][GetBestFromTopHits(newnode,NJ,nActive,tophits[newnode],m)];
  ResetVisible(NJ, nActive, tophits[newnode], m, /*IN/OUT*/visible);

  uniqueList = myfree(uniqueList, 2*m*sizeof(besthit_t));

  /* Forget top-hit list of children */
  int c;
  for(c = 0; c < NJ->child[newnode].nChild; c++) {
    int child = NJ->child[newnode].child[c];
    tophits[child] = myfree(tophits[child], m*sizeof(besthit_t));
  }
}

void ResetVisible(NJ_t *NJ, int nActive,
		   /*IN*/besthit_t *tophits,
		   int nTopHits,
		   /*IN/UPDATE*/besthit_t *visible) {
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
      if(verbose>5) fprintf(stderr,"NewVisible %d %d %f %f\n",
			    hit->j,visible[hit->j].j,visible[hit->j].dist,visible[hit->j].criterion);
    } else {
      /* see if this is a better hit -- if it is, "reset" */
      SetCriterion(/*IN/OUT*/NJ, nActive, /*IN/OUT*/&visible[hit->j]);
      if (hit->criterion < visible[hit->j].criterion) {
	visible[hit->j] = *hit;
	visible[hit->j].j = visible[hit->j].i;
	visible[hit->j].i = hit->j;
	if(verbose>5) fprintf(stderr,"ResetVisible %d %d %f %f\n",
			      hit->j,visible[hit->j].j,visible[hit->j].dist,visible[hit->j].criterion);
	nVisibleReset++;
      }
    }
  } /* end loop over hits */
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
		    /*IN/OUT*/besthit_t **tophits,
		    /*OUT*/besthit_t *join) {
  besthit_t *visibleSorted = mymalloc(sizeof(besthit_t)*nActive);
  int nVisible = 0;		/* #entries in visibleSorted */
  int iNode;
  for (iNode = 0; iNode < NJ->maxnode; iNode++) {
    /* skip joins involving stale nodes or joins we've already saved */
    if (NJ->parent[iNode] >= 0) continue;
    assert(visible[iNode].i == iNode);
    int j = visible[iNode].j;
    if(NJ->parent[j] >= 0) continue;
    if (j < iNode && visible[j].j == iNode) continue;

    assert(nVisible < nActive);
    visibleSorted[nVisible++] = visible[iNode];
  }
  assert(nVisible > 0);

  qsort(/*IN/OUT*/visibleSorted,nVisible,sizeof(besthit_t),CompareHitsByCriterion);

  /* Only keep the top m items */
  if(nVisible > m) nVisible = m;

  int iBest;
  for (iBest = 0; iBest < nVisible; iBest++)
    SetCriterion(/*IN/OUT*/NJ, nActive, /*IN/OUT*/&visibleSorted[iBest]);

  qsort(/*IN/OUT*/visibleSorted,nVisible,sizeof(besthit_t),CompareHitsByCriterion);

  *join = visibleSorted[0];
  visibleSorted = myfree(visibleSorted, sizeof(besthit_t)*nActive);

  if(fastest) return;

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

hashstrings_t *DeleteHashtable(hashstrings_t* hash) {
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
    
