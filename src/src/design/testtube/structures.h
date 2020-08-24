
#ifndef TESTTUBE_STRUCTURES_H
#define TESTTUBE_STRUCTURES_H

#include <shared.h>

typedef struct OPTIONS_T_ {
  DBL_TYPE temperature; // temperature (K)
  DBL_TYPE sodium;      // sodium concentration (M)
  DBL_TYPE magnesium;   // magnesium concentration (M)
  DBL_TYPE min_ppair_saved;  // minimum pair probability saved
  DBL_TYPE M_unfavorable;   // scaled number of unfavorable leaf mutations to allow
  DBL_TYPE M_leafopt;   // number of failed leaf reoptimizations to allow
  DBL_TYPE M_reseed;   // number of base pairs to reseed
  DBL_TYPE f_split;  // minimum pair prob of helix for ppair decomp
  DBL_TYPE f_passive;     // fraction that concentrations are deflated by
  DBL_TYPE f_stringent; // Margin for relaxation during tree optimization
  DBL_TYPE f_redecomp; 
  DBL_TYPE f_refocus;
  DBL_TYPE bonus_per_split;
  DBL_TYPE gc_init_prob;  // Probability of choosing GC vs AU on initialization
  // undesired structures are added or redecomposition is repeated until
  // at least include_frac of the excess defect of the parent is captured
  // in the child
  int material;    // DNA1998/RNA1995/RNA1999 type
  int dangle_type;
  int use_long_helix;
  unsigned int seed;
  int H_split;
  int N_split;
  int allow_wobble;
  int allow_mismatch;
  int print_leaves;
  int print_steps;
  int single_decomp;
  int max_print_steps;
  int disable_defect_weights;
  int disable_focus;
  int forbid_splits;
  int include_all;
  int fake_dummies;
  int include_dummies;
  int designing; // We are designing! (used for output)

  char * file_prefix;

  DBL_TYPE start_time;
  DBL_TYPE allowed_opt_time;
} options_t;
/*****************************************************************************
 * Sequence properties
 * **************************************************************************/
typedef struct SEQUENCE_T_ {
  int * nucs;
  int n;
} sequence_t;

typedef struct SEQSPEC_LIST_T_ {
  sequence_t * specs;
  char ** names;
  int n;
  int cap;
} seqspec_list_t;


typedef struct SEQUENCE_SPEC_T_ {
  sequence_t nucs;
  sequence_t comp_map;
  seqspec_list_t domains;
  seqspec_list_t strands;
} sequence_spec_t;

typedef struct SEQUENCE_LIST_T_ {
  sequence_t * seqs;
  int n;
} sequence_list_t;

typedef struct SEQSTATE_T_ {
  sequence_t nucspec;
  sequence_t dumspec;
  sequence_list_t strands;
  sequence_list_t domains;
} seqstate_t;

typedef struct MUTATION_T_ {
  int * ids;
  int * nucs;
  int * dum;
  int n;
  int cap;
} mutation_t;

typedef struct MUTATION_LIST_T_ {
  mutation_t * muts;
  int n;
  int cap;
} mutation_list_t;

typedef struct SPLIT_EL_T {
  int lsplit;// proposed left split point
  int rsplit;// proposed right split point
  struct SPLIT_EL_T * next;
} split_el_t;

typedef struct SPLIT_LIST_T_ {
  double cost;
  double ppair;
  split_el_t * head;
  split_el_t * tail;
  struct SPLIT_LIST_T_ * next;
} split_list_t;

typedef struct SPLIT_TRACKER_T_ {
  int n;
  split_list_t * head;
  split_list_t * tail;
} split_tracker_t;

// I really should make this a DAG, but simplicity first
typedef struct STRUC_TREE_T_ {
  int * seg_start;   
  int * seg_stop;   
  int n_segments;

  int * assumed_i;
  int * assumed_j;
  int n_assumed;

  split_tracker_t * forbidden;
  struct STRUC_TREE_T_ ** children; // If multiple decompositions exist, see the next one
  int n_children;
} struc_tree_t;


typedef struct DESIGN_STRUC_T_ {
  int * strands;      // index of strands in the current structure.
  int n_strands;      // number of strands in the current structure
  int * struc; // pairing of each domain
  int n_nucs;      // the number of domains
  int * split_forbidden;
  int symmetry;

  int modifiable;

  struc_tree_t * tree;

  struct DESIGN_STRUC_T_ * next;
} design_struc_t ;

typedef struct PPAIRS_T_ {
  DBL_TYPE ** ppairs;
  int ** ppairs_i;
  int ** ppairs_j;
  int * ppairs_n;
  int n;
} ppairs_t;

typedef struct STRUC_STATE_T_ {
  design_struc_t * struc;
  int n; // number of leaves of the decomposition
  struc_tree_t * tree; 
} struc_state_t;

typedef struct design_state_t {
  struc_state_t * states;
  int n_strucs;
  int cap_strucs;
} design_state_t ;

typedef struct MOD_STATE_T_ {
  int ** mod;
  int * n_nodes;
  int n_strucs;
} mod_state_t;

typedef struct DESIGN_TUBE_T_ {
  DBL_TYPE * target_x; // The desired concentrations for the 
  int * target;         // 1 if structure is a positive target structure.
  int * included;       // 1 if full calc is performed for structure
  int * included_ind;   // The index of the full calculation
  int * generated_ind;  // Index of strand ordering (for unincluded strucs)
  int maxsize;
  DBL_TYPE stop;
  int n_strucs;         // The number of structures
  int cap_strucs;
} design_tube_t;

// Should actually separate out sequence spec from sequence result
typedef struct DESIGN_SPEC_T_ {
  sequence_spec_t seqs;
  options_t opts;

  design_struc_t * strucs;
  char ** struc_names;
  int n_strucs;
  int cap_strucs;

  design_tube_t * tubes;
  char ** tube_names;
  int n_tubes;
  int cap_tubes;

  int ** orderings;
  int * symmetry;
  int * struc_map;    // index of current ordering in strucs
  int * n_strands;
  int n_orderings;
  int cap_orderings;
} design_spec_t;

typedef struct STATE_TUBE_T_ {
  design_tube_t * tube;
  int * included;
} state_tube_t ;

// Any of the fields aside from children may be blank for interior nodes.
// They might also be filled out, I'm not sure which is preferable yet.
typedef struct RESULT_TREE_T_ {
  int * sequence;         // Sequence of node
  int * native_map;       // Map from nucleotide indices to native indices
  int * dummy_flag;       // 1 if nuc is dummy 0 otherwise
  DBL_TYPE * nuc_defects; // Defects within the node
  int n_nucs;             // Number of node nucleotides

  DBL_TYPE pfunc;         // reduced partition function of the node
  DBL_TYPE eval_time;     // Evaluation time for the node (or total of children)
  int n_children;         // number of children the node has
  struct RESULT_TREE_T_ ** children; // map to the children of the node

  DBL_TYPE * ppairs;
  int * ppairs_i;
  int * ppairs_j;
  int ppairs_n;
  int ppairs_cap;
} result_tree_t;

typedef struct RESULT_STRUC_T_ {
  int * sequence;       // The nucleotide sequence
  int * f_sequence;  // The nucleotide sequence for fake bases
  int * nuc_ids;
  int * structure;      // The base-pairing structure
  int * modifiable;     // 1 if base is modifiable, 0 otherwise
  DBL_TYPE * defects;   // nucleotide defects
  DBL_TYPE * f_defects; // nucleotide defects for fake bases
  int n_nucs;           // The number of nucleotides
  int cap_nucs;         // the capacity of the nucleotide-indexed results

  int * breaks;         // Native strand breaks in nucleotides
  int n_breaks;         // The number of strand breaks

  DBL_TYPE pfunc; 
  DBL_TYPE time;

  DBL_TYPE * ppairs;
  int * ppairs_i;
  int * ppairs_j;
  int ppairs_n;
  int ppairs_cap;

  DBL_TYPE defect;

  result_tree_t * tree;
} result_struc_t;

typedef struct RESULT_TUBE_T_ {
  DBL_TYPE * x;
  DBL_TYPE * target_x;
  int * target;         // 1 if structure is a positive target structure.
  int * included;       // 1 if full calc is performed for structure
  int * included_ind;   // The index of the structure in result_struc_t * array
  int * generated_ind;  // Index of strand ordering

  DBL_TYPE defect;
  int n_strucs;
  int cap_strucs;
} result_tube_t;

typedef struct RESULT_T_ {
  result_struc_t * strucs;
  int n_strucs;
  int cap_strucs;

  result_tube_t * tubes;
  int n_tubes;
  int cap_tubes;

  int ** orderings;
  int * n_strands;
  DBL_TYPE * dG;
  DBL_TYPE * eval_time;
  int * included;
  int * order_len;
  int * struc_map;    // index of current ordering in strucs
  int n_orderings;
  int cap_orderings;

  DBL_TYPE total_defect;

  DBL_TYPE elapsed_time;
  DBL_TYPE root_time;

  seqstate_t * offtarget_seqs;
  int tot_n_strands;
  int id;
} result_t;

enum recursiontype {C_U, C_E, C_B, C_S, C_P, C_G, C_GL, C_GR, C_GRS};

typedef struct RECURSION_T_ {
  enum recursiontype type; /* The type of recursion. */
  int lower; /* Lower, middle, and upper indices. */
  int midl;
  int midu;
  int upper;
  double score;
} recursion_el_t;

#endif // TESTTUBE_STRUCTURES_H
