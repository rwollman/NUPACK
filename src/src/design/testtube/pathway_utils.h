#ifndef __PATHWAY_UTILS_H__
#define __PATHWAY_UTILS_H__
/*******************************************************************
 *******************************************************************
 * pathway_design.c
 * These functions are used for calculating defects and optimizing
 * leaves in a way very specific to pathway design.
 *
 * Also included are brief unit tests used for development and test
 * purposes (and called in test_main.c)
 ******************************************************************/

/*
 * parse_design
 * parse and fill out the structures for design
 */
int parse_design(char * filename, 
    design_spec_t * spec);

/*
 * optimize_tubes
 * Optimize the all test tubes included in the design specification
 * subject to the options specified in opts
 */
int optimize_tubes(result_t * res, seqstate_t * res_seqs, 
    design_spec_t * spec);

int copy_dummy(seqstate_t * res);

/*
 * optimize_forest
 * Optimize the full decomposition trees.
 */
int optimize_forest(result_t * res, seqstate_t * seqs, 
    design_spec_t * spec);

int optimize_leaves(result_t * result, seqstate_t * seqs,
    design_state_t * state, DBL_TYPE goal_frac, 
    design_spec_t * spec);

int mutate_leaves(result_t * result, seqstate_t * seqs,
    design_state_t * state, DBL_TYPE goal_frac,
    design_spec_t * spec);

int get_nuc_ids(int * nuc_ids, design_struc_t * struc, 
    design_spec_t * spec) ;
/*
 * Get the remapped, thresholded, and normalized ensemble defect of the 
 * test tubes.
 */
DBL_TYPE get_defect(result_t * res);

int reseed_next_leaf(seqstate_t * seqs, result_t * res,
    design_spec_t * spec) ;
int reseed_leaf(seqstate_t * seqs, result_t * res, 
      int r_tube, int r_struc, int r_node, design_spec_t * spec) ;


int pick_mutation(mutation_t * mut, result_t * res, seqstate_t * seqs, design_spec_t * spec);

int pick_mutation_location(int * nuc_id, int * dummy, 
    result_t * res, design_spec_t * spec);

/*
 * Get the true sum of test tube ensemble defect fractions
 */
DBL_TYPE get_true_defect(result_t * res, 
    design_spec_t * spec);

/*******************************************************************
 * pathway_utils.c
 * These functions are used for decomposition, calculation of 
 * various intermediary quantities and manipulating the data 
 * structures used in pathway design. This is meant to be a 
 * collection of functions only really used by the functions in
 * pathway_design.c
 *
 * Many of these are tested directly in the function test_A
 ******************************************************************/

/*
 * init_sequence_table
 * Initialize the sequence table (used for holding length, domain, and 
 * strand specifications and state)
 */
int init_sequence_spec(sequence_spec_t * seqs);
int init_design_spec(design_spec_t * spec);

int init_seqspec_list(seqspec_list_t * seqs);
int delete_sequence_seqspec(seqspec_list_t * seqs, int i_str);
int append_sequence_seqspec(seqspec_list_t * dest, sequence_t * src,
    const char * name);
int free_seqspec_list(seqspec_list_t * seqs);

int init_seqstate(seqstate_t * seqs, design_spec_t * spec );
int copy_seqstate(seqstate_t * dest,
    seqstate_t * src);
int refresh_seqstate(seqstate_t * seqs, design_spec_t * seqspec);

void free_seqstate(seqstate_t * seqs);

int init_sequence(sequence_t * seq, int n_nucs);
int init_copy_sequence(sequence_t * dest, sequence_t * src);
int resize_sequence(sequence_t * seq, int n_nucs);
int copy_sequence(sequence_t * dest, sequence_t * src);
void free_sequence(sequence_t * seq);

int get_n_mutable_nucleotides(design_spec_t * spec);

int init_sequence_list(sequence_list_t * seqlist, seqspec_list_t * spec);
int init_copy_sequence_list(sequence_list_t * dest, sequence_list_t * src);
int copy_sequence_list(sequence_list_t * dest, sequence_list_t * src);
int append_sequence(sequence_list_t * dest, sequence_t * src);
void free_sequence_list(sequence_list_t * seqs);

/* add_*
 * spec: the design specification to append to
 * *: properties for the given object to add
 */

int add_strand(design_spec_t * spec, char * name,
      char ** domain_names, int n_domains);

int add_strand_basic(design_spec_t * spec, int * ind, char * name, int len);
int add_domain(design_spec_t * spec,
      char * name, char * constraint);
int add_structure_basic(design_spec_t * des, char * name, char * struc);
int add_structure_str(design_spec_t * spec, char * name, char ** strands,
    int n_strands, char * domain_str);
int add_structure(design_spec_t * spec, char * name, 
    int * strands, int n_strands, int * dom_struc, int n_domains);

int add_altstruc_basic(design_spec_t * spec, char * name, char * struc);
int add_altstruc(design_spec_t * spec, int c_struc, int * pairing, int n_nucs);

int define_structure_strands(design_spec_t * spec, char * strucname, 
    char ** strandnames, int n_strands);
int set_concentration(design_spec_t * spec,
    const char * tubename, const char * strucname, DBL_TYPE conc);
int add_tube_basic( design_spec_t * spec,
      char * name, char ** strucnames, int n_strucs) ;
int add_tube_str(design_spec_t * spec, char * name, char ** strucs,
    DBL_TYPE * concs, int n_strucs, DBL_TYPE stop, int maxsize);
int add_tube(design_spec_t * spec, char * name, int * strucs,
    DBL_TYPE * concs, int n_strucs, DBL_TYPE stop, int maxsize);
int make_off_targets(design_spec_t * spec);
int init_struc(design_struc_t* struc, sequence_spec_t * seqs, 
    int * strands, int n_strands,
    int * dom_struc, int n_domains);

int init_opts(options_t * opts);
int init_seqs_random(seqstate_t * seqs, sequence_spec_t * seqspec, 
    options_t * opts);
int init_constraint_random(int * f_seq, int * c_seq, 
    int * f_con, int * c_con, int n_nucs, options_t * opts);

int init_result_tube(result_tube_t * tube);
int init_result_tube_spec(result_tube_t * tube, design_tube_t * spec);
int init_result(result_t * res, design_spec_t * spec);
int copy_result(result_t * dest, result_t * src);

int index_of(const char * name, char ** const names, int n_names);

int print_design_result(FILE * f, result_t * res, seqstate_t * seqs,
    int indent, design_spec_t * spec, char * executable) ;

int print_leafplot(enum NUPACK_DESIGN_STEP step, 
    seqstate_t * seqs,
    design_spec_t * spec);
/*
 * get_sequence_plus
 * Put the sequence string into str with the strand breaks indicated
 * with a + sign
 */
int get_sequence_plus(char * str, result_struc_t * res, 
    options_t * opts);

/*
 * get_pairing
 * convert dot-paren-plus to pairing and strand breaks
 */
int get_pairing(int * pairing, int * breaks, int * n_nucs, int * n_breaks, 
    char * dpp, int n_pos);
/*
 * get_dpp
 * convert a a structure result to dot-paren-plus
 */
int get_dpp(char * str, result_struc_t * res);
int pairs_to_dpp(char * str, int * struc, int n_nucs) ;

int get_spec_dpp(char * str, design_struc_t * des, 
    design_spec_t * spec) ;

int ppair_structure_at_leaf(result_struc_t * res_struc,
    design_struc_t * des_struc, seqstate_t * seqs, 
    design_spec_t * spec, int i_node, DBL_TYPE fsplit);

int ppair_decompose_at_node(result_tree_t * restree, 
    result_struc_t * res, seqstate_t * seqs, struc_tree_t * tree, 
    design_struc_t * struc, design_spec_t * spec, 
    DBL_TYPE fsplit) ;
/*
 * ppair_structure
 * convert a a structure result to dot-paren-plus
 */
int ppair_structure(
    result_struc_t * res,
    design_struc_t * struc,
    struc_state_t * state,
    seqstate_t * seqs,
    design_spec_t * spec);

/* convert_nucs_to_str
 * Convert integer nucs to a string representation
 */
int convert_nucs_to_str(
    char * str,
    int * nucs,
    int n_nucs,
    options_t * opts);

/*
 * convert_str_to_nucs
 * Convert string nucleotides to the integer representation
 */
int convert_str_to_nucs(
    int * nucs,
    char * str,
    int n_nucs);

/*
 * get_struc_length
 * get the length of the structure 
 */
int get_struc_length(design_struc_t * struc, sequence_spec_t * seqspec);

/* 
 * get_max_node_lengths
 * get the maximum lengths of the nodes.
 */
int get_max_node_lengths(int * lengths, struc_state_t * state, 
    sequence_spec_t * seqspec, options_t * opts);

int add_nucleotides(design_spec_t * spec, int * start, int n);
/*
 * init_struc_state
 * Initialize the decomposition state for the given structure
 */
int init_struc_state(struc_state_t * state, design_struc_t * des);
int copy_struc_state(struc_state_t * dest, struc_state_t * src);
int init_states(design_state_t * states, design_spec_t * spec) ;
int copy_states(design_state_t * dest, design_state_t * src);

int init_decomposition(design_spec_t * spec);
int get_decomposition(design_state_t * states, int maxdepth,
    design_spec_t * spec);

/*
 * decompose_greedy_guided
 * Decompose in the Zadeh fashion, greedily splitting children.
 * must include all split points specified in the guide states
 * if guide_state is NULL, it is ignored
 */

int decompose_greedy(
      struc_state_t * state,
      struc_state_t * guide_state,
      design_struc_t * des,
      int max_depth,
      design_spec_t * spec);


/*
 * init_result_struc
 * Initialize the result struc for a given structure/state
 */
int init_result_struc(result_struc_t * res);
int copy_result_struc(result_struc_t * dest, result_struc_t * src);

/*
 * ensure_capacity
 * Ensure, even for after length changes, that there is enough room
 * to hold all the results
 */

int ensure_capacity(result_struc_t * result, struc_state_t * struc);
int ensure_capacity_vals(result_struc_t * result, int n_nucs,
    int n_strands);

int ensure_ppair_cap(
      result_struc_t * dest, int ppaircount
    );
/*
 * fill_in_sequences
 * Fill in the sequence information to the result structure
 */
int fill_in_sequences(result_struc_t * result, design_struc_t * struc,
      seqstate_t * seqs, design_spec_t * spec);

int fill_out_blank_seqs(result_t * res, design_state_t * state, 
    design_spec_t * spec);

int update_structure(result_struc_t * result,
    struc_state_t * struc_state, design_spec_t * spec);

/*
 * update result
 */
int update_result(result_t * result, design_state_t * state, 
    seqstate_t * seqs, design_spec_t * spec);

/*
 * update_structure_result
 * Update the structure result, including moving unmodified node's defects
 * to the correct nucleotides and re-evaluation of the remaining nodes.
 */
int update_structure_result(result_struc_t * result,
    struc_state_t * struc_state, seqstate_t * seqs,
    design_spec_t * spec);

/*
 * evaluate_node
 * Evaluate a single node's partition function and nucleotide defect
 * information.
 */
int evaluate_node(result_struc_t * res, options_t * opts, int node);
int evaluate_undesired(result_t * res, seqstate_t * seqs,
    design_spec_t * spec);

int calculate_concentrations(result_t * res, design_spec_t * seqs);

int print_full_decomposition(FILE * f, 
    design_struc_t * struc,
    design_spec_t * spec,
    int indent);
int print_decomposition_state(FILE * f, struc_state_t * state, 
    design_spec_t * spec, int indent) ;

int print_leaves(FILE * f, 
    design_struc_t * struc,
    design_spec_t * spec,
    int indent);
int print_decompositions(FILE * f,
    design_state_t * states,
    design_spec_t * spec,
    int indent) ;
/*
 * get_nodes_changed_struc
 * Find and return all changed nodes in a particular structure by
 * comparing two different results
 */
int get_nodes_changed_struc(int ** nodes_changed_p, int * n_nodes_p,
    int * cap_nodes_p, result_struc_t * res1, result_struc_t * res2);

/*
 * get_leaf_defect(res, i_tube, i_struc, i_node)
 */
DBL_TYPE get_nodal_defect(result_t * res, int i_tube, int i_struc, int i_node,
    design_spec_t * spec);

/*
 * get_n_leaves
 */
int get_n_leaves(result_struc_t * res);

/*
 * get_nodes_changed
 * Find and return all changed nodes in a list of structures by
 * comparing two sets of results (must be identically ordered)
 */
int get_nodes_changed(int ** strucs_p,
    int ** nodes_p, int * n_nodes_p, int * cap_nodes_p, 
    result_struc_t * res1, result_struc_t * res2, 
    int n_strucs);

int print_domains(FILE * f, 
    seqstate_t * seqs, int indent, design_spec_t * spec);
int print_strands(FILE * f, 
    seqstate_t * seqs, int indent, design_spec_t * spec);
int print_objectives(FILE * f, result_t * res,
    int indent, design_spec_t * spec);
int print_struc_result(FILE * f, result_struc_t * res, int indent,
    options_t * opts);
int print_struc_spec(FILE * f,
    design_struc_t * spec, 
    seqstate_t * seqs,
    int indent);

int print_sequences(FILE * f, seqstate_t * seqs, int indent,
    options_t * opts);

/*
 * General translation utils
 */

int get_comp_code(int nuc_code, int allow_wobble);
int get_nuc_code(char nuc);

/*
 * Mutation utilities
 */
int init_mutation(mutation_t * mut);
int resize_mutation(mutation_t * mut, int size);
int add_mutation_nuc(mutation_t * mut, int id, int nuc, int dummy);
int clear_mutation(mutation_t * mut);
int copy_mutation(mutation_t * dest, mutation_t * src);
int apply_mutation(seqstate_t * seqs, 
    mutation_t * mut);
int sample_mutation(mutation_t * mut, int nuc_id, int dummy, 
      seqstate_t * seqs, design_spec_t * spec);
void print_mutation(FILE * f, mutation_t * mut, int indent);
void print_mutation_weights(FILE * f, result_t * res, 
    mutation_t * mut, int indent) ;
void free_mutation(mutation_t * mut);

/*
 * Mutation list utilities
 */
int init_mutation_list(mutation_list_t * muts);
int resize_mutation_list(mutation_list_t * muts, int size) ;
int find_mutation(mutation_list_t * muts, mutation_t * mut);
int append_mutation(mutation_list_t * muts, mutation_t * mut);
int clear_mutation_list(mutation_list_t * muts);
void free_mutation_list(mutation_list_t * muts);


/*
 * split tracker utilities (for decomposition)
 */
split_tracker_t * alloc_split_tracker(void);
int init_split_tracker(split_tracker_t * tracker);
int append_split_tracker(split_tracker_t * tracker, split_list_t * prev_els, 
    split_el_t * new_el);
int reduce_splits(split_tracker_t * tracker, int * split_markers, int n_nucs,
    int * to_full, int n_full_nucs,
    design_spec_t * spec);

int copy_split_tracker(split_tracker_t * dest, split_tracker_t * src);
int get_unique_id(split_tracker_t * tracker, split_list_t * prev_els);
int free_split_tracker(split_tracker_t * tracker);
int destroy_split_tracker(split_tracker_t * tracker);

/*
 * Split list utilities
 */

split_list_t * alloc_split_list(void);
int init_split_list(split_list_t * l);
int append_split_list(split_list_t * l, split_el_t * el);
int copy_split_list(split_list_t * dest, split_list_t * src);
int free_split_list(split_list_t * l);
int destroy_split_list(split_list_t * l);


struc_tree_t * alloc_struc_tree(int n_segs, int * seg_start, int * seg_stop,
    int n_assume, int * i_assume, int * j_assume);
int free_struc_tree(struc_tree_t * t);
int tree_is_decomp_for(struc_tree_t * tree, design_struc_t * struc);

int is_modifiable(result_t *, int i_struc, int i_node);

/* 
 * Result tree utilities
 */
result_tree_t * get_result_node(result_t * res, int struc, int node,
    DBL_TYPE * prob, design_spec_t * spec);
result_tree_t * alloc_empty_result_tree(void) ;
result_tree_t * alloc_result_tree(struc_tree_t * tree, 
    design_struc_t * struc, design_spec_t * spec);
int print_child_counts(FILE * f, result_tree_t * restree);

int copy_result_tree(result_tree_t * dest, result_tree_t * src);

int map_tree_defects(result_tree_t * res, result_struc_t * struc, 
    design_spec_t * spec);

int sync_result_tree(result_tree_t * res,
    struc_tree_t * tree, design_struc_t * struc, 
    design_spec_t * spec);

int eval_tree(result_tree_t * res, result_struc_t * res_struc, 
    struc_tree_t * tree, design_struc_t * des_struc, design_spec_t * spec);

int eval_leaf(result_tree_t * res, result_struc_t * struc,
    struc_tree_t * tree, design_spec_t * spec);

DBL_TYPE compute_pfunc_merge(result_tree_t * left, result_tree_t * right,
    result_struc_t * res, design_spec_t * spec) ;
int free_result_tree(result_tree_t * tree);
int destroy_result_tree(result_tree_t * tree);
int update_properties(result_struc_t * res, design_struc_t * struc, design_spec_t * spec);

/* free_*(*)
 * free the given struct and all stored descendants
 */
int free_result_tube(result_tube_t * res);
int free_result_struc(result_struc_t * res);
int free_struc_state(struc_state_t * state);
int free_states(design_state_t * state);
int free_struc(design_struc_t * struc);
int free_sequence_spec(sequence_spec_t * seqs);
// int free_level(design_level_t * cur_level);
int free_tube(design_tube_t * tube);
int free_design_spec(design_spec_t * spec);
int free_result(result_t * res);


/*
 * Unused junk
 */

int get_strand_breaks(int * breaks, design_struc_t * struc,
    design_spec_t * spec);

#endif
