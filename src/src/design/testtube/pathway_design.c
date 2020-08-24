#include "pathway_design.h"

static int perturb_sequences(
      seqstate_t * seqs,
      result_t * res,
      design_spec_t * spec
    );
static int update_decomposition(
    result_t * results,
    design_state_t * level_states, 
    seqstate_t * seqs,
    int reset_level, 
    int n_levels, 
    design_spec_t * spec);

static int recalculate_defect(
    result_t * res,
    design_spec_t * spec);

static int get_insertion_point_helper(
      result_tree_t ** parent_p, result_tree_t ** child_p,
      result_tree_t *** parent_ins_point, result_tree_t * parent_tree,
      result_tree_t * child_tree, int i_node, int * c_node, 
      design_spec_t * spec
    ) {
  int i_c;
  if (parent_tree->n_children == 0) {
    *c_node += 1;
  } else {
    check(parent_tree->n_children == child_tree->n_children, 
        "Parent and child trees are out of sync");
    for (i_c = 0; i_c < parent_tree->n_children; i_c++) {
      if (i_node == *c_node && parent_tree->children[i_c]->n_children == 0) {
        *parent_p = parent_tree->children[i_c];
        *child_p = child_tree->children[i_c];
        *parent_ins_point = &(parent_tree->children[i_c]);
        *c_node += 1;
      } else {
        check(ERR_OK == get_insertion_point_helper(parent_p, child_p, 
              parent_ins_point, parent_tree->children[i_c], 
              child_tree->children[i_c], i_node, c_node, spec),
            "Error finding insertion point");
      }
    }
  }
  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

static int get_insertion_point(result_tree_t ** parent_p, 
    result_tree_t ** child_p, result_tree_t *** parent_ins_point,
    result_struc_t * parent, result_struc_t * child,
    int i_node, design_spec_t * spec) {

  int c_node = 0;
  *parent_p = NULL;
  *child_p = NULL;
  *parent_ins_point = NULL;
  if (parent->tree->n_children == 0) {
    check(i_node == 0, "Invalid node in insertion point");
    *parent_p = parent->tree;
    *child_p = child->tree;
    *parent_ins_point = &(parent->tree);
  } else {
    check(ERR_OK == get_insertion_point_helper(parent_p, child_p,
          parent_ins_point, parent->tree, child->tree, i_node, &c_node, spec),
        "Error finding insertion point");
  }
  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

static DBL_TYPE calculate_merge_defect(result_t * parent, result_t * child, 
    int i_struc, int i_node, design_spec_t * spec) {

  result_tree_t * parent_node;
  result_tree_t * child_node;
  result_tree_t ** insert_location;
  DBL_TYPE res = 1e100;

  check(ERR_OK == get_insertion_point(&parent_node, &child_node, &insert_location,
      &(parent->strucs[i_struc]), &(child->strucs[i_struc]),
      i_node, spec),
    "Error in insertion point");

  if (child_node->n_children > 0) {
    check(insert_location, "Insert location came back NULL");
    check(*insert_location == parent_node, "Insert location is incorrect");
    check(parent_node != NULL && child_node != NULL && insert_location != NULL,
      "Failed to find insertion point");

    *insert_location = child_node;
    update_properties(&(parent->strucs[i_struc]), &(spec->strucs[i_struc]), spec);
    recalculate_defect(parent, spec);
    res = get_defect(parent);

    *insert_location = parent_node;
    update_properties(&(parent->strucs[i_struc]), &(spec->strucs[i_struc]), spec);
    recalculate_defect(parent, spec);
    debug("Decomp defect: %Lf %Lf %i %i", res, get_defect(parent), parent_node->n_children, child_node->n_children);
  }

  return res;
error:
  return -1;
}

// static DBL_TYPE calculate_min_merge_ppair(result_t * parent, result_t * child,
//     int i_struc, int i_node, design_spec_t * spec) {
//   DBL_TYPE min_ppair = 1.0;
//   DBL_TYPE res = 1.0;
// 
//   result_tree_t * parent_node;
//   result_tree_t * child_node;
//   result_tree_t ** insert_location;
//   result_tree_t * ch_child;
//   int i_ch;
//   int i_nuc, j_nuc, d_nuc, e_nuc;
//   int found;
// 
//   check(ERR_OK == get_insertion_point(&parent_node, &child_node, &insert_location,
//       &(parent->strucs[i_struc]), &(child->strucs[i_struc]),
//       i_node, spec),
//     "Error in insertion point");
// 
//   check(parent_node != NULL && child_node != NULL && insert_location != NULL,
//     "Failed to find insertion point");
// 
//   int i_pos;
// 
//   res = 0.0;
// 
//   for (i_ch = 0; i_ch < child_node->n_children / 2; i_ch++) {
//     min_ppair = 1.0;
//     ch_child = child_node->children[2 * i_ch + 1];
//     i_nuc = ch_child->native_map[spec->opts.H_split];
//     j_nuc = ch_child->native_map[ch_child->n_nucs - spec->opts.H_split - 1];
//     found = 0;
//     for (i_pos = 0; i_pos < parent_node->ppairs_n; i_pos++) {
//       d_nuc = parent_node->ppairs_i[i_pos];
//       e_nuc = parent_node->ppairs_j[i_pos];
//       if ((i_nuc == d_nuc && j_nuc == e_nuc) 
//         ||(i_nuc - 1 == d_nuc && j_nuc + 1 == e_nuc)) {
//         debug("i %i j %i ppair %Lf", d_nuc, e_nuc, parent_node->ppairs[i_pos]);
//         if (parent_node->ppairs[i_pos] < min_ppair) {
//           min_ppair = parent_node->ppairs[i_pos];
//           found = 1;
//         }
//       }
//     }
//     if (!found) {
//       min_ppair = 0;
//     }
//     res += min_ppair;
//   }
// 
// 
//   return res;
// error:
//   return 1.0;
// }

static int get_spec_and_res_helper(
      result_tree_t ** res_p, struc_tree_t ** struc_tree_p,
      result_tree_t * res, struc_tree_t * struc_tree, 
      int i_node, int * c_node, 
      design_spec_t * spec
    ) {
  int i_c;
  if (res->n_children == 0) {
    if (*c_node == i_node) {
      *res_p = res;
      *struc_tree_p = struc_tree;
    }
    *c_node += 1;
  } else {
    // I should have taken all or none of the children
    check(struc_tree->n_children == res->n_children, 
        "Result and struc trees are out of sync: # child %i != %i",
        struc_tree->n_children, res->n_children);
    for (i_c = 0; i_c < struc_tree->n_children; i_c++) {
      check(ERR_OK == get_spec_and_res_helper(res_p, struc_tree_p, 
            res->children[i_c], struc_tree->children[i_c], 
            i_node, c_node, spec),
          "Error finding insertion point");
    }
  }
  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

static int get_spec_and_res(result_tree_t ** res_p, 
    struc_tree_t ** struc_tree_p, result_struc_t * res, 
    design_struc_t * struc, 
    int i_node, design_spec_t * spec) {

  int c_node = 0;
  *res_p = NULL;
  *struc_tree_p = NULL;

  if (res->tree->n_children == 0) {
    check(i_node == 0, "Invalid node in insertion point");
    *res_p = res->tree;
    *struc_tree_p = struc->tree;
  } else {
    check(ERR_OK == get_spec_and_res_helper(res_p, struc_tree_p,
          res->tree, struc->tree, i_node, &c_node, spec),
        "Error finding insertion point");
  }
  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

int ppair_structure_at_leaf(result_struc_t * res,
        design_struc_t * struc, seqstate_t * seqs, design_spec_t * spec, 
        int i_node, DBL_TYPE fsplit) {
  
  int i_c;
  int i_nuc, j_nuc;
  int i_assume, j_assume, d_nuc, e_nuc, foundassume;
  result_tree_t * restree = NULL;
  struc_tree_t * structree = NULL;

  split_list_t temp_split;
  split_el_t temp_el;

  get_spec_and_res(&restree, &structree, res, struc, i_node, spec);
  check(restree != NULL && structree != NULL, "Error finding spec and res");

  // debug("Adding forbid");
  init_split_list(&temp_split);
  for (i_c = 0; i_c < structree->n_children / 2; i_c++) {
    // Split point used
    for (i_assume = 0; 
        i_assume < structree->children[i_c * 2 + 1]->n_assumed; 
        i_assume++) {

      i_nuc = structree->children[i_c * 2 + 1]->assumed_i[i_assume];
      j_nuc = structree->children[i_c * 2 + 1]->assumed_j[i_assume];
      foundassume = 0;
      for (j_assume = 0; j_assume < structree->n_assumed; j_assume++) {
        d_nuc = structree->assumed_i[i_assume];
        e_nuc = structree->assumed_j[i_assume];
        if ((i_nuc == d_nuc && e_nuc == j_nuc) ||
            (i_nuc == e_nuc && j_nuc == d_nuc)) {
          foundassume = 1;
        }
      }
      if (!foundassume) {
        temp_el.lsplit = i_nuc;
        temp_el.rsplit = j_nuc;
        temp_el.next = NULL;
        append_split_list(&temp_split, &temp_el);
        append_split_tracker(struc->tree->forbidden, &temp_split, NULL);
        free_split_list(&temp_split);
      }
    }
  }
  split_list_t * tl;
  split_el_t * te;
  debug("forbidden:");
  tl = struc->tree->forbidden->head;
  while (tl) {
    debug("Next:");
    te = tl->head;
    while (te) {
      debug("%i %i", te->lsplit, te->rsplit);
      te = te->next;
    }
    tl = tl->next;
  }
  // append_split_tracker(struc->tree->forbidden, &temp_split, NULL);

  free_split_list(&temp_split);

#ifndef NDEBUG
  fprintf(stderr, "Before\n");
  print_full_decomposition(stderr, struc, spec, 4);
#endif
  check(ERR_OK == ppair_decompose_at_node(restree, res, seqs, structree, struc,
      spec, fsplit), "Error decomposing node %i", i_node);

#ifndef NDEBUG
  debug("Leaves");
  print_leaves(stderr, struc, spec, 0);
  print_full_decomposition(stderr, struc, spec, 4);
#endif

  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

static int update_decomposition(
      result_t * results,
      design_state_t * level_states,
      seqstate_t * seqs,
      int reset_level,
      int n_level,
      design_spec_t * spec
    ) {


  result_t * parent_res = results + reset_level;
  result_t * child_res = results + reset_level + 1;

  DBL_TYPE parent_defect = get_defect(parent_res);
  DBL_TYPE child_defect = get_defect(child_res);
  DBL_TYPE accept_defect = 
    (parent_defect - (parent_defect - child_defect / spec->opts.f_stringent) 
                      * spec->opts.f_redecomp) 
      * spec->opts.f_stringent;
  DBL_TYPE temp_defect;
  DBL_TYPE min_defect;
  DBL_TYPE splitting_prob;
  int i_struc, n_strucs, i_node, n_nodes;
  int min_struc, min_node;
  int i_level;
  int i_cutoff = 0;
  n_strucs = spec->n_strucs;

  result_t parent_temp;
  result_t child_temp;
  result_t parent_eval;
  result_t child_eval;

  design_state_t parent_state;
  design_state_t child_state;


  init_result(&parent_eval, spec);
  init_result(&child_eval, spec);

  init_result(&parent_temp, spec);
  init_result(&child_temp, spec);

  copy_result(&parent_eval, results + reset_level);
  copy_result(&child_eval, results + reset_level + 1);

  while (child_defect < accept_defect && i_cutoff < 20) {
    init_states(&parent_state, spec);
    check(ERR_OK == get_decomposition(&parent_state, reset_level, spec),
        "Error decomposing children");
    debug("Parent eval");
    update_result(&parent_eval, &parent_state, seqs, spec);
    copy_result(&parent_temp, &parent_eval);

    init_states(&child_state, spec);
    check(ERR_OK == get_decomposition(&child_state, reset_level + 1, spec),
        "Error decomposing children");
    debug("Child eval");
    update_result(&child_eval, &child_state, seqs, spec);
    copy_result(&child_temp, &child_eval);

    update_result(&parent_temp, &parent_state, seqs, spec);
    update_result(&child_temp, &child_state, seqs, spec);

    min_node = 0;
    min_struc = 0;
    min_defect = 1e10;
    for (i_struc = 0; i_struc < n_strucs ; i_struc++) {
      n_nodes = get_n_leaves(&(parent_temp.strucs[i_struc]));
      for (i_node = 0; i_node < n_nodes; i_node++) {
        debug("Checking merge defect. Struc: %i  node: %i", i_struc, i_node);
        temp_defect = calculate_merge_defect(&parent_temp, &child_temp, i_struc, i_node, spec);
        // debug("Merged defect %i %i = %Lf  child = %Lf", i_struc, i_node, temp_defect, child_defect);
        check(!(temp_defect < 0),"Error merging defect");
        if (temp_defect < min_defect) {
          debug("Setting best defect: %Lf", temp_defect);
          min_defect = temp_defect;
          min_struc = i_struc;
          min_node = i_node;
        }
      }
    }

    debug("min node: %i min struc: %i min defect: %Lf", min_node, min_struc, min_defect);

    splitting_prob = spec->opts.f_split;

    debug("Current level: %i", reset_level);
    check(ERR_OK == ppair_structure_at_leaf(
          &(parent_temp.strucs[min_struc]),
          &(spec->strucs[min_struc]), 
          seqs, spec, min_node, splitting_prob),
        "error performing ppair decomposition on structure %i", i_struc);

    free_states(&child_state);

    init_states(&child_state, spec);

    get_decomposition(&child_state, reset_level + 1, spec);
    update_result(&child_eval, &child_state, seqs, spec);

    child_defect = get_defect(&child_eval);
    debug("Child defect: %Lf / %Lf = %Lf | %Lf", child_defect, 
        accept_defect, child_defect / accept_defect, parent_defect);

    free_states(&parent_state);
    free_states(&child_state);

    free_result(&parent_temp);
    free_result(&child_temp);
    i_cutoff ++;
  }

  free_result(&parent_eval);
  free_result(&child_eval);

  for (i_level = reset_level + 1; i_level < n_level; i_level++) {
    free_states(level_states + i_level);
    init_states(level_states + i_level, spec);
    get_decomposition(level_states + i_level, 
        i_level, spec);

    free_result(results + i_level);
    init_result(results + i_level, spec);
  }

  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

DBL_TYPE get_defect(result_t * res) {
  return res->total_defect;
}

int reseed_next_leaf(seqstate_t * seqs, result_t * res, 
    design_spec_t * spec) {

  int i_tube = 0;

  int i_struc = 0;
  int c_struc;

  int r_struc = -1;
  int r_node = -1;
  int n_strucs;
  int n_nodes;
  int i_node; 
  int r_tube;
  DBL_TYPE total_weight = 0;
  DBL_TYPE cur_weight = 0;
  DBL_TYPE cur_tot_weight = 0.0;
  DBL_TYPE stop = 0.0;

  n_strucs = res->tubes[i_tube].n_strucs;
  for (i_struc = 0; i_struc < n_strucs; i_struc++) {
    if (res->tubes[i_tube].included[i_struc] 
        && res->tubes[i_tube].target[i_struc]) {
      c_struc = res->tubes[i_tube].included_ind[i_struc];
      n_nodes = get_n_leaves(&(res->strucs[c_struc]));

      for (i_node = 0; i_node < n_nodes; i_node++) {
        total_weight += 
            get_nodal_defect(res, i_tube, i_struc, i_node, spec);
      }
    }
  }
  stop = total_weight * genrand_real1();

  for (i_struc = 0; i_struc < n_strucs; i_struc++) {
    if (res->tubes[i_tube].included[i_struc] 
        && res->tubes[i_tube].target[i_struc]) {
      c_struc = res->tubes[i_tube].included_ind[i_struc];
      n_nodes = get_n_leaves(&(res->strucs[c_struc]));

      for (i_node = 0; i_node < n_nodes; i_node++) {
        cur_weight = get_nodal_defect(res, i_tube, i_struc, i_node, spec);
        cur_tot_weight += cur_weight;
        if (cur_tot_weight > stop) {
          r_tube = i_tube;
          r_struc = i_struc;
          r_node = i_node;
          break;
        }
      }
      if (cur_tot_weight > stop) {
        break;
      }
    }
  }

  if (r_struc >= 0 && r_node >= 0) {
    debug("Reseeding leaf at leaf level: %i %i %i %Le %Le", r_tube, 
        r_struc, r_node, cur_weight, total_weight);
    check(ERR_OK == reseed_leaf(seqs, res, r_tube, r_struc, r_node, spec),
        "Error reseeding leaf");
  }
  return ERR_OK;
error:
  return ERR_OOM;
}


int reseed_leaf(seqstate_t * seqs, result_t * res, 
      int r_tube, int r_struc, int r_node, design_spec_t * spec) {
  int c_struc;

  result_tree_t * node = NULL;

  int * f_cons =  NULL; 
  int * c_cons =  NULL; 
  int * f_seq =   NULL; 
  int * c_seq =   NULL; 
  int * nuc_ids = NULL; 

  check(res, "NULL result pointed to");

  c_struc = res->tubes[r_tube].included_ind[r_struc];
  node = get_result_node(res, c_struc, r_node, NULL, spec);
  f_cons =  (int*) malloc(sizeof(int) * node->n_nucs);
  c_cons =  (int*) malloc(sizeof(int) * node->n_nucs);
  f_seq =   (int*) malloc(sizeof(int) * node->n_nucs);
  c_seq =   (int*) malloc(sizeof(int) * node->n_nucs);
  nuc_ids = (int*) malloc(sizeof(int) * spec->strucs[c_struc].n_nucs);

  check_mem(f_cons);
  check_mem(c_cons);
  check_mem(f_seq);
  check_mem(c_seq);


  int i_nuc, c_nuc, c_id;

  check(ERR_OK == get_nuc_ids(nuc_ids, &(spec->strucs[c_struc]), spec),
    "Error getting nucleotide ids");

  for (i_nuc = 0; i_nuc < node->n_nucs; i_nuc++) {
    c_nuc = node->native_map[i_nuc];
    c_id = nuc_ids[c_nuc];
    f_seq[i_nuc] = BASE_N;
    c_seq[i_nuc] = BASE_N;
    if (node->dummy_flag[i_nuc]) {
      f_cons[i_nuc] = BASE_N;
      c_cons[i_nuc] = BASE_N;
    } else {
      f_cons[i_nuc] = spec->seqs.nucs.nucs[c_id];
      if (spec->seqs.comp_map.nucs[c_id] >= 0) {
        c_cons[i_nuc] = spec->seqs.nucs.nucs[spec->seqs.comp_map.nucs[c_id]];
      } else {
        c_cons[i_nuc] = BASE_N;
      }
    }
  }

  check(ERR_OK == init_constraint_random(f_seq, c_seq, f_cons, c_cons, 
        node->n_nucs, &(spec->opts)),
      "Error reseeding leaf");

  for (i_nuc = 0; i_nuc < node->n_nucs; i_nuc++) {
    c_nuc = node->native_map[i_nuc];
    c_id = nuc_ids[c_nuc];
    if (node->dummy_flag[i_nuc]) {
      seqs->dumspec.nucs[c_id] = f_seq[i_nuc];
      if (spec->seqs.comp_map.nucs[c_id] >= 0) {
        seqs->dumspec.nucs[spec->seqs.comp_map.nucs[c_id]] = c_seq[i_nuc];
      }
    } else {
      seqs->nucspec.nucs[c_id] = f_seq[i_nuc];
      if (spec->seqs.comp_map.nucs[c_id] >= 0) {
        seqs->nucspec.nucs[spec->seqs.comp_map.nucs[c_id]] = c_seq[i_nuc];
      }
    }
  }

  free(f_cons);
  free(c_cons);
  free(f_seq);
  free(c_seq);
  free(nuc_ids);



  return ERR_OK;
error:
  free(f_cons);
  free(c_cons);
  free(f_seq);
  free(c_seq);
  free(nuc_ids);
  return ERR_OOM;
}

int evaluate_undesired(
    result_t * res, 
    seqstate_t * seqs,
    design_spec_t * spec) {

  int * old_seqs_same = NULL;
  int n_strs = spec->seqs.strands.n;
  int m_strs;
  int i_str;
  int c_str;
  int * cur_seq = NULL;
  int n_nucs;
  int i_nuc;
  int j_nuc;
  int strs_eq;

  int n_comps = res->n_orderings;
  int i_comp;
  int maxsize = 0;
  int cursize = 0;
  int n_eval = 0;

  struct timeval starttime;
  struct timeval endtime;
  DBL_TYPE cur_pfunc = 0;
  design_state_t state;
  seqstate_t * old_seqs = res->offtarget_seqs;

  old_seqs_same = (int*) calloc(n_strs, sizeof(int));

  refresh_seqstate(seqs, spec);

  init_states(&state, spec);
  // Construct strands

  for (i_str = 0; i_str < n_strs; i_str++) {
    old_seqs_same[i_str] = 0;
  }


  if (old_seqs) {
    for (i_str = 0; i_str < n_strs; i_str++) {
      n_nucs = old_seqs->strands.seqs[i_str].n;
      if (n_nucs == seqs->strands.seqs[i_str].n) {
        strs_eq = 0;
        for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
          if (seqs->strands.seqs[i_str].nucs[i_nuc] != 
              old_seqs->strands.seqs[i_str].nucs[i_nuc]) {
            strs_eq += 1;
          }
        }
        if (strs_eq == 0) {
          old_seqs_same[i_str] = 1;
        } else {
          old_seqs_same[i_str] = 0;
        }
      } else {
        old_seqs_same[i_str] = 0;
      }
    }
  }

  for (i_comp = 0; i_comp < n_comps; i_comp++) {
    cursize = 0;
    m_strs = res->n_strands[i_comp];
    for (i_str = 0; i_str < m_strs; i_str++) {
      c_str = res->orderings[i_comp][i_str];
      cursize += seqs->strands.seqs[c_str].n + 1;
    }
    if (cursize > maxsize) {
      maxsize = cursize;
    }
  }
  check(maxsize > 0, "maxsize is invalid %i", maxsize);
  cur_seq = (int*) malloc(maxsize * sizeof(int));


  for (i_comp = 0; i_comp < n_comps; i_comp++) {
    // For each undesired complex (ordering in struc_map == -1)
    m_strs = res->n_strands[i_comp];
    j_nuc = 0;
    strs_eq = 1;
    // Construct the full sequence
    for (i_str = 0; i_str < m_strs; i_str++) {
      c_str = res->orderings[i_comp][i_str];
      n_nucs = seqs->strands.seqs[c_str].n;
      for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
        cur_seq[i_nuc + j_nuc] = seqs->strands.seqs[c_str].nucs[i_nuc];
      }

      strs_eq = strs_eq && old_seqs_same[c_str];

      cur_seq[i_nuc + j_nuc] = STRAND_PLUS;
      j_nuc += n_nucs + 1;
    }
    cur_seq[j_nuc - 1] = -1;
    res->order_len[i_comp] = j_nuc - m_strs;

    if ((!strs_eq) && res->struc_map[i_comp] < 0) {
      // Evaluate the partition function
      gettimeofday(&starttime,NULL);

      cur_pfunc = pfuncFull(cur_seq, 3, spec->opts.material, 
          spec->opts.dangle_type, 
          spec->opts.temperature - ZERO_C_IN_KELVIN, 0, 
          spec->opts.sodium, 
          spec->opts.magnesium, spec->opts.use_long_helix);
      n_eval ++;
      cur_pfunc /= spec->symmetry[i_comp];

      gettimeofday(&endtime,NULL);

      res->eval_time[i_comp] = (endtime.tv_sec - starttime.tv_sec) 
        + 1e-6 * (endtime.tv_usec - starttime.tv_usec);

      // Set included = 1
      res->dG[i_comp] = - LOG_FUNC(cur_pfunc);
      res->included[i_comp] = 1;
    }
  }

  if (!res->offtarget_seqs) {
    res->offtarget_seqs = (seqstate_t *) malloc(sizeof(seqstate_t));
    check_mem(res->offtarget_seqs);
    init_seqstate(res->offtarget_seqs, spec);
  }

  copy_seqstate(res->offtarget_seqs, seqs);

  // debug("Evaluated %i undesired structures", n_eval);
  // Update all defects, including desired complex partition functions
  // (those should all be evaluated by this point anyway)

  free_states(&state);

  free(cur_seq);
  free(old_seqs_same);
  return ERR_OK;
error:
  free(cur_seq);
  free(old_seqs_same);
  return ERR_INVALID_STATE;
}

static int get_tubes_satisfied(
    result_t * res, 
    result_t * res_deflated,
    design_spec_t * spec) {

  int n_tubes = res->n_tubes;
  int i_tube;
  
  int satisfied = 1;

  for (i_tube = 0; i_tube < n_tubes; i_tube++) {
    if (res->tubes[i_tube].defect > spec->tubes[i_tube].stop 
        && (1 - 1e-7) * res->tubes[i_tube].defect > res_deflated->tubes[i_tube].defect) {
      satisfied = 0;
    }
  }
  return satisfied;
}

static int add_undesired_ordering(
      design_spec_t * spec, 
      seqstate_t * seqs,
      int i_comp
    ) {
  int n_strs = spec->n_strands[i_comp];
  int i_str;
  int c_str;
  int n_doms = 0;

  int * cur_strs = spec->orderings[i_comp];
  int * dom_struc = NULL;
  struc_state_t state;
  int i_struc = spec->n_strucs;
  int i_dom;
  result_struc_t res_struc;
  init_result_struc(&res_struc);

  for (i_str = 0; i_str < n_strs; i_str++) {
    c_str = cur_strs[i_str];
    n_doms += spec->seqs.strands.specs[c_str].n;
  }
  check(n_doms > 0, "N doms must be > 0 %i", n_doms);

  dom_struc = (int*) malloc(n_doms * sizeof(int));

  for (i_dom = 0; i_dom < n_doms; i_dom++) {
    dom_struc[i_dom] = -1;
  }

  add_structure(spec, "::gen", spec->orderings[i_comp],
      spec->n_strands[i_comp], dom_struc, n_doms);
  spec->strucs[i_struc].modifiable = 1;
  init_struc_state(&state, spec->strucs + i_struc);
  debug("Evaluating structure");
  // Update the structure result
  check(ERR_OK == update_structure_result(&res_struc, &state, seqs, spec),
      "Error updating structure result");

  debug("Decomposing structure");
  ppair_structure(&res_struc, spec->strucs + i_struc, &state, seqs, spec);
#ifndef NDEBUG
  // print_full_decomposition(stderr, spec->strucs + i_struc, spec, 4);
#endif

  free_struc_state(&state);
  free_result_struc(&res_struc);
  free(dom_struc);

  return ERR_OK;
error:
  free_result_struc(&res_struc);
  free(dom_struc);
  return ERR_OOM;
}

static int recalculate_defect(
    result_t * res,
    design_spec_t * spec
    ) {
  int n_tubes, n_strucs, i_tube, i_struc;
  int c_struc;
  n_tubes = res->n_tubes;
  DBL_TYPE deficiency, desired_x, nuc_conc, tube_defect, total_defect;
  int n_nucs;

  calculate_concentrations(res, spec);

  total_defect = 0;
  for (i_tube = 0; i_tube < n_tubes; i_tube++) {
    n_strucs = res->tubes[i_tube].n_strucs;
    tube_defect = 0;
    nuc_conc = 0;
    for (i_struc = 0; i_struc < n_strucs; i_struc++) {
      if (res->tubes[i_tube].included[i_struc]) {
        c_struc = res->tubes[i_tube].included_ind[i_struc];
        deficiency = res->tubes[i_tube].target_x[i_struc] - res->tubes[i_tube].x[i_struc];
        desired_x = res->tubes[i_tube].x[i_struc];
        if (deficiency < 0) {
          deficiency = 0;
          desired_x = res->tubes[i_tube].target_x[i_struc];
        }

        n_nucs = res->strucs[c_struc].n_nucs;

        nuc_conc += n_nucs * res->tubes[i_tube].target_x[i_struc];

        tube_defect += res->strucs[c_struc].defect * desired_x + n_nucs * deficiency;
      }
    }
    res->tubes[i_tube].defect = tube_defect / nuc_conc;
    total_defect += tube_defect / nuc_conc;
  }

  res->total_defect = total_defect;

  return ERR_OK;
}

static int add_undesired_complexes(
    design_spec_t * spec, 
    result_t * res, 
    result_t * res_deflated,
    seqstate_t * seqs) {
  int n_strucs;
  int i_struc;
  int c_struc;

  int n_tubes = res->n_tubes;
  int i_tube;
  
  DBL_TYPE c_conc;
  DBL_TYPE * total_concs = 
    (DBL_TYPE *) malloc(res->tot_n_strands * sizeof(DBL_TYPE));
  DBL_TYPE * max_contribs = 
    (DBL_TYPE *) malloc(res->tot_n_strands * sizeof(DBL_TYPE));
  DBL_TYPE * total_contribs = 
    (DBL_TYPE *) malloc(res->tot_n_strands * sizeof(DBL_TYPE));
  int * max_contrib_inds = 
    (int *) malloc(res->tot_n_strands * sizeof(int));
  DBL_TYPE * cur_contribs = 
    (DBL_TYPE *) malloc(res->tot_n_strands * sizeof(DBL_TYPE));
  DBL_TYPE max_conc = 0;
  DBL_TYPE * orig_deflated = 
    (DBL_TYPE *) malloc(res->n_tubes * sizeof(DBL_TYPE));
  DBL_TYPE goal_defect;

  int max_conc_struc = 0;
  int n_undesired = 0;

  for (i_tube = 0; i_tube < n_tubes; i_tube++) {
    orig_deflated[i_tube] = res_deflated->tubes[i_tube].defect;
  }


  check_mem(total_concs);
  check_mem(total_contribs);
  check_mem(max_contribs);
  check_mem(max_contrib_inds);

  for (i_tube = 0; i_tube < n_tubes; i_tube++) {
    debug("Tube: %i Defect: %Lf / %Lf", i_tube, res->tubes[i_tube].defect, 
        spec->tubes[i_tube].stop);
    debug("Deflated: %i Defect: %Lf / %Lf", i_tube, orig_deflated[0],
        spec->tubes[i_tube].stop);
    if (res->tubes[i_tube].defect > spec->tubes[i_tube].stop
        && (1 - 1e-7) * res->tubes[i_tube].defect > res_deflated->tubes[i_tube].defect) {

      debug("Actual defect: %Lf", res->tubes[i_tube].defect);
      debug("include frac:  %Lf", spec->opts.f_refocus);
      debug("Deflated defect: %Lf", orig_deflated[i_tube]);
      goal_defect = res->tubes[i_tube].defect - 
        (spec->opts.f_refocus * (res->tubes[i_tube].defect - orig_deflated[i_tube]));

      while (goal_defect > res_deflated->tubes[i_tube].defect) {
        n_strucs = res->tubes[i_tube].n_strucs;

        max_conc = 0;
        max_conc_struc = -1;
        for (i_struc = 0; i_struc < n_strucs; i_struc ++) {
          c_struc = res->tubes[i_tube].generated_ind[i_struc];
          if (spec->struc_map[c_struc] == -1 && res_deflated->included[c_struc] == 0) {
            c_conc = res->tubes[i_tube].x[i_struc];
            if (c_conc > max_conc) {
              max_conc = c_conc;
              max_conc_struc = c_struc;
            }
          }
        }
        check(max_conc_struc >= 0, "Invalid structure set to be added");
        res_deflated->dG[max_conc_struc] = res->dG[max_conc_struc];
        res_deflated->included[max_conc_struc] = 1;
        recalculate_defect(res_deflated, spec); 
        debug("Added structure %i", max_conc_struc);
        debug("New Deflated: %i Defect: %Lf : %Lf", i_tube, 
            res_deflated->tubes[i_tube].defect, goal_defect);
      }
    }
  }
  n_undesired = 0;
  for (i_struc = 0; i_struc < res_deflated->n_orderings; i_struc++) {
    if (res_deflated->included[i_struc] && res->struc_map[i_struc] == -1) {
      debug("Adding ordering %i", i_struc);
      add_undesired_ordering(spec, seqs, i_struc);
      n_undesired ++;
    }
  }

  check(n_undesired > 0, "Unsatisfied tubes, but no complexes added");
  free(orig_deflated);
  free(max_contribs);
  free(max_contrib_inds);
  free(total_concs);
  free(total_contribs);
  free(cur_contribs);

  return ERR_OK;
error:

  free(orig_deflated);
  free(max_contribs);
  free(max_contrib_inds);
  free(total_concs);
  free(total_contribs);
  free(cur_contribs);
  return ERR_INVALID_STATE;
}

static DBL_TYPE get_tube_defects(
    result_t * res,
    design_spec_t * spec) {
  DBL_TYPE cur_defect = res->tubes[0].defect;
  (void)spec;
  // if (cur_defect < spec->tubes[0].stop) {
  //   cur_defect = spec->tubes[0].stop;
  // }
  return cur_defect;
}

int optimize_tubes(
      result_t * res,
      seqstate_t * res_seqs,
      design_spec_t * spec
    ) {
  
  result_t active_res;
  result_t current_res;
  result_t temp_res;

  design_state_t tempstate;

  seqstate_t best_seqs;
  seqstate_t active_seqs;
  seqstate_t current_seqs;

  struct timeval starttime;
  struct timeval endtime;
  DBL_TYPE elapsed;

  DBL_TYPE temp_defect = DBL_MAX;
  DBL_TYPE best_defect = DBL_MAX;
  DBL_TYPE tube_time = 0;
  DBL_TYPE add_time = 0;

  struct timeval curtime;
  DBL_TYPE dbltime;

  gettimeofday(&starttime, NULL);
  init_seqstate(&best_seqs, spec);
  check(ERR_OK == copy_seqstate(&best_seqs, res_seqs),
      "Error copying sequences");

  init_seqstate(&active_seqs, spec);
  check(ERR_OK == copy_seqstate(&active_seqs, res_seqs),
      "Error copying sequences");

  init_seqstate(&current_seqs, spec);

  if (spec->opts.print_steps) {
    print_leafplot(NP_STEP_NOOPT, &best_seqs, spec);
  }

  init_result(&active_res, spec);
  init_result(&current_res, spec);
  init_result(&temp_res, spec);

  // Decompose and optimize forest
  check(ERR_OK == optimize_forest(&active_res, &active_seqs, spec),
        "Error optimizing trees");
  check(ERR_OK == copy_result(res, &active_res),
      "Error copying result");

  gettimeofday(&curtime, NULL);
  dbltime = curtime.tv_sec + (1e-6 * curtime.tv_usec);
  init_states(&tempstate, spec);

  // Evaluate tube defect
  check(ERR_OK == evaluate_undesired(res, &active_seqs, spec), 
    "Error evaluating undesired complexes");
  check(ERR_OK == update_result(res, &tempstate, &active_seqs, spec),
      "Error updating result");

  free_states(&tempstate);

  best_defect = get_tube_defects(res, spec);

  gettimeofday(&curtime, NULL);
  dbltime -= curtime.tv_sec + (1e-6 * curtime.tv_usec);
  tube_time -= dbltime;
  debug("Setting best defect %Lf", best_defect);
  // Accept as best defect
  check(ERR_OK == copy_seqstate(&best_seqs, &active_seqs),
      "Error copying sequences");

  // Set the current sequences/result
  check(ERR_OK == copy_result(&current_res, res),
      "Error copying result");

  check(ERR_OK == copy_seqstate(&current_seqs, &active_seqs),
      "Error copying sequences");

  if (spec->opts.print_steps) {
    print_leafplot(NP_STEP_TREEOPT, &best_seqs, spec);
  }

  gettimeofday(&curtime, NULL);
  dbltime = curtime.tv_sec + (1e-6 * curtime.tv_usec);

  while (!get_tubes_satisfied(&current_res, &active_res, spec)
          && (dbltime - spec->opts.start_time) < spec->opts.allowed_opt_time
      ) {

    gettimeofday(&curtime, NULL);
    dbltime = curtime.tv_sec + (1e-6 * curtime.tv_usec);
    check(ERR_OK == copy_dummy(&active_seqs),
        "Error copying true sequences to dummy sequences");
    check(ERR_OK == add_undesired_complexes(spec, &current_res, 
          &active_res, &active_seqs),
          "Error adding undesired complexes");

    gettimeofday(&curtime, NULL);
    dbltime -= curtime.tv_sec + (1e-6 * curtime.tv_usec);
    add_time -= dbltime;

    check(ERR_OK == optimize_forest(&active_res, &active_seqs, spec), 
      "Error optimizing trees");
    check(ERR_OK == copy_result(&temp_res, &active_res),
        "Error copying result");

    gettimeofday(&curtime, NULL);
    dbltime = curtime.tv_sec + (1e-6 * curtime.tv_usec);

    init_states(&tempstate, spec);

    check(ERR_OK == evaluate_undesired(&temp_res, &active_seqs, 
          spec),
        "Error evaluating undesired complexes");
    check(ERR_OK == update_result(&temp_res, &tempstate, 
          &active_seqs, spec),
      "Error updating result");

    free_states(&tempstate);

    gettimeofday(&curtime, NULL);
    dbltime -= curtime.tv_sec + (1e-6 * curtime.tv_usec);
    tube_time -= dbltime;

    temp_defect = get_tube_defects(&temp_res, spec);
    check(ERR_OK == copy_result(&current_res, &temp_res),
        "Error copying result");
    check(ERR_OK == copy_seqstate(&current_seqs, &active_seqs),
        "Error copying sequences");

    if (temp_defect < best_defect) {
      best_defect = temp_defect;
      debug("Setting best defect %Lf", best_defect);
      check(ERR_OK == copy_result(res, &current_res),
          "Error copying result");
      check(ERR_OK == copy_seqstate(&best_seqs, &active_seqs),
          "Error copying sequence spec");
    }

    if (spec->opts.print_steps) {
      print_leafplot(NP_STEP_TREEOPT, &best_seqs, spec);
    }

    gettimeofday(&curtime, NULL);
    dbltime = curtime.tv_sec + (1e-6 * curtime.tv_usec);
  }


  if (spec->opts.print_steps) {
    print_leafplot(NP_STEP_TUBEOPT, &best_seqs, spec);
  }
  debug("Final defect: %Lf", best_defect);

  check(ERR_OK == copy_seqstate(res_seqs, &best_seqs),
      "Error copying sequences");

  free_result(&current_res);
  free_result(&temp_res);
  free_result(&active_res);

  free_seqstate(&current_seqs);
  free_seqstate(&best_seqs);
  free_seqstate(&active_seqs);

  gettimeofday(&endtime,NULL);

  elapsed = endtime.tv_sec - starttime.tv_sec 
    + 1e-6 * (endtime.tv_usec - starttime.tv_usec);
  res->elapsed_time = elapsed;
  refresh_seqstate(res_seqs, spec);

  debug("Add time: %Lf", add_time);
  debug("Tube time: %Lf", tube_time);
  debug("Tot time: %Lf", elapsed);
  // Copy the result into the main result structure.
  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}


int optimize_forest(
      result_t * res,
      seqstate_t * seqs,
      design_spec_t * spec
    ) {

  int n_strucs = spec->n_strucs;
  int i_struc;
  int max_levels =  20;// 4000 / 20 = 200 leaves. 10 levels should be plenty
  int n_levels = 0;
  int i_level;
  int j_level;
  int decomp_success;
  int all_satisfied;
  int err = ERR_OK;

  design_state_t temp_state;
  design_state_t * level_states;

  result_t * level_results;
  result_t * child_results;
  result_t temp_res;
  result_t temp_res2;

  seqstate_t * level_seqs;
  seqstate_t temp_seqs;

  DBL_TYPE * best_defects;
  DBL_TYPE * stop_conditions;
  DBL_TYPE cur_defect;
  int reset_level = 0;

  struct timeval curtime;
  DBL_TYPE dbltime;
  DBL_TYPE elapsetime;
  static DBL_TYPE leaftime = 0;
  static DBL_TYPE mergetime = 0;
  static DBL_TYPE redecomptime = 0;

  if (spec->opts.include_all) {
    max_levels = 1;
  }

  /***************************************
   * Initialize decomposition holders    *
   ***************************************/
  level_states = (design_state_t *) malloc(sizeof(design_state_t) * 
      max_levels);
  level_results = (result_t *) malloc(sizeof(result_t) *
      max_levels);
  child_results = (result_t *) malloc(sizeof(result_t) * 
      max_levels);
  level_seqs = (seqstate_t *) malloc(sizeof(seqstate_t) * 
      max_levels);

  best_defects = (DBL_TYPE *) malloc(sizeof(DBL_TYPE) * max_levels);
  stop_conditions = (DBL_TYPE *) malloc(sizeof(DBL_TYPE) * max_levels);

  /***************************************
   * Initialize decomposition tree
   ***************************************/
  init_seqstate(&temp_seqs, spec);
  init_states(&temp_state, spec);
  init_result(&temp_res, spec);

  stop_conditions[0] = spec->tubes[0].stop;
  for (i_level = 1; i_level < max_levels; i_level++) {
    stop_conditions[i_level] = stop_conditions[i_level - 1] * 
      spec->opts.f_stringent;
  }

  for (i_level = 0; i_level < max_levels; i_level++) {
    // Decompose all the structures
    n_levels ++;
    decomp_success = 0;


    // Copy the original sequences down
    init_seqstate(level_seqs + i_level, spec);
    err = (ERR_OK == copy_seqstate(level_seqs + i_level, seqs)) 
      ? err : ERR_OOM;

    // Initialize structure results / best results
    init_result(level_results + i_level, spec);
    init_result(child_results + i_level, spec);
    best_defects[i_level] = DBL_MAX;


    // Initialize decomposition states and copy last decomposition
    init_states(level_states + i_level, spec);
    copy_states(level_states + i_level, &temp_state);

    fill_out_blank_seqs(level_results + i_level, level_states + i_level,
        spec);

    free_states(&temp_state);

    init_states(&temp_state, spec);
    // Decompose the temp state to see if there are any nodes at the 
    // next level
    check(ERR_OK == get_decomposition(&temp_state, i_level + 1, spec),
        "Error decomposing forest");

    for (i_struc = 0; i_struc < n_strucs; i_struc++) {
      if (level_states[i_level].states[i_struc].n < 
          temp_state.states[i_struc].n) {
        decomp_success = 1;
      }
    }

    if (!decomp_success) {
      break;
    }
  }
  check(err == ERR_OK, "Error initializing trees");

  /***************************************
   * Actually optimize the tree          *
   ***************************************/
  all_satisfied = 0;
  // gettimeofday(&curtime, NULL);
  // dbltime = curtime.tv_sec + (1e-6 * curtime.tv_usec);

  while (!all_satisfied ) {
    all_satisfied = 1;

    gettimeofday(&curtime,NULL);
    elapsetime = curtime.tv_sec + (1e-6 * curtime.tv_usec);
    check(ERR_OK == optimize_leaves(level_results + (n_levels - 1), 
        level_seqs + (n_levels - 1), level_states + (n_levels-1),
        stop_conditions[n_levels - 1], 
        spec), "Error in optimize leaves loop");
    gettimeofday(&curtime,NULL);
    elapsetime -= curtime.tv_sec + (1e-6 * curtime.tv_usec);
    leaftime -= elapsetime;

    elapsetime = curtime.tv_sec + (1e-6 * curtime.tv_usec);

    free_result(&temp_res);
    init_result(&temp_res, spec);
    check(ERR_OK == update_result(&temp_res, level_states + (n_levels - 1),
        level_seqs + (n_levels - 1), spec),
        "error evaluating level");

    debug("Actual leaf defect: %Lf", get_defect(&temp_res));

    best_defects[n_levels - 1] = get_defect(&temp_res);

    // Merge up the decomposition tree
    reset_level = -1;
    for (i_level = n_levels - 2; i_level >= 0; i_level --) {
      // Copy results and sequences up
      check(ERR_OK == copy_seqstate(&temp_seqs, &(level_seqs[i_level + 1])),
          "Error copying sequences");
      copy_result(&temp_res, level_results + i_level);

      // Calculate result for current level
      check(ERR_OK == update_result(&temp_res, 
            level_states + i_level, &temp_seqs, spec),
          "Error updating level %i", i_level);

      cur_defect = get_defect(&temp_res);

      // Keep the best sequence seen so far
      if (cur_defect < best_defects[i_level]) {
        check(ERR_OK == copy_seqstate(&(level_seqs[i_level]), &temp_seqs),
            "Error copying sequences");
        copy_result(level_results + i_level, &temp_res);
        copy_result(child_results + i_level, level_results + i_level + 1);
        debug("Level %i new defect: %Lf", i_level, 
            cur_defect);
        best_defects[i_level] = cur_defect;
      } 

      gettimeofday(&curtime, NULL);
      dbltime = curtime.tv_sec + (1e-6 * curtime.tv_usec);

      if (cur_defect * spec->opts.f_stringent > 
            best_defects[i_level + 1]
          && cur_defect > stop_conditions[i_level]
          && (dbltime - spec->opts.start_time) < spec->opts.allowed_opt_time
          ) {

        gettimeofday(&curtime,NULL);
        elapsetime -= curtime.tv_sec + (1e-6 * curtime.tv_usec);
        mergetime -= elapsetime;

        elapsetime = curtime.tv_sec + (1e-6 * curtime.tv_usec);

        init_result(&temp_res2, spec);
        copy_result(&temp_res2, level_results + i_level);
        copy_result(level_results + i_level, &temp_res);
        // Use the latest sequences to do the re-decomposition

        reset_level = i_level;
        all_satisfied = 0;

        debug("Reoptimizing from: %i", reset_level);
        // Keep using the current sequences
        // check(ERR_OK == copy_seqstate(&temp_seqs, level_seqs + reset_level),
        //     "Error copying sequences");

        check(ERR_OK == update_decomposition(level_results,
              level_states, &temp_seqs,
              reset_level, n_levels, spec),
            "Error updating decomposition");

        copy_result(level_results + i_level, &temp_res2);
        free_result(&temp_res2);

        // Reseed descendants
        for (j_level = reset_level + 1; j_level < n_levels; j_level++) {
          // Copy the reseeded sequences to all descendents
          check(ERR_OK == copy_seqstate(level_seqs + j_level, &temp_seqs),
              "Error copying sequences");
          // Reset the best defects in all children levels
          best_defects[j_level] = DBL_MAX;
        }
        gettimeofday(&curtime,NULL);
        elapsetime -= curtime.tv_sec + (1e-6 * curtime.tv_usec);
        redecomptime -= elapsetime;

        break;
      }
    }
  }

  gettimeofday(&curtime,NULL);
  elapsetime -= curtime.tv_sec + (1e-6 * curtime.tv_usec);
  mergetime -= elapsetime;
  debug("Leaf time : %Lf", leaftime);
  debug("Merge time: %Lf", mergetime);
  debug("Decom time: %Lf", redecomptime);
  debug("Final Level 0 defect: %Lf", get_defect(level_results + 0));
  /***************************************
   * Copy the results over               *
   ***************************************/
  check(ERR_OK == copy_seqstate(seqs, &(level_seqs[0])),
      "Error copying sequences");
  // Copy over the results
  copy_result(res, level_results);

  for (i_level = 0; i_level < n_levels; i_level++) {
    free_states(level_states + i_level);
    free_result(level_results + i_level);
    free_result(child_results + i_level);
    free_seqstate(level_seqs + i_level);
  }
  free(level_seqs);
  free(level_states);
  free(level_results);
  free(child_results);
  free(best_defects);
  free(stop_conditions);

  free_result(&temp_res);
  free_states(&temp_state);
  free_seqstate(&temp_seqs);

  return ERR_OK;
error:
  if (level_states) {
    for (i_level = 0; i_level < n_levels; i_level++) {
      free_states(level_states + i_level);
      free_result(level_results + i_level);
      free_seqstate(level_seqs + i_level);
    }
  }
  free(level_seqs);
  free(level_states);
  free(level_results);
  free(best_defects);

  free_result(&temp_res);
  free_states(&temp_state);
  free_seqstate(&temp_seqs);

  return ERR_OTHER;
}

int optimize_leaves(
      result_t * res,
      seqstate_t * seqs,
      design_state_t * state,
      DBL_TYPE goal_frac,
      design_spec_t * spec
    ) {
  seqstate_t temp_seqs;
  result_t temp_res;
  DBL_TYPE best_defect;
  DBL_TYPE cur_defect;
  int fail_count = 0;
  struct timeval curtime;
  DBL_TYPE dbltime;

  init_seqstate(&temp_seqs, spec);
  init_result(&temp_res, spec);

  check(ERR_OK == update_result(res, state, seqs, spec),
      "Error updating result");
  copy_result(&temp_res, res);

  check(ERR_OK == mutate_leaves(res, seqs, state, 
        goal_frac, spec), "Error optimizing leaves");

  check(ERR_OK == update_result(res, state, seqs, spec),
      "Error updating result");

  best_defect = get_defect(res);

  debug("LO defect: %Lf     %Lf", best_defect, 
      get_defect(res));

  debug("defect fraction: %Lf / %Lf", best_defect, goal_frac);

  gettimeofday(&curtime, NULL);
  dbltime = curtime.tv_sec + (1e-6 * curtime.tv_usec);

  while (fail_count < spec->opts.M_leafopt && best_defect > goal_frac
      && dbltime - spec->opts.start_time < spec->opts.allowed_opt_time) {
    // Reseed
    check(ERR_OK == copy_seqstate(&temp_seqs, seqs),
        "Error copying sequences");
    copy_result(&temp_res, res);

    // debug("Not satisfied. Reseeding a leaf");
    // reseed_next_leaf(&temp_seqs, &temp_res, spec);
    debug("Not satisfied. Perturbing...");
    perturb_sequences(&temp_seqs, &temp_res, spec);

    check(ERR_OK == mutate_leaves(&temp_res, &temp_seqs, state, 
          goal_frac, spec), "Error optimizing leaves");

    check(ERR_OK == update_result(&temp_res, state, &temp_seqs, spec),
        "Error updating result");

    cur_defect = get_defect(&temp_res);
    debug("LO defect: %Lf", cur_defect);

    // Copy over on improved defect. Reset otherwise
    if (cur_defect < best_defect) {
      debug("Accepted %Lf < %Lf", cur_defect, best_defect);
      best_defect = cur_defect;
      check(ERR_OK == copy_seqstate(seqs, &temp_seqs),
          "Error copying sequences");
      copy_result(res, &temp_res);
      fail_count = 0;
    } else {
      debug("Rejected %Lf > %Lf", cur_defect, best_defect);
      check(ERR_OK == copy_seqstate(&temp_seqs, seqs),
          "Error copying sequences");
      copy_result(&temp_res, res);
      fail_count += 1;
    }

    gettimeofday(&curtime, NULL);
    dbltime = curtime.tv_sec + (1e-6 * curtime.tv_usec);
  }
  if (spec->opts.print_steps) {
    print_leafplot(NP_STEP_LEAFOPT, seqs, spec);
  }


  free_seqstate(&temp_seqs);
  free_result(&temp_res);

  return ERR_OK;
error:
  free_seqstate(&temp_seqs);
  free_result(&temp_res);
  return ERR_INVALID_STATE;
}

static int get_n_mutable(
      int * n_mut_positions_p,
      result_t * result,
      design_spec_t * spec
    ) {
  int i_str, n_strucs;
  int i_nuc, n_nucs;
  int n_mut_nucs;
  int * on_targets = NULL;

  n_strucs = spec->n_strucs;
  check(n_strucs == result->n_strucs,
      "result and spec are not the same");
  on_targets = (int*) malloc(n_strucs * sizeof(int));
  check_mem(on_targets);
  n_mut_nucs = 0;
  for (i_str = 0; i_str < n_strucs; i_str++) {
    on_targets[i_str] = 0;
  }

  for (i_str = 0; i_str < result->tubes[0].n_strucs; i_str++) {
    if (result->tubes[0].target[i_str]) {
      on_targets[result->tubes[0].included_ind[i_str]] = 1;
    }
  }

  for (i_str = 0; i_str < n_strucs; i_str++) {
    n_nucs = result->strucs[i_str].n_nucs;
    for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
      if (on_targets[i_str]) {
        n_mut_nucs += 1;
      }
    }
  }
  free(on_targets);
  *n_mut_positions_p = n_mut_nucs;
  return ERR_OK;
error:
  return ERR_OOM;
}

int pick_mutation(
    mutation_t * mut,
    result_t * res,
    seqstate_t * seqs,
    design_spec_t * spec
    ) {
  int nuc_id = -1;
  int dummy = -1;
  check(ERR_OK == pick_mutation_location(&nuc_id, &dummy, res, spec),
      "Error picking mutation location");
  check(ERR_OK == sample_mutation(mut, nuc_id, dummy, seqs, spec),
      "Error generating mutation");
  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}


int perturb_sequences(
      seqstate_t * seqs,
      result_t * res,
      design_spec_t * spec
    ) {

  int * native = (int *) calloc(spec->seqs.nucs.n, sizeof(int));
  int * dummies = (int *) calloc(spec->seqs.nucs.n, sizeof(int));
  int i_mutnuc = 0;

  int n_strucs, i_struc, n_nucs, i_nuc;
  int c_struc;
  int c_id;
  DBL_TYPE tot_weight, cur_tot_weight,
           stop_weight;

  DBL_TYPE target_x, x, deficiency;

  n_strucs = spec->tubes[0].n_strucs;

  mutation_t mut;
  init_mutation(&mut);
  check_mem(native);
  check_mem(dummies);
  check(n_strucs == res->tubes[0].n_strucs,
      "res and spec are out of sync");

  for (i_mutnuc = 0; i_mutnuc < spec->opts.M_reseed; i_mutnuc++) {
    tot_weight = 0;
    for (i_struc = 0; i_struc < n_strucs; i_struc++) {
      if (res->tubes[0].target[i_struc]) {
        c_struc = res->tubes[0].included_ind[i_struc];
        n_nucs = res->strucs[c_struc].n_nucs;
        target_x = res->tubes[0].target_x[i_struc];
        x = res->tubes[0].x[i_struc];
        if (target_x > x) {
          deficiency = target_x - x;
        } else {
          x = target_x;
          deficiency = 0;
        }
 
        for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
          c_id = res->strucs[c_struc].nuc_ids[i_nuc];
          if (!native[c_id]) {
            if (spec->opts.disable_defect_weights) {
              tot_weight += target_x;
            } else {
              tot_weight += deficiency + 
                  x * res->strucs[c_struc].defects[i_nuc];
            }
          }

          if (!dummies[c_id]) {
            if (spec->opts.fake_dummies && res->strucs[c_struc].f_defects[i_nuc] > 0) {
              if (spec->opts.disable_defect_weights) {
                tot_weight += target_x;
              } else {
                tot_weight += deficiency + 
                    x * res->strucs[c_struc].f_defects[i_nuc];
              }
            }
          }
        }
      }
    }

    stop_weight = tot_weight * genrand_real1();

    cur_tot_weight = 0;
    for (i_struc = 0; i_struc < n_strucs; i_struc++) {
      if (res->tubes[0].target[i_struc]) {
        c_struc = res->tubes[0].included_ind[i_struc];
        n_nucs = res->strucs[c_struc].n_nucs;
        target_x = res->tubes[0].target_x[i_struc];
        x = res->tubes[0].x[i_struc];
        if (target_x > x) {
          deficiency = target_x - x;
        } else {
          x = target_x;
          deficiency = 0;
        }
 
        for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
          c_id = res->strucs[c_struc].nuc_ids[i_nuc];

          if (!native[c_id]) {
            if (spec->opts.disable_defect_weights) {
              cur_tot_weight += target_x;
            } else {
              cur_tot_weight += deficiency + 
                x * res->strucs[c_struc].defects[i_nuc];
            }

            if (stop_weight < cur_tot_weight) {
              native[res->strucs[c_struc].nuc_ids[i_nuc]] = 1;
              break;
            }
          }

          if (!dummies[c_id]) {
            if(spec->opts.fake_dummies && res->strucs[c_struc].f_defects[i_nuc] > 0) {
              if (spec->opts.disable_defect_weights) {
                cur_tot_weight += target_x;
              } else {
                cur_tot_weight += deficiency + 
                  x * res->strucs[c_struc].f_defects[i_nuc];
              }

              if (stop_weight < cur_tot_weight) {
                dummies[res->strucs[c_struc].nuc_ids[i_nuc]] = 1;
                break;
              }
            }
          }
        }

        if (stop_weight < cur_tot_weight) {
          break;
        }
      }
    }
  }

  for (i_nuc = 0; i_nuc < spec->seqs.nucs.n; i_nuc++) {
    if (native[i_nuc]) {
      check(ERR_OK == sample_mutation(&mut, i_nuc, 0, seqs, spec),
        "Error generating mutation");
      apply_mutation(seqs, &mut);
    } else if (dummies[i_nuc]) {
      check(ERR_OK == sample_mutation(&mut, i_nuc, 1, seqs, spec),
        "Error generating mutation");
      apply_mutation(seqs, &mut);
    }
  }

  free(native);
  free(dummies);
  return ERR_OK;
error:
  free(native);
  free(dummies);
  return ERR_OOM;
}

int pick_mutation_location(
      int * nuc_id,
      int * dummy,
      result_t * res,
      design_spec_t * spec
    ) {
  int n_strucs, i_struc, n_nucs, i_nuc;
  int c_struc;
  DBL_TYPE tot_weight, cur_tot_weight,
           stop_weight;

  DBL_TYPE target_x, x, deficiency;

  n_strucs = spec->tubes[0].n_strucs;
  check(n_strucs == res->tubes[0].n_strucs,
      "res and spec are out of sync");

  *nuc_id = -1;
  *dummy = -1;

  tot_weight = 0;
  for (i_struc = 0; i_struc < n_strucs; i_struc++) {
    if (res->tubes[0].target[i_struc]) {
      c_struc = res->tubes[0].included_ind[i_struc];
      n_nucs = res->strucs[c_struc].n_nucs;
      target_x = res->tubes[0].target_x[i_struc];
      x = res->tubes[0].x[i_struc];
      if (target_x > x) {
        deficiency = target_x - x;
      } else {
        x = target_x;
        deficiency = 0;
      }
 
      for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
        if (spec->opts.disable_defect_weights) {
          tot_weight += target_x;

          if (spec->opts.fake_dummies && res->strucs[c_struc].f_defects > 0) {
            tot_weight += target_x;
          }

        } else {
          tot_weight += deficiency + 
            x * res->strucs[c_struc].defects[i_nuc];

          if (spec->opts.fake_dummies) {
            tot_weight += deficiency + 
                x * res->strucs[c_struc].f_defects[i_nuc];
          }
        }
      }
    }
  }

  stop_weight = tot_weight * genrand_real1();

  cur_tot_weight = 0;
  for (i_struc = 0; i_struc < n_strucs; i_struc++) {
    if (res->tubes[0].target[i_struc]) {
      c_struc = res->tubes[0].included_ind[i_struc];
      n_nucs = res->strucs[c_struc].n_nucs;
      target_x = res->tubes[0].target_x[i_struc];
      x = res->tubes[0].x[i_struc];
      if (target_x > x) {
        deficiency = target_x - x;
      } else {
        x = target_x;
        deficiency = 0;
      }
 
      for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {

        if (spec->opts.disable_defect_weights) {
          cur_tot_weight += target_x;
        } else {
          cur_tot_weight += deficiency + 
            x * res->strucs[c_struc].defects[i_nuc];
        }

        if (stop_weight < cur_tot_weight) {
          *nuc_id = res->strucs[c_struc].nuc_ids[i_nuc];
          *dummy = 0;
          break;
        }

        if(spec->opts.fake_dummies && 
              res->strucs[c_struc].f_defects[i_nuc] > 0) {
          if (spec->opts.disable_defect_weights) {
            cur_tot_weight += target_x;
          } else {
            cur_tot_weight += deficiency + 
              x * res->strucs[c_struc].f_defects[i_nuc];
          }

          if (stop_weight < cur_tot_weight) {
            *nuc_id = res->strucs[c_struc].nuc_ids[i_nuc];
            *dummy = 1;
            break;
          }
        }
      }

      if (stop_weight < cur_tot_weight) {
        break;
      }
    }
  }

  if (*nuc_id == -1 && *dummy == -1) {
    *nuc_id = 0;
    dummy = 0;
  }

  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

int mutate_leaves(
      result_t * res,
      seqstate_t * seqs,
      design_state_t * state,
      DBL_TYPE goal_frac,
      design_spec_t * spec
    ) {

  seqstate_t temp_seqs;
  result_t temp_res;
  DBL_TYPE best_defect ;
  DBL_TYPE cur_defect = 0;
  DBL_TYPE dbltime;
  struct timeval curtime;
  int mut_failed = 0;
  int n_mutable_nucs = 0;
  int n_failed = 0;
  mutation_t mut;             // xi
  mutation_list_t muts_tried; // gamma

  init_seqstate(&temp_seqs, spec);
  check(ERR_OK == copy_seqstate(&temp_seqs, seqs),
      "Error copying sequences");

  check(ERR_OK == init_mutation(&mut), 
      "Error initializing mutation");
  check(ERR_OK == init_mutation_list(&muts_tried), 
      "Error initializing mutation list");

  init_result(&temp_res, spec);
  copy_result(&temp_res, res);

  // Compute the initial defect
  check(ERR_OK == update_result(res, state, seqs, spec),
      "Error in update_structure_result");

  best_defect = get_defect(res);
  check(ERR_OK == get_n_mutable(&n_mutable_nucs, res, spec),
      "Error getting mutable count");

  gettimeofday(&curtime, NULL);
  dbltime = curtime.tv_sec + (1e-6 * curtime.tv_usec);
  while (n_failed < spec->opts.M_unfavorable 
      && best_defect > goal_frac 
      && dbltime - spec->opts.start_time < spec->opts.allowed_opt_time) {
    mut_failed = 0;

    check(ERR_OK == pick_mutation(&mut, res, &temp_seqs, spec),
        "Error picking mutation");

    // Check if the mutation has been tried before
    if (mut.n == 0 || -1 != find_mutation(&muts_tried, &mut)) {
      mut_failed = 1;
    } else {
      // Evaluate with the new mutation
      apply_mutation(&temp_seqs, &mut);
      check(ERR_OK == update_result(&temp_res, state, &temp_seqs, spec),
          "Error updating result");

      cur_defect = get_defect(&temp_res);
      // Get the new defect
      if (cur_defect < best_defect) {
        // If the defect is improved, save the sequences/properties

        check(ERR_OK == copy_seqstate(seqs, &temp_seqs),
            "Error copying sequences");

        copy_result(res, &temp_res);
        best_defect = cur_defect;
        debug("Success: %Lf %4i / %4i", best_defect, n_failed, 
            (int)(spec->opts.M_unfavorable));

        clear_mutation_list(&muts_tried);
        check(ERR_OK == get_n_mutable(&n_mutable_nucs, res, spec),
            "Error getting mutable count");
        
      } else {
        // Mutation failed. Copy back the old properties.
        check(ERR_OK == copy_seqstate(&temp_seqs, seqs),
            "Error copying sequences");
        copy_result(&temp_res, res);
        append_mutation(&muts_tried, &mut);
        mut_failed = 1;
      }
    }

    if (mut_failed) {
      n_failed += 1;
    } else {
      n_failed = 0;
    }
    gettimeofday(&curtime, NULL);
    dbltime = curtime.tv_sec + (1e-6 * curtime.tv_usec);
  }

  free_mutation(&mut);
  free_mutation_list(&muts_tried);
  free_result(&temp_res);
  free_seqstate(&temp_seqs);
  return ERR_OK;
error:
  free_mutation(&mut);
  free_mutation_list(&muts_tried);
  free_result(&temp_res);
  free_seqstate(&temp_seqs);
  return ERR_OTHER;
}

int print_leafplot(
      enum NUPACK_DESIGN_STEP step,
      seqstate_t * seqs,
      design_spec_t * spec
    ) {


  static char token[NP_STEP_COUNT][20] = {
    "noopt__",
    "leafopt",
    "treeopt",
    "tubeopt"
  };
  static int counts[NP_STEP_COUNT] = {
    0, 
    0,
    0,
    0
  };
  
  char * cur_tok = token[step];
  char * fn = NULL;
  int n_char = 0;
  int fn_len = 0;
  int tok_len = 0;
  int cur_count = counts[step];

  FILE * outf = NULL;
  result_t temp_res;
  design_state_t temp_state;

  refresh_seqstate(seqs, spec);
  counts[step] ++;
  if (counts[step] < spec->opts.max_print_steps) {

    tok_len = strlen(cur_tok);
    if (spec->opts.file_prefix) {
      fn_len = strlen(spec->opts.file_prefix);
    }

    n_char += fn_len + tok_len + 20;

    fn = (char *) malloc(n_char * sizeof(char));
    fn[0] = '\0';

    if (spec->opts.file_prefix) {
      snprintf(fn , fn_len + tok_len + 20, "%s_%s_%03i.strands", 
          spec->opts.file_prefix, cur_tok, cur_count);
    } else {
      snprintf(fn, tok_len + 20, "%s_%03i.strands", cur_tok, cur_count);
    }

    debug("Printing %s", fn);
    outf = fopen(fn, "w");
    print_strands(outf, seqs, 0, spec);
    fclose(outf);

    if (spec->opts.print_leaves) {
      if (spec->opts.file_prefix) {
        snprintf(fn , fn_len + tok_len + 20, "%s_%s_%03i_root.npo", 
            spec->opts.file_prefix, cur_tok, cur_count);
      } else {
        snprintf(fn, tok_len + 20, "%s_%03i_root.npo", cur_tok, cur_count);
      }

      debug("Evaluating for %s", fn);
      outf = fopen(fn, "w");
      
      init_result(&temp_res, spec);
      init_states(&temp_state, spec);

      get_decomposition(&temp_state, 0, spec);

      evaluate_undesired(&temp_res, seqs, spec);
      update_result(&temp_res, &temp_state, seqs, spec);

      print_design_result(outf, &temp_res, seqs, 0, spec, "tubedesign steps");
      fclose(outf);

      free_states(&temp_state);
      free_result(&temp_res);

      if (spec->opts.file_prefix) {
        snprintf(fn , fn_len + tok_len + 20, "%s_%s_%03i_leaf.npo", 
            spec->opts.file_prefix, cur_tok, cur_count);
      } else {
        snprintf(fn, tok_len + 20, "%s_%03i_leaf.npo", cur_tok, cur_count);
      }

      debug("Evaluating for %s", fn);
      outf = fopen(fn, "w");
 
      init_result(&temp_res, spec);
      init_states(&temp_state, spec);

      get_decomposition(&temp_state, 20, spec);

      update_result(&temp_res, &temp_state, seqs, spec);

      print_design_result(outf, &temp_res, seqs, 0, spec, "tubedesign steps");
      fclose(outf);

      free_states(&temp_state);
      free_result(&temp_res);
    }

    free(fn);
  }
  return ERR_OK;
}

