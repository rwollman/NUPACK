#include "pathway_design.h"

/***********************************************************/
int convert_nucs_to_str(char * str, 
      int * nucs, 
      int n_nucs, 
      options_t * opts
    ) {
  char nucs_rna[] = "NACGURMSWKYVHDB+";
  char nucs_dna[] = "NACGTRMSWKYVHDB+";
  int i_nuc;
  if (opts->material == DNA) {
    for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
      check(nucs[i_nuc] <= 15, "Invalid nucleotide %i encountered", 
          nucs[i_nuc]);
      str[i_nuc] = nucs_dna[nucs[i_nuc]];
    }
  } else {
    for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
      check(nucs[i_nuc] <= 15, "Invalid nucleotide %i encountered", 
          nucs[i_nuc]);
      str[i_nuc] = nucs_rna[nucs[i_nuc]];
    }
  }
  str[n_nucs] = '\0';
  return ERR_OK;
error:
  return ERR_INVALID_INPUT;
}
/***********************************************************/

/***********************************************************/
int convert_str_to_nucs(int * nucs, char * str, int n_nucs) {
  int i_nuc;
  int nuc;
  for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    switch (str[i_nuc]) {
      case 'N':
        nuc = BASE_N;
        break;
      case 'A':
        nuc = BASE_A;
        break;
      case 'C':
        nuc = BASE_C;
        break;
      case 'G':
        nuc = BASE_G;
        break;
      case 'T':
        nuc = BASE_U;
        break;
      case 'U':
        nuc = BASE_U;
        break;
      case 'R':
        nuc = BASE_R;
        break;
      case 'Y':
        nuc = BASE_Y;
        break;
      case 'M':
        nuc = BASE_M;
        break;
      case 'K':
        nuc = BASE_K;
        break;
      case 'S':
        nuc = BASE_S;
        break;
      case 'W':
        nuc = BASE_W;
        break;
      case 'V':
        nuc = BASE_V;
        break;
      case 'B':
        nuc = BASE_B;
        break;
      case 'H':
        nuc = BASE_H;
        break;
      case 'D':
        nuc = BASE_D;
        break;
      default:
        sentinel("Invalid nucleotide encountered");
        break;
    }
    nucs[i_nuc] = nuc;
  }
  return ERR_OK;
error:
  return ERR_INVALID_INPUT;
}
/***********************************************************/

/***********************************************************/
static int print_tree(FILE * f, struc_tree_t * tree, 
    design_struc_t * struc, int * ind, design_spec_t * spec, int indent,
    int print_interior) {
  int i_seg, n_segs;
  int i_nuc, j_nuc, k_nuc;

  int * breaks = (int*) malloc(struc->n_strands * sizeof(int));

  int i_br, n_br;
  int i_c;
  int i_struc;
  int n_print = 0;

  char * struc_str = (char*) malloc(struc->n_nucs + tree->n_segments + 
      struc->n_strands);

  char * indstr = (char*) malloc((indent + 1) * sizeof(char));

  check_mem(indstr);
  check_mem(struc_str);
  check_mem(breaks);

  if (print_interior || tree->n_children == 0) {

    for (i_c = 0; i_c < indent; i_c++) {
      indstr[i_c] = ' ';
    }
    indstr[indent] = '\0';

    design_struc_t * cur_struc = struc;

    n_br = struc->n_strands - 1;

    i_nuc = 0;
    for (i_br = 0; i_br < n_br; i_br++) {
      i_nuc += spec->seqs.strands.specs[struc->strands[i_br]].n;
      breaks[i_br] = i_nuc;
    }

    n_segs = tree->n_segments;

    i_struc = 0;
    while (cur_struc) {
      j_nuc = 0;
      k_nuc = 0;
      if (tree_is_decomp_for(tree, cur_struc)) {
        for (i_seg = 0; i_seg < n_segs; i_seg++) {
          for ( ; k_nuc < tree->seg_start[i_seg]; k_nuc++) {
            struc_str[j_nuc] = ' ';
            j_nuc++;
          }
          for (i_nuc = tree->seg_start[i_seg]; i_nuc < tree->seg_stop[i_seg]; 
              i_nuc++) {
            while (i_br < n_br && breaks[i_br] < i_nuc) {
              i_br ++;
              struc_str[j_nuc] = ' ';
              j_nuc ++;
            }

            if (i_br < n_br && breaks[i_br] == i_nuc) {
              struc_str[j_nuc] = '+';
              j_nuc ++;
              i_br ++;
            }
            if (cur_struc->struc[i_nuc] < 0) {
              struc_str[j_nuc] = '.';
            } else if (cur_struc->struc[i_nuc] > i_nuc) {
              struc_str[j_nuc] = '(';
            } else if (cur_struc->struc[i_nuc] < i_nuc) {
              struc_str[j_nuc] = ')';
            } else {
              sentinel("Invalid structure present");
            }
            j_nuc ++;
            k_nuc ++;
          }
        }
        struc_str[j_nuc] = '\0';

        fprintf(f, "%4i %3i %3i %s\n", *ind, tree->n_segments, i_struc, struc_str);
        n_print ++;
      }
      i_struc ++;
      cur_struc = cur_struc->next;
    }

    if (n_print == 0) {
      for (i_seg = 0; i_seg < n_segs; i_seg++) {
        for ( ; k_nuc < tree->seg_start[i_seg]; k_nuc++) {
          struc_str[j_nuc] = ' ';
          j_nuc++;
        }
        for (i_nuc = tree->seg_start[i_seg]; i_nuc < tree->seg_stop[i_seg]; 
            i_nuc++) {
          while (i_br < n_br && breaks[i_br] < i_nuc) {
            i_br ++;
            struc_str[j_nuc] = ' ';
            j_nuc ++;
          }

          if (i_br < n_br && breaks[i_br] == i_nuc) {
            struc_str[j_nuc] = '+';
            j_nuc ++;
            i_br ++;
          }
          struc_str[j_nuc] = '-';

          j_nuc ++;
          k_nuc ++;
        }
      }
      struc_str[j_nuc] = '\0';

      fprintf(f, "%4i %3i %3i %s\n", *ind, tree->n_segments, i_struc, struc_str);
     
    }

    (*ind) ++;
  }

  for (i_c = 0; i_c < tree->n_children; i_c++) {
    print_tree(f, tree->children[i_c], struc, ind, spec, indent, print_interior);
  }

  free(breaks);
  free(struc_str);
  free(indstr);
  return ERR_OK;
error:
  free(breaks);
  free(struc_str);
  free(indstr);
  return ERR_INVALID_STATE;
}
/***********************************************************/

int print_leaves(FILE * f, 
    design_struc_t * struc,
    design_spec_t * spec,
    int indent) {
  int node_id = 0;
  print_tree(f, struc->tree, struc, &node_id, spec, indent, 0);
  return ERR_OK;
}

int print_full_decomposition(FILE * f, 
    design_struc_t * struc,
    design_spec_t * spec,
    int indent) {
  int node_id = 0;
  print_tree(f, struc->tree, struc, &node_id, spec, indent, 1);
  return ERR_OK;
}

int print_decompositions(FILE * f,
    design_state_t * states,
    design_spec_t * spec,
    int indent) {
  int i_struc;
  int node_id;
  for (i_struc = 0; i_struc < spec->n_strucs; i_struc++) {
    node_id = 0;
    print_tree(f, states->states[i_struc].tree, states->states[i_struc].struc, 
        &node_id, spec, indent, 1);
  }
  return ERR_OK;
}


/***********************************************************/
int print_decomposition_state(FILE * f, 
    struc_state_t * state, 
    design_spec_t * spec,
    int indent) {

  int node_id = 0;
  check(state->tree != NULL, "Printing null decomposition state");

  print_tree(f, state->tree, state->struc, &node_id, spec, indent, 1);

  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}
/***********************************************************/

/***********************************************************/
int print_design_result(
    FILE * f,
    result_t * res,
    seqstate_t * seqs,
    int indent,
    design_spec_t * spec,
    char * executable) {

  char * istr = (char *) malloc(indent + 1);
  options_t * opts = &(spec->opts);
  memset(istr, ' ', indent);
  istr[indent] = '\0';

  fprintf(f, "%sNUPACK version: %s\n", istr, NUPACK_VERSION);
  fprintf(f, "%sBackend: %s\n", istr, executable);
  fprintf(f, "%sDesign index: %i\n", istr, res->id);
  fprintf(f, "%sOptions: \n", istr);
  // Physical options
  fprintf(f, "%s    temperature: %lf\n", istr, (double) opts->temperature);
  if (opts->material == RNA) {
    fprintf(f, "%s    material: rna1995\n", istr);
  } else if (opts->material == RNA37) {
    fprintf(f, "%s    material: rna1999\n", istr);
  } else if (opts->material == DNA) {
    fprintf(f, "%s    material: dna1998\n", istr);
  } else if (opts->material == USE_SPECIFIED_PARAMETERS_FILE) {
    fprintf(f, "%s    material: %s\n", istr, PARAM_FILE);
  } else {
    fprintf(f, "%s    material: undefined\n", istr);
  }
  // DNA1998/RNA1995/RNA1999 type
  fprintf(f, "%s    sodium: %lf\n", istr, (double) opts->sodium);      
  fprintf(f, "%s    magnesium: %lf\n", istr, (double) opts->magnesium);   
  fprintf(f, "%s    dangle_type: %i\n", istr,  opts->dangle_type);
  fprintf(f, "%s    use_long_helix: %i\n", istr,  opts->use_long_helix);

  if (opts->designing) {
    // Optimization options
    fprintf(f, "%s    seed: %u\n", istr, opts->seed);
    fprintf(f, "%s    H_split: %i\n", istr, opts->H_split);
    fprintf(f, "%s    N_split: %i\n", istr,  opts->N_split);
    fprintf(f, "%s    M_bad: %Lf\n", istr,  opts->M_unfavorable);
    fprintf(f, "%s    M_reopt: %Lf\n", istr,  opts->M_leafopt);
    fprintf(f, "%s    M_reseed: %Lf\n", istr,  opts->M_reseed);
    fprintf(f, "%s    f_split: %Lf\n", istr, opts->f_split);
    fprintf(f, "%s    f_passive: %Lf\n", istr, opts->f_passive);
    fprintf(f, "%s    f_stringent: %Lf\n", istr, opts->f_stringent);
    fprintf(f, "%s    f_refocus: %Lf\n", istr, opts->f_refocus);
    fprintf(f, "%s    f_redecomp: %Lf\n", istr, opts->f_redecomp);
    fprintf(f, "%s    initgc: %Lf\n", istr, opts->gc_init_prob);
    fprintf(f, "%s    dgclamp: %Lf\n", istr, opts->bonus_per_split);
  }

  // Optimization modes
  fprintf(f, "%s    allow_wobble: %i\n", istr,  opts->allow_wobble);
#ifndef NDEBUG
  fprintf(f, "%s    allow_mismatch: %i\n", istr,  opts->allow_mismatch);
  fprintf(f, "%s    disable_mutation_weights: %i\n", istr, opts->disable_defect_weights);
  fprintf(f, "%s    include_all: %i\n", istr, opts->include_all);

  // Result options
  fprintf(f, "%s    min_ppair_saved: %lf\n", istr, (double) opts->min_ppair_saved);  
  fprintf(f, "%s    output_format: out\n", istr);
#endif

  if (opts->designing) {
    fprintf(f, "%sRoot time: %Lf\n", istr, res->root_time);
  }
  fprintf(f, "%sElapsed time: %Lf\n", istr, res->elapsed_time);
  // fprintf(f, "%sTotal defect: %lf\n", istr, (double) res->total_defect);
  // fprintf(f, "%sNormalized defect: %lf\n", 
  //     istr, (double) res->total_defect / res->n_tubes);

  print_domains(f, seqs, indent, spec);
  print_strands(f, seqs, indent, spec);
  print_objectives(f, res, indent, spec);
  free(istr);
  return ERR_OK;
}
/***********************************************************/

/***********************************************************/
int print_domains(
    FILE * f, 
    seqstate_t * seqs, 
    int indent, 
    design_spec_t * spec) {
  int i_dom;
  char * istr = (char *) malloc(indent + 1);
  char * dom_str = (char *) malloc(1000 * sizeof(char));

  memset(istr, ' ', indent);
  istr[indent] = '\0';

  fprintf(f, "%sDomains:\n", istr);
  for (i_dom = 0; i_dom < seqs->domains.n; i_dom++) {
    convert_nucs_to_str(dom_str, seqs->domains.seqs[i_dom].nucs, 
        seqs->domains.seqs[i_dom].n, &(spec->opts));

    fprintf(f, "%s  - Name: %s\n", istr, 
        spec->seqs.domains.names[i_dom]);
    fprintf(f, "%s    Sequence: %s\n", istr, dom_str);
  }
  free(dom_str);
  free(istr);
  return ERR_OK;
}
/***********************************************************/


/***********************************************************/
int print_strands(
    FILE * f, 
    seqstate_t * seqs, 
    int indent, 
    design_spec_t * spec) {

  int i_str;
  char * istr = (char *) malloc(indent + 1);
  char * str_str = (char *) malloc(1000 * sizeof(char));

  memset(istr, ' ', indent);
  istr[indent] = '\0';

  refresh_seqstate(seqs, spec);

  fprintf(f, "%sStrands:\n", istr);
  for (i_str = 0; i_str < seqs->strands.n; i_str++) {
    convert_nucs_to_str(str_str, seqs->strands.seqs[i_str].nucs, 
        seqs->strands.seqs[i_str].n, &(spec->opts));

    if (spec->seqs.strands.names[i_str][0] == '-') {
      fprintf(f, "%s  - Name: %s\n", istr, 
          spec->seqs.strands.names[i_str] + 1);
    } else {
      fprintf(f, "%s  - Name: %s\n", istr, 
          spec->seqs.strands.names[i_str]);
    }
    fprintf(f, "%s    Sequence: %s\n", istr, str_str);
  }
  free(str_str);
  free(istr);

  return ERR_OK;
}
/***********************************************************/

/***********************************************************/
int print_objectives(
    FILE * f,
    result_t * res,
    int indent,
    design_spec_t * spec) {

  int n_tubes = res->n_tubes;
  int i_tube;
  int n_strucs;
  int i_struc;
  int c_struc;
  char * ind_str = (char *) malloc(indent + 1);
  char * structure = NULL;
  char * sequence = NULL;
  int n_nucs;
  int n_breaks;
  int n_chars = 0;
  int cap_chars = 0;
  int g_struc;
  int i_str;
#ifndef NDEBUG
  design_struc_t * cur_struc;
#endif

  DBL_TYPE pfunc;
  DBL_TYPE dG;
  DBL_TYPE conc;
  DBL_TYPE t_conc;
  DBL_TYPE defect;
  DBL_TYPE water_conc = water_density(
      (double)(spec->opts.temperature - ZERO_C_IN_KELVIN));
  DBL_TYPE * s_conc = 
    (DBL_TYPE *) malloc(spec->seqs.strands.n * sizeof(DBL_TYPE));
  DBL_TYPE * i_conc = 
    (DBL_TYPE *) malloc(spec->seqs.strands.n * sizeof(DBL_TYPE));
  DBL_TYPE * u_conc = 
    (DBL_TYPE *) malloc(spec->seqs.strands.n * sizeof(DBL_TYPE));

  memset(ind_str, ' ', indent);
  ind_str[indent] = '\0';

  fprintf(f, "%sObjectives:\n", ind_str);
  for (i_tube = 0; i_tube < n_tubes; i_tube++) {
    for (i_str = 0; i_str < spec->seqs.strands.n; i_str++) {
      s_conc[i_str] = 0;
      i_conc[i_str] = 0;
      u_conc[i_str] = 0;
    }
    fprintf(f, "%s  - Name: %s\n", ind_str, spec->tube_names[i_tube]);

    fprintf(f, "%s    Normalized defect: %Lf\n", ind_str, res->tubes[i_tube].defect);
    if (spec->opts.designing) {
      fprintf(f, "%s    Stop fraction: %Lf\n", ind_str, spec->tubes[i_tube].stop);
    }
    fprintf(f, "%s    Max size: %i\n", ind_str, spec->tubes[i_tube].maxsize);
    fprintf(f, "%s    Structures:\n", ind_str);
    n_strucs = res->tubes[i_tube].n_strucs;
    for (i_struc = 0; i_struc < n_strucs; i_struc++) {
      if (res->tubes[i_tube].included[i_struc]) {
        c_struc = res->tubes[i_tube].included_ind[i_struc];
        g_struc = res->tubes[i_tube].generated_ind[i_struc];
        n_nucs = res->strucs[c_struc].n_nucs;
        n_breaks = res->strucs[c_struc].n_breaks;
        n_chars = n_nucs + n_breaks + 1;
        if (n_chars > cap_chars) {
          cap_chars = n_chars * 2;
          structure = (char *)realloc(structure, cap_chars * sizeof(char));
          sequence = (char *)realloc(sequence, cap_chars * sizeof(char));
        }
        get_dpp(structure, res->strucs + c_struc);
        get_sequence_plus(sequence, res->strucs + c_struc, &(spec->opts));
        pfunc = res->strucs[c_struc].pfunc;
        dG = -kB * spec->opts.temperature * LOG_FUNC(pfunc);
        conc = res->tubes[i_tube].x[i_struc] * water_conc;
        t_conc = spec->tubes[i_tube].target_x[i_struc] * water_conc;
        defect = res->strucs[c_struc].defect;
        for (i_str = 0; i_str < spec->n_strands[g_struc]; i_str++) {
          s_conc[spec->orderings[g_struc][i_str]] += conc;
        }
        for (i_str = 0; i_str < spec->n_strands[g_struc]; i_str++) {
          i_conc[spec->orderings[g_struc][i_str]] += conc;
        }

        if (res->tubes[i_tube].target[i_struc]) {
          fprintf(f, "%s      - Name: %s\n", 
              ind_str, spec->struc_names[c_struc]);
          fprintf(f, "%s        Structure: %s\n", ind_str, structure);
        } else {
          if (spec->seqs.strands.names[spec->orderings[g_struc][0]][0] == '-'){
            fprintf(f, "%s      - Name: gen:%s", ind_str, 
                spec->seqs.strands.names[spec->orderings[g_struc][0]]+1);
            for (i_str = 1; i_str < spec->n_strands[g_struc]; i_str++) {
              fprintf(f, "%s", 
                  spec->seqs.strands.names[spec->orderings[g_struc][i_str]]);
            }
            fprintf(f, "\n");
          } else {
            fprintf(f, "%s      - Name: gen:%s", ind_str, 
                spec->seqs.strands.names[spec->orderings[g_struc][0]]);
            for (i_str = 1; i_str < spec->n_strands[g_struc]; i_str++) {
              fprintf(f, "-%s", 
                  spec->seqs.strands.names[spec->orderings[g_struc][i_str]]);
            }
            fprintf(f, "\n");
          }
        }
        fprintf(f, "%s        Sequence : %s\n", ind_str, sequence);
#ifndef NDEBUG
        fprintf(f, "%s        Childcounts  : ", ind_str);
        print_child_counts(f, res->strucs[c_struc].tree);
        fprintf(f, "\n");
        cur_struc = spec->strucs[c_struc].next;
        while (cur_struc) {
          get_spec_dpp(structure, cur_struc, spec);
          fprintf(f, "%s        Altstruc : %s\n", ind_str, structure);
          cur_struc = cur_struc->next;
        }
        fprintf(f, "%s        I struc    : %i\n", ind_str, i_struc);
        fprintf(f, "%s        C struc    : %i\n", ind_str, c_struc);
        fprintf(f, "%s        G struc    : %i\n", ind_str, g_struc);
        fprintf(f, "%s        Eval time: %Lf\n", ind_str, res->eval_time[g_struc]);
#endif
        fprintf(f, "%s        Free energy: %Lf\n", ind_str, dG);
        fprintf(f, "%s        Defect: %Lf\n", ind_str, defect);
        fprintf(f, "%s        Normalized defect: %Lf\n", 
            ind_str, defect / n_nucs);
        fprintf(f, "%s        Concentration: %Le\n", ind_str, conc);
        fprintf(f, "%s        Target concentration: %Le\n", ind_str, t_conc);
 
      } else {
#ifndef NDEBUG
        g_struc = res->tubes[i_tube].generated_ind[i_struc];
        dG = kB * spec->opts.temperature * res->dG[g_struc];
        conc = res->tubes[i_tube].x[i_struc] * water_conc;
        t_conc = spec->tubes[i_tube].target_x[i_struc] * water_conc;
        for (i_str = 0; i_str < spec->n_strands[g_struc]; i_str++) {
          s_conc[spec->orderings[g_struc][i_str]] += conc;
        }
        for (i_str = 0; i_str < spec->n_strands[g_struc]; i_str++) {
          u_conc[spec->orderings[g_struc][i_str]] += conc;
        }

        fprintf(f, "%s      - Name: gen:%s", ind_str,
                spec->seqs.strands.names[spec->orderings[g_struc][0]]);
        for (i_str = 1; i_str < spec->n_strands[g_struc]; i_str++) {
          fprintf(f, "-%s", 
                spec->seqs.strands.names[spec->orderings[g_struc][i_str]]);
        }
        fprintf(f, "\n");
        fprintf(f, "%s        I struc    : %i\n", ind_str, i_struc);
        fprintf(f, "%s        G struc    : %i\n", ind_str, g_struc);
        fprintf(f, "%s        Eval time: %Lf\n", ind_str, res->eval_time[g_struc]);
        fprintf(f, "%s        Free energy: %Lf\n", ind_str, dG);
        fprintf(f, "%s        Concentration: %Le\n", ind_str, conc);
        fprintf(f, "%s        Target concentration: %Le\n", ind_str, t_conc);
#endif
      }
    }
  }
  free(s_conc);
  free(u_conc);
  free(i_conc);

  free(ind_str);
  free(structure);
  free(sequence);
  return ERR_OK;
}
/***********************************************************/

/***********************************************************/
int print_struc_result(
      FILE * f, 
      result_struc_t * res, 
      int indent,
      options_t * opts
    ) {
  // Print a summary of the structure properties

  char * seq_string = NULL;
  char * struc_string = NULL;
  char * ind = (char *) malloc((indent + 1)* sizeof(char));
  int n_chars;
  DBL_TYPE defect = 0;
  int n_nucs = res->n_nucs;
  int i_nuc = 0;
  int i_break = 0;
  int i_ch;

  for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    defect += res->defects[i_nuc];
  }
  for (i_ch = 0; i_ch < indent; i_ch++) {
    ind[i_ch] = ' ';
  }
  ind[indent] = '\0';

  // int n_nodes = res->n_nodes;
  // int i_node;

  n_chars = res->n_nucs + res->n_breaks + 1;

  seq_string = (char *) malloc(n_chars * sizeof(char));
  struc_string = (char *) malloc(n_chars * sizeof(char));
  get_sequence_plus(seq_string, res, opts);
  get_dpp(struc_string, res);


  fprintf(f, "Sequence : %s\n", seq_string);
  fprintf(f, "Structure: %s\n", struc_string);

  i_break = 0;
  fprintf(f, "IntSeq   : ");
  for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    if (i_break < res->n_breaks && res->breaks[i_break] == i_nuc) {
      fprintf(f, "+");
      i_break ++;
    }
    fprintf(f, "%i", res->sequence[i_nuc]);
  }
  fprintf(f, "\n");

  fprintf(f, "Pfunc : %Lf\n", res->pfunc);
  fprintf(f, "defect : %Lf\n", defect);


  free(seq_string);
  free(struc_string);
  free(ind);
  return ERR_OK;
}
/***********************************************************/

int pairs_to_dpp(char * str, int * struc, int n_nucs) {
  int i_nuc, j_nuc;
  int i_ch;
  j_nuc = 0;
  i_ch = 0;
  for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    if (struc[j_nuc] < 0) {
      str[i_ch] = '.';
    } else if (struc[j_nuc] < j_nuc) {
      str[i_ch] = ')';
    } else {
      str[i_ch] = '(';
    }
    i_ch ++;
    j_nuc ++;
  }
  str[i_ch] = '\0';
  return ERR_OK;
}
/***********************************************************/
int get_spec_dpp(char * str, design_struc_t * des, 
    design_spec_t * spec) {
  int n_nucs;
  int i_nuc, j_nuc;
  int i_ch;
  int i_str;
  int c_str;
  j_nuc = 0;
  i_ch = 0;
  for (i_str = 0; i_str < des->n_strands; i_str++) {
    c_str = des->strands[i_str];
    n_nucs = spec->seqs.strands.specs[c_str].n;
    for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
      if (des->struc[j_nuc] < 0) {
        str[i_ch] = '.';
      } else if (des->struc[j_nuc] < j_nuc) {
        str[i_ch] = ')';
      } else {
        str[i_ch] = '(';
      }
      i_ch ++;
      j_nuc ++;
    }
    str[i_ch] = '+';
    i_ch ++;
  }
  str[i_ch-1] = '\0';
  return ERR_OK;
}
/***********************************************************/

/***********************************************************/
int get_dpp(char * str, result_struc_t * res) {
  int n_nucs = res->n_nucs;
  int i_nuc;
  int j_nuc;
  int i_break = 0;
  int k_nuc = 0;
  for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    if (i_break < res->n_breaks && res->breaks[i_break] == i_nuc ) {
      str[k_nuc] = '+';
      i_break++;
      k_nuc ++;
    }
    j_nuc = res->structure[i_nuc];
    if (j_nuc == -1) {
      str[k_nuc] = '.';
    } else if (j_nuc < i_nuc) {
      str[k_nuc] = ')';
    } else {
      str[k_nuc] = '(';
    }
    k_nuc ++;
  }
  str[k_nuc] = '\0';

  return ERR_OK;
}
/***********************************************************/

/***********************************************************/
int get_sequence_plus(char * str, result_struc_t * res,
    options_t * opts) {
  int n_nucs = res->n_nucs;
  int i_nuc;
  int n_breaks = res->n_breaks;
  int i_break;
  int j_nuc;

  convert_nucs_to_str(str, res->sequence, n_nucs, opts);
  i_break = n_breaks - 1;
  j_nuc = n_nucs + n_breaks - 1;

  for (i_nuc = n_nucs - 1; i_nuc >= 0; i_nuc--) {
    str[j_nuc] = str[i_nuc];
    j_nuc --;
    if (i_break >= 0 && res->breaks[i_break] == i_nuc) {
      str[j_nuc] = '+';
      j_nuc --;
      i_break --;
    }
  }
  str[n_nucs + n_breaks] = '\0';
  return ERR_OK;
}
/***********************************************************/

void print_mutation_weights(FILE * f, result_t * res, 
    mutation_t * mut, int indent) {
  int i;

  char * ind = (char *) malloc((indent + 1)* sizeof(char));
  for (i = 0; i < indent; i++) {
    ind[i] = ' ';
  }
  ind[indent] = '\0';

  fprintf(f, "%sloc: ", ind);
  for (i = 0; i < mut->n; i++) {
    fprintf(f, "%4i ", mut->ids[i]);
  }
  fprintf(f, "\n");
  fprintf(f, "%snuc: ", ind);
  for (i = 0; i < mut->n; i++) {
    fprintf(f, "%4i ", mut->nucs[i]);
  }
  fprintf(f, "\n");
  fprintf(f, "%sdum: ", ind);
  for (i = 0; i < mut->n; i++) {
    fprintf(f, "%4i ", mut->dum[i]);
  }
  fprintf(f, "\n");
  fprintf(f, "%sweight: ", ind);
  int si, ni;
  int csi;
  DBL_TYPE cur_defect, cur_conc, cur_def, cur_des;
  for (i = 0; i < mut->n; i++) {
    cur_defect = 0;
    for (si = 0; si < res->tubes[0].n_strucs; si++) {
      if (res->tubes[0].included[si] && res->tubes[0].target[si]) {
        csi = res->tubes[0].included_ind[si];
        cur_conc = res->tubes[0].x[si];
        cur_des = res->tubes[0].target_x[si];
        cur_def = cur_des - cur_conc;
        if (cur_conc > cur_des) {
          cur_conc = cur_des;
          cur_def = 0;
        }
        for (ni = 0; ni < res->strucs[csi].n_nucs; ni++) {
          if (res->strucs[si].nuc_ids[ni] == mut->ids[i]) {
            cur_defect += 
              cur_conc * res->strucs[si].defects[ni] + 
              cur_def;
          }
        }
      }
    }
    fprintf(f, "%8.3Le ", cur_defect);
  }
  fprintf(f, "\n");
  free(ind);
}

void print_mutation(FILE * f, mutation_t * mut, int indent) {
  int i;

  char * ind = (char *) malloc((indent + 1)* sizeof(char));
  for (i = 0; i < indent; i++) {
    ind[i] = ' ';
  }
  ind[indent] = '\0';

  fprintf(f, "%sloc: ", ind);
  for (i = 0; i < mut->n; i++) {
    fprintf(f, "%4i ", mut->ids[i]);
  }
  fprintf(f, "\n");
  fprintf(f, "%snuc: ", ind);
  for (i = 0; i < mut->n; i++) {
    fprintf(f, "%4i ", mut->nucs[i]);
  }
  fprintf(f, "\n");
  fprintf(f, "%sdum: ", ind);
  for (i = 0; i < mut->n; i++) {
    fprintf(f, "%4i ", mut->dum[i]);
  }
  fprintf(f, "\n");
  free(ind);
}


