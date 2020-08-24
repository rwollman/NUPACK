#include "pathway_design.h"

static void print_assumed_bp(result_tree_t * restree, struc_tree_t * tree) ;
static int pairs_cross(int i_nuc, int j_nuc, int d_nuc, int e_nuc);
static DBL_TYPE get_helix_min_ppair(DBL_TYPE * ppairs, 
    int i_nuc, int j_nuc, int n_nucs, int h_split);
static int get_helix_comp(int i_nuc, int j_nuc, int n_nucs, 
    int h_split, int * nucs, int * comp_map) ;


int init_opts(options_t * opts) {
  struct timeval curtime;

  opts->temperature = ZERO_C_IN_KELVIN + 37.0;
  opts->sodium = NUPACK_DEF_SODIUM;
  opts->magnesium = NUPACK_DEF_MAGNESIUM;
  opts->min_ppair_saved = 0.00001;
  opts->f_split = NUPACK_DEF_F_SPLIT;
  opts->gc_init_prob = NUPACK_DEF_GC_INIT_PROB;
  opts->M_unfavorable = NUPACK_DEF_M_BAD;
  opts->M_reseed = NUPACK_DEF_M_RESEED;
  opts->N_split = NUPACK_DEF_N_SPLIT;
  opts->H_split = NUPACK_DEF_H_SPLIT;
  opts->M_leafopt = NUPACK_DEF_M_REOPT;
  opts->allow_mismatch = 0;
  opts->allow_wobble = 0;
  opts->material = NUPACK_DEF_MATERIAL;
  opts->dangle_type = 1;
  opts->use_long_helix = 0;
  opts->print_leaves = 0;
  opts->print_steps = 0;
  opts->max_print_steps = 10;
  opts->file_prefix = NULL;
  opts->disable_defect_weights = 0;
  opts->disable_focus = 0;
  opts->forbid_splits = 0;
  opts->include_all = 0;
  opts->f_passive = NUPACK_DEF_F_PASSIVE; // Deflate concentrations by 
  opts->fake_dummies = 0;
  opts->include_dummies = 0;
  opts->single_decomp = 0;
                          //  f_passive * f_stop 
  opts->f_stringent = NUPACK_DEF_F_STRINGENT;   // Deflate stop condition by f_relax
  opts->f_refocus =   NUPACK_DEF_F_REFOCUS;  // Only include up to f_include of badness on refocusing
  opts->f_redecomp =  NUPACK_DEF_F_REDECOMP; // Only include up to f_redecomp of badness on redecomposition
  opts->bonus_per_split = NUPACK_DEF_BONUS_PER_SPLIT;
  opts->seed = 0;

  opts->designing = 1;

  opts->allowed_opt_time = 31536000; 
  // Approximate seconds in a year. Doubtful anyone will exceed this.

  gettimeofday(&curtime, NULL);
  opts->start_time = curtime.tv_sec + (1e-6 * curtime.tv_usec);
  return ERR_OK;
}

int get_n_mutable_nucleotides(design_spec_t * spec) {
  int n_strucs = spec->n_strucs;
  int n_strands, n_nucs;
  int i_struc, i_strand, c_strand, i_nuc;
  int n_free = 0;
  int cur_nuc_id;
  int cur_nuc;

  for (i_struc = 0; i_struc < n_strucs; i_struc++) {
    n_strands = spec->strucs[i_struc].n_strands;
    for (i_strand = 0; i_strand < n_strands; i_strand++) {
      c_strand = spec->strucs[i_struc].strands[i_strand];
      n_nucs = spec->seqs.strands.specs[c_strand].n;
      for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
        cur_nuc_id = spec->seqs.strands.specs[c_strand].nucs[i_nuc];
        cur_nuc = spec->seqs.nucs.nucs[cur_nuc_id];

        if (!(cur_nuc == BASE_A
              || cur_nuc == BASE_C
              || cur_nuc == BASE_G
              || cur_nuc == BASE_U
              || cur_nuc == BASE_T)) {
          n_free ++;
        }
      }
    }
  }

  return n_free;
}

int init_sequence_spec(
      sequence_spec_t * seqs
    ) {
  
  seqs->nucs.nucs = NULL;
  seqs->nucs.n = 0;

  seqs->comp_map.nucs = NULL;
  seqs->comp_map.n = 0;

  seqs->domains.specs = NULL;
  seqs->domains.names = NULL;
  seqs->domains.n = 0;
  seqs->domains.cap = 0;

  seqs->strands.specs = NULL;
  seqs->strands.names = NULL;
  seqs->strands.n = 0;
  seqs->strands.cap = 0;
  return ERR_OK;
}

int init_design_spec (
      design_spec_t * spec
    ) {

  int cap_strucs = 10;
  int cap_tubes = 1;
  init_sequence_spec(&(spec->seqs));
  init_opts(&(spec->opts));

  spec->strucs = (design_struc_t *) malloc(sizeof(design_struc_t) * 
      cap_strucs);
  spec->tubes = (design_tube_t *) malloc(sizeof(design_tube_t) *
      cap_tubes);

  spec->struc_names = (char **) malloc(sizeof(char*) * cap_strucs);
  spec->tube_names = (char **) malloc(sizeof(char*) * cap_tubes);

  spec->n_strucs = 0;
  spec->n_tubes = 0;

  spec->cap_strucs = cap_strucs;
  spec->cap_tubes = cap_tubes;

  spec->n_orderings = 0;
  spec->cap_orderings = cap_strucs;
  spec->orderings = (int **) malloc(sizeof(int*) * cap_strucs);
  spec->symmetry = (int*) malloc(sizeof(int) * cap_strucs);
  spec->struc_map = (int*) malloc(sizeof(int) * cap_strucs);
  spec->n_strands = (int*) malloc(sizeof(int) * cap_strucs);

  return ERR_OK;
}

int init_sequence(sequence_t * seq, int n_nucs) {
  int i_nuc;
  seq->n = n_nucs;
  seq->nucs = (int*) malloc(n_nucs * sizeof(int));
  check_mem(seq->nucs);
  for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    seq->nucs[i_nuc] = 0;
  }
  return ERR_OK;
error:
  return ERR_OOM;
}

int resize_sequence(sequence_t * seq, int n_nucs) {
  int i_nuc;
  int * tempint;
  int oldn = seq->n;
  seq->n = n_nucs;
  tempint = (int*) realloc(seq->nucs, n_nucs * sizeof(int));
  check_mem(tempint);
  seq->nucs = tempint;
  for (i_nuc = oldn; i_nuc < n_nucs; i_nuc++) {
    seq->nucs[i_nuc] = 0;
  }
  return ERR_OK;
error:
  return ERR_OOM;
}

int init_seqspec_list(seqspec_list_t * seqs) {
  seqs->specs = NULL;
  seqs->names = NULL;
  seqs->n = 0;
  seqs->cap = 0;
  return ERR_OK;
}

int free_seqspec_list(seqspec_list_t * seqs) {
  int i;
  for (i = 0; i < seqs->n; i++) {
    free(seqs->names[i]);
    free_sequence(seqs->specs + i);
  }
  free(seqs->names);
  free(seqs->specs);
  seqs->names = NULL;
  seqs->specs = NULL;
  seqs->n = 0;
  seqs->cap = 0;
  return ERR_OK;
}

int append_sequence_seqspec(
    seqspec_list_t * seqs, 
    sequence_t * seq,
    const char * name) {
  int oldn = seqs->n;
  int newcap;
  sequence_t * temp_seq;
  char ** temp_chars;
  int i;
  int namelen = strlen(name);
  if (oldn >= seqs->cap) {
    newcap = (oldn + 1) * 2;
    temp_seq = (sequence_t*) realloc(seqs->specs, newcap * sizeof(sequence_t));
    check_mem(temp_seq);
    seqs->specs = temp_seq;

    temp_chars = (char**) realloc(seqs->names, newcap * sizeof(char*));
    check_mem(temp_chars);
    seqs->names = temp_chars;
    for (i = oldn; i < newcap; i++) {
      seqs->names[i] = NULL;
      seqs->specs[i].nucs = NULL;
      seqs->specs[i].n = 0;
    }
    seqs->cap = newcap;
  }

  check(ERR_OK == init_sequence(seqs->specs + oldn, seq->n),
      "Error initializing sequence");
  check(ERR_OK == copy_sequence(seqs->specs + oldn, seq),
      "Error copying sequence");

  seqs->names[oldn] = (char *) malloc((namelen + 1) * sizeof(char));
  check_mem(seqs->names[oldn]);

  strncpy(seqs->names[oldn], name, namelen + 1);
  seqs->names[oldn][namelen] = '\0';

  seqs->n += 1;
  return ERR_OK;
error:
  return ERR_OOM;
}

int delete_sequence_seqspec(seqspec_list_t * seqs, int i_str) {
  int oldn = seqs->n;
  int i;
  free_sequence(seqs->specs + i_str);
  free(seqs->names[i_str]);
  for (i = i_str; i < oldn - 1; i++) {
    seqs->specs[i] = seqs->specs[i+1];
    seqs->names[i] = seqs->names[i+1];
  }

  seqs->names[oldn-1] = NULL;
  seqs->specs[i].nucs = NULL;
  seqs->specs[i].n = 0;
  seqs->n --;

  return ERR_OK;
}

int init_sequence_list(sequence_list_t * seqs, seqspec_list_t * spec) {
  int i_str, n_strs = spec->n;
  int n_nucs ;
  seqs->n = spec->n;
  seqs->seqs = (sequence_t *) malloc(n_strs * sizeof(sequence_t));
  for (i_str = 0; i_str < n_strs; i_str++) {
    n_nucs = spec->specs[i_str].n;
    init_sequence(seqs->seqs + i_str, n_nucs);
  }
  return ERR_OK;
}

static int init_copy_seqstate(seqstate_t * dest, seqstate_t * src) {
  check(ERR_OK == init_copy_sequence(&(dest->nucspec), &(src->nucspec)),
      "Error initializing nucspec");
  check(ERR_OK == init_copy_sequence(&(dest->dumspec), &(src->dumspec)),
      "Error initializing dumspec");
  check(ERR_OK == init_copy_sequence_list(&(dest->domains), &(src->domains)),
      "Error initializing domains");
  check(ERR_OK == init_copy_sequence_list(&(dest->strands), &(src->strands)),
      "Error initializing strands");
  return ERR_OK;
error:
  return ERR_OOM;
}

int init_seqstate(seqstate_t * seqs, design_spec_t * spec) {
  check(ERR_OK == init_sequence(&(seqs->nucspec), spec->seqs.nucs.n),
      "Error initializing nucspec");
  check(ERR_OK == init_sequence(&(seqs->dumspec), spec->seqs.nucs.n),
      "Error initializing dumspec");
  check(ERR_OK == init_sequence_list(&(seqs->domains), &(spec->seqs.domains)),
      "Error initializing domains");
  check(ERR_OK == init_sequence_list(&(seqs->strands), &(spec->seqs.strands)),
      "Error initializing strands");
  return ERR_OK;
error:
  return ERR_OOM;
}

int init_copy_sequence(sequence_t * dest, sequence_t * src) {
  init_sequence(dest, src->n);
  copy_sequence(dest, src);
  return ERR_OK;
}

int copy_sequence(sequence_t * dest, sequence_t * src) {
  check(dest->n == src->n, "sequence lengths are unequal dest: %i src: %i", dest->n, src->n);
  int i;
  for (i = 0; i < dest->n; i++) {
    dest->nucs[i] = src->nucs[i];
  }
  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

int init_copy_sequence_list(sequence_list_t * dest, sequence_list_t * src) {
  int i_str, n_strs = src->n;
  dest->n = src->n;
  dest->seqs = (sequence_t *) malloc(n_strs * sizeof(sequence_t));
  check_mem(dest->seqs);
  for (i_str = 0; i_str < n_strs; i_str++) {
    init_copy_sequence(dest->seqs + i_str, src->seqs + i_str);
  }
  return ERR_OK;
error:
  return ERR_OOM;
}

int copy_sequence_list(sequence_list_t * dest, sequence_list_t * src) {
  check(dest->n == src->n, "number of sequences are unequal dest: %i src: %i", dest->n, src->n);
  int i;
  for (i = 0; i < dest->n; i++) {
    check(ERR_OK == copy_sequence(dest->seqs + i, src->seqs + i),
        "Error copying sequence %i", i);
  }
  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

int copy_dummy(seqstate_t * st) {
  check(ERR_OK == copy_sequence(&(st->dumspec), &(st->nucspec)),
      "Error copying sequence");
  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

int copy_seqstate(seqstate_t * dest, seqstate_t * src) {
  check(ERR_OK == copy_sequence(&(dest->nucspec), &(src->nucspec)),
      "Error copying nucspec");
  check(ERR_OK == copy_sequence(&(dest->dumspec), &(src->dumspec)),
      "Error copying dumspec");
  check(ERR_OK == copy_sequence_list(&(dest->domains), &(src->domains)),
      "Error copying domains");
  check(ERR_OK == copy_sequence_list(&(dest->strands), &(src->strands)),
      "Error copying strands");
  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}


int refresh_seqstate(
      seqstate_t * seqs,
      design_spec_t * spec
    ) {

  int i_nuc, c_nuc, c_nt, n_nucs, i, n;
  n = spec->seqs.strands.n;
  for (i = 0; i < n; i++) {
    n_nucs = spec->seqs.strands.specs[i].n;
    for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
      c_nuc = spec->seqs.strands.specs[i].nucs[i_nuc];
      c_nt = seqs->nucspec.nucs[c_nuc];
      seqs->strands.seqs[i].nucs[i_nuc] = c_nt;
    }
  }
  n = spec->seqs.domains.n;
  for (i = 0; i < n; i++) {
    n_nucs = spec->seqs.domains.specs[i].n;
    for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
      c_nuc = spec->seqs.domains.specs[i].nucs[i_nuc];
      c_nt = seqs->nucspec.nucs[c_nuc];
      seqs->domains.seqs[i].nucs[i_nuc] = c_nt;
    }
  }
  return ERR_OK;
}

int apply_mutation(seqstate_t * seqs, mutation_t * mut) {
  int i_nuc;
  int c_id;
  int c_nuc;
  int c_dum;

  for (i_nuc = 0; i_nuc < mut->n; i_nuc++) {
    c_id = mut->ids[i_nuc];
    c_nuc = mut->nucs[i_nuc];
    c_dum = mut->dum[i_nuc];
    if (c_dum) {
      seqs->dumspec.nucs[c_id] = c_nuc;
    } else {
      seqs->nucspec.nucs[c_id] = c_nuc;
    }
  }
  return ERR_OK;
}

void free_sequence(sequence_t * seq) {
  free(seq->nucs);
  seq->nucs = NULL;
  seq->n = 0;
}

void free_sequence_list(sequence_list_t * seqs) {
  int i;
  for (i = 0; i < seqs->n; i++) {
    free_sequence(seqs->seqs + i);
  }
  free(seqs->seqs);
  seqs->n = 0;
  seqs->seqs = NULL;
}

void free_seqstate(seqstate_t * seqs) {
  free_sequence(&(seqs->nucspec));
  free_sequence(&(seqs->dumspec));
  free_sequence_list(&(seqs->domains));
  free_sequence_list(&(seqs->strands));
}


int init_mutation(mutation_t * mut) {
  mut->n = 0;
  mut->cap = 0;
  mut->ids = NULL;
  mut->nucs = NULL;
  mut->dum = NULL;
  return ERR_OK;
}

int resize_mutation(mutation_t * mut, int size) {
  int * temp_int;
  temp_int = (int*) realloc(mut->ids, size * sizeof(int));
  check_mem(temp_int);
  mut->ids = temp_int;

  temp_int = (int*) realloc(mut->nucs, size * sizeof(int));
  check_mem(temp_int);
  mut->nucs = temp_int;

  temp_int = (int*) realloc(mut->dum, size * sizeof(int));
  check_mem(temp_int);
  mut->dum = temp_int;

  mut->cap = 2 * (mut->n + 1);
  return ERR_OK;
error:
  return ERR_OOM;
}

int copy_mutation(mutation_t * dest, mutation_t * src) {
  int i;
  if (dest->cap < src->n) {
    check(ERR_OK == resize_mutation(dest, src->n),
        "Error resizing mutation");
  }
  
  for (i = 0; i < src->n; i++) {
    dest->ids[i] = src->ids[i];
    dest->nucs[i] = src->nucs[i];
    dest->dum[i] = src->dum[i];
  }
  dest->n = src->n;
  return ERR_OK;
error:
  return ERR_OOM;
}

int clear_mutation(mutation_t * mut) {
  mut->n = 0;
  return ERR_OK;
}

int add_mutation_nuc(mutation_t * mut, int id, int nuc, int dummy) {
  int i, cid, cnuc, cdummy, tid, tnuc, tdummy;
  if (mut->n >= mut->cap) {
    check(ERR_OK == resize_mutation(mut, 2 * (mut->n + 1)),
          "Error resizing mutation");
  }
  cid = id;
  cnuc = nuc;
  cdummy = dummy;
  for (i = 0; i < mut->n; i++) {
    if (cid < mut->ids[i] || (cid == mut->ids[i] && cnuc < mut->nucs[i])
        || (cid == mut->ids[i] && cnuc == mut->nucs[i] 
          && cdummy < mut->dum[i])) {
      tid = mut->ids[i];
      tnuc = mut->nucs[i];
      tdummy = mut->dum[i];
      mut->ids[i] = cid;
      mut->nucs[i] = cnuc;
      mut->dum[i] = cdummy;
      cid = tid;
      cnuc = tnuc;
      cdummy = tdummy;
    }
  }
  mut->ids[mut->n] = cid;
  mut->nucs[mut->n] = cnuc;
  mut->dum[mut->n] = cdummy;
  mut->n += 1;

  return ERR_OK;
error:
  return ERR_OOM;
}

int get_comp_code(int nuc_code, int allow_wobble) {
  int codes[] = {
    BASE_N,
    BASE_A,
    BASE_C,
    BASE_G,
    BASE_U,
    BASE_R,
    BASE_M,
    BASE_S,
    BASE_W,
    BASE_K,
    BASE_Y,
    BASE_V,
    BASE_H,
    BASE_D,
    BASE_B
  };
  int comps[] = {
    BASE_N, // N -> N
    BASE_U, // A -> U
    BASE_G, // C -> G
    BASE_Y, // G -> Y
    BASE_R, // U -> R
    BASE_Y, // R -> Y
    BASE_K, // M -> K
    BASE_S, // S -> S
    BASE_W, // W -> W
    BASE_N, // K -> N
    BASE_R, // Y -> R
    BASE_B, // V -> B
    BASE_D, // H -> D
    BASE_N, // D -> N
    BASE_N  // B -> N
  };

  int comps_no_wobble[] = {
    BASE_N, // N -> N
    BASE_U, // A -> U
    BASE_G, // C -> G
    BASE_C, // G -> C
    BASE_A, // U -> A
    BASE_Y, // R -> Y
    BASE_K, // M -> K
    BASE_S, // S -> S
    BASE_W, // W -> W
    BASE_M, // K -> M
    BASE_R, // Y -> R
    BASE_B, // V -> B
    BASE_D, // H -> D
    BASE_H, // D -> H
    BASE_V  // B -> V
  };

  int n = 14;
  int i;
  for (i = 0; i < n; i++) {
    if (nuc_code == codes[i]) {
      break;
    }
  }
  if (allow_wobble) {
    return comps[i];
  } else {
    return comps_no_wobble[i];
  }
}

int get_nuc_code(char nuc_char) {
  char * chars = "NACGUTMRWSYKVHDB";
  int codes[] = {
    BASE_N, BASE_A, BASE_C, BASE_G, BASE_U, BASE_U, 
    BASE_M, BASE_R, BASE_W, BASE_S, BASE_Y, BASE_K, 
    BASE_V, BASE_H, BASE_D, BASE_B
  };
  int i;
  int n_codes = 16;
  for (i = 0; i < n_codes; i++) {
    if (nuc_char == chars[i]) {
      return codes[i];
    }
  }
  return codes[0];
}

int sample_mutation(
      mutation_t * mut, 
      int nuc_id, 
      int dummy, 
      seqstate_t * seqs, 
      design_spec_t * spec
    ) {
  int f_con = spec->seqs.nucs.nucs[nuc_id];
  int c_id = spec->seqs.comp_map.nucs[nuc_id];

  int f_cur = seqs->nucspec.nucs[nuc_id];
  int c_con;
  int c_cur;

  int f_nuc = 0;
  int c_nuc = 0;

  int n_choices;
  int can_wobble = spec->opts.allow_wobble;
  int can_mismatch = spec->opts.allow_mismatch;

  int nuc, i_con;
  DBL_TYPE tot_weight, cur_weight, r_num;


  int n_allowed[15] = {
    4, 
    1, 1, 1, 1, 
    2, 2, 2, 2, 2, 2,
    3, 3, 3, 3
  };
  int allowed_seqs[15][4] = {
    {1, 2, 3, 4}, // N
    {1,-1,-1,-1}, // A
    {2,-1,-1,-1}, // C
    {3,-1,-1,-1}, // G
    {4,-1,-1,-1}, // U / T
    {1, 3,-1,-1}, // AG
    {1, 2,-1,-1}, // AC
    {2, 3,-1,-1}, // CG
    {1, 4,-1,-1}, // AU
    {3, 4,-1,-1}, // GU
    {2, 4,-1,-1}, // CU
    {1, 2, 3,-1}, // ACG
    {1, 2, 4,-1}, // ACU
    {1, 3, 4,-1}, // AGU
    {2, 3, 4,-1}, // CGU
  };
  int nuc_comps[5] = {
    0, 4, 3, 2, 1
  };
  int nuc_wobble[5] = {
    0, -1, -1, 4, 3
  };


  double def_weight[] = {
    0.0, 1.0, 1.0, 1.0, 1.0
  };
  double weight[] = {
    0.0, 1.0, 1.0, 1.0, 1.0
  };
  if (c_id >= 0) {
    c_con = spec->seqs.nucs.nucs[c_id];
    c_cur = seqs->nucspec.nucs[c_id];
  } else {
    c_con = 0;
    c_cur = 1;
  }
  clear_mutation(mut);

  check(f_con >= 0 && f_con <= 15, 
      "Invalid constraint at %d: %d", nuc_id, f_con);
  check(c_con >= 0 && c_con <= 15,
      "Invalid comp constraint at %d.%d: %d", nuc_id, c_id, c_con);
  n_choices = n_allowed[f_con];

  for (i_con = 1; i_con < 5; i_con++) {
    weight[i_con] = 0;
  }
  for (i_con = 0; i_con < n_choices; i_con++) {
    weight[allowed_seqs[f_con][i_con]] = def_weight[allowed_seqs[f_con][i_con]];
  }
  weight[f_cur] = 0;

  tot_weight = 0;
  for (nuc = 1; nuc < 5; nuc++) {
    tot_weight += weight[nuc];
  }
  if (tot_weight > 0) {
    nuc = 5;
    while (nuc == 5) {
      r_num = ((DBL_TYPE) genrand_real1()) * tot_weight;
      cur_weight = 0;
      for (nuc = 1; nuc < 4; nuc++) {
        cur_weight += weight[nuc];
        if (cur_weight > r_num) {
          break;
        }
      }
    }
    check(nuc >= 1 && nuc <= 4, "Invalid nucleotide code: %i (%Lf / %Lf)", nuc,
        cur_weight, tot_weight);

    f_nuc = nuc;

    for (i_con = 1; i_con < 5; i_con++) {
      weight[i_con] = 0;
    }

    weight[nuc_comps[f_nuc]] = def_weight[nuc_comps[f_nuc]];

    if (can_wobble) {
      weight[nuc_wobble[f_nuc]] = def_weight[nuc_wobble[f_nuc]];
    }

    if (can_mismatch) {
      weight[c_cur] = def_weight[c_cur];
    }

    tot_weight = 0;
    for (nuc = 1; nuc < 5; nuc++) {
      tot_weight += weight[nuc];
    }

    if (tot_weight > 0) {
      nuc = 5;
      while (nuc == 5) {
        r_num = ((DBL_TYPE) genrand_real1()) * tot_weight;
        cur_weight = 0;
        for (nuc = 1; nuc < 5; nuc++) {
          cur_weight += weight[nuc];
          if (cur_weight > r_num) {
            break;
          }
        }
      }

      check(nuc >= 1 && nuc <= 4, "Invalid nucleotide code: %i", nuc);

      c_nuc = nuc;

      add_mutation_nuc(mut, nuc_id, f_nuc, dummy);
      if (c_id >= 0) {
        add_mutation_nuc(mut, c_id, c_nuc, dummy);
      }
    }
  }
  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

void free_mutation(mutation_t * mut) {
  mut->n = 0;
  mut->cap = 0;
  free(mut->ids);
  free(mut->nucs);
  free(mut->dum);
  mut->dum = NULL;
  mut->ids = NULL;
  mut->nucs = NULL;
}

int init_mutation_list(mutation_list_t * muts) {
  muts->n = 0;
  muts->cap = 0;
  muts->muts = NULL;
  return ERR_OK;
}

int clear_mutation_list(mutation_list_t * muts) {
  int i;
  for (i = 0; i < muts->n; i++) {
    free_mutation(muts->muts + i);
  }
  muts->n = 0;
  return ERR_OK;
}

int find_mutation(mutation_list_t * muts, mutation_t * mut) {
  int i, j, equal, ind;
  ind = -1;
  for (i = 0; i < muts->n; i++) {
    equal = 0;
    if (muts->muts[i].n == mut->n) {
      equal = 1;
      for (j = 0; j < muts->muts[i].n; j++) {
        if (muts->muts[i].ids[j] != mut->ids[j]
            || muts->muts[i].nucs[j] != mut->nucs[j]
            || muts->muts[i].dum[j] != mut->dum[j]) {
          equal = 0;
          break;
        }
      }
    } 
    if (equal) {
      ind = i;
      break;
    }
  }
  return ind;
}

int resize_mutation_list(mutation_list_t * muts, int size) {
  int i;
  mutation_t * tmpmuts;
  if (size < muts->n) {
    for (i = size; i < muts->n; i++) {
      free_mutation(muts->muts + i);
    }
    muts->n = size;
  }
  tmpmuts = (mutation_t *) realloc(muts->muts, size * sizeof(mutation_t));
  check_mem(tmpmuts);
  muts->muts = tmpmuts;
  muts->cap = size;
  return ERR_OK;
error:
  return ERR_OOM;
}

int append_mutation(mutation_list_t * muts, mutation_t * mut) {
  if (muts->n >= muts->cap) {
    check(ERR_OK == resize_mutation_list(muts, (muts->n + 1) * 2),
        "Error resizing mutation list");
  }
  init_mutation(muts->muts + muts->n);
  copy_mutation(muts->muts + muts->n, mut);
  muts->n += 1;
  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

void free_mutation_list(mutation_list_t * muts) {
  int i;
  for (i = 0; i < muts->n; i++) {
    free_mutation(muts->muts + i);
  }
  free(muts->muts);
  muts->n = 0;
  muts->cap = 0;
  muts->muts = NULL;
}

int init_result_tube(
      result_tube_t * r_tube
    ) {
  r_tube->n_strucs = 0;
  r_tube->cap_strucs = 0;
  r_tube->x = NULL;
  r_tube->target_x = NULL;
  r_tube->target = NULL;
  r_tube->included = NULL;
  r_tube->included_ind = NULL;
  r_tube->generated_ind = NULL;
  r_tube->defect = DBL_MAX;

  return ERR_OK;
}

static int copy_result_tube (
      result_tube_t * dest,
      result_tube_t * src
    ) {
  int cap_strucs = src->cap_strucs;
  int n_strucs = src->n_strucs;
  int i_struc;

  dest->n_strucs = n_strucs;
  dest->cap_strucs = cap_strucs;
  dest->x = (DBL_TYPE*) realloc(dest->x, 
      cap_strucs * sizeof(DBL_TYPE));
  dest->target_x = (DBL_TYPE*) realloc(dest->target_x,
      cap_strucs * sizeof(DBL_TYPE));
  dest->target = (int*) realloc(dest->target, 
      cap_strucs * sizeof(int));
  dest->included = (int*) realloc(dest->included, 
      cap_strucs * sizeof(int));
  dest->included_ind = (int*) realloc(dest->included_ind, 
      cap_strucs * sizeof(int));
  dest->generated_ind = (int*) realloc(dest->generated_ind, 
      cap_strucs * sizeof(int));
  dest->defect = src->defect;

  for (i_struc = 0; i_struc < n_strucs; i_struc++) {
    dest->x[i_struc] = src->x[i_struc];
    dest->target_x[i_struc] = src->target_x[i_struc];
    dest->target[i_struc] = src->target[i_struc];
    dest->included[i_struc] = src->included[i_struc];
    dest->included_ind[i_struc] = src->included_ind[i_struc];
    dest->generated_ind[i_struc] = src->generated_ind[i_struc];
  }
  return ERR_OK;
}

int init_result_tube_spec(
    result_tube_t * r_tube,  
    design_tube_t * spec_tube ) {
  int n_strucs = spec_tube->n_strucs;
  int cap_strucs = spec_tube->n_strucs;
  int i_struc;

  r_tube->n_strucs = n_strucs;
  r_tube->cap_strucs = cap_strucs;
  r_tube->x = (DBL_TYPE*) malloc(cap_strucs * sizeof(DBL_TYPE));
  r_tube->target_x = (DBL_TYPE*) malloc(cap_strucs * sizeof(DBL_TYPE));
  r_tube->target = (int*) malloc(cap_strucs * sizeof(int));
  r_tube->included = (int*) malloc(cap_strucs * sizeof(int));
  r_tube->included_ind = (int*) malloc(cap_strucs * sizeof(int));
  r_tube->generated_ind = (int*) malloc(cap_strucs * sizeof(int));
  r_tube->defect = DBL_MAX;
  
  for (i_struc = 0; i_struc < n_strucs; i_struc++) {
    r_tube->x[i_struc] = spec_tube->target_x[i_struc];
    r_tube->target_x[i_struc] = spec_tube->target_x[i_struc];
    r_tube->target[i_struc] = spec_tube->target[i_struc];
    r_tube->included[i_struc] = spec_tube->included[i_struc];
    r_tube->included_ind[i_struc] = spec_tube->included_ind[i_struc];
    r_tube->generated_ind[i_struc] = spec_tube->generated_ind[i_struc];
  }

  return ERR_OK;
}

int init_result(
      result_t * res,
      design_spec_t * spec) {

  int n_strucs = spec->n_strucs;
  int n_ords = spec->n_orderings;

  int i_struc;
  int i_ord;

  int n_strs;
  int i_str;

  res->strucs = (result_struc_t *) 
    malloc(sizeof(result_struc_t) * n_strucs);
  res->tubes = (result_tube_t *) 
    malloc(sizeof(result_tube_t));
  res->orderings = (int**) malloc(sizeof(int*) * n_ords);
  res->n_strands = (int*) malloc(sizeof(int) * n_ords);
  res->struc_map = (int*) malloc(sizeof(int) * n_ords);
  res->order_len = (int*) malloc(sizeof(int) * n_ords);
  res->included = (int*) malloc(sizeof(int) * n_ords);
  res->dG = (DBL_TYPE*) malloc(sizeof(DBL_TYPE) * n_ords);
  res->eval_time = (DBL_TYPE*) malloc(sizeof(DBL_TYPE) * n_ords);
  res->offtarget_seqs = NULL;

  for (i_struc = 0; i_struc < n_strucs; i_struc++) {
    init_result_struc(res->strucs + i_struc);
  }

  init_result_tube_spec(res->tubes, spec->tubes);

  for (i_ord = 0; i_ord < n_ords; i_ord++) {
    n_strs = spec->n_strands[i_ord];
    res->orderings[i_ord] = (int*) malloc(sizeof(int) * n_strs);
    res->n_strands[i_ord] = n_strs;
    for (i_str = 0; i_str < n_strs; i_str++) {
      res->orderings[i_ord][i_str] = spec->orderings[i_ord][i_str];
    }
    res->struc_map[i_ord] = spec->struc_map[i_ord];
    res->order_len[i_ord] = 0;
    res->included[i_ord] = spec->struc_map[i_ord] == -1 ? 0 : 1;
    res->dG[i_ord] = 100;
    res->eval_time[i_ord] = 0;
  }

  res->n_strucs = n_strucs;
  res->cap_strucs = n_strucs;

  res->n_tubes = 1;
  res->cap_tubes = 1;

  res->n_orderings = n_ords;
  res->cap_orderings = n_ords;
  res->tot_n_strands = spec->seqs.strands.n;

  res->total_defect = DBL_MAX;
  res->id = 0;

  res->elapsed_time = 0;
  res->root_time = 0;

  return ERR_OK;
}

int fill_out_blank_seqs(
      result_t * res,
      design_state_t * state, 
      design_spec_t * spec
    ) {

  int n_strucs = spec->n_strucs;
  int i_struc;
  seqstate_t tseq;
  init_seqstate(&tseq, spec);

  for (i_struc = 0; i_struc < n_strucs; i_struc++) {
    check(ERR_OK == ensure_capacity(res->strucs + i_struc, 
        state->states + i_struc), 
        "Error ensuring capacity");
    check(ERR_OK == update_structure(res->strucs + i_struc, 
          state->states + i_struc, spec), "Error updating structure");
  }

  free_seqstate(&tseq);
  return ERR_OK;
error:
  free_seqstate(&tseq);
  return ERR_OK;
}

int copy_result(
    result_t * dest,
    result_t * src) {
  int n_strucs = src->n_strucs;
  int cap_strucs = src->cap_strucs;
  int i_struc;
  int ** temp_ords;
  int * temp_strs;
  int * temp_strucs;
  DBL_TYPE * temp_dG;
  DBL_TYPE * temp_time;

  int n_tubes = 1;
  int cap_tubes = src->cap_tubes;
  int i_tube;

  int n_ords;
  int i_ord;
  int n_strs;
  int i_str;

  if (cap_strucs > dest->cap_strucs) {
    dest->strucs = (result_struc_t *) realloc(dest->strucs,
        sizeof(result_struc_t) * cap_strucs);
    for (i_struc = dest->n_strucs; i_struc < cap_strucs; i_struc++) {
      init_result_struc(dest->strucs + i_struc);
    }
    dest->cap_strucs = cap_strucs;
  } else if (src->n_strucs < dest->n_strucs) {
    for (i_struc = src->n_strucs; i_struc < dest->n_strucs; i_struc++) {
      free_result_struc(dest->strucs + i_struc);
    }
  }


  if (cap_tubes > dest->cap_tubes) {
    dest->tubes = (result_tube_t *) realloc(dest->tubes,
        sizeof(result_tube_t) * cap_tubes);
    init_result_tube(dest->tubes);
    dest->cap_tubes = cap_tubes;
  }

  for (i_struc = 0; i_struc < n_strucs; i_struc++) {
    copy_result_struc(dest->strucs + i_struc, src->strucs + i_struc);
  }
  dest->n_strucs = n_strucs;

  for (i_tube = 0; i_tube < n_tubes; i_tube++) {
    copy_result_tube(dest->tubes + i_tube, src->tubes + i_tube);
  }
  dest->n_tubes = n_tubes;

  n_ords = src->n_orderings;
  if (dest->cap_orderings < n_ords) {
    temp_dG = (DBL_TYPE*) realloc(dest->dG, sizeof(DBL_TYPE) * n_ords);
    check_mem(temp_dG);
    dest->dG = temp_dG;

    temp_time = (DBL_TYPE*) realloc(dest->eval_time, sizeof(DBL_TYPE) * n_ords);
    check_mem(temp_time);
    dest->eval_time = temp_time;

    temp_ords = (int**) realloc(dest->orderings,
        sizeof(int*) * n_ords);
    check_mem(temp_ords);
    dest->orderings = temp_ords;

    temp_strs = (int*) realloc(dest->n_strands, sizeof(int) * n_ords);
    check_mem(temp_strs);
    dest->n_strands = temp_strs;

    temp_strucs = (int*) realloc(dest->struc_map, sizeof(int) * n_ords);
    check_mem(temp_strucs);
    dest->struc_map = temp_strucs;

    temp_strucs = (int*) realloc(dest->included, sizeof(int) * n_ords);
    check_mem(temp_strucs);
    dest->included = temp_strucs;

    temp_strucs = (int*) realloc(dest->order_len, sizeof(int) * n_ords);
    check_mem(temp_strucs);
    dest->order_len = temp_strucs;

    dest->cap_orderings = n_ords;
  }

  if (dest->n_orderings > n_ords) {
    for (i_ord = n_ords; i_ord < dest->n_orderings; i_ord++) {
      free(dest->orderings[i_ord]);
      dest->orderings[i_ord] = NULL;
    }
  }

  for (i_ord = dest->n_orderings; i_ord < dest->cap_orderings; i_ord++) {
    dest->dG[i_ord] = 0;
    dest->eval_time[i_ord] = 0;
    dest->orderings[i_ord] = NULL;
    dest->n_strands[i_ord] = 0;
    dest->struc_map[i_ord] = 0;
    dest->included[i_ord] = 0;
    dest->order_len[i_ord] = 0;
  }

  for (i_ord = 0; i_ord < n_ords; i_ord++) {
    n_strs = src->n_strands[i_ord];
    temp_strs = (int*) realloc(dest->orderings[i_ord], 
        sizeof(int) * n_strs);
    check_mem(temp_strs);
    dest->orderings[i_ord] = temp_strs;
    dest->n_strands[i_ord] = n_strs;
    for (i_str = 0; i_str < n_strs; i_str++) {
      dest->orderings[i_ord][i_str] = src->orderings[i_ord][i_str];
    }
    dest->struc_map[i_ord] = src->struc_map[i_ord];
    dest->order_len[i_ord] = src->order_len[i_ord];
    dest->included[i_ord] = src->included[i_ord];
  }

  for (i_ord = 0; i_ord < n_ords; i_ord++) {
    dest->dG[i_ord] = src->dG[i_ord];
    dest->eval_time[i_ord] = src->eval_time[i_ord];
  }

  dest->n_orderings = n_ords;

  dest->total_defect = src->total_defect;
  dest->elapsed_time = src->elapsed_time;
  dest->root_time = src->root_time;
  dest->id = src->id;
  dest->n_orderings = src->n_orderings;

  if (src->offtarget_seqs) {
    if (!dest->offtarget_seqs) {
      dest->offtarget_seqs = (seqstate_t *) malloc(sizeof(seqstate_t));
      check_mem(dest->offtarget_seqs);
      init_copy_seqstate(dest->offtarget_seqs, src->offtarget_seqs);
    } else {
      copy_seqstate(dest->offtarget_seqs, src->offtarget_seqs);
    }
  }

  return ERR_OK;
error:
  return ERR_OOM;
}

int free_result(
    result_t * res
    ) {
  int n_strucs = res->n_strucs;
  int n_tubes = res->n_tubes;
  int n_ords = res->n_orderings;

  int i_struc;
  int i_tube;
  int i_ord;

  for (i_struc = 0; i_struc < n_strucs; i_struc++) {
    free_result_struc(res->strucs + i_struc);
  }
  for (i_tube = 0; i_tube < n_tubes; i_tube++) {
    free_result_tube(res->tubes + i_tube);
  }
  for (i_ord = 0; i_ord < n_ords; i_ord++) {
    free(res->orderings[i_ord]);
  }

  if (res->offtarget_seqs) {
    free_seqstate(res->offtarget_seqs);
    free(res->offtarget_seqs);
  }

  free(res->strucs);
  free(res->tubes);
  free(res->orderings);
  free(res->n_strands);
  free(res->struc_map);
  free(res->order_len);
  free(res->included);
  free(res->dG);
  free(res->eval_time);

  res->tubes = NULL;
  res->strucs = NULL;
  res->n_strucs = 0;
  res->cap_strucs = 0;
  res->n_tubes = 0;
  res->cap_tubes = 0;
  res->total_defect = 0;
  res->id = 0;

  res->orderings = NULL;
  res->n_strands = NULL;
  res->struc_map = NULL;
  res->order_len = NULL;
  res->included = NULL;
  res->dG = NULL;
  res->eval_time = NULL;
  res->n_orderings = 0;
  res->cap_orderings = 0;
  return ERR_OK;
}

int free_design_spec (
      design_spec_t * spec
    ) {
  int i;
  free_sequence_spec(&(spec->seqs));
  for (i = 0; i < spec->n_strucs; i++) {
    free_struc(spec->strucs + i);
    free(spec->struc_names[i]);
  }

  for (i = 0; i < spec->n_tubes; i++) {
    free_tube(spec->tubes + i);
    free(spec->tube_names[i]);
  }

  for (i = 0; i < spec->n_orderings; i++) {
    free(spec->orderings[i]);
  }

  free(spec->struc_names);
  free(spec->strucs);
  free(spec->tube_names);
  free(spec->tubes);
  free(spec->orderings);
  free(spec->symmetry);
  free(spec->struc_map);
  free(spec->n_strands);
  spec->struc_names = NULL;
  spec->strucs = NULL;
  spec->tube_names = NULL;
  spec->tubes = NULL;
  spec->orderings = NULL;
  spec->symmetry = NULL;
  spec->n_strands = NULL;

  spec->cap_strucs = 0;
  spec->cap_tubes = 0;
  spec->cap_orderings = 0;

  spec->n_strucs = 0;
  spec->n_tubes = 0;
  spec->n_orderings = 0;
  return ERR_OK;
}

int free_tube (
      design_tube_t * tube
    ) {

  free(tube->target_x);
  free(tube->target);
  free(tube->included);
  free(tube->included_ind);
  free(tube->generated_ind);

  tube->target_x = NULL;
  tube->target = NULL;
  tube->included = NULL;
  tube->included_ind = NULL;
  tube->generated_ind = NULL;

  tube->maxsize = 0;
  tube->stop = 0;
  tube->n_strucs = 0;
  tube->cap_strucs = 0;
  return ERR_OK;
}

int copy_result_struc(
      result_struc_t * dest,
      result_struc_t * src
    ) {
  int n_nucs, n_strands;
  int i;
  int j;
  DBL_TYPE * temp_dbl;
  int * temp_int;

  n_nucs = src->n_nucs;
  n_strands = src->n_breaks + 1;

  check(ERR_OK == 
      ensure_capacity_vals(dest, n_nucs, n_strands),
      "Error allocating memory for struc result");
  for (i = 0; i < n_nucs; i ++) {
    dest->sequence[i] = src->sequence[i];
    dest->f_sequence[i] = src->f_sequence[i];
    dest->nuc_ids[i] = src->nuc_ids[i];
    dest->structure[i] = src->structure[i];
    dest->modifiable[i] = src->modifiable[i];
    dest->defects[i] = src->defects[i];
    dest->f_defects[i] = src->f_defects[i];
  }

  for (i = 0; i < n_strands - 1; i++) {
    dest->breaks[i] = src->breaks[i];
  }

  dest->pfunc = src->pfunc;
  dest->time = src->time;

  if (src->ppairs_n > dest->ppairs_cap) {
    temp_dbl = (DBL_TYPE*) 
      realloc(dest->ppairs, sizeof(DBL_TYPE) * src->ppairs_n);
    check_mem(temp_dbl);
    dest->ppairs = temp_dbl;

    temp_int = (int*) 
      realloc(dest->ppairs_i, sizeof(int) * src->ppairs_n);
    check_mem(temp_int);
    dest->ppairs_i = temp_int;

    temp_int = (int*) 
      realloc(dest->ppairs_j, sizeof(int) * src->ppairs_n);
    check_mem(temp_int);
    dest->ppairs_j = temp_int;
  }

  for (j = 0; j < src->ppairs_n; j++) {
    dest->ppairs[j] = src->ppairs[j];
    dest->ppairs_i[j] = src->ppairs_i[j];
    dest->ppairs_j[j] = src->ppairs_j[j];
  }
  dest->ppairs_n = src->ppairs_n;
  dest->ppairs_cap = src->ppairs_n;

  dest->defect = src->defect;

  if (src->tree) {
    if (!dest->tree) {
      dest->tree = alloc_empty_result_tree();
    }
    copy_result_tree(dest->tree, src->tree);
  } else {
    if (dest->tree) {
      destroy_result_tree(dest->tree);
      dest->tree = NULL;
    }
  }

  return ERR_OK;
error:
  return ERR_OOM;
}

int ensure_ppair_cap(
      result_struc_t * dest, int ppaircount
    ) {
  DBL_TYPE * temp_dbl;
  int * temp_int;

  if (dest->ppairs_cap < ppaircount) {
    temp_dbl = (DBL_TYPE*) realloc(dest->ppairs, 
        ppaircount * sizeof(DBL_TYPE));
    check_mem(temp_dbl);
    dest->ppairs = temp_dbl;

    temp_int = (int*) realloc(dest->ppairs_i,
        ppaircount * sizeof(int));
    check_mem(temp_int);
    dest->ppairs_i = temp_int;
    
    temp_int = (int*) realloc(dest->ppairs_j,
        ppaircount * sizeof(int));
    check_mem(temp_int);
    dest->ppairs_j = temp_int;
  }

  return ERR_OK;
error:
  return ERR_OOM;
}

int ensure_capacity_vals(
      result_struc_t * dest,
      int n_nucs,
      int n_strands
    ) {
  // Allocate memory for everything
  int i_nuc;
  
  if (n_nucs > dest->cap_nucs) {
    dest->sequence = (int*) realloc(dest->sequence, 
        n_nucs*sizeof(int));
    dest->nuc_ids = 
      (int*) realloc(dest->nuc_ids, n_nucs * sizeof(int));
    dest->f_sequence = 
      (int*) realloc(dest->f_sequence, n_nucs * sizeof(int));
    dest->modifiable = 
      (int*) realloc(dest->modifiable, n_nucs * sizeof(int));
    dest->structure = 
      (int*) realloc(dest->structure, n_nucs * sizeof(int));
    dest->defects = 
      (DBL_TYPE*) realloc(dest->defects, 
          n_nucs * sizeof(DBL_TYPE));
    dest->f_defects = 
      (DBL_TYPE*) realloc(dest->f_defects, 
          n_nucs * sizeof(DBL_TYPE));

    check_mem(dest->sequence);
    check_mem(dest->f_sequence);
    check_mem(dest->nuc_ids);
    check_mem(dest->structure);
    check_mem(dest->modifiable);
    check_mem(dest->defects);
    check_mem(dest->f_defects);

    for (i_nuc = dest->n_nucs; i_nuc < n_nucs; i_nuc++) {
      dest->sequence[i_nuc] = 0;
      dest->f_sequence[i_nuc] = 0;
      dest->nuc_ids[i_nuc] = 0;
      dest->structure[i_nuc] = 0;
      dest->modifiable[i_nuc] = 0;
      dest->defects[i_nuc] = 0;
      dest->f_defects[i_nuc] = 0;
    }
    dest->cap_nucs = n_nucs;
  }

  dest->n_nucs = n_nucs;

  if (n_strands - 1 > dest->n_breaks 
      || !dest->breaks) {
    dest->breaks = (int*) realloc(dest->breaks, 
          (n_strands - 1) * sizeof(int));
  }


  dest->n_breaks = n_strands - 1;
  return ERR_OK;
error:
  return ERR_OOM;
}

int add_nucleotides(
    design_spec_t * spec,
    int * start,
    int n) {
  int oldstart = spec->seqs.nucs.n;
  check(ERR_OK == resize_sequence(&(spec->seqs.nucs), n + oldstart),
      "Error resizing nucs");
  check(ERR_OK == resize_sequence(&(spec->seqs.comp_map), n + oldstart),
      "Error resizing comp_map");
  int i;
  for (i = oldstart; i < oldstart + n; i++) {
    spec->seqs.comp_map.nucs[i] = -1;
  }
  *start = oldstart;
  return ERR_OK;
error:
  return ERR_OOM;
}

int add_strand_basic(
      design_spec_t * spec,
      int * ind,
      char * name,
      int len
    ) {
  int start = 0;
  sequence_t seq;
  int i;
  int cur_n = spec->seqs.strands.n;
  // Add the requisite nucleotides
  add_nucleotides(spec, &start, len);
  // Add the strand

  *ind = 0;
  init_sequence(&seq, len);
  for (i = 0; i < len; i++) {
    seq.nucs[i] = start + i;
  }
  check(ERR_OK == append_sequence_seqspec(&(spec->seqs.strands), &seq, name),
      "Error appending sequence specification %i", cur_n);
  *ind = cur_n;
  free_sequence(&seq);
  return ERR_OK;
error:
  free_sequence(&seq);
  return ERR_INVALID_STATE;
}

int add_strand(
      design_spec_t * spec,
      char * name,
      char ** domain_names,
      int n_domains
    ) {

  int i_dom;
  int c_dom;
  int i_nuc;
  int j_nuc;
  int n_nucs = 0;
  sequence_t seq;

  for (i_dom = 0; i_dom < n_domains; i_dom++) {
    c_dom = index_of(domain_names[i_dom], 
        spec->seqs.domains.names, spec->seqs.domains.n);

    check(c_dom >= 0, "invalid domain name: %s", domain_names[i_dom]);
    n_nucs += spec->seqs.domains.specs[c_dom].n;
  }

  init_sequence(&seq, n_nucs);

  i_nuc = 0;
  for (i_dom = 0; i_dom < n_domains; i_dom++) {
    c_dom = index_of(domain_names[i_dom], 
        spec->seqs.domains.names, spec->seqs.domains.n);
    check(c_dom >= 0, "invalid domain name: %s", domain_names[i_dom]);
    n_nucs = spec->seqs.domains.specs[c_dom].n;

    for (j_nuc = 0; j_nuc < n_nucs; j_nuc++) {
      seq.nucs[i_nuc] =
        spec->seqs.domains.specs[c_dom].nucs[j_nuc];
      i_nuc += 1;
    }
  }

  check(ERR_OK == append_sequence_seqspec(&(spec->seqs.strands), &seq, name),
      "Error appending sequence specification %s", name);
  free_sequence(&seq);

  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

static int get_symmetry(
      int * ids, int n_ids
    ) {
  int i, j;
  int cur_match;
  for (i = 1; i < n_ids; i++) {
    cur_match = 1;
    for (j = 0; j < n_ids; j++) {
      if (ids[j] != ids[(i + j) % n_ids]) {
        cur_match = 0;
      }
    }

    if (cur_match) {
      break;
    }
  }

  return n_ids / i;
}

int define_structure_strands(
      design_spec_t * spec,
      char * strucname,
      char ** strandnames,
      int n_strands
    ) {
  
  int g_struc;
  int c_struc;
  int i_str;
  int c_str;
  int o_str;

  c_struc = index_of(strucname, spec->struc_names, spec->n_strucs);
  check(c_struc >= 0, "Unknown structure %s", strucname);
  check(spec->strucs[c_struc].n_strands == n_strands,
      "Mismatch in strand count %i in structure %i specified",
      spec->strucs[c_struc].n_strands, n_strands);
  for (i_str = 0; i_str < n_strands; i_str++) {
    o_str = spec->strucs[c_struc].strands[i_str];
    check(spec->seqs.strands.names[o_str][0] == '-', 
        "sequences for structure %s defined twice.", strucname);
    c_str = index_of(strandnames[i_str], 
        spec->seqs.strands.names, spec->seqs.strands.n);
    check(c_str >= 0, "strand %s not defined", strandnames[i_str]);
    spec->strucs[c_struc].strands[i_str] = c_str;
  }
  for (g_struc = 0; g_struc < spec->n_orderings; g_struc++) {
    if (c_struc == spec->struc_map[g_struc]) {
      break;
    }
  }

  check(g_struc < spec->n_orderings, "Ordering for structure %s not found",
      strucname);
  check(spec->n_strands[g_struc] == n_strands,
      "Ordering for structure %s has the incorrect number of strands %i %i",
      strucname, n_strands, spec->n_strands[g_struc]);
  for (i_str = 0; i_str < n_strands; i_str++) {
    spec->orderings[g_struc][i_str] = 
      spec->strucs[c_struc].strands[i_str];
  }

  spec->symmetry[g_struc] = get_symmetry(spec->strucs[c_struc].strands, 
      n_strands);
  spec->strucs[c_struc].symmetry = spec->symmetry[g_struc];

  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

int add_domain(
      design_spec_t * spec,
      char * name, 
      char * constraint
    ) {
  sequence_t seq;
  int start = 0;
  int n_nucs = strlen(constraint);
  int i;
  int cur_n = spec->seqs.domains.n;
  int namelen = strlen(name);
  char * compname = NULL;

  int i_dom = index_of(name, spec->seqs.domains.names, spec->seqs.domains.n);
  // int j_dom = -1;

  compname = (char *) malloc((namelen + 2) * sizeof(char));
  strncpy(compname, name, namelen + 1);
  if (name[namelen - 1] == '*') {
    compname[namelen - 1] = '\0';
  } else {
    compname[namelen] = '*';
    compname[namelen + 1] = '\0';
  }
  seq.nucs = NULL;
  seq.n = 0;

  // j_dom = index_of(compname, spec->seqs.domains.names, spec->seqs.domains.n);

  if (i_dom >= 0) {
    sentinel("Complement or duplicate assignments not currently allowed");
    // for (i = 0; i < n_nucs; i++) {

    // }
  } else {
    check(ERR_OK == add_nucleotides(spec, &start, n_nucs * 2),
        "Error adding nucleotides for domain");

    for (i = 0; i < n_nucs; i++) {
      spec->seqs.nucs.nucs[start + i] = get_nuc_code(constraint[i]);
      spec->seqs.comp_map.nucs[start + i] = start + 2 * n_nucs - i - 1;
    }
    check(ERR_OK == init_sequence(&seq, n_nucs),
        "Error initializing domain");
    for (i = 0; i < n_nucs; i++) {
      seq.nucs[i] = start + i;
    }

    check(ERR_OK == append_sequence_seqspec(&(spec->seqs.domains), 
          &seq, name),
      "Error appending domain specification %i", cur_n);

    for (i = 0; i < n_nucs; i++) {
      spec->seqs.nucs.nucs[start + n_nucs + i] = 
        get_comp_code(get_nuc_code(constraint[n_nucs - i - 1]), spec->opts.allow_wobble);
      spec->seqs.comp_map.nucs[start + n_nucs + i] = start + n_nucs - i - 1;
    }

    for (i = 0; i < n_nucs; i++) {
      seq.nucs[i] = start + n_nucs + i;
    }
    check(ERR_OK == append_sequence_seqspec(&(spec->seqs.domains), 
          &seq, compname),
      "Error appending domain specification %i", cur_n);

  }

  free(compname);
  free_sequence(&seq);
  return ERR_OK;
error:
  free(compname);
  free_sequence(&seq);
  return ERR_INVALID_STATE;
}

int index_of(const char * name, char ** names, int n_names) {
  int i;
  check(name != NULL, "Invalid null string in name");
  check(names != NULL || n_names == 0, "Invalid null list for names");
  for (i = 0; i < n_names; i++) {
    if (strcmp(name, names[i]) == 0) {
      return i;
    }
  }
error:
  return -1;
}

int set_concentration(
      design_spec_t * spec, 
      const char * tube_name,
      const char * struc_name,
      DBL_TYPE conc
    ) {
  int c_struc = 0;
  int i_struc = 0;
  c_struc = index_of(struc_name, spec->struc_names, spec->n_strucs);
  for (i_struc = 0; i_struc < spec->tubes[0].n_strucs; i_struc++) {
    if (spec->tubes[0].included_ind[i_struc] == c_struc) {
      spec->tubes[0].target_x[i_struc] = conc;
      break;
    }
  }
  check(i_struc < spec->tubes[0].n_strucs, 
      "Tube %s does not contain structure %s", tube_name, struc_name);
  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

int add_altstruc_basic(
      design_spec_t * spec,
      char * name, 
      char * struc
    ) {
  int n_chars = strlen(struc);
  int namelen = strlen(name);
  int * pairing = (int*) malloc(n_chars * sizeof(int));
  int * breaks = (int*) malloc(n_chars * sizeof(int));
  int n_nucs = 0;
  int n_breaks = 0;
  int i_str;

  check_mem(pairing);
  check_mem(breaks);

  get_pairing(pairing, breaks, &n_nucs, &n_breaks, struc, n_chars);

  for (i_str = 0; i_str < spec->n_strucs; i_str++) {
    if (0 == strncmp(name, spec->struc_names[i_str], namelen + 1)) {
      break;
    }
  }

  check(i_str < spec->n_strucs, "Structure %s not found in altstruc assignment",
     name); 

  check(n_breaks == spec->strucs[i_str].n_strands - 1, 
      "Invalid number of breaks in altstruc for structure %s", name);
  
  check(n_nucs == spec->strucs[i_str].n_nucs,
      "Invalid length of altstruc for %s, new: %i old: %i", name, 
      n_nucs, spec->strucs[i_str].n_nucs);

  check(ERR_OK == add_altstruc(spec, i_str, pairing, n_nucs),
      "Error adding altstruc for %s", name);

  return ERR_OK;
error:
  free(pairing);
  free(breaks);
  return ERR_INVALID_STATE;
}

int add_altstruc(
      design_spec_t * spec,
      int c_struc,
      int * pairing,
      int n_nucs
    ) {
  
  design_struc_t * head = &(spec->strucs[c_struc]);
  design_struc_t ** next_p = &(spec->strucs[c_struc].next);

  while (*next_p) {
    next_p = &((*next_p)->next);
  }

  design_struc_t * newstruc = (design_struc_t*) malloc(sizeof(design_struc_t));

  check(ERR_OK == init_struc(newstruc, &(spec->seqs), head->strands, 
        head->n_strands, pairing, n_nucs),
      "Error initializing structure");


  newstruc->symmetry = head->symmetry;

  *next_p = newstruc;
  return ERR_OK;
error:
  if (newstruc) free_struc(newstruc);
  free(newstruc);
  return ERR_INVALID_STATE;
}

int add_structure_basic(
      design_spec_t * spec,
      char * name,
      char * struc
    ) {
  int n_chars = strlen(struc);
  int * pairing = (int*) malloc(n_chars * sizeof(int));
  int * breaks = (int*) malloc(n_chars * sizeof(int));
  int * strands = NULL;
  int n_nucs = 0;
  int n_breaks = 0;
  int i_str = 0;
  int c_str;
  int i_nuc;
  int lastbreak = 0;
  int namelen = strlen(name);
  char * strandname = (char*) malloc((namelen + 20) * sizeof(char));

  check_mem(pairing);
  check_mem(breaks);
  check_mem(strandname);

  get_pairing(pairing, breaks, &n_nucs, &n_breaks, struc, n_chars);
  strands = (int*) malloc((n_breaks + 1) * sizeof(int));
  check_mem(strands);
  check(n_breaks < 10000, "Too many breaks");
  

  i_nuc = 0;
  for (i_str = 0; i_str < n_breaks; i_str++) {
    snprintf(strandname, namelen + 19, "-%s_strand_%i", name, i_str + 1);
    strandname[namelen + 19] = '\0';
    add_strand_basic(spec, &c_str, strandname, breaks[i_str] - i_nuc);
    strands[i_str] = c_str;
    i_nuc = breaks[i_str];
  }
  if (n_breaks > 0) {
    lastbreak = breaks[n_breaks - 1];
  }

  snprintf(strandname, namelen + 19, "-%s_strand_%i", name, n_breaks + 1);
  strandname[namelen + 19] = '\0';
  add_strand_basic(spec, &c_str, strandname, n_nucs - lastbreak);
  strands[n_breaks] = c_str;

  check(ERR_OK == 
      add_structure(spec, name, strands, n_breaks+1, pairing, n_nucs),
      "Error adding structure: %s", name);

  free(strandname);
  free(pairing);
  free(breaks);
  free(strands);

  return ERR_OK;
error:
  free(strandname);
  free(pairing);
  free(breaks);
  free(strands);
  return ERR_OOM;
}


int add_structure_str(
      design_spec_t * spec,
      char * name,
      char ** strand_names, 
      int n_strands,
      char * domain_struc
    ) {

  int * strands = NULL;
  int * domain_pairing = NULL;
  int * breaks = NULL;
  int n_tot_strands = spec->seqs.strands.n;
  int i_str;
  int j_str;
  int n_doms = 0;
  int n_pos = strlen(domain_struc);
  int n_breaks;
  char ** name_lookup = spec->seqs.strands.names;

  strands = (int*) malloc(n_strands * sizeof(int));
  domain_pairing = (int*) malloc(n_pos * sizeof(int));
  breaks = (int*) malloc(n_pos * sizeof(int));

  get_pairing(domain_pairing, breaks, &n_doms, &n_breaks, 
      domain_struc, n_pos);

  // Look up strand indices
  for (i_str = 0; i_str < n_strands; i_str++) {
    for (j_str = 0; j_str < n_tot_strands; j_str++) {
      if (strcmp(strand_names[i_str], name_lookup[j_str]) == 0) {
        strands[i_str] = j_str;
        break;
      }
    }
    check(j_str < n_tot_strands, 
        "Error in structure %s. Strand %s not found",
        name, strand_names[i_str]);
  }

  add_structure(spec, name, strands, n_strands, domain_pairing, n_doms);

  free(strands);
  free(domain_pairing);
  free(breaks);

  return ERR_OK;
error:
  free(strands);
  free(domain_pairing);
  free(breaks);

  return ERR_INVALID_STATE;
}

int add_structure(
      design_spec_t * spec,
      char * name,
      int * strands,
      int n_strands,
      int * dom_struc,
      int n_domains
    ) {

  int n_strucs = spec->n_strucs;
  int i_strand;
  int j_strand;
  int name_len = strlen(name) + 1;
  int cap_strucs = spec->cap_strucs;
  int cap_ords;
  int n_ords;
  int i_ord;
  int c_ord = -1;
  int * temp_n_strands = NULL;
  int * temp_struc_map = NULL;
  int * temp_symmetry = NULL;
  int ** temp_orderings = NULL;
  char ** temp_names = NULL;
  int order_present = 0;
  int orders_equal = 0;
  int i_tube;
  int n_tubes;

  design_struc_t * temp_strucs = NULL;

  if (n_strucs + 1 >= cap_strucs) {
    cap_strucs = 2 * (n_strucs + 1);
    temp_strucs = (design_struc_t *) realloc(spec->strucs,
        cap_strucs * sizeof(design_struc_t));
    check_mem(temp_strucs);
    spec->strucs = temp_strucs;

    temp_names = (char **) realloc(spec->struc_names,
        cap_strucs * sizeof(char*));
    check_mem(temp_names);
    spec->struc_names = temp_names;

    spec->cap_strucs = cap_strucs;
  }

  spec->struc_names[n_strucs] = (char *) malloc(name_len * sizeof(char));
  strncpy(spec->struc_names[n_strucs], name, name_len);
  
  init_struc(spec->strucs + n_strucs, &(spec->seqs), strands, n_strands,
      dom_struc, n_domains);

  cap_ords = spec->cap_orderings;
  n_ords = spec->n_orderings;
  // Search for the structure in the current set of orderings
  order_present = 0;
  for (i_ord = 0; i_ord < n_ords; i_ord++) {
    if (n_strands == spec->n_strands[i_ord]) {
      orders_equal = 1;
      for (i_strand = 0; i_strand < n_strands; i_strand++) {
        orders_equal = 1;
        for (j_strand = 0; j_strand < n_strands; j_strand++) {
          if (strands[j_strand] != 
              spec->orderings[i_ord][(i_strand + j_strand) % n_strands]) {
            orders_equal = 0;
            break;
          }
        }
        if (orders_equal) {
          order_present = 1;
          c_ord = i_ord;
          break;
        }
      }
      if (orders_equal) {
        break;
      }
    }
  }
  if (order_present) {
    check(spec->struc_map[c_ord] == -1, 
        "Structure %s is being added twice", name);
    spec->struc_map[c_ord] = n_strucs;
    n_tubes = spec->n_tubes;
    for (i_tube = 0; i_tube < n_tubes; i_tube++) {
      n_ords = spec->tubes[i_tube].n_strucs;
      for (i_ord = 0; i_ord < n_ords; i_ord++) {
        if (spec->tubes[i_tube].generated_ind[i_ord] == c_ord) {
          spec->tubes[i_tube].included[i_ord] = 1;
          spec->tubes[i_tube].included_ind[i_ord] = n_strucs;
        }
      }
    }
    spec->strucs[n_strucs].symmetry = spec->symmetry[c_ord];

  } else {
    // Add the on-target orderings
    if (n_ords + 1 >= cap_ords) {
      cap_ords = 2*(spec->n_orderings + 1);
      temp_orderings = (int**) realloc(spec->orderings, 
          sizeof(int*) * cap_ords);
      check_mem(temp_orderings);
      spec->orderings = temp_orderings;

      temp_n_strands = (int*) realloc(spec->n_strands,
          sizeof(int) * cap_ords);
      check_mem(temp_n_strands);
      spec->n_strands = temp_n_strands;

      temp_struc_map = (int*) realloc(spec->struc_map,
          sizeof(int) * cap_ords);
      check_mem(temp_struc_map);
      spec->struc_map = temp_struc_map;

      temp_symmetry = (int*) realloc(spec->symmetry,
          sizeof(int) * cap_ords);
      check_mem(temp_symmetry);
      spec->symmetry = temp_symmetry;

      for (i_ord = n_ords; i_ord < cap_ords; i_ord++) {
        spec->orderings[i_ord] = NULL;
        spec->n_strands[i_ord] = 0;
        spec->struc_map[i_ord] = -1;
        spec->symmetry[i_ord] = 1;
      }
      spec->cap_orderings = cap_ords;
    }

    spec->orderings[n_ords] = (int*) malloc(n_strands * sizeof(int));


    for (i_strand = 0; i_strand < n_strands; i_strand++) {
      spec->orderings[n_ords][i_strand] = strands[i_strand];
    }

    for (i_strand = 1; i_strand < n_strands; i_strand++) {
      order_present = 1;
      for (j_strand = 0; j_strand < n_strands; j_strand++) {
        if (strands[(i_strand + j_strand) % n_strands] != strands[j_strand]) {
          order_present = 0;
          break;
        }
      }
      if (order_present) {
        break;
      }
    }


    spec->n_strands[n_ords] = n_strands;
    spec->struc_map[n_ords] = n_strucs;
    spec->symmetry[n_ords] = n_strands / i_strand;

    spec->strucs[n_strucs].symmetry = spec->symmetry[n_ords];

    spec->n_orderings = n_ords + 1;
  }

  spec->n_strucs ++;
  
  return ERR_OK;
error:
  return ERR_OOM;
}

int init_struc(
      design_struc_t * struc,
      sequence_spec_t * seqs,
      int * strands,
      int n_strands,
      int * nucstruc,
      int n_s_nucs
    ) {
  int i_str;
  int i_nuc;
  int n_nucs;

  struc->struc = NULL;
  struc->strands = NULL;
  struc->split_forbidden = NULL;
  struc->next = NULL;
  struc->tree = NULL;

  check(n_strands > 0, "Number of strands %i <= 0", n_strands);
  check(n_s_nucs > 0, "Number of domains %i <= 0", n_s_nucs);

  struc->strands = (int*) malloc(n_strands * sizeof(int));
  check_mem(struc->strands);


  n_nucs = 0;
  for (i_str = 0; i_str < n_strands; i_str++) {
    struc->strands[i_str] = strands[i_str];
    n_nucs += seqs->strands.specs[strands[i_str]].n;
  }

  check(n_s_nucs == n_nucs,
      "# str: %i # nucs in struc: %i != # nucs in strand(s): %i", 
      n_strands, n_s_nucs, n_nucs);

  struc->split_forbidden = (int*) malloc(n_nucs * sizeof(int));
  struc->struc = (int*) malloc(n_nucs * sizeof(int));
  check_mem(struc->split_forbidden);
  check_mem(struc->struc);
  for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    struc->struc[i_nuc] = nucstruc[i_nuc];
    struc->split_forbidden[i_nuc] = 0;
  }

  struc->n_strands = n_strands;
  struc->n_nucs = n_nucs;
  struc->symmetry = 1;
  struc->modifiable = 0;

  i_nuc = 0;
  struc->tree = alloc_struc_tree(1, &i_nuc, &n_nucs, 0, NULL, NULL); 

  return ERR_OK;
error:
  if (struc->strands) free(struc->strands);
  if (struc->struc) free(struc->struc);
  if (struc->split_forbidden) free(struc->split_forbidden);
  if (struc->tree) free_struc_tree(struc->tree);

  return ERR_OOM;
}

int add_tube_str(
      design_spec_t * spec,
      char * name, 
      char ** strucs,
      DBL_TYPE * concs,
      int n_strucs,
      DBL_TYPE stop,
      int maxsize
    ) {

  int tot_n_strucs = spec->n_strucs;
  int i_str;
  int j_str;
  int namelen;
  
  int * struc_ints = (int*) malloc(n_strucs * sizeof(int));
  check_mem(struc_ints);

  for (i_str = 0; i_str < n_strucs; i_str++) {
    namelen = strlen(strucs[i_str]) + 1;
    for (j_str = 0; j_str < tot_n_strucs; j_str++) {
      if (strncmp(strucs[i_str], spec->struc_names[j_str], namelen) == 0) {
        struc_ints[i_str] = j_str;
        break;
      }
    }
    check(j_str < tot_n_strucs, "Tube %s: Couldn't find structure %s",
        name, strucs[i_str]);
  }

  add_tube(spec, name, struc_ints, concs, n_strucs, stop, maxsize);

  free(struc_ints);
  return ERR_OK;
error:
  free(struc_ints);
  return ERR_INVALID_STATE;
}

static int generate_complex_orders(
      int *** result_p,
      int ** n_comp_strands,
      int * n_results_p,
      int * cap_results_p,
      int * strand_ids,
      int n_strands,
      int maxsize
    ) {
  
  int n;
  int i;
  int j;
  int k;
  // int l;
  // int same;
  int ** res = *result_p;
  int * n_comp_strs = *n_comp_strands;
  int cap = *cap_results_p;
  int n_res = *n_results_p;

  int ** temp_res = NULL;
  int * temp_comp_strs = NULL;

  int * a = (int *) malloc(maxsize * sizeof(int));
  int * tempa = NULL;

  for (n = 1; n <= maxsize; n++) {
    for (i = 0; i < n; i++) {
      a[i] = 0;
    }
    if (n_res + 1 > cap) {
      cap = ((n_res + 1) * 2);
      temp_res = (int**) realloc(res, sizeof(int*) * cap);
      check_mem(temp_res);
      res = temp_res;

      temp_comp_strs = (int*) realloc(n_comp_strs, sizeof(int) * cap);
      check_mem(temp_comp_strs);
      n_comp_strs = temp_comp_strs;
    }

    tempa = (int*) malloc(n * sizeof(int));
    for (k = 0; k < n; k++) {
      check(a[k] < n_strands, "a[%i] = %i n_strands %i", 
          k, a[k], n_strands);
      tempa[k] = strand_ids[a[k]];
    }
    res[n_res] = tempa;
    tempa = NULL;

    n_comp_strs[n_res] = n;

    n_res++;
    
    i = n;
    if (n_strands == 1) {
      i = 0;
    }
    while (i > 0) {
      a[i-1] = a[i-1] + 1;
      for (j = 1; j < n - i + 1; j++) {
        a[i + j - 1] = a[j - 1];
      }

      if (n % i == 0) {
        if (n_res + 1 > cap) {
          cap = ((n_res + 1) * 2);
          temp_res = (int**) realloc(res, sizeof(int*) * cap);
          check_mem(temp_res);
          res = temp_res;

          temp_comp_strs = (int*) realloc(n_comp_strs, sizeof(int) * cap);
          check_mem(temp_comp_strs);
          n_comp_strs = temp_comp_strs;
        }

        tempa = (int*) malloc(n * sizeof(int));
        for (k = 0; k < n; k++) {
          check(a[k] < n_strands, "a[%i] = %i n_strands %i", 
              k, a[k], n_strands);
          tempa[k] = strand_ids[a[k]];
        }
        res[n_res] = tempa;
        tempa = NULL;

        n_comp_strs[n_res] = n;
        n_res++;
      }

      for (i = n; i > 0; i --) {
        if (a[i-1] != n_strands - 1) {
          break;
        }
      }
    }
  }

  // n = 1;
  // for (i = 1; i < n_res; i++) {
  //   for (j = 0; j < n; j++) {
  //     same = 0;
  //     if (n_comp_strs[i] == n_comp_strs[j]) {
  //       for (k = 0; k < n_comp_strs[i]; k++) {
  //         same = 1;
  //         for (l = 0; l < n_comp_strs[j]; l++) {
  //           if (res[i][(k + l) % n_comp_strs[i]] != res[j][l]) {
  //             same = 0;
  //           }
  //         }
  //         if (same) {
  //           break;
  //         }
  //       }
  //     }
  //     if (same) {
  //       break;
  //     }
  //   }
  //   if (same) {
  //     debug("Duplicate found: %i", i);
  //     free(res[i]);
  //     res[i] = NULL;
  //     n_comp_strs[i] = 0;
  //   } else {
  //     res[n] = res[i];
  //     n_comp_strs[n] = n_comp_strs[i];
  //     n += 1;
  //   }

  // }

  free(a);

  *result_p = res;
  *n_comp_strands = n_comp_strs;
  *cap_results_p = cap;
  *n_results_p = n_res;

  return ERR_OK;
error:
  free(tempa);
  free(a);
  *result_p = res;
  *n_comp_strands = n_comp_strs;
  *cap_results_p = cap;
  *n_results_p = n_res;
  return ERR_OOM;
}

int add_tube_basic(
      design_spec_t * spec,
      char * name,
      char ** strucnames,
      int n_strucs
    ) {
  DBL_TYPE * concs = (DBL_TYPE *) malloc(n_strucs * sizeof(DBL_TYPE));
  DBL_TYPE water_conc = water_density(
      (double)(spec->opts.temperature - ZERO_C_IN_KELVIN));
  int i;
  for (i = 0; i < n_strucs; i++) {
    concs[i] = DEFAULT_CONCENTRATION / water_conc;
  }

  add_tube_str(spec, name, strucnames, concs, n_strucs, DEFAULT_STOP_CONDITION,
      DEFAULT_MAXSIZE);

  free(concs);

  return ERR_OK;
}


int add_tube(
      design_spec_t * spec,
      char * name,
      int * strucs,
      DBL_TYPE * concs,
      int n_strucs,
      DBL_TYPE stop,
      int maxsize
    ) {
  int i_tube = spec->n_tubes;
  int namelen = strlen(name);
  int i_struc;
  int c_struc;
  int g_struc;
  int n_ords;
  int i_ord;
  debug("Adding tube %s", name);

  if (spec->cap_tubes == 0) {
    spec->cap_tubes = 1;
    spec->tubes = (design_tube_t *) realloc(spec->tubes, 
        sizeof(design_tube_t) * spec->cap_tubes);
    spec->tube_names = (char**) realloc(spec->tube_names,
        sizeof(char *) * spec->cap_tubes);
    check_mem(spec->tubes);
    check_mem(spec->tube_names);
  }

  spec->tube_names[i_tube] = (char*) malloc((namelen + 1) * sizeof(char));
  strncpy(spec->tube_names[i_tube], name, namelen + 1);
  
  spec->tubes[i_tube].target_x = (DBL_TYPE *) 
    malloc(n_strucs * sizeof(DBL_TYPE));

  spec->tubes[i_tube].target = (int*) malloc(n_strucs * sizeof(int));
  spec->tubes[i_tube].included = (int*) malloc(n_strucs * sizeof(int));
  spec->tubes[i_tube].included_ind = (int*) malloc(n_strucs * sizeof(int));
  spec->tubes[i_tube].generated_ind = (int*) malloc(n_strucs * sizeof(int));

  spec->tubes[i_tube].maxsize = maxsize;
  spec->tubes[i_tube].stop = stop;
  spec->tubes[i_tube].n_strucs = n_strucs;
  spec->tubes[i_tube].cap_strucs = n_strucs;
  
  n_ords = spec->n_orderings;

  for (i_struc = 0; i_struc < n_strucs; i_struc++) {
    c_struc = strucs[i_struc];
    spec->tubes[i_tube].target_x[i_struc] = concs[i_struc];
    spec->tubes[i_tube].target[i_struc] = concs[i_struc] > 0;
    spec->tubes[i_tube].included[i_struc] = 1;
    spec->tubes[i_tube].included_ind[i_struc] = c_struc;
    spec->tubes[i_tube].generated_ind[i_struc] = -1; 
    g_struc = -1;

    for (i_ord = 0; i_ord < n_ords; i_ord++) {
      if (spec->struc_map[i_ord] == c_struc) {
        g_struc = i_ord;
        break;
      }
    }
    spec->tubes[i_tube].generated_ind[i_struc] = g_struc;

    check(spec->tubes[i_tube].generated_ind[i_struc] >= 0, 
        "Couldn't find ordering");
  }

  spec->tubes[i_tube].maxsize = maxsize;
  spec->tubes[i_tube].stop = stop;
  spec->tubes[i_tube].n_strucs = n_strucs;
  spec->tubes[i_tube].cap_strucs = n_strucs;
  debug("Total complexes in tube %i", n_strucs);
  spec->n_tubes = i_tube + 1;
  return ERR_OK;
error:

  return ERR_OOM;
}

int make_off_targets(design_spec_t * spec) {
  int i_struc;
  int c_struc;
  int n_ords = spec->n_orderings;
  int i_ord;
  int i_str;
  int j_str;
  int n_strs;
  int n_strucs ;
  int tot_n_strucs;
  int * strands_included = (int*) malloc(sizeof(int) * spec->seqs.strands.n);
  int cap_orders = 100;
  int ** new_orders = (int**) malloc(sizeof(int*) * cap_orders);
  int * new_n_strands = (int*) malloc(sizeof(int) * cap_orders);
  DBL_TYPE * temp_dbl;
  int * temp_int;

  int ** temp_orderings;
  int * temp_strs;
  int * temp_map;
  int * temp_sym;

  int n_orders = 0;
  int n_included = 0;
  int order_present;
  int repl_ind;

  n_strucs = spec->tubes[0].n_strucs;

  for (i_struc = 0; i_struc < n_strucs; i_struc++) {
    if (spec->tubes[0].included[i_struc] 
        && spec->tubes[0].target_x[i_struc] > 0) {
      c_struc = spec->tubes[0].included_ind[i_struc];
      n_strs = spec->strucs[c_struc].n_strands;
      for (i_str = 0; i_str < n_strs; i_str++) {
        for (j_str = 0; j_str < n_included; j_str++) {
          if (strands_included[j_str] == spec->strucs[c_struc].strands[i_str]){
            break;
          }
        }
        if (j_str == n_included) {
          strands_included[n_included] = spec->strucs[c_struc].strands[i_str];
          n_included ++;
        }
      }
    }
  }

  generate_complex_orders(&new_orders, &new_n_strands, &n_orders, &cap_orders, 
      strands_included, n_included, spec->tubes[0].maxsize);

  tot_n_strucs = n_strucs + n_orders;

  temp_dbl = (DBL_TYPE *) realloc(spec->tubes[0].target_x, 
      tot_n_strucs * sizeof(DBL_TYPE));
  check_mem(temp_dbl);
  spec->tubes[0].target_x = temp_dbl;
  
  temp_int = (int*) realloc(spec->tubes[0].target, 
      tot_n_strucs * sizeof(int));
  check_mem(temp_int);
  spec->tubes[0].target = temp_int;

  temp_int = (int*) realloc(spec->tubes[0].included, 
      tot_n_strucs * sizeof(int));
  check_mem(temp_int);
  spec->tubes[0].included = temp_int;


  temp_int = (int*) realloc(spec->tubes[0].included_ind, 
      tot_n_strucs * sizeof(int));
  check_mem(temp_int);
  spec->tubes[0].included_ind = temp_int;

  temp_int = (int*) realloc(spec->tubes[0].generated_ind, 
      tot_n_strucs * sizeof(int));
  check_mem(temp_int);
  spec->tubes[0].generated_ind = temp_int;

 
  for (i_struc = 0; i_struc < n_orders; i_struc++) {
    n_strs = new_n_strands[i_struc];
    check(n_strs > 0, "Invalid number of strands %i", 
        new_n_strands[i_struc]);
    for (i_ord = 0; i_ord < n_ords; i_ord++) {
      if (spec->n_strands[i_ord] == n_strs) {
        for (i_str = 0; i_str < n_strs; i_str++) {
          for (j_str = 0; j_str < n_strs; j_str++) {
            if (spec->orderings[i_ord][j_str] != 
                new_orders[i_struc][(i_str + j_str) % n_strs]) {
              break;
            }
          }
          if (j_str == n_strs) {
            break;
          }
        }

        if (i_str != n_strs) {
          break;
        }
      }
    }

    if (i_ord == n_ords) {
      if (spec->cap_orderings < n_ords + 1) {
        spec->cap_orderings = (n_ords + 1) * 2;
        temp_orderings = (int**) realloc(spec->orderings, spec->cap_orderings 
            * sizeof(int*));
        check_mem(temp_orderings);
        spec->orderings = temp_orderings;

        temp_strs = (int*) realloc(spec->n_strands, spec->cap_orderings 
            * sizeof(int));
        check_mem(temp_strs);
        spec->n_strands = temp_strs;

        temp_map = (int*) realloc(spec->struc_map, spec->cap_orderings 
            * sizeof(int));
        check_mem(temp_map);
        spec->struc_map = temp_map;
        
        temp_sym = (int*) realloc(spec->symmetry, spec->cap_orderings 
            * sizeof(int));
        check_mem(temp_sym);
        spec->symmetry = temp_sym;
      }
      repl_ind = n_strs;

      for (i_str = 1; i_str < n_strs; i_str++) {
        order_present = 1;
        for (j_str = 0; j_str < n_strs; j_str++) {
          if (new_orders[i_struc][(i_str + j_str)% n_strs] != 
              new_orders[i_struc][j_str]) {
            order_present = 0;
            break;
          } 
        }
        if (order_present) {
          repl_ind = i_str;
          break;
        }
      }

      spec->symmetry[n_ords] = n_strs / repl_ind;
      spec->orderings[n_ords] = new_orders[i_struc];
      spec->n_strands[n_ords] = n_strs;
      spec->struc_map[n_ords] = -1;

      spec->tubes[0].target_x[n_strucs] = 0.0;
      spec->tubes[0].target[n_strucs] = 0;
      spec->tubes[0].included[n_strucs] = 0;
      spec->tubes[0].included_ind[n_strucs] = -1;
      spec->tubes[0].generated_ind[n_strucs] = n_ords;

      n_strucs ++;
      n_ords ++;
      spec->n_orderings = n_ords;
    } else {
      for (i_str = 0; i_str < n_strucs; i_str++) {
        if (spec->tubes[0].generated_ind[i_str] == i_ord) {
          break;
        }
      } 
      if (i_str == n_strucs) {
        spec->tubes[0].target_x[n_strucs] = 0.0;
        spec->tubes[0].target[n_strucs] = 0;
        spec->tubes[0].included[n_strucs] = 0;
        spec->tubes[0].included_ind[n_strucs] = -1;
        spec->tubes[0].generated_ind[n_strucs] = i_ord;
        n_strucs ++;
      }
      free(new_orders[i_struc]);
    }
  }
  spec->tubes[0].n_strucs = n_strucs;

  free(new_orders);
  free(new_n_strands);
  free(strands_included);

  return ERR_OK;
error:
  free(new_orders);
  free(new_n_strands);
  free(strands_included);
  return ERR_INVALID_STATE;
}


int init_seqs_random(
      seqstate_t * seqs,
      sequence_spec_t * seqspec,
      options_t * opts
    ) {

  int n_nucs, i_nuc, c_nuc;
  int * f_seq, *f_con, *c_seq, *c_con;

  n_nucs = seqspec->nucs.n;
  f_seq = NULL;
  f_con = NULL;
  c_seq = NULL;
  c_con = NULL;

  f_seq = (int*) malloc(n_nucs * sizeof(int));
  f_con = (int*) malloc(n_nucs * sizeof(int));
  c_seq = (int*) malloc(n_nucs * sizeof(int));
  c_con = (int*) malloc(n_nucs * sizeof(int));

  for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    f_seq[i_nuc] = 0;
    c_seq[i_nuc] = 0;

    f_con[i_nuc] = seqspec->nucs.nucs[i_nuc];
    if (seqspec->comp_map.nucs[i_nuc] < 0) {
      c_con[i_nuc] = BASE_N;
    } else {
      c_con[i_nuc] = seqspec->nucs.nucs[seqspec->comp_map.nucs[i_nuc]];
    }
  }

  check(ERR_OK == 
      init_constraint_random(f_seq, c_seq, f_con, c_con, n_nucs, opts),
      "Error initializing constraints");

  for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    seqs->nucspec.nucs[i_nuc] = f_seq[i_nuc];
    seqs->dumspec.nucs[i_nuc] = f_seq[i_nuc];
    check(f_seq[i_nuc] >= 1 && f_seq[i_nuc] <= 4, 
        "Invalid nucleotide id: %i nucleotide code: %i", i_nuc, f_seq[i_nuc]);
    c_nuc = seqspec->comp_map.nucs[i_nuc];
    if (c_nuc >= 0) {
      seqs->nucspec.nucs[c_nuc] = c_seq[i_nuc];
      seqs->dumspec.nucs[c_nuc] = c_seq[i_nuc];
      check(c_seq[i_nuc] >= 1 && c_seq[i_nuc] <= 4, 
        "Invalid nucleotide id: %i nucleotide code: %i", i_nuc, c_seq[i_nuc]);
    }
  }

  free(f_seq);
  free(f_con);
  free(c_seq);
  free(c_con);
  return ERR_OK;
error:
  free(f_seq);
  free(f_con);
  free(c_seq);
  free(c_con);
  return ERR_OOM;
}

int init_constraint_random(
      int * f_seq,
      int * c_seq,
      int * f_con,
      int * c_con,
      int n_nucs,
      options_t * opts
    ) {

  int n_allowed[15] = {
    4, 
    1, 1, 1, 1, 
    2, 2, 2, 2, 2, 2,
    3, 3, 3, 3
  };
  int allowed_seqs[15][4] = {
    {1, 2, 3, 4}, // N
    {1,-1,-1,-1}, // A
    {2,-1,-1,-1}, // C
    {3,-1,-1,-1}, // G
    {4,-1,-1,-1}, // U / T
    {1, 3,-1,-1}, // AG
    {1, 2,-1,-1}, // AC
    {2, 3,-1,-1}, // CG
    {1, 4,-1,-1}, // AU
    {3, 4,-1,-1}, // GU
    {2, 4,-1,-1}, // CU
    {1, 2, 3,-1}, // ACG
    {1, 2, 4,-1}, // ACU
    {1, 3, 4,-1}, // AGU
    {2, 3, 4,-1}, // CGU
  };

  int nuc_comps[5] = {
    0, 4, 3, 2, 1
  };
  int nuc_wobble[5] = {
    0, -1, -1, 4, 3
  };
  double def_weight[5] = {
    0.0, 1.0, 1.0, 1.0, 1.0
  };
  double weight[5] = {
    0.0, 1.0, 1.0, 1.0, 1.0
  };
  double r_num;
  double temp_num;
  double tot_weight;
  unsigned int i_nuc;
  int c;
  int n;
  int n_comp;
  int n_wobble;
  int i_con;
  int wc_allowed;
  int wobble_allowed;
  int n_choices;

  def_weight[BASE_A] = (double) (1.0 - opts->gc_init_prob);
  def_weight[BASE_U] = (double) (1.0 - opts->gc_init_prob);
  def_weight[BASE_G] = (double) opts->gc_init_prob;
  def_weight[BASE_C] = (double) opts->gc_init_prob;

  for ( i_nuc = 0; i_nuc < n_nucs; i_nuc++ ) {
    c = f_con[i_nuc];
    check(c >= 0 && c <= 15, "Invalid constraint at %u: %i", i_nuc, c);
    n_choices = n_allowed[c];
    for (i_con = 1; i_con < 5; i_con++) {
      weight[i_con] = 0;
    }
    for (i_con = 0; i_con < n_choices; i_con++) {
      weight[allowed_seqs[c][i_con]] = def_weight[allowed_seqs[c][i_con]];
    }
    tot_weight = 0;
    for (i_con = 1; i_con < 5; i_con++) {
      tot_weight += weight[i_con];
    }
    r_num = (double) tot_weight * genrand_real1();
    temp_num = 0;
    for (n = 1; n <= 4; n++) {
      temp_num += weight[n];
      if (temp_num > r_num) {
        break;
      }
    }
    check_debug(n >= 1 && n <= 4, "Invalid nucleotide selected: %i", n);
    f_seq[i_nuc] = n;
    // Set the WC complement and wobble complement
    n_comp = nuc_comps[n];
    n_wobble = nuc_wobble[n];
    
    // Get the constraints for the complement
    c = c_con[i_nuc];
    check(c >= 0 && c <= 14, "Invalid constraint for nucleotide %u of %i: %i",
       i_nuc, n_nucs, c);
    n_choices = n_allowed[c];
    check(n_choices > 0, "Invalid number of choices");

    wc_allowed = 0;
    wobble_allowed = 0;
    // Check if the wc or wobble pairs are allowed
    for (i_con = 0; i_con < n_choices; i_con++) {
      if (allowed_seqs[c][i_con] == n_comp) {
        wc_allowed = 1;
      }
      if (allowed_seqs[c][i_con] == n_wobble) {
        wobble_allowed = 1;
      }
    }
    if (wc_allowed) {
      c_seq[i_nuc] = n_comp;
    } else if (wobble_allowed && opts->allow_wobble) {
      c_seq[i_nuc] = n_wobble;
    } else if (opts->allow_mismatch) {
      for (i_con = 0; i_con < 5; i_con++) {
        weight[i_con] = 0;
      }
      for (i_con = 0; i_con < n_choices; i_con++) {
        weight[allowed_seqs[c][i_con]] = def_weight[allowed_seqs[c][i_con]];
      }
      tot_weight = 0;
      for (i_con = 0; i_con < 5; i_con++) {
        tot_weight += weight[i_con];
      }
      r_num = tot_weight * genrand_real1();

      temp_num = 0;
      for (n = 1; n < 4; n++) {
        temp_num += weight[n];
        if (temp_num > r_num) {
          break;
        }
      }
      c_seq[i_nuc] = n;
    } else {
      sentinel("Invalid constraint %u f_con: %i c_con: %i f_seq: %i", 
          i_nuc, f_con[i_nuc], c_con[i_nuc], f_seq[i_nuc]);
    }
  }
  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

/***********************************************************/
int get_struc_length(design_struc_t * struc, 
    sequence_spec_t * seq_table) {
  int n_strands = struc->n_strands;
  int i_strand;
  int n_nucs = 0;

  // int tok_len = 1; // This will eventually be variable, not yet though
  int c_strand;
  for (i_strand = 0; i_strand < n_strands; i_strand++) {
    c_strand = struc->strands[i_strand];
    n_nucs += seq_table->strands.specs[c_strand].n;
  }
  return n_nucs;
}
/***********************************************************/

struc_tree_t * alloc_struc_tree(int n_segs, int * seg_start,
    int * seg_stop, int n_assumed, int * assumed_i, int * assumed_j) {
  struc_tree_t * res = (struc_tree_t *) malloc(sizeof(struc_tree_t));
  int i;
  int * in_use = (int*) malloc(seg_stop[n_segs-1] * sizeof(int));
  int i_seg, s_start, s_stop, j_nuc;
  res->n_segments = n_segs;
  res->children = NULL;
  res->n_children = 0;
  res->seg_start = (int*) malloc(n_segs * sizeof(int));
  res->seg_stop = (int*) malloc(n_segs * sizeof(int));
  res->forbidden = NULL;
  for (i = 0; i < n_segs; i++) {
    res->seg_start[i] = seg_start[i];
    res->seg_stop[i] = seg_stop[i];
  }

  res->assumed_i = (int*) malloc(n_assumed * sizeof(int));
  res->assumed_j = (int*) malloc(n_assumed * sizeof(int));
  for (i = 0; i < n_assumed; i++) {
    res->assumed_i[i] = assumed_i[i];
    res->assumed_j[i] = assumed_j[i];
  }


  for (i_seg = 0; i_seg < n_segs; i_seg++) {
    s_start = seg_start[i_seg];
    s_stop = seg_stop[i_seg];

    for (j_nuc = s_start; j_nuc < s_stop; j_nuc++) {
      in_use[j_nuc] = 1;
    }
  }
  for (i = 0; i < n_assumed; i++) {
    // debug("assumed %i: %i %i max = %i", i, res->assumed_i[i], res->assumed_j[i], seg_stop[n_segs - 1]);
    check(res->assumed_i[i] < seg_stop[n_segs-1] 
        && in_use[res->assumed_i[i]], "%i marked as assumed, but not used", 
        res->assumed_i[i]);
    check(res->assumed_j[i] < seg_stop[n_segs-1] 
        && in_use[res->assumed_j[i]], "%i marked as assumed, but not used", 
        res->assumed_j[i]);
  }

  res->forbidden = alloc_split_tracker();

  res->n_assumed = n_assumed;

  free(in_use);
  return res;
error:
  free_struc_tree(res);
  free(in_use);
  return NULL;
}

/***********************************************************/
int init_struc_state(struc_state_t * state, 
    design_struc_t * struc) {
  int seg1, seg2;
  state->struc = struc;
  state->n = 1;
  seg1 = 0;
  seg2 = struc->n_nucs;
  state->tree = alloc_struc_tree(1, &seg1, &seg2, 0, NULL, NULL);
  return ERR_OK;
}
/***********************************************************/

/***********************************************************/
int free_struc_tree(struc_tree_t * t) {
  int i_c;
  if (t) {
    for (i_c = 0; i_c < t->n_children; i_c++) {
      free_struc_tree(t->children[i_c]);
    }
    free(t->seg_start);
    free(t->seg_stop);
    free(t->children);
    free(t->assumed_i);
    free(t->assumed_j);

    t->assumed_i = NULL;
    t->assumed_j = NULL;
    t->n_assumed = 0;

    t->seg_start = NULL;
    t->seg_stop = NULL;
    t->children = NULL;
    t->n_children = 0;
    destroy_split_tracker(t->forbidden);
    t->forbidden = NULL;
    free(t);
  }
  return ERR_OK;
}
/***********************************************************/

/***********************************************************/
static int free_struc_tree_children(struc_tree_t * t) {
  int i_c;
  for (i_c = 0; i_c < t->n_children; i_c++) {
    free_struc_tree(t->children[i_c]);
    t->children[i_c] = NULL;
  }
  free(t->children);
  t->children = NULL;
  t->n_children = 0;
  return ERR_OK;
}
/***********************************************************/

/***********************************************************/
int free_struc_state(struc_state_t * s) {
  free_struc_tree(s->tree);
  s->n = 0;

  return ERR_OK;
}
/***********************************************************/

/***********************************************************/
static int copy_struc_tree(struc_tree_t * dest, struc_tree_t * src) {
  int i_c;

  if (src->n_children != dest->n_children) {
    for (i_c = 0; i_c < dest->n_children; i_c++) {
      free_struc_tree(dest->children[i_c]);
    }

    free(dest->children);

    dest->children = (struc_tree_t **) malloc(sizeof(struc_tree_t*) * src->n_children);
    for (i_c = 0; i_c < src->n_children; i_c++) {
      dest->children[i_c] = alloc_struc_tree(src->children[i_c]->n_segments,
          src->children[i_c]->seg_start, src->children[i_c]->seg_stop,
          src->children[i_c]->n_assumed, src->children[i_c]->assumed_i,
          src->children[i_c]->assumed_j);
    }
  }

  if (dest->n_assumed != src->n_assumed) {
    free(dest->assumed_i);
    free(dest->assumed_j);
    dest->assumed_i = (int*) malloc(src->n_assumed);
    dest->assumed_j = (int*) malloc(src->n_assumed);
  }

  for (i_c = 0; i_c < src->n_assumed; i_c++) {
    dest->assumed_i[i_c] = src->assumed_i[i_c];
    dest->assumed_j[i_c] = src->assumed_j[i_c];
  }

  if (dest->n_segments != src->n_segments) {
    free(dest->seg_start);
    free(dest->seg_stop);
    dest->seg_start = (int*) malloc(src->n_segments);
    dest->seg_stop = (int*) malloc(src->n_segments);
  }

  for (i_c = 0; i_c < src->n_segments; i_c++) {
    dest->seg_start[i_c] = src->seg_start[i_c];
    dest->seg_stop[i_c] = src->seg_stop[i_c];
  }
  dest->n_assumed = src->n_assumed;
  dest->n_segments = src->n_segments;

  for (i_c = 0; i_c < src->n_children; i_c++) {
    copy_struc_tree(dest->children[i_c], src->children[i_c]);
  }
  dest->n_children = src->n_children;

  return ERR_OK;
}
/***********************************************************/

/***********************************************************/
int copy_struc_state(struc_state_t * dest, struc_state_t * src) {
  dest->struc = src->struc;
  dest->n = src->n;
  if (!dest->tree) {
    dest->tree = alloc_struc_tree(src->tree->n_segments,
        src->tree->seg_start, src->tree->seg_stop,
        src->tree->n_assumed, src->tree->assumed_i,
        src->tree->assumed_j);
  }
  copy_struc_tree(dest->tree, src->tree);
  return ERR_OK;
}
/***********************************************************/

static int mark_split_points(int * split_marker, int * struc,
    int n_nucs, int * breaks, int n_br, design_spec_t * spec) {
  int i_nuc, j_nuc;
  int i_h = 0;
  int h_spl = spec->opts.H_split;
  int i_br = 0;
  int l1, l2;
  int k_nuc, l_nuc;

  for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    split_marker[i_nuc] = -1;
  };

  for (i_nuc = 0; i_nuc < n_nucs - 1; i_nuc++) {
    if (struc[i_nuc] >= 0) {
      j_nuc = struc[i_nuc];
      if (j_nuc >= 0 && struc[i_nuc + 1] == j_nuc - 1) {
        if (i_h == 0) {
          i_h = 1;
        }
        i_h ++;
      } else {
        i_h = 0;
      }

      for (i_br = 0; i_br < n_br; i_br++) {
        if (breaks[i_br] == i_nuc + 1 || breaks[i_br] == j_nuc) {
          i_h = 0;
        }
      }

      l1 = (i_nuc + n_nucs - j_nuc - 1) ;
      l2 = (j_nuc - i_nuc + 1) ;

      if (i_nuc < j_nuc && i_h >= 2*h_spl 
          && l1 >= spec->opts.N_split && l2 >= spec->opts.N_split) {
        debug("i: %i j: %i i_h: %i l1: %i l2: %i", i_nuc, j_nuc, i_h, l1, l2);
        k_nuc = i_nuc - spec->opts.H_split + 2;
        l_nuc = struc[k_nuc];
        split_marker[k_nuc] = l_nuc;
        split_marker[l_nuc] = k_nuc;
      }
    } else {
      i_h = 0;
    }
  }
  return ERR_OK;
}

static int create_children(struc_tree_t * tree, design_struc_t * struc, 
    split_list_t * cur, design_spec_t * spec) {
  split_el_t * c_split = cur->head;
  int i_seg, j_seg, n_segs;
  int i_nuc, j_nuc; // start and end of segment in full indices
  int * seg_start = (int*) malloc(2 * tree->n_segments * sizeof(int));
  int * seg_stop = (int*) malloc(2 * tree->n_segments * sizeof(int));
  int * assumed_i = (int*) malloc((2 * tree->n_assumed + 1) * sizeof(int));
  int * assumed_j = (int*) malloc((2 * tree->n_assumed + 1) * sizeof(int));
  int i_assume;
  int j_assume;
  int i_c;
  int n_c;
  n_segs = tree->n_segments;
  (void)spec;
  (void)struc;

  if (tree->children) {
    for (i_c = 0; i_c < tree->n_children; i_c++) {
      free_struc_tree(tree->children[i_c]);
      tree->children[i_c] = NULL;
    }
  }
  free(tree->children);
  tree->children = NULL;
  tree->n_children = 0;

  n_c = 0;
  while (c_split) {
    n_c += 2;
    c_split = c_split->next;
  }
  i_c = 0;
  check(n_c > 0, "No children requested");

  tree->children = (struc_tree_t **) malloc(sizeof(struc_tree_t *) * n_c);

  c_split = cur->head;
  while (c_split) {
    // Make left child
    j_seg = 0;
    j_assume = 0;

    debug("Split %i: %i %i", i_c/2, c_split->lsplit, c_split->rsplit);
    for (i_seg = 0; i_seg < n_segs; i_seg++) {
      i_nuc = tree->seg_start[i_seg];
      j_nuc = tree->seg_stop[i_seg];
      if (i_nuc < c_split->lsplit) {
        seg_start[j_seg] = i_nuc;

        if (j_nuc < c_split->lsplit) {
          seg_stop[j_seg] = j_nuc;
        } else {
          seg_stop[j_seg] = c_split->lsplit;
        }
        j_seg ++;
      } 
      
      if (j_nuc > c_split->rsplit) {
        seg_stop[j_seg] = j_nuc;
        if (i_nuc < c_split->rsplit) {
          seg_start[j_seg] = c_split->rsplit + 1;
        } else {
          seg_start[j_seg] = i_nuc;
        }
        j_seg ++;
      }
    }

    for (i_assume = 0; i_assume < tree->n_assumed; i_assume++) {
      if (tree->assumed_i[i_assume] < c_split->lsplit || 
          tree->assumed_i[i_assume] > c_split->rsplit) {

        check(tree->assumed_j[i_assume] < c_split->lsplit ||
            tree->assumed_j[i_assume] > c_split->rsplit,
            "Invalid assumption, split to two nodes %i %i | %i %i",
            tree->assumed_i[i_assume], tree->assumed_j[i_assume],
            c_split->lsplit, c_split->rsplit);
        assumed_i[j_assume] = tree->assumed_i[i_assume];
        assumed_j[j_assume] = tree->assumed_j[i_assume];
        j_assume ++;
      }
    }

    assumed_i[j_assume] = c_split->lsplit - 1;
    assumed_j[j_assume] = c_split->rsplit + 1;
    j_assume ++;

    tree->children[i_c] = alloc_struc_tree(j_seg, seg_start, seg_stop,
        j_assume, assumed_i, assumed_j);
    i_c ++;

    j_seg = 0;
    for (i_seg = 0; i_seg < n_segs; i_seg++) {
      i_nuc = tree->seg_start[i_seg];
      j_nuc = tree->seg_stop[i_seg];
      if (j_nuc > c_split->lsplit && i_nuc < c_split->rsplit) {
        if (i_nuc > c_split->lsplit) {
          seg_start[j_seg] = i_nuc;
        } else {
          seg_start[j_seg] = c_split->lsplit;
        }
        if (j_nuc < c_split->rsplit) {
          seg_stop[j_seg] = j_nuc;
        } else {
          seg_stop[j_seg] = c_split->rsplit + 1;
        }
        j_seg ++;
      }
    }

    j_assume = 0;
    for (i_assume = 0; i_assume < tree->n_assumed; i_assume++) {
      if (tree->assumed_i[i_assume] > c_split->lsplit && 
          tree->assumed_i[i_assume] < c_split->rsplit) {
  
        check(tree->assumed_j[i_assume] > c_split->lsplit &&
            tree->assumed_j[i_assume] < c_split->rsplit,
            "Invalid assumption, split to two nodes %i %i | %i %i",
            tree->assumed_i[i_assume], tree->assumed_j[i_assume],
            c_split->lsplit, c_split->rsplit);

        assumed_i[j_assume] = tree->assumed_i[i_assume];
        assumed_j[j_assume] = tree->assumed_j[i_assume];
        j_assume ++;
      }
    }
    

    assumed_i[j_assume] = c_split->rsplit;
    assumed_j[j_assume] = c_split->lsplit;
    j_assume ++;

    tree->children[i_c] = alloc_struc_tree(j_seg, seg_start, seg_stop,
        j_assume, assumed_i, assumed_j);
    i_c ++;
    tree->n_children = i_c;

    c_split = c_split->next;
  }
  free(assumed_i);
  free(assumed_j);
  free(seg_start);
  free(seg_stop);

  tree->n_children = i_c;
  return ERR_OK;
error:
  free_struc_tree_children(tree);
  free(assumed_i);
  free(assumed_j);
  free(seg_start);
  free(seg_stop);
  return ERR_INVALID_STATE;
}

static int get_decomp_strucs(int ** struc_list, 
    struc_tree_t * tree, design_struc_t * struc, 
    design_spec_t * spec) {
  design_struc_t * cur_struc = struc;
  int i_struc = 0;
  (void)spec;

  while (cur_struc) {
    if (tree_is_decomp_for(tree, cur_struc)) {
      struc_list[i_struc] = cur_struc->struc;
      i_struc ++;
    }
    cur_struc = cur_struc->next;
  }

  return i_struc;
}

static int get_native_nuc_maps(int * to_full, int * to_node, 
    int * n_node_nucs, int * breaks, int * n_breaks, 
    struc_tree_t * tree, design_struc_t * struc,
    design_spec_t * spec) {

  int i_nuc, j_nuc, n_nucs, i_seg, n_segs, i_br;
  int i_str, c_str, n_strs, n_str_nucs;
  n_segs = tree->n_segments;
  n_nucs = struc->n_nucs;

  for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    to_full[i_nuc] = -1;
    to_node[i_nuc] = -1;
  }

  j_nuc = 0; // cur nuc within node

  n_strs = struc->n_strands; // n of strands
  i_str = 0; // cur strand index
  c_str = struc->strands[i_str]; // cur strand id

  n_str_nucs = spec->seqs.strands.specs[c_str].n;
  // index of the end of the current strand

  i_br = 0;

  for (i_seg = 0; i_seg < n_segs; i_seg++) {
    debug("Segment: %i  start: %i stop: %i", i_seg, tree->seg_start[i_seg], tree->seg_stop[i_seg]);
    for (i_nuc = tree->seg_start[i_seg]; 
        i_nuc < tree->seg_stop[i_seg]; i_nuc++) {
      // need to increment strand
      while (i_nuc >= n_str_nucs) {
        if (j_nuc > 0) {
          // if inside of node
          if (i_br == 0 || breaks[i_br-1] != j_nuc) {
            breaks[i_br] = j_nuc;
            i_br ++;
          }
        }
        i_str ++;

        check(i_str < n_strs, 
            "N strands exceeded in get_native_nuc_maps");
        c_str = struc->strands[i_str];
        n_str_nucs += spec->seqs.strands.specs[c_str].n;
      }
      to_full[j_nuc] = i_nuc;
      to_node[i_nuc] = j_nuc;
      j_nuc ++;
    }
    if (i_seg != n_segs - 1) {
      breaks[i_br] = j_nuc;
      i_br ++;
    }
  }
  *n_node_nucs = j_nuc;
  *n_breaks = i_br;

  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

int get_strand_breaks(int * breaks, design_struc_t * struc,
    design_spec_t * spec) {
  int i_str;
  int i_nuc = 0;
  int c_str;
  for (i_str = 0; i_str < struc->n_strands - 1; i_str++) {
    c_str = struc->strands[i_str];
    check(c_str < spec->seqs.strands.n, 
        "Invalid strand %i encountered", c_str);
    i_nuc += spec->seqs.strands.specs[c_str].n;
    breaks[i_str] = i_nuc;
  }
  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

int tree_is_decomp_for(struc_tree_t * tree, design_struc_t * struc) {
  int i_assume;
  int i_nuc, j_nuc;
  int rval = 1;
  if (struc->modifiable) {
    rval = 0;
  }
  for (i_assume = 0; i_assume < tree->n_assumed; i_assume++) {
    i_nuc = tree->assumed_i[i_assume];
    j_nuc = tree->assumed_j[i_assume];
    if (struc->struc[i_nuc] != j_nuc) {
      rval = 0;
      break;
    }
  }
  return rval;
}


static int init_result_tree(result_tree_t * res, struc_tree_t * tree,
    design_struc_t * struc, design_spec_t * spec) {
  int i_c;
  int i_seg, n_segs;
  int i_nuc, j_nuc, s_start, s_stop, n_nucs;
  res->sequence = NULL;
  res->native_map = NULL;
  res->nuc_defects = NULL;
  res->n_nucs = 0;

  res->eval_time = 0.0;
  res->pfunc = 0.0;

  res->ppairs = NULL;
  res->ppairs_i = NULL;
  res->ppairs_j = NULL;
  res->ppairs_n = 0;
  res->ppairs_cap = 0;

  res->n_children = tree->n_children;
  res->children = NULL;
  if (res->n_children > 0) {
    res->children = (result_tree_t **) 
      malloc(res->n_children * sizeof(result_tree_t *));


    n_segs = tree->n_segments;
    n_nucs = 0;
    for (i_seg = 0; i_seg < n_segs; i_seg++) {
      if (tree->seg_start[i_seg] != 0 && spec->opts.include_dummies) {
        n_nucs += spec->opts.H_split;
      }
      n_nucs += tree->seg_stop[i_seg] - tree->seg_start[i_seg];
      if (tree->seg_stop[i_seg] != struc->n_nucs && spec->opts.include_dummies) {
        n_nucs += spec->opts.H_split;
      }
    }
    res->n_nucs = n_nucs;

    res->sequence = (int*) malloc(n_nucs * sizeof(int));
    res->native_map = (int*) malloc(n_nucs * sizeof(int));
    res->dummy_flag = (int*) malloc(n_nucs * sizeof(int));
    res->nuc_defects = (DBL_TYPE *) malloc(n_nucs * sizeof(DBL_TYPE));

    check_mem(res->sequence);
    check_mem(res->native_map);
    check_mem(res->dummy_flag);
    check_mem(res->nuc_defects);

    i_nuc = 0;
    for (i_seg = 0; i_seg < n_segs; i_seg++) {
      s_start = tree->seg_start[i_seg];
      if (tree->seg_start[i_seg] != 0 && spec->opts.include_dummies) {
        s_start -= spec->opts.H_split;
      }
      s_stop = tree->seg_stop[i_seg];
      if (tree->seg_stop[i_seg] != struc->n_nucs && spec->opts.include_dummies) {
        s_stop += spec->opts.H_split;
      }
      for (j_nuc = s_start; j_nuc < s_stop; j_nuc++) {
        res->native_map[i_nuc] = j_nuc;
        i_nuc ++;
      }
    }

    for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
      res->sequence[i_nuc] = 0;
    }

    i_nuc = 0;
    for (i_seg = 0; i_seg < n_segs; i_seg++) {
      s_start = tree->seg_start[i_seg];
      if (tree->seg_start[i_seg] != 0 && spec->opts.include_dummies) {
        s_start -= spec->opts.H_split;
      }
      s_stop = tree->seg_stop[i_seg];
      if (tree->seg_stop[i_seg] != struc->n_nucs && spec->opts.include_dummies) {
        s_stop += spec->opts.H_split;
      }
      for (j_nuc = s_start; j_nuc < s_stop; j_nuc++) {
        res->native_map[i_nuc] = j_nuc;
        if (j_nuc < tree->seg_start[i_seg] || j_nuc >= tree->seg_stop[i_seg]) {
          res->dummy_flag[i_nuc] = 1;
        } else {
          res->dummy_flag[i_nuc] = 0;
        }
        res->nuc_defects[i_nuc] = 0;
        i_nuc ++;
      }
    }

    for (i_c = 0; i_c < res->n_children; i_c++) {
      res->children[i_c] = alloc_result_tree(tree->children[i_c], struc, spec);
    }
  } else {
    n_segs = tree->n_segments;
    n_nucs = 0;
    for (i_seg = 0; i_seg < n_segs; i_seg++) {
      if (tree->seg_start[i_seg] != 0 && spec->opts.include_dummies) {
        n_nucs += spec->opts.H_split;
      }
      n_nucs += tree->seg_stop[i_seg] - tree->seg_start[i_seg];
      if (tree->seg_stop[i_seg] != struc->n_nucs && spec->opts.include_dummies) {
        n_nucs += spec->opts.H_split;
      }
    }

    check(n_nucs > 0, "Invalid result init. length = %i", n_nucs);

    res->sequence = (int*) malloc(n_nucs * sizeof(int));
    res->native_map = (int*) malloc(n_nucs * sizeof(int));
    res->dummy_flag = (int*) malloc(n_nucs * sizeof(int));
    res->nuc_defects = (DBL_TYPE *) malloc(n_nucs * sizeof(DBL_TYPE));
    check_mem(res->sequence);
    check_mem(res->native_map);
    check_mem(res->dummy_flag);
    check_mem(res->nuc_defects);

    for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
      res->sequence[i_nuc] = 0;
    }

    i_nuc = 0;
    for (i_seg = 0; i_seg < n_segs; i_seg++) {
      s_start = tree->seg_start[i_seg];
      s_stop = tree->seg_stop[i_seg];
      if (spec->opts.include_dummies) {
        if (tree->seg_start[i_seg] != 0 && spec->opts.include_dummies) {
          s_start -= spec->opts.H_split;
        }
        if (tree->seg_stop[i_seg] != struc->n_nucs && spec->opts.include_dummies) {
          s_stop += spec->opts.H_split;
        }
      }
      for (j_nuc = s_start; j_nuc < s_stop; j_nuc++) {
        res->native_map[i_nuc] = j_nuc;
        if (j_nuc < tree->seg_start[i_seg] || j_nuc >= tree->seg_stop[i_seg]) {
          res->dummy_flag[i_nuc] = 1;
        } else {
          res->dummy_flag[i_nuc] = 0;
        }
        res->nuc_defects[i_nuc] = 0;
        i_nuc ++;
      }
    }
    check(n_nucs == i_nuc, "Invalid nucleotide counts %i != %i", n_nucs, i_nuc);
    res->n_nucs = n_nucs;
  }
  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

static int init_empty_result_tree(result_tree_t * res) {
  res->sequence = NULL;
  res->native_map = NULL;
  res->dummy_flag = NULL;
  res->nuc_defects = NULL;
  res->n_nucs = 0;

  res->pfunc = 0;
  res->eval_time = 0;
  res->n_children = 0;
  res->children = NULL;

  res->ppairs = NULL;
  res->ppairs_i = NULL;
  res->ppairs_j = NULL;
  res->ppairs_n = 0;
  res->ppairs_cap = 0;
  return ERR_OK;
}

result_tree_t * alloc_empty_result_tree() {
  result_tree_t * res = (result_tree_t*) malloc(sizeof(result_tree_t));
  init_empty_result_tree(res);
  return res;
}

static int print_child_counts_helper(FILE * f, result_tree_t * tree) {
  int i_c;
  fprintf(f, ", %i", tree->n_children);
  for (i_c = 0; i_c < tree->n_children; i_c++) {
    print_child_counts_helper(f, tree->children[i_c]);
  }
  return ERR_OK;

}

int print_child_counts(FILE * f, result_tree_t * tree) {
  int i_c;
  fprintf(f, "[%i", tree->n_children);
  for (i_c = 0; i_c < tree->n_children; i_c++) {
    print_child_counts_helper(f, tree->children[i_c]);
  }
  fprintf(f, "]");
  return ERR_OK;
}


result_tree_t * alloc_result_tree(struc_tree_t * tree, 
    design_struc_t * struc, design_spec_t * spec) {

  result_tree_t * res = (result_tree_t*) malloc(sizeof(result_tree_t));
  init_result_tree(res, tree, struc, spec);

  return res;
}

int copy_result_tree(result_tree_t * dest, result_tree_t * src) {
  int i_nuc;
  int n_nucs;
  int i_ch;
  int * temp_int;
  DBL_TYPE * temp_dbl;
  n_nucs = src->n_nucs;

  if (n_nucs > dest->n_nucs) {
    temp_int = (int*) realloc(dest->sequence, n_nucs * sizeof(int));
    check_mem(temp_int);
    dest->sequence = temp_int;
    temp_int = (int*) realloc(dest->native_map, n_nucs * sizeof(int));
    check_mem(temp_int);
    dest->native_map = temp_int;
    temp_int = (int*) realloc(dest->dummy_flag, n_nucs * sizeof(int));
    check_mem(temp_int);
    dest->dummy_flag = temp_int;
    temp_dbl = (DBL_TYPE*) realloc(dest->nuc_defects, n_nucs * sizeof(DBL_TYPE));
    check_mem(temp_dbl);
    dest->nuc_defects = temp_dbl;
  }

  for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    dest->sequence[i_nuc] = src->sequence[i_nuc];
    dest->native_map[i_nuc] = src->native_map[i_nuc];
    dest->dummy_flag[i_nuc] = src->dummy_flag[i_nuc];
    dest->nuc_defects[i_nuc] = src->nuc_defects[i_nuc];
  }
  dest->n_nucs = src->n_nucs;

  dest->pfunc = src->pfunc;
  dest->eval_time = src->eval_time;

  if (src->n_children != dest->n_children) {
    for (i_ch = 0; i_ch < dest->n_children; i_ch++) {
      destroy_result_tree(dest->children[i_ch]);
    }
    free(dest->children);
    dest->children = (result_tree_t **) 
      malloc(sizeof(result_tree_t *) * src->n_children);

    for (i_ch = 0; i_ch < src->n_children; i_ch++) {
      dest->children[i_ch] = alloc_empty_result_tree();
    }
  }
  dest->n_children = src->n_children;

  for (i_ch = 0; i_ch < src->n_children; i_ch++) {
    copy_result_tree(dest->children[i_ch], src->children[i_ch]);
  }

  if (src->ppairs_n > dest->ppairs_cap) {
    temp_dbl = (DBL_TYPE*) realloc(dest->ppairs, src->ppairs_n * sizeof(DBL_TYPE));
    check_mem(temp_dbl);
    dest->ppairs = temp_dbl;
    temp_int = (int*) realloc(dest->ppairs_i, src->ppairs_n * sizeof(int));
    check_mem(temp_int);
    dest->ppairs_i = temp_int;
    temp_int = (int*) realloc(dest->ppairs_j, src->ppairs_n * sizeof(int));
    check_mem(temp_int);
    dest->ppairs_j = temp_int;

    dest->ppairs_cap = src->ppairs_n;
  }

  dest->ppairs_n = src->ppairs_n;
  for (i_nuc = 0; i_nuc < src->ppairs_n; i_nuc++) {
    dest->ppairs[i_nuc] = src->ppairs[i_nuc];
    dest->ppairs_i[i_nuc] = src->ppairs_i[i_nuc];
    dest->ppairs_j[i_nuc] = src->ppairs_j[i_nuc];
  }
  
  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

static int get_eval_sequence(int * evalseq,
    result_tree_t * res, result_struc_t * struc,
    struc_tree_t * tree, design_spec_t * spec) {

  int n_brs = struc->n_breaks;
  int i_br;
  int n_segs = tree->n_segments;
  int i_seg;
  int i_nuc, j_nuc, seg_start, seg_stop;


  i_nuc = 0;
  i_br = 0;
  for (i_seg = 0; i_seg < n_segs; i_seg++ ) {
    seg_start = tree->seg_start[i_seg];
    seg_stop = tree->seg_stop[i_seg];

    if (seg_start > 0 && spec->opts.include_dummies) {
      seg_start -= spec->opts.H_split;
    }
    if (seg_stop < struc->n_nucs && spec->opts.include_dummies) {
      seg_stop += spec->opts.H_split;
    }
    while (i_br < n_brs && struc->breaks[i_br] < seg_start) {
      i_br ++;
    }
    for (j_nuc = seg_start; j_nuc < seg_stop; j_nuc++) {
      if (i_br < n_brs && struc->breaks[i_br] == j_nuc) {
        if (i_nuc > 0 && evalseq[i_nuc-1] != STRAND_PLUS) {
          evalseq[i_nuc] = STRAND_PLUS;
          i_nuc ++;
        }
        i_br ++;
      } 
      check_debug(j_nuc >= 0 && j_nuc < struc->n_nucs,
          "j_nuc is not a valid: %i length: %i", j_nuc, struc->n_nucs);
      if (spec->opts.fake_dummies && 
          (j_nuc < tree->seg_start[i_seg] || j_nuc >= tree->seg_stop[i_seg])) {
        evalseq[i_nuc] = struc->f_sequence[j_nuc];
      } else {
        evalseq[i_nuc] = struc->sequence[j_nuc];
      }
      i_nuc ++;
    }
    evalseq[i_nuc] = STRAND_PLUS;
    i_nuc ++;
  }

  evalseq[i_nuc-1] = -1;

  check(i_nuc <= res->n_nucs + tree->n_segments + struc->n_breaks + 2,
      "number of elements in eval exceeds assumption: %i %i",
      i_nuc, res->n_nucs + tree->n_segments + struc->n_breaks + 2);

  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

static DBL_TYPE evaluate_dummy_pfuncs(result_struc_t * struc,
    struc_tree_t * tree, design_spec_t * spec) {

  int i_assume, i_nuc, j_nuc;
  int nuc_start, nuc_end;
  int n_nucs = 2 * (spec->opts.H_split + 1);

  // One for break, one for termination
  int * curseq = (int*) malloc((n_nucs + 2)* sizeof(int));
  DBL_TYPE res = 1;

  // TODO make the dummy sequences, should be from assumed_i to assumed_i + H_split then
  // strand break, then assumed_j - H_split + 1 through assumed_j evaluate the partition function
  // and take QB of the upper right hand corner
  EXTERN_QB = (DBL_TYPE*) malloc(n_nucs * n_nucs * sizeof(DBL_TYPE));

  check_mem(curseq);
  check_mem(EXTERN_QB);

  if (spec->opts.include_dummies) {
    for (i_assume = 0; i_assume < tree->n_assumed; i_assume++) {
      j_nuc = 0;
      nuc_start = tree->assumed_i[i_assume];
      nuc_end = nuc_start + spec->opts.H_split + 1;
      curseq[j_nuc] = struc->sequence[nuc_start];
      j_nuc ++;

      for (i_nuc = nuc_start + 1; i_nuc < nuc_end; i_nuc++) {
        if (spec->opts.fake_dummies) {
          curseq[j_nuc] = struc->f_sequence[i_nuc];
        } else {
          curseq[j_nuc] = struc->sequence[i_nuc];
        }
        j_nuc ++;
      }

      curseq[j_nuc] = STRAND_PLUS;
      j_nuc ++;

      nuc_end = tree->assumed_j[i_assume] + 1;
      nuc_start = nuc_end - spec->opts.H_split - 1;

      for (i_nuc = nuc_start; i_nuc < nuc_end - 1; i_nuc++) {
        if (spec->opts.fake_dummies) {
          curseq[j_nuc] = struc->f_sequence[i_nuc];
        } else {
          curseq[j_nuc] = struc->sequence[i_nuc];
        }
        j_nuc ++;
      }
      curseq[j_nuc] = struc->sequence[nuc_end-1];
      j_nuc ++;
      curseq[j_nuc] = -1;
      j_nuc ++;

      check(j_nuc == (n_nucs + 2), "Invalid sequence length");
      // char * tmpseq = (char*) malloc((n_nucs + 2) * sizeof(char));
      // convert_nucs_to_str(tmpseq, curseq, n_nucs + 1, &(spec->opts));
      // debug("Dummy: %s", tmpseq);
      // free(tmpseq);

      pfuncFull(curseq, 3, spec->opts.material, spec->opts.dangle_type, 
          spec->opts.temperature - ZERO_C_IN_KELVIN, 1, 
          spec->opts.sodium, spec->opts.magnesium,
          spec->opts.use_long_helix);

      // debug("dummy: %Lf", -kB * spec->opts.temperature * LOG_FUNC(EXTERN_QB[pf_index(0, 2*(spec->opts.H_split + 1) - 1, 
      //     2*(spec->opts.H_split + 1))]));
      res *= EXTERN_QB[pf_index(0, 2*(spec->opts.H_split + 1) - 1, 
          2*(spec->opts.H_split + 1))];
    }
  } else {
    // for (i_assume = 0; i_assume < tree->n_assumed; i_assume++) {
    //   i_nuc = tree->assumed_i[i_assume];
    //   j_nuc = tree->assumed_j[i_assume];
    //   if (struc->sequence[i_nuc] != BASE_C && 
    //       struc->sequence[j_nuc] != BASE_C) {
    //     res *= EXP_FUNC(-AT_PENALTY / (kB * spec->opts.temperature));
    //   }
    // }
  }
  free(curseq);
  curseq = NULL;
  free(EXTERN_QB);
  EXTERN_QB = NULL;

  return res;
error:
  free(curseq);
  free(EXTERN_QB);
  EXTERN_QB = NULL;
  return 1;
}

int eval_leaf(result_tree_t * res, result_struc_t * struc,
    struc_tree_t * tree, design_spec_t * spec) {

  int i_nuc, j_nuc, k_nuc, m_i_nuc, m_j_nuc, n_nucs;
  int i_assume;
  int * tempseq = (int*) malloc((res->n_nucs + tree->n_segments 
        + struc->n_breaks + 2) * sizeof(int));
  int * tmp_int = NULL;
  DBL_TYPE * tmp_dbl = NULL;
  int changed = 0;
  DBL_TYPE pfunc_corrected, cur_defect, cur_prob;
  DBL_TYPE min_ppair;
  DBL_TYPE bonus_per_split = spec->opts.bonus_per_split;
  DBL_TYPE total_bonus = 1;
  int n_bonuses = 0;
  int * rev_map = NULL;
  struct timeval start_time;
  struct timeval end_time;


  // Get updated sequence
  n_nucs = res->n_nucs;
  for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    j_nuc = res->native_map[i_nuc];
    check_debug(j_nuc < struc->n_nucs, "Invalid nucleotide map i: %i j: %i dum: %i ",
        i_nuc, j_nuc, res->dummy_flag[i_nuc]);
    if (res->dummy_flag[i_nuc] && spec->opts.fake_dummies) {
      tempseq[i_nuc] = struc->f_sequence[j_nuc];
    } else {
      tempseq[i_nuc] = struc->sequence[j_nuc];
    }
  }


  // Check if sequence has changed
  for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    if (tempseq[i_nuc] != res->sequence[i_nuc]) {
      changed = 1;
      break;
    }
  }

// #define NUPACK_TTDESIGN_NOCACHING
#ifdef NUPACK_TTDESIGN_NOCACHING
  changed = 1;
#endif
  DBL_TYPE * bonuses = NULL; 
  int i_seg = 0;

  // If sequence has changed
  if (changed) {
    gettimeofday(&start_time, NULL);
    rev_map = (int*) malloc((struc->n_nucs * sizeof(int)));

    for (i_nuc = 0; i_nuc < struc->n_nucs; i_nuc++) {
      rev_map[i_nuc] = -1;
    }

    for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
      rev_map[res->native_map[i_nuc]] = i_nuc;
    }

    // Fill out sequence with breaks included
    for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
      res->sequence[i_nuc] = tempseq[i_nuc];
    }
   
    check(ERR_OK == get_eval_sequence(tempseq, res, struc, tree, spec),
        "Error constructing sequence");

    char * printstring = (char *) malloc((n_nucs*2 + 1) * sizeof(char));

    for (i_nuc = 0; i_nuc < 2*n_nucs; i_nuc++) {
      if (tempseq[i_nuc] == -1) {
        break;
      }
    }
    convert_nucs_to_str(printstring, tempseq, i_nuc, &(spec->opts));

    debug("Sequence: %s", printstring);
    free(printstring);

    // evaluate pfunc and pair probs
    EXTERN_Q = (DBL_TYPE*) calloc(res->n_nucs * res->n_nucs, sizeof(DBL_TYPE));
    pairPr = (DBL_TYPE*) calloc(res->n_nucs * (res->n_nucs + 1), sizeof(DBL_TYPE));

    check_mem(EXTERN_Q);
    check_mem(pairPr);

    TEMP_K = spec->opts.temperature;
    SODIUM_CONC = spec->opts.sodium;
    MAGNESIUM_CONC = spec->opts.magnesium;
    USE_LONG_HELIX_FOR_SALT_CORRECTION = spec->opts.use_long_helix;
    DNARNACOUNT = spec->opts.material;
    DBL_TYPE Q_val = 1;

    if (! spec->opts.include_dummies) {
 
      bonuses = (DBL_TYPE*) calloc(struc->n_nucs * struc->n_nucs, sizeof(DBL_TYPE));
      i_nuc = 0;
      for (i_nuc = 0; i_nuc < struc->n_nucs * struc->n_nucs; i_nuc++) {
        bonuses[i_nuc] = 1;
      }
      i_nuc = 0;
      for (i_seg = 0; i_seg < tree->n_segments - 1; i_seg++) {
        i_nuc += tree->seg_stop[i_seg] - tree->seg_start[i_seg];
        total_bonus *= EXP_FUNC(-bonus_per_split / (kB * spec->opts.temperature));
        bonuses[(i_nuc - 1) * n_nucs + i_nuc] = EXP_FUNC(-bonus_per_split / (kB * spec->opts.temperature));
        n_bonuses ++;
      }


      if (tree->seg_start[0] != 0) {
        total_bonus *= EXP_FUNC(-bonus_per_split / (kB * spec->opts.temperature));
        bonuses[0 * n_nucs + n_nucs - 1] = EXP_FUNC(-bonus_per_split / (kB * spec->opts.temperature));
        n_bonuses ++;
      }

      Q_val = pfuncFullWithBonuses(tempseq, 3, spec->opts.material, spec->opts.dangle_type,
          spec->opts.temperature - ZERO_C_IN_KELVIN, 1, 1,
          spec->opts.sodium, spec->opts.magnesium,
          spec->opts.use_long_helix, bonuses);

      Q_val /= total_bonus;

      free(bonuses);
      bonuses = NULL;
    } else {
      pfuncFull(tempseq, 3, spec->opts.material, spec->opts.dangle_type, 
          spec->opts.temperature - ZERO_C_IN_KELVIN, 1, 
          spec->opts.sodium, spec->opts.magnesium,
          spec->opts.use_long_helix);
    }


    // save the nucleotide defects (only for the top structure)
    for (i_nuc = 0; i_nuc < res->n_nucs; i_nuc++) {
      m_i_nuc = res->native_map[i_nuc];
      m_j_nuc = struc->structure[m_i_nuc];
      if (m_j_nuc < 0) {
        j_nuc = res->n_nucs;
      } else if (rev_map[m_j_nuc] == -1) {
        j_nuc = i_nuc;
      } else {
        j_nuc = rev_map[m_j_nuc];
      }

      if (i_nuc < j_nuc) {
        k_nuc = i_nuc;
      } else {
        k_nuc = j_nuc;
        j_nuc = i_nuc;
      }

      check(k_nuc >= 0 && k_nuc < res->n_nucs,
          "Invalid nucleotide after map: %i", k_nuc);
      check(j_nuc >= 0 && j_nuc <= res->n_nucs,
          "Invalid nucleotide after reverse map: %i", j_nuc); 

      cur_defect = 1.0 - pairPr[k_nuc * (res->n_nucs + 1) + j_nuc];

      res->nuc_defects[i_nuc] = cur_defect;
    }

    res->ppairs_n = 0;

    // Save only native ppairs with prob > min_prob 
    for (i_nuc = 0; i_nuc < res->n_nucs; i_nuc++) {
      if (!res->dummy_flag[i_nuc]) {
        if (res->ppairs_n + res->n_nucs > res->ppairs_cap) {
          res->ppairs_cap = (res->ppairs_n + res->n_nucs) * 2;
          tmp_dbl = (DBL_TYPE*) realloc(res->ppairs, 
              res->ppairs_cap * sizeof(DBL_TYPE));
          check_mem(tmp_dbl);
          res->ppairs = tmp_dbl;

          tmp_int = (int*) realloc(res->ppairs_i,
              res->ppairs_cap * sizeof(int));
          check_mem(tmp_int);
          res->ppairs_i = tmp_int;

          tmp_int = (int*) realloc(res->ppairs_j,
              res->ppairs_cap * sizeof(int));
          check_mem(tmp_int);
          res->ppairs_j = tmp_int;
        }
        m_i_nuc = res->native_map[i_nuc];
        for (j_nuc = i_nuc; j_nuc < res->n_nucs; j_nuc++) {
          if (!res->dummy_flag[j_nuc]) {
            m_j_nuc = res->native_map[j_nuc];
            cur_prob = pairPr[i_nuc * (res->n_nucs + 1) + j_nuc];
            if (cur_prob > spec->opts.min_ppair_saved) { 
              res->ppairs_i[res->ppairs_n] = m_i_nuc;
              res->ppairs_j[res->ppairs_n] = m_j_nuc;
              res->ppairs[res->ppairs_n] = cur_prob;
              res->ppairs_n ++;
            }
          }
        }
        cur_prob = pairPr[i_nuc * (res->n_nucs + 1) + res->n_nucs];
        if (cur_prob > spec->opts.min_ppair_saved) {
          res->ppairs_i[res->ppairs_n] = m_i_nuc;
          res->ppairs_j[res->ppairs_n] = -1;
          res->ppairs[res->ppairs_n] = cur_prob;
          res->ppairs_n ++;
        }
      }
    }

    // correct the partition function
    pfunc_corrected = EXTERN_Q[pf_index(0, res->n_nucs - 1, res->n_nucs)];

    if (spec->opts.include_dummies) {
      min_ppair = 1.0;
      for (i_assume = 0; i_assume < tree->n_assumed; i_assume++) {
        m_i_nuc = tree->assumed_i[i_assume];
        m_j_nuc = tree->assumed_j[i_assume];

        if (m_i_nuc > m_j_nuc) {
          m_i_nuc = m_j_nuc;
          m_j_nuc = tree->assumed_i[i_assume];
        }

        i_nuc = rev_map[m_i_nuc];
        j_nuc = rev_map[m_j_nuc];

        check(i_nuc >= 0 && j_nuc >= 0, "Invalid assumed nucleotides %i %i %i %i %i",
            i_nuc, j_nuc, m_i_nuc, m_j_nuc, i_assume);
        check(!res->dummy_flag[i_nuc], "Dummy flag set for assumed nuc %i", i_nuc);
        check(!res->dummy_flag[j_nuc], "Dummy flag set for assumed nuc %i", j_nuc);

        // if (pairPr[i_nuc * (res->n_nucs + 1) + j_nuc] < 0.99) {
        //   debug("Nuc pair %i %i -> %i %i pair %Lf", m_i_nuc, m_j_nuc, i_nuc, j_nuc, 
        //       pairPr[i_nuc * (res->n_nucs + 1) + j_nuc]);
        // }

        if (pairPr[i_nuc * (res->n_nucs + 1) + j_nuc] < min_ppair) {
          min_ppair = pairPr[i_nuc * (res->n_nucs + 1) + j_nuc];
        }
      }
      min_ppair = 1.0;

      if (min_ppair < 0.000001) {
        min_ppair =   0.000001;
        debug("Min ppair ~= 0 !!, expanding");
      }
      pfunc_corrected *= min_ppair;
    } else {
      pfunc_corrected /= total_bonus;
    }


    pfunc_corrected /= evaluate_dummy_pfuncs(struc, tree, spec);

    free(rev_map);
    free(EXTERN_Q);
    free(pairPr);
    EXTERN_Q = NULL;
    pairPr = NULL;


    res->pfunc = pfunc_corrected;
    gettimeofday(&end_time, NULL);
    res->eval_time = end_time.tv_sec - start_time.tv_sec 
      + 1.0e-6 * (end_time.tv_usec - start_time.tv_usec);
  }


  free(tempseq);
  return ERR_OK;
error:
  free(tempseq);
  free(rev_map);
  return ERR_INVALID_STATE;
}

int destroy_result_tree(result_tree_t * res) {
  free_result_tree(res);
  free(res);
  return ERR_OK;
}

int free_result_tree(result_tree_t * res) {
  int i_c;
  for (i_c = 0; i_c < res->n_children; i_c++) {
    destroy_result_tree(res->children[i_c]);
  }
  free(res->sequence);
  free(res->native_map);
  free(res->dummy_flag);
  free(res->nuc_defects);
  free(res->children);
  free(res->ppairs);
  free(res->ppairs_i);
  free(res->ppairs_j);

  res->sequence = NULL;
  res->native_map = NULL;
  res->dummy_flag = NULL;
  res->nuc_defects = NULL;
  res->children = NULL;
  res->ppairs = NULL;
  res->ppairs_i = NULL;
  res->ppairs_j = NULL;

  res->n_nucs = 0;
  res->n_children = 0;
  res->pfunc = 0;
  res->eval_time = 0;
  res->ppairs_n = 0;
  res->ppairs_cap = 0;
  return ERR_OK;
}

static int merge_tree_pfuncs(result_tree_t * restree, result_struc_t * res, 
    design_spec_t * spec) {
  int i_c;
  DBL_TYPE pfunc = 0;
  DBL_TYPE eval_time = 0;
  if (restree->n_children > 0) {
    for (i_c = 0; i_c < restree->n_children; i_c++) {
      merge_tree_pfuncs(restree->children[i_c], res, spec);
    }
    for (i_c = 0; i_c < restree->n_children / 2; i_c++) {
      pfunc += compute_pfunc_merge(restree->children[i_c * 2],
          restree->children[i_c * 2 + 1], res, spec);
      eval_time += restree->children[i_c * 2]->eval_time;
      eval_time += restree->children[i_c * 2 + 1]->eval_time;
    }
    restree->eval_time = eval_time;
    restree->pfunc = pfunc;
    // debug("Parent %Lf", - LOG_FUNC(restree->pfunc));
  } else {
    // debug("Leaf %Lf", - LOG_FUNC(restree->pfunc));
  }
  return ERR_OK;
}

int update_properties(result_struc_t * res, design_struc_t * struc,
    design_spec_t * spec) {
  int i_nuc, n_nucs;

  merge_tree_pfuncs(res->tree, res, spec);

  res->pfunc = res->tree->pfunc;

  res->pfunc *= EXP_FUNC(-((BIMOLECULAR + SALT_CORRECTION) * 
        (struc->n_strands - 1)) 
      / (kB * spec->opts.temperature));

  res->pfunc /= struc->symmetry;
  res->time = res->tree->eval_time;

  n_nucs = res->n_nucs;
  for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    res->defects[i_nuc] = 0;
    res->f_defects[i_nuc] = 0;
  }

  // Map nucleotide defects from the rest of the tree to the full thing
  map_tree_defects(res->tree, res, spec);

  res->defect = 0.0;
  for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    res->defect += res->defects[i_nuc];
  }

  return ERR_OK;
}

int eval_tree(result_tree_t * res, result_struc_t * res_struc, 
    struc_tree_t * tree, design_struc_t * des_struc, design_spec_t * spec) {

  int i_c;
  // Check that this node actually represents the tree given
  int s_start, s_stop;
  int i_nuc, j_nuc, n_nucs;
  int i_seg;
  int same = 1;
  int * native_map = (int*) malloc(res_struc->n_nucs * sizeof(int));


  i_nuc = 0;
  for (i_seg = 0; i_seg < tree->n_segments; i_seg++) {
    s_start = tree->seg_start[i_seg];
    if (tree->seg_start[i_seg] != 0 && spec->opts.include_dummies) {
      s_start -= spec->opts.H_split;
    }
    s_stop = tree->seg_stop[i_seg];
    if (tree->seg_stop[i_seg] != des_struc->n_nucs && spec->opts.include_dummies) {
      s_stop += spec->opts.H_split;
    }
    for (j_nuc = s_start; j_nuc < s_stop; j_nuc++) {
      native_map[i_nuc] = j_nuc;
      i_nuc ++;
    }
  }
  n_nucs = i_nuc;
  if (n_nucs != res->n_nucs) {
    same = 0;
  }

  for (i_nuc = 0; i_nuc < res->n_nucs && same; i_nuc++) {
    if (native_map[i_nuc] != res->native_map[i_nuc]) {
      same = 0;
    }
  }

  if (tree->n_children != res->n_children) {
    same = 0;
  }
  free(native_map);

  if (same == 0) {
    debug("Reallocating: %i %i %i %i", tree->n_children, res->n_children, n_nucs, res->n_nucs);
    free_result_tree(res);
    init_result_tree(res, tree, des_struc, spec);
  }

  if (tree->n_children > 0) {
    if (res->n_children != tree->n_children) {
      for (i_c = 0; i_c < res->n_children; i_c++) {
        destroy_result_tree(res->children[i_c]);
      }
      free(res->children);

      res->children = (result_tree_t **) 
        malloc(sizeof(result_tree_t *) * tree->n_children);

      for (i_c = 0; i_c < tree->n_children; i_c++) {
        res->children[i_c] = alloc_result_tree(tree->children[i_c],
            des_struc, spec);
      }

      res->n_children = tree->n_children;
    }
    for (i_c = 0; i_c < tree->n_children; i_c++) {
      check(ERR_OK == eval_tree(res->children[i_c], res_struc, 
          tree->children[i_c], des_struc, spec),
          "Error evaluating tree");
    }
    res->pfunc = 0;
    res->eval_time = 0;
    for (i_c = 0; i_c < tree->n_children / 2; i_c++) {
      res->pfunc += compute_pfunc_merge(res->children[2*i_c],
          res->children[2*i_c + 1], res_struc, spec);
      res->eval_time += res->children[2*i_c]->eval_time;
      res->eval_time += res->children[2*i_c + 1]->eval_time;
    }
  } else {
    check(ERR_OK == eval_leaf(res, res_struc, tree, spec),
        "Error evaluating leaf node");
  }
  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

static int reduce_struc_splits(split_tracker_t * tracker, 
    int * struc, int * breaks, int n_br,
    int * to_full, int n_full_nucs, int * to_node, 
    int n_node_nucs, design_spec_t * spec) {

  int * short_struc = NULL; 
  int * split_markers = NULL; 
  int i_nuc;
  int i_br;
  int s_nuc, e_nuc;
  check(n_node_nucs > 0, "N node nucs not valid: %i", n_node_nucs);
  short_struc = (int*) malloc(n_node_nucs * sizeof(int));
  split_markers = (int*) malloc(n_node_nucs * sizeof(int));
  // Find all potential split points
  for (i_nuc = 0; i_nuc < n_node_nucs; i_nuc++) {
    if (struc[to_full[i_nuc]] >= 0) {
      short_struc[i_nuc] = to_node[struc[to_full[i_nuc]]];
    } else {
      short_struc[i_nuc] = -1;
    }
  }

  for (i_br = 0; i_br < n_br; i_br++) {
    s_nuc = breaks[i_br] - 1;
    e_nuc = breaks[i_br] + 1;
    if (s_nuc < 0) {
      s_nuc = 0;
    } 
    if (e_nuc >= n_node_nucs) {
      e_nuc = n_node_nucs;
    }
    for (i_nuc = s_nuc; i_nuc < e_nuc; i_nuc++) {
      short_struc[i_nuc] = -1;
    }
    debug("Forbidding %i to %i around native %i", to_full[s_nuc], to_full[e_nuc - 1] + 1, to_full[breaks[i_br]]);
  }
  mark_split_points(split_markers, short_struc,
      n_node_nucs, breaks, n_br, spec);

  // Reduce the current split points into the given list
  reduce_splits(tracker, split_markers, n_node_nucs, to_full, n_full_nucs,
      spec); 

  free(split_markers);
  free(short_struc);

  return ERR_OK;
error:

  free(split_markers);
  free(short_struc);
  return ERR_INVALID_STATE;
}

static int get_best_split_list(split_list_t ** best_p, 
    split_tracker_t * tracker, split_tracker_t * forbidden, int * to_full) {

  double best_cost = DBL_MAX;
  split_list_t * best = NULL;
  split_list_t * cur = tracker->head;
  split_list_t * cur_forbidden;
  split_el_t * cur_el;
  split_el_t * cur_for_el;
  int forbidden_flag = 0;
  int match_all;

  // debug("Choosing best of %i", tracker->n);

  // Loop through the set of compatible splits, take the minimal cost
  // set, up to a cost of twice the current node. 
  while (cur) {
    forbidden_flag = 0;
    cur_forbidden = forbidden->head;
    while (cur_forbidden) {
      cur_el = cur->head;

      match_all = 1;
      while (cur_el) {
        cur_for_el = cur_forbidden->head;

        while (cur_for_el) {
          if ((to_full[cur_el->lsplit] == cur_for_el->lsplit
                && to_full[cur_el->rsplit] == cur_for_el->rsplit)
              || (to_full[cur_el->lsplit] == cur_for_el->rsplit 
                && to_full[cur_el->rsplit] == cur_for_el->lsplit)) {
            // match one of the forbidden splits
            break;
          }
          cur_for_el = cur_for_el->next;
        }

        if (!cur_for_el) {
          // This wasn't forbidden by current line
          match_all = 0;
          break;
        }
        cur_el = cur_el->next;
      }
      if (match_all) {
        forbidden_flag = 1;
        // debug("Split was forbidden");
        // debug("%i %i...", to_full[cur->head->lsplit], to_full[cur->head->rsplit]);
        break;
      }
      cur_forbidden = cur_forbidden->next;
    }

    if (!forbidden_flag && cur->cost < best_cost) {
      best = cur;
      best_cost = cur->cost;
    }
    cur = cur->next;
  }

  *best_p = best;
  return ERR_OK;
}
          
static int ispresent_split_tracker(split_tracker_t * forbidden, split_list_t * test_l, 
    int * to_full) {

  split_list_t * cur_forbidden = forbidden->head;
  split_el_t * cur_el;
  split_el_t * cur_for_el;
  int match_any = 0;

  while (cur_forbidden && !match_any) {
    cur_for_el = cur_forbidden->head;
    while (cur_for_el && !match_any) {
      cur_el = test_l->head;
      while (cur_el && !match_any) {

        if ((to_full[cur_el->lsplit] == cur_for_el->lsplit
              && to_full[cur_el->rsplit] == cur_for_el->rsplit)
            || (to_full[cur_el->lsplit] == cur_for_el->rsplit 
              && to_full[cur_el->rsplit] == cur_for_el->lsplit)) {
          // match one of the forbidden splits
          match_any = 1;
          break;
        }
        cur_el = cur_el->next;
      }
      cur_for_el = cur_for_el->next;
    }
    cur_forbidden = cur_forbidden->next;
  }

  // while (cur_forbidden) {
  //   cur_el = test_l->head;

  //   match_all = 1;
  //   while (cur_el) {
  //     cur_for_el = cur_forbidden->head;

  //     while (cur_for_el) {
  //       if ((to_full[cur_el->lsplit] == cur_for_el->lsplit
  //             && to_full[cur_el->rsplit] == cur_for_el->rsplit)
  //           || (to_full[cur_el->lsplit] == cur_for_el->rsplit 
  //             && to_full[cur_el->rsplit] == cur_for_el->lsplit)) {
  //         // match one of the forbidden splits
  //         break;
  //       }
  //       cur_for_el = cur_for_el->next;
  //     }

  //     if (!cur_for_el) {
  //       // This wasn't forbidden by current line
  //       match_all = 0;
  //       break;
  //     }
  //     cur_el = cur_el->next;
  //   }
  //   if (match_all) {
  //     forbidden_flag = 1;
  //     // debug("Split is forbidden");
  //     // debug("%i %i...", to_full[test_l->head->lsplit], to_full[test_l->head->rsplit]);
  //     break;
  //   }
  //   cur_forbidden = cur_forbidden->next;
  // }
  return match_any;
}

static int reduce_ppair_splits(split_tracker_t * tracker, DBL_TYPE * ppair_full,
    int * ppair_i_full, int * ppair_j_full, int ppair_n_full, int * breaks, int n_breaks,
    int * to_full, int n_full, int * to_node, int n_node, 
    design_struc_t * struc, design_spec_t * spec, DBL_TYPE fsplit) {

  int i_pos, i_nuc, j_nuc, d_nuc;
  int i_br;
  // ignore unpaired for this
  DBL_TYPE * ppairmat = (DBL_TYPE*) calloc((n_full + 1) * n_full, sizeof(DBL_TYPE));
  DBL_TYPE split_prob;
  int m_i_nuc, m_j_nuc;
  int j_pos;

  int l1, l2;

  int ppair_n;
  split_el_t temp_el;
  split_list_t * l = NULL;
  split_list_t * best_l = NULL;
  split_list_t * test_l = NULL;
  split_el_t * el = NULL;
  int max_depth = 6;
  int cur_depth = 0;
  int j_stack;
  double max_ppair;
  int found_forbidden;
  double tmp_dbl1, tmp_dbl2;

  int * can_break = (int*) malloc(n_full * sizeof(int));
  int * nuc_ids = (int*) malloc(n_full * sizeof(int));
  // int * to_fnode = (int*) malloc(n_full * sizeof(int));
  int * ppair_i = (int*) malloc(ppair_n_full * sizeof(int));
  int * ppair_j = (int*) malloc(ppair_n_full * sizeof(int));
  int * chosen_pos_stack = (int*) malloc((max_depth + 1) * sizeof(int));
  double * cost_stack = (double *) malloc((max_depth + 1) * sizeof(double));
  double * ppair_stack = (double *) malloc((max_depth + 1) * sizeof(double));
  DBL_TYPE * ppair = (DBL_TYPE*) malloc(ppair_n_full * sizeof(DBL_TYPE));

  int ** allowed_pos_stack = (int**) malloc((max_depth + 1) * sizeof(int*));

  split_tracker_t * old_tracker = alloc_split_tracker();

  double * split_lookup = NULL;
  double * cost_lookup = NULL; 
  double best_cost = DBL_MAX;

  check_mem(can_break);
  check_mem(nuc_ids);
  // check_mem(to_fnode);
  check_mem(ppair);
  check_mem(ppair_i);
  check_mem(ppair_j);
  check_mem(ppairmat);
  check_mem(cost_stack);
  check_mem(ppair_stack);
  check_mem(chosen_pos_stack);
  check_mem(allowed_pos_stack);


  if (spec->opts.single_decomp) {
    max_depth = 1;
  }
  for (i_nuc = 0; i_nuc < max_depth + 1; i_nuc++) {
    allowed_pos_stack[i_nuc] = NULL;
  }

  for (i_nuc = 0; i_nuc < n_full; i_nuc++) {
    can_break[i_nuc] = 1;
  }

  // if (to_full[0] == 0) {
  //   j_nuc = 0;
  // } else {
  //   j_nuc = spec->opts.H_split;
  // }
  // to_fnode[0] = j_nuc;

  // for (i_nuc = 1; i_nuc < n_node; i_nuc++) {
  //   if (to_full[i_nuc - 1] == to_full[i_nuc] - 1) {
  //     j_nuc += 1;
  //   } else {
  //     j_nuc += 2 * spec->opts.H_split;
  //   }
  //   to_fnode[i_nuc] = j_nuc;
  // }
  // n_fnode_nucs = j_nuc;
  // if (to_full[0] == 0) {
  //   n_fnode_nucs += spec->opts.H_split;
  // }

  check(ERR_OK == get_nuc_ids(nuc_ids, struc, spec),
      "Error getting nucleotide ids");

  for (i_br = 0; i_br < n_breaks; i_br++) {
    i_nuc = breaks[i_br];
    check(i_nuc >= spec->opts.H_split && i_nuc < n_full - spec->opts.H_split, 
        "Invalid split point, too close to break");
    for (d_nuc = i_nuc - spec->opts.H_split; d_nuc < i_nuc + spec->opts.H_split;
        d_nuc ++) {
      can_break[d_nuc] = 0;
    }
  }

  for (d_nuc = 0; d_nuc < spec->opts.H_split; d_nuc++) {
    can_break[d_nuc] = 0;
    can_break[n_node - d_nuc - 1] = 0;
  }

  j_pos = 0;
  for (i_pos = 0; i_pos < ppair_n_full; i_pos++) {
    i_nuc = ppair_i_full[i_pos];
    j_nuc = ppair_j_full[i_pos];
    split_prob = ppair_full[i_pos];
    if (j_nuc >= 0 && nuc_ids[i_nuc] == spec->seqs.comp_map.nucs[nuc_ids[j_nuc]]) {
      ppair[j_pos] = split_prob;
      ppair_i[j_pos] = i_nuc;
      ppair_j[j_pos] = j_nuc;
      j_pos ++;
    }
  }
  ppair_n = j_pos;

  if (ppair_n > 0) {
    split_lookup = (double *) malloc(ppair_n * sizeof(double));
    cost_lookup = (double *) malloc(ppair_n * sizeof(double));
  }
  for (i_pos = 0; i_pos < ppair_n; i_pos++) {
    i_nuc = ppair_i[i_pos];
    j_nuc = ppair_j[i_pos];
    if (j_nuc >= 0) {
      if (to_node[i_nuc] >= 0 && to_node[j_nuc] >= 0) {
        ppairmat[i_nuc * (n_full + 1) + j_nuc] = ppair[i_pos];
        ppairmat[j_nuc * (n_full + 1) + i_nuc] = ppair[i_pos];
      }
    }
  }


  l = tracker->head;
  while (l) {
    l->ppair = 0;
    el = l->head;
    while (el) {

      m_i_nuc = el->lsplit;
      m_j_nuc = el->rsplit;

      i_nuc = to_full[m_i_nuc];
      j_nuc = to_full[m_j_nuc];


      l->ppair += (double)get_helix_min_ppair(ppairmat, i_nuc, j_nuc, 
          n_full, spec->opts.H_split);

      el = el->next;
    }
    l = l->next;
  }

  copy_split_tracker(old_tracker, tracker);
  free_split_tracker(tracker);

  best_l = alloc_split_list();

  l = old_tracker->head;
  while (l) {
    // debug("%lf", l->ppair);
    el = l->head;
    while (el) {
      // debug(" %i %i", el->lsplit, el->rsplit);
      el = el->next;
    }
    if (l->ppair > fsplit && l->cost < best_cost) {
      free_split_list(best_l);
      copy_split_list(best_l, l);
      best_cost = l->cost;
    }
    l = l->next;
  }
  for (cur_depth = 0; cur_depth < max_depth; cur_depth++) {
    if (ppair_n > 0) {
      allowed_pos_stack[cur_depth] = (int*) malloc(ppair_n * sizeof(int));
      for (i_pos = 0; i_pos < ppair_n; i_pos++) {
        allowed_pos_stack[cur_depth][i_pos] = 1;
      }
    } else {
      allowed_pos_stack[cur_depth] = NULL;
    }
  }

  l = old_tracker->head;

  while (l) {

    cur_depth = 0;
    for (i_pos = 0; i_pos < ppair_n; i_pos++) {
      allowed_pos_stack[cur_depth][i_pos] = 1;
    }
    for (i_pos = 0; i_pos < ppair_n; i_pos++) {
      el = l->head;
      i_nuc = ppair_i[i_pos];
      j_nuc = ppair_j[i_pos];
      if (i_nuc >= 0 && j_nuc >= 0) {
        m_i_nuc = to_node[i_nuc];
        m_j_nuc = to_node[j_nuc];
        // l1 = to_fnode[m_i_nuc] + n_fnode_nucs - to_fnode[m_j_nuc] - 1;
        // l2 = to_fnode[m_j_nuc] - to_fnode[m_i_nuc] + 1;
        
        l1 = m_i_nuc + n_node - m_j_nuc - 1;
        l2 = m_j_nuc - m_i_nuc + 1;

        tmp_dbl1 = pow((double)l1, 3);
        tmp_dbl2 = pow((double)l2, 3);
        cost_lookup[i_pos] = tmp_dbl1 + tmp_dbl2;

        allowed_pos_stack[cur_depth][i_pos] =
          allowed_pos_stack[cur_depth][i_pos] && 
          (l1 >= spec->opts.N_split) && (l2 >= spec->opts.N_split)
          && can_break[m_i_nuc] && can_break[m_j_nuc] &&
          get_helix_comp(i_nuc, j_nuc, n_full, spec->opts.H_split, nuc_ids,
              spec->seqs.comp_map.nucs);

        while (el && allowed_pos_stack[cur_depth][i_pos]) {
          allowed_pos_stack[cur_depth][i_pos] = 
              allowed_pos_stack[cur_depth][i_pos] 
            && pairs_cross(m_i_nuc, m_j_nuc, el->lsplit, el->rsplit);
          el = el->next;
        }

        if (allowed_pos_stack[cur_depth][i_pos]) {
          split_lookup[i_pos] = (double)get_helix_min_ppair(ppairmat, 
              i_nuc, j_nuc, n_full, spec->opts.H_split);
        } else {
          split_lookup[i_pos] = 0;
        }
        // debug("Cost %i (n_nucs %i): %i %i | %i %i = %lf %lf %Lf", i_pos, n_node, i_nuc, j_nuc, m_i_nuc, m_j_nuc,
        //     cost_lookup[i_pos], split_lookup[i_pos], ppair[i_pos]);
      } else {
        cost_lookup[i_pos] = DBL_MAX;
        split_lookup[i_pos] = 0;
        allowed_pos_stack[cur_depth][i_pos] = 0;
      }
    }

    cost_stack[cur_depth] = l->cost;
    ppair_stack[cur_depth] = l->ppair;
    while (cur_depth >= 0) {
      // Branch
      max_ppair = 0;
      j_pos = -1;
      for (i_pos = 0; i_pos < ppair_n; i_pos++) {
        if (max_ppair < split_lookup[i_pos] && allowed_pos_stack[cur_depth][i_pos] &&
            cost_lookup[i_pos] + cost_stack[cur_depth] < best_cost) {
          j_pos = i_pos;
          max_ppair = split_lookup[i_pos];
        }
      }

      // Found something that is below the bound
      if (j_pos >= 0) {
        i_nuc = ppair_i[j_pos];
        j_nuc = ppair_j[j_pos];
        chosen_pos_stack[cur_depth] = j_pos;

        if (split_lookup[j_pos] + ppair_stack[cur_depth] > fsplit) {
          test_l = alloc_split_list();
          copy_split_list(test_l, l);

          for (j_stack = 0; j_stack <= cur_depth; j_stack++) {
            i_pos = chosen_pos_stack[j_stack];
            temp_el.lsplit = to_node[ppair_i[i_pos]];
            temp_el.rsplit = to_node[ppair_j[i_pos]];
            temp_el.next = NULL;
            append_split_list(test_l, &temp_el);
          }

          found_forbidden = ispresent_split_tracker(struc->tree->forbidden, test_l, to_full);

          if (!found_forbidden) {
            free_split_list(best_l);
            copy_split_list(best_l, test_l);
            best_l->cost = cost_stack[cur_depth] + cost_lookup[j_pos];
            best_l->ppair = split_lookup[j_pos] + ppair_stack[cur_depth];
            best_cost = best_l->cost;
          }

          destroy_split_list(test_l);
          allowed_pos_stack[cur_depth][j_pos] = 0;

        } else if (cur_depth < max_depth - 1) {
          allowed_pos_stack[cur_depth][j_pos] = 0;

          cost_stack[cur_depth + 1] = cost_stack[cur_depth] + cost_lookup[j_pos];
          ppair_stack[cur_depth + 1] = ppair_stack[cur_depth] + split_lookup[j_pos];
          for (i_pos = 0; i_pos < ppair_n; i_pos++) {
            allowed_pos_stack[cur_depth+1][i_pos] = allowed_pos_stack[cur_depth][i_pos];
          }

          for (i_pos = 0; i_pos < ppair_n; i_pos++) {
            allowed_pos_stack[cur_depth+1][i_pos] = 
              allowed_pos_stack[cur_depth + 1][i_pos] && pairs_cross(i_nuc, j_nuc, 
                  ppair_i[i_pos], ppair_j[i_pos]);
          }

          cur_depth ++;
        } else {
          cur_depth --;
        }
      } else {
        cur_depth --;
      }
    }
    l = l->next;
  }

  if (best_l->head) {
    append_split_tracker(tracker, best_l, NULL);
  }

  l = tracker->head;
  while (l) {
    debug("%lf %lf", l->ppair, l->cost);
    el = l->head;
    while (el) {
      debug(" %i %i", to_full[el->lsplit], to_full[el->rsplit]);
      el = el->next;
    }
    l = l->next;
  }

  destroy_split_list(best_l);
  destroy_split_tracker(old_tracker);

  free(cost_stack);
  free(ppair_stack);
  free(chosen_pos_stack);

  for (cur_depth = 0; cur_depth < max_depth; cur_depth++) {
    free(allowed_pos_stack[cur_depth]);
  }
  free(allowed_pos_stack);
  free(split_lookup);
  free(cost_lookup);

  free(ppairmat);
  free(can_break);
  free(nuc_ids);
  // free(to_fnode);

  free(ppair);
  free(ppair_i);
  free(ppair_j);

  return ERR_OK;
error:
  if (best_l) destroy_split_list(best_l);
  if (old_tracker) destroy_split_tracker(old_tracker);

  free(cost_stack);
  free(ppair_stack);
  free(chosen_pos_stack);

  for (cur_depth = 0; cur_depth < max_depth; cur_depth++) {
    free(allowed_pos_stack[cur_depth]);
  }
  free(allowed_pos_stack);
  free(split_lookup);
  free(cost_lookup);


  free(ppairmat);
  free(can_break);
  free(nuc_ids);

  free(ppair);
  free(ppair_i);
  free(ppair_j);
  return ERR_INVALID_STATE;
}

/***********************************************************/
static int decompose_tree(struc_tree_t * tree, design_struc_t * struc,
    design_spec_t * spec) {
  // Find all structures that this can be a component for
  int n_strucs, i_struc;
  int n_nucs, n_node_nucs;
  int i_c;
  design_struc_t * cur_struc = struc;
  int ** strucs = NULL;
  int * to_full = NULL;
  int * to_node = NULL;
  int * breaks = NULL;
  int n_br;
  int n_leaves = 0;
  split_list_t * best;
  split_el_t * el;

  n_nucs = struc->n_nucs;

  // return if maximum depth has been reached
  n_strucs = 0;
  while (cur_struc) {
    n_strucs ++;
    cur_struc = cur_struc->next;
  }

  strucs = (int**) malloc(n_strucs * sizeof(int*));

  // look at all structures that are compatible with the current node
  n_strucs = get_decomp_strucs(strucs, tree, struc, spec);

  breaks = (int*) malloc((struc->n_strands + tree->n_segments) * sizeof(int));
  to_full = (int*) malloc(n_nucs * sizeof(int));
  to_node = (int*) malloc(n_nucs * sizeof(int));

  get_native_nuc_maps(to_full, to_node, &n_node_nucs, breaks, 
      &n_br, tree, struc, spec);

  // Figure out the best split states
  // Loop through each structure, building a compatible split list:
  // a list of splits that are mutually compatible along with their
  // total costs.
  split_tracker_t * tracker = (split_tracker_t *) alloc_split_tracker();

  split_list_t empty;
  empty.head = NULL;
  empty.tail = NULL;
  empty.cost = 0;
  empty.ppair = 0;
  empty.next = NULL;
  append_split_tracker(tracker, &empty, NULL);

  for (i_struc = 0; i_struc < n_strucs; i_struc++) {
    check(ERR_OK == reduce_struc_splits(tracker, strucs[i_struc], 
          breaks, n_br, to_full, struc->n_nucs, to_node, 
          n_node_nucs, spec),
        "error reducing split points");
  }

  check(ERR_OK == get_best_split_list(&best, tracker, struc->tree->forbidden, to_full),
    "Error getting best split list");

  // debug("Best: %p", best);

  if (best) {
    // debug("Best cost: %lf", best->cost);
    // debug("Best ppair: %lf", best->ppair);
    // map the split points back to full indices.
    el = best->head;
    while (el) {
      el->lsplit = to_full[el->lsplit];
      el->rsplit = to_full[el->rsplit];
      debug("Split: %i %i", el->lsplit, el->rsplit);
      el = el->next;
    }
    // Create the leaves and attempt to decompose them
    check(ERR_OK == create_children(tree, struc, best, spec), 
        "Error creating children");

    for (i_c = 0; i_c < tree->n_children; i_c++) {
      n_leaves += decompose_tree(tree->children[i_c], struc,
          spec);
    }
  } else {
    n_leaves += 1;
  }

  destroy_split_tracker(tracker);
  free(to_node);
  free(to_full);
  free(breaks);
  free(strucs);

  return n_leaves;
error:

  if (tracker) destroy_split_tracker(tracker);
  free(to_node);
  free(to_full);
  free(breaks);
  free(strucs);
  return 0;
}
/***********************************************************/

int init_decomposition(design_spec_t * spec) {
  int i_struc, n_strucs;
  design_struc_t * cur;
  int n_leaves = 0;

  n_strucs = spec->n_strucs;

  for (i_struc = 0; i_struc < n_strucs; i_struc++) {
    cur = spec->strucs + i_struc;
    n_leaves = decompose_tree(cur->tree, cur, spec);
#ifndef NDEBUG
    fprintf(stderr, "Decomposition:\n");
    print_full_decomposition(stderr, cur, spec, 4);
#endif
    check(n_leaves > 0, "initializing decomposition failed");
  }
  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

/***********************************************************/
split_tracker_t * alloc_split_tracker() {
  split_tracker_t * res = 
    (split_tracker_t *) malloc(sizeof(split_tracker_t));

  check_mem(res);
  init_split_tracker(res);
  return res;
error:
  return NULL;
}
/***********************************************************/

/***********************************************************/
int init_split_tracker(split_tracker_t * t) {
  t->head = NULL;
  t->tail = NULL;
  t->n = 0;
  return ERR_OK;
}
/***********************************************************/
static int pairs_cross(int i_nuc, int j_nuc, int d_nuc, int e_nuc) {
  int a1, a2, b1, b2;
  int rval;
  if (i_nuc < j_nuc) {
    a1 = i_nuc;
    a2 = j_nuc;
  } else {
    a1 = j_nuc;
    a2 = i_nuc;
  }

  if (d_nuc < e_nuc) {
    b1 = d_nuc;
    b2 = e_nuc;
  } else {
    b1 = e_nuc;
    b2 = d_nuc;
  }


  rval = (a1 < b1 && a2 > b1 && a2 < b2) ||
         (a1 > b1 && a1 < b2 && a2 > b2) ||
         (a1 == b1 && a2 != b2) || 
         (a2 == b2 && a1 != b1) ||
         (a1 == b2 && a2 != b1) ||
         (a2 == b1 && a1 != b2);
  
  return rval;
}

// static int same_node(int i_nuc, int j_nuc, int lsplit, int rsplit) {
//   return !(((i_nuc <  lsplit || i_nuc > rsplit) 
//             && (j_nuc >= lsplit && j_nuc <= rsplit))
//         ||((i_nuc >= lsplit && i_nuc <= rsplit)
//             && (j_nuc < lsplit || j_nuc > rsplit)));
// }

/***********************************************************/
int reduce_splits(split_tracker_t * tracker, int * split_markers,
    int n_nucs, int * to_full, int n_full_nucs, design_spec_t * spec) {
  int i_nuc;
  int j_nuc;
  split_list_t * l;
  split_el_t * el;
  split_el_t tmp_el;
  split_tracker_t * new_t = alloc_split_tracker();
  split_tracker_t temp_t;
  int found_match = 0;
  int found_cross = 0;
  double cur_cost;
  // int * to_fnode = (int*) malloc(n_full_nucs * sizeof(int));
  int l1, l2;
  check_mem(new_t);
  (void)to_full;
  (void)n_full_nucs;
  (void)spec;

  // check_mem(to_fnode);

  // if (to_full[0] == 0) {
  //   j_nuc = 0;
  // } else {
  //   j_nuc = spec->opts.H_split;
  // }
  // to_fnode[0] = j_nuc;

  // for (i_nuc = 1; i_nuc < n_nucs; i_nuc++) {
  //   if (to_full[i_nuc - 1] == to_full[i_nuc] - 1) {
  //     j_nuc += 1;
  //   } else {
  //     j_nuc += 2 * spec->opts.H_split;
  //   }
  //   to_fnode[i_nuc] = j_nuc;
  // }
  // n_fnode_nucs = j_nuc;
  // if (to_full[0] == 0) {
  //   n_fnode_nucs += spec->opts.H_split;
  // }

  l = tracker->head;
  // debug("Reducing");
  while (l) {
    el = l->head;
    found_match = 0;

    while (el && !found_match) {
      if (split_markers[el->lsplit] == el->rsplit) {
        found_match = 1;
      }
      el = el->next;
    }

    if (found_match) {
      append_split_tracker(new_t, l, NULL);
    } else {
      for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
        if (split_markers[i_nuc] > i_nuc) {
          el = l->head;
          j_nuc = split_markers[i_nuc];
          found_cross = 1;
          while (el && found_cross) {
            if (!pairs_cross(i_nuc, j_nuc, el->lsplit, el->rsplit)) {
              found_cross = 0;
            } 
            el = el->next;
          }
          if (found_cross) {
            cur_cost = l->cost;
            l1 = i_nuc + n_nucs - j_nuc - 1;  // to_fnode[i_nuc] + n_fnode_nucs - to_fnode[j_nuc] - 1;
            l2 = j_nuc - i_nuc + 1;           // to_fnode[j_nuc] - to_fnode[i_nuc] + 1;
            // cur_cost += pow((double)(2 * spec->opts.H_split + l1), 3.0);
            // cur_cost += pow((double)(2 * spec->opts.H_split + l2), 3.0);
            cur_cost += pow((double)(l1), 3.0);
            cur_cost += pow((double)(l2), 3.0);
            tmp_el.lsplit = i_nuc;
            tmp_el.rsplit = j_nuc;
            tmp_el.next = NULL;
            append_split_tracker(new_t, l, &tmp_el);
            new_t->tail->cost = cur_cost;
          }
        }
      }
    }

    l = l->next;
  }

  temp_t = *tracker;
  *tracker = *new_t;
  *new_t = temp_t;
  
  destroy_split_tracker(new_t);
  // free(to_fnode);

  return ERR_OK;
error:
  return ERR_OOM;
}
/***********************************************************/

/***********************************************************/
int free_split_tracker(split_tracker_t * t) {
  split_list_t * cur = t->head;
  split_list_t * prev ;

  while (cur) {
    prev = cur;
    cur = cur->next;
    destroy_split_list(prev);
  }

  t->head = NULL;
  t->tail = NULL;
  t->n = 0;
  return ERR_OK;
}
/***********************************************************/

/***********************************************************/
int destroy_split_tracker(split_tracker_t * t) {
  if (t) {
    free_split_tracker(t);
    free(t);
  }
  return ERR_OK;
}
/***********************************************************/

/***********************************************************/
int append_split_tracker(split_tracker_t * tracker, split_list_t * l,
    split_el_t * el) {

  split_list_t * tmp_list;

  tmp_list = alloc_split_list();
  copy_split_list(tmp_list, l);
  if (el) {
    append_split_list(tmp_list, el);
  }

  if (tracker->head == NULL || tracker->tail == NULL) {
    check(tracker->head == NULL, "Invalid state in split tracker");
    check(tracker->tail == NULL, "Invalid state in split tracker");
    tracker->head = tmp_list;
    tracker->tail = tmp_list;
  } else {
    tracker->tail->next = tmp_list;
    tracker->tail = tmp_list;
    tmp_list->next = NULL;
  }
  tracker->n ++;
  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}
/***********************************************************/

/***********************************************************/
int init_split_list(split_list_t * l) {
  l->cost = 0;
  l->ppair = 0;
  l->head = NULL;
  l->tail = NULL;
  l->next = NULL;
  return ERR_OK;
}

split_list_t * alloc_split_list() {
  split_list_t * res = (split_list_t *) malloc(sizeof(split_list_t));
  check_mem(res);
  check(ERR_OK == init_split_list(res), "Error initializing split list");

  return res;
error:
  return NULL;
}
/***********************************************************/

/***********************************************************/
int append_split_list(split_list_t * l, split_el_t * el) {
  split_el_t * new_el = (split_el_t *) malloc(sizeof(split_el_t));
  check(el, "Null pointer passed in for element");
  check_mem(new_el);
  new_el->lsplit = el->lsplit;
  new_el->rsplit = el->rsplit;
  new_el->next = NULL;

  if (l->head == NULL || l->tail == NULL) {
    check(l->head == NULL && l->tail == NULL,
        "split head and tail out of sync");
    l->head = new_el;
    l->tail = new_el;
  } else {
    l->tail->next = new_el;
    l->tail = new_el;
  }
  return ERR_OK;
error:
  free(new_el);
  return ERR_INVALID_STATE;
}
/***********************************************************/

/***********************************************************/
int copy_split_tracker(split_tracker_t * dest, split_tracker_t * src) {
  free_split_tracker(dest);
  split_list_t * l = src->head;
  while (l) {
    append_split_tracker(dest, l, NULL);
    l = l->next;
  }

  return ERR_OK;
}
/***********************************************************/

/***********************************************************/
int copy_split_list(split_list_t * dest, split_list_t * src) {
  split_el_t * cur_el = NULL;

  check(src && dest, "Error copying %p -> %p", src, dest);

  cur_el = src->head;
  free_split_list(dest);

  while (cur_el) {
    append_split_list(dest, cur_el);
    cur_el = cur_el->next;
  }
  dest->ppair = src->ppair;
  dest->cost = src->cost;
  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}
/***********************************************************/

/***********************************************************/
int free_split_list(split_list_t * l) {
  split_el_t * cur_el = l->head;
  split_el_t * prev_el ;

  while (cur_el) {
    prev_el = cur_el;
    cur_el = cur_el->next;
    free(prev_el);
  }
  l->head = NULL;
  l->tail = NULL;
  l->cost = 0;
  l->ppair = 0;
  return ERR_OK;
}
/***********************************************************/

/***********************************************************/
int destroy_split_list(split_list_t * l) {
  int rval = free_split_list(l);
  free(l);
  return rval;
}
/***********************************************************/

/***********************************************************/
static int get_struc_decomposition(struc_state_t * state, 
    design_struc_t * struc, 
    int maxdepth, design_spec_t * spec) {

  int i_stack;
  int n_child;
  int i_child;
  int i_ch;

  struc_tree_t ** nodestack = (struc_tree_t**) malloc((maxdepth + 1)*
      sizeof(struc_tree_t*));

  int * childstack = (int*) malloc((maxdepth  + 1)* sizeof(int));

  (void) spec;

  check(ERR_OK == copy_struc_tree(state->tree, struc->tree),
      "Error copying tree");
  check_mem(nodestack);
  check_mem(childstack);

  nodestack[0] = state->tree;
  childstack[0] = 0;
  i_stack = 0;
  state->n = 0;

  while (i_stack >= 0) {
    i_child = childstack[i_stack];
    n_child = nodestack[i_stack]->n_children;

    if (i_stack == maxdepth) {
      for (i_ch = 0; i_ch < nodestack[i_stack]->n_children; i_ch++) {
        free_struc_tree(nodestack[i_stack]->children[i_ch]);
      }
      free(nodestack[i_stack]->children);
      nodestack[i_stack]->children = NULL;
      nodestack[i_stack]->n_children = 0;

      state->n += 1;
    } else if (n_child == 0) {
      state->n += 1;
    }

    if (i_stack < maxdepth && i_child < n_child) {
      nodestack[i_stack + 1] = 
        nodestack[i_stack]->children[childstack[i_stack]];
      childstack[i_stack] ++;
      i_stack++;
      childstack[i_stack] = 0;
    } else {
      i_stack --;
    }
  }

  check(state->n > 0, "Error creating decomposition tree");

  free(nodestack);
  free(childstack);

  return ERR_OK;
error:
  free(nodestack);
  free(childstack);
  return ERR_INVALID_STATE;
}
/***********************************************************/

int get_decomposition(design_state_t * states, int maxdepth,
    design_spec_t * spec) {
  int n_strucs = states->n_strucs;
  int i_struc;

  for (i_struc = 0; i_struc < n_strucs; i_struc++) {
    check(ERR_OK == get_struc_decomposition(
          states->states + i_struc, spec->strucs + i_struc,
          maxdepth, spec), "Error getting decomposition");
  }
  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

/***********************************************************/
/*
static int fill_in_split_points(design_struc_t * cur_struc, 
    result_tree_t * restree, struc_tree_t * tree, split_el_t * split, 
    design_spec_t * spec) {
  int i_nuc, i_assume, j_nuc, start_i, stop_i, start_j, stop_j;
  int m_i_nuc;

  // clear out all base pairs contained by this node
  for (i_nuc = 0; i_nuc < restree->n_nucs; i_nuc++) {
    m_i_nuc = restree->native_map[i_nuc];
    cur_struc->struc[m_i_nuc] = -1;
  }

  // put in helices/base pairs assumed by this node
  for (i_assume = 0; i_assume < tree->n_assumed; i_assume++) {
    start_i = tree->assumed_i[i_assume] + spec->opts.H_split;
    stop_i = tree->assumed_i[i_assume] - spec->opts.H_split;

    start_j = tree->assumed_j[i_assume] - spec->opts.H_split;
    stop_j = tree->assumed_j[i_assume] + spec->opts.H_split ;

    for (i_nuc = start_i, j_nuc = start_j ; 
        i_nuc > stop_i && j_nuc < stop_j; 
        i_nuc--, j_nuc ++) {
      cur_struc->struc[i_nuc] = j_nuc;
      cur_struc->struc[j_nuc] = i_nuc;
    }
  }

  if (split) {
    // put in helices/base pairs assumed for child
    start_i = split->lsplit + spec->opts.H_split;
    stop_i = split->lsplit - spec->opts.H_split;

    start_j = split->rsplit - spec->opts.H_split;
    stop_j = split->rsplit + spec->opts.H_split;
    
    for (i_nuc = start_i, j_nuc = start_j;
        i_nuc > stop_i && j_nuc < stop_j;
        i_nuc--, j_nuc ++) {
      cur_struc->struc[i_nuc] = j_nuc;
      cur_struc->struc[j_nuc] = i_nuc;
    }
  }

  return ERR_OK;
}
*/
/***********************************************************/


/***********************************************************/
static int get_helix_comp(int i_nuc, int j_nuc, int n_nucs, 
    int h_split, int * nucs, int * comp_map) {

  int rcomp = 1;
  int i_h;
  int k_nuc, l_nuc;
  int k_id, l_id;
  for (i_h = h_split;  i_h > -h_split; i_h--) {
    k_nuc = i_nuc + i_h;
    l_nuc = j_nuc - i_h;
    if (k_nuc >= 0 && k_nuc < n_nucs && l_nuc >= 0 && l_nuc < n_nucs) {
      k_id = nucs[k_nuc];
      l_id = nucs[l_nuc];
      if (comp_map[k_id] != l_id) {
        rcomp = 0;
      }
    } else {
      rcomp = 0;
    }
  }
  return rcomp;
}
/***********************************************************/


/***********************************************************/
static DBL_TYPE get_helix_min_ppair(DBL_TYPE * ppairs, 
    int i_nuc, int j_nuc, int n_nucs, int h_split) {

  DBL_TYPE ppair = 1.0;
  int i_h;
  int k_nuc, l_nuc;
  for (i_h = h_split;  i_h > -h_split; i_h--) {
    k_nuc = i_nuc + i_h;
    l_nuc = j_nuc - i_h;
    if (l_nuc < 0 || l_nuc >= n_nucs || k_nuc < 0 || k_nuc >= n_nucs) {
      ppair = 0;
    } else if (ppairs[k_nuc * (n_nucs + 1) + l_nuc] < ppair) {
      ppair = ppairs[k_nuc * (n_nucs + 1) + l_nuc];
    }
  }
  return ppair;
}
/***********************************************************/

/***********************************************************/
int get_nuc_ids(int * nuc_ids, design_struc_t * struc, 
    design_spec_t * spec) {
  int i_nuc, j_nuc, i_str, c_str;
  sequence_t * curseq;

  j_nuc = 0;
  for (i_str = 0; i_str < struc->n_strands; i_str++) {
    c_str = struc->strands[i_str];
    curseq = spec->seqs.strands.specs + c_str;
    check(struc->n_nucs >= j_nuc + curseq->n,
        "Invalid nucleotide count in get_nuc_ids %i > %i + %i",
        struc->n_nucs, j_nuc, curseq->n);
    for (i_nuc = 0; i_nuc < curseq->n; i_nuc++) {
      nuc_ids[j_nuc] = curseq->nucs[i_nuc];
      j_nuc ++;
    }
  }

  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}
/***********************************************************/

static void print_assumed_bp(result_tree_t * restree, struc_tree_t * tree) {
  int i, j;
  int i_nuc, j_nuc;
  for (i = 0; i < tree->n_assumed; i++) {
    i_nuc = tree->assumed_i[i];
    j_nuc = tree->assumed_j[i];
    for (j = 0; j < restree->ppairs_n; j++) {
      if (restree->ppairs_i[j] == i_nuc && restree->ppairs_j[j] == j_nuc) {
        debug("%3i %3i %Lf", i_nuc, j_nuc, (long double) restree->ppairs[j]);
      }
    }
  }
}

// static int check_assumed_bp(result_tree_t * child_res, 
//     struc_tree_t * parent, struc_tree_t * child, design_struc_t * struc,
//     design_spec_t * spec) {
//   int i, j;
//   int i_nuc, j_nuc;
//   int found_assumed;
//   int success = 1;
//   split_list_t temp_split;
//   init_split_list(&temp_split);
// 
//   for (i = 0; i < child->n_assumed; i++) {
//     i_nuc = child->assumed_i[i];
//     j_nuc = child->assumed_j[i];
//     found_assumed = 0;
//     for (j = 0; j < parent->n_assumed; j++) {
//       if (parent->assumed_i[j] == i_nuc && parent->assumed_j[j] == j_nuc) {
//         found_assumed = 1;
//       }
//     }
//     if (!found_assumed) {
//       for (j = 0; j < child_res->ppairs_n; j++) {
//         if (child_res->ppairs_i[j] == i_nuc && child_res->ppairs_j[j] == j_nuc){
//           if (child_res->ppairs[j] < spec->opts.f_split) {
//             debug("Split failed: %i %i %Lf < %Lf", i_nuc, j_nuc, child_res->ppairs[j], spec->opts.f_split);
//             success = 0;
//           }
//         }
//       }
//     }
//   }
// 
//   return success;
// }

/***********************************************************/
int ppair_decompose_at_node(result_tree_t * restree, 
    result_struc_t * res, seqstate_t * seqs, struc_tree_t * tree, 
    design_struc_t * struc, design_spec_t * spec, 
    DBL_TYPE fsplit) {

  // Find all structures that this can be a component for
  int n_strucs, i_struc;
  int n_nucs, n_node_nucs;
  int i_c;
  design_struc_t * cur_struc = struc;
  int ** strucs = NULL;
  int * to_full = NULL;
  int * to_node = NULL;
  int * breaks = NULL;
  int n_br;
  split_list_t * best;
  split_el_t * el;
  split_tracker_t * tracker = NULL;
  int success = 0;
  

  n_nucs = struc->n_nucs;

  check(ERR_OK == eval_leaf(restree, res, tree, spec),
      "Error evaluating leaf");
  // return if maximum depth has been reached

  print_assumed_bp(restree, tree);
  n_strucs = 0;
  while (cur_struc) {
    n_strucs ++;
    cur_struc = cur_struc->next;
  }

  strucs = (int**) malloc(n_strucs * sizeof(int*));

  // look at all structures that are compatible with the current node
  n_strucs = get_decomp_strucs(strucs, tree, struc, spec);

  breaks = (int*) malloc((struc->n_strands + tree->n_segments) * sizeof(int));
  to_full = (int*) malloc(n_nucs * sizeof(int));
  to_node = (int*) malloc(n_nucs * sizeof(int));

  get_native_nuc_maps(to_full, to_node, &n_node_nucs, breaks, 
      &n_br, tree, struc, spec);

  // Figure out the best split states
  // Loop through each structure, building a compatible split list:
  // a list of splits that are mutually compatible along with their
  // total costs.

  while (!success) {
    success = 1;

    tracker = (split_tracker_t *) alloc_split_tracker();
    check_mem(tracker);
    split_list_t empty;
    empty.head = NULL;
    empty.tail = NULL;
    empty.cost = 0;
    empty.ppair = 0;
    empty.next = NULL;
    append_split_tracker(tracker, &empty, NULL);

    for (i_struc = 0; i_struc < n_strucs; i_struc++) {
      check(ERR_OK == reduce_struc_splits(tracker, strucs[i_struc], 
            breaks, n_br, to_full, struc->n_nucs, to_node, 
            n_node_nucs, spec),
          "error reducing split points");
    }

    check(ERR_OK == reduce_ppair_splits(tracker, 
          restree->ppairs, restree->ppairs_i, restree->ppairs_j, 
          restree->ppairs_n, breaks, n_br, to_full, struc->n_nucs,
          to_node, n_node_nucs, struc, spec, fsplit),
        "Error reducing ppair splits");

    check(ERR_OK == get_best_split_list(&best, tracker, struc->tree->forbidden, to_full),
      "Error getting best split list");

    free_struc_tree_children(tree);
    if (best) {
      // If there is a valid decomposition of low enough cost
      
      debug("Decomposing cost: %lf  ppair: %lf", best->cost, best->ppair);
      // map the split points back to full indices.
      el = best->head;
      while (el) {
        el->lsplit = to_full[el->lsplit];
        el->rsplit = to_full[el->rsplit];
        debug("i: %i  j: %i", el->lsplit, el->rsplit);

        el = el->next;
      }
      // Create the leaves and attempt to decompose them
      check(ERR_OK == create_children(tree, struc, best, spec), 
          "Error creating children");

      for (i_c = 0; i_c < restree->n_children; i_c++) {
        destroy_result_tree(restree->children[i_c]);
      }
      free(restree->children);

      // debug("N children: %i", tree->n_children);

      restree->children = (result_tree_t **) malloc(
          sizeof(result_tree_t*) * tree->n_children);

      for (i_c = 0; i_c < tree->n_children; i_c++) {
        restree->children[i_c] = alloc_result_tree(tree->children[i_c], struc,
            spec);

        eval_leaf(restree->children[i_c], res, tree->children[i_c], spec);
        // check that leaf split point satisfies
        // success = success && 
        //   check_assumed_bp(restree->children[i_c], tree, 
        //     tree->children[i_c], struc, spec);
      }
      restree->n_children = tree->n_children;

      if (success) {
        // debug("== Split successful %Lf", fsplit);
        for (i_c = 0; i_c < tree->n_children; i_c++) {
          ppair_decompose_at_node(restree->children[i_c], res,
              seqs, tree->children[i_c], struc, spec, spec->opts.f_split);
        }
      } else {
        fsplit = (best->ppair + 1.0) / 2.0;
        // debug("== Setting fsplit to %Lf", fsplit);
      }
    } else {
      if (tree->children) {
        for (i_c = 0; i_c < tree->n_children; i_c++) {
          free_struc_tree(tree->children[i_c]);
          tree->children[i_c] = NULL;
        }
        free(tree->children);
        tree->children = NULL;
        tree->n_children = 0;
      }

      if (restree->children) {
        for (i_c = 0; i_c < restree->n_children; i_c++) {
          destroy_result_tree(restree->children[i_c]);
        }
        free(restree->children);
        restree->children = NULL;
        restree->n_children = 0;
      }
      // debug("== No split found %Lf", fsplit);
    }
  }

  destroy_split_tracker(tracker);
  free(to_node);
  free(to_full);
  free(breaks);
  free(strucs);

  return ERR_OK;
error:
  if (tracker) destroy_split_tracker(tracker);
  free(to_node);
  free(to_full);
  free(breaks);
  free(strucs);
  return ERR_INVALID_STATE;
}
/***********************************************************/

/***********************************************************/
int ppair_structure(
    result_struc_t * res,
    design_struc_t * struc,
    struc_state_t * state,
    seqstate_t * seqs,
    design_spec_t * spec) {

  check(ERR_OK == ppair_decompose_at_node(res->tree, res, seqs,
      struc->tree, struc,  spec, 
      spec->opts.f_split), "Error generating ppair structure");
  (void)state;

#ifndef NDEBUG
  debug("Leaves");
  print_leaves(stderr, struc, spec, 0);
#endif
  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}
/***********************************************************/

/***********************************************************/
int init_states(
    design_state_t * states, 
    design_spec_t * spec) {
  int n_strucs = spec->n_strucs;
  int i_struc;

  states->states = (struc_state_t *) 
    malloc(sizeof(struc_state_t) * n_strucs);
  states->n_strucs = n_strucs;
  states->cap_strucs = n_strucs;

  check_mem(states->states);

  for (i_struc = 0; i_struc < n_strucs; i_struc++) {
    check( ERR_OK == init_struc_state(states->states + i_struc, 
          spec->strucs + i_struc), 
        "Error initializing state structure %i", i_struc);
  }

  return ERR_OK;
error:
  return ERR_OOM;
}
/***********************************************************/

/***********************************************************/
int free_states(
      design_state_t * states
    ) {
  int n_strucs = states->n_strucs;
  int i_struc;

  if (states) {
    for (i_struc = 0; i_struc < n_strucs; i_struc++) {
      free_struc_state(states->states + i_struc);
    }
  }
  free(states->states);

  states->states = NULL;
  states->n_strucs = 0;
  states->cap_strucs = 0;
  return ERR_OK;
}
/***********************************************************/

/***********************************************************/
int copy_states(
    design_state_t * dest,
    design_state_t * src) {
  int n_strucs = src->n_strucs;
  int i_struc;

  if (dest->cap_strucs < n_strucs) {
    dest->states = (struc_state_t *)
      realloc(dest->states, sizeof(struc_state_t) * n_strucs);
    for (i_struc = dest->cap_strucs; i_struc < n_strucs; i_struc++) {
      init_struc_state(dest->states + i_struc, src->states[i_struc].struc);
    }
    dest->cap_strucs = n_strucs;
  }

  check_mem(dest->states);

  for (i_struc = 0; i_struc < n_strucs; i_struc++) {
    check(ERR_OK == copy_struc_state(
          dest->states + i_struc, src->states + i_struc),
        "Error copying structure state %i", i_struc);
  }

  dest->n_strucs = n_strucs;

  return ERR_OK;
error:
  return ERR_OOM;
}

/***********************************************************/
int init_result_struc(
      result_struc_t * result
    ) {

  result->sequence = NULL;
  result->f_sequence = NULL;
  result->nuc_ids = NULL;
  result->structure = NULL;
  result->modifiable = NULL;
  result->defects = NULL;
  result->f_defects = NULL;

  result->n_nucs = 0;
  result->cap_nucs = 0;
  result->breaks = NULL;
  result->n_breaks = 0;
  result->pfunc = 0;
  result->time = 0;
  result->ppairs = NULL;
  result->ppairs_i = NULL;
  result->ppairs_j = NULL;
  result->ppairs_n = 0;
  result->ppairs_cap = 0;
  result->defect = 0;

  result->tree = NULL;
  return ERR_OK;
}
/***********************************************************/

/***********************************************************/
int ensure_capacity(
      result_struc_t * result,
      struc_state_t * state
    ) {
  design_struc_t * struc = state->struc;
  int n_tot_nucs = struc->n_nucs;
  int n_strands = state->struc->n_strands;

  check(n_tot_nucs > 0, "Structure length must be > 0");
  check(ERR_OK == ensure_capacity_vals(result, n_tot_nucs, 
        n_strands), "Error allocating memory for struc result");

  return ERR_OK;
error:
  return ERR_OOM;
}
/***********************************************************/

/***********************************************************/
int fill_in_sequences(
      result_struc_t * result,
      design_struc_t * struc,
      seqstate_t * seqs,
      design_spec_t * spec
    ) {

  int n_strs, n_nucs, i_str, c_str, i_nuc, j_nuc, c_nuc;
  n_strs = struc->n_strands;

  j_nuc = 0;
  for (i_str = 0; i_str < n_strs; i_str++) {
    c_str = struc->strands[i_str];
    n_nucs = spec->seqs.strands.specs[c_str].n;
    for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
      c_nuc = spec->seqs.strands.specs[c_str].nucs[i_nuc];
      result->nuc_ids[j_nuc] = c_nuc;
      result->sequence[j_nuc] = seqs->nucspec.nucs[c_nuc];
      result->f_sequence[j_nuc] = seqs->dumspec.nucs[c_nuc];

      check(result->sequence[j_nuc] >= 1 && result->sequence[j_nuc] <= 4,
          "Invalid nucleotide code %i", result->sequence[j_nuc]);
      check(result->f_sequence[j_nuc] >= 1 && result->f_sequence[j_nuc] <= 4,
          "Invalid fake nucleotide code %i", result->f_sequence[j_nuc]);
      j_nuc += 1;
    }
  }
  // j_nuc = 0;
  // for (i_str = 0; i_str < n_strs; i_str++) {
  //   c_str = struc->strands[i_str];
  //   n_nucs = spec->seqs.strands.specs[c_str].n;
  //   for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
  //     c_nuc = spec->seqs.strands.specs[c_str].nucs[i_nuc];
  //     c_nt = seqs->nucspec.nucs[c_nuc];
  //     j_nuc += 1;
  //   }
  // }
  check(j_nuc == struc->n_nucs, "Invalid sequence fill %i nucs out of %i",
      j_nuc, struc->n_nucs);
  result->n_nucs = struc->n_nucs;

  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}
/***********************************************************/

/***********************************************************/
int update_result(
      result_t * res,
      design_state_t * state,
      seqstate_t * seqs,
      design_spec_t * spec
    ) {
  int n_strucs = state->n_strucs;
  int i_struc;
  int c_struc;
  int n_tubes = res->n_tubes;
  int i_tube;

  DBL_TYPE total_defect = 0;
  DBL_TYPE tube_defect = 0;
  DBL_TYPE deficiency = 0;
  DBL_TYPE desired_x = 0;
  DBL_TYPE nuc_conc = 0;
  DBL_TYPE eval_time = 0;
  int n_nucs = 0;

  for (i_struc = 0; i_struc < n_strucs; i_struc++) {
    check(ERR_OK == update_structure_result(res->strucs + i_struc, 
        state->states + i_struc, seqs, spec),
        "Error updating the resulting structure");
  }

  if (spec->opts.include_all) {
    evaluate_undesired(res, seqs, spec);
  }

  eval_time = 0;
  for (i_struc = 0; i_struc < res->n_orderings; i_struc++) {
    if (res->struc_map[i_struc] >= 0) {
      c_struc = res->struc_map[i_struc];
      res->eval_time[i_struc] = res->strucs[c_struc].time;
    }
    eval_time += res->eval_time[i_struc];
  }

  res->root_time = eval_time;

  calculate_concentrations(res, spec);

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
error:
  return ERR_OOM;
}
/***********************************************************/

/***********************************************************/
int calculate_concentrations(
      result_t * res,
      design_spec_t * spec
    ) {
  int n_strucs;
  int m_strucs;
  int m_strucs_pruned;
  int n_tubes = res->n_tubes;
  int i_struc;
  int j_struc;
  int k_struc;
  int c_struc;
  int g_struc;
  int i_tube;
  int i_str;
  int c_str;

  double * dG = NULL;
  double * A = NULL;
  double * A_pruned = NULL;
  double * x = NULL;
  double * x0 = NULL;
  int * excluded = NULL;

  int cap_strucs = 0;
  int n_strands = spec->seqs.strands.n;
  int cur_n_strands = 0;
  int n_points = 1;
  int max_iters = 10000;
  double tol = 1e-8;
  double delta_bar = TRUST_REGION_DELTABAR;
  double eta = TRUST_REGION_ETA;
  double min_delta = 1e-12;
  int max_trial = 1000;
  int perturb_scale = 100;
  int quiet = 1;
  int write_log_file = 0;
  char *log_file = NULL;
  int seed = genrand_int32();
  int unincluded = 0;
  DBL_TYPE inc_frac;

  DBL_TYPE pfunc = 0;

  for (i_tube = 0; i_tube < n_tubes; i_tube++) {

    n_strucs = res->tubes[i_tube].n_strucs;
    m_strucs = 0;
    unincluded = 0;
    for (i_struc = 0; i_struc < n_strucs; i_struc ++) {
      g_struc = res->tubes[i_tube].generated_ind[i_struc];
      // debug("res->tubes[i_tube].included[%i]", i_struc);
      // debug("%i", res->tubes[i_tube].included[i_struc]);
      // debug("res->included[%i]", g_struc);
      // debug("%i", res->included[g_struc]);
      if (res->included[g_struc] || res->tubes[i_tube].included[i_struc]) {
        m_strucs ++;
      }
    }
    if (m_strucs > cap_strucs) {
      dG = (double *) realloc(dG, sizeof(double) * m_strucs);
      A = (double *) realloc(A, sizeof(double) * m_strucs * n_strands);
      A_pruned = (double *) realloc(A_pruned, sizeof(double) * m_strucs * n_strands);
      x = (double *) realloc(x, sizeof(double) * m_strucs);
      x0 = (double *) realloc(x0, sizeof(double) * m_strucs);
      excluded = (int * ) realloc(excluded, sizeof(int) * m_strucs);

      cap_strucs = m_strucs;
    }

    inc_frac = 1.0; 
    for (i_struc = 0; i_struc < n_strucs; i_struc++) {
      g_struc = res->tubes[i_tube].generated_ind[i_struc];
      if (!res->tubes[i_tube].included[i_struc] 
          && !res->included[g_struc]) {
        inc_frac = 1.0 - (spec->opts.f_passive * spec->tubes[i_tube].stop);
      } 
    }
    j_struc = 0;
    for (i_struc = 0; i_struc < n_strucs; i_struc++) {

      g_struc = res->tubes[i_tube].generated_ind[i_struc];

      if (res->tubes[i_tube].included[i_struc] || res->included[g_struc]) {
        for (i_str = 0; i_str < n_strands; i_str++) {
          A[i_str * m_strucs + j_struc] = 0;
        }
        cur_n_strands = res->n_strands[g_struc];
        for (i_str = 0; i_str < cur_n_strands; i_str++) {
          c_str = res->orderings[g_struc][i_str];
          A[c_str * m_strucs + j_struc] += 1;
        }
      }


      if( res->tubes[i_tube].included[i_struc]) {
        c_struc = res->tubes[i_tube].included_ind[i_struc];
        pfunc = res->strucs[c_struc].pfunc;
        dG[j_struc] = - (double) LOG_FUNC(pfunc);
        x[j_struc] = 0;
        x0[j_struc] = (double) (res->tubes[i_tube].target_x[i_struc]
          * inc_frac);
        // debug("dG[%i] = %lf", j_struc, (double) dG[j_struc]);
        j_struc ++;
      } else if (res->included[g_struc]) {
        dG[j_struc] = (double) res->dG[g_struc];
        x[j_struc] = 0;
        x0[j_struc] = (double) (res->tubes[i_tube].target_x[i_struc]
          * inc_frac);
        unincluded ++;
        j_struc ++;
      }
    }

    m_strucs_pruned = m_strucs;
    for (i_struc = 0; i_struc < m_strucs; i_struc++) {
      if (isinf(dG[i_struc])) {
        excluded[i_struc] = 1;
        dG[i_struc] = NAD_INFINITY;
        m_strucs_pruned --;
      } else {
        excluded[i_struc] = 0;
      }
    }

    j_struc = 0;
    for (i_struc = 0 ; i_struc < m_strucs; i_struc++) {
      if (!excluded[i_struc]) {
        for (i_str = 0; i_str < n_strands; i_str++) {
          A_pruned[i_str * m_strucs_pruned + j_struc] = 
            A[i_str * m_strucs + i_struc];
        }
        dG[j_struc] = dG[i_struc];
        x0[j_struc] = x0[i_struc];
        x[j_struc] = 0;
        j_struc ++;
      }
    }
  
    // Calculate the latest concentrations
    check(ERR_OK == calc_conc_from_free_energies(x, A_pruned, dG, x0, 
          n_strands, m_strucs_pruned, 
        n_points, max_iters, tol, delta_bar, eta, min_delta, max_trial,
        perturb_scale, quiet, write_log_file, log_file, seed, NULL), 
        "Concentrations did not converge");

    j_struc = 0;
    k_struc = 0;
    for (i_struc = 0; i_struc < n_strucs; i_struc++) {
      g_struc = res->tubes[i_tube].generated_ind[i_struc];
      if (res->included[g_struc] || res->tubes[i_tube].included[i_struc]) {
        if (!excluded[j_struc]) {
          res->tubes[i_tube].x[i_struc] = (DBL_TYPE) x[k_struc];
          if (isnan(res->tubes[i_tube].x[i_struc])) {
            debug("Found nan result");
          }
          k_struc++;
        } else {
          res->tubes[i_tube].x[i_struc] = 0;
        }
        j_struc ++;
      } else {
        res->tubes[i_tube].x[i_struc] = 0;
      }
    }
  }
  free(dG);
  free(A);
  free(x);
  free(x0);
  free(A_pruned);
  free(excluded);

  return ERR_OK;
error:
  free(dG);
  free(A);
  free(x);
  free(x0);
  free(A_pruned);
  free(excluded);
  return ERR_OOM;
}

/***********************************************************/
int update_structure_result(
      result_struc_t * result,
      struc_state_t * struc_state,
      seqstate_t * seqs,
      design_spec_t * spec
    ) {

  int i_nuc, n_nucs;

  // Ensure capacity for everything in the structure result 
  check(ERR_OK == ensure_capacity(result, struc_state), 
      "Error in ensure_capacity");

  // Update sequences
  check(ERR_OK == fill_in_sequences(result, struc_state->struc,
        seqs, spec),
      "Error in fill_in_sequences");

  // Update the structure
  check(ERR_OK == update_structure(result, struc_state,
        spec),
       "Error in update_structure");

  if (result->tree == NULL) {
    result->tree = alloc_result_tree(struc_state->tree, struc_state->struc, spec);
    check(result->tree, "error allocating result tree");
  }


  check(ERR_OK == eval_tree(result->tree, result, struc_state->tree, 
        struc_state->struc, spec),
      "Error evaluating tree");

  result->pfunc = result->tree->pfunc;
  result->time = result->tree->eval_time;

  result->pfunc *= EXP_FUNC(-((BIMOLECULAR + SALT_CORRECTION) * 
        (struc_state->struc->n_strands - 1)) 
      / (kB * spec->opts.temperature));

  result->pfunc /= struc_state->struc->symmetry;

  n_nucs = result->n_nucs;
  for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    result->defects[i_nuc] = 0;
    result->f_defects[i_nuc] = 0;
  }

  // Map nucleotide defects from the rest of the tree to the full thing
  map_tree_defects(result->tree, result, spec);

  result->defect = 0.0;
  for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    result->defect += result->defects[i_nuc];
  }

  return ERR_OK;
error:
  return ERR_OTHER;
}
/***********************************************************/

/***********************************************************/
DBL_TYPE compute_pfunc_merge(result_tree_t * left, result_tree_t * right,
    result_struc_t * res, design_spec_t * spec) {
  DBL_TYPE merge_pfunc = left->pfunc * right->pfunc;
  int i, j, d, e;
  int a0, a1, b0, b1;
  if (spec->opts.include_dummies) {
    i = right->native_map[spec->opts.H_split];
    j = right->native_map[right->n_nucs - spec->opts.H_split - 1];
  } else {
    i = right->native_map[0];
    j = right->native_map[right->n_nucs - 1];
  }
  d = i - 1;
  e = j + 1;
  a0 = res->sequence[d];
  a1 = res->sequence[i];
  b0 = res->sequence[j];
  b1 = res->sequence[e];

  DBL_TYPE int_ene = HelixEnergy(a0, b1, a1, b0);

  if (!spec->opts.include_dummies) {
    if (a1 != BASE_C && b0 != BASE_C) {
      int_ene -= AT_PENALTY;
    }

    if (a0 != BASE_C && b1 != BASE_C) {
      int_ene -= AT_PENALTY;
    }
  }

  merge_pfunc *= EXP_FUNC(-int_ene / (kB * spec->opts.temperature));

  // DBL_TYPE ene1, ene2, ene_res;

  // ene1 =- kB * spec->opts.temperature * LOG_FUNC(left->pfunc);
  // ene2 =- kB * spec->opts.temperature * LOG_FUNC(right->pfunc);
  // ene_res =- kB * spec->opts.temperature * LOG_FUNC(merge_pfunc);

  // debug("Energies. Left: %Lf, Right:  %Lf,  Merge: %Lf", 
  //     ene1, ene2, ene_res);

  return merge_pfunc;
}
/***********************************************************/

/***********************************************************/
static int map_tree_defects_helper(DBL_TYPE * defects, DBL_TYPE * f_defects,
    result_tree_t * res, result_struc_t * struc, design_spec_t * spec) {
  DBL_TYPE * ds = NULL; 
  DBL_TYPE * fds = NULL; 
  DBL_TYPE tot_pfunc = res->pfunc;
  DBL_TYPE temp_pfunc;
  int i_c;
  int i_nuc;
  int j_nuc;
  int n_nucs = struc->n_nucs;
  DBL_TYPE temp_defect;

  if (res->n_children > 0) {
    ds = (DBL_TYPE*) malloc(struc->n_nucs * sizeof(DBL_TYPE));
    fds = (DBL_TYPE*) malloc(struc->n_nucs * sizeof(DBL_TYPE));

    for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
      defects[i_nuc] = 0;
      f_defects[i_nuc] = 0;
    }

    for (i_c = 0; i_c < res->n_children / 2; i_c++) {
      temp_pfunc = compute_pfunc_merge(res->children[2*i_c], 
          res->children[2*i_c + 1], struc, spec);

      check(ERR_OK == map_tree_defects_helper(ds, fds, 
            res->children[2*i_c], struc, spec),
          "Error mapping tree defects");

      for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
        defects[i_nuc] += ds[i_nuc] * temp_pfunc / tot_pfunc;
        f_defects[i_nuc] += fds[i_nuc] * temp_pfunc / tot_pfunc;
      }

      check(ERR_OK == map_tree_defects_helper(ds, fds, 
            res->children[2*i_c + 1], struc, spec),
          "Error mapping tree defects");

      temp_defect = 0;
      for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
        defects[i_nuc] += ds[i_nuc] * temp_pfunc / tot_pfunc;
        f_defects[i_nuc] += fds[i_nuc] * temp_pfunc / tot_pfunc;
        temp_defect += defects[i_nuc];
      }
    }

    for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
      if (defects[i_nuc] > 1.0000000001) {
        debug("Probability at %i exceeds 1!!  %Lf", i_nuc, defects[i_nuc]);
      }
    }

  } else {
    for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
      f_defects[i_nuc] = 0;
      defects[i_nuc] = 0;
    }
    for (i_nuc = 0; i_nuc < res->n_nucs; i_nuc++) {
      j_nuc = res->native_map[i_nuc];
      if (res->dummy_flag[i_nuc]) {
        f_defects[j_nuc] = res->nuc_defects[i_nuc];
      } else {
        defects[j_nuc] = res->nuc_defects[i_nuc];
      }
    }
  }

  free(fds);
  free(ds);
  return ERR_OK;
error:
  free(fds);
  free(ds);
  return ERR_INVALID_STATE;
}
/***********************************************************/

/***********************************************************/
int map_tree_defects(result_tree_t * res, result_struc_t * struc, 
    design_spec_t * spec) {
  check(ERR_OK == map_tree_defects_helper(struc->defects, struc->f_defects, 
      res, struc, spec), 
    "Error mapping tree defects");
  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}
/***********************************************************/

/***********************************************************/
int update_structure(
    result_struc_t * res,
    struc_state_t * state,
    design_spec_t * spec) {

  int i_nuc, n_nucs;
  int i_str, c_str, n_strands;

  n_nucs = state->struc->n_nucs;
  for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    res->structure[i_nuc] = state->struc->struc[i_nuc];
  }

  i_nuc = 0;
  n_strands = state->struc->n_strands;
  for (i_str = 0; i_str < n_strands - 1; i_str++) {
    c_str = state->struc->strands[i_str];
    i_nuc += spec->seqs.strands.specs[c_str].n;
    res->breaks[i_str] = i_nuc;
  }
  

  return ERR_OK;
}
/***********************************************************/

/***********************************************************/
static int get_n_leaves_helper(result_tree_t * tree) {
  int n_d = 0;
  int i;
  if (tree->n_children == 0) {
    n_d = 1;
  } else {
    for (i = 0; i < tree->n_children; i++) {
      n_d += get_n_leaves_helper(tree->children[i]);
    }
  }
  return n_d;
}
/***********************************************************/

/***********************************************************/
int get_n_leaves(result_struc_t * res) {
  return get_n_leaves_helper(res->tree);
}
/***********************************************************/

static result_tree_t * get_result_node_helper(result_tree_t * res, int node,
    DBL_TYPE * prob, int * marker, result_struc_t * struc, design_spec_t * spec) {
  result_tree_t * ret = NULL;
  int i;
  DBL_TYPE start_prob = *prob;
  DBL_TYPE cur_prob = start_prob;
  if (res->n_children == 0) {
    if (node == *marker) {
      ret = res;
    }
    *marker += 1;
  } else {
    for (i = 0; i < res->n_children; i += 2) {

      cur_prob = start_prob * compute_pfunc_merge(res->children[i], res->children[i+1],
          struc, spec) / res->pfunc;

      ret = get_result_node_helper(res->children[i], node, &cur_prob, marker,
          struc, spec);
      if (ret) {
        break;
      }
      ret = get_result_node_helper(res->children[i+1], node, &cur_prob, marker,
          struc, spec);
      if (ret) {
        break;
      }
    }
  }

  if (ret) {
    *prob = cur_prob;
  }
  return ret;
}

/***********************************************************/
result_tree_t * get_result_node(result_t * res, int struc, int node,
    DBL_TYPE * prob, design_spec_t * spec) {
  int marker = 0;
  DBL_TYPE prob_temp = 1.0;
  result_tree_t * ret = get_result_node_helper(res->strucs[struc].tree, node, 
      &prob_temp, &marker, &(res->strucs[struc]), spec);

  if (prob) {
    *prob = prob_temp;
  }
  return ret;
}
/***********************************************************/

/***********************************************************/
DBL_TYPE get_nodal_defect(result_t * res, int i_tube, int i_struc, int i_node,
    design_spec_t * spec) {

  DBL_TYPE deficiency, desired_x, act_x;

  int i_nuc, n_nucs;
  int c_struc = res->tubes[i_tube].included_ind[i_struc];
  DBL_TYPE prob;
  result_tree_t * node = get_result_node(res, c_struc, i_node,
      &prob, spec);

  DBL_TYPE struc_defect = 0;
  int n_native = 0;

  n_nucs = node->n_nucs;
  for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    if (!node->dummy_flag[i_nuc]) {
      struc_defect += node->nuc_defects[i_nuc];
      n_native ++;
    }
  }
  desired_x = res->tubes[i_tube].target_x[i_struc];
  act_x = res->tubes[i_tube].x[i_struc];

  deficiency = act_x - desired_x;

  if (deficiency < 0) {
    act_x = desired_x;
    deficiency = 0;
  }

  if (spec->opts.disable_defect_weights) {
    struc_defect = n_native;
  }

  return prob * (deficiency * n_native + act_x * struc_defect);
}
/***********************************************************/

/***********************************************************/
int get_pairing(int * pairing, int * breaks, int * n_nucs, 
    int * n_breaks, char * dpp, int n_pos) {
  int * stack = (int*) malloc(n_pos * sizeof(int));
  int i_break = 0;
  int i_nuc = 0;
  int j_nuc = 0;
  int i_pos = 0;
  int i_stack = 0;
  for (i_pos = 0; i_pos < n_pos; i_pos++) {
    if (dpp[i_pos] == '(') {
      pairing[i_nuc] = -2;
      stack[i_stack] = i_nuc;
      i_stack ++;
      i_nuc ++;
    } else if (dpp[i_pos] == ')') {
      i_stack --;
      check(i_stack >= 0, "Invalid structure %s", dpp);
      j_nuc = stack[i_stack];
      pairing[i_nuc] = j_nuc;
      pairing[j_nuc] = i_nuc;
      i_nuc ++;
    } else if (dpp[i_pos] == '.') {
      pairing[i_nuc] = -1;
      i_nuc ++;
    } else {
      breaks[i_break] = i_nuc;
      i_break ++;
    }
  }

  *n_nucs = i_nuc;
  *n_breaks = i_break;

  free(stack);
  return ERR_OK;
error:
  free(stack);
  return ERR_INVALID_STATE;
}
/***********************************************************/

/***********************************************************/
int free_result_tube(result_tube_t * res) {
  free(res->x);
  free(res->target_x);
  free(res->target);
  free(res->included);
  free(res->included_ind);
  free(res->generated_ind);

  res->x = NULL;
  res->target_x = NULL;
  res->target = NULL;
  res->included = NULL;
  res->included_ind = NULL;
  res->generated_ind = NULL;

  res->defect = 0;
  res->n_strucs = 0;
  res->cap_strucs = 0;
  return ERR_OK;
}
/***********************************************************/

/***********************************************************/
int free_result_struc(result_struc_t * res) {
  free(res->sequence);
  free(res->f_sequence);
  free(res->nuc_ids);
  free(res->structure);
  free(res->modifiable);
  free(res->defects);
  free(res->f_defects);

  free(res->breaks);

  free(res->ppairs);
  free(res->ppairs_i);
  free(res->ppairs_j);

  if (res->tree) {
    destroy_result_tree(res->tree);
    res->tree = NULL;
  }

  res->sequence = NULL;
  res->f_sequence = NULL;
  res->nuc_ids = NULL;
  res->structure = NULL;
  res->modifiable = NULL;
  res->defects = NULL;
  res->f_defects = NULL;

  res->breaks = NULL;
  res->pfunc = 0;
  res->time = 0;
  res->ppairs = NULL;
  res->ppairs_i = NULL;
  res->ppairs_j = NULL;
  res->n_nucs = 0;
  res->cap_nucs = 0;
  res->n_breaks = 0;
  res->ppairs_n = 0;
  res->ppairs_cap = 0;
  return ERR_OK;
}
/***********************************************************/

/***********************************************************/
int free_struc(design_struc_t * s) {

  if (s->strands) free(s->strands);
  if (s->struc) free(s->struc);

  free(s->split_forbidden);
  s->split_forbidden = NULL;

  s->strands = NULL;
  s->n_strands = 0;
  s->struc = NULL;
  s->n_nucs = 0;
  s->symmetry = 1;

  if (s->next) {
    free_struc(s->next);
    free(s->next);
    s->next = NULL;
  }

  if (s->tree) {
    free_struc_tree(s->tree);
  }
  s->tree = NULL;

  return ERR_OK;
}
/***********************************************************/

/***********************************************************/
int free_sequence_spec(sequence_spec_t * s) {
  free_sequence(&(s->nucs));
  free_sequence(&(s->comp_map));
  free_seqspec_list(&(s->domains));
  free_seqspec_list(&(s->strands));
  return ERR_OK;
}
/***********************************************************/
