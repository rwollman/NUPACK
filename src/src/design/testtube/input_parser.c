#include "pathway_design.h"


static int add_domains(sequence_spec_t * seqs, cJSON * json_domains);
static int add_strands(sequence_spec_t * seqs, cJSON * json_strands);
static int add_structures(design_spec_t * spec, cJSON * json_strucs);
static int add_tubes(design_spec_t * spec, cJSON * json_tubes, 
    DBL_TYPE def_temp);
static int add_orig_domains(sequence_spec_t * seqs, cJSON * json_orig_domains);

int parse_design(char * filename, 
    design_spec_t * des,
    options_t * opts
  ) {
  char * buf = NULL;
  FILE * f = fopen(filename, "r");
  int b_to_read = BUF_SIZE;
  int b_read = BUF_SIZE;
  int tot_b_read = 0;
  cJSON * json = NULL;
  cJSON * base = NULL;
  cJSON * lengths = NULL;
  cJSON * domains = NULL;
  cJSON * strands = NULL;
  cJSON * structures = NULL;
  cJSON * tubes = NULL;
  cJSON * orig_domains = NULL;

  check(f != NULL, "Error opening file %s", filename);

  while (b_read > 0 && b_to_read == b_read) {
    buf = (char *) realloc(buf, sizeof(char) * (b_to_read + tot_b_read));
    check_mem(buf);
    b_read = fread(buf + tot_b_read, sizeof(char), b_to_read, f);
    tot_b_read += b_read;
  }
  check(feof(f), "Error reading from file");
  fclose(f);

  json = cJSON_Parse(buf);

  base = json->child;
  init_opts(opts);
  init_design_spec(des);

  while (base != NULL) {
    if (0 == strncmp(base->string, "Material", 10)) {
      if (0 == strncmp(base->valuestring, "rna", 5) ||
          0 == strncmp(base->valuestring, "rna1995", 8)) {
        opts->material = RNA;
      } else if (0 == strncmp(base->valuestring, "rna37", 5) ||
          0 == strncmp(base->valuestring, "rna1999", 8)) {
        opts->material = RNA37;
      } else if (0 == strncmp(base->valuestring, "dna", 4) ||
          0 == strncmp(base->valuestring, "dna1998", 8)) {
        opts->material = DNA;
      } else {
        opts->material = USE_SPECIFIED_PARAMETERS_FILE;
        strncpy(PARAM_FILE, base->valuestring, 100);
      }
    }

    if (0 == strncmp(base->string, "Temperature", 14)) {
      opts->temperature = base->valuedouble + ZERO_C_IN_KELVIN;
    }
    if (0 == strncmp(base->string, "RelaxationFrac", 14)) {
      opts->relaxation_frac = (DBL_TYPE) base->valuedouble;
    }
    if (0 == strncmp(base->string, "IncludeFrac", 14)) {
      debug("Include frac: %lf", base->valuedouble);
      opts->include_frac = (DBL_TYPE) base->valuedouble;
    }
    if (0 == strncmp(base->string, "Sodium", 14)) {
      opts->sodium = base->valuedouble ;
    }
    if (0 == strncmp(base->string, "Magnesium", 14)) {
      opts->magnesium = base->valuedouble ;
    }
    if (0 == strncmp(base->string, "Dangles", 14)) {
      if (0 == strncmp(base->valuestring, "all", 4)) {
        opts->dangle_type = 2;
      } else if (0 == strncmp(base->valuestring, "some", 4)) {
        opts->dangle_type = 1;
      } else if (0 == strncmp(base->valuestring, "none", 4)) {
        opts->dangle_type = 0;
      }
    }
    if (0 == strncmp(base->string, "MinSplitHelix", 14)) {
      opts->H_split = base->valueint / 2;
    }
    if (0 == strncmp(base->string, "MinLeafSize", 14)) {
      debug("Found N_split %i", base->valueint);
      opts->N_split = base->valueint;
    }
    if (0 == strncmp(base->string, "MUnfavorable", 14)) {
      opts->M_unfavorable = base->valuedouble;
    }
    if (0 == strncmp(base->string, "MLeafopt", 14)) {
      opts->M_leafopt = (DBL_TYPE) base->valuedouble;
    }
    if (0 == strncmp(base->string, "TrueConcentrations", 20)) {
      opts->use_true_concentrations = base->valueint;
    }
    if (0 == strncmp(base->string, "AllowMismatch", 14)) {
      opts->allow_mismatch = base->valueint;
    }
    if (0 == strncmp(base->string, "AllowWobble", 14)) {
      opts->allow_mismatch = base->valueint;
    }
    if (0 == strncmp(base->string, "Seed", 14)) {
      opts->seed = base->valueint;
    }
    if (0 == strncmp(base->string, "PrintLeaves", 14)) {
      opts->print_leaves = base->valueint;
    }
    if (0 == strncmp(base->string, "PrintSteps", 14)) {
      opts->print_steps = base->valueint;
    }
    if (0 == strncmp(base->string, "Lengths", 14)) {
      lengths = base->child;
    }
    if (0 == strncmp(base->string, "Domains", 14)) {
      domains = base->child;
    }
    if (0 == strncmp(base->string, "Strands", 14)) {
      strands = base->child;
    }
    if (0 == strncmp(base->string, "Structures", 14)) {
      structures = base->child;
    }
    if (0 == strncmp(base->string, "Tubes", 14)) {
      tubes = base->child;
    }
    if (0 == strncmp(base->string, "Original domains", 17)) {
      orig_domains = base->child;
    }
    if (0 == strncmp(base->string, "PPairDecomp", 12)) {
      opts->decomp_ppair = base->valuedouble;
    }
    if (0 == strncmp(base->string, "FUnincluded", 12)) {
      opts->max_unincluded = base->valuedouble;
    }
    if (0 == strncmp(base->string, "DisableDefectWeights", 20)) {
      opts->disable_defect_weights = base->valueint;
    }
    if (0 == strncmp(base->string, "DisableFocus", 13)) {
      opts->disable_focus = base->valueint;
    } 
    if (0 == strncmp(base->string, "GCInitProb", 10)) {
      opts->gc_init_prob = base->valuedouble;
    }
    if (0 == strncmp(base->string, "ForbidSplits", 13)) {
      opts->forbid_splits = base->valueint;
    }
    base = base->next;
  }

  check(lengths != NULL, "No lengths found");
  check(domains != NULL, "No domains found");
  check(strands != NULL, "No strands found");
  check(structures != NULL, "No structures found");
  check(tubes != NULL, "No structures found");

  debug("Adding domains");
  add_domains(des->seqs, domains);
  debug("Adding strands");
  add_strands(des->seqs, strands);
  debug("Adding structures");
  add_structures(des, structures);
  debug("Adding tubes");
  add_tubes(des, tubes, opts->temperature);
  debug("Design loaded");
  check(ERR_OK == add_orig_domains(des->seqs, orig_domains),
      "Failed to add original domains");

  debug("Options");
  debug("temperature: %lf", (double) opts->temperature);
  debug("sodium: %lf", (double) opts->sodium);      // sodium concentration (M)
  debug("magnesium: %lf", (double) opts->magnesium);   // magnesium concentration (M)
  debug("min_ppair_saved: %lf", (double) opts->min_ppair_saved);  // minimum pair probability saved
  debug("H_split: %i", opts->H_split);
  debug("N_split: %i",  opts->N_split);
  debug("M_unfavorable: %Lf",  opts->M_unfavorable);
  debug("M_leafopt: %Lf",  opts->M_leafopt);
  debug("use_true_concentrations: %i",  opts->use_true_concentrations);
  debug("allow_mismatch: %i",  opts->allow_mismatch);
  debug("allow_wobble: %i",  opts->allow_wobble);
  if (opts->material == RNA) {
    debug("material: rna1995");
  } else if (opts->material == RNA37) {
    debug("material: rna1999");
  } else if (opts->material == DNA) {
    debug("material: dna1998");
  } else if (opts->material == USE_SPECIFIED_PARAMETERS_FILE) {
    debug("material: %s", PARAM_FILE);
  } else {
    debug("material: undefined");
  }
  debug("seed: %u", opts->seed);
  debug("H_split: %i", opts->H_split);
  debug("N_split: %i", opts->N_split);
  debug("M_unfavorable: %Lf", opts->M_unfavorable);
  debug("M_leafopt: %Lf", opts->M_leafopt);
  debug("F_unfavorable: %Lf", opts->max_unincluded);
  debug("decomp_ppair: %Lf", opts->decomp_ppair);
  debug("gc_init_prob: %Lf", opts->gc_init_prob);
  debug("F_relax: %Lf", opts->relaxation_frac);
  debug("F_include: %Lf", opts->include_frac);



  if (opts->N_split < opts->H_split * 8) {
    log_warn("N_split < 8 * H_split found. N_split: %i H_split: %i", 
      opts->N_split, opts->H_split);
    opts->N_split = opts->H_split * 8;
    log_warn("N_split set to %i", opts->N_split);
  }

  cJSON_Delete(json);
  free(buf);

  return ERR_OK;
error:
  free(buf);
  return ERR_OTHER;
}

int add_domains(sequence_spec_t * seqs, cJSON * domains) {
  cJSON * cur = domains;

  while (cur != NULL) {
    // Throw away the dummy constant length
    add_domain_str(seqs, cur->string, cur->valuestring, 
        strlen(cur->valuestring));
    cur = cur->next;
  }

  return ERR_OK;
}

int add_strands(sequence_spec_t * seqs, cJSON * strands) {
  cJSON * cur = strands;
  cJSON * doms;
  cJSON * cur_d;
  char ** names = NULL;
  int cap_doms = 0;
  int i_dom;
  int n_doms;
  while (cur != NULL) {
    doms = cur->child;
    cur_d = doms->child;
    n_doms = 0;
    while (cur_d != NULL) {
      n_doms ++;
      cur_d = cur_d->next;
    }
    if (n_doms > cap_doms) {
      names = (char **) realloc(names, sizeof(char *) * n_doms * 2);
    }

    cur_d = doms->child;
    i_dom = 0;
    while (cur_d != NULL) {
      names[i_dom] = cur_d->valuestring;
      i_dom ++;
      cur_d = cur_d->next;
    }

    add_strand_strs(seqs, cur->string, names, n_doms);
    cur = cur->next;
  }
  free(names);
  return ERR_OK;
}

int add_structures(design_spec_t * des, cJSON * strucs) {
  cJSON * cur = strucs;
  cJSON * el;
  cJSON * strands;
  int i_strand;
  int n_strands;
  int cap_strands = 0;

  char ** strand_names = NULL;
  char * domain_struc;

  while (cur != NULL) {
    strands = NULL;
    domain_struc = NULL;

    el = cur->child;
    while (el != NULL) {
      if (strncmp(el->string, "Strands", 8) == 0) {
        strands = el->child;
      } else if (strncmp(el->string, "Structure", 10) == 0) {
        domain_struc = el->valuestring;
      }
      el = el->next;
    }

    check(strands != NULL, "No strands found");
    check(domain_struc != NULL, "No structure found");

    el = strands;
    n_strands = 0;
    while (el != NULL) {
      n_strands ++;
      el = el->next;
    }

    if (n_strands > cap_strands) {
      cap_strands = n_strands * 2;
      strand_names = (char **) realloc(strand_names, 
          sizeof(char*) * cap_strands);
    }
    el = strands;
    i_strand = 0;
    while (el != NULL) {
      strand_names[i_strand] = el->valuestring;
      i_strand++;
      el = el->next;
    }
    add_structure_str(des, cur->string, strand_names, n_strands, 
        domain_struc);
    cur = cur->next;
  }

  free(strand_names);

  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

int add_tubes(design_spec_t * des, cJSON * tubes, DBL_TYPE def_temp) {
  cJSON * cur = tubes;
  cJSON * el;
  cJSON * strucs;
  int maxsize;
  DBL_TYPE stop;

  int n_strucs = 0;
  int i_struc;
  int cap_strucs = 0;
  DBL_TYPE * struc_concs = NULL;
  char ** struc_names = NULL;
  DBL_TYPE water_conc = 
    water_density((double) (def_temp - ZERO_C_IN_KELVIN));

  if (cur != NULL) {
    el = cur->child;
    stop = -1;
    maxsize = -1;
    strucs = NULL;

    while (el != NULL) {
      if (strncmp(el->string, "Structures", 12) == 0) {
        strucs = el->child;
      } else if (strncmp(el->string, "Stop", 5) == 0) {
        stop = (DBL_TYPE) el->valuedouble;
      } else if (strncmp(el->string, "Size", 5) == 0) {
        maxsize = el->valueint;
      }
      el = el->next;
    }
    check(strucs != NULL, "No structures found");
    check(stop > 0, "No valid stop condition found");
    check(maxsize >= 0, "No valid max size found");

    el = strucs;
    n_strucs = 0;
    while (el != NULL) {
      n_strucs ++;
      el = el->next;
    }
    
    if (n_strucs > cap_strucs) {
      cap_strucs = n_strucs * 2;
      struc_concs = (DBL_TYPE *) realloc(struc_concs, 
          sizeof(DBL_TYPE) * cap_strucs);
      struc_names = (char **) realloc(struc_names, 
          sizeof(char*) * cap_strucs);
    }

    el = strucs;
    i_struc = 0;
    while (el != NULL) {
      struc_names[i_struc] = el->child->valuestring;
      struc_concs[i_struc] = (DBL_TYPE) el->child->next->valuedouble / 
        water_conc;
      i_struc++;
      el = el->next;
    }

    add_tube_str(des, cur->string, struc_names, struc_concs, n_strucs, 
        stop, maxsize);
  }

  free(struc_concs);
  free(struc_names);

  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

int add_orig_domains(sequence_spec_t * seqs, cJSON * json_orig_domains) {
  cJSON * cur = json_orig_domains;
  cJSON * rdom = NULL;

  int * real_doms = NULL;
  int * temp_real_doms = NULL;
  int * real_comp = NULL;
  int * temp_real_comp = NULL;
  int cap_r_doms = 0;
  int n_r_doms = 0;
  char * cur_name = NULL;
  char * real_name = NULL;
  int n_len = 0;
  int i_dom;
  int cur_comp = 0;

  while (cur != NULL) {
    rdom = cur->child;
    cur_name = cur->string;
    n_len = strlen(cur_name);
    if (cur_name[n_len - 1] == '*') {
      rdom = NULL;
    }
    n_r_doms = 0;
    while (rdom != NULL) {
      real_name = rdom->valuestring;
      rdom = rdom->next;
      n_len = strlen(real_name);
      cur_comp = real_name[n_len - 1] == '*';
      if (cur_comp) {
        real_name[n_len - 1] = '\0';
        n_len -= 1;
      }

      if (n_r_doms + 1 > cap_r_doms) {
        temp_real_doms = (int*) realloc(real_doms, 
            sizeof(int) * (n_r_doms + 1) * 2);
        check_mem(temp_real_doms);
        real_doms = temp_real_doms;

        temp_real_comp = (int*) realloc(real_comp, 
            sizeof(int) * (n_r_doms + 1) * 2);
        check_mem(temp_real_comp);
        real_comp = temp_real_comp;
        cap_r_doms = (n_r_doms + 1) * 2;
      }
      for (i_dom = 0; i_dom < seqs->n_domains; i_dom++) {
        if (0 == strncmp(real_name, seqs->domain_names[i_dom], n_len + 1)) {
          real_doms[n_r_doms] = i_dom;
          break;
        }
      }
      check(i_dom < seqs->n_domains, "Domain not found %s", real_name);
      real_comp[n_r_doms] = cur_comp;
      n_r_doms += 1;
    }
    if (n_r_doms > 0) {
      check(ERR_OK == add_orig_domain(seqs, cur_name, real_doms, 
            real_comp, n_r_doms),
          "Error adding original domain");
    }
    cur = cur->next;
  }
  free(real_doms);
  free(real_comp);
  return ERR_OK;
error:
  free(real_comp);
  free(real_doms);
  return ERR_OOM;
}


