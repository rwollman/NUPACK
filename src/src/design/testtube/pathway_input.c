
#include <stdio.h>
#include "parsestruc.h"
#include "npparser.h"
#include "nplexer.h"


#include "pathway_design.h"
static value_struc_t * make_val_struc(void) ;
static int collect_string_list(char *** strlist_p, int * n_strs_p,
    value_struc_t * linkedlist);
static DBL_TYPE process_stop_units(value_struc_t * cur, 
      DBL_TYPE baseval);
static DBL_TYPE process_conc_units(design_spec_t * spec, 
      value_struc_t * cur, DBL_TYPE baseval);
static int process_design(design_spec_t * spec, value_struc_t * root);
static int min_ident(int i_nuc, int * iden_map) ;
static int comp_union(int i_nuc, int j_nuc, int * iden_map, int * comp_map, int n_nucs) ;
static int constraint_union(int con1, int con2) ;
static int remap_nucleotides(design_spec_t * spec);
static int print_strand_constraints(design_spec_t * spec);

static value_struc_t * make_val_struc(void) {
  value_struc_t * val = malloc(sizeof(value_struc_t));
  val->type = NP_TT_INTEGER;
  val->intval = 0;
  val->doubleval = 0.0;
  val->strval = NULL;
  val->strlen = 0;
  val->strcap = 0;
  val->stype = NULL;
  val->id = NULL;
  val->def = NULL;
  val->next = NULL;
  return val;
}

char* np_tt_strdup(const char* str)
{
      size_t len;
      char* copy;

      len = strlen(str) + 1;
      copy = (char*)malloc(len);
      if (!copy)  return 0;
      memcpy(copy, str, len);
      return copy;
}



value_struc_t * np_tt_make_int(int input) {
  value_struc_t * val = make_val_struc();
  val->type = NP_TT_INTEGER;
  val->intval = input;
  val->doubleval = (double) input;
  return val;
}

value_struc_t * np_tt_make_double(double input) {
  value_struc_t * val = make_val_struc();
  val->type = NP_TT_FLOAT;
  val->doubleval = input;
  val->intval = (int) input;
  return val;
}

value_struc_t * np_tt_make_string(char * input) {
  value_struc_t * val = make_val_struc();
  int l = strlen(input);
  val->strval = (char*) malloc((l + 1) * sizeof(char));
  strncpy(val->strval, input, l + 1);
  val->type = NP_TT_STRING;
  val->strlen = l + 1;
  val->strcap = l + 1;
  return val;
}

value_struc_t * np_tt_append_string(value_struc_t * cur, 
    char * next) {
  int l = strlen(next);
  char * tempchar = cur->strval;
  char * nextdup = np_tt_strdup(next);
  if (l + cur->strlen + 1 > cur->strcap) {
    tempchar = (char*) realloc(cur->strval, 
          2 * (l + cur->strlen + 1));
    check_mem(tempchar);
    cur->strval = tempchar;
    cur->strcap = 2 * (l + cur->strlen + 1);
  }
  strncat(tempchar, nextdup, l+1);
  cur->strlen = l + cur->strlen + 1;
  free(nextdup);

  return cur;
error:
  free(nextdup);
  cur->strval = NULL;
  return NULL;
}

value_struc_t * np_tt_make_definition(
    value_struc_t * stype,
    value_struc_t * id,
    value_struc_t * def) {
  value_struc_t * val = make_val_struc();
  val->type = NP_TT_DEFINITION;
  val->stype = stype;
  val->id = id;
  val->def = def;
  return val;
}

value_struc_t * np_tt_append_value(value_struc_t * cur,
    value_struc_t * app) {
  value_struc_t * cv = cur;
  while (cv->next) {
    cv = cv->next;
  }
  cv->next = app;
  return cur;
}

void np_tt_destroy_value_struc(value_struc_t * val) {
  value_struc_t * cur = val;
  value_struc_t * prev ;

  while (cur) {
    np_tt_destroy_value_struc(cur->stype);
    np_tt_destroy_value_struc(cur->id);
    np_tt_destroy_value_struc(cur->def);
    if (cur->strval && cur->strcap > 0) {
      free(cur->strval);
    }
    prev = cur;
    cur = cur->next;
    free(prev);
  }
}

void np_tt_print_ast(value_struc_t * root, int indent) {
  value_struc_t * cur = root;
  char * indent_string = (char *) malloc((indent + 1) 
      * sizeof(char));
  int i;
  for (i = 0; i < indent; i++) {
    indent_string[i] = ' ';
  }
  indent_string[indent] = '\0';
  while (cur) {
    if (cur->type == NP_TT_DEFINITION) {
      printf("%sDEFINITION:\n", indent_string);
      np_tt_print_ast(cur->stype, indent + 4);
      np_tt_print_ast(cur->id, indent + 4);
      np_tt_print_ast(cur->def, indent + 4);
    } else if (cur->type == NP_TT_STRING) {
      printf("%s%s\n", indent_string, cur->strval);
    } else if (cur->type == NP_TT_FLOAT) {
      printf("%s%lf\n", indent_string, cur->doubleval);
    } else if (cur->type == NP_TT_INTEGER) {
      printf("%s%i\n", indent_string, cur->intval);
    }
    cur = cur->next;
  }
  free(indent_string);
}


static value_struc_t * get_ast(const char * parse_string) {
  value_struc_t * vals;
  yyscan_t scanner;
  YY_BUFFER_STATE state;

  if (yylex_init(&scanner)) {
    return NULL;
  }

  state = yy_scan_string(parse_string, scanner);

  vals = NULL;
  if (yyparse(&vals, scanner)) {
    return NULL;
  }

  yy_delete_buffer(state, scanner);

  yylex_destroy(scanner);

  return vals;
}

static int collect_string_list(char *** strlist_p, int * n_strs_p,
    value_struc_t * linkedlist) {
  char ** strlist = NULL;
  int cap_strs = *n_strs_p;
  int n_strs = 0;
  char ** templist = NULL;
  value_struc_t * cur = linkedlist;

  check(strlist_p, "NULL passed to collect string list");
  check(n_strs_p, "NULL passed to collect string list");

  strlist = *strlist_p;
  while (cur) {
    if (n_strs >= cap_strs) {
      cap_strs = (n_strs + 1) * 2;
      templist = (char**)realloc(strlist, cap_strs * sizeof(char*));
      check_mem(templist);
      strlist = templist;
    }
    strlist[n_strs] = cur->strval;
    n_strs ++;
    cur = cur->next;
  }
  *strlist_p = strlist;
  *n_strs_p = n_strs;

  return ERR_OK;
error:
  if(strlist) free(strlist);
  if (strlist_p) *strlist_p = NULL;
  if (n_strs_p) *n_strs_p = 0;

  return ERR_INVALID_STATE;
}

static DBL_TYPE process_stop_units(
      value_struc_t * cur, 
      DBL_TYPE baseval
    ) {
  DBL_TYPE multiplier = 1.0;
  DBL_TYPE multipliers[] = {
    1.0,
    0.01
  };

  const char units[][20] = {
    "frac",
    "%"
  };

  int i;
  int n_units = 2;

  if (cur) {
    for (i = 0; i < n_units; i++) {
      if (0 == strncmp(cur->strval, units[i], 5)) {
        break;
      }
    }
    check(i < n_units, "Invalid stop units specified: %s", cur->strval);
    multiplier = multipliers[i];
  }

  return baseval * multiplier;
error:
  return 0.01;
}


static DBL_TYPE process_conc_units(
      design_spec_t * spec, 
      value_struc_t * cur, 
      DBL_TYPE baseval
    ) {
  DBL_TYPE multiplier = 1.0;
  DBL_TYPE multipliers[] = {
    1.0,
    1.0,
    1.0e-1,
    1.0e-2,
    1.0e-3,
    1.0e-6,
    1.0e-9,
    1.0e-12,
    1.0e-15,
    1.0e-18,
    1.0e-21,
    1.0e-24,
  };

  const char units[][20] = {
    "frac",
    "M",
    "dM",
    "cM",
    "mM",
    "uM",
    "nM",
    "pM",
    "fM",
    "aM",
    "zM",
    "yM"
  };

  int n_units = 12;
  DBL_TYPE water_conc = 0.0;
  int i;
  water_conc = water_density(
      (double)(spec->opts.temperature - ZERO_C_IN_KELVIN));
  multiplier = 1 / water_conc;

  if (cur) {
    for (i = 1; i < n_units; i++) {
      multipliers[i] /= water_conc;
    }
    for (i = 0; i < n_units; i++) {
      if (0 == strncmp(cur->strval, units[i], 5)) {
        break;
      }
    }
    multiplier = multipliers[i];
  }

  return baseval * multiplier;
}

static int process_design(design_spec_t * spec, value_struc_t * root) {
  int i_line = 1;
  value_struc_t * cur = root;
  // value_struc_t * tmp;
  char ** strlist = NULL;
  int nstrs = 0;
  DBL_TYPE curdbl = 0;
  int hsplit_set = 0;

  while (cur) {
    if (cur->type == NP_TT_DEFINITION) {
      switch(cur->stype->intval) {
        /* Physical parameters */
        case TOK_TEMPERATURE:
          if (!cur->stype->next 
              || strncmp(cur->stype->next->strval, "C", 2) == 0) {
            spec->opts.temperature = (DBL_TYPE) cur->def->doubleval + 
              ZERO_C_IN_KELVIN;
          } else if (strncmp(cur->stype->next->strval, "K", 2) == 0) {
            spec->opts.temperature = (DBL_TYPE) cur->def->doubleval;
          } else {
            sentinel("Invalid units for temperature.");
          }
          break;
        case TOK_MATERIAL:
          if (0 == strncmp(cur->def->strval, "rna", 5) ||
              0 == strncmp(cur->def->strval, "rna1995", 8)) {
            spec->opts.material = RNA;
            if (!hsplit_set) {
              spec->opts.H_split = 2;
            }
          } else if (0 == strncmp(cur->def->strval, "rna37", 5) ||
              0 == strncmp(cur->def->strval, "rna1999", 8)) {
            spec->opts.material = RNA37;
            if (!hsplit_set) {
              spec->opts.H_split = 2;
            }
          } else if (0 == strncmp(cur->def->strval, "dna", 4) ||
              0 == strncmp(cur->def->strval, "dna1998", 8)) {
            spec->opts.material = DNA;
            if (!hsplit_set) {
              spec->opts.H_split = 3;
            }
          } else {
            spec->opts.material = USE_SPECIFIED_PARAMETERS_FILE;
            strncpy(PARAM_FILE, cur->def->strval, 99);
          }
          break;
        case TOK_DGCLAMP:
          spec->opts.bonus_per_split = (DBL_TYPE) cur->def->doubleval;
          break;
        case TOK_SODIUM:
          spec->opts.sodium = (DBL_TYPE) cur->def->doubleval;
          break;
        case TOK_MAGNESIUM:
          spec->opts.magnesium = (DBL_TYPE) cur->def->doubleval;
          break;
        case TOK_DANGLES:
          if (0 == strncmp(cur->def->strval, "all", 4)) {
            spec->opts.dangle_type = 2;
          } else if (0 == strncmp(cur->def->strval, "some", 4)) {
            spec->opts.dangle_type = 1;
          } else if (0 == strncmp(cur->def->strval, "none", 4)) {
            spec->opts.dangle_type = 0;
          }
          break;

        /* Optimization parameters */
        case TOK_SEED:
          spec->opts.seed = (unsigned int) cur->def->intval;
          break;
        case TOK_HSPLIT:
          hsplit_set = 1;
          spec->opts.H_split = (int) cur->def->intval;
          break;
        case TOK_NSPLIT:
          spec->opts.N_split = (int) cur->def->intval;
          break;
        case TOK_MUNFAVORABLE:
          spec->opts.M_unfavorable = (DBL_TYPE) cur->def->doubleval;
          break;
        case TOK_MLEAFOPT:
          spec->opts.M_leafopt = (DBL_TYPE) cur->def->doubleval;
          break;
        case TOK_MRESEED:
          spec->opts.M_reseed = (DBL_TYPE) cur->def->doubleval;
          break;
        case TOK_FSPLIT:
          spec->opts.f_split = (DBL_TYPE) cur->def->doubleval;
          break;
        case TOK_FREDECOMP:
          spec->opts.f_redecomp = (DBL_TYPE) cur->def->doubleval;
          break;
        case TOK_FREFOCUS:
          spec->opts.f_refocus = (DBL_TYPE) cur->def->doubleval;
          break;
        case TOK_FSTRINGENT:
          spec->opts.f_stringent = (DBL_TYPE) cur->def->doubleval;
          break;
        case TOK_FPASSIVE:
          spec->opts.f_passive = (DBL_TYPE) cur->def->doubleval;
          break;
        case TOK_GC_INIT:
          spec->opts.gc_init_prob = (DBL_TYPE) cur->def->doubleval;
          break;
        /* Optimization modes */
        case TOK_ALLOWWOBBLE:
          spec->opts.allow_wobble = (int) cur->def->intval;
          break;
        case TOK_ALLOWMISMATCH:
          spec->opts.allow_mismatch = (int) cur->def->intval;
          break;
        case TOK_DISABLEMUTWEIGHTS:
          spec->opts.disable_defect_weights = (int) cur->def->intval;
          break;
        case TOK_SINGLE_DECOMP:
          spec->opts.single_decomp = (int) cur->def->intval;
          break;
        case TOK_INCLUDE_ALL:
          spec->opts.include_all = (int) cur->def->intval;
          break;
        /* Output parameters */
        case TOK_PRINTLEAVES:
          spec->opts.print_leaves = (int) cur->def->intval;
          break;
        case TOK_PRINTSTEPS:
          spec->opts.print_steps = (int) cur->def->intval;
          break;
        case TOK_MINPAIR:
          spec->opts.min_ppair_saved = (DBL_TYPE) cur->def->doubleval;
          break;
        case TOK_STRUCTURE:
          check(ERR_OK == 
              add_structure_basic(spec, cur->id->strval, cur->def->strval),
              "Error adding structure %s on line %i", cur->id->strval, i_line);
          break;
        case TOK_DOMAIN:
          check(ERR_OK == 
              add_domain(spec, cur->id->strval, cur->def->strval),
              "Error adding domain %s", cur->id->strval);
          break;
        case TOK_STRAND:
          check(ERR_OK ==
              collect_string_list(&strlist, &nstrs, cur->def),
              "Error preprocessing strand %s on line %i", 
              cur->id->strval, i_line);
          check(ERR_OK == add_strand(spec, cur->id->strval,
                strlist, nstrs),
              "Error adding strand %s", cur->id->strval);
          break;
        case TOK_TUBE:
          check(ERR_OK ==
              collect_string_list(&strlist, &nstrs, cur->def),
              "Error preprocessing tube %s on line %i", cur->id->strval, i_line);
          check(ERR_OK == 
              add_tube_basic(spec, cur->id->strval, strlist, nstrs),
              "Error adding tube %s on line %i", cur->id->strval, i_line);
          check(spec->n_tubes == 1, "Only one tube allowed %s isn't allowed",
              cur->id->strval);

          break;
        case TOK_CONCDEF:
          curdbl = (DBL_TYPE)cur->def->doubleval;
          curdbl = process_conc_units(spec, cur->stype->next, 
              curdbl);

          check(curdbl >= -1e-50, "Invalid concentration/units specified on"
             " line %i", i_line);

          check(ERR_OK == 
              set_concentration(spec, cur->id->strval, 
                cur->id->next->strval, curdbl),
              "Error setting concentration on line %i", i_line);

          break;
        case TOK_SEQ:
          check(ERR_OK ==
              collect_string_list(&strlist, &nstrs, cur->def),
              "Error preprocessing strand %s on line %i", 
              cur->id->strval, i_line);
          check(ERR_OK == define_structure_strands(spec, cur->id->strval,
                strlist, nstrs),
              "Error defining sequence for %s on line %i", 
              cur->id->strval, i_line);
          break;
        case TOK_STOPDEF:
          curdbl = (DBL_TYPE)cur->def->doubleval;
          curdbl = process_stop_units(cur->stype->next,
              curdbl);

          check(curdbl >= -1e-50, "Invalid stop condition specified on"
              "line %i", i_line);
          // i_tube = index_of(cur->id->strval, spec->tube_names, spec->n_tubes)
          spec->tubes[0].stop = curdbl;
          break;
        case TOK_MAXSIZE:
          spec->tubes[0].maxsize = cur->def->intval;
          break;
        case TOK_TRIALS:
          break;
        case TOK_OPTTIME:
          spec->opts.allowed_opt_time = (DBL_TYPE) cur->def->doubleval;
          break;
        default:
          sentinel("Invalid definition code %i on line %i", 
              cur->stype->intval, i_line);
      }
    }
    i_line += 1;
    cur = cur->next;
  }

  if (spec->n_tubes == 0) {
    debug("Adding structure tube");
    add_tube_basic(spec, "T", spec->struc_names, spec->n_strucs);
  }

  free(strlist);

  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}

static int min_ident(int i_nuc, int * iden_map) {
  while (iden_map[i_nuc] < i_nuc) {
    i_nuc = iden_map[i_nuc];
  }
  if (iden_map[i_nuc] != i_nuc) {
    log_err("Invalid identical map %i -> %i", i_nuc, iden_map[i_nuc]);
  }
  return i_nuc;
}

static int comp_union(int i_nuc, int j_nuc, int * iden_map, 
    int * comp_map, int n_nucs) {
  int id1 = min_ident(i_nuc, iden_map);
  int id2 = min_ident(j_nuc, iden_map);
  int id3 = comp_map[id1];
  int id4 = comp_map[id2];

  int id5 = id1;
  int id6 = id2;

  if (id4 < n_nucs) {
    if (id1 < id4) {
      iden_map[id4] = id1;
      id5 = id1;
    } else {
      iden_map[id1] = id4;
      id5 = id4;
    }
  }
  if (id3 < n_nucs) {
    if (id2 < id3) {
      iden_map[id3] = id2;
      id6 = id2;
    } else {
      iden_map[id2] = id3;
      id6 = id3;
    }
  }

  comp_map[id5] = id6;
  comp_map[id6] = id5;

  return ERR_OK;
}

static int constraint_union(int con1, int con2) {
  // first bit A, second bit C, third bit G, fourth bit U
  int codes[] = {
    BASE_N,
    BASE_A,
    BASE_C,
    BASE_G,
    BASE_U,
    BASE_M,
    BASE_R,
    BASE_W,
    BASE_S,
    BASE_Y,
    BASE_K,
    BASE_V,
    BASE_H,
    BASE_D,
    BASE_B
  };
  int bit_masks[] = {
    0xF, // 0 N : A | C | G | U
    0x1, // 1 A : A
    0x2, // 2 C : C
    0x4, // 3 G : G
    0x8, // 4 U : U
    0x3, // 5 M : A | C
    0x5, // 6 R : A | G
    0x9, // 7 W : A | U
    0x6, // 8 S : C | G
    0xA, // 9 Y : C | U
    0xC, // a K : G | U
    0x7, // b V : A | C | G
    0xB, // c H : A | C | U
    0xD, // d D : A | G | U
    0xE, // e B : C | G | U
    0x0,
  };
  int rev_bit_masks[] = {
    -1, // invalid sequence, nothing remaining
    BASE_A,
    BASE_C,
    BASE_M,
    BASE_G,
    BASE_R,
    BASE_S,
    BASE_V,
    BASE_U,
    BASE_W,
    BASE_Y,
    BASE_H,
    BASE_K,
    BASE_D,
    BASE_B,
    BASE_N
  };
  int i;
  int bm1, bm2;
  for (i = 0; i < 14; i++) {
    if (con1 == codes[i]) {
      break;
    }
  }
  bm1 = bit_masks[i];
  for (i = 0; i < 14; i++) {
    if (con2 == codes[i]) {
      break;
    }
  }
  bm2 = bit_masks[i];

  int mn = rev_bit_masks[bm1 & bm2];

  if (mn < 0 || mn > 14) {
    log_err("Inconsistent constraints specified %i and %i on identical nucs", 
        con1, con2);
  }
  return mn;
}

static int remap_nucleotides(design_spec_t * spec) {
  int n_nucs;
  // post-processing, iden_map[i] will be the minimum nucleotide 
  // guaranteed to be identical to nucleotide i
  int * iden_map = NULL; 
  // post-processing, comp_map[i] will be the minimum nucleotide 
  // guaranteed to be the complement of nucleotide i 
  // (unless mismatch-allowed is specified)
  int * comp_map = NULL;
  int * used_flag = NULL;
  int * final_map = NULL;
  int i;
  int i_nuc, j_nuc, k_nuc;
  int i_dom;
  int j_dom;
  int c_n_nucs;

  int * struc_nucs = NULL;
  int * strand_used = NULL;
  int name_len;
  char * comp_name = NULL; 
  int * new_cons = NULL;
  int * new_comps = NULL;
  int n_new_nucs = 0;
  int i_struc, n_strs, i_str, j_str, c_str;
  design_struc_t * cur_struc;

  n_nucs = spec->seqs.nucs.n;

  iden_map = (int*) malloc(n_nucs * sizeof(int));
  comp_map = (int*) malloc(n_nucs * sizeof(int));
  used_flag = (int*) malloc(n_nucs * sizeof(int));
  strand_used = (int*) malloc(spec->seqs.strands.n * sizeof(int));
  final_map = (int*) malloc(n_nucs * sizeof(int));

  for (i = 0; i < n_nucs; i++) {
    iden_map[i] = i;
    comp_map[i] = n_nucs + 1;
    used_flag[i] = 0;
    final_map[i] = -1;
  }
  
  /* Map each domain nucleotide to its complement */
  for (i_dom = 0; i_dom < spec->seqs.domains.n; i_dom++) {
    c_n_nucs = spec->seqs.domains.specs[i_dom].n;
    name_len = strlen(spec->seqs.domains.names[i_dom]);
    comp_name = (char *) malloc((name_len + 2) * sizeof(char));
    strncpy(comp_name, spec->seqs.domains.names[i_dom], name_len + 1);
    if (comp_name[name_len - 1] == '*') {
      comp_name[name_len - 1] = '\0';
    } else {
      comp_name[name_len] = '*';
      comp_name[name_len + 1] = '\0';
    }
    j_dom = index_of(comp_name, spec->seqs.domains.names, 
        spec->seqs.domains.n);
    check(j_dom >= 0, "Domain %s was not created", comp_name);

    if (j_dom > i_dom) {
      check(spec->seqs.domains.specs[i_dom].n == 
          spec->seqs.domains.specs[j_dom].n,
          "Domain %s and its complement have different nucleotide counts",
          comp_name);

      for (i = 0; i < c_n_nucs; i++) {
        i_nuc = spec->seqs.domains.specs[i_dom].nucs[i];
        j_nuc = spec->seqs.domains.specs[j_dom].nucs[c_n_nucs - i - 1];
        comp_map[i_nuc] = j_nuc < comp_map[i_nuc] ? j_nuc: comp_map[i_nuc];
        comp_map[j_nuc] = i_nuc < comp_map[j_nuc] ? i_nuc: comp_map[j_nuc];
      }
    }
    free(comp_name);
  }

  /* Map each nucleotide in a duplex structure to its complement */
  for (i_struc = 0; i_struc < spec->n_strucs; i_struc++) {
    struc_nucs = (int *) malloc(spec->strucs[i_struc].n_nucs * sizeof(int));
    n_strs = spec->strucs[i_struc].n_strands;
    i_nuc = 0;

    for (i_str = 0; i_str < n_strs; i_str++) {
      c_str = spec->strucs[i_struc].strands[i_str];
      c_n_nucs = spec->seqs.strands.specs[c_str].n;
      for (i = 0; i < c_n_nucs; i++) {
        struc_nucs[i_nuc] = spec->seqs.strands.specs[c_str].nucs[i];
        i_nuc ++;
      }
    }
    for (i_nuc = 0; i_nuc < spec->strucs[i_struc].n_nucs; i_nuc++) {
      if (spec->strucs[i_struc].struc[i_nuc] >= 0) {
        check(ERR_OK == comp_union(
              struc_nucs[i_nuc], 
              struc_nucs[spec->strucs[i_struc].struc[i_nuc]], 
            iden_map, comp_map, n_nucs),
            "Invalid comp union attempted");
      }
    }
    free(struc_nucs);

  }

  for (i_str = 0; i_str < spec->seqs.strands.n; i_str++) {
    strand_used[i_str] = 0;
  }

  for (i_struc = 0; i_struc < spec->n_strucs; i_struc++) {
    struc_nucs = (int *) malloc(spec->strucs[i_struc].n_nucs * sizeof(int));
    n_strs = spec->strucs[i_struc].n_strands;
    i_nuc = 0;

    for (i_str = 0; i_str < n_strs; i_str++) {
      c_str = spec->strucs[i_struc].strands[i_str];
      strand_used[c_str] = 1;
      c_n_nucs = spec->seqs.strands.specs[c_str].n;
      for (i = 0; i < c_n_nucs; i++) {
        struc_nucs[i_nuc] = spec->seqs.strands.specs[c_str].nucs[i];
        i_nuc ++;
      }
    }
    for (i_nuc = 0; i_nuc < spec->strucs[i_struc].n_nucs; i_nuc++) {
      j_nuc = min_ident(struc_nucs[i_nuc], iden_map);
      used_flag[j_nuc] = 1;
    }
    free(struc_nucs);
  }

  for (i_dom = 0; i_dom < spec->seqs.domains.n; i_dom++) {
    for (i_nuc= 0; i_nuc < spec->seqs.domains.specs[i_dom].n; i_nuc++) {
      j_nuc = spec->seqs.domains.specs[i_dom].nucs[i_nuc];
      j_nuc = min_ident(j_nuc, iden_map);
      used_flag[j_nuc] = 1;
    }
  }

  for (i_str = 0; i_str < spec->seqs.strands.n; i_str++) {
    if (strand_used[i_str] == 0 || spec->seqs.strands.names[i_str][0] != '-') {
      for (i_nuc= 0; i_nuc < spec->seqs.strands.specs[i_str].n; i_nuc++) {
        j_nuc = spec->seqs.strands.specs[i_str].nucs[i_nuc];
        j_nuc = min_ident(j_nuc, iden_map);
        used_flag[j_nuc] = 1;
      }
    }
  }

  for (i_str = spec->seqs.strands.n - 1; i_str >= 0 ; i_str--) {
    check(i_str < spec->seqs.strands.n && i_str >= 0, "Invalid loop check (for clang analyzer)");
    if (strand_used[i_str] == 0 && spec->seqs.strands.names[i_str][0] == '-') {
      delete_sequence_seqspec(&(spec->seqs.strands), i_str);
      debug("Deleting strand %i", i_str);
      for (i_struc = 0; i_struc < spec->n_strucs; i_struc++) {
        cur_struc = spec->strucs + i_struc;
        while (cur_struc) {
          for (j_str = 0; j_str < cur_struc->n_strands; j_str++) {
            c_str = cur_struc->strands[j_str];
            // debug("Structure strand %i", c_str);
            if (c_str > i_str) {
              // debug("Moved strand %i in struc %i", c_str, i_struc);
              cur_struc->strands[j_str] -= 1;
            }
          }
          cur_struc = cur_struc->next;
        }
      }
      for (i_struc = 0; i_struc < spec->n_orderings; i_struc++) {
        for (j_str = 0; j_str < spec->n_strands[i_struc]; j_str++) {
          c_str = spec->orderings[i_struc][j_str];
          if (c_str > i_str) {
            // debug("Moved strand %i in ordering %i", c_str, i_struc);
            spec->orderings[i_struc][j_str] -= 1;
          }
        }
      }
    }
  }

  j_nuc = 0;
  for (i = 0; i < n_nucs; i++) {
    i_nuc = min_ident(i, iden_map);
    if (used_flag[i_nuc]) {
      if (i_nuc == i) {
        final_map[i] = j_nuc;
        j_nuc += 1;
      } else {
        final_map[i] = final_map[i_nuc];
      }
    } else {
      final_map[i] = -1;
    }
  }
  n_new_nucs = j_nuc;
  check(n_new_nucs > 0, "Error simplifying constraints. No nucleotides left!");
  new_cons = (int*) malloc(j_nuc * sizeof(int));
  new_comps = (int*) malloc(j_nuc * sizeof(int));

  for (i_dom = 0; i_dom < spec->seqs.domains.n; i_dom++) {
    for (i_nuc= 0; i_nuc < spec->seqs.domains.specs[i_dom].n; i_nuc++) {
      j_nuc = spec->seqs.domains.specs[i_dom].nucs[i_nuc];
      spec->seqs.domains.specs[i_dom].nucs[i_nuc] = final_map[j_nuc];
    }
  }

  for (i_str = 0; i_str < spec->seqs.strands.n; i_str++) {
    for (i_nuc = 0; i_nuc < spec->seqs.strands.specs[i_str].n; i_nuc++) {
      j_nuc = spec->seqs.strands.specs[i_str].nucs[i_nuc];
      spec->seqs.strands.specs[i_str].nucs[i_nuc] = final_map[j_nuc];
    }
  }
  // Map the sequence constraints down

  for (i_nuc = 0; i_nuc < n_new_nucs; i_nuc++) {
    new_cons[i_nuc] = 0;
    new_comps[i_nuc] = -1;
  }

  for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
    j_nuc = final_map[i_nuc];
    k_nuc = min_ident(i_nuc, iden_map); // id in iden_map and comp_map
    if (j_nuc >= 0) {
      new_cons[j_nuc] = 
        constraint_union(new_cons[j_nuc], 
            spec->seqs.nucs.nucs[i_nuc]);
      if (comp_map[k_nuc] >= 0 && comp_map[k_nuc] < n_nucs) {
        if (new_comps[j_nuc] < 0) {
          new_comps[j_nuc] = final_map[comp_map[k_nuc]];
          new_comps[final_map[comp_map[k_nuc]]] = j_nuc;
        }
      }
    }
  }

  spec->seqs.nucs.n = n_new_nucs;
  spec->seqs.comp_map.n = n_new_nucs;
  free(spec->seqs.nucs.nucs);
  free(spec->seqs.comp_map.nucs);

  spec->seqs.nucs.nucs = new_cons;
  spec->seqs.comp_map.nucs = new_comps;

  free(iden_map);
  free(comp_map);
  free(used_flag);
  free(strand_used);
  free(final_map);

  return ERR_OK;
error:
  free(iden_map);
  free(comp_map);
  free(used_flag);
  free(strand_used);
  free(final_map);
  return ERR_INVALID_STATE;
}

static int print_strand_constraints(design_spec_t * spec) {
  int i_str, n_strs, i_nuc, n_nucs, c_nuc, c_nt;
  char * curseq = NULL;
  int * cur_cons = NULL;
  n_strs = spec->seqs.strands.n;
  for (i_str = 0; i_str < n_strs; i_str++) {
    n_nucs = spec->seqs.strands.specs[i_str].n;
    cur_cons = (int*) malloc(n_nucs * sizeof(int));
    curseq = (char*) malloc((n_nucs + 1) * sizeof(char)); 

    for (i_nuc = 0; i_nuc < n_nucs; i_nuc++) {
      c_nuc = spec->seqs.strands.specs[i_str].nucs[i_nuc];
      c_nt = spec->seqs.nucs.nucs[c_nuc];
      cur_cons[i_nuc] = c_nt;
    }
    convert_nucs_to_str(curseq, cur_cons,
        n_nucs, &(spec->opts));
    debug("Strand %s: %s", spec->seqs.strands.names[i_str], curseq);
    free(curseq);
    free(cur_cons);
  }
  return ERR_OK;
}

int parse_design(char * filename, design_spec_t * spec)
{
  value_struc_t * val = NULL;
  char * read_file = NULL;
  int read_chunk =      100000;
  int max_file_size = 10000000;
  int curstart = 0;
  FILE * npfile = NULL;
  int bytes_read = 0;

  npfile = fopen(filename, "r");

  check(npfile != NULL, "Error reading file");
  read_file = (char *) malloc((read_chunk + 2) 
      * sizeof(char));
  bytes_read = fread(read_file, sizeof(char), read_chunk,
      npfile);

  while (!feof(npfile) && curstart < max_file_size && 
      bytes_read == read_chunk) {
    curstart += bytes_read;
    read_file = (char *) realloc(read_file, 
        (curstart + read_chunk + 2) * sizeof(char));
    check_mem(read_file);
    bytes_read = fread(read_file + curstart, 
        sizeof(char), read_chunk, npfile);
  }
  fclose(npfile);
  read_file[curstart + bytes_read] = '\n';
  read_file[curstart + bytes_read + 1] = '\0';
  init_design_spec(spec);
  
  val = get_ast(read_file);
  check(val != NULL, "Error parsing file");
  // np_tt_print_ast(val, 4);
  check(ERR_OK == process_design(spec, val),
      "Error processing design");
  check(ERR_OK == remap_nucleotides(spec),
      "Error remapping nucleotides");
  check(ERR_OK == make_off_targets(spec),
      "Error making off targets");

  np_tt_destroy_value_struc(val);

  print_strand_constraints(spec);

  free(read_file);
  return ERR_OK;
error:
  return ERR_INVALID_STATE;
}
