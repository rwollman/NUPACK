#include "pathway_design.h"

int main (int argc, char** argv) {
  design_spec_t spec;
  result_t result;
  design_state_t states;
  result_t leafres;
  seqstate_t seqs;
  
  struct timeval starttime;
  struct timeval endtime;
  char * input_fn = NULL;
  char * output_fn = NULL;
  char * output_fn_leaves = NULL;
  int input_len;
  FILE * outf;
  FILE * leafoutf;
  char * input_prefix = NULL;
  int rank = 0;
  int depth = 0;
  int unspecified = 0;

  if (argc != 2) {
    fprintf(stderr, "Expected 1 argument. Received %i\n", argc - 1);
    fprintf(stderr, "Usage: %s <design_prefix>\n", argv[0]);
    goto error;
  }
  
  if (rank == 0) {
    input_len = strlen(argv[1]);

    input_prefix = (char *) malloc((input_len + 2) * sizeof(char));
    input_fn = (char *) malloc((input_len + 6) * sizeof(char));
    output_fn = (char *) malloc((input_len + 10) * sizeof(char));

    check_mem(input_prefix);
    check_mem(input_fn);
    check_mem(output_fn);
    strncpy(input_prefix, argv[1], input_len + 2);
    strncpy(input_fn, argv[1], input_len + 2);
    strncpy(output_fn, argv[1], input_len + 2);

    input_prefix[input_len] = '\0';

    strncat(input_fn, ".np", 6);
    strncat(output_fn, ".out", 10);

    check(ERR_OK == parse_design(input_fn, &spec), 
        "Error parsing file");

    spec.opts.seed = 0;

    gettimeofday(&starttime, NULL);

    spec.opts.file_prefix = input_prefix;
    spec.opts.designing = 0; // We aren't designing :(

    spec.opts.start_time = starttime.tv_sec + (1e-6 * starttime.tv_usec);
    check(ERR_OK == init_result(&result, &spec),
        "Error initializing result structure");

    init_decomposition(&spec);
    init_seqstate(&seqs, &spec);

    init_states(&states, &spec);

    get_decomposition(&states, depth, &spec);

    unspecified = get_n_mutable_nucleotides(&spec);
    check(unspecified == 0,
        "Defect can only be calculated over fully defined sequences, %i generalized nucleotide found",
        unspecified);

    // Ensure that the sequence is fully filled out, this should just copy
    // the sequence over for a fully determined sequence.
    check(ERR_OK == init_seqs_random(&seqs, &spec.seqs, &spec.opts),
        "Error initializing sequences");

    // Evaluate the defect
    evaluate_undesired(&result, &seqs, &spec);
    update_result(&result, &states, &seqs, &spec);

    gettimeofday(&endtime,NULL);

    result.elapsed_time = endtime.tv_sec - starttime.tv_sec 
      + 1e-6 * (endtime.tv_usec - starttime.tv_usec);

    outf = fopen(output_fn, "w");
    print_design_result(outf, &result, &seqs, 0, &spec, "tubedefect");
    fclose(outf);
    outf = NULL;

    if (spec.opts.print_leaves) {
      output_fn_leaves = (char *) malloc((input_len + 16) * sizeof(char));
      strncpy(output_fn_leaves, argv[1], input_len + 1);
      strncat(output_fn_leaves, "_leaves.npo", 16);

      leafoutf = fopen(output_fn_leaves, "w");

      init_states(&states, &spec);
      init_result(&leafres, &spec);

      get_decomposition(&states, 20, &spec);
      update_result(&leafres, &states, &seqs, &spec);

      print_design_result(leafoutf, &leafres, &seqs, 0, &spec, "tubedefect");

      free_result(&leafres);
      free_states(&states);
      fclose(leafoutf);
      leafoutf = NULL;
    }
    free_design_spec(&spec);
    free_seqstate(&seqs);
    free_result(&result);
  } else {
  }

  free(input_prefix);
  free(input_fn);
  free(output_fn);
  free(output_fn_leaves);

  return 0;
error:
  log_err("Error occurred in %s. Exiting", argv[0]);
  free(input_prefix);
  free(input_fn);
  free(output_fn);
  free(output_fn_leaves);
  return 2;
}
