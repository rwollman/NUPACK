
#include <tclap/CmdLine.h>

#include "designer.h"
#include "design_spec.h"
#include "design_result.h"
#include "pathway_input.h"


#include <fstream>
#include <list>
#include <sys/time.h>

int main(int argc, char ** argv) {
  int rval = 0;
  int myrank = 0;


  SetExecutionPath(argc, argv);

  if (myrank == 0) {
    try {
      std::string filename;
      bool _drawjson = false;
      bool _prettyjson = false;
      bool _ppairsopt = false;
      bool _jsonopt = false;

      try {
        TCLAP::CmdLine cmds(
            "multitubedefect is the multi-tube defect calculation from NUPACK. "
            "For a consistent, fully-specified design, it calculates the defect "
            "in the multi-tube ensemble for the design."
            "See the user guide for further information", 
            ' ', NUPACK_VERSION);

        TCLAP::SwitchArg ppairsopt("p", "ppairs", 
            "print ppairs file for each structure", false);
        TCLAP::SwitchArg jsonopt("", "strucjson",
            "print json file for each structure", false);
        TCLAP::UnlabeledValueArg<std::string> file_arg("input", "script input file (.np)",
            true, "", 
            "input script");
#ifdef JSONCPP_FOUND
        TCLAP::SwitchArg drawjson("j", "drawjson", 
            "Output a json summary file for the specification", false);
        TCLAP::SwitchArg prettyjson("", "prettyjson", 
            "Output a pretty-printed json summary of the specification "
            "(overrides drawjson)", false);

        cmds.add(drawjson);
        cmds.add(prettyjson);
#endif
        cmds.add(ppairsopt);
        cmds.add(jsonopt);
        cmds.add(file_arg);

        cmds.parse(argc, argv);

        _ppairsopt = ppairsopt.getValue();
        _jsonopt = jsonopt.getValue();
#ifdef JSONCPP_FOUND
        _drawjson = drawjson.getValue();
        _prettyjson = prettyjson.getValue();
#endif
        filename = file_arg.getValue();
      } catch (TCLAP::ArgException &e) {
        NUPACK_ERROR(e.error() + " " + e.argId());
      }

      nupack::DesignSpec fullspec;
      
      auto & invars = fullspec.eval.options;
      invars.print_json = _jsonopt;
      invars.print_ppairs = _ppairsopt;
      // default seed to prevent diffs indicating seed is different (value is
      // arbitrary)
      invars.seed = 93;
      // avoid needing a "stop = ..." line for multitubedefect since there is
      // no optimization
      invars.add_global_stop = true;

      nupack::ScriptProcessor(filename, fullspec).parse_design();
      auto problem_domain = fullspec.eval.sequences.partially_specified_domain();
      NUPACK_CHECK(problem_domain == "", 
          "Unspecified nucleotides in domain: " + problem_domain 
          + "; cannot evaluate defect.");

      std::string output_prefix = fullspec.eval.options.file_prefix;

#ifdef JSONCPP_FOUND
      if (_drawjson || _prettyjson) {
        std::string json_filename = output_prefix + ".json";
        std::fstream jsonfile(json_filename, std::fstream::out | std::fstream::trunc);
        fullspec.serialize_json(jsonfile, _prettyjson);
      }
#endif

      timeval starttime;
      gettimeofday(&starttime, NULL);
      
      init_genrand(fullspec.eval.options.seed);
      
      nupack::Designer designer(fullspec);
      
      designer.evaluate_defect();
      const auto & res = designer.get_results();
      
      timeval endtime;
      gettimeofday(&endtime, NULL);
      DBL_TYPE elapsed = (DBL_TYPE) (endtime.tv_sec - starttime.tv_sec);
      elapsed += 1e-6 * (endtime.tv_usec - starttime.tv_usec);
      fullspec.eval.options.elapsed_time = elapsed;
      
      auto it = res.begin();
      int i = 0;
      for ( ; it != res.end(); ++it, ++i) {
        std::string full_output_prefix = output_prefix + "_" + nupack::to_string(i);
        std::string npo_filename = full_output_prefix + ".npo";
        
        std::fstream npofile(npo_filename, std::fstream::out | std::fstream::trunc);
        
        fullspec.eval.options.serialize(npofile, 0, "");
        it->eval.sequences.serialize(fullspec.eval.sequences, fullspec.eval.options, 
          npofile, 0, "");
        it->eval.physical.serialize(fullspec.eval.physical, fullspec.eval.options,
          npofile, 0, "");
        it->eval.physical.print_res_files(fullspec.eval.physical, fullspec.eval.options,
          full_output_prefix);
        fullspec.objectives.serialize(fullspec.eval, it->eval, npofile, 0, "");
      }
    } catch (nupack::NupackException & np_exc) {
      np_exc.print_message();
      rval = 1;
    }
  } else {
  }


  return rval;
}
