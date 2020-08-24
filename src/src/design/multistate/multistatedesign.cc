
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
      bool _nodesign = false;
      bool _drawjson = false;
      bool _prettyjson = false;
      bool _default_stops = false;
      bool _ppairsopt = false;
      bool _jsonopt = false;

      try {
        TCLAP::CmdLine cmds(
            "multitubedesign is the multi-tube design algorithm from NUPACK.  "
            "It allows linked optimization of test tube and structural "
            "ensemble defects under a variety of sequence constraints.  "
            "See the user guide for further information", 
            ' ', NUPACK_VERSION);

        TCLAP::SwitchArg nodesign("n", "nodesign", 
            "Don't design (just parse) and output any preliminary results", false);
        TCLAP::SwitchArg ppairsopt("p", "ppairs", 
            "print ppairs file for each structure", false);
        TCLAP::SwitchArg jsonopt("", "strucjson",
            "print json file for each structure", false);
        TCLAP::SwitchArg default_stops("", "addstops",
            "Automatically set default stop conditions for"
            " structures and tubes without a stop condition. "
            "(typically 1%)", false);
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
        cmds.add(nodesign);
        cmds.add(file_arg);
        cmds.add(default_stops);

        cmds.parse(argc, argv);

        _nodesign = nodesign.getValue();
        _ppairsopt = ppairsopt.getValue();
        _jsonopt = jsonopt.getValue();
#ifdef JSONCPP_FOUND
        _drawjson = drawjson.getValue();
        _prettyjson = prettyjson.getValue();
#endif
        _default_stops = default_stops.getValue();
        filename = file_arg.getValue();
      } catch (TCLAP::ArgException &e) {
        NUPACK_ERROR(e.error() + " " + e.argId());
      }

      nupack::DesignSpec fullspec;
      
      auto & invars = fullspec.eval.options;
      invars.add_global_stop = _default_stops;
      invars.print_json = _jsonopt;
      invars.print_ppairs = _ppairsopt;

      nupack::ScriptProcessor(filename, fullspec).parse_design();

      std::string output_prefix = invars.file_prefix;

#ifdef JSONCPP_FOUND
      if (_drawjson || _prettyjson) {
        std::string json_filename = output_prefix + ".json";
        std::fstream jsonfile(json_filename, std::fstream::out | std::fstream::trunc);
        fullspec.serialize_json(jsonfile, _prettyjson);
      }
#endif

      if (! _nodesign) {
        timeval starttime;
        gettimeofday(&starttime, NULL);
        
        NUPACK_DEBUG("seed: " + std::to_string(invars.seed));
        init_genrand(invars.seed);

        nupack::Designer maindes(fullspec);
        maindes.optimize();
        const auto & res = maindes.get_results();
        
        timeval endtime;
        gettimeofday(&endtime, NULL);
        DBL_TYPE elapsed = (DBL_TYPE) (endtime.tv_sec - starttime.tv_sec);
        elapsed += 1e-6 * (endtime.tv_usec - starttime.tv_usec);
        invars.elapsed_time = elapsed;

        auto it = res.begin();
        for ( ; it != res.end(); ++it) {
          auto i = it - res.begin();
          std::string full_output_prefix = output_prefix + "_" + nupack::to_string(i);
          std::string npo_filename = full_output_prefix + ".npo";

          std::fstream npofile(npo_filename, std::fstream::out | std::fstream::trunc);

          invars.serialize(npofile, 0, "");
          it->eval.sequences.serialize(fullspec.eval.sequences, invars, 
              npofile, 0, "");
          it->eval.physical.serialize(fullspec.eval.physical, invars,
              npofile, 0, "");
          it->eval.physical.print_res_files(fullspec.eval.physical, invars,
              full_output_prefix);
          fullspec.objectives.serialize(fullspec.eval, it->eval, npofile, 0, "");
        }
      }
    } catch (nupack::NupackException & np_exc) {
      np_exc.print_message();
      rval = 1;
    }
  } else {
  }


  return rval;
}
