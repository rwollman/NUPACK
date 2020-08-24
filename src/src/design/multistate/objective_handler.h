#pragma once

#include "weight_spec.h"
#include "eval_spec.h"
#include "eval_result.h"

#include "types.h"

#include <memory>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <functional>

namespace nupack {
  // generic implementation 
  inline int get_nuc_id(const std::vector<int> & ind) { return ind.back(); }
  Map weight_defects(const Map, WeightSpec &, std::vector<int>);
  
  void get_variable_mut_weights(const std::vector<int> &, WeightSpec &, 
      std::vector<DBL_TYPE> &, std::vector<int>, const Map &, DBL_TYPE);
  
  class Objective {
    public:
      Objective(std::string name, std::string ob_type, DBL_TYPE stop) :
          name(name), objective_type(ob_type), stop(stop) {}
      virtual DBL_TYPE get_defect(const EvalSpec & spec,
          const EvalResult & res, WeightSpec & weightspec) = 0;

      virtual void get_variable_mut_weights(const EvalSpec & spec,
          const EvalResult & res, WeightSpec & weightspec,
          std::vector<DBL_TYPE> & weights) = 0;

      virtual std::shared_ptr<Objective> clone() const = 0;

      virtual int get_physical_id() const { return -1; };

      virtual std::vector<int> get_struc_id() const { return {}; }

      virtual std::vector<int> get_tube_id() const { return {}; }

      bool satisfied(const EvalSpec & spec, const EvalResult & res, 
          WeightSpec & weightspec) {
        return !(get_defect(spec, res, weightspec) > get_stop(spec, res));
      }

      virtual std::string get_objective_type(const EvalSpec & spec, const EvalResult & res) {
        return objective_type;
      }

      DBL_TYPE get_stop(const EvalSpec & spec, const EvalResult & res) const {
        return stop * pow(spec.options.f_stringent, res.physical.get_max_depth());
      }

      const std::string & get_name() const { return name; }

      virtual void serialize(const EvalSpec & spec, const EvalResult & res, 
          WeightSpec & weightspec, std::ostream & out, int indent=4, 
          std::string prefix = "") {
        std::string ind_string = prefix + std::string(indent, ' ');
        std::string name_ind = ind_string;
        if (indent > 2) name_ind = prefix + std::string(indent - 2, ' ') + "- ";
        
        out << name_ind << "name: " << name << std::endl;
        out << ind_string << "defect: " << get_defect(spec, res, weightspec) << std::endl;
        out << ind_string << "satisfied: " << 
            (satisfied(spec, res, weightspec) ? "true" : "false") << std::endl;
        out << ind_string << "stop: " << get_stop(spec, res) << std::endl;
        out << ind_string << "type: " << get_objective_type(spec, res) << std::endl;
      };

      virtual ~Objective() = default;

    protected:
      
      std::string name;
      std::string objective_type;
      DBL_TYPE stop;
  };

  template <typename ObjectiveType> 
  class Objective_CRTP : public Objective {
    public:
      Objective_CRTP(std::string name, std::string ob_type, DBL_TYPE stop) :
          Objective(name, ob_type, stop) {}
      virtual std::shared_ptr<Objective> clone() const {
        return std::make_shared<ObjectiveType>(static_cast<ObjectiveType const &>(*this));
      };
  };

  class TubeObjective : public Objective_CRTP<TubeObjective> {
    public:
      TubeObjective(std::string name, DBL_TYPE stop) :
          Objective_CRTP(name, "tube", stop), id(-1) {}
      ~TubeObjective() {}

      DBL_TYPE get_defect(const EvalSpec & spec, const EvalResult & res, 
            WeightSpec & weightspec);
      void get_variable_mut_weights(const EvalSpec & spec, const EvalResult & res, 
          WeightSpec & weightspec, std::vector<DBL_TYPE> & weights);

      std::vector<int> get_tube_id() const;
      int get_physical_id() const { return 0; }

      DBL_TYPE get_stop(const EvalSpec & spec, const EvalResult & res) const;

      bool satisfied(const EvalSpec & spec, const EvalResult & res,
          WeightSpec & weightspec);

      /**
       * \brief Weight the given map of defects according to
       *    the given weight specification. This is used to
       *    perform user-specified weighting of defects at the
       *    nucleotide level.
       */
      Map weight_defects(Map nuc_defects, WeightSpec & weights);

    private:
      int id;
      // Add in ignored nucleotide ids/domain names/strand names here.
      // Not sure precisely how to organize it.
  };

  class StructureObjective : public Objective_CRTP<StructureObjective> {
    public:
      StructureObjective(std::string name, DBL_TYPE stop) :
          Objective_CRTP(name, "structure", stop), id(-1), targ_id(-1) {}

      DBL_TYPE get_defect(const EvalSpec & spec,
          const EvalResult & res, WeightSpec & weightspec);
      DBL_TYPE get_stop(const EvalSpec & spec, const EvalResult & res) const;
      void get_variable_mut_weights(const EvalSpec & spec, 
          const EvalResult & res, WeightSpec & weightspec,
          std::vector<DBL_TYPE> & weights);

      std::vector<int> get_struc_id() const { return {this->id}; }
      
      int get_physical_id() const { return 0; }

      bool satisfied(const EvalSpec & spec, const EvalResult & res,
          WeightSpec & weightspec);

      Map weight_defects(Map nuc_defects, WeightSpec & weights);

    private:
      int id;
      int targ_id;
      // Add in ignored nucleotide ids/domain names/strand names here. 
      // Not sure how exactly to organize it
  };

  class CombinedObjective : public Objective_CRTP<CombinedObjective> {
    public:
      CombinedObjective(std::vector<std::string> names, DBL_TYPE stop) :
          Objective_CRTP("combined_objective", "combined", stop), element_names(names) {}

      DBL_TYPE get_defect(const EvalSpec & spec,
          const EvalResult & res, WeightSpec & weightspec);
      DBL_TYPE get_stop(const EvalSpec & spec, const EvalResult & res) const;
      void get_variable_mut_weights(const EvalSpec & spec, const EvalResult & res,
          WeightSpec & weightspec, std::vector<DBL_TYPE> & weights);

      int get_physical_id() const { return 0; }
      bool satisfied(const EvalSpec & spec, const EvalResult & res,
          WeightSpec & weightspec);

      std::vector<int> get_struc_id() const { return struc_ids; }
      std::vector<int> get_tube_id() const { return tube_ids; }

      std::vector<Map> weight_defects(std::vector<Map> nuc_defects, WeightSpec & weights);
      
      virtual void serialize(const EvalSpec & spec, const EvalResult & res, 
          WeightSpec & weightspec, std::ostream & out, int indent=4, 
          std::string prefix = "");

    private:
      void resolve_names(const EvalSpec & spec);

      std::vector<std::string> element_names;
      std::vector<int> struc_ids;
      std::vector<int> target_ids;
      std::vector<int> tube_ids;
  };

  class ObjectiveHandler {
    public:
      ObjectiveHandler() {}
      ObjectiveHandler(const ObjectiveHandler & other) { this->copy(other); }
      ObjectiveHandler & operator=(const ObjectiveHandler & other);
      ~ObjectiveHandler() {}

      void add_objective(const Objective & ob) { this->objectives.push_back(ob.clone()); }
      void set_weights(const WeightSpec & weights) { this->weights = weights; }

      std::vector<DBL_TYPE> get_defects(EvalSpec & spec, EvalResult & res);

      void get_variable_mut_weights(const EvalSpec & spec, 
          const EvalResult & res, std::vector<DBL_TYPE> & weights);

      void serialize(const EvalSpec & spec, const EvalResult & res,
          std::ostream & out, int indent=0, std::string prefix = "");

      void short_serial(const EvalSpec & spec, const EvalResult & res,
          std::ostream & out, int indent=0, std::string prefix = "");
      
      int get_physical_id(int ob_i) { return objectives.at(ob_i)->get_physical_id(); }
      std::vector<int> get_struc_id(int ob_i) { return objectives.at(ob_i)->get_struc_id(); }
      std::vector<int> get_tube_id(int ob_i) { return objectives.at(ob_i)->get_tube_id(); }

      bool all_satisfied(const EvalSpec & spec, const EvalResult & res);
      std::vector<bool> satisfied(const EvalSpec & spec, const EvalResult & res);

      std::vector<DBL_TYPE> get_stops(const EvalSpec & spec, const EvalResult & res);

      std::vector<std::string> get_names();

    private:
      void copy(const ObjectiveHandler & other);
      void clear() { objectives.clear(); }

      std::vector<std::shared_ptr<Objective>> objectives;
      WeightSpec weights;
  };
  
  
}
