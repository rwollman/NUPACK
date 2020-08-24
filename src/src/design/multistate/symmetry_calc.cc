#include "symmetry_calc.h"
#include "symmetry_objective.h"
#include "sequence_spec.h"
#include "weight_spec.h"
#include "design_debug.h"

namespace nupack {
const DBL_TYPE SingleSymmetrySpec::default_weight = 5;

SingleSymmetrySpec::SingleSymmetrySpec(const std::string & name,
    const std::vector<std::string> & domains) : 
    NamedSpec(name, "SSM"), names(domains), word_len{NUPACK_DEF_WORD_LEN}, 
    weights{default_weight}, weightscale(default_weight) {}

void SingleSymmetrySpec::resolve(const SequenceSpec & spec) {
  int i;
  for (i = 0; i < this->names.size(); i++) {
    const SingleSequenceSpec & cur_spec = spec.get_element(this->names[i]);
    const std::vector<int> & nuc_ids = cur_spec.get_nuc_ids();
    this->nuc_ids.push_back(nuc_ids);
  }
}

void SingleSymmetrySpec::set_word_len(std::vector<int> len) {
  this->word_len = len;
  int i;
  this->weights.clear();
  for (i = 0; i < this->word_len.size(); i++) {
    NUPACK_CHECK(this->word_len[i] > 0,
        "Word lengths in SSM objective must be > 0");
    int wl = this->word_len[i];
    this->weights.push_back(pow(this->weightscale, wl - this->word_len[0]));
  }
}

void SingleSymmetrySpec::set_weights(std::vector<DBL_TYPE> weights) {
  NUPACK_CHECK(weights.size() == this->word_len.size(),
      "number of weights must be equal to the number of word sizes");
  this->weights = weights;
}

void SingleSymmetrySpec::set_weights(DBL_TYPE weight) {
  this->weightscale = weight;
  int i;
  this->weights.clear();
  for (i = 0; i < this->word_len.size(); i++) {
    int wl = this->word_len[i];
    this->weights.push_back(pow(this->weightscale, wl - 1));
  }
}

void SingleSymmetryResult::evaluate(
    const SequenceState & seqs,
    const SingleSymmetrySpec & spec
    ) {
  const std::vector<std::vector<int> > & nuc_ids = spec.get_nuc_ids();

  DBL_TYPE total = 0;
  DBL_TYPE total_poss = 0;
  const std::vector<int> & word_len = spec.get_word_len();
  const std::vector<DBL_TYPE> & weights = spec.get_weights();
  int j;

  if (nuc_ids.size() != this->last_seqs.size()) {
    last_seqs.resize(nuc_ids.size());
  }

  std::vector<bool> eval(nuc_ids.size(), false);

  int i;
  for (i = 0; i < nuc_ids.size(); i++) {
    std::vector<int> curseq = seqs.get_sequence(nuc_ids[i]);
    if (curseq != last_seqs[i]) {
      eval[i] = true;
    }
  }

  this->nuc_defects.clear();

  for (i = 0; i < nuc_ids.size(); i++) {
    for (j = 0; j < nuc_ids[i].size(); j++) {
      this->nuc_defects[nuc_ids[i][j]] = 0;
    }
  }

  this->symmetries.resize(nuc_ids.size());

  NupackInvariants tmp_opts;
  
  for (i = 0; i < nuc_ids.size(); i++) {
    this->symmetries[i].resize(i + 1);
    for (j = 0; j <= i; j++) {
      std::vector<int> seq1 = seqs.get_sequence(nuc_ids[i]);
      std::vector<int> seq2 = seqs.get_sequence(nuc_ids[j]);

      std::vector<int> match_length(nuc_ids[j].size() + 1, 0);
      std::vector<int> old_match_length(nuc_ids[j].size() + 1, 0);
      int k, l, d;
      int cur_len;
      DBL_TYPE defect = 0;

      if (i != j) {
        for (k = 0; k < seq1.size(); k++) {
          for (l = 0; l < seq2.size(); l++) {
            if (seq1[k] == seq2[l]) {
              cur_len = old_match_length[l] + 1;
              match_length[l + 1] = cur_len;
              int i_wl;
              for (i_wl = 0; i_wl < word_len.size(); i_wl++) {
                int wl = word_len[i_wl];
                DBL_TYPE weight = weights[i_wl];
                if (cur_len >= wl) {
                  for (d = 0; d < wl; d++) {
                    this->nuc_defects[nuc_ids[i][k - d]] += weight / wl;
                    this->nuc_defects[nuc_ids[j][l - d]] += weight / wl;
                  }
                  defect += weight;
                }  else {
                  break;
                }
              }
            } else {
              match_length[l + 1] = 0;
            }
          }

          swap(match_length, old_match_length);
        }

        int i_wl, wl;
        for (i_wl = 0; i_wl < word_len.size(); i_wl++) {
          wl = word_len[i_wl];
          DBL_TYPE weight = weights[i_wl];
          int n_words1 = seq1.size() - wl + 1;
          int n_words2 = seq2.size() - wl + 1;
          total_poss += weight * n_words1 * n_words2;
        }
      } else {
        for (k = 0; k < seq1.size(); k++) {
          for (l = k + 1; l < seq2.size(); l++) {
            if (seq1[k] == seq2[l]) {
              cur_len = old_match_length[l] + 1;
              match_length[l + 1] = cur_len;
              int i_wl;
              for (i_wl = 0; i_wl < word_len.size(); i_wl++) {
                int wl = word_len[i_wl];
                DBL_TYPE weight = weights[i_wl];
                if (cur_len >= wl) {
                  for (d = 0; d < wl; d++) {
                    this->nuc_defects[nuc_ids[i][k - d]] += weight / wl;
                    this->nuc_defects[nuc_ids[j][l - d]] += weight / wl;
                  }
                  defect += weight;
                }  else {
                  break;
                }
              }
            } else {
              match_length[l + 1] = 0;
            }
          }
          swap(match_length, old_match_length);
        }
        int i_wl, wl;
        for (i_wl = 0; i_wl < word_len.size(); i_wl++) {
          wl = word_len[i_wl];
          DBL_TYPE weight = weights[i_wl];
          int n_words = seq1.size() - wl + 1;
          total_poss += weight * n_words * (n_words - 1) / 2;
        }
      }


      this->symmetries[i][j] = defect;
      total += defect;
    }
  }

  std::map<int, DBL_TYPE>::iterator it;
  for (it = this->nuc_defects.begin(); it != this->nuc_defects.end(); it++) {
    it->second /= total_poss;
  }

  this->defect = total / total_poss;
}

void SymmetryResult::evaluate(
    const SequenceState & seqs,
    const SymmetrySpec & spec 
  ) {
  int i;
  const std::vector<SingleSymmetrySpec> & specs = spec.get_specs();

  this->results.resize(specs.size());

  for (i = 0; i < results.size(); i++) {
    this->results[i].evaluate(seqs, specs[i]);
  }
}

DBL_TYPE SymmetryObjective::get_defect(const EvalSpec & spec, 
    const EvalResult & res, WeightSpec & weightspec) {
  if (this->id < 0) {
    int i;
    const std::vector<SingleSymmetrySpec> & specs = spec.symmetries.get_specs();
    for (i = 0; i < specs.size(); i++) {
      const std::string & curname = specs[i].get_name();
      if (curname == name) {
        this->id = i;
        break;
      }
    }
    NUPACK_CHECK(this->id >= 0, "Failed to find objective to link to");
  }

  const std::vector<SingleSymmetryResult> & results = res.symmetries.get_results();
  NUPACK_CHECK(this->id < results.size(),
      "Invalid ID or bad evaluation in SSM::get_defect. name: " + this->name + 
      " id: " +  to_string(this->id));

  return results[this->id].get_defect();
}

void SymmetryObjective::get_variable_mut_weights(const EvalSpec & spec,
    const EvalResult & res, WeightSpec & weightspec, std::vector<DBL_TYPE> & weights) {
  if (! this->satisfied(spec, res, weightspec)) {
    const std::vector<SingleSymmetryResult> & results = res.symmetries.get_results();
    const std::map<int, DBL_TYPE> & tmp = results[this->id].get_nuc_defects();

    std::map<int, DBL_TYPE>::const_iterator it;
    for (it = tmp.begin(); it != tmp.end(); it++) {
      NUPACK_CHECK(it->first >= 0 && it->first < weights.size(),
          "Invalid nucleotide id " + to_string(it->first));
      NUPACK_CHECK(!(it->second < 0),  
          "Invalid negative defect " + to_string(it->second));
      // NUPACK_DEBUG(it->first << " " << it->second);
      weights[it->first] += it->second;
    }
  }

  return;
}

bool SymmetryObjective::satisfied(const EvalSpec & spec, const EvalResult & res,
    WeightSpec & weightspec) {
  DBL_TYPE defect = this->get_defect(spec, res, weightspec);
  bool ret = true;
  if (defect > this->stop) {
    ret = false;
  }
  return ret;
}
}
