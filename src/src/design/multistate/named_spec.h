#pragma once

#include <string>

namespace nupack {
  class NamedSpec {
    public:
      NamedSpec() : name("--no-name--"), type("--no-type--"), id(-1) {}
      NamedSpec(std::string name, std::string type, int id = -1) :
          name(name), type(type), id(id) {}
      virtual ~NamedSpec() {};

      const std::string & get_name() const { return name; }
      const std::string & get_type() const { return type; }
      int get_id() const { return id; } 
      
      void set_name(std::string name) { this->name = name; } 
      void set_type(std::string type) { this->type = type; } 
      virtual void set_id(int id) { this->id = id; } 

    protected:
      void clone(const NamedSpec & other);
      std::string name;
      std::string type;
      int id;
  };

  inline void NamedSpec::clone(const NamedSpec & other) {
    this->name = other.name;
    this->type = other.type;
    this->id = other.id;
  }
  
  inline bool operator< (const NamedSpec & s1, const NamedSpec & s2) {
    return s1.get_name() < s2.get_name();
  }
  inline bool operator== (const NamedSpec & s1, const NamedSpec & s2) {
    return s1.get_name() == s2.get_name();
  }
  inline bool operator> (const NamedSpec & s1, const NamedSpec & s2) {
    return s1.get_name() > s2.get_name();
  }
  inline bool operator!= (const NamedSpec & s1, const NamedSpec & s2) {
    return s1.get_name() != s2.get_name();
  }
}