//#############################################################################
//#
//# Copyright 2008-2025, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################

#ifndef TEMPLATE_MODEL_H
#define TEMPLATE_MODEL_H

#include <memory>
#include <variant>
#include <string>
#include <sstream>
#include <map>
#include <vector>

namespace Loci {

  class TemplateValue ;
  class SimpleTemplateValue ;
  class DictionaryTemplateValue ;
  class ArrayTemplateValue ;

  typedef std::unique_ptr<TemplateValue> TemplateValueP ;

  typedef std::unique_ptr<TemplateValue> const_TemplateValueP ;

  class TemplateValue {
  public:
    typedef std::variant<int, double, std::string> value_t ;
    typedef std::map<std::string, TemplateValueP> dictionary_t ;
    typedef std::pair<std::string, TemplateValueP> dictionary_item_t ;
    typedef std::vector<TemplateValueP> array_t ;

  public:
    class const_Accessor ;

    class Accessor {
      friend class const_Accessor ;

    private:
      TemplateValueP & ptr_ ;

    public:
      Accessor(TemplateValueP & ptr) : ptr_(ptr) { }

      TemplateValue & operator=(int value) ;

      TemplateValue & operator=(double value) ;

      TemplateValue & operator=(std::string const & value) ;

      TemplateValue & operator=(dictionary_t const & value) ;

      TemplateValue & operator=(DictionaryTemplateValue const & value) ;

      TemplateValue & operator=(array_t const & value) ;

      TemplateValue & operator=(ArrayTemplateValue const & value) ;

      TemplateValueP clone() const ;

      std::string to_string() const ;

      bool is_value() const ;

      bool is_object() const ;

      bool is_array() const ;

      Accessor operator[](std::string const & key) ;

      const_Accessor operator[](std::string const & key) const ;

      Accessor operator[](size_t index) ;

      const_Accessor operator[](size_t index) const ;

      TemplateValue & insert(std::string const & key, int value) ;

      TemplateValue & insert(std::string const & key, double value) ;

      TemplateValue & insert(std::string const & key, std::string const & value) ;

      TemplateValue & insert(std::string const & key, DictionaryTemplateValue const & value) ;

      TemplateValue & insert(std::string const & key, ArrayTemplateValue const & value) ;

      TemplateValue & append(int value) ;

      TemplateValue & append(double value) ;

      TemplateValue & append(std::string const & value) ;

      TemplateValue & append(DictionaryTemplateValue const & value) ;

      TemplateValue & append(ArrayTemplateValue const & value) ;
    } ;

    class const_Accessor {
      friend class Accessor ;

    private:
      const_TemplateValueP const & ptr_ ;

    public:
      const_Accessor(const_TemplateValueP const & ptr) : ptr_(ptr) { }

      const_Accessor(Accessor const & src) : ptr_(src.ptr_) { }

      TemplateValueP clone() const ;

      std::string to_string() const ;

      bool is_value() const ;

      bool is_object() const ;

      bool is_array() const ;

      const_Accessor operator[](std::string const & key) const ;

      const_Accessor operator[](size_t index) const ;
    } ;

  public:
    TemplateValue() = default ;

    virtual ~TemplateValue() = default ;

    virtual std::unique_ptr<TemplateValue> clone() const = 0 ;

    virtual std::string to_string() const = 0 ;

    virtual bool is_value() const = 0 ;

    virtual bool is_object() const = 0 ;

    virtual bool is_array() const = 0 ;

    virtual Accessor operator[](std::string const & key) = 0 ;

    virtual const_Accessor operator[](std::string const & key) const = 0 ;

    virtual Accessor operator[](size_t index)  = 0 ;

    virtual const_Accessor operator[](size_t index) const = 0 ;

    virtual TemplateValue & insert(std::string const & key, int value) = 0 ;

    virtual TemplateValue & insert(std::string const & key, double value) = 0 ;

    virtual TemplateValue & insert(
      std::string const & key, std::string const & value
    ) = 0 ;

    virtual TemplateValue & insert(
      std::string const & key, DictionaryTemplateValue const & value
    ) = 0 ;

    virtual TemplateValue & insert(
      std::string const & key, ArrayTemplateValue const & value
    ) = 0 ;

    virtual TemplateValue & append(int value) = 0 ;

    virtual TemplateValue & append(double value) = 0 ;

    virtual TemplateValue & append(std::string const & value) = 0 ;

    virtual TemplateValue & append(DictionaryTemplateValue const & value) = 0 ;

    virtual TemplateValue & append(ArrayTemplateValue const & value) = 0 ;
  } ;

  std::ostream & operator<<(std::ostream & s, TemplateValue::Accessor const & obj) ;

  std::ostream & operator<<(std::ostream & s, TemplateValue::const_Accessor const & obj) ;

  std::ostream & operator<<(std::ostream & s, TemplateValue const & obj) ;

  TemplateValue const * lookup_path(
    TemplateValue const & root, std::string const & path
  ) ;

  bool is_truthy(TemplateValue const * value) ;

  class SimpleTemplateValue : public TemplateValue {
  public:
    typedef std::variant<int, double, std::string> value_t ;

  private:
    value_t value_ ;

  public:
    SimpleTemplateValue() = default ;

    virtual ~SimpleTemplateValue() = default ;

    SimpleTemplateValue(SimpleTemplateValue const & src) : value_(src.value_) {}

    SimpleTemplateValue & operator=(SimpleTemplateValue const & src) {
      value_ = src.value_ ;
      return *this ;
    }

    SimpleTemplateValue(int value) : value_(value) { }

    SimpleTemplateValue(double value) : value_(value) { }

    SimpleTemplateValue(std::string const & value) : value_(value) { }

    SimpleTemplateValue & operator=(int value) {
      value_ = value ;
      return *this ;
    }

    SimpleTemplateValue & operator=(double value) {
      value_ = value ;
      return *this ;
    }

    SimpleTemplateValue & operator=(std::string const & value) {
      value_ = value ;
      return *this ;
    }

    std::unique_ptr<TemplateValue> clone() const override ;

    std::string to_string() const override ;

    bool is_value() const override { return true ; }

    bool is_object() const override { return false ; }

    bool is_array() const override { return false ; }

    Accessor operator[](std::string const & key) override ;

    const_Accessor operator[](std::string const & key) const override ;

    Accessor operator[](size_t index) override ;

    const_Accessor operator[](size_t index) const override ;

    TemplateValue & insert(std::string const & key, int value) override ;

    TemplateValue & insert(std::string const & key, double value) override ;

    TemplateValue & insert(
      std::string const & key, std::string const & value
    ) override ;

    TemplateValue & insert(
      std::string const & key, DictionaryTemplateValue const & value
    ) override ;

    TemplateValue & insert(
      std::string const & key, ArrayTemplateValue const & value
    ) override ;

    TemplateValue & append(int value) override ;

    TemplateValue & append(double value) override ;

    TemplateValue & append(std::string const & value) override ;

    TemplateValue & append(DictionaryTemplateValue const & value) override ;

    TemplateValue & append(ArrayTemplateValue const & value) override ;
  } ;

  class DictionaryTemplateValue : public TemplateValue {
  private:
    dictionary_t members_ ;

  public:
    DictionaryTemplateValue() = default ;

    virtual ~DictionaryTemplateValue() = default ;

    DictionaryTemplateValue(DictionaryTemplateValue const & src) ;

    DictionaryTemplateValue(dictionary_t const & value) ;

    DictionaryTemplateValue & operator=(DictionaryTemplateValue const & src) ;

    DictionaryTemplateValue & operator=(dictionary_t const & value) ;

    TemplateValueP clone() const override ;

    std::string to_string() const override ;

    bool is_value() const override { return false ; }

    bool is_object() const override { return true ; }

    bool is_array() const override { return false ; }

    Accessor operator[](std::string const & key) override ;

    const_Accessor operator[](std::string const & key) const override ;

    Accessor operator[](size_t index) override ;

    const_Accessor operator[](size_t index) const override ;

    TemplateValue & insert(std::string const & key, int value) override ;

    TemplateValue & insert(std::string const & key, double value) override ;

    TemplateValue & insert(
      std::string const & key, std::string const & value
    ) override ;

    TemplateValue & insert(
      std::string const & key, DictionaryTemplateValue const & value
    ) override ;

    TemplateValue & insert(
      std::string const & key, ArrayTemplateValue const & value
    ) override ;

    TemplateValue & append(int value) ;

    TemplateValue & append(double value) ;

    TemplateValue & append(std::string const & value) ;

    TemplateValue & append(DictionaryTemplateValue const & value) ;

    TemplateValue & append(ArrayTemplateValue const & value) ;

    TemplateValue const * get(std::string const & key) const ;
  } ;

  class ArrayTemplateValue : public TemplateValue {
  private:
    array_t members_ ;

  public:
    ArrayTemplateValue() = default ;

    virtual ~ArrayTemplateValue() = default ;

    ArrayTemplateValue(ArrayTemplateValue const & src) ;

    ArrayTemplateValue(array_t const & value) ;

    ArrayTemplateValue & operator=(ArrayTemplateValue const & src) ;

    ArrayTemplateValue & operator=(array_t const & value) ;

    TemplateValueP clone() const override ;

    std::string to_string() const override ;

    bool is_value() const override { return false ; }

    bool is_object() const override { return false ; }

    bool is_array() const override { return true ; }

    size_t size() const ;

    TemplateValue const * get(size_t index) const ;

    Accessor operator[](std::string const & key) override ;

    const_Accessor operator[](std::string const & key) const override ;

    Accessor operator[](size_t index) override ;

    const_Accessor operator[](size_t index) const override ;

    TemplateValue & insert(std::string const & key, int value) override ;

    TemplateValue & insert(std::string const & key, double value) override ;

    TemplateValue & insert(
      std::string const & key, std::string const & value
    ) override ;

    TemplateValue & insert(
      std::string const & key, DictionaryTemplateValue const & value
    ) override ;

    TemplateValue & insert(
      std::string const & key, ArrayTemplateValue const & value
    ) override ;

    TemplateValue & append(int value) override ;

    TemplateValue & append(double value) override ;

    TemplateValue & append(std::string const & value) override ;

    TemplateValue & append(DictionaryTemplateValue const & value) override ;

    TemplateValue & append(ArrayTemplateValue const & value) override ;
  } ;

} // end: namespace Loci

#endif
