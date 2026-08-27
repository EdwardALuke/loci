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

#include "template_model.h"

#include <iostream>
#include <stdexcept>

namespace Loci {

  TemplateValue const * lookup_path(
    TemplateValue const & root, std::string const & path
  ) {
    if(path.empty()) {
      return &root ;
    }

    TemplateValue const * cur = &root ;
    size_t start = 0 ;
    while(start <= path.size()) {
      size_t dot = path.find('.', start) ;
      std::string key = (dot == std::string::npos)
        ? path.substr(start)
        : path.substr(start, dot - start) ;

      if(key.empty() || !cur->is_object()) {
        return nullptr ;
      }

      auto const * dict = dynamic_cast<DictionaryTemplateValue const *>(cur) ;
      if(!dict) {
        return nullptr ;
      }

      cur = dict->get(key) ;
      if(!cur) {
        return nullptr ;
      }

      if(dot == std::string::npos) {
        break ;
      }
      start = dot + 1 ;
    }

    return cur ;
  }

  bool is_truthy(TemplateValue const * value) {
    if(!value) {
      return false ;
    }
    if(value->is_value()) {
      std::string s = value->to_string() ;
      return !(s.empty() || s == "0" || s == "false") ;
    }
    if(value->is_array()) {
      auto const * arr = dynamic_cast<ArrayTemplateValue const *>(value) ;
      return arr && arr->size() > 0 ;
    }
    if(value->is_object()) {
      return true ;
    }
    return false ;
  }


  TemplateValue & TemplateValue::Accessor::operator=(int value) {
    ptr_.reset(new SimpleTemplateValue(value)) ;
    return *ptr_ ;
  }

  TemplateValue & TemplateValue::Accessor::operator=(double value) {
    ptr_.reset(new SimpleTemplateValue(value)) ;
    return *ptr_ ;
  }

  TemplateValue & TemplateValue::Accessor::operator=(std::string const & value) {
    ptr_.reset(new SimpleTemplateValue(value)) ;
    return *ptr_ ;
  }

  TemplateValue & TemplateValue::Accessor::operator=(
    TemplateValue::dictionary_t const & value
  ) {
    ptr_.reset(new DictionaryTemplateValue(value)) ;
    return *ptr_ ;
  }

  TemplateValue & TemplateValue::Accessor::operator=(
    DictionaryTemplateValue const & value
  ) {
    ptr_.reset(new DictionaryTemplateValue(value)) ;
    return *ptr_ ;
  }

  TemplateValue & TemplateValue::Accessor::operator=(
    TemplateValue::array_t const & value
  ) {
    ptr_.reset(new ArrayTemplateValue(value)) ;
    return *ptr_ ;
  }

  TemplateValue & TemplateValue::Accessor::operator=(
    ArrayTemplateValue const & value
  ) {
    ptr_.reset(new ArrayTemplateValue(value)) ;
    return *ptr_ ;
  }

  TemplateValueP TemplateValue::Accessor::clone() const {
    if(!ptr_)
      throw std::runtime_error("uninitialized TemplateValue::Accessor") ;
    return ptr_->clone() ;
  }

  std::string TemplateValue::Accessor::to_string() const {
    if(!ptr_)
      throw std::runtime_error("uninitialized TemplateValue::Accessor") ;
    return ptr_->to_string() ;
  }

  bool TemplateValue::Accessor::is_value() const {
    if(!ptr_)
      throw std::runtime_error("uninitialized TemplateValue::Accessor") ;
    return ptr_->is_value() ;
  }

  bool TemplateValue::Accessor::is_object() const {
    if(!ptr_)
      throw std::runtime_error("uninitialized TemplateValue::Accessor") ;
    return ptr_->is_object() ;
  }

  bool TemplateValue::Accessor::is_array() const {
    if(!ptr_) {
      throw std::runtime_error("uninitialized TemplateValue::Accessor") ;
    }
    return ptr_->is_array() ;
  }

  TemplateValue::Accessor TemplateValue::Accessor::operator[](
    std::string const & key
  ) {
    if(!ptr_)
      throw std::runtime_error("uninitialized TemplateValue::Accessor") ;
    return ptr_->operator[](key) ;
  }

  TemplateValue::const_Accessor TemplateValue::Accessor::operator[](
    std::string const & key
  ) const {
    if(!ptr_)
      throw std::runtime_error("uninitialized TemplateValue::Accessor") ;
    return ptr_->operator[](key) ;
  }

  TemplateValue::Accessor TemplateValue::Accessor::operator[](size_t index) {
    if(!ptr_)
      throw std::runtime_error("uninitialized TemplateValue::Accessor") ;
    return ptr_->operator[](index) ;
  }

  TemplateValue::const_Accessor TemplateValue::Accessor::operator[](size_t index) const {
    if(!ptr_)
      throw std::runtime_error("uninitialized TemplateValue::Accessor") ;
    return ptr_->operator[](index) ;
  }

  TemplateValue & TemplateValue::Accessor::insert(
    std::string const & key, int value
  ) {
    if(!ptr_)
      throw std::runtime_error("uninitialized TemplateValue::Accessor") ;
    return ptr_->insert(key, value) ;
  }

  TemplateValue & TemplateValue::Accessor::insert(
    std::string const & key, double value
  ) {
    if(!ptr_)
      throw std::runtime_error("uninitialized TemplateValue::Accessor") ;
    return ptr_->insert(key, value) ;
  }

  TemplateValue & TemplateValue::Accessor::insert(
    std::string const & key, std::string const & value
  ) {
    if(!ptr_)
      throw std::runtime_error("uninitialized TemplateValue::Accessor") ;
    return ptr_->insert(key, value) ;
  }

  TemplateValue & TemplateValue::Accessor::insert(
    std::string const & key, DictionaryTemplateValue const & value
  ) {
    if(!ptr_)
      throw std::runtime_error("uninitialized TemplateValue::Accessor") ;
    return ptr_->insert(key, value) ;
  }

  TemplateValue & TemplateValue::Accessor::insert(
    std::string const & key, ArrayTemplateValue const & value
  ) {
    if(!ptr_)
      throw std::runtime_error("uninitialized TemplateValue::Accessor") ;
    return ptr_->insert(key, value) ;
  }

  TemplateValue & TemplateValue::Accessor::append(int value) {
    if(!ptr_)
      throw std::runtime_error("uninitialized TemplateValue::Accessor") ;
    return ptr_->append(value) ;
  }

  TemplateValue & TemplateValue::Accessor::append(double value) {
    if(!ptr_)
      throw std::runtime_error("uninitialized TemplateValue::Accessor") ;
    return ptr_->append(value) ;
  }

  TemplateValue & TemplateValue::Accessor::append(std::string const & value) {
    if(!ptr_)
      throw std::runtime_error("uninitialized TemplateValue::Accessor") ;
    return ptr_->append(value) ;
  }

  TemplateValue & TemplateValue::Accessor::append(
    DictionaryTemplateValue const & value
  ) {
    if(!ptr_)
      throw std::runtime_error("uninitialized TemplateValue::Accessor") ;
    return ptr_->append(value) ;
  }

  TemplateValue & TemplateValue::Accessor::append(ArrayTemplateValue const & value) {
    if(!ptr_)
      throw std::runtime_error("uninitialized TemplateValue::Accessor") ;
    return ptr_->append(value) ;
  }

  TemplateValueP TemplateValue::const_Accessor::clone() const {
    if(!ptr_)
      throw std::runtime_error("uninitialized TemplateValue::const_Accessor") ;
    return ptr_->clone() ;
  }

  std::string TemplateValue::const_Accessor::to_string() const {
    if(!ptr_)
      throw std::runtime_error("uninitialized TemplateValue::const_Accessor") ;
    return ptr_->to_string() ;
  }

  bool TemplateValue::const_Accessor::is_value() const {
    if(!ptr_)
      throw std::runtime_error("uninitialized TemplateValue::const_Accessor") ;
    return ptr_->is_value() ;
  }

  bool TemplateValue::const_Accessor::is_object() const {
    if(!ptr_)
      throw std::runtime_error("uninitialized TemplateValue::const_Accessor") ;
    return ptr_->is_object() ;
  }

  bool TemplateValue::const_Accessor::is_array() const {
    if(!ptr_)
      throw std::runtime_error("uninitialized TemplateValue::const_Accessor") ;
    return ptr_->is_array() ;
  }

  TemplateValue::const_Accessor TemplateValue::const_Accessor::operator[](
    std::string const & key
  ) const {
    if(!ptr_)
      throw std::runtime_error("uninitilized TemplateValue::const_Accessor") ;
    return ptr_->operator[](key) ;
  }

  TemplateValue::const_Accessor TemplateValue::const_Accessor::operator[](
    size_t index
  ) const {
    if(!ptr_)
      throw std::runtime_error("unintialized TemplateValue::const_Accessor") ;
    return ptr_->operator[](index) ;
  }

  std::ostream & operator<<(std::ostream & s, TemplateValue::Accessor const & obj) {
    s << obj.to_string() ;
    return s ;
  }

  std::ostream & operator<<(std::ostream & s, TemplateValue::const_Accessor const & obj) {
    s << obj.to_string() ;
    return s ;
  }

  std::ostream & operator<<(std::ostream & s, TemplateValue const & obj) {
    s << obj.to_string() ;
    return s ;
  }

  std::unique_ptr<TemplateValue> SimpleTemplateValue::clone() const {
    return std::make_unique<SimpleTemplateValue>(*this) ;
  }

  std::string SimpleTemplateValue::to_string() const {
    std::string s = std::visit([](const auto & v) {
      std::ostringstream ss ;
      ss << v ;
      return ss.str() ;
    }, value_) ;
    return s ;
  }

  TemplateValue::Accessor SimpleTemplateValue::operator[](std::string const & key) {
    throw std::runtime_error("SimpleTemplateValue::operator[](string) unsupported") ;
  }

  TemplateValue::const_Accessor SimpleTemplateValue::operator[](
    std::string const & key
  ) const {
    throw std::runtime_error("SimpleTemplateValue::operator[](string) unsupported") ;
  }

  TemplateValue::Accessor SimpleTemplateValue::operator[](size_t index) {
    throw std::runtime_error("SimpleTemplateValue::operator[](size_t) unsupported") ;
  }

  TemplateValue::const_Accessor SimpleTemplateValue::operator[](size_t index) const {
    throw std::runtime_error("SimpleTemplateValue::operator[](size_t) unsupported") ;
  }

  TemplateValue & SimpleTemplateValue::insert(
    std::string const & key, int value
  ) {
    throw std::runtime_error(
      "SimpleTemplateValue::insert(string, int) unsupported"
    ) ;
    return *this ;
  }

  TemplateValue & SimpleTemplateValue::insert(
    std::string const & key, double value
  ) {
    throw std::runtime_error(
      "SimpleTemplateValue::insert(string, double) unsupported"
    ) ;
    return *this ;
  }

  TemplateValue & SimpleTemplateValue::insert(
    std::string const & key, std::string const & value
  ) {
    throw std::runtime_error(
      "SimpleTemplateValue::insert(string, string) unsupported"
    ) ;
    return *this ;
  }

  TemplateValue & SimpleTemplateValue::insert(
    std::string const & key, DictionaryTemplateValue const & value
  ) {
    throw std::runtime_error(
      "SimpleTemplateValue::insert(string, DictionaryTemplateValue) unsupported"
    ) ;
    return *this ;
  }

  TemplateValue & SimpleTemplateValue::insert(
    std::string const & key, ArrayTemplateValue const & value
  ) {
    throw std::runtime_error(
      "SimpleTemplateValue::insert(string, ArrayTemplateValue) unsupported"
    ) ;
    return *this ;
  }

  TemplateValue & SimpleTemplateValue::append(int value) {
    throw std::runtime_error("SimpleTemplateValue::append(int) unsupported") ;
    return *this ;
  }

  TemplateValue & SimpleTemplateValue::append(double value) {
    throw std::runtime_error("SimpleTemplateValue::append(double) unsupported") ;
    return *this ;
  }
  TemplateValue & SimpleTemplateValue::append(std::string const & value) {
    throw std::runtime_error("SimpleTemplateValue::append(string) unsupported") ;
    return *this ;
  }

  TemplateValue & SimpleTemplateValue::append(DictionaryTemplateValue const & value) {
    throw std::runtime_error(
      "SimpleTemplateValue::append(DictionaryTemplateValue) unsupported"
    ) ;
    return *this ;
  }

  TemplateValue & SimpleTemplateValue::append(ArrayTemplateValue const & value) {
    throw std::runtime_error(
      "SimpleTemplateValue::append(ArrayTemplateValue) unsupported"
    ) ;
    return *this ;
  }

  DictionaryTemplateValue::DictionaryTemplateValue(
    DictionaryTemplateValue const & src
  ) {
    dictionary_t::const_iterator iter ;
    iter = src.members_.begin() ;
    while(iter != src.members_.end()) {
      members_.emplace(std::make_pair(iter->first, iter->second->clone())) ;
      ++iter ;
    }
  }

  DictionaryTemplateValue::DictionaryTemplateValue(dictionary_t const & value) {
    dictionary_t::const_iterator iter ;
    iter = value.begin() ;
    while(iter != value.end()) {
      members_.emplace(std::make_pair(iter->first, iter->second->clone())) ;
      ++iter ;
    }
  }

  DictionaryTemplateValue & DictionaryTemplateValue::operator=(
    DictionaryTemplateValue const & src
  ) {
    members_.clear() ;
    dictionary_t::const_iterator iter ;
    iter = src.members_.begin() ;
    while(iter != src.members_.end()) {
      members_.emplace(std::make_pair(iter->first, iter->second->clone())) ;
      ++iter ;
    }
    return *this ;
  }

  DictionaryTemplateValue & DictionaryTemplateValue::operator=(
    dictionary_t const & value
  ) {
    members_.clear() ;
    dictionary_t::const_iterator iter ;
    iter = value.begin() ;
    while(iter != value.end()) {
      members_.insert(std::make_pair(iter->first, iter->second->clone())) ;
      ++iter ;
    }
    return *this ;
  }

  std::unique_ptr<TemplateValue> DictionaryTemplateValue::clone() const {
    return std::make_unique<DictionaryTemplateValue>(*this) ;
  }

  std::string DictionaryTemplateValue::to_string() const {
    std::ostringstream ss ;
    ss << "{" ;
    dictionary_t::const_iterator iter ;
    iter = members_.begin() ;
    while(iter != members_.end()) {
      ss << iter->first << "=" << iter->second->to_string() ;
      if(std::next(iter) != members_.end())
        ss << "," ;
      ++iter ;
    }
    ss << "}" ;
    return ss.str() ;
  }

  TemplateValue::Accessor DictionaryTemplateValue::operator[](
    std::string const & key
  ) {
    size_t dotpos = key.find('.') ;
    if(dotpos != std::string::npos)
      throw std::runtime_error("DictionaryTemplateValue key cannot contain '.'") ;
    return TemplateValue::Accessor(members_[key]) ;
  }

  TemplateValue::const_Accessor DictionaryTemplateValue::operator[](
    std::string const & key
  ) const {
    return TemplateValue::const_Accessor(members_.at(key)) ;
  }

  TemplateValue::Accessor DictionaryTemplateValue::operator[](size_t index) {
    throw std::runtime_error(
      "DictionaryTemplateValue::operator[](size_t) unsupported"
    ) ;
  }

  TemplateValue::const_Accessor DictionaryTemplateValue::operator[](
    size_t index
  ) const {
    throw std::runtime_error(
      "DictionaryTemplateValue::operator[](size_t) unsupported"
    ) ;
  }

  TemplateValue & DictionaryTemplateValue::insert(
    std::string const & key, int value
  ) {
    members_.emplace(std::make_pair(
      key, new SimpleTemplateValue(value)
    )) ;
    return *this ;
  }

  TemplateValue & DictionaryTemplateValue::insert(
    std::string const & key, double value
  ) {
    members_.emplace(std::make_pair(
      key, new SimpleTemplateValue(value)
    )) ;
    return *this ;
  }

  TemplateValue & DictionaryTemplateValue::insert(
    std::string const & key, std::string const & value
  ) {
    members_.emplace(std::make_pair(
      key, new SimpleTemplateValue(value)
    )) ;
    return *this ;
  }

  TemplateValue & DictionaryTemplateValue::insert(
    std::string const & key, DictionaryTemplateValue const & value
  ) {
    members_.emplace(std::make_pair(key, value.clone())) ;
    return *this ;
  }

  TemplateValue & DictionaryTemplateValue::insert(
    std::string const & key, ArrayTemplateValue const & value
  ) {
    members_.emplace(std::make_pair(key, value.clone())) ;
    return *this ;
  }

  TemplateValue & DictionaryTemplateValue::append(int value) {
    throw std::runtime_error("DictionaryTemplateValue::append(int) unsupported") ;
    return *this ;
  }

  TemplateValue & DictionaryTemplateValue::append(double value) {
    throw std::runtime_error("DictionaryTemplateValue::append(double) unsupported") ;
    return *this ;
  }

  TemplateValue & DictionaryTemplateValue::append(std::string const & value) {
    throw std::runtime_error("DictionaryTemplateValue::append(std::string) unsupported") ;
    return *this ;
  }

  TemplateValue & DictionaryTemplateValue::append(DictionaryTemplateValue const & value) {
    throw std::runtime_error("DictionaryTemplateValue::append(DictionaryTemplateValue) unsupported") ;
    return *this ;
  }

  TemplateValue & DictionaryTemplateValue::append(ArrayTemplateValue const & value) {
    throw std::runtime_error("DictionaryTemplateValue::append(ArrayTemplateValue) unsupported") ;
    return *this ;
  }

  TemplateValue const * DictionaryTemplateValue::get(std::string const & key) const {
    dictionary_t::const_iterator iter = members_.find(key) ;
    if(iter == members_.end() || !iter->second) {
      return nullptr ;
    }
    return iter->second.get() ;
  }

  ArrayTemplateValue::ArrayTemplateValue(ArrayTemplateValue const & src) {
    array_t::const_iterator iter = src.members_.begin() ;
    while(iter != src.members_.end()) {
      members_.emplace_back((*iter)->clone()) ;
      ++iter ;
    }
  }

  ArrayTemplateValue::ArrayTemplateValue(TemplateValue::array_t const & value) {
    array_t::const_iterator iter = value.begin() ;
    while(iter != value.end()) {
      members_.emplace_back((*iter)->clone()) ;
      ++iter ;
    }
  }

  ArrayTemplateValue & ArrayTemplateValue::operator=(ArrayTemplateValue const & src) {
    members_.clear() ;
    array_t::const_iterator iter = src.members_.begin() ;
    while(iter != src.members_.end()) {
      members_.emplace_back((*iter)->clone()) ;
      ++iter ;
    }
    return *this ;
  }

  ArrayTemplateValue & ArrayTemplateValue::operator=(TemplateValue::array_t const & value) {
    members_.clear() ;
    array_t::const_iterator iter = value.begin() ;
    while(iter != value.end()) {
      members_.emplace_back((*iter)->clone()) ;
      ++iter ;
    }
    return *this ;
  }

  TemplateValueP ArrayTemplateValue::clone() const {
    return std::make_unique<ArrayTemplateValue>(*this) ;
  }

  size_t ArrayTemplateValue::size() const {
    return members_.size() ;
  }

  TemplateValue const * ArrayTemplateValue::get(size_t index) const {
    if(index >= members_.size() || !members_[index]) {
      return nullptr ;
    }
    return members_[index].get() ;
  }

  std::string ArrayTemplateValue::to_string() const {
    std::ostringstream ss ;
    ss << "[" ;
    array_t::const_iterator iter = members_.begin() ;
    while(iter != members_.end()) {
      ss << *iter ;
      if(std::next(iter) != members_.end())
        ss << ',' ;
      ++iter ;
    }
    ss << "]" ;
    return ss.str() ;
  }

  TemplateValue::Accessor ArrayTemplateValue::operator[](
    std::string const & key
  ) {
    throw std::runtime_error(
      "ArrayTemplateValue::operator[](string) unsupported"
    ) ;
  }

  TemplateValue::const_Accessor ArrayTemplateValue::operator[](
    std::string const & key
  ) const {
    throw std::runtime_error(
      "ArrayTemplateValue::operator[](string) unsupported"
    ) ;
  }

  TemplateValue::Accessor ArrayTemplateValue::operator[](size_t index) {
    if(index >= members_.size()) {
      members_.resize(index+1) ;
    }
    return Accessor(members_[index]) ;
  }

  TemplateValue::const_Accessor ArrayTemplateValue::operator[](size_t index) const {
    return const_Accessor(members_[index]) ;
  }

  TemplateValue & ArrayTemplateValue::insert(std::string const & key, int value) {
    throw std::runtime_error("ArrayTemplateValue::insert(string, int) unsupported") ;
    return *this ;
  }

  TemplateValue & ArrayTemplateValue::insert(std::string const & key, double value) {
    throw std::runtime_error("ArrayTemplateValue::insert(string, double) unsupported") ;
    return *this ;
  }

  TemplateValue & ArrayTemplateValue::insert(
    std::string const & key, std::string const & value
  ) {
    throw std::runtime_error("ArrayTemplateValue::insert(string, string) unsupported") ;
    return *this ;
  }

  TemplateValue & ArrayTemplateValue::insert(
    std::string const & key, DictionaryTemplateValue const & value
  ) {
    throw std::runtime_error(
      "ArrayTemplateValue::insert(string, DictionaryTemplateValue) unsupported"
    ) ;
    return *this ;
  }

  TemplateValue & ArrayTemplateValue::insert(
    std::string const & key, ArrayTemplateValue const & value
  ) {
    throw std::runtime_error(
      "ArrayTemplateValue::insert(string, ArrayTemplateValue) unsupported"
    ) ;
    return *this ;
  }

  TemplateValue & ArrayTemplateValue::append(int value) {
    members_.emplace_back(new SimpleTemplateValue(value)) ;
    return *this ;
  }

  TemplateValue & ArrayTemplateValue::append(double value) {
    members_.emplace_back(new SimpleTemplateValue(value)) ;
    return *this ;
  }

  TemplateValue & ArrayTemplateValue::append(std::string const & value) {
    members_.emplace_back(new SimpleTemplateValue(value)) ;
    return *this ;
  }

  TemplateValue & ArrayTemplateValue::append(DictionaryTemplateValue const & value) {
    members_.emplace_back(value.clone()) ;
    return *this ;
  }

  TemplateValue & ArrayTemplateValue::append(ArrayTemplateValue const & value) {
    members_.emplace_back(value.clone()) ;
    return *this ;
  }
}

