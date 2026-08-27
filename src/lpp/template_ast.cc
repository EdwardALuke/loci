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

#include "template_ast.h"
#include "template_engine.h"

#include <cctype>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace Loci {

  TemplateScope::TemplateScope(
    TemplateValue const & v, TemplateScope const * p
  ) : value(&v),
      parent(p),
      engine(p ? p->engine : nullptr),
      active_partial(nullptr),
      partial_depth(p ? p->partial_depth : 0),
      newline_handler(p ? p->newline_handler : nullptr) { }

  TemplateScope::TemplateScope(
    TemplateValue const & v, TemplateEngine const & eng
  ) : value(&v),
      parent(nullptr),
      engine(&eng),
      active_partial(nullptr),
      partial_depth(0),
      newline_handler(nullptr) { }

  TemplateScope::TemplateScope(
    TemplateScope const & outer, char const * partial_name
  ) : value(outer.value),
      parent(&outer),
      engine(outer.engine),
      active_partial(partial_name),
      partial_depth(outer.partial_depth + 1),
      newline_handler(outer.newline_handler) { }

  TemplateValue const * TemplateScope::lookup(std::string const & path) const {
    for(TemplateScope const * s = this ; s != nullptr ; s = s->parent) {
      if(!s->value) {
        continue ;
      }
      TemplateValue const * found = lookup_path(*s->value, path) ;
      if(found) {
        return found ;
      }
    }
    return nullptr ;
  }

  ExpressionTemplateNode::ExpressionTemplateNode() { }

  ExpressionTemplateNode::ExpressionTemplateNode(
    ExpressionTemplateNode const & src
  ) {
    TemplateNodeList::const_iterator iter = src.children_.begin() ;
    while(iter != src.children_.end()) {
      children_.emplace_back((*iter)->clone()) ;
      ++iter ;
    }
  }

  ExpressionTemplateNode::~ExpressionTemplateNode() { }

  ExpressionTemplateNode & ExpressionTemplateNode::operator=(
    ExpressionTemplateNode const & src
  ) {
    children_.clear() ;
    TemplateNodeList::const_iterator iter = src.children_.begin() ;
    while(iter != src.children_.end()) {
      children_.emplace_back((*iter)->clone()) ;
      ++iter ;
    }
    return *this ;
  }

  size_t ExpressionTemplateNode::size() const {
    return children_.size() ;
  }

  TemplateNodeList::iterator ExpressionTemplateNode::begin_children() {
    return children_.begin() ;
  }

  TemplateNodeList::const_iterator ExpressionTemplateNode::end_children() const {
    return children_.end() ;
  }

  void ExpressionTemplateNode::push_back(TemplateNode const & node) {
    children_.emplace_back(node.clone()) ;
  }

  TemplateNode::Kind ExpressionTemplateNode::kind() const {
    return Kind::Expression ;
  }

  TemplateNodeP ExpressionTemplateNode::clone() const {
    return std::make_unique<ExpressionTemplateNode>(*this) ;
  }

  std::ostream & ExpressionTemplateNode::evaluate(
    std::ostream & s, TemplateScope const & scope
  ) const {
    for(TemplateNodeList::const_iterator iter = children_.begin() ;
        iter != children_.end() ; ++iter) {
      (*iter)->evaluate(s, scope) ;
    }
    return s ;
  }


  TextTemplateNode::TextTemplateNode() { }

  TextTemplateNode::TextTemplateNode(std::string const & text)
    : text_(text) { }

  TextTemplateNode::TextTemplateNode(TextTemplateNode const & src)
    : text_(src.text_) { }

  TextTemplateNode::~TextTemplateNode() { }

  TextTemplateNode & TextTemplateNode::operator=(std::string const & text) {
    text_ = text ;
    return *this ;
  }

  TextTemplateNode & TextTemplateNode::operator=(TextTemplateNode const & src) {
    text_ = src.text_ ;
    return *this ;
  }

  std::string const & TextTemplateNode::get_text() const {
    return text_ ;
  }

  void TextTemplateNode::set_text(std::string const & text) {
    text_ = text ;
  }

  TemplateNode::Kind TextTemplateNode::kind() const {
    return Kind::Text ;
  }

  TemplateNodeP TextTemplateNode::clone() const {
    return std::make_unique<TextTemplateNode>(*this) ;
  }

  std::ostream & TextTemplateNode::evaluate(
    std::ostream & s, TemplateScope const & scope
  ) const {
    (void)scope ;
    s << text_ ;
    return s ;
  }


  NewlineTemplateNode::NewlineTemplateNode() { }

  NewlineTemplateNode::NewlineTemplateNode(NewlineTemplateNode const & src) { }

  NewlineTemplateNode::~NewlineTemplateNode() { }

  NewlineTemplateNode & NewlineTemplateNode::operator=(
    NewlineTemplateNode const & src) { return *this ; }

  TemplateNode::Kind NewlineTemplateNode::kind() const {
    return Kind::Newline ;
  }

  TemplateNodeP NewlineTemplateNode::clone() const {
    return std::make_unique<NewlineTemplateNode>(*this) ;
  }

  std::ostream & NewlineTemplateNode::evaluate(
    std::ostream & s, TemplateScope const & scope
  ) const {
    if(scope.newline_handler && *scope.newline_handler) {
      (*scope.newline_handler)(s, scope.active_partial) ;
      return s ;
    }
    s << std::endl ;
    return s ;
  }


  SpaceTemplateNode::SpaceTemplateNode() { }

  SpaceTemplateNode::SpaceTemplateNode(std::string const & text) {
    for(auto c : text) {
      if(std::isspace(c)) {
        text_ += c ;
      }
    }
  }

  SpaceTemplateNode::SpaceTemplateNode(SpaceTemplateNode const & src)
    : text_(src.text_) { }

  SpaceTemplateNode::~SpaceTemplateNode() { }

  SpaceTemplateNode & SpaceTemplateNode::operator=(std::string const & text) {
    text_.clear() ;
    for(auto c : text) {
      if(std::isspace(c)) {
        text_ += c ;
      }
    }
    return *this ;
  }

  SpaceTemplateNode & SpaceTemplateNode::operator=(SpaceTemplateNode const & src) {
    text_ = src.text_ ;
    return *this ;
  }

  TemplateNode::Kind SpaceTemplateNode::kind() const {
    return Kind::Space ;
  }

  TemplateNodeP SpaceTemplateNode::clone() const {
    return std::make_unique<SpaceTemplateNode>(*this) ;
  }

  std::ostream & SpaceTemplateNode::evaluate(
    std::ostream & s, TemplateScope const & scope
  ) const {
    (void)scope ;
    s << text_ ;
    return s ;
  }


  VariableTemplateNode::VariableTemplateNode() { }

  VariableTemplateNode::VariableTemplateNode(std::string const & name)
    : name_(name) { }

  VariableTemplateNode::VariableTemplateNode(VariableTemplateNode const & src)
    : name_(src.name_) { }

  VariableTemplateNode::~VariableTemplateNode() { }

  VariableTemplateNode & VariableTemplateNode::operator=(
    std::string const & name
  ) {
    name_ = name ;
    return *this ;
  }

  VariableTemplateNode & VariableTemplateNode::operator=(
    VariableTemplateNode const & src
  ) {
    name_ = src.name_ ;
    return *this ;
  }

  std::string const & VariableTemplateNode::get_name() const {
    return name_ ;
  }

  void VariableTemplateNode::set_name(std::string const & name) {
    name_ = name ;
  }

  TemplateNode::Kind VariableTemplateNode::kind() const {
    return Kind::Variable ;
  }

  TemplateNodeP VariableTemplateNode::clone() const {
    return std::make_unique<VariableTemplateNode>(*this) ;
  }

  std::ostream & VariableTemplateNode::evaluate(
    std::ostream & s, TemplateScope const & scope
  ) const {
    TemplateValue const * value = scope.lookup(name_) ;
    if(value) {
      s << value->to_string() ;
    }
    return s ;
  }


  IfTemplateNode::IfTemplateNode() { }

  IfTemplateNode::IfTemplateNode(IfTemplateNode const & src) {
    TemplateNodeList::const_iterator iter ;

    iter = src.conditions_.begin() ;
    while(iter != src.conditions_.end()) {
      conditions_.emplace_back((*iter)->clone()) ;
      ++iter ;
    }

    iter = src.branches_.begin() ;
    while(iter != src.branches_.end()) {
      branches_.emplace_back((*iter)->clone()) ;
      ++iter ;
    }

    if(src.fallback_) {
      fallback_ = src.fallback_->clone() ;
    } else {
      fallback_.reset() ;
    }
  }

  IfTemplateNode::~IfTemplateNode() { }

  IfTemplateNode & IfTemplateNode::operator=(IfTemplateNode const & src) {
    TemplateNodeList::const_iterator iter ;

    conditions_.clear() ;
    iter = src.conditions_.begin() ;
    while(iter != src.conditions_.end()) {
      conditions_.emplace_back((*iter)->clone()) ;
      ++iter ;
    }

    branches_.clear() ;
    iter = src.branches_.begin() ;
    while(iter != src.branches_.end()) {
      branches_.emplace_back((*iter)->clone()) ;
      ++iter ;
    }

    if(src.fallback_) {
      fallback_ = src.fallback_->clone() ;
    } else {
      fallback_.reset() ;
    }

    return *this ;
  }

  void IfTemplateNode::append_branch(
    TemplateNodeP const & condition,
    TemplateNodeP const & branch
  ) {
    conditions_.emplace_back(condition->clone()) ;
    branches_.emplace_back(branch->clone()) ;
  }

  void IfTemplateNode::set_fallback(TemplateNodeP const & fallback) {
    fallback_ = fallback->clone() ;
  }

  TemplateNode::Kind IfTemplateNode::kind() const {
    return Kind::If ;
  }

  TemplateNodeP IfTemplateNode::clone() const {
    return std::make_unique<IfTemplateNode>(*this) ;
  }

  std::ostream & IfTemplateNode::evaluate(
    std::ostream & s, TemplateScope const & scope
  ) const {
    TemplateNodeList::const_iterator citer = conditions_.begin() ;
    TemplateNodeList::const_iterator biter = branches_.begin() ;
    while(citer != conditions_.end() && biter != branches_.end()) {
      bool matched = false ;
      if((*citer)->kind() == Kind::Variable) {
        auto const * var =
          static_cast<VariableTemplateNode const *>(citer->get()) ;
        matched = is_truthy(scope.lookup(var->get_name())) ;
      } else {
        std::ostringstream ss ;
        (*citer)->evaluate(ss, scope) ;
        std::string text = ss.str() ;
        matched = !(text.empty() || text == "0" || text == "false") ;
      }

      if(matched) {
        (*biter)->evaluate(s, scope) ;
        return s ;
      }
      ++citer ;
      ++biter ;
    }

    if(fallback_) {
      fallback_->evaluate(s, scope) ;
    }
    return s ;
  }


  EachTemplateNode::EachTemplateNode() { }

  EachTemplateNode::EachTemplateNode(
    std::string const & iterate, TemplateNodeP body
  ) : iterate_(iterate), body_(std::move(body)) { }

  EachTemplateNode::EachTemplateNode(EachTemplateNode const & src)
    : iterate_(src.iterate_) {
    if(src.body_) {
      body_ = src.body_->clone() ;
    }
  }

  EachTemplateNode::~EachTemplateNode() { }

  EachTemplateNode & EachTemplateNode::operator=(EachTemplateNode const & src) {
    iterate_ = src.iterate_ ;
    if(src.body_) {
      body_ = src.body_->clone() ;
    } else {
      body_.reset() ;
    }
    return *this ;
  }

  std::string const & EachTemplateNode::get_iterate() const {
    return iterate_ ;
  }

  TemplateNode::Kind EachTemplateNode::kind() const {
    return Kind::Each ;
  }

  TemplateNodeP EachTemplateNode::clone() const {
    return std::make_unique<EachTemplateNode>(*this) ;
  }

  std::ostream & EachTemplateNode::evaluate(
    std::ostream & s, TemplateScope const & scope
  ) const {
    if(!body_) {
      return s ;
    }

    TemplateValue const * collection = scope.lookup(iterate_) ;
    if(!collection || !collection->is_array()) {
      return s ;
    }

    auto const * arr = dynamic_cast<ArrayTemplateValue const *>(collection) ;
    if(!arr) {
      return s ;
    }

    size_t const n = arr->size() ;
    for(size_t i = 0 ; i < n ; ++i) {
      TemplateValue const * item = arr->get(i) ;
      if(!item) {
        continue ;
      }

      DictionaryTemplateValue meta ;
      meta["@index"] = static_cast<int>(i) ;
      meta["@count"] = static_cast<int>(n) ;
      meta["@first"] = (i == 0) ? 1 : 0 ;
      meta["@last"] = (i + 1 == n) ? 1 : 0 ;

      TemplateScope item_scope(*item, &scope) ;
      TemplateScope child(meta, &item_scope) ;
      body_->evaluate(s, child) ;
    }
    return s ;
  }

  PartialTemplateNode::PartialTemplateNode() { }

  PartialTemplateNode::PartialTemplateNode(std::string const & name)
    : name_(name) { }

  PartialTemplateNode::PartialTemplateNode(PartialTemplateNode const & src)
    : name_(src.name_) { }

  PartialTemplateNode::~PartialTemplateNode() { }

  PartialTemplateNode & PartialTemplateNode::operator=(
    PartialTemplateNode const & src
  ) {
    name_ = src.name_ ;
    return *this ;
  }

  std::string const & PartialTemplateNode::get_name() const {
    return name_ ;
  }

  TemplateNode::Kind PartialTemplateNode::kind() const {
    return Kind::Partial ;
  }

  TemplateNodeP PartialTemplateNode::clone() const {
    return std::make_unique<PartialTemplateNode>(*this) ;
  }

  std::ostream & PartialTemplateNode::evaluate(
    std::ostream & s, TemplateScope const & scope
  ) const {
    if(!scope.engine) {
      throw std::runtime_error("partial evaluate without render engine") ;
    }

    for(TemplateScope const * fr = &scope ; fr != nullptr ; fr = fr->parent) {
      if(fr->active_partial && name_ == fr->active_partial) {
        throw std::runtime_error("recursive partial: " + name_) ;
      }
    }

    if(scope.partial_depth >= TemplateEngine::kMaxPartialDepth) {
      throw std::runtime_error("partial recursion depth exceeded") ;
    }

    TemplateNodeList const * nodes = scope.engine->find(name_) ;
    if(!nodes) {
      throw std::runtime_error("undefined partial: " + name_) ;
    }

    TemplateScope frame(scope, name_.c_str()) ;
    TemplateNodeList::const_iterator iter = nodes->begin() ;
    while(iter != nodes->end()) {
      (*iter)->evaluate(s, frame) ;
      ++iter ;
    }
    return s ;
  }
}
