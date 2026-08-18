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

#ifndef TEMPLATE_AST_H
#define TEMPLATE_AST_H

#include "template_model.h"

#include <functional>
#include <iosfwd>
#include <memory>
#include <vector>
#include <string>
#include <list>

namespace Loci {
  class TemplateNode ;
  typedef std::unique_ptr<TemplateNode> TemplateNodeP ;
  typedef std::list<TemplateNodeP> TemplateNodeList ;

  class ExpressionTemplateNode ;
  typedef std::unique_ptr<ExpressionTemplateNode> ExpressionTemplateNodeP ;

  class TextTemplateNode ;
  typedef std::unique_ptr<TextTemplateNode> TextTemplateNodeP ;

  class NewlineTemplateNode ;
  typedef std::unique_ptr<NewlineTemplateNode> NewlineTemplateNodeP ;

  class SpaceTemplateNode ;
  typedef std::unique_ptr<SpaceTemplateNode> SpaceTemplateNodeP ;

  class VariableTemplateNode ;
  typedef std::unique_ptr<VariableTemplateNode> VariableTemplateNodeP ;

  class IfTemplateNode ;
  typedef std::unique_ptr<IfTemplateNode> IfTemplateNodeP ;

  class EachTemplateNode ;
  typedef std::unique_ptr<EachTemplateNode> EachTemplateNodeP ;

  class PartialTemplateNode ;
  typedef std::unique_ptr<PartialTemplateNode> PartialTemplateNodeP ;

  class TemplateEngine ;

  using TemplateNewlineHandler =
    std::function<void(std::ostream & s, char const * active_partial)> ;

  struct TemplateScope {
    TemplateValue const * value ;
    TemplateScope const * parent ;
    TemplateEngine const * engine ;
    char const * active_partial ;
    int partial_depth ;
    TemplateNewlineHandler const * newline_handler ;

    TemplateScope(
      TemplateValue const & v, TemplateScope const * p = nullptr
    ) ;
    TemplateScope(
      TemplateValue const & v, TemplateEngine const & eng
    ) ;
    // Same data value as outer; marks this frame as entering partial `name`.
    TemplateScope(TemplateScope const & outer, char const * partial_name) ;

    TemplateValue const * lookup(std::string const & path) const ;
  } ;

  class TemplateNode {
  public:
    enum class Kind {
      Expression,
      Text,
      Newline,
      Space,
      Variable,
      If,
      Each,
      Partial
    } ;

  public:
    virtual ~TemplateNode() { }
    virtual Kind kind() const = 0 ;
    virtual TemplateNodeP clone() const = 0 ;
    virtual std::ostream & evaluate(
      std::ostream & s, TemplateScope const & scope
    ) const = 0 ;
  } ;

  class ExpressionTemplateNode : public TemplateNode {
    TemplateNodeList children_ ;
  public:
    ExpressionTemplateNode() ;
    ExpressionTemplateNode(ExpressionTemplateNode const & src) ;
    virtual ~ExpressionTemplateNode() ;
    ExpressionTemplateNode & operator=(ExpressionTemplateNode const & src) ;
    size_t size() const ;
    TemplateNodeList::iterator begin_children() ;
    TemplateNodeList::const_iterator end_children() const ;
    void push_back(TemplateNode const & node) ;
    Kind kind() const override ;
    TemplateNodeP clone() const override ;
    std::ostream & evaluate(
      std::ostream & s, TemplateScope const & scope
    ) const override ;
  } ;

  class TextTemplateNode : public TemplateNode {
    std::string text_ ;
  public:
    TextTemplateNode() ;
    TextTemplateNode(std::string const & text) ;
    TextTemplateNode(TextTemplateNode const & src) ;
    virtual ~TextTemplateNode() ;
    TextTemplateNode & operator=(std::string const & text) ;
    TextTemplateNode & operator=(TextTemplateNode const & src) ;
    std::string const & get_text() const ;
    void set_text(std::string const & text) ;
    Kind kind() const override ;
    TemplateNodeP clone() const override ;
    std::ostream & evaluate(
      std::ostream & s, TemplateScope const & scope
    ) const override ;
  } ;

  class NewlineTemplateNode : public TemplateNode {
  public:
    NewlineTemplateNode() ;
    NewlineTemplateNode(NewlineTemplateNode const & src) ;
    virtual ~NewlineTemplateNode() ;
    NewlineTemplateNode & operator=(NewlineTemplateNode const & src) ;
    Kind kind() const override ;
    TemplateNodeP clone() const override ;
    std::ostream & evaluate(
      std::ostream & s, TemplateScope const & scope
    ) const override ;
  } ;

  class SpaceTemplateNode : public TemplateNode {
    std::string text_ ;
  public:
    SpaceTemplateNode() ;
    SpaceTemplateNode(std::string const & text) ;
    SpaceTemplateNode(SpaceTemplateNode const & src) ;
    virtual ~SpaceTemplateNode() ;
    SpaceTemplateNode & operator=(std::string const & text) ;
    SpaceTemplateNode & operator=(SpaceTemplateNode const & src) ;
    Kind kind() const override ;
    TemplateNodeP clone() const override ;
    std::ostream & evaluate(
      std::ostream & s, TemplateScope const & scope
    ) const override ;
  } ;

  class VariableTemplateNode : public TemplateNode {
    std::string name_ ;
  public:
    VariableTemplateNode() ;
    VariableTemplateNode(std::string const & name) ;
    VariableTemplateNode(VariableTemplateNode const & src) ;
    virtual ~VariableTemplateNode() ;
    VariableTemplateNode & operator=(std::string const & name) ;
    VariableTemplateNode & operator=(VariableTemplateNode const & src) ;
    std::string const & get_name() const ;
    void set_name(std::string const & name) ;
    Kind kind() const override ;
    TemplateNodeP clone() const override ;
    std::ostream & evaluate(
      std::ostream & s, TemplateScope const & scope
    ) const override ;
  } ;

  class IfTemplateNode : public TemplateNode {
    TemplateNodeList conditions_ ;
    TemplateNodeList branches_ ;
    TemplateNodeP fallback_ ;

  public:
    IfTemplateNode() ;
    IfTemplateNode(IfTemplateNode const & src) ;
    virtual ~IfTemplateNode() ;
    IfTemplateNode & operator=(IfTemplateNode const & src) ;
    void append_branch(
      TemplateNodeP const & condition,
      TemplateNodeP const & branch
    ) ;
    void set_fallback(TemplateNodeP const & fallback) ;
    Kind kind() const override ;
    TemplateNodeP clone() const override ;
    std::ostream & evaluate(
      std::ostream & s, TemplateScope const & scope
    ) const override ;
  } ;

  class EachTemplateNode : public TemplateNode {
    std::string iterate_ ;
    TemplateNodeP body_ ;

  public:
    EachTemplateNode() ;
    EachTemplateNode(std::string const & iterate, TemplateNodeP body) ;
    EachTemplateNode(EachTemplateNode const & src) ;
    virtual ~EachTemplateNode() ;
    EachTemplateNode & operator=(EachTemplateNode const & src) ;
    std::string const & get_iterate() const ;
    Kind kind() const override ;
    TemplateNodeP clone() const override ;
    std::ostream & evaluate(
      std::ostream & s, TemplateScope const & scope
    ) const override ;
  } ;

  class PartialTemplateNode : public TemplateNode {
    std::string name_ ;

  public:
    PartialTemplateNode() ;
    PartialTemplateNode(std::string const & name) ;
    PartialTemplateNode(PartialTemplateNode const & src) ;
    virtual ~PartialTemplateNode() ;
    PartialTemplateNode & operator=(PartialTemplateNode const & src) ;
    std::string const & get_name() const ;
    Kind kind() const override ;
    TemplateNodeP clone() const override ;
    std::ostream & evaluate(
      std::ostream & s, TemplateScope const & scope
    ) const override ;
  } ;
}

#endif
