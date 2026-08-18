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

#include "template_parser.h"

#include <cctype>
#include <stdexcept>
#include <sstream>
#include <iostream>

namespace Loci {
  TemplateSourceLocation::TemplateSourceLocation() : line(-1), column(-1) { }

  TemplateSourceLocation::TemplateSourceLocation(size_t line, size_t column)
    : line(line), column(column) { }


  std::ostream & operator<<(
    std::ostream & s, TemplateSourceLocation const & loc
  ) {
    s << '(' << loc.line << ',' << loc.column << ')' ;
    return s ;
  }

  std::string TemplateToken::to_string(TemplateToken::Kind kind) {
    switch (kind) {
    case Kind::Text:
      return "Text" ;
    case Kind::Space:
      return "Space" ;
    case Kind::Newline:
      return "Newline" ;
    case Kind::BlockOpen:
      return "BlockOpen" ;
    case Kind::BlockClose:
      return "BlockClose" ;
    case Kind::BlockIf:
      return "BlockIf" ;
    case Kind::BlockElseIf:
      return "BlockElseIf" ;
    case Kind::BlockElse:
      return "BlockElse" ;
    case Kind::BlockEndIf:
      return "BlockEndIf" ;
    case Kind::BlockEach:
      return "BlockEach" ;
    case Kind::BlockEndEach:
      return "BlockEndEach" ;
    case Kind::BlockPartial:
      return "BlockPartial" ;
    case Kind::InlineOpen:
      return "InlineOpen" ;
    case Kind::InlineClose:
      return "InlineClose" ;
    case Kind::InlineIf:
      return "InlineIf" ;
    case Kind::InlineElseIf:
      return "InlineElseIf" ;
    case Kind::InlineElse:
      return "InlineElse" ;
    case Kind::InlineEndIf:
      return "InlineEndIf" ;
    case Kind::InlineEach:
      return "InlineEach" ;
    case Kind::InlineEndEach:
      return "InlineEndEach" ;
    case Kind::InlinePartial:
      return "InlinePartial" ;
    case Kind::End:
      return "End" ;
    }
    return "Unknown" ;
  }

  TemplateToken::TemplateToken() : kind(Kind::Unknown) { }

  TemplateToken::TemplateToken(
    TemplateToken::Kind kind, std::string value, size_t line, size_t column
  ) : kind(kind), value(value), location{line, column} { }

  bool TemplateToken::is_text() const {
    return kind == Kind::Text ;
  }

  bool TemplateToken::is_space() const {
    return kind == Kind::Space ;
  }

  bool TemplateToken::is_newline() const {
    return kind == Kind::Newline ;
  }

  bool TemplateToken::is_block_open() const {
    return kind == Kind::BlockOpen ;
  }

  bool TemplateToken::is_block_close() const {
    return kind == Kind::BlockClose ;
  }

  bool TemplateToken::is_block_if() const {
    return kind == Kind::BlockIf ;
  }

  bool TemplateToken::is_block_elseif() const {
    return kind == Kind::BlockElseIf ;
  }

  bool TemplateToken::is_block_else() const {
    return kind == Kind::BlockElse ;
  }

  bool TemplateToken::is_block_endif() const {
    return kind == Kind::BlockEndIf ;
  }

  bool TemplateToken::is_block_each() const {
    return kind == Kind::BlockEach ;
  }

  bool TemplateToken::is_block_endeach() const {
    return kind == Kind::BlockEndEach ;
  }

  bool TemplateToken::is_block_partial() const {
    return kind == Kind::BlockPartial ;
  }

  bool TemplateToken::is_block_control() const {
    return kind == Kind::BlockOpen ||
      kind == Kind::BlockIf ||
      kind == Kind::BlockElseIf ||
      kind == Kind::BlockElse ||
      kind == Kind::BlockEndIf ||
      kind == Kind::BlockEach ||
      kind == Kind::BlockEndEach ||
      kind == Kind::BlockPartial ;
  }

  bool TemplateToken::is_inline_open() const {
    return kind == Kind::InlineOpen ;
  }

  bool TemplateToken::is_inline_close() const {
    return kind == Kind::InlineClose ;
  }

  bool TemplateToken::is_inline_if() const {
    return kind == Kind::InlineIf ;
  }

  bool TemplateToken::is_inline_elseif() const {
    return kind == Kind::InlineElseIf ;
  }

  bool TemplateToken::is_inline_else() const {
    return kind == Kind::InlineElse ;
  }

  bool TemplateToken::is_inline_endif() const {
    return kind == Kind::InlineEndIf ;
  }

  bool TemplateToken::is_inline_each() const {
    return kind == Kind::InlineEach ;
  }

  bool TemplateToken::is_inline_endeach() const {
    return kind == Kind::InlineEndEach ;
  }

  bool TemplateToken::is_inline_partial() const {
    return kind == Kind::InlinePartial ;
  }

  bool TemplateToken::is_inline_control() const {
    return kind == Kind::InlineOpen ||
      kind == Kind::InlineIf ||
      kind == Kind::InlineElseIf ||
      kind == Kind::InlineElse ||
      kind == Kind::InlineEndIf ||
      kind == Kind::InlineEach ||
      kind == Kind::InlineEndEach ||
      kind == Kind::InlinePartial ;
  }

  bool TemplateToken::is_end() const {
    return kind == Kind::End ;
  }

  std::ostream & operator<<(std::ostream & s, TemplateToken const & tok) {
    s << '[' << tok.TemplateToken::to_string(tok.kind)
      << ' ' << tok.value << '@' << tok.location << ']' ;
    return s ;
  }


  TemplateLexerStream::TemplateLexerStream(std::string const & input)
    : input_(input), pos_(0), line_(0), column_(0) { }

  TemplateLexerChar TemplateLexerStream::get() {
    if(buffer_.empty()) {
      if(pos_ >= input_.size()) {
        char value = '\0' ;
        buffer_.emplace_back(TemplateLexerChar{value, line_, column_}) ;
      } else {
        char value = input_[pos_] ;
        buffer_.emplace_back(TemplateLexerChar{value, line_, column_}) ;
        ++pos_ ;
        if(value == '\n') {
          ++line_ ;
          column_ = 0 ;
        } else {
          ++column_ ;
        }
      }
    }

    TemplateLexerChar c = buffer_.back() ;
    buffer_.pop_back() ;
    return c ;
  }

  void TemplateLexerStream::unget(TemplateLexerChar const & c) {
    if(c.value != '\0') {
      buffer_.emplace_back(c) ;
    }
  }

  TemplateLexerChar TemplateLexerStream::peek() {
    TemplateLexerChar c = get() ;
    unget(c) ;
    return c ;
  }

  TemplateLexer::TemplateLexer(std::string const & input) : stream_(input) { }

  TemplateToken TemplateLexer::get() {
    if(tokens_.empty()) {
      tokens_.emplace_back(next_token()) ;
    }
    TemplateToken tok = tokens_.back() ;
    tokens_.pop_back() ;
    return tok ;
  }

  void TemplateLexer::unget(TemplateToken const & tok) {
    if(tok.kind != TemplateToken::Kind::End) {
      tokens_.emplace_back(tok) ;
    }
  }

  TemplateToken TemplateLexer::next_token() {
    TemplateLexerChar c = stream_.get() ;

    if(c.value == '\0') {
      return TemplateToken(
        TemplateToken::Kind::End, "", c.line, c.column
      ) ;
    }

    if(c.value == '\r' || c.value == '\n') {
      if(c.value == '\r') {
        TemplateLexerChar n = stream_.peek() ;
        if(n.value == '\n') {
          stream_.get() ;
        }
      }
      return TemplateToken(
        TemplateToken::Kind::Newline, "\n", c.line, c.column
      ) ;
    }

    if(std::isspace(c.value)) {
      TemplateToken tok ;
      tok.kind = TemplateToken::Kind::Space ;
      tok.value += c.value ;
      tok.location.line = c.line ;
      tok.location.column = c.column ;
      c = stream_.get() ;
      while(c.value != '\0' && c.value != '\n' && c.value != '\r' &&
            std::isspace(c.value)) {
        tok.value += c.value ;
        c = stream_.get() ;
      }
      stream_.unget(c) ;
      return tok ;
    }

    if(c.value == '$') {
      TemplateLexerChar c1 = stream_.get() ;
      if(c1.value == '{') {
        return scan_control(c, true) ;
      } else if(c1.value == '[') {
        return scan_control(c, false) ;
      }
      stream_.unget(c1) ;
    }

    if(c.value == '}') {
      TemplateLexerChar c1 = stream_.get() ;
      if(c1.value == '$') {
        return TemplateToken(
          TemplateToken::Kind::BlockClose, "}$", c.line, c.column
        ) ;
      }
      stream_.unget(c1) ;
    }

    if(c.value == ']') {
      TemplateLexerChar c1 = stream_.get() ;
      if(c1.value == '$') {
        return TemplateToken(
          TemplateToken::Kind::InlineClose, "]$", c.line, c.column
        ) ;
      }
      stream_.unget(c1) ;
    }

    return scan_text(c) ;
  }

  TemplateToken TemplateLexer::scan_text(TemplateLexerChar c) {
    TemplateToken tok ;
    tok.kind = TemplateToken::Kind::Text ;
    tok.value += c.value ;
    tok.location.line = c.line ;
    tok.location.column = c.column ;
    c = stream_.get() ;
    while(c.value != '\0' && !std::isspace(c.value) &&
          c.value != '$' && c.value != '}' && c.value != ']') {
      tok.value += c.value ;
      c = stream_.get() ;
    }
    stream_.unget(c) ;
    return tok ;
  }

  TemplateToken TemplateLexer::scan_control(TemplateLexerChar c, bool block) {
    TemplateLexerChar c2 = stream_.get() ;
    if(c2.value == 'i') {
      TemplateLexerChar c3 = stream_.get() ;
      if(c3.value == 'f') {
        TemplateLexerChar c4 = stream_.get() ;
        if(std::isspace(c4.value)) {
          return TemplateToken(
            block ? TemplateToken::Kind::BlockIf : TemplateToken::Kind::InlineIf,
            block ? "${if " : "$[if ",
            c.line, c.column
          ) ;
        }
        stream_.unget(c4) ;
      }
      stream_.unget(c3) ;
    } else if(c2.value == 'e') {
      TemplateLexerChar c3 = stream_.get() ;
      if(c3.value == 'a') {
        TemplateLexerChar c4 = stream_.get() ;
        if(c4.value == 'c') {
          TemplateLexerChar c5 = stream_.get() ;
          if(c5.value == 'h') {
            TemplateLexerChar c6 = stream_.get() ;
            if(std::isspace(c6.value)) {
              return TemplateToken(
                block ? TemplateToken::Kind::BlockEach : TemplateToken::Kind::InlineEach,
                block ? "${each " : "$[each ",
                c.line, c.column
              ) ;
            }
            stream_.unget(c6) ;
          }
          stream_.unget(c5) ;
        }
        stream_.unget(c4) ;
      } else if(c3.value == 'n') {
        TemplateLexerChar c4 = stream_.get() ;
        if(c4.value == 'd') {
          TemplateLexerChar c5 = stream_.get() ;
          if(c5.value == 'i') {
            TemplateLexerChar c6 = stream_.get() ;
            if(c6.value == 'f') {
              return TemplateToken(
                block ? TemplateToken::Kind::BlockEndIf : TemplateToken::Kind::InlineEndIf,
                block ? "${endif" : "$[endif",
                c.line, c.column
              ) ;
            }
            stream_.unget(c6) ;
          } else if(c5.value == 'e') {
            TemplateLexerChar c6 = stream_.get() ;
            if(c6.value == 'a') {
              TemplateLexerChar c7 = stream_.get() ;
              if(c7.value == 'c') {
                TemplateLexerChar c8 = stream_.get() ;
                if(c8.value == 'h') {
                  return TemplateToken(
                    block ? TemplateToken::Kind::BlockEndEach : TemplateToken::Kind::InlineEndEach,
                    block ? "${endeach" : "$[endeach",
                    c.line, c.column
                  ) ;
                }
                stream_.unget(c8) ;
              }
              stream_.unget(c7) ;
            }
            stream_.unget(c6) ;
          }
          stream_.unget(c5) ;
        }
        stream_.unget(c4) ;
      } else if(c3.value == 'l') {
        TemplateLexerChar c4 = stream_.get() ;
        if(c4.value == 's') {
          TemplateLexerChar c5 = stream_.get() ;
          if(c5.value == 'e') {
            TemplateLexerChar c6 = stream_.get() ;
            if(c6.value == 'i') {
              TemplateLexerChar c7 = stream_.get() ;
              if(c7.value == 'f') {
                TemplateLexerChar c8 = stream_.get() ;
                if(std::isspace(c8.value)) {
                  return TemplateToken(
                    block ? TemplateToken::Kind::BlockElseIf : TemplateToken::Kind::InlineElseIf,
                    block ? "${elseif " : "$[elseif ",
                    c.line, c.column
                  ) ;
                }
                stream_.unget(c8) ;
              }
              stream_.unget(c7) ;
            }
            stream_.unget(c6) ;

            return TemplateToken(
              block ? TemplateToken::Kind::BlockElse : TemplateToken::Kind::InlineElse,
              block ? "${else" : "$[else",
              c.line, c.column
            ) ;
          }
          stream_.unget(c5) ;
        }
        stream_.unget(c4) ;
      }
      stream_.unget(c3) ;
    } else if(c2.value == '>') {
      TemplateLexerChar c3 = stream_.get() ;
      if(std::isspace(c3.value)) {
        return TemplateToken(
          block ? TemplateToken::Kind::BlockPartial : TemplateToken::Kind::InlinePartial,
          block ? "${> " : "$[> ",
          c.line, c.column
        ) ;
      }
      stream_.unget(c3) ;
    }
    stream_.unget(c2) ;

    return TemplateToken(
      block ? TemplateToken::Kind::BlockOpen : TemplateToken::Kind::InlineOpen,
      block ? "${" : "$[",
      c.line, c.column
    ) ;
  }


  TemplateParser::TemplateParser(std::string const & input) : lexer_(input) { }

  TemplateNodeList TemplateParser::parse() {
    TemplateNodeList nodes ;
    tok_ = lexer_.get() ;
    while(!tok_.is_end()) {
      nodes.emplace_back(parse_expression()) ;
    }
    return nodes ;
  }

  bool TemplateParser::is_block_starter() const {
    return tok_.is_block_open() ||
      tok_.is_block_if() ||
      tok_.is_block_each() ||
      tok_.is_block_partial() ;
  }

  void TemplateParser::parse_block_expression_close() {
    while(tok_.is_space() || tok_.is_newline()) {
      tok_ = lexer_.get() ;
    }

    if(tok_.is_block_close()) {
      tok_ = lexer_.get() ;
    } else {
      std::ostringstream ss ;
      ss << "unclosed block expression " << tok_ ;
      throw std::runtime_error(ss.str()) ;
    }
  }

  void TemplateParser::eat_block_trailer(bool allow_end) {
    while(tok_.is_space()) {
      tok_ = lexer_.get() ;
    }
    if(tok_.is_newline()) {
      tok_ = lexer_.get() ;
    } else if(allow_end && tok_.is_end()) {
      // Complete block form may end at EOF.
    } else {
      std::ostringstream ss ;
      ss << "block statement must end with newline, got " << tok_ ;
      throw std::runtime_error(ss.str()) ;
    }
  }

  void TemplateParser::parse_if_clause_close(bool block, bool allow_end) {
    if(block) {
      parse_block_expression_close() ;
      eat_block_trailer(allow_end) ;
    } else {
      parse_inline_expression_close() ;
    }
  }

  TemplateNodeP TemplateParser::parse_expression() {
    TemplateNodeP expr ;

    // Leading Space* before a block starter is consumed (not emitted).
    if(tok_.is_space()) {
      std::string spaces = tok_.value ;
      tok_ = lexer_.get() ;
      while(tok_.is_space()) {
        spaces += tok_.value ;
        tok_ = lexer_.get() ;
      }
      if(!is_block_starter()) {
        return std::make_unique<SpaceTemplateNode>(spaces) ;
      }
      // else fall through with tok_ on the block starter
    }

    if(tok_.is_inline_open()) {
      tok_ = lexer_.get() ;
      expr = parse_inline_expression() ;
    } else if(tok_.is_inline_if()) {
      // Lexer already folded "$[if " into InlineIf.
      tok_ = lexer_.get() ;
      expr = parse_if_block(false) ;
    } else if(tok_.is_block_open()) {
      tok_ = lexer_.get() ;
      expr = parse_variable() ;
      parse_block_expression_close() ;
      eat_block_trailer(true) ;
    } else if(tok_.is_block_if()) {
      // Lexer already folded "${if " into BlockIf.
      tok_ = lexer_.get() ;
      expr = parse_if_block(true) ;
    } else if(tok_.is_inline_each()) {
      tok_ = lexer_.get() ;
      expr = parse_each(false) ;
    } else if(tok_.is_block_each()) {
      tok_ = lexer_.get() ;
      expr = parse_each(true) ;
    } else if(tok_.is_inline_partial()) {
      tok_ = lexer_.get() ;
      expr = parse_partial(false) ;
    } else if(tok_.is_block_partial()) {
      tok_ = lexer_.get() ;
      expr = parse_partial(true) ;
    } else if(tok_.is_text()) {
      expr = std::make_unique<TextTemplateNode>(tok_.value) ;
      tok_ = lexer_.get() ;
    } else if(tok_.is_newline()) {
      expr = std::make_unique<NewlineTemplateNode>() ;
      tok_ = lexer_.get() ;
    } else {
      std::ostringstream ss ;
      ss << "unexpected token " << tok_ ;
      throw std::runtime_error(ss.str()) ;
    }
    return expr ;
  }

  TemplateNodeP TemplateParser::parse_inline_expression() {
    // Called after InlineOpen: only variables use the bare open token.
    // Keyword forms (if/each/...) arrive as dedicated tokens in parse_expression.
    while(tok_.is_space() || tok_.is_newline()) {
      tok_ = lexer_.get() ;
    }

    TemplateNodeP expr = parse_variable() ;
    parse_inline_expression_close() ;
    return expr ;
  }

  TemplateNodeP TemplateParser::parse_variable() {
    std::string variable_name ;

    while(tok_.is_space() || tok_.is_newline()) {
      tok_ = lexer_.get() ;
    }

    if(tok_.is_text()) {
      variable_name = tok_.value ;
      tok_ = lexer_.get() ;
    } else {
      std::ostringstream ss ;
      ss << "unexpected token instead of variable name " << tok_ ;
      throw std::runtime_error(ss.str()) ;
    }

    if(variable_name.empty()) {
      std::ostringstream ss ;
      ss << "empty variable name " << tok_ ;
      throw std::runtime_error(ss.str()) ;
    }

    return std::make_unique<VariableTemplateNode>(variable_name) ;
  }

  bool TemplateParser::is_if_terminator(bool block) const {
    if(block) {
      return tok_.is_block_elseif() ||
        tok_.is_block_else() ||
        tok_.is_block_endif() ;
    }
    return tok_.is_inline_elseif() ||
      tok_.is_inline_else() ||
      tok_.is_inline_endif() ;
  }

  TemplateNodeP TemplateParser::parse_if_body(bool block) {
    auto expr = std::make_unique<ExpressionTemplateNode>() ;
    while(!tok_.is_end() && !is_if_terminator(block)) {
      // For block if, Space* before elseif/else/endif belongs to the clause,
      // not the body — discard it and stop. Space* before a nested block
      // starter is absorbed (not emitted), same as document-level Item.
      if(block && tok_.is_space()) {
        std::string spaces = tok_.value ;
        tok_ = lexer_.get() ;
        while(tok_.is_space()) {
          spaces += tok_.value ;
          tok_ = lexer_.get() ;
        }
        if(is_if_terminator(block)) {
          break ;
        }
        if(is_block_starter()) {
          // indent before nested block tag — do not emit
        } else {
          expr->push_back(SpaceTemplateNode(spaces)) ;
          continue ;
        }
      }
      expr->push_back(*parse_expression()) ;
    }
    return expr ;
  }

  TemplateNodeP TemplateParser::parse_if_block(bool block) {
    // Called after InlineIf / BlockIf has already been consumed.
    IfTemplateNode if_node ;

    TemplateNodeP condition = parse_if_conditional() ;
    parse_if_clause_close(block, false) ;
    TemplateNodeP branch = parse_if_body(block) ;
    if_node.append_branch(condition, branch) ;

    while((block && tok_.is_block_elseif()) ||
          (!block && tok_.is_inline_elseif())) {
      tok_ = lexer_.get() ;
      condition = parse_if_conditional() ;
      parse_if_clause_close(block, false) ;
      branch = parse_if_body(block) ;
      if_node.append_branch(condition, branch) ;
    }

    if((block && tok_.is_block_else()) ||
       (!block && tok_.is_inline_else())) {
      tok_ = lexer_.get() ;
      parse_if_clause_close(block, false) ;
      TemplateNodeP fallback = parse_if_body(block) ;
      if_node.set_fallback(fallback) ;
    }

    if((block && tok_.is_block_endif()) ||
       (!block && tok_.is_inline_endif())) {
      tok_ = lexer_.get() ;
      parse_if_clause_close(block, true) ;
    } else {
      std::ostringstream ss ;
      ss << "expected endif, got " << tok_ ;
      throw std::runtime_error(ss.str()) ;
    }

    return if_node.clone() ;
  }

  void TemplateParser::parse_inline_expression_close() {
    while(tok_.is_space() || tok_.is_newline()) {
      tok_ = lexer_.get() ;
    }

    if(tok_.is_inline_close()) {
      tok_ = lexer_.get() ;
    } else {
      std::ostringstream ss ;
      ss << "unclosed inline expression " << tok_ ;
      throw std::runtime_error(ss.str()) ;
    }
  }

  TemplateNodeP TemplateParser::parse_if_conditional() {
    while(tok_.is_space() || tok_.is_newline()) {
      tok_ = lexer_.get() ;
    }

    TemplateNodeP condition = parse_variable() ;

    return condition ;
  }

  bool TemplateParser::is_each_terminator(bool block) const {
    return block ? tok_.is_block_endeach() : tok_.is_inline_endeach() ;
  }

  TemplateNodeP TemplateParser::parse_each_body(bool block) {
    auto expr = std::make_unique<ExpressionTemplateNode>() ;
    while(!tok_.is_end() && !is_each_terminator(block)) {
      // Space* before endeach belongs to the closer. Space* before a nested
      // block starter is absorbed (not emitted), same as document-level Item.
      if(block && tok_.is_space()) {
        std::string spaces = tok_.value ;
        tok_ = lexer_.get() ;
        while(tok_.is_space()) {
          spaces += tok_.value ;
          tok_ = lexer_.get() ;
        }
        if(is_each_terminator(block)) {
          break ;
        }
        if(is_block_starter()) {
          // indent before nested block tag — do not emit
        } else {
          expr->push_back(SpaceTemplateNode(spaces)) ;
          continue ;
        }
      }
      expr->push_back(*parse_expression()) ;
    }
    return expr ;
  }

  TemplateNodeP TemplateParser::parse_each(bool block) {
    // Called after InlineEach / BlockEach has already been consumed.
    std::string iterate ;

    while(tok_.is_space() || tok_.is_newline()) {
      tok_ = lexer_.get() ;
    }
    if(!tok_.is_text()) {
      std::ostringstream ss ;
      ss << "unexpected token instead of each path " << tok_ ;
      throw std::runtime_error(ss.str()) ;
    }
    iterate = tok_.value ;
    tok_ = lexer_.get() ;

    parse_if_clause_close(block, false) ;
    TemplateNodeP body = parse_each_body(block) ;

    if(is_each_terminator(block)) {
      tok_ = lexer_.get() ;
      parse_if_clause_close(block, true) ;
    } else {
      std::ostringstream ss ;
      ss << "expected endeach, got " << tok_ ;
      throw std::runtime_error(ss.str()) ;
    }

    return std::make_unique<EachTemplateNode>(iterate, std::move(body)) ;
  }

  TemplateNodeP TemplateParser::parse_partial(bool block) {
    // Called after InlinePartial / BlockPartial has already been consumed.
    while(tok_.is_space() || tok_.is_newline()) {
      tok_ = lexer_.get() ;
    }
    if(!tok_.is_text()) {
      std::ostringstream ss ;
      ss << "unexpected token instead of partial name " << tok_ ;
      throw std::runtime_error(ss.str()) ;
    }
    std::string name = tok_.value ;
    tok_ = lexer_.get() ;

    if(block) {
      parse_block_expression_close() ;
      eat_block_trailer(true) ;
    } else {
      parse_inline_expression_close() ;
    }

    return std::make_unique<PartialTemplateNode>(name) ;
  }

}
