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

#ifndef TEMPLATE_PARSER_H
#define TEMPLATE_PARSER_H

#include "template_ast.h"

#include <string>
#include <vector>
#include <memory>

namespace Loci {
  struct TemplateSourceLocation {
    size_t line ;
    size_t column ;

    TemplateSourceLocation() ;
    TemplateSourceLocation(size_t line, size_t column) ;
  } ;

  std::ostream & operator<<(
    std::ostream & s, TemplateSourceLocation const & loc
  ) ;

  struct TemplateToken {
  public:
    enum class Kind {
      Text,
      Space,
      Newline,
      BlockOpen,
      BlockClose,
      BlockIf,
      BlockElseIf,
      BlockElse,
      BlockEndIf,
      BlockEach,
      BlockEndEach,
      BlockPartial,
      InlineOpen,
      InlineClose,
      InlineIf,
      InlineElseIf,
      InlineElse,
      InlineEndIf,
      InlineEach,
      InlineEndEach,
      InlinePartial,
      End,
      Unknown
    } ;

    static std::string to_string(Kind kind) ;

  public:
    Kind kind ;
    std::string value ;
    TemplateSourceLocation location ;

  public:
    TemplateToken() ;
    TemplateToken(Kind kind, std::string value, size_t line, size_t column) ;
    bool is_text() const ;
    bool is_space() const ;
    bool is_newline() const ;
    bool is_block_open() const ;
    bool is_block_close() const ;
    bool is_block_if() const ;
    bool is_block_elseif() const ;
    bool is_block_else() const ;
    bool is_block_endif() const ;
    bool is_block_each() const ;
    bool is_block_endeach() const ;
    bool is_block_partial() const ;
    bool is_block_control() const ;
    bool is_inline_open() const ;
    bool is_inline_close() const ;
    bool is_inline_if() const ;
    bool is_inline_elseif() const ;
    bool is_inline_else() const ;
    bool is_inline_endif() const ;
    bool is_inline_each() const ;
    bool is_inline_endeach() const ;
    bool is_inline_partial() const ;
    bool is_inline_control() const ;
    bool is_end() const ;
  } ;

  std::ostream & operator<<(std::ostream & s, TemplateToken const & tok) ;

  struct TemplateLexerChar {
    char value ;
    size_t line ;
    size_t column ;
  } ;

  class TemplateLexerStream {
    std::string input_ ;
    size_t pos_ ;
    size_t line_ ;
    size_t column_ ;
    std::vector<TemplateLexerChar> buffer_ ;

  public:
    TemplateLexerStream(std::string const & input) ;
    TemplateLexerChar get() ;
    void unget(TemplateLexerChar const & c) ;
    TemplateLexerChar peek() ;
  } ;

  class TemplateLexer {
    TemplateLexerStream stream_ ;
    std::vector<TemplateToken> tokens_ ;

    TemplateToken next_token() ;
    TemplateToken scan_text(TemplateLexerChar c) ;
    TemplateToken scan_control(TemplateLexerChar c, bool block) ;

  public:
    TemplateLexer(std::string const & input) ;
    TemplateToken get() ;
    void unget(TemplateToken const & tok) ;
  };

  /** Grammar (matches TemplateParser / TemplateLexer)

      // Lexer tokens
      Text
      Space              // isspace, except CR/LF (those are Newline)
      Newline            // '\n' | '\r' | '\r' '\n'
      End
      BlockOpen          // "${"
      BlockClose         // "}$"
      InlineOpen         // "$["
      InlineClose        // "]$"
      BlockIfToken       // "${if "     (space after keyword required)
      BlockElseIfToken   // "${elseif "
      BlockElseToken     // "${else"    (no extra space required)
      BlockEndIfToken    // "${endif"
      BlockEachToken     // "${each "
      BlockEndEachToken  // "${endeach"
      BlockPartialToken  // "${> "
      InlineIfToken      // "$[if "
      InlineElseIfToken  // "$[elseif "
      InlineElseToken    // "$[else"
      InlineEndIfToken   // "$[endif"
      InlineEachToken    // "$[each "
      InlineEndEachToken // "$[endeach"
      InlinePartialToken // "$[> "

      Ws := (Space | Newline)*

      // parse()
      Document := Item*

      // parse_expression() — ordered choice. Leading Space* is absorbed into
      // a block form when the next token is a block starter (BlockOpen,
      // BlockIfToken, BlockEachToken, BlockPartialToken); otherwise Space is
      // an Item.
      Item := BlockVariable
            | BlockIf
            | BlockEach
            | BlockPartial
            | InlineVariable
            | InlineIf
            | InlineEach
            | InlinePartial
            | Text
            | Space
            | Newline

      BlockVariable := Space*
                       BlockOpen
                       Ws Text Ws
                       BlockClose
                       Space* (Newline | End)

      InlineVariable := InlineOpen
                        Ws Text Ws
                        InlineClose

      BlockIf := Space*
                 BlockIfToken
                 Ws Text Ws
                 BlockClose
                 Space* Newline
                 IfBody
                 (Space*
                  BlockElseIfToken
                  Ws Text Ws
                  BlockClose
                  Space* Newline
                  IfBody)*
                 (Space*
                  BlockElseToken
                  Ws
                  BlockClose
                  Space* Newline
                  IfBody)?
                 Space*
                 BlockEndIfToken
                 Ws
                 BlockClose
                 Space* (Newline | End)

      InlineIf := InlineIfToken
                  Ws Text Ws
                  InlineClose
                  IfBody
                  (InlineElseIfToken
                   Ws Text Ws
                   InlineClose
                   IfBody)*
                  (InlineElseToken
                   Ws
                   InlineClose
                   IfBody)?
                  InlineEndIfToken
                  Ws
                  InlineClose

      BlockEach := Space*
                   BlockEachToken
                   Ws Text Ws
                   BlockClose
                   Space* Newline
                   EachBody
                   Space*
                   BlockEndEachToken
                   Ws
                   BlockClose
                   Space* (Newline | End)

      InlineEach := InlineEachToken
                    Ws Text Ws
                    InlineClose
                    EachBody
                    InlineEndEachToken
                    Ws
                    InlineClose

      InlinePartial := InlinePartialToken
                       Ws Text Ws
                       InlineClose

      BlockPartial := Space*
                      BlockPartialToken
                      Ws Text Ws
                      BlockClose
                      Space* (Newline | End)

      // IfBody / EachBody are Item* that stop before a same-style closer
      // (inline vs block). Opposite-style constructs nest as Item.
      // For block bodies, Space* immediately before:
      //   - elseif / else / endif / endeach  → belongs to that clause (not the body)
      //   - a nested block starter (BlockOpen / BlockIf / BlockEach / BlockPartial)
      //     → absorbed like document-level Item (not emitted)
      // Other Space* in the body (e.g. indent before Text) remains an Item.
      IfBody := Item*
      EachBody := Item*

      // End is valid only on complete block forms (BlockVariable, BlockEndIf,
      // BlockEndEach, BlockPartial). End after if/elseif/else/each open
      // trailers is an error (those require Newline).

   */

  class TemplateParser {
    TemplateLexer lexer_ ;
    TemplateToken tok_ ;

    TemplateNodeP parse_expression() ;
    TemplateNodeP parse_inline_expression() ;
    void parse_inline_expression_close() ;
    void parse_block_expression_close() ;
    void eat_block_trailer(bool allow_end) ;
    bool is_block_starter() const ;
    TemplateNodeP parse_variable() ;
    bool is_if_terminator(bool block) const ;
    TemplateNodeP parse_if_body(bool block) ;
    TemplateNodeP parse_if_block(bool block) ;
    void parse_if_clause_close(bool block, bool allow_end) ;
    TemplateNodeP parse_if_conditional() ;
    bool is_each_terminator(bool block) const ;
    TemplateNodeP parse_each_body(bool block) ;
    TemplateNodeP parse_each(bool block) ;
    TemplateNodeP parse_partial(bool block) ;

  public:
    TemplateParser(std::string const & input) ;
    TemplateNodeList parse() ;
  } ;
}

#endif
