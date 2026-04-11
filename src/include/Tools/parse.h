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
#ifndef PARSE_H
#define PARSE_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>


#include <istream>
#include <string>

namespace Loci {
namespace parse {

    /// @brief Advances an input stream past leading white space and line
    /// comments beginning with `//`.
    ///
    /// @param [s] input stream to be positioned at the next non-whitespace,
    /// non-comment character
    void kill_white_space(std::istream &s) ;

    /// @brief Checks whether the next token in the stream begins with a valid
    /// identifier lead character.
    ///
    /// @param [s] input stream to inspect
    /// @return true if the next non-whitespace character is alphabetic or
    /// underscore
    bool is_name(std::istream &s) ;

    /// @brief Extracts an identifier token from the input stream.
    ///
    /// @param [s] input stream to consume from
    /// @return identifier made of alphanumeric and underscore characters, or
    /// an empty string if the next token is not a valid name
    std::string get_name(std::istream &s) ;

    /// @brief Checks whether the next token in the stream can begin a valid
    /// integer literal.
    ///
    /// @param [s] input stream to inspect
    /// @return true if the next non-whitespace token begins with digits or
    /// with a sign followed by digits
    bool is_int(std::istream &s) ;

    /// @brief Extracts an integer value from the input stream.
    ///
    /// @param [s] input stream to consume from
    /// @return parsed integer value, or zero if the next token does not begin
    /// an integer literal
    long get_int(std::istream &s) ;

    /// @brief Checks whether the next token in the stream can begin a valid
    /// floating-point literal.
    ///
    /// @param [s] input stream to inspect
    /// @return true if the next non-whitespace token begins with digits, a
    /// decimal point, or a sign followed by digits or a decimal point
    bool is_real(std::istream &s) ;

    /// @brief Extracts a floating-point value from the input stream.
    ///
    /// @param [s] input stream to consume from
    /// @return parsed floating-point value, or zero if the next token does not
    /// begin a real literal
    double get_real(std::istream &s) ;

    /// @brief Checks whether the next token in the stream is a quoted string
    /// literal.
    ///
    /// @param [s] input stream to inspect
    /// @return true if the next non-whitespace character is a double quote
    bool is_string(std::istream &s) ;

    /// @brief Extracts the contents of a quoted string literal.
    ///
    /// @param [s] input stream to consume from
    /// @return characters between the opening and closing double quotes, or an
    /// empty string if the next token is not a string literal
    std::string get_string(std::istream &s) ;

    /// @brief Checks whether the next token in the stream matches a specific
    /// literal without consuming it.
    ///
    /// @param [s] input stream to inspect
    /// @param [token] literal token to compare against the input
    /// @return true if the next characters match token exactly
    bool is_token(std::istream &s, const std::string &token) ;

    /// @brief Consumes a specific literal token from the input stream if it
    /// matches exactly.
    ///
    /// @param [s] input stream to consume from
    /// @param [token] literal token expected at the current position
    /// @return true if token was matched and consumed, false otherwise
    bool get_token(std::istream &s, const std::string &token) ;
}
}

#endif
