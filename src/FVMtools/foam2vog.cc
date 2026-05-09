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

#include <Loci.h>
#include "vogtools.h"

#include <sys/stat.h>

#include <stdint.h>

#include <ctype.h>
#include <string.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using std::cerr ;
using std::cout ;
using std::endl ;
using std::ifstream ;
using std::ios ;
using std::istringstream ;
using std::pair ;
using std::string ;
using std::vector ;

// -----------------------------------------------------------------------
/// @brief Stores the OpenFOAM boundary metadata needed to build VOG
/// surface ids.
struct FoamPatchInfo {
  string name ;
  string type ;
  int nFaces ;
  int startFace ;
  int id ;
  // -----------------------------------------------------------------------
  /// @brief Initializes a patch record with unset indices and the default
  /// OpenFOAM patch type.
  FoamPatchInfo() {
    type = "patch" ;
    nFaces = -1 ;
    startFace = -1 ;
    id = -1 ;
  }
} ;

// -----------------------------------------------------------------------
/// @brief Stores one OpenFOAM cell zone and the cells assigned to it.
struct FoamCellZoneInfo {
  string name ;
  vector<int> cellLabels ;
} ;

// -----------------------------------------------------------------------
/// @brief Collects the mesh arrays read from an OpenFOAM polyMesh
/// directory.
struct FoamMeshData {
  vector<vector3d<double> > points ;
  vector<vector<int> > faces ;
  vector<int> owner ;
  vector<int> neighbour ;
  vector<FoamPatchInfo> patches ;
  vector<FoamCellZoneInfo> cellZones ;
} ;

// -----------------------------------------------------------------------
/// @brief Returns true when @p path exists and is a directory.
///
/// @param [path] Filesystem path to test.
/// @return True when @p path is a directory.
bool isDirectory(const string &path) {
  struct stat st ;
  if(stat(path.c_str(),&st) != 0)
    return false ;
  return S_ISDIR(st.st_mode) ;
}

// -----------------------------------------------------------------------
/// @brief Returns true when @p path exists.
///
/// @param [path] Filesystem path to test.
/// @return True when @p path exists.
bool fileExists(const string &path) {
  struct stat st ;
  return stat(path.c_str(),&st) == 0 ;
}

// -----------------------------------------------------------------------
/// @brief Removes trailing slashes while preserving the root directory.
///
/// @param [path] Input path that may include trailing '/' characters.
/// @return A normalized path without redundant trailing slashes.
string trimTrailingSlashes(string path) {
  while(path.size() > 1 && path[path.size()-1] == '/')
    path.erase(path.size()-1) ;
  return path ;
}

// -----------------------------------------------------------------------
/// @brief Returns the final path component of @p path.
///
/// @param [path] Input filesystem path.
/// @return Basename of @p path after trimming trailing slashes.
string baseName(string path) {
  path = trimTrailingSlashes(path) ;
  string::size_type pos = path.find_last_of('/') ;
  if(pos == string::npos)
    return path ;
  return path.substr(pos+1) ;
}

// -----------------------------------------------------------------------
/// @brief Returns the parent directory of @p path.
///
/// @param [path] Input filesystem path.
/// @return Parent directory, "." for relative single-component paths, and
/// "/" for root-owned children.
string parentPath(string path) {
  path = trimTrailingSlashes(path) ;
  string::size_type pos = path.find_last_of('/') ;
  if(pos == string::npos)
    return "." ;
  if(pos == 0)
    return "/" ;
  return path.substr(0,pos) ;
}

// -----------------------------------------------------------------------
/// @brief Concatenates two path components with one separating slash.
///
/// @param [lhs] Left-hand path component.
/// @param [rhs] Right-hand path component.
/// @return Joined filesystem path.
string joinPath(const string &lhs, const string &rhs) {
  if(lhs == "")
    return rhs ;
  if(lhs[lhs.size()-1] == '/')
    return lhs + rhs ;
  return lhs + "/" + rhs ;
}

// -----------------------------------------------------------------------
/// @brief Checks whether a directory contains the minimum polyMesh files
/// required by the converter.
///
/// @param [mesh_dir] Candidate polyMesh directory.
/// @return True when the required points, faces, owner, neighbour, and
/// boundary files exist.
bool hasPolyMeshFiles(const string &mesh_dir) {
  return fileExists(joinPath(mesh_dir,"points")) &&
         fileExists(joinPath(mesh_dir,"faces")) &&
         fileExists(joinPath(mesh_dir,"owner")) &&
         fileExists(joinPath(mesh_dir,"neighbour")) &&
         fileExists(joinPath(mesh_dir,"boundary")) ;
}

// -----------------------------------------------------------------------
/// @brief Resolves a user-supplied path into a concrete polyMesh
/// directory and output case name.
///
/// @param [input_name] Case directory, constant directory, or polyMesh
/// directory supplied on the command line.
/// @param [mesh_dir] Resolved polyMesh directory on success.
/// @param [case_name] Basename used to form the output `.vog` filename.
/// @return True when a valid polyMesh location is found.
bool resolveMeshDirectory(const string &input_name,
                          string &mesh_dir,
                          string &case_name) {
  string input = trimTrailingSlashes(input_name) ;
  string direct_poly = joinPath(input,"constant/polyMesh") ;
  string child_poly = joinPath(input,"polyMesh") ;

  if(isDirectory(input) && hasPolyMeshFiles(input)) {
    mesh_dir = input ;
    string mesh_base = baseName(input) ;
    if(mesh_base == "polyMesh") {
      case_name = baseName(parentPath(parentPath(input))) ;
    } else {
      case_name = mesh_base ;
    }
    return true ;
  }

  if(isDirectory(direct_poly) && hasPolyMeshFiles(direct_poly)) {
    mesh_dir = direct_poly ;
    case_name = baseName(input) ;
    return true ;
  }

  if(isDirectory(child_poly) && hasPolyMeshFiles(child_poly)) {
    mesh_dir = child_poly ;
    string input_base = baseName(input) ;
    if(input_base == "constant")
      case_name = baseName(parentPath(input)) ;
    else
      case_name = input_base ;
    return true ;
  }

  return false ;
}

// -----------------------------------------------------------------------
/// @brief Identifies punctuation tokens used by OpenFOAM dictionaries and
/// lists.
///
/// @param [c] Character being classified.
/// @return True when @p c should be emitted as a standalone token.
bool isDelimiter(char c) {
  return c == '{' || c == '}' || c == '(' || c == ')' ||
         c == ';' || c == '[' || c == ']' ;
}

// -----------------------------------------------------------------------
/// @brief Tokenizes an ASCII OpenFOAM file, skipping comments and keeping
/// structural delimiters as separate tokens.
///
/// @param [filename] OpenFOAM dictionary or list file to read.
/// @param [tokens] Parsed tokens appended by the tokenizer.
/// @return True when the file is read and tokenized successfully.
bool tokenizeFoamFile(const string &filename, vector<string> &tokens) {
  ifstream file(filename.c_str(),ios::in) ;
  if(file.fail()) {
    cerr << "unable to open '" << filename << "'" << endl ;
    return false ;
  }

  while(file.peek() != EOF) {
    int next = file.peek() ;
    if(isspace(next)) {
      file.get() ;
      continue ;
    }

    if(next == '/') {
      file.get() ;
      int c2 = file.peek() ;
      if(c2 == '/') {
        while(file.peek() != EOF && file.get() != '\n')
          ;
        continue ;
      }
      if(c2 == '*') {
        file.get() ;
        int prev = 0 ;
        bool closed = false ;
        while(file.peek() != EOF) {
          int c = file.get() ;
          if(prev == '*' && c == '/') {
            closed = true ;
            break ;
          }
          prev = c ;
        }
        if(!closed) {
          cerr << "unterminated block comment in '" << filename << "'" << endl ;
          return false ;
        }
        continue ;
      }
      tokens.push_back("/") ;
      continue ;
    }

    if(next == '"') {
      file.get() ;
      string tok ;
      while(file.peek() != EOF) {
        int c = file.get() ;
        if(c == '"')
          break ;
        tok += char(c) ;
      }
      tokens.push_back(tok) ;
      continue ;
    }

    if(isDelimiter(char(next))) {
      string tok ;
      tok += char(file.get()) ;
      tokens.push_back(tok) ;
      continue ;
    }

    string tok ;
    while(file.peek() != EOF) {
      int c = file.peek() ;
      if(isspace(c) || isDelimiter(char(c)))
        break ;
      if(c == '/') {
        file.get() ;
        int c2 = file.peek() ;
        file.unget() ;
        if(c2 == '/' || c2 == '*')
          break ;
      }
      tok += char(file.get()) ;
    }
    if(tok != "")
      tokens.push_back(tok) ;
  }

  return true ;
}

// -----------------------------------------------------------------------
/// @brief Converts an ASCII token to a signed integer.
///
/// @param [tok] Token text to convert.
/// @param [val] Parsed integer value on success.
/// @return True when the token is a valid base-10 integer.
bool tokenToInt(const string &tok, int &val) {
  char *endptr = 0 ;
  long tmp = strtol(tok.c_str(),&endptr,10) ;
  if(endptr == tok.c_str() || *endptr != '\0')
    return false ;
  val = int(tmp) ;
  return true ;
}

// -----------------------------------------------------------------------
/// @brief Converts an ASCII token to a floating-point value.
///
/// @param [tok] Token text to convert.
/// @param [val] Parsed floating-point value on success.
/// @return True when the token is a valid scalar.
bool tokenToDouble(const string &tok, double &val) {
  char *endptr = 0 ;
  val = strtod(tok.c_str(),&endptr) ;
  if(endptr == tok.c_str() || *endptr != '\0')
    return false ;
  return true ;
}

// -----------------------------------------------------------------------
/// @brief Reads an entire file into a byte buffer.
///
/// @param [filename] File to load.
/// @param [data] Raw file contents on success.
/// @return True when the file is read completely.
bool readFileBytes(const string &filename, vector<char> &data) {
  ifstream file(filename.c_str(),ios::in | ios::binary) ;
  if(file.fail()) {
    cerr << "unable to open '" << filename << "'" << endl ;
    return false ;
  }

  file.seekg(0,ios::end) ;
  std::streamoff sz = file.tellg() ;
  file.seekg(0,ios::beg) ;
  if(sz < 0) {
    cerr << "unable to determine file size for '" << filename << "'" << endl ;
    return false ;
  }

  data.resize(size_t(sz)) ;
  if(sz != 0)
    file.read(&data[0],sz) ;
  if(file.fail()) {
    cerr << "error reading '" << filename << "'" << endl ;
    return false ;
  }
  return true ;
}

// -----------------------------------------------------------------------
/// @brief Verifies that the next ASCII token matches an expected string.
///
/// @param [tokens] Token stream being parsed.
/// @param [idx] Current token index, advanced on success.
/// @param [expected] Required token text.
/// @param [filename] Source filename used in diagnostics.
/// @return True when the expected token is present.
bool expectToken(const vector<string> &tokens,
                 size_t &idx,
                 const string &expected,
                 const string &filename) {
  if(idx >= tokens.size()) {
    cerr << "unexpected end of file while reading '" << filename
         << "'" << endl ;
    return false ;
  }
  if(tokens[idx] != expected) {
    cerr << "expected '" << expected << "' in '" << filename
         << "', got '" << tokens[idx] << "'" << endl ;
    return false ;
  }
  ++idx ;
  return true ;
}

// -----------------------------------------------------------------------
/// @brief Verifies that the next raw byte matches an expected character.
///
/// @param [data] Raw file contents being parsed.
/// @param [idx] Current byte index, advanced on success.
/// @param [expected] Required byte value.
/// @param [filename] Source filename used in diagnostics.
/// @return True when the expected character is present.
bool expectChar(const vector<char> &data,
                size_t &idx,
                char expected,
                const string &filename) {
  if(idx >= data.size()) {
    cerr << "unexpected end of file while reading '" << filename
         << "'" << endl ;
    return false ;
  }
  if(data[idx] != expected) {
    cerr << "expected '" << expected << "' in '" << filename
         << "', got '" << data[idx] << "'" << endl ;
    return false ;
  }
  ++idx ;
  return true ;
}

// -----------------------------------------------------------------------
/// @brief Skips one ASCII OpenFOAM entry value, including nested list,
/// dictionary, and bracketed expressions.
///
/// @param [tokens] Token stream being parsed.
/// @param [idx] Current token index, advanced to the terminating ';'.
void skipEntryValue(const vector<string> &tokens, size_t &idx) {
  int nparen = 0 ;
  int nbrace = 0 ;
  int nbracket = 0 ;

  while(idx < tokens.size()) {
    string tok = tokens[idx++] ;
    if(tok == "(")
      ++nparen ;
    else if(tok == ")")
      --nparen ;
    else if(tok == "{")
      ++nbrace ;
    else if(tok == "}")
      --nbrace ;
    else if(tok == "[")
      ++nbracket ;
    else if(tok == "]")
      --nbracket ;
    else if(tok == ";" &&
            nparen == 0 && nbrace == 0 && nbracket == 0)
      return ;
  }
}

// -----------------------------------------------------------------------
/// @brief Parses the ASCII `FoamFile` header and validates the file
/// format.
///
/// @param [tokens] Tokenized file contents.
/// @param [idx] Current token index, advanced past the header.
/// @param [filename] Source filename used in diagnostics.
/// @param [allow_binary] True when a binary format header is acceptable.
/// @return True when the header is valid or absent and the format is
/// supported.
bool parseFoamHeader(const vector<string> &tokens,
                     size_t &idx,
                     const string &filename,
                     bool allow_binary = false) {
  if(idx >= tokens.size() || tokens[idx] != "FoamFile")
    return true ;

  ++idx ;
  if(!expectToken(tokens,idx,"{",filename))
    return false ;

  string format = "ascii" ;
  while(idx < tokens.size() && tokens[idx] != "}") {
    string key = tokens[idx++] ;
    if(key == "format") {
      if(idx >= tokens.size()) {
        cerr << "missing format value in '" << filename << "'" << endl ;
        return false ;
      }
      format = tokens[idx++] ;
      if(!expectToken(tokens,idx,";",filename))
        return false ;
    } else {
      skipEntryValue(tokens,idx) ;
    }
  }

  if(!expectToken(tokens,idx,"}",filename))
    return false ;

  if(format != "ascii" && !(allow_binary && format == "binary")) {
    cerr << "unsupported OpenFOAM file format '" << format
         << "' in '" << filename << "'" << endl ;
    return false ;
  }
  return true ;
}

// -----------------------------------------------------------------------
/// @brief Extracts the next token from a raw OpenFOAM byte buffer.
///
/// @param [data] Raw file contents.
/// @param [idx] Current byte index, advanced to the end of the token.
/// @param [tok] Parsed token text on success.
/// @return True when a token is found before end of file.
bool nextFoamToken(const vector<char> &data, size_t &idx, string &tok) {
  while(idx < data.size()) {
    unsigned char c = unsigned(data[idx]) ;
    if(isspace(c)) {
      ++idx ;
      continue ;
    }
    if(data[idx] == '/' && idx+1 < data.size()) {
      if(data[idx+1] == '/') {
        idx += 2 ;
        while(idx < data.size() && data[idx] != '\n')
          ++idx ;
        continue ;
      }
      if(data[idx+1] == '*') {
        idx += 2 ;
        while(idx+1 < data.size()) {
          if(data[idx] == '*' && data[idx+1] == '/') {
            idx += 2 ;
            break ;
          }
          ++idx ;
        }
        continue ;
      }
    }
    break ;
  }

  if(idx >= data.size()) {
    tok = "" ;
    return false ;
  }

  if(isDelimiter(data[idx])) {
    tok = string(1,data[idx]) ;
    ++idx ;
    return true ;
  }

  tok = "" ;
  while(idx < data.size()) {
    unsigned char c = unsigned(data[idx]) ;
    if(isspace(c) || isDelimiter(data[idx]))
      break ;
    if(data[idx] == '/' && idx+1 < data.size() &&
       (data[idx+1] == '/' || data[idx+1] == '*'))
      break ;
    tok += data[idx++] ;
  }

  return tok != "" ;
}

// -----------------------------------------------------------------------
/// @brief Skips one entry value in a raw OpenFOAM byte buffer.
///
/// @param [data] Raw file contents being parsed.
/// @param [idx] Current byte index, advanced past the entry terminator.
void skipEntryValue(const vector<char> &data, size_t &idx) {
  int nparen = 0 ;
  int nbrace = 0 ;
  int nbracket = 0 ;
  string tok ;

  while(nextFoamToken(data,idx,tok)) {
    if(tok == "(")
      ++nparen ;
    else if(tok == ")")
      --nparen ;
    else if(tok == "{")
      ++nbrace ;
    else if(tok == "}")
      --nbrace ;
    else if(tok == "[")
      ++nbracket ;
    else if(tok == "]")
      --nbracket ;
    else if(tok == ";" &&
            nparen == 0 && nbrace == 0 && nbracket == 0)
      return ;
  }
}

// -----------------------------------------------------------------------
/// @brief Parses the binary-friendly `FoamFile` header from raw bytes.
///
/// @param [data] Raw file contents.
/// @param [idx] Byte index reset and advanced past the header.
/// @param [filename] Source filename used in diagnostics.
/// @param [format] Declared OpenFOAM file format.
/// @param [class_name] Declared OpenFOAM class name when present.
/// @return True when the header is parsed successfully.
bool parseFoamHeader(const vector<char> &data,
                     size_t &idx,
                     const string &filename,
                     string &format,
                     string &class_name) {
  idx = 0 ;
  string tok ;
  if(!nextFoamToken(data,idx,tok) || tok != "FoamFile") {
    cerr << "unable to find FoamFile header in '" << filename << "'"
         << endl ;
    return false ;
  }
  if(!nextFoamToken(data,idx,tok) || tok != "{") {
    cerr << "unable to parse FoamFile header in '" << filename << "'"
         << endl ;
    return false ;
  }

  format = "ascii" ;
  class_name = "" ;
  while(nextFoamToken(data,idx,tok) && tok != "}") {
    string key = tok ;
    if(key == "format") {
      if(!nextFoamToken(data,idx,format)) {
        cerr << "missing format value in '" << filename << "'" << endl ;
        return false ;
      }
      if(!nextFoamToken(data,idx,tok) || tok != ";") {
        cerr << "missing ';' after format in '" << filename << "'"
             << endl ;
        return false ;
      }
    } else if(key == "class") {
      if(!nextFoamToken(data,idx,class_name)) {
        cerr << "missing class value in '" << filename << "'" << endl ;
        return false ;
      }
      if(!nextFoamToken(data,idx,tok) || tok != ";") {
        cerr << "missing ';' after class in '" << filename << "'"
             << endl ;
        return false ;
      }
    } else {
      skipEntryValue(data,idx) ;
    }
  }

  if(tok != "}") {
    cerr << "unterminated FoamFile header in '" << filename << "'"
         << endl ;
    return false ;
  }
  return true ;
}

// -----------------------------------------------------------------------
/// @brief Reads the leading count and opening parenthesis of a binary
/// OpenFOAM list.
///
/// @param [data] Raw file contents.
/// @param [idx] Current byte index, advanced past the list prefix.
/// @param [filename] Source filename used in diagnostics.
/// @param [count] Parsed number of list entries.
/// @return True when the list prefix is valid.
bool parseBinaryListPrefix(const vector<char> &data,
                           size_t &idx,
                           const string &filename,
                           size_t &count) {
  string tok ;
  if(!nextFoamToken(data,idx,tok)) {
    cerr << "missing binary list size in '" << filename << "'" << endl ;
    return false ;
  }
  int tmp = 0 ;
  if(!tokenToInt(tok,tmp) || tmp < 0) {
    cerr << "invalid binary list size '" << tok << "' in '" << filename
         << "'" << endl ;
    return false ;
  }
  count = tmp ;
  if(!nextFoamToken(data,idx,tok) || tok != "(") {
    cerr << "expected '(' after list size in '" << filename << "'"
         << endl ;
    return false ;
  }
  return true ;
}

// -----------------------------------------------------------------------
/// @brief Decodes a binary OpenFOAM label stored in 4-byte or 8-byte
/// form.
///
/// @param [ptr] Pointer to the encoded label bytes.
/// @param [width] Number of bytes in the encoded label.
/// @param [val] Decoded integer value on success.
/// @return True when the label width is supported and fits in `int`.
bool parseBinaryLabel(const char *ptr, int width, int &val) {
  if(width == 4) {
    int32_t tmp = 0 ;
    memcpy(&tmp,ptr,4) ;
    val = tmp ;
    return true ;
  }
  if(width == 8) {
    int64_t tmp = 0 ;
    memcpy(&tmp,ptr,8) ;
    if(tmp < std::numeric_limits<int>::min() ||
       tmp > std::numeric_limits<int>::max())
      return false ;
    val = int(tmp) ;
    return true ;
  }
  return false ;
}

// -----------------------------------------------------------------------
/// @brief Decodes a binary OpenFOAM scalar stored in 4-byte or 8-byte
/// form.
///
/// @param [ptr] Pointer to the encoded scalar bytes.
/// @param [width] Number of bytes in the encoded scalar.
/// @return Scalar value converted to `double`.
double parseBinaryScalar(const char *ptr, int width) {
  if(width == 4) {
    float tmp = 0.f ;
    memcpy(&tmp,ptr,4) ;
    return tmp ;
  }
  double tmp = 0.0 ;
  memcpy(&tmp,ptr,8) ;
  return tmp ;
}

// -----------------------------------------------------------------------
/// @brief Locates the closing parenthesis that terminates a binary list
/// payload.
///
/// @param [data] Raw file contents.
/// @param [filename] Source filename used in diagnostics.
/// @return Byte offset of the closing ')' or `data.size()` on failure.
size_t findBinaryFooter(const vector<char> &data, const string &filename) {
  for(size_t i=data.size();i>4;--i) {
    size_t p = i-1 ;
    if(data[p] == ')' && p+4 < data.size() &&
       data[p+1] == '\n' && data[p+2] == '\n' &&
       data[p+3] == '/' && data[p+4] == '/')
      return p ;
  }
  cerr << "unable to locate binary list footer in '" << filename << "'"
       << endl ;
  return data.size() ;
}

// -----------------------------------------------------------------------
/// @brief Reads a binary OpenFOAM points file into Cartesian coordinates.
///
/// @param [filename] Binary `points` file.
/// @param [points] Parsed point coordinates on success.
/// @return True when the binary point list is supported and parsed.
bool readBinaryPointsFile(const string &filename,
                          vector<vector3d<double> > &points) {
  vector<char> data ;
  if(!readFileBytes(filename,data))
    return false ;

  size_t idx = 0 ;
  string format, class_name ;
  if(!parseFoamHeader(data,idx,filename,format,class_name))
    return false ;
  if(format != "binary") {
    cerr << "binary reader called for non-binary file '" << filename
         << "'" << endl ;
    return false ;
  }

  size_t count = 0 ;
  if(!parseBinaryListPrefix(data,idx,filename,count))
    return false ;

  size_t footer = findBinaryFooter(data,filename) ;
  if(footer <= idx) {
    cerr << "corrupt binary point list in '" << filename << "'" << endl ;
    return false ;
  }
  size_t body_len = footer - idx ;
  if(count == 0) {
    points.clear() ;
    return true ;
  }
  if(body_len % (count*3) != 0) {
    cerr << "binary point data size mismatch in '" << filename << "'"
         << endl ;
    return false ;
  }
  int scalar_bytes = body_len/(count*3) ;
  if(scalar_bytes != 4 && scalar_bytes != 8) {
    cerr << "unsupported binary scalar width in '" << filename << "'"
         << endl ;
    return false ;
  }

  points.resize(count) ;
  for(size_t i=0;i<count;++i) {
    const char *ptr = &data[idx + i*3*scalar_bytes] ;
    points[i].x = parseBinaryScalar(ptr,scalar_bytes) ;
    points[i].y = parseBinaryScalar(ptr+scalar_bytes,scalar_bytes) ;
    points[i].z = parseBinaryScalar(ptr+2*scalar_bytes,scalar_bytes) ;
  }
  return true ;
}

// -----------------------------------------------------------------------
/// @brief Reads a binary OpenFOAM label list such as `owner` or
/// `neighbour`.
///
/// @param [filename] Binary label-list file.
/// @param [vals] Parsed labels on success.
/// @return True when the binary label list is supported and parsed.
bool readBinaryLabelFile(const string &filename, vector<int> &vals) {
  vector<char> data ;
  if(!readFileBytes(filename,data))
    return false ;

  size_t idx = 0 ;
  string format, class_name ;
  if(!parseFoamHeader(data,idx,filename,format,class_name))
    return false ;
  if(format != "binary") {
    cerr << "binary reader called for non-binary file '" << filename
         << "'" << endl ;
    return false ;
  }

  size_t count = 0 ;
  if(!parseBinaryListPrefix(data,idx,filename,count))
    return false ;

  size_t footer = findBinaryFooter(data,filename) ;
  if(footer < idx) {
    cerr << "corrupt binary label list in '" << filename << "'"
         << endl ;
    return false ;
  }
  size_t body_len = footer - idx ;
  if(count == 0) {
    vals.clear() ;
    return true ;
  }
  if(body_len % count != 0) {
    cerr << "binary label data size mismatch in '" << filename << "'"
         << endl ;
    return false ;
  }
  int label_bytes = body_len/count ;
  if(label_bytes != 4 && label_bytes != 8) {
    cerr << "unsupported binary label width in '" << filename << "'"
         << endl ;
    return false ;
  }

  vals.resize(count) ;
  for(size_t i=0;i<count;++i) {
    if(!parseBinaryLabel(&data[idx + i*label_bytes],label_bytes,vals[i])) {
      cerr << "binary label value out of range in '" << filename << "'"
           << endl ;
      return false ;
    }
  }
  return true ;
}

// -----------------------------------------------------------------------
/// @brief Reads a binary OpenFOAM `faceCompactList` into per-face node
/// connectivity.
///
/// @param [filename] Binary faces file.
/// @param [faces] Parsed face-to-node connectivity on success.
/// @return True when the binary face list is supported and parsed.
bool readBinaryFacesFile(const string &filename, vector<vector<int> > &faces) {
  vector<char> data ;
  if(!readFileBytes(filename,data))
    return false ;

  size_t idx = 0 ;
  string format, class_name ;
  if(!parseFoamHeader(data,idx,filename,format,class_name))
    return false ;
  if(format != "binary") {
    cerr << "binary reader called for non-binary file '" << filename
         << "'" << endl ;
    return false ;
  }
  if(class_name != "faceCompactList") {
    cerr << "unsupported binary face class '" << class_name
         << "' in '" << filename << "'" << endl ;
    return false ;
  }

  size_t offset_count = 0 ;
  if(!parseBinaryListPrefix(data,idx,filename,offset_count))
    return false ;

  int label_bytes = 0 ;
  size_t entry_count = 0 ;
  size_t entries_start = 0 ;
  size_t offsets_end = 0 ;
  for(int candidate=4;candidate<=8;candidate+=4) {
    size_t end1 = idx + offset_count*candidate ;
    if(end1 >= data.size() || data[end1] != ')')
      continue ;
    size_t tmp_idx = end1 + 1 ;
    if(!parseBinaryListPrefix(data,tmp_idx,filename,entry_count))
      continue ;
    size_t footer = findBinaryFooter(data,filename) ;
    if(footer <= tmp_idx)
      continue ;
    size_t body2_len = footer - tmp_idx ;
    if(entry_count == 0 || body2_len != entry_count*size_t(candidate))
      continue ;
    label_bytes = candidate ;
    entries_start = tmp_idx ;
    offsets_end = end1 ;
    break ;
  }

  if(label_bytes == 0) {
    cerr << "unable to determine binary face label width in '"
         << filename << "'" << endl ;
    return false ;
  }
  if(offset_count == 0) {
    faces.clear() ;
    return true ;
  }

  vector<int> offsets(offset_count) ;
  for(size_t i=0;i<offset_count;++i) {
    if(!parseBinaryLabel(&data[idx + i*label_bytes],label_bytes,
                         offsets[i])) {
      cerr << "binary face offset out of range in '" << filename << "'"
           << endl ;
      return false ;
    }
  }

  vector<int> entries(entry_count) ;
  for(size_t i=0;i<entry_count;++i) {
    if(!parseBinaryLabel(&data[entries_start + i*label_bytes],
                         label_bytes,entries[i])) {
      cerr << "binary face entry out of range in '" << filename << "'"
           << endl ;
      return false ;
    }
  }

  if(offsets_end != idx + offset_count*label_bytes) {
    cerr << "binary face offset parsing failed in '" << filename << "'"
         << endl ;
    return false ;
  }
  if(offsets[offset_count-1] != int(entry_count)) {
    cerr << "binary face offsets do not match packed entry count in '"
         << filename << "'" << endl ;
    return false ;
  }

  faces.resize(offset_count-1) ;
  for(size_t i=0;i+1<offset_count;++i) {
    int begin = offsets[i] ;
    int end = offsets[i+1] ;
    if(begin < 0 || end < begin || size_t(end) > entry_count) {
      cerr << "invalid face offsets in '" << filename << "'" << endl ;
      return false ;
    }
    faces[i].resize(end-begin) ;
    for(int j=begin;j<end;++j)
      faces[i][j-begin] = entries[j] ;
  }

  return true ;
}

// -----------------------------------------------------------------------
/// @brief Parses the size prefix of an ASCII OpenFOAM list.
///
/// @param [tokens] Token stream being parsed.
/// @param [idx] Current token index, advanced past the count.
/// @param [filename] Source filename used in diagnostics.
/// @param [count] Parsed list size.
/// @return True when the list size token is present and non-negative.
bool parseListSize(const vector<string> &tokens,
                   size_t &idx,
                   const string &filename,
                   int &count) {
  if(idx >= tokens.size()) {
    cerr << "missing list size in '" << filename << "'" << endl ;
    return false ;
  }
  if(!tokenToInt(tokens[idx],count)) {
    cerr << "expected list size in '" << filename
         << "', got '" << tokens[idx] << "'" << endl ;
    return false ;
  }
  ++idx ;
  if(count < 0) {
    cerr << "negative list size in '" << filename << "'" << endl ;
    return false ;
  }
  return true ;
}

// -----------------------------------------------------------------------
/// @brief Rewrites names into VOG-safe alphanumeric identifiers.
///
/// @param [input_name] Raw patch or zone name from OpenFOAM.
/// @return Sanitized identifier suitable for VOG tags and boundary names.
string sanitizeName(const string &input_name) {
  string name = input_name ;
  if(name == "")
    name = "unnamed_patch" ;

  string tmp = name ;
  if(!(isalpha(unsigned(name[0]))))
    name[0] = '_' ;
  for(size_t i=1;i<name.size();++i) {
    if(!(isalpha(unsigned(name[i]))) &&
       !(isdigit(unsigned(name[i]))))
      name[i] = '_' ;
  }

  if(tmp != name) {
    cerr << "Renaming tag '" << tmp << "' to '" << name << "'!" << endl ;
  }
  return name ;
}

// -----------------------------------------------------------------------
/// @brief Converts a string to lowercase using byte-wise ASCII rules.
///
/// @param [s] Input string.
/// @return Lowercase copy of @p s.
string lowerCase(string s) {
  for(size_t i=0;i<s.size();++i)
    s[i] = tolower(unsigned(s[i])) ;
  return s ;
}

// -----------------------------------------------------------------------
/// @brief Identifies OpenFOAM patch types that are periodic-like and are
/// imported as ordinary boundaries.
///
/// @param [patch_type] Raw OpenFOAM patch type.
/// @return True when the patch type is `cyclic`, `cyclicAMI`, or
/// `cyclicACMI`.
bool periodicPatchType(const string &patch_type) {
  string t = lowerCase(patch_type) ;
  if(t == "cyclic" || t == "cyclicami" || t == "cyclicacmi")
    return true ;
  return false ;
}

// -----------------------------------------------------------------------
/// @brief Identifies OpenFOAM mapped/interface patch types that are
/// preserved only as boundary names.
///
/// @param [patch_type] Raw OpenFOAM patch type.
/// @return True when the patch is a mapped-style interface.
bool mappedPatchType(const string &patch_type) {
  string t = lowerCase(patch_type) ;
  if(t == "mapped" || t == "mappedwall" || t == "mappedpatch")
    return true ;
  return false ;
}

// -----------------------------------------------------------------------
/// @brief Identifies OpenFOAM patch types that the converter refuses to
/// import.
///
/// @param [patch_type] Raw OpenFOAM patch type.
/// @return True when the patch requires unsupported coupled/parallel
/// behavior.
bool rejectPatchType(const string &patch_type) {
  string t = lowerCase(patch_type) ;
  if(t == "processor" || t == "processorcyclic")
    return true ;
  if(t == "regioncoupled")
    return true ;
  return false ;
}

// -----------------------------------------------------------------------
/// @brief Verifies that sanitized patch names remain unique.
///
/// @param [patches] Boundary patch records to validate.
/// @return True when no two sanitized names collide.
bool uniquePatchNames(const vector<FoamPatchInfo> &patches) {
  for(size_t i=0;i<patches.size();++i) {
    for(size_t j=i+1;j<patches.size();++j) {
      if(patches[i].name == patches[j].name) {
        cerr << "boundary name collision after sanitizing patch names: '"
             << patches[i].name << "'" << endl ;
        return false ;
      }
    }
  }
  return true ;
}

// -----------------------------------------------------------------------
/// @brief Verifies that sanitized cell-zone names remain unique.
///
/// @param [cellZones] Cell-zone records to validate.
/// @return True when no two sanitized names collide.
bool uniqueCellZoneNames(const vector<FoamCellZoneInfo> &cellZones) {
  for(size_t i=0;i<cellZones.size();++i) {
    for(size_t j=i+1;j<cellZones.size();++j) {
      if(cellZones[i].name == cellZones[j].name) {
        cerr << "cell zone name collision after sanitizing zone names: '"
             << cellZones[i].name << "'" << endl ;
        return false ;
      }
    }
  }
  return true ;
}

// -----------------------------------------------------------------------
/// @brief Prints a readable summary of imported boundary patches.
///
/// @param [patches] Boundary patches to summarize.
void printPatchSummary(const vector<FoamPatchInfo> &patches) {
  cout << "boundary patches:" << endl ;
  for(size_t i=0;i<patches.size();++i) {
    const FoamPatchInfo &patch = patches[i] ;
    cout << "  [" << patch.id << "] " << patch.name
         << " type=" << patch.type
         << " nFaces=" << patch.nFaces
         << " startFace=" << patch.startFace << endl ;
  }
}

// -----------------------------------------------------------------------
/// @brief Prints a readable summary of imported cell zones.
///
/// @param [cellZones] Cell zones to summarize.
void printCellZoneSummary(const vector<FoamCellZoneInfo> &cellZones) {
  cout << "cell zones:" << endl ;
  if(cellZones.empty()) {
    cout << "  none" << endl ;
    return ;
  }
  for(size_t i=0;i<cellZones.size();++i) {
    const FoamCellZoneInfo &zone = cellZones[i] ;
    cout << "  " << zone.name
         << " nCells=" << zone.cellLabels.size() << endl ;
  }
}

// -----------------------------------------------------------------------
/// @brief Prints a readable summary of VOG volume tags prepared for
/// output.
///
/// @param [volTags] Volume tags to summarize.
void printVolumeTagSummary(const vector<pair<string,entitySet> > &volTags) {
  cout << "volume tags:" << endl ;
  if(volTags.empty()) {
    cout << "  none" << endl ;
    return ;
  }
  for(size_t i=0;i<volTags.size();++i) {
    cout << "  " << volTags[i].first
         << " intervals=" << volTags[i].second.num_intervals() << endl ;
  }
}

// -----------------------------------------------------------------------
/// @brief Parses an ASCII OpenFOAM points file.
///
/// @param [filename] ASCII `points` file.
/// @param [points] Parsed point coordinates on success.
/// @return True when the file is parsed successfully.
bool parsePointsFile(const string &filename,
                     vector<vector3d<double> > &points) {
  vector<string> tokens ;
  if(!tokenizeFoamFile(filename,tokens))
    return false ;

  size_t idx = 0 ;
  if(!parseFoamHeader(tokens,idx,filename))
    return false ;

  int npnts = 0 ;
  if(!parseListSize(tokens,idx,filename,npnts))
    return false ;
  if(!expectToken(tokens,idx,"(",filename))
    return false ;

  points.resize(npnts) ;
  for(int i=0;i<npnts;++i) {
    if(!expectToken(tokens,idx,"(",filename))
      return false ;
    if(idx+2 >= tokens.size()) {
      cerr << "point " << i << " is incomplete in '" << filename << "'"
           << endl ;
      return false ;
    }
    if(!tokenToDouble(tokens[idx++],points[i].x) ||
       !tokenToDouble(tokens[idx++],points[i].y) ||
       !tokenToDouble(tokens[idx++],points[i].z)) {
      cerr << "unable to read point " << i << " in '" << filename << "'"
           << endl ;
      return false ;
    }
    if(!expectToken(tokens,idx,")",filename))
      return false ;
  }

  if(!expectToken(tokens,idx,")",filename))
    return false ;
  return true ;
}

// -----------------------------------------------------------------------
/// @brief Parses an ASCII OpenFOAM faces file.
///
/// @param [filename] ASCII `faces` file.
/// @param [faces] Parsed face-to-node connectivity on success.
/// @return True when the file is parsed successfully.
bool parseFacesFile(const string &filename, vector<vector<int> > &faces) {
  vector<string> tokens ;
  if(!tokenizeFoamFile(filename,tokens))
    return false ;

  size_t idx = 0 ;
  if(!parseFoamHeader(tokens,idx,filename))
    return false ;

  int nfaces = 0 ;
  if(!parseListSize(tokens,idx,filename,nfaces))
    return false ;
  if(!expectToken(tokens,idx,"(",filename))
    return false ;

  faces.resize(nfaces) ;
  for(int i=0;i<nfaces;++i) {
    int fsz = 0 ;
    if(idx >= tokens.size() || !tokenToInt(tokens[idx],fsz)) {
      cerr << "unable to read face size for face " << i
           << " in '" << filename << "'" << endl ;
      return false ;
    }
    ++idx ;
    if(fsz < 0) {
      cerr << "negative face size in '" << filename << "'" << endl ;
      return false ;
    }
    if(!expectToken(tokens,idx,"(",filename))
      return false ;
    faces[i].resize(fsz) ;
    for(int j=0;j<fsz;++j) {
      if(idx >= tokens.size() || !tokenToInt(tokens[idx],faces[i][j])) {
        cerr << "unable to read node index for face " << i
             << " in '" << filename << "'" << endl ;
        return false ;
      }
      ++idx ;
    }
    if(!expectToken(tokens,idx,")",filename))
      return false ;
  }

  if(!expectToken(tokens,idx,")",filename))
    return false ;
  return true ;
}

// -----------------------------------------------------------------------
/// @brief Parses an ASCII OpenFOAM label list.
///
/// @param [filename] ASCII label-list file.
/// @param [vals] Parsed labels on success.
/// @return True when the file is parsed successfully.
bool parseLabelFile(const string &filename, vector<int> &vals) {
  vector<string> tokens ;
  if(!tokenizeFoamFile(filename,tokens))
    return false ;

  size_t idx = 0 ;
  if(!parseFoamHeader(tokens,idx,filename))
    return false ;

  int nvals = 0 ;
  if(!parseListSize(tokens,idx,filename,nvals))
    return false ;
  if(!expectToken(tokens,idx,"(",filename))
    return false ;

  vals.resize(nvals) ;
  for(int i=0;i<nvals;++i) {
    if(idx >= tokens.size() || !tokenToInt(tokens[idx],vals[i])) {
      cerr << "unable to read label " << i << " in '" << filename << "'"
           << endl ;
      return false ;
    }
    ++idx ;
  }

  if(!expectToken(tokens,idx,")",filename))
    return false ;
  return true ;
}

// -----------------------------------------------------------------------
/// @brief Parses a named `List<label>` entry embedded inside an ASCII
/// dictionary.
///
/// @param [tokens] Token stream being parsed.
/// @param [idx] Current token index, advanced past the list entry.
/// @param [filename] Source filename used in diagnostics.
/// @param [entry_name] Human-readable entry name for diagnostics.
/// @param [vals] Parsed label values on success.
/// @return True when the entry is parsed successfully.
bool parseLabelListValue(const vector<string> &tokens,
                         size_t &idx,
                         const string &filename,
                         const string &entry_name,
                         vector<int> &vals) {
  int nvals = -1 ;
  while(idx < tokens.size()) {
    if(tokenToInt(tokens[idx],nvals))
      break ;
    if(tokens[idx] == ";" || tokens[idx] == "(" || tokens[idx] == ")" ||
       tokens[idx] == "{" || tokens[idx] == "}") {
      cerr << "unable to determine list size for " << entry_name
           << " in '" << filename << "'" << endl ;
      return false ;
    }
    ++idx ;
  }
  if(idx >= tokens.size()) {
    cerr << "missing list size for " << entry_name
         << " in '" << filename << "'" << endl ;
    return false ;
  }
  ++idx ;
  if(nvals < 0) {
    cerr << "negative list size for " << entry_name
         << " in '" << filename << "'" << endl ;
    return false ;
  }
  if(!expectToken(tokens,idx,"(",filename))
    return false ;

  vals.resize(nvals) ;
  for(int i=0;i<nvals;++i) {
    if(idx >= tokens.size() || !tokenToInt(tokens[idx],vals[i])) {
      cerr << "unable to read label " << i << " for " << entry_name
           << " in '" << filename << "'" << endl ;
      return false ;
    }
    ++idx ;
  }

  if(!expectToken(tokens,idx,")",filename))
    return false ;
  if(!expectToken(tokens,idx,";",filename))
    return false ;
  return true ;
}

// -----------------------------------------------------------------------
/// @brief Parses an ASCII OpenFOAM boundary file and applies converter
/// policy for supported patch types.
///
/// @param [filename] `boundary` file to parse.
/// @param [patches] Parsed boundary patch records on success.
/// @return True when the file is parsed and all patch types are accepted.
bool parseBoundaryFile(const string &filename,
                       vector<FoamPatchInfo> &patches) {
  vector<string> tokens ;
  if(!tokenizeFoamFile(filename,tokens))
    return false ;

  size_t idx = 0 ;
  if(!parseFoamHeader(tokens,idx,filename,true))
    return false ;

  int npatches = 0 ;
  if(!parseListSize(tokens,idx,filename,npatches))
    return false ;
  if(!expectToken(tokens,idx,"(",filename))
    return false ;

  patches.resize(npatches) ;
  for(int i=0;i<npatches;++i) {
    if(idx >= tokens.size()) {
      cerr << "missing patch name in '" << filename << "'" << endl ;
      return false ;
    }
    patches[i].name = sanitizeName(tokens[idx++]) ;
    if(!expectToken(tokens,idx,"{",filename))
      return false ;

    while(idx < tokens.size() && tokens[idx] != "}") {
      string key = tokens[idx++] ;
      if(key == "type") {
        if(idx >= tokens.size()) {
          cerr << "missing patch type for " << patches[i].name
               << " in '" << filename << "'" << endl ;
          return false ;
        }
        patches[i].type = tokens[idx++] ;
        if(!expectToken(tokens,idx,";",filename))
          return false ;
      } else if(key == "nFaces") {
        if(idx >= tokens.size() ||
           !tokenToInt(tokens[idx],patches[i].nFaces)) {
          cerr << "unable to read nFaces for patch " << patches[i].name
               << " in '" << filename << "'" << endl ;
          return false ;
        }
        ++idx ;
        if(!expectToken(tokens,idx,";",filename))
          return false ;
      } else if(key == "startFace") {
        if(idx >= tokens.size() ||
           !tokenToInt(tokens[idx],patches[i].startFace)) {
          cerr << "unable to read startFace for patch "
               << patches[i].name << " in '" << filename << "'" << endl ;
          return false ;
        }
        ++idx ;
        if(!expectToken(tokens,idx,";",filename))
          return false ;
      } else {
        skipEntryValue(tokens,idx) ;
      }
    }

    if(!expectToken(tokens,idx,"}",filename))
      return false ;

    if(patches[i].nFaces < 0 || patches[i].startFace < 0) {
      cerr << "patch " << patches[i].name
           << " is missing nFaces/startFace in '" << filename << "'"
           << endl ;
      return false ;
    }
    if(rejectPatchType(patches[i].type)) {
      cerr << "patch " << patches[i].name << " has unsupported type '"
           << patches[i].type << "'" << endl ;
      return false ;
    }
    if(periodicPatchType(patches[i].type)) {
      cerr << "WARNING: patch " << patches[i].name << " has type '"
           << patches[i].type << "' and will be imported as a "
           << "boundary patch only; periodic pairing is expected to be "
           << "set by the solver boundary condition." << endl ;
    } else if(mappedPatchType(patches[i].type)) {
      cerr << "WARNING: patch " << patches[i].name << " has type '"
           << patches[i].type << "' and will be imported as a "
           << "boundary patch only; mapped/interface behavior is expected "
           << "to be set by the solver boundary condition." << endl ;
    }
    patches[i].id = i+1 ;
  }

  if(!expectToken(tokens,idx,")",filename))
    return false ;
  if(!uniquePatchNames(patches))
    return false ;
  return true ;
}

// -----------------------------------------------------------------------
/// @brief Parses an ASCII OpenFOAM `cellZones` file.
///
/// @param [filename] `cellZones` file to parse.
/// @param [cellZones] Parsed cell-zone records on success.
/// @return True when the file is parsed successfully.
bool parseCellZonesFile(const string &filename,
                        vector<FoamCellZoneInfo> &cellZones) {
  vector<string> tokens ;
  if(!tokenizeFoamFile(filename,tokens))
    return false ;

  size_t idx = 0 ;
  if(!parseFoamHeader(tokens,idx,filename))
    return false ;

  int nzones = 0 ;
  if(!parseListSize(tokens,idx,filename,nzones))
    return false ;
  if(!expectToken(tokens,idx,"(",filename))
    return false ;

  cellZones.resize(nzones) ;
  for(int i=0;i<nzones;++i) {
    if(idx >= tokens.size()) {
      cerr << "missing cell zone name in '" << filename << "'" << endl ;
      return false ;
    }
    cellZones[i].name = sanitizeName(tokens[idx++]) ;
    if(!expectToken(tokens,idx,"{",filename))
      return false ;

    bool haveCellLabels = false ;
    while(idx < tokens.size() && tokens[idx] != "}") {
      string key = tokens[idx++] ;
      if(key == "cellLabels") {
        if(haveCellLabels) {
          cerr << "cell zone " << cellZones[i].name
               << " defines cellLabels more than once in '" << filename
               << "'" << endl ;
          return false ;
        }
        if(!parseLabelListValue(tokens,idx,filename,
                                string("cell zone '") + cellZones[i].name +
                                "' cellLabels",
                                cellZones[i].cellLabels))
          return false ;
        haveCellLabels = true ;
      } else {
        skipEntryValue(tokens,idx) ;
      }
    }

    if(!expectToken(tokens,idx,"}",filename))
      return false ;
    if(!haveCellLabels) {
      cerr << "cell zone " << cellZones[i].name
           << " is missing cellLabels in '" << filename << "'" << endl ;
      return false ;
    }
  }

  if(!expectToken(tokens,idx,")",filename))
    return false ;
  if(!uniqueCellZoneNames(cellZones))
    return false ;
  return true ;
}

// -----------------------------------------------------------------------
/// @brief Reads optional OpenFOAM cell-zone metadata when present.
///
/// @param [mesh_dir] polyMesh directory that may contain `cellZones`.
/// @param [cellZones] Parsed cell zones on success; cleared when absent.
/// @return True when the optional file is absent or successfully handled.
bool readOptionalCellZones(const string &mesh_dir,
                           vector<FoamCellZoneInfo> &cellZones) {
  string filename = joinPath(mesh_dir,"cellZones") ;
  cellZones.clear() ;
  if(!fileExists(filename))
    return true ;

  vector<char> data ;
  if(!readFileBytes(filename,data))
    return false ;

  size_t idx = 0 ;
  string format, class_name ;
  if(!parseFoamHeader(data,idx,filename,format,class_name))
    return false ;
  if(class_name != "" && class_name != "cellZoneList") {
    cerr << "unsupported cellZones class '" << class_name
         << "' in '" << filename << "'" << endl ;
    return false ;
  }
  if(format == "binary") {
    cerr << "WARNING: binary cellZones file '" << filename
         << "' is not yet supported; volume tags will be omitted." << endl ;
    return true ;
  }

  return parseCellZonesFile(filename,cellZones) ;
}

// -----------------------------------------------------------------------
/// @brief Derives the number of cells implied by the owner and neighbour
/// arrays.
///
/// @param [mesh] OpenFOAM mesh arrays.
/// @param [ncells] Number of cells inferred from connectivity.
/// @return True when the cell count can be determined and indices are
/// valid.
bool determineCellCount(const FoamMeshData &mesh, int &ncells) {
  ncells = 0 ;
  for(size_t i=0;i<mesh.owner.size();++i) {
    if(mesh.owner[i] < 0) {
      cerr << "owner list contains a negative cell index" << endl ;
      return false ;
    }
    if(mesh.owner[i] >= ncells)
      ncells = mesh.owner[i] + 1 ;
  }
  for(size_t i=0;i<mesh.neighbour.size();++i) {
    if(mesh.neighbour[i] < 0) {
      cerr << "neighbour list contains a negative cell index" << endl ;
      return false ;
    }
    if(mesh.neighbour[i] >= ncells)
      ncells = mesh.neighbour[i] + 1 ;
  }
  if(ncells == 0) {
    cerr << "unable to determine number of cells from owner/neighbour"
         << endl ;
    return false ;
  }
  return true ;
}

// -----------------------------------------------------------------------
/// @brief Reads the OpenFOAM mesh arrays needed by the converter from a
/// resolved polyMesh directory.
///
/// @param [mesh_dir] Resolved polyMesh directory.
/// @param [mesh] Aggregate mesh data structure filled on success.
/// @return True when all required mesh files are parsed successfully.
bool readFoamMesh(const string &mesh_dir, FoamMeshData &mesh) {
  string points_file = joinPath(mesh_dir,"points") ;
  string faces_file = joinPath(mesh_dir,"faces") ;
  string owner_file = joinPath(mesh_dir,"owner") ;
  string neighbour_file = joinPath(mesh_dir,"neighbour") ;

  vector<char> points_data ;
  if(!readFileBytes(points_file,points_data))
    return false ;
  size_t idx = 0 ;
  string format, class_name ;
  if(!parseFoamHeader(points_data,idx,points_file,format,class_name))
    return false ;
  if(format == "binary") {
    if(!readBinaryPointsFile(points_file,mesh.points))
      return false ;
  } else {
    if(!parsePointsFile(points_file,mesh.points))
      return false ;
  }

  vector<char> faces_data ;
  if(!readFileBytes(faces_file,faces_data))
    return false ;
  if(!parseFoamHeader(faces_data,idx,faces_file,format,class_name))
    return false ;
  if(format == "binary") {
    if(!readBinaryFacesFile(faces_file,mesh.faces))
      return false ;
  } else {
    if(!parseFacesFile(faces_file,mesh.faces))
      return false ;
  }

  vector<char> owner_data ;
  if(!readFileBytes(owner_file,owner_data))
    return false ;
  if(!parseFoamHeader(owner_data,idx,owner_file,format,class_name))
    return false ;
  if(format == "binary") {
    if(!readBinaryLabelFile(owner_file,mesh.owner))
      return false ;
  } else {
    if(!parseLabelFile(owner_file,mesh.owner))
      return false ;
  }

  vector<char> neighbour_data ;
  if(!readFileBytes(neighbour_file,neighbour_data))
    return false ;
  if(!parseFoamHeader(neighbour_data,idx,neighbour_file,
                      format,class_name))
    return false ;
  if(format == "binary") {
    if(!readBinaryLabelFile(neighbour_file,mesh.neighbour))
      return false ;
  } else {
    if(!parseLabelFile(neighbour_file,mesh.neighbour))
      return false ;
  }

  if(!parseBoundaryFile(joinPath(mesh_dir,"boundary"),mesh.patches))
    return false ;
  if(!readOptionalCellZones(mesh_dir,mesh.cellZones))
    return false ;
  return true ;
}

// -----------------------------------------------------------------------
/// @brief Converts OpenFOAM connectivity into the VOG node, face, and
/// surface structures.
///
/// @param [mesh] Parsed OpenFOAM mesh data.
/// @param [posScale] Unit conversion factor applied to point coordinates.
/// @param [pos] Output VOG node coordinates.
/// @param [cl] Output left-cell map for each face.
/// @param [cr] Output right-cell or boundary-id map for each face.
/// @param [face2node] Output face-to-node connectivity.
/// @param [surf_ids] Output VOG surface-id to name mapping.
/// @return True when the mesh connectivity is valid and converted.
bool buildVOGMesh(const FoamMeshData &mesh,
                  double posScale,
                  store<vector3d<double> > &pos,
                  Map &cl, Map &cr,
                  multiMap &face2node,
                  vector<pair<int,string> > &surf_ids) {
  const int nfaces = mesh.faces.size() ;
  const int npnts = mesh.points.size() ;
  const int ninternal = mesh.neighbour.size() ;

  if(int(mesh.owner.size()) != nfaces) {
    cerr << "owner list size does not match number of faces" << endl ;
    return false ;
  }
  if(ninternal > nfaces) {
    cerr << "neighbour list is longer than face list" << endl ;
    return false ;
  }

  int ncells = 0 ;
  if(!determineCellCount(mesh,ncells))
    return false ;

  vector<int> facePatch(nfaces,0) ;
  surf_ids.clear() ;
  for(size_t i=0;i<mesh.patches.size();++i) {
    const FoamPatchInfo &patch = mesh.patches[i] ;
    if(patch.startFace < ninternal && patch.nFaces > 0) {
      cerr << "patch " << patch.name
           << " overlaps internal faces; coupled patches are not "
           << "supported in this converter" << endl ;
      return false ;
    }
    if(patch.startFace + patch.nFaces > nfaces) {
      cerr << "patch " << patch.name
           << " extends past the end of the face list" << endl ;
      return false ;
    }
    for(int f=patch.startFace;f<patch.startFace+patch.nFaces;++f) {
      if(facePatch[f] != 0) {
        cerr << "face " << f << " appears in more than one patch" << endl ;
        return false ;
      }
      facePatch[f] = patch.id ;
    }
    surf_ids.push_back(pair<int,string>(patch.id,patch.name)) ;
  }

  for(int i=0;i<ninternal;++i) {
    if(facePatch[i] != 0) {
      cerr << "internal face " << i << " appears in boundary data" << endl ;
      return false ;
    }
  }
  for(int i=ninternal;i<nfaces;++i) {
    if(facePatch[i] == 0) {
      cerr << "boundary face " << i << " is not covered by any patch"
           << endl ;
      return false ;
    }
  }

  entitySet pdom = interval(0,npnts-1) ;
  pos.allocate(pdom) ;
  for(int i=0;i<npnts;++i)
    pos[i] = mesh.points[i]*posScale ;

  entitySet fdom = interval(0,nfaces-1) ;
  store<int> count ;
  count.allocate(fdom) ;
  cl.allocate(fdom) ;
  cr.allocate(fdom) ;
  for(int i=0;i<nfaces;++i) {
    count[i] = mesh.faces[i].size() ;
    if(count[i] < 3) {
      cerr << "face " << i << " has fewer than 3 nodes" << endl ;
      return false ;
    }
  }

  face2node.allocate(count) ;
  for(int i=0;i<nfaces;++i) {
    for(int j=0;j<count[i];++j) {
      int node = mesh.faces[i][j] ;
      if(node < 0 || node >= npnts) {
        cerr << "face " << i << " references invalid node " << node
             << endl ;
        return false ;
      }
      face2node[i][j] = node ;
    }
    cl[i] = mesh.owner[i] ;
    if(i < ninternal)
      cr[i] = mesh.neighbour[i] ;
    else
      cr[i] = -facePatch[i] ;
  }

  return true ;
}

// -----------------------------------------------------------------------
/// @brief Converts OpenFOAM cell zones into VOG volume tags.
///
/// @param [mesh] Parsed OpenFOAM mesh data.
/// @param [volTags] Output VOG volume tags built from non-empty cell zones.
/// @return True when the cell-zone data is valid and converted.
bool buildVolumeTags(const FoamMeshData &mesh,
                     vector<pair<string,entitySet> > &volTags) {
  volTags.clear() ;
  if(mesh.cellZones.empty())
    return true ;

  int ncells = 0 ;
  if(!determineCellCount(mesh,ncells))
    return false ;

  for(size_t i=0;i<mesh.cellZones.size();++i) {
    const FoamCellZoneInfo &zone = mesh.cellZones[i] ;
    if(zone.cellLabels.empty())
      continue ;

    vector<int> sortedCells = zone.cellLabels ;
    std::sort(sortedCells.begin(),sortedCells.end()) ;

    entitySet tagSet ;
    int start = -1 ;
    int prev = -1 ;
    for(size_t j=0;j<sortedCells.size();++j) {
      const int cell = sortedCells[j] ;
      if(cell < 0 || cell >= ncells) {
        cerr << "cell zone " << zone.name
             << " references invalid cell " << cell << endl ;
        return false ;
      }
      if(j == 0) {
        start = cell ;
        prev = cell ;
        continue ;
      }
      if(cell == prev)
        continue ;
      if(cell == prev + 1) {
        prev = cell ;
        continue ;
      }
      tagSet += interval(start,prev) ;
      start = cell ;
      prev = cell ;
    }
    tagSet += interval(start,prev) ;
    volTags.push_back(pair<string,entitySet>(zone.name,tagSet)) ;
  }

  return true ;
}

// -----------------------------------------------------------------------
/// @brief Writes the command-line help text for the converter.
///
/// @param [s] Output stream that receives the help text.
/// @return The same stream to support chained output.
std::ostream &showUsage(std::ostream &s) {
  s << "foam2vog converts an OpenFOAM polyMesh into VOG format."
    << endl ;
  s << "It accepts either an OpenFOAM case directory containing"
    << endl ;
  s << "\"constant/polyMesh\", the \"constant\" directory, or the "
    << "polyMesh directory itself." << endl ;
  s << endl ;
  s << "Current support:" << endl ;
  s << "  ASCII and binary points, faces, owner, neighbour, and "
    << "boundary files"
    << endl ;
  s << "  ASCII cellZones files imported as volume tags when present"
    << endl ;
  s << "  Serial execution only" << endl ;
  s << "  OpenFOAM cyclic patches imported as boundary patches; periodic"
    << endl ;
  s << "  pairing is expected to be specified by solver boundary"
    << endl ;
  s << "  conditions" << endl ;
  s << "  OpenFOAM mapped patches imported as boundary patches; interface"
    << endl ;
  s << "  behavior is expected to be specified by solver boundary"
    << endl ;
  s << "  conditions" << endl ;
  s << endl ;
  s << "Usage: foam2vog <options> <case directory | polyMesh directory>"
    << endl ;
  s << "flags:" << endl ;
  s << "  -help : print this help text" << endl ;
  s << "  -v    : display version" << endl ;
  s << "  -o    : disable optimization that reorders nodes and faces"
    << endl ;
  s << "  -in   : input grid is in inches" << endl ;
  s << "  -ft   : input grid is in feet" << endl ;
  s << "  -cm   : input grid is in centimeters" << endl ;
  s << "  -m    : input grid is in meters" << endl ;
  s << "  -mm   : input grid is in millimeters" << endl ;
  s << "  -Lref <units> : 1 unit in input grid is <units> long" << endl ;
  s << endl ;
  return s ;
}

// -----------------------------------------------------------------------
/// @brief Prints command-line usage information to standard error.
void Usage() {
  showUsage(cerr) ;
}

// -----------------------------------------------------------------------
/// @brief Entry point for the OpenFOAM-to-VOG converter.
///
/// @param [ac] Command-line argument count.
/// @param [av] Command-line argument vector.
/// @return Zero on success and nonzero on failure.
int main(int ac, char *av[]) {
  using namespace Loci ;

  bool optimize = true ;
  string Lref = "NOSCALE" ;

  for(int i=1;i<ac;++i) {
    if(!strcmp(av[i],"-help") || !strcmp(av[i],"-h")) {
      showUsage(cout) ;
      return 0 ;
    }
  }

  Loci::Init(&ac,&av) ;

  if(Loci::MPI_processes != 1) {
    cerr << "foam2vog converter can only run serial at present" << endl ;
    Loci::Abort() ;
  }

  while(ac >= 2 && av[1][0] == '-') {
    if(ac >= 3 && !strcmp(av[1],"-Lref")) {
      Lref = av[2] ;
      ac -= 2 ;
      av += 2 ;
    } else if(ac >= 2 && !strcmp(av[1],"-v")) {
      cout << "Loci version: " << Loci::version() << endl ;
      if(ac == 2) {
        Loci::Finalize() ;
        return 0 ;
      }
      --ac ;
      ++av ;
    } else if(ac >= 2 && !strcmp(av[1],"-o")) {
      optimize = false ;
      --ac ;
      ++av ;
    } else if(ac >= 2 && !strcmp(av[1],"-in")) {
      Lref = "1 inch" ;
      --ac ;
      ++av ;
    } else if(ac >= 2 && !strcmp(av[1],"-ft")) {
      Lref = "1 foot" ;
      --ac ;
      ++av ;
    } else if(ac >= 2 && !strcmp(av[1],"-cm")) {
      Lref = "1 centimeter" ;
      --ac ;
      ++av ;
    } else if(ac >= 2 && !strcmp(av[1],"-m")) {
      Lref = "1 meter" ;
      --ac ;
      ++av ;
    } else if(ac >= 2 && !strcmp(av[1],"-mm")) {
      Lref = "1 millimeter" ;
      --ac ;
      ++av ;
    } else {
      cerr << "argument " << av[1] << " is not understood." << endl ;
      Usage() ;
      Loci::Finalize() ;
      return -1 ;
    }
  }

  if(Lref == "NOSCALE") {
    cerr << "Must set grid units!" << endl
         << "Use options -in, -ft, -cm, -m, or -Lref to set grid units."
         << endl ;
    Loci::Finalize() ;
    return -1 ;
  }

  if(ac != 2) {
    Usage() ;
    Loci::Finalize() ;
    return -1 ;
  }

  if(Lref == "")
    Lref = "1 meter" ;
  if(!isdigit(Lref[0]))
    Lref = string("1") + Lref ;

  Loci::UNIT_type tp ;
  istringstream iss(Lref) ;
  iss >> tp ;
  double posScale = tp.get_value_in("meter") ;

  if(Loci::MPI_rank == 0) {
    cout << "input grid file units = " << tp ;
    if(posScale != 1.0)
      cout << " = " << posScale << " meters " ;
    cout << endl ;
  }

  string mesh_dir ;
  string case_name ;
  if(!resolveMeshDirectory(av[1],mesh_dir,case_name)) {
    cerr << "unable to locate OpenFOAM mesh files from '" << av[1] << "'"
         << endl
         << "expected either a case directory containing constant/polyMesh"
         << endl
         << "or the polyMesh directory itself" << endl ;
    Loci::Finalize() ;
    return -1 ;
  }

  if(case_name == "" || case_name == "." || case_name == "..")
    case_name = "foam" ;

  string outfile = case_name + ".vog" ;
  FoamMeshData mesh ;
  if(!readFoamMesh(mesh_dir,mesh)) {
    Loci::Finalize() ;
    return -1 ;
  }

  cout << "resolved polyMesh directory = " << mesh_dir << endl
       << "points = " << mesh.points.size() << endl
       << "faces = " << mesh.faces.size() << endl
       << "owner entries = " << mesh.owner.size() << endl
       << "neighbour entries = " << mesh.neighbour.size() << endl
       << "patches = " << mesh.patches.size() << endl
       << "cell zones = " << mesh.cellZones.size() << endl
       << "output vog file = " << outfile << endl ;
  printPatchSummary(mesh.patches) ;
  printCellZoneSummary(mesh.cellZones) ;

  store<vector3d<double> > pos ;
  Map cl, cr ;
  multiMap face2node ;
  vector<pair<int,string> > surf_ids ;
  if(!buildVOGMesh(mesh,posScale,pos,cl,cr,face2node,surf_ids)) {
    Loci::Finalize() ;
    return -1 ;
  }

  vector<pair<string,entitySet> > volTags ;
  if(!buildVolumeTags(mesh,volTags)) {
    Loci::Finalize() ;
    return -1 ;
  }
  printVolumeTagSummary(volTags) ;

  if(MPI_rank == 0)
    cerr << "orienting faces" << endl ;
  VOG::orientFaces(pos,cl,cr,face2node) ;

  if(MPI_rank == 0)
    cerr << "coloring matrix" << endl ;
  VOG::colorMatrix(pos,cl,cr,face2node) ;

  if(optimize) {
    if(MPI_rank == 0)
      cerr << "optimizing mesh layout" << endl ;
    VOG::optimizeMesh(pos,cl,cr,face2node) ;
  }

  if(MPI_rank == 0)
    cerr << "writing VOG file" << endl ;
  if(volTags.empty())
    Loci::writeVOG(outfile,pos,cl,cr,face2node,surf_ids) ;
  else
    Loci::writeVOG(outfile,pos,cl,cr,face2node,surf_ids,volTags) ;

  Loci::Finalize() ;
  return 0 ;
}
