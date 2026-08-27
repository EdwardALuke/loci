//#############################################################################
//#
//# Copyright 2008-2026, Mississippi State University
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

#include "fvmadapt_core.h"

#include <cerrno>
#include <cstring>
#include <sys/stat.h>

/// Format split-plan char codes as a readable integer list.
std::string plan_string(const std::vector<char>& plan) {
  std::ostringstream out ;
  out << "[" ;
  for(size_t i = 0; i < plan.size(); ++i) {
    if(i != 0) {
      out << "," ;
    }
    out << int(plan[i]) ;
  }
  out << "]" ;
  return out.str() ;
}

/// Create an output directory if it does not already exist.
void make_directory(const std::string& path) {
  if(path.empty()) {
    return ;
  }
  if(mkdir(path.c_str(), 0775) != 0 && errno != EEXIST) {
    std::ostringstream msg ;
    msg << "could not create directory '" << path << "': " << std::strerror(errno) ;
    throw std::runtime_error(msg.str()) ;
  }
}

/// Join a directory and file name with one path separator.
std::string join_path(const std::string& dir, const std::string& file) {
  if(dir.empty()) {
    return file ;
  }
  if(dir[dir.size() - 1] == '/') {
    return dir + file ;
  }
  return dir + "/" + file ;
}
