 #pragma once
#include "ThirdPartyHeadersBegin.h"
#  include <sstream>
#  include <stdio.h>
#  include <string>
#include "ThirdPartyHeadersEnd.h"
 #ifdef _MSC_VER
 #define snprintf std_snprintf
namespace{ inline int std_vsnprintf(char* str, size_t size, const char* format, va_list ap) { int count = -1; if (size != 0) count = _vsnprintf_s(str, size, _TRUNCATE, format, ap); if (count == -1) count = _vscprintf(format, ap); return count; } inline int std_snprintf(char* str, size_t size, const char* format, ...) { int count; va_list ap; va_start(ap, format); count = std_vsnprintf(str, size, format, ap); va_end(ap); return count; } }
 #endif 
template <typename T> std::string ___4187(T ___4314) { std::ostringstream o; o << ___4314; return o.str(); } namespace tecplot { template <typename T> std::string toString(T ___4314) { std::ostringstream o; o << ___4314; return o.str(); } }
