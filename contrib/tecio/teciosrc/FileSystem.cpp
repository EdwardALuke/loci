#include "FileSystem.h"
#include "CodeContract.h"
 #if defined MSWIN
#   include "UnicodeStringUtils.h"
 #endif
namespace tecplot { namespace filesystem { FILE* fileOpen(std::string const& ___1394, std::string const& ___2504) { REQUIRE(!___1394.empty()); REQUIRE(!___2504.empty());
 #if defined MSWIN
return _wfopen(tecplot::utf8ToWideString(___1394).c_str(), tecplot::utf8ToWideString(___2504).c_str());
 #else
return fopen(___1394.c_str(), ___2504.c_str());
 #endif
} FILE* fileReopen(std::string const& ___1394, std::string const& ___2504, FILE* file) { REQUIRE(!___1394.empty()); REQUIRE(!___2504.empty()); REQUIRE(VALID_REF(file));
 #if defined MSWIN
return _wfreopen(tecplot::utf8ToWideString(___1394).c_str(), tecplot::utf8ToWideString(___2504).c_str(), file);
 #else
return freopen(___1394.c_str(), ___2504.c_str(), file);
 #endif
} }}
