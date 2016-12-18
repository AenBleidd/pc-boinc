#ifndef _FILESYS_H_
#define _FILESYS_H_

#include <string.h>
#include <string>

/*#ifdef __cplusplus
extern "C" {
#endif*/

extern int boinc_resolve_filename_s(const char*, std::string&);
extern FILE* boinc_fopen(const char* path, const char* mode);


inline int boinc_resolve_filename_s(const char* path, std::string& resolved)
{
	resolved = path;
	return false;
}

inline FILE* boinc_fopen(const char* path, const char* mode)
{
	return fopen(path, mode);
}


/*#ifdef __cplusplus
}
#endif*/

#endif // _FILESYS_H_
