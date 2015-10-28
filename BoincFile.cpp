/**
    Copyright (c) 2013, All Right Reserved
    
    This software is in the public domain, furnished "as is", without technical
    support, and with no warranty, express or implied, as to its usefulness for
    any purpose.
    
    BoincFile.cpp
    Library for managing I/O in files within the BOINC infrastructure.
    
    University of Trento,
    Department of Information Engineering and Computer Science
    
    Trento, fall 2013 - spring 2014

    Authors: (alphabetically ordered) Francesco Asnicar, Luca Masera,
             Paolo Morettin, Nadir Sella, Thomas Tolio.
*/

#include <string>
#include <iostream>
#ifndef _BOINCAPI
#include "boinc_api.h"
#endif
#ifndef _FILESYS
#include "filesys.h"
#endif
#ifndef _BOINCFILE
#include "BoincFile.hpp"
#endif

using namespace std;

/** Opens a file within BOINC.
	
	@param string path
	The path of the file to (resolve and) open.

	@param string mode
	The modality in which the file will be opened.

	@return TRUE if the file is correctly opened, FALSE otherwise.
*/
bool BoincFile::open(string path, string mode) {
	if (boinc_is_standalone()) {
		// no need to resolve the logical path
		wrappedFile = fopen(path.c_str(), mode.c_str());
	} else {
		string resolvedName;
		bool fail = boinc_resolve_filename_s(path.c_str(), resolvedName);
		
		if (fail) {
			cerr << "[E] Cannot resolve \"" << path.c_str() << "\"" << endl;
			return false;
		}

		wrappedFile = boinc_fopen(resolvedName.c_str(), mode.c_str());
	}

	return (wrappedFile != NULL);
}

/** Closes a previously opened file.

	@return TRUE il the file is correctly closed, FALSE otherwise.
*/
bool BoincFile::close() {
	return !fclose(wrappedFile);
}

/** Gets the next line of the opened file.

	@param string& out
	The read line.

	@return TRUE if it's returned a valid (non empty) string, FALSE otherwise.
*/
bool BoincFile::getLine(string& out) {
	string str = "";
	char c = fgetc(wrappedFile);
	
	while ((c != '\n') && (c != EOF)) {
		str.append(1, c);
		c = fgetc(wrappedFile);
	}

	out = str;

	return (strcmp(str.c_str(), "") != 0) ? true : false;
}

/** Writes the given string into a file within BOINC.

	@param string str
	Teh string to write.

	@return TRUE if the writing is succefully done, FALSE otherwise.
*/
bool BoincFile::write(string str) {
	return fputs(str.c_str(), wrappedFile);
}
