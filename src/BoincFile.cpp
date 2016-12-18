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
#include "boinc_api.h"
#include "filesys.h"
#include "BoincFile.hpp"

using namespace std;

/** Opens a file within BOINC.
	
	@param string path
	The path of the file to (resolve and) open.

	@param string mode
	The modality in which the file will be opened.

	@return TRUE if the file is correctly opened, FALSE otherwise.
*/
bool BoincFile::open(const string& path, const string& mode) {
	bufPos = 0;
	bufSize = 0;
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
	out.clear();
	
	while(1) {
		if (bufPos == bufSize) {
			bufPos = 0;
			bufSize = fread(buffer, 1, BufferSize - 1, wrappedFile);
			buffer[bufSize] = '\0';
			if (0 == bufSize) // End of file or error
				break;
		}
		
		char* ptr = strchr(buffer + bufPos, '\n');
		if (NULL == ptr) {
			out.append(buffer + bufPos, bufSize - bufPos);
			bufPos = bufSize;
		} else {
			size_t len = ptr - buffer - bufPos;
			out.append(buffer + bufPos, len);
			bufPos += len + 1;
			break;
		}
	}
	
	return !out.empty();
}

/** Writes the given string into a file within BOINC.

	@param string str
	Teh string to write.

	@return TRUE if the writing is succefully done, FALSE otherwise.
*/
bool BoincFile::write(const string& str) {
	return fputs(str.c_str(), wrappedFile);
}
