/**
	Copyright (c) 2013, All Right Reserved
	
	This software is in the public domain, furnished "as is", without technical
	support, and with no warranty, express or implied, as to its usefulness for
	any purpose.
	
	utlity.cpp
	File with utility functions for an efficient computation of the PC
	algorithm for gene networks expansion.
	
	University of Trento,
	Department of Information Engineering and Computer Science
	
	Trento, fall 2013 - spring 2014

	Authors: (alphabetically ordered) Francesco Asnicar, Luca Masera,
			 Paolo Morettin, Nadir Sella, Thomas Tolio.
*/

#define M_SQRT1_2 0.70710678118654757273731092936941422522068023681640625 // definition for M_SQRT1_2
#include <cmath>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <map>
#ifndef _UTILITY
#include "utility.hpp"
#endif
#ifndef _BOINCFILE
#include "BoincFile.hpp"
#endif
#ifndef _ERF
#include "erf.h"
#endif

using namespace std;

/** Compares two pairs of the type <probe identifier, lookup index>.
	
	@param const intpair &l
	First intpair (formally <p1,l1>).
	
	@param const intpair &r
	Second intpair (formally <p1,l1>).

	@return
	TRUE if p1 is less than p2 (formally p1 < p2). Otherwise FALSE.
*/
bool comparator (const intpair &l, const intpair &r) {
	return l.first < r.first; 
}

/** Reads the file TILE modifying both sizes and data.

	@param const string tilePath
	Path of the file containing the selected tiles.

	@param int* &tilesDim
	Reference of the array containing the lengths of the selected tiles.

	@param intpair** &tiles
	Reference of the matrix of pairs <probe identifier, lookup index>.

	@param int &tileRows
	Reference of the number of selected tiles.

	@return nothing.
*/
void readTile(const string tilePath, int* &tilesDim, intpair** &tiles, int &tileRows, string** &experiments, const int experimentsDim) {
	BoincFile tile;
	string line;
	tileRows = 0;

	// count the numbers of rows
	if (!tile.open(tilePath, "r")) {
		cerr << "[E] Failed to open \"" << tilePath << "\" file" << endl;
		exit(1);
	}

	// read the input experiments (first line)
	tile.getLine(line);
	readExperiments(line, experiments, experimentsDim);

	while (tile.getLine(line)) {
		tileRows++;
	}

	tile.close();

	// create tileRows arrays
	tilesDim = new int[tileRows];

	// intialize tilesDim structure
	for (int i = 0; i < tileRows; i++) {
		tilesDim[i] = 0;
	}

	tiles = new intpair*[tileRows];

	if (!tile.open(tilePath, "r")) {
		cerr << "[E] Failed to open \"" << tilePath << "\" file" << endl;
		exit(1);
	}

	// trash the first line, it contains the experiments
	tile.getLine(line);

	// For each line (i.e., tile), counts the number of words in order to know how many nodes there are in the subgraph.
	for (int  i = 0; tile.getLine(line); i++) {
		stringstream stream(line);

		while (getline(stream, line, ' ')) {
			tilesDim[i]++;
		}

		// Instantiate a new array of intpair as long as the previous read.
		tiles[i] = new intpair[tilesDim[i]];
	}

	tile.close();

	if (!tile.open(tilePath, "r")) {
		cerr << "[E] Failed to open \"" << tilePath << "\" file" << endl;
		exit(1);
	}

	// trash the first line, it contains the experiments
	tile.getLine(line);

	// For each tile, extracts the index of involved probes and associates it to a (increasing) lookup index.
	for (int r = 0; r < tileRows; r++) {
		tile.getLine(line);
		stringstream stream(line);

		for (int c = 0; c < tilesDim[r]; c++) {
			getline(stream, line, ' ');
			int temp = atoi(line.c_str());
			tiles[r][c] = make_pair(temp, c);
		}
	}

	tile.close();

	// sort the indexes of the subgraphs
	for (int i = 0; i < tileRows; i++) {
		sort(tiles[i], tiles[i] + tilesDim[i]);
	}
 }
 
 /**
  *
  */
 void readExperiments(string line1, string** &experiments, const int experimentsDim) {
	// BoincFile expp;
	string line;

	// expp.open(experimentsFile, "r");
	// expp.getLine(line);
	// expp.close();

	stringstream stream(line1);
	
	for (int i = 0; (getline(stream, line, ' ') && (i < experimentsDim)); i++) {
		experiments[i] = new string(line);
	}
 }

/** Reads the file CGN saving the biological data that will be used to compute the correlation coefficients.

	@param const string cgnPath
	Path of the file containing the complete gene network.

	@param const intpair* nodesPermutation
	Permutation of the nodes taken into consideration. 

	@param Graph* &g
	The reference of the Graph object representing the gene network.

	@return nothing.
*/
void readCGN(const intpair* nodesPermutation, string** experiments, const int experimentsDim, Graph* &g, const int hibridizationDim) {
	BoincFile cgn;
	string line;
	int column = 0;

	// Initializes the bioData matrix
	g->nCols = hibridizationDim;
	g->bioData = new double*[g->nRows];

	for (int i = 0; i < g->nRows; i++) {
		g->bioData[i] = new double[g->nCols];
	}

	for (int expp = 0; expp < experimentsDim; expp++) {
		// read the file and save the values in bioData
		if (!cgn.open(*(experiments[expp]), "r")) {
			cerr << "[E] Failed to open \"" << *(experiments[expp]) << "\" file" << endl;
			exit(1);
		}

		// counter for the subset graph (aka, tile)
		int c = 0, j = 0;

		// trash the header
		cgn.getLine(line);

		// I start from 0, as the number of rows shows in the text editors
		for (int i = 0; (c < g->nRows) && cgn.getLine(line); i++) {
			// check if it is the rows chosen in readTile()
			if (i == nodesPermutation[c].first) {
				stringstream stream(line);

				for (j = 0; ((j - 1 + column) < hibridizationDim) && getline(stream, line, ','); j++) {
					if (j == 0) { // I'm reading the first token that is the probe id
						if (expp == 0) {
							g->probeIDs[nodesPermutation[c].second] = line;
						}
					} else {
						g->bioData[nodesPermutation[c].second][j - 1 + column] = atof(line.c_str());
						g->means[nodesPermutation[c].second] += g->bioData[nodesPermutation[c].second][j - 1 + column];
					}
				}

				c++;
			}
		}

		if (j > 0) {
			column += j - 1;
		}

		cgn.close();
	}

	// compute means
	for (int i = 0; i < g->nRows; i++) {
		g->means[nodesPermutation[i].second] /= (double) g->nCols;
	}
}

/** Computes the continous density function.
	M_SQRT1_2 takes value 1/sqrt(2).
	
	@param const double value
	Value for which it will be computed its cumulative normal distribution.

	@return The cumulative normal distribution decimal value for the passed parameter.
*/
double comulativeNormalDistribution(const double value) {
	// return 0.5 * erfc(-value * M_SQRT1_2);
	return 0.5 * a_erfc(-value * M_SQRT1_2);
}

/** Finds the correlation coefficient when l is greater than 1 (formally, when l > 1).
	
	@param const int a
	Index of the selected row, that also represents the "departure" node i of an edge i->j.

	@param const int b
	Index of the selected column, that also represents the "arrival" node j of an edge i->j.

	@param const Graph* g
	The Graph object representing the gene network.

	@param const int* neighbours
	Set of the nearby nodes.

	@param const int* subset
	Subset of the lookup indexes of the nearby nodes.
	Note that this subset has cardinality l.

	@param const int l
	Reached dimension actually taken into account for the subset cardinality.
	In this case l is greater than 1 (formally, l > 1).

	@param double ** p
	Rho value.
	
	@return The decimal value of the computed correlation for the edge a->b depending on the given neighbours' subset.
*/
double correlations(const int a, const int b, const Graph* g, const int* neighbours, const int* subset, const int l, double ** p) {
	int dim = l + 2;

	// initialization of p (looks like rho)
	for (int i = 0; i < dim - 1; i++) {
		for (int j = i + 1; j < dim; j++) {
			int first, second;

			if (i == 0) {
				first = a;
			} else if (i == 1) {
				first = b;
			} else {
				first = neighbours[subset[i - 2]];
			}

			if (j == 1) {
				second = b;
			} else {
				second = neighbours[subset[j - 2]];
			}

			p[i][j] = p[j][i] = g->rho[first][second];
		}
	}

	for (int k = 1; k <= l; k++) {
		for (int i = 0; i <= (l - k); i++) {
			for (int j = i + 1; j < (dim - k); j++) {
				p[i][j] = p[j][i] = (p[i][j] - p[i][dim - k] * p[j][dim - k]) / (sqrt((1 - pow(p[i][dim - k], 2)) * (1 - pow(p[j][dim - k], 2))));
			}
		}
	}

	return p[0][1];
}

/** Checks if a given string (of the form array of chars) whether representing a float number or not.
	
	@param const char* number
	String (or, rather, array of characters) that should represent a decimal number.

	@return TRUE if the string follows the correct format for representing a float number. FALSE otherwise.
*/
bool isFloat(const char* number) {
	bool noMorePoint = false;

	for (int i = 0; number[i] != '\0'; i++) {
		if ((number[i] < '0') || (number[i] > '9')) {
			if (number[i] == '.') {
				if (!noMorePoint) {
					noMorePoint = true;
				} else {
					return false;
				}
			}
		}
	}

	return true;
}

/** Prints the uncutted edges of the graph in a .csv file.

	@param Graph* g
	The Graph object representing the gene network.

	@param const intpair* nodesPermutation
	Permutation of the nodes taken into consideration. 

	@param const string fileName
	Name of the output file.

	@return nothing.
*/
void fprintEdgesCSV(Graph* g, const string resultsNotSaved, const intpair* nodesPermutation, const string fileName) {
	BoincFile out;
	ostringstream strs;

	strs.str(std::string()); // cleans the stringstream

	if (out.open(fileName, "ab")) {
		// check if there are results not saved...
		if (strcmp(resultsNotSaved.c_str(), "") != 0) {
			//... if there are then save them!
			out.write(resultsNotSaved);
		}

		// print the edges
		for (int i = 0; i < g->nRows - 1; ++i) {
			for (int j = i + 1; j < g->nRows; ++j) {
				if (g->matrix[nodesPermutation[i].second][nodesPermutation[j].second]) {
					strs << g->probeIDs[nodesPermutation[i].second] << "," << g->probeIDs[nodesPermutation[j].second] << "\n";
				}
			}
		}

		out.write(strs.str());
		out.close();
	} else {
		cerr << "[E] Cannot open \"" << fileName << "\"" << endl;
	}
 }

/** Counts the number of (uncutted) edges in the graph.

	@param bool** matrix
	Matrix of booleans representing a tabular form of the presence and absence of all the edges in the graph.
	The boolean located in the cell of i-th row and j-th column represents the presence/absence of the edge i->j.

	@param const int rows
	The number of rows in the matrix.

	@param const int cols
	The number of rows in the matrix.

	@return The number of uncutted edges.
*/
int countArcs(bool** matrix, const int rows, const int cols) {
	int counter = 0;

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			if (matrix[i][j]) counter++;
		}
	}

	return counter;
}

/**
 *
 */
void calculatePostProcessing(const string filename) {
	map<string, int> postProcessing;
	BoincFile file;
	string line;
	ostringstream strs;

	file.open(filename, "rb");
	
	// count the numbers of rows
	while (file.getLine(line)) {
		if(postProcessing[line] == 0) {
			postProcessing[line] = 1;
		} else {
			postProcessing[line]++;
		}
	}

	file.close();
	
	// write the output file in the format "arc<TAB>count"
	file.open(filename, "wb");
	
	for (map<string, int>::iterator it = postProcessing.begin(); it != postProcessing.end(); it++) {
		strs << it->first << "\t" << it->second << "\n";
	}

	file.write(strs.str());
	file.close();
}
