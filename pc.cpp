/**
    Copyright (c) 2013, All Right Reserved
    
    This software is in the public domain, furnished "as is", without technical
    support, and with no warranty, express or implied, as to its usefulness for
    any purpose.
    
    pc.cpp
    [TODO]
    
    University of Trento,
    Department of Information Engineering and Computer Science
    
    Trento, fall 2013 - spring 2014

    Authors: (alphabetically ordered) Francesco Asnicar, Luca Masera,
             Paolo Morettin, Nadir Sella, Thomas Tolio.
*/

#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#ifndef _GRAPH
#include "Graph.hpp"
#endif
#ifndef _UTILITY
#include "utility.hpp"
#endif
#ifndef _BOINCFUNCTIONS
#include "boinc_functions.hpp"
#endif
#ifndef _PC
#include "pc.hpp"
#endif
#ifndef _BOINCAPI
#include "boinc_api.h"
#endif
#ifndef _BOINCFILE
#include "BoincFile.hpp"
#endif

using namespace std;

/** Tests the causality and, when appropriate, cuts both the edges r->c and c->r.
	
	Tests the causality w.r.t. the correlation coefficient calculating the probability value of causal dependency.
	Whether this value is greater than the alpha value (complement of confidence interval), both the edges r->c
	and c->r are cutted. 
	
	@param double correlationCoefficient
	Correlation coefficient of the edge r->c.

	@param Graph* &g
	The reference of the Graph object representing the gene network.

	@param const int r
	Index of the selected row, that also represents the node r of an edge r->c.

	@param const int c
	Index of the selected column, that also represents the node c of an edge r->c.

	@param const int l
	Reached dimension actually taken into account for the subset cardinality.
	In this case l is greater than 1 (formally, l > 1).	

	@param const double alpha
	Complement of the confidence interval.

	@return nothing.	
*/
void testAndRemove(double correlationCoefficient, Graph* &g, const int r, const int c, const int l, const double alpha) {
	double pVal;
	double const cutAt = 0.9999999;
	bool NAdelete = true;

	// increment score only every 1000 calls of testAndRemove
	if (g->count_score >= 1000) {
		g->score++;
		g->count_score = 1;
	} else {
		g->count_score++;
	}

#ifdef _MSC_VER
	if (_isnan(correlationCoefficient)) {
#else
	if (std::isnan(correlationCoefficient)) {
#endif
		correlationCoefficient = 0.0;
	}
	

	correlationCoefficient = min(cutAt, max(-cutAt, correlationCoefficient));
	pVal = sqrt(g->nCols - l - 3.0) * 0.5 * log((1.0 + correlationCoefficient) / (1.0 - correlationCoefficient));

#ifdef _MSC_VER
	if (_isnan(pVal)) {
#else
	if (std::isnan(pVal)) {
#endif
		pVal = 0.0;
	}

	pVal = 2.0 * (1.0 - comulativeNormalDistribution(abs(pVal)));

#ifdef _MSC_VER
	if (_isnan(pVal)) {
#else
	if (std::isnan(pVal)) {
#endif
		pVal = NAdelete ? 1.0 : 0.0;
	}

	// test d-separation
	if (pVal >= alpha) {
		// remove arcs
		g->matrix[r][c] = g->matrix[c][r] = false;

		// decrement neighbours
		g->numNeighbours[r]--;
		g->numNeighbours[c]--;
	}
}

/** Computes the correlation coefficient related to the edge r->c considering
	subsets of cardinality l of neaighbour edges.
	
	@param int* neighbours
	Set of the indexes of the nearby nodes of the given edge i->j

	@param const int* subset
	Subset of the lookup indexes of the nearby nodes.
	Note that this subset has cardinality l.

	@param const int l
	Reached dimension actually taken into account for the subset cardinality.
	In this case l is greater than 1 (formally, l > 1).

	@param Graph* &g
	The reference of the Graph object representing the gene network.

	@param const int r
	Index of the selected row, that also represents the node r of an edge r->c.

	@param const int c
	Index of the selected column, that also represents the node c of an edge r->c.

	@param double** p
	Rho.

	@return The coefficient of correlation.	
*/
double getCorrelationCoefficient(const int* neighbours, const int* subset, const int l, Graph* &g, const int r, const int c, double** p) {
	double correlationCoefficient;

	if (l == 2) {
		double rhoRC, rhoRH1, rhoCH1;
		int h1 = neighbours[subset[0]];
		int h2 = neighbours[subset[1]];

		rhoRC = (g->rho[r][c] - (g->rho[r][h2] * g->rho[c][h2])) / sqrt((1 - pow(g->rho[r][h2], 2)) * (1 - pow(g->rho[c][h2], 2)));
		rhoRH1 = (g->rho[r][h1] - (g->rho[r][h2] * g->rho[h1][h2])) / sqrt((1 - pow(g->rho[r][h2], 2)) * (1 - pow(g->rho[h1][h2], 2)));
		rhoCH1 = (g->rho[c][h1] - (g->rho[c][h2] * g->rho[h1][h2])) / sqrt((1 - pow(g->rho[c][h2], 2)) * (1 - pow(g->rho[h1][h2], 2)));

		correlationCoefficient = (rhoRC - (rhoRH1 * rhoCH1)) / sqrt((1 - pow(rhoRH1, 2)) * (1 - pow(rhoCH1, 2)));
	} else if (l == 3) {
		double rhoRC_H2H3, rhoRH1_H2H3, rhoCH1_H2H3, rhoRC_H3, rhoRH1_H3, rhoRH2_H3, rhoCH1_H3, rhoCH2_H3, rhoH1H2_H3;
		int h1 = neighbours[subset[0]];
		int h2 = neighbours[subset[1]];
		int h3 = neighbours[subset[2]];

		rhoRC_H3 = (g->rho[r][c] - (g->rho[r][h3] * g->rho[c][h3])) / sqrt((1 - pow(g->rho[r][h3], 2)) * (1 - pow(g->rho[c][h3], 2)));
		rhoRH1_H3 = (g->rho[r][h1] - (g->rho[r][h3] * g->rho[h1][h3])) / sqrt((1 - pow(g->rho[r][h3], 2)) * (1 - pow(g->rho[h1][h3], 2)));
		rhoRH2_H3 = (g->rho[r][h2] - (g->rho[r][h3] * g->rho[h2][h3])) / sqrt((1 - pow(g->rho[r][h3], 2)) * (1 - pow(g->rho[h2][h3], 2)));
		rhoCH1_H3 = (g->rho[c][h1] - (g->rho[c][h3] * g->rho[h1][h3])) / sqrt((1 - pow(g->rho[c][h3], 2)) * (1 - pow(g->rho[h1][h3], 2)));
		rhoCH2_H3 = (g->rho[c][h2] - (g->rho[c][h3] * g->rho[h2][h3])) / sqrt((1 - pow(g->rho[c][h3], 2)) * (1 - pow(g->rho[h2][h3], 2)));
		rhoH1H2_H3 = (g->rho[h1][h2] - (g->rho[h1][h3] * g->rho[h2][h3])) / sqrt((1 - pow(g->rho[h1][h3], 2)) * (1 - pow(g->rho[h2][h3], 2)));

		rhoRC_H2H3 = (rhoRC_H3 - (rhoRH2_H3 * rhoCH2_H3)) / sqrt((1 - pow(rhoRH2_H3, 2)) * (1 - pow(rhoCH2_H3, 2)));
		rhoRH1_H2H3 = (rhoRH1_H3 - (rhoRH2_H3 * rhoH1H2_H3)) / sqrt((1 - pow(rhoRH2_H3, 2)) * (1 - pow(rhoH1H2_H3, 2)));
		rhoCH1_H2H3 = (rhoCH1_H3 - (rhoCH2_H3 * rhoH1H2_H3)) / sqrt((1 - pow(rhoCH2_H3, 2)) * (1 - pow(rhoH1H2_H3, 2)));

		correlationCoefficient = (rhoRC_H2H3 - (rhoRH1_H2H3 * rhoCH1_H2H3)) / sqrt((1 - pow(rhoRH1_H2H3, 2)) * (1 - pow(rhoCH1_H2H3, 2)));
	} else {
		correlationCoefficient = correlations(r, c, g, neighbours, subset, l, p);
	}

	return correlationCoefficient;
}

/** [TODO]

	@param int* neighbours
	Set of the indexes of the nearby nodes of the given edge i->j.

	@param const int neighboursDim
	The length of the array of nwighbours.

	@param const int l
	Reached dimension actually taken into account for the subset cardinality.
	In this case l is greater than 1 (formally, l > 1).

	@param Graph* &g
	The reference of the Graph object representing the gene network.

	@param const int r
	Index of the selected row, that also represents the node r of an edge r->c.

	@param const int c
	Index of the selected column, that also represents the node c of an edge r->c.

	@param const double alpha
	Complement of the confidence interval.

	@param double** p
	Rho.

	@return nothing. 
	
	@see Thanks to Ilya Bursov, http://stackoverflow.com/questions/19327847/n-choose-k-for-large-n-and-k
*/
void iterativeComb(int* neighbours, const int neighboursDim, const int l, Graph* &g, const int r, const int c, const double alpha, double** p) {
	double coeff;
	int* currentCombination = new int[neighboursDim];

	for (int i = 0; i < neighboursDim; i++) {
		currentCombination[i] = i;
	}

	currentCombination[l - 1] = l - 1 - 1;

	do {
		if (currentCombination[l - 1] == (neighboursDim - 1)) {
			int i = l - 1 - 1;

			while (currentCombination[i] == (neighboursDim - l + i)) {
				i--;
			}

			currentCombination[i]++;

			for (int j = i + 1; j < l; j++) {
				currentCombination[j] = currentCombination[i] + j - i;
			}
		} else {
			currentCombination[l - 1]++;
		}

		coeff = getCorrelationCoefficient(neighbours, currentCombination, l, g, r, c, p);
		testAndRemove(coeff, g, r, c, l, alpha);
	} while (g->matrix[r][c] && !((currentCombination[0] == (neighboursDim - 1 - l + 1)) && (currentCombination[l - 1] == (neighboursDim - 1))));

	delete(currentCombination);
}

/** Finds all the subsets adj(i)\{j} with cardinality equals to l (formally, |adj(i)\{j}| = l).
	
	@param Graph* &g
	The reference of the Graph object representing the gene network.

	@param const int i
	Index of the selected row, that also represents the "departure" node i of an edge i->j.

	@param const int j
	Index of the selected column, that also represents the "arrival" node j of an edge i->j.

	@param const int l
	Reached dimension actually taken into account for the subset cardinality.
	In this case l is greater than 1 (formally, l > 1).

	@param bool &hasWorked


	@param const double alpha
	Complement of the confidence interval

	@param int* neighbours
	Set of the indexes of the nearby nodes of the given edge i->j

	@param double** p
	Rhos.

	@return nothing.
*/
void findAllSubsets(Graph* &g, const int i, const int j, const int l, const double alpha, int* neighbours, double** p) {
	int neighboursDim = 0;

	if (l == 0) {
		testAndRemove(g->rho[i][j], g, i, j, l, alpha);
	} else if (l == 1) {
		// find neighbours (when l > 0) of i without considering j
		for (int k = 0; (g->matrix[i][j]) && (k < g->nRows) && (neighboursDim < g->numNeighbours[i]); k++) {
			if (g->matrix[i][k] && (k != j)) {
				neighboursDim++;
				double correlationCoefficient = (g->rho[i][j] - (g->rho[i][k] * g->rho[j][k])) / sqrt((1 - pow(g->rho[i][k], 2)) * (1 - pow(g->rho[j][k], 2)));
				testAndRemove(correlationCoefficient, g, i, j, l, alpha);
			}
		}
	} else {
		// find neighbours (when l > 0) of i without considering j
		for (int k = 0; (k < g->nRows) && (neighboursDim < g->numNeighbours[i]); k++) {
			if ((g->matrix[i][k]) && (k != j)) {
				neighbours[neighboursDim++] = k;
			}
		}

		// look for all subset of length l
		iterativeComb(neighbours, neighboursDim, l, g, i, j, alpha, p);
	}
}

/** Body of the PC algorithm.
	
	@see [TODO reference to the paper]
*/
int main(int argc, char* argv[]) {
	int l;
	int experimentsDim; // number of experiment files to be read
	int hibridizationDim; // total number of columns of all experiment files
	int* tilesDim = NULL; // array that contains the length of each tiles read
	int tileRows; // number of rows of the tile.txt file 
	bool hasWorked; // boolean to see that there is at least an arc i,j s.t. |adj(C, i)\{j}| >= l TO CHECK FROM THE TEXT
	intpair** tiles = NULL; // array that contains the numbers of the rows that will be extracted and represent the subgraph
	Graph* g = NULL; // the graph
	double alpha;
	int* neighbours = NULL;
	double** p = NULL;
	string* tilePath = NULL;
	string* experimentsFile = NULL;
	string* outputFile = NULL;
	string** experiments = NULL;
	string str;
	int retval;
	char buf[256];
	ostringstream strs;
	string resultsNotSaved = ""; // carries the results that were not allowed to write
	time_t rawtime;
	struct tm* timeinfo;
	char* tmp;
	long long unsigned int score = 0;

	if (argc != 6) {
		cerr << "[E] wrong input parameters" << endl;
		cerr << "Usage: ./pc.exe tiles output_file alpha n_experiments n_hibridizations" << endl;
		return -1;
	}

	retval = boinc_init();
	
	if (retval) {
		cerr << boinc_msg_prefix(buf, sizeof(buf)) << " boinc_init() returned: " << retval;
		exit(retval);
	}

	// print a timestamp
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	tmp = asctime(timeinfo);
	cerr << "Start  @ " << tmp;
	
	tilePath = new string(argv[1]);
	outputFile = new string(argv[2]);
	experimentsDim = atoi(argv[4]);
	hibridizationDim = atoi(argv[5]);

	if (isFloat(argv[3])) {
		alpha = atof(argv[3]);
	} else {
		alpha = 0.05;
		cerr << "[E] alpha not a float. Use the default value " << alpha << endl;
	}

	// istantiate the array of string
	experiments = new string*[experimentsDim];

	// read all the subgragraph that will be computed
	readTile(*tilePath, tilesDim, tiles, tileRows, experiments, experimentsDim);

	for (int c = readCheckpoint(*(outputFile), score); c < tileRows; c++) {
		// start with l = 0 & hasWorked = true
		l = -1;
		hasWorked = true;

		// update progress-bar
		progressBar(c, tileRows);

		// create the graph
		g = new Graph(tilesDim[c]);

		// extract the subgraph from the complete genes file
		readCGN(tiles[c], experiments, experimentsDim, g, hibridizationDim);

		// compute the standard deviations
		g->computeStandardDeviations();
		
		// compute the correlations coefficients
		g->computeCorrelations();

		// alloc the rho for correlations() (save time)
		p = new double*[g->nRows];

		for (int i = 0; i < g->nRows; i++) {
			p[i] = new double[g->nRows];
		}

		// alloc the neighbours array (save time)
		neighbours = new int[g->nRows];

		// PC-algorithm
		while ((hasWorked) && (l < g->nRows)) {
			l++;
			hasWorked = false;

			for (int i = 0; i < g->nRows; i++) {
				for (int j = 0; j < g->nRows; j++) {
					// check if exists the arc between i and j
					if (g->matrix[i][j] && (g->numNeighbours[i] > l)) {
						hasWorked = true;
						findAllSubsets(g, i, j, l, alpha, neighbours, p);
					}
				}
			}
		}

		/*** START CRITICAL SECTION ***/
		// tries to checkpoint
		if (writeCheckpoint(c, score)) {
			// print output file
			fprintEdgesCSV(g, resultsNotSaved, tiles[c], *(outputFile));

			// chiude la critical section
			boinc_checkpoint_completed();
		/*** END of CRITICAL SECTION ***/

			resultsNotSaved = "";
		} else {
			// save the results in g and, hopefully, write the next time
			strs.str(std::string()); // cleans the stringstream

			for (int i = 0; i < g->nRows - 1; ++i) {
				for (int j = i + 1; j < g->nRows; ++j) {
					if (g->matrix[tiles[c][i].second][tiles[c][j].second]) {
						strs << g->probeIDs[tiles[c][i].second] << "," << g->probeIDs[tiles[c][j].second] << "\n";
					}
				}
			}

			resultsNotSaved.append(strs.str());
		}

		// save the partial score
		score += g->score;
		
		// free the memory		
		delete[] neighbours;

		for (int i = 0; i < g->nRows; i++) {
			delete[] p[i];
		}
		delete[] p;

		delete g;
	}

	/*** START CRITICAL SECTION ***/
	if (!boinc_is_standalone()) {
		boinc_begin_critical_section();

		// write the last results not saved
		if (strcmp(resultsNotSaved.c_str(), "") != 0) {
			BoincFile out;

			if (out.open(*(outputFile), "ab")) {
				//... if there are then save them!
				out.write(resultsNotSaved);
				out.close();
			} else {
				cerr << "[E] Cannot open \"" << *(outputFile) << "\"" << endl;
			}
		}

		//calculate the post processing to return a more compressed file
		calculatePostProcessing(*(outputFile));

		// print score
		cerr << "|" << score << "|" << endl;

		boinc_end_critical_section();
	/*** END of CRITICAL SECTION ***/
	}

	// print a timestamp
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	tmp = asctime(timeinfo);
	cerr << "Finish @ " << tmp;
	
	// free the memory
	for (int i = 0; i < tileRows; i++) {
		delete[] tiles[i];
	}
	delete[] tiles;
	
	for (int i = 0; i < experimentsDim; i++) {
		delete experiments[i];
	}
	delete[] experiments;

	delete[] tilesDim;
	delete experimentsFile;
	delete tilePath;
	delete outputFile;

	return boinc_finish(0);
}
