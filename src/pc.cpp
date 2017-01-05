/**
	Copyright (c) 2013,2016, All Rights Reserved.
	
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

	Optimizations by Daniel Fruzynski
*/

#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include "Graph.hpp"
#include "utility.hpp"
#include "boinc_functions.hpp"
#include "pc.hpp"
#include "boinc_api.h"
#include "BoincFile.hpp"
#include "DoubleVector.hpp"

#ifdef _MSC_VER
#define __attribute__(x)
#endif

using namespace std;

#ifdef __SSE2__
static const __m128d vd2_1_0 = { 1.0, 1.0 };
#endif
#ifdef __AVX__
static const __m256d vd4_1_0 = { 1.0, 1.0, 1.0, 1.0 };
#endif

__attribute__((hot, always_inline)) static inline
double pow2(double val)
{
	return val * val;
}

__attribute__((hot, always_inline)) static inline
double getCorrelationCoefficient1(double a, double b, double c)
{
	return (a - (b * c)) / sqrt((1 - pow2(b)) * (1 - pow2(c)));
}

__attribute__((hot, always_inline)) static inline
void getCorrelationCoefficient2(double a1, double b1, double c1, double a2, double b2, double c2, double* __restrict__ result)
{
#ifdef __SSE2__
	__m128d va = _mm_set_pd(a2, a1);
	__m128d vb = _mm_set_pd(b2, b1);
	__m128d vc = _mm_set_pd(c2, c1);
	__m128d vres = (va - vb * vc) / _mm_sqrt_pd((vd2_1_0 - vb * vb) * (vd2_1_0 - vc * vc));
	_mm_store_pd(result, vres);
#else
	result[0] = getCorrelationCoefficient1(a1, b1, c1);
	result[1] = getCorrelationCoefficient1(a2, b2, c2);
#endif
}

__attribute__((hot, always_inline)) static inline
void getCorrelationCoefficient4(double a1, double b1, double c1, double a2, double b2, double c2,
		double a3, double b3, double c3, double a4, double b4, double c4, double* __restrict__ result)
{
#ifdef __AVX__
	__m256d va = _mm256_set_pd(a4, a3, a2, a1);
	__m256d vb = _mm256_set_pd(b4, b3, b2, b1);
	__m256d vc = _mm256_set_pd(c4, c3, c2, c1);
	__m256d vres = (va - vb * vc) / _mm256_sqrt_pd((vd4_1_0 - vb * vb) * (vd4_1_0 - vc * vc));
	_mm256_store_pd(result, vres);
#else
	getCorrelationCoefficient2(a1, b1, c1, a2, b2, c2, result);
	getCorrelationCoefficient2(a3, b3, c3, a4, b4, c4, result + 2);
#endif
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
__attribute__((hot, always_inline)) static inline
double correlations(const int a, const int b, const Graph* __restrict__ g, const int* __restrict__ neighbours,
	const int* __restrict__ subset, const int l, double * __restrict__ * __restrict__ p) {
	int dim = l + 2;

	// initialization of p (looks like rho)
	for (int i = 0; i < dim - 1; i++) {
		int first;
		if (i == 0) {
			first = a;
		} else if (i == 1) {
			first = b;
		} else {
			first = neighbours[subset[i - 2]];
		}

		for (int j = i + 1; j < dim; j++) {
			int second;

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
			const double p_i_dim_SUB_k = p[i][dim - k];
			const double _1_SUB_p_i_dim_SUB_k_sqr = 1 - pow2(p_i_dim_SUB_k);
#ifdef __AVX__
			const __m256d v4_p_i_dim_SUB_k = _mm256_set1_pd(p_i_dim_SUB_k);
			const __m256d v4_1_SUB_p_i_dim_SUB_k_sqr = _mm256_set1_pd(_1_SUB_p_i_dim_SUB_k_sqr);
#endif
#ifdef __SSE2__
			const __m128d v2_p_i_dim_SUB_k = _mm_set1_pd(p_i_dim_SUB_k);
			const __m128d v2_1_SUB_p_i_dim_SUB_k_sqr = _mm_set1_pd(_1_SUB_p_i_dim_SUB_k_sqr);
#endif
			
			int j = i + 1;
#ifdef __AVX__
			for (; j < (dim - k) - 3; j += 4) {
				__m256d v4_val = _mm256_loadu_pd(&p[dim - k][j]);
				__m256d v4_2 = _mm256_loadu_pd(&p[i][j]);
				__m256d v4;
				v4 = v4_1_SUB_p_i_dim_SUB_k_sqr * (vd4_1_0 - v4_val * v4_val);
				v4 = _mm256_sqrt_pd(v4);
				v4 = (v4_2 - v4_p_i_dim_SUB_k * v4_val) / v4;
				_mm256_storeu_pd(&p[i][j], v4);
				
				p[j+0][i] = p[i][j+0];
				p[j+1][i] = p[i][j+1];
				p[j+2][i] = p[i][j+2];
				p[j+3][i] = p[i][j+3];
			}
#endif

#if defined(__AVX__) || defined(__SSE2__)
#ifdef __AVX__
			if (j < (dim - k) - 1) {
#else
			for (; j < (dim - k) - 1; j += 2) {
#endif
				__m128d v2_val = _mm_loadu_pd(&p[dim - k][j]);
				__m128d v2_2 = _mm_loadu_pd(&p[i][j]);
				__m128d v2;
				v2 = v2_1_SUB_p_i_dim_SUB_k_sqr * (vd2_1_0 - v2_val * v2_val);
				v2 = _mm_sqrt_pd(v2);
				v2 = (v2_2 - v2_p_i_dim_SUB_k * v2_val) / v2;
				_mm_storeu_pd(&p[i][j], v2);
				
				p[j+0][i] = p[i][j+0];
				p[j+1][i] = p[i][j+1];
#ifdef __AVX__
				j += 2;
#endif
			}
#endif

#if defined(__AVX__) || defined(__SSE2__)
			if (j < (dim - k)) {
#else
			for (; j < (dim - k); j++) {
#endif
				p[i][j] = p[j][i] = 
					(p[i][j] - p_i_dim_SUB_k * p[j][dim - k]) /
					(sqrt( (_1_SUB_p_i_dim_SUB_k_sqr) * (1 - pow2(p[j][dim - k])) ));
			}
		}
	}

	return p[0][1];
}

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

__attribute__((hot, always_inline)) static inline
bool IsNan(double val) {
#ifdef _MSC_VER
	return _isnan(val);
#else
	return std::isnan(val);
#endif
}

__attribute__((hot, always_inline)) static inline
void testAndRemove(double correlationCoefficient, Graph* __restrict__ g, const int r, const int c, const double ccRangeMin, const double ccRangeMax) {
	++g->score;
	
	if (IsNan(correlationCoefficient) || ((ccRangeMin <= correlationCoefficient) && (correlationCoefficient <= ccRangeMax))) {
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
__attribute__((hot, always_inline)) static inline
double getCorrelationCoefficient(const int* __restrict__ neighbours, const int* __restrict__ subset, const int l,
	Graph* __restrict__ g, const int r, const int c, double* __restrict__ * __restrict__ p) {
	double correlationCoefficient;
	
	if (l == 2) {
		double rhoRC, rhoRH1, rhoCH1;
		int h1 = neighbours[subset[0]];
		int h2 = neighbours[subset[1]];
		
		double result[2] __attribute__((aligned(16)));
		getCorrelationCoefficient2(g->rho[r][c], g->rho[r][h2], g->rho[c][h2], g->rho[r][h1], g->rho[r][h2], g->rho[h1][h2], result);
		rhoRC = result[0];
		rhoRH1 = result[1];
		rhoCH1 = getCorrelationCoefficient1(g->rho[c][h1], g->rho[c][h2], g->rho[h1][h2]);
		
		correlationCoefficient = getCorrelationCoefficient1(rhoRC, rhoRH1, rhoCH1);
	} else if (l == 3) {
		double rhoRC_H2H3, rhoRH1_H2H3, rhoCH1_H2H3, rhoRC_H3, rhoRH1_H3, rhoRH2_H3, rhoCH1_H3, rhoCH2_H3, rhoH1H2_H3;
		int h1 = neighbours[subset[0]];
		int h2 = neighbours[subset[1]];
		int h3 = neighbours[subset[2]];

		double result[4] __attribute__((aligned(32)));
		getCorrelationCoefficient4(g->rho[r][c], g->rho[r][h3], g->rho[c][h3], g->rho[r][h1], g->rho[r][h3], g->rho[h1][h3],
				g->rho[r][h2], g->rho[r][h3], g->rho[h2][h3], g->rho[c][h1], g->rho[c][h3], g->rho[h1][h3], result);
		rhoRC_H3 = result[0];
		rhoRH1_H3 = result[1];
		rhoRH2_H3 = result[2];
		rhoCH1_H3 = result[3];
		getCorrelationCoefficient2(g->rho[c][h2], g->rho[c][h3], g->rho[h2][h3], g->rho[h1][h2], g->rho[h1][h3], g->rho[h2][h3], result);
		rhoCH2_H3 = result[0];
		rhoH1H2_H3 = result[1];

		getCorrelationCoefficient2(rhoRC_H3, rhoRH2_H3, rhoCH2_H3, rhoRH1_H3, rhoRH2_H3, rhoH1H2_H3, result);
		rhoRC_H2H3 = result[0];
		rhoRH1_H2H3 = result[1];
		rhoCH1_H2H3 = getCorrelationCoefficient1(rhoCH1_H3, rhoCH2_H3, rhoH1H2_H3);
		
		correlationCoefficient = getCorrelationCoefficient1(rhoRC_H2H3, rhoRH1_H2H3, rhoCH1_H2H3);
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
__attribute__((hot, always_inline)) static inline
void iterativeComb(int* __restrict__ neighbours, const int neighboursDim, const int l, Graph* __restrict__ g, const int r,
	const int c, const double ccRangeMin, const double ccRangeMax, double* __restrict__ * __restrict__ p, int* __restrict__ currentCombination) {
	double coeff;

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
		testAndRemove(coeff, g, r, c, ccRangeMin, ccRangeMax);
	} while (g->matrix[r][c] && !((currentCombination[0] == (neighboursDim - 1 - l + 1)) && (currentCombination[l - 1] == (neighboursDim - 1))));
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
__attribute__((hot, always_inline)) static inline
void findAllSubsets(Graph* __restrict__ g, const int i, const int j, const int l, const double ccRangeMin, const double ccRangeMax,
		int* __restrict__ neighbours, double* __restrict__ * __restrict__ p, int* __restrict__ currentCombination) {
	int neighboursDim = 0;

	if (l == 0) {
		testAndRemove(g->rho[i][j], g, i, j, ccRangeMin, ccRangeMax);
	} else if (l == 1) {
		// find neighbours (when l > 0) of i without considering j
		double rho_i_j = g->rho[i][j];
		for (int k = 0; (k < g->nRows) && (neighboursDim < g->numNeighbours[i]); k++) {
			if (!g->matrix[i][j])
				break;
			if (g->matrix[i][k] && (k != j)) {
				neighboursDim++;
				double correlationCoefficient = getCorrelationCoefficient1(rho_i_j, g->rho[i][k], g->rho[j][k]);
				
				testAndRemove(correlationCoefficient, g, i, j, ccRangeMin, ccRangeMax);
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
		iterativeComb(neighbours, neighboursDim, l, g, i, j, ccRangeMin, ccRangeMax, p, currentCombination);
	}
}

__attribute__((hot, noinline)) static
void pcAlgorithm(Graph* __restrict__ g, double inverseCNDForAlpha, int* __restrict__ neighbours, double* __restrict__ * __restrict__ p,
		int* __restrict__ currentCombination)
{
	int l;
	bool hasWorked; // boolean to see that there is at least an arc i,j s.t. |adj(C, i)\{j}| >= l TO CHECK FROM THE TEXT
	
	// start with l = 0 & hasWorked = true
	l = -1;
	hasWorked = true;
	
	while ((hasWorked) && (l < g->nRows)) {
		l++;
		hasWorked = false;
		
		double expVal = exp(inverseCNDForAlpha / (sqrt(g->nCols - l - 3.0) * 0.5));
		
		double ccRangeMin = (1.0 - expVal) / (1.0 + expVal);
		double ccRangeMax = -ccRangeMin;
		
		for (int i = 0; i < g->nRows; i++) {
			if (g->numNeighbours[i] <= l)
				continue;
			for (int j = 0; j < g->nRows; j++) {
				// check if exists the arc between i and j
				if (g->matrix[i][j] && (g->numNeighbours[i] > l)) {
					hasWorked = true;
					findAllSubsets(g, i, j, l, ccRangeMin, ccRangeMax, neighbours, p, currentCombination);
				}
			}
		}
	}
}

/** Body of the PC algorithm.
	
	@see [TODO reference to the paper]
*/
int main(int argc, char* argv[]) {
	//int l;
	int experimentsDim; // number of experiment files to be read
	int hibridizationDim; // total number of columns of all experiment files
	int* tilesDim = NULL; // array that contains the length of each tiles read
	int tileRows; // number of rows of the tile.txt file 
	//bool hasWorked; // boolean to see that there is at least an arc i,j s.t. |adj(C, i)\{j}| >= l TO CHECK FROM THE TEXT
	intpair** tiles = NULL; // array that contains the numbers of the rows that will be extracted and represent the subgraph
	Graph* __restrict__ g = NULL; // the graph
	double alpha;
	int* __restrict__ neighbours = NULL;
	double* __restrict__ * __restrict__ p = NULL;
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
	double inverseCNDForAlpha;
	int* __restrict__ currentCombination;
	
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
	
	inverseCNDForAlpha = inverseComulativeNormalDistribution(1.0 - (alpha / 2.0));

	// istantiate the array of string
	experiments = new string*[experimentsDim];

	// read all the subgragraph that will be computed
	readTile(*tilePath, tilesDim, tiles, tileRows, experiments, experimentsDim);

	int numChkpt = readCheckpoint(*(outputFile), score);
	for (int c = numChkpt; c < tileRows; c++) {
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
			p[i] = DoubleVector::alloc(g->nRowsAligned);
			if (g->nRowsAligned > g->nRows) {
				std::memset(p[i] + g->nRows, 0, (g->nRowsAligned - g->nRows) * sizeof(double));
			}
		}

		currentCombination = new int[g->nRows];

		// alloc the neighbours array (save time)
		neighbours = new int[g->nRows];

		// PC-algorithm
		pcAlgorithm(g, inverseCNDForAlpha, neighbours, p, currentCombination);

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
		// originally score was incremented only every 1000 calls of testAndRemove
		score += g->score / 1000;
		
		// free the memory		
		delete[] neighbours;

		for (int i = 0; i < g->nRows; i++) {
			DoubleVector::free(p[i]);
		}
		delete[] p;
		
		delete[] currentCombination;

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
