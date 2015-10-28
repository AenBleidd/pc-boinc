/**
	Copyright (c) 2013, All Right Reserved
	
	This software is in the public domain, furnished "as is", without technical
	support, and with no warranty, express or implied, as to its usefulness for
	any purpose.
	
	Graph.cpp
	File defining the class Graph, representing the graph/network in which we
	will run the PC algorithm for causality discovering
	
	University of Trento,
	Department of Information Engineering and Computer Science
	
	Trento, fall 2013 - spring 2014

	Authors: (alphabetically ordered) Francesco Asnicar, Luca Masera,
	         Paolo Morettin, Nadir Sella, Thomas Tolio.
*/

#include <cmath>
#include <cstdlib>
#include <algorithm>
#ifndef _GRAPH
#include "Graph.hpp"
#endif

/**	Constructor method.
	
	The graph is implemented as a matrix of edges. Each cell has coordinates
	(x,y), representing the edge x->y. If the edge exists within the graph, the
	correspondent cell contains the value TRUE, otherwise FALSE.
	
	@param const int dim
	The dimension with which the matrix will be built.
*/
Graph::Graph(const int dim) {
	nRows = dim;

	numNeighbours = new int[nRows];

	// create the bool matrix
	matrix = new bool*[nRows];

	for (int i = 0; i < nRows; i++) {
		matrix[i] = new bool[nRows];
	}

	means = new double[nRows];
	probeIDs = new std::string[nRows];
	standardDeviations = new double[nRows];

	// create rho matrix
	rho = new double*[nRows];

	for (int i = 0; i < nRows; i++) {
		rho[i] = new double[nRows];
	}

	score = 0;
	count_score = 1;

	// initialize matrix, numNeighbours and l
	initializeMatrix(matrix, nRows);
	initializeNeighbours(numNeighbours, nRows);
	initializeZero(means, nRows);
	initializeZero(standardDeviations, nRows);
}

/** Destructor method.
	
	The method simply deallocates all the previously allocated data.
*/
Graph::~Graph(void) {
	// empty the memory for matrix
	for (int i = 0; i < nRows; i++) {
		delete[] matrix[i];
	}
	delete[] matrix;

	// empty the memory for bioData
	for (int i = 0; i < nRows; i++) {
		delete[] bioData[i];
	}
	delete[] bioData;

	// empty the memory for means
	delete[] means;

	// empty the memory for probeIDs
	delete[] probeIDs;

	// empty the memory for standardDeviations
	delete[] standardDeviations;

	// empty the memory for numNeighbours
	delete[] numNeighbours;

	// empty the memory for rho
	for (int i = 0; i < nRows; i++) {
		delete[] rho[i];
	}
	delete[] rho;
}

/** Shuffles the given probes.

	@param const time_t seed
	Seed for shuffling univocally the probes passed as input.

	@return nothing.
*/
void Graph::shuffleInputProbes(const time_t seed) {

	std::srand(seed);
	std::random_shuffle(&bioData[0], &bioData[nRows]);

	std::srand(seed);
	std::random_shuffle(&means[0], &means[nRows]);

	std::srand(seed);
	std::random_shuffle(&probeIDs[0], &probeIDs[nRows]);
}

/** Compute the standard deviations for each node in the graph.

	@return nothing.
*/
void Graph::computeStandardDeviations(void) {
	for (int r = 0; r < nRows; r++) {
		for (int c = 0; c < nCols; c++) {
			standardDeviations[r] += pow((bioData[r][c] - means[r]), 2);
		}

		standardDeviations[r] /= (double) nCols;
		standardDeviations[r] = sqrt(standardDeviations[r]);
	}
}

/** Compute the correlation coefficient of the base case, and store it in rho.

	@return nothing.
*/
void Graph::computeCorrelations(void) {
	double covariance = 0.0;

	for (int i = 0; i < nRows; i++) {
		for (int j = 0; j < nRows; j++) {
			covariance = 0.0;

			for (int k = 0; k < nCols; k++) {
				covariance += (bioData[i][k] - means[i]) * (bioData[j][k] - means[j]);
			}

			// divide covariance by nCols
			covariance /= nCols;

			// covariance(i, j) / sigma(i) * sigma(j)
			rho[i][j] = covariance / (standardDeviations[i] * standardDeviations[j]);
		}
	}
}

/** Initialize the boolean matrix to TRUE, but the diagonal, set to FALSE.
	
	@param bool** matrix
	Matrix to initialize.

	@param const int dim
	Dimension (number of rows, and then also number of columns) of the matrix.

	@return nothing.
*/
void Graph::initializeMatrix(bool** matrix, const int dim) {
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			if (i != j) {
				matrix[i][j] = true;
			} else {
				matrix[i][j] = false;
			}
		}
	}
}

/** Initialize the array numNeighbours with the value dim - 1, since the
	initial graph is connected.
	
	@param int* numNeighbours
	List of the number of neighbours for each node of the graph.
	
	@param const int dim
	The dimension of the list of neighbours.
	
	@return nothing.
*/
void Graph::initializeNeighbours(int* numNeighbours, const int dim) {
	for (int i = 0; i < dim; i++) {
		numNeighbours[i] = dim - 1;
	}
}

/** Initialize the given array till dimension to 0.0.

	@param double* array
	[TODO]

	@param const int dim
	The dimension of the array of probes.

	@return nothing.
*/
void Graph::initializeZero(double* array, const int dim) {
	for (int i = 0; i < dim; i++) {
		array[i] = 0.0;
	}
}
