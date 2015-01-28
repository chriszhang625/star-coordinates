#pragma once

#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <complex>
#include <string>
#include <map>
#include <QString>
#include <QFile>
#include <qpoint.h>
#include <ctime>
#include <assert.h>
#include <QTextStream>
#include <QStringList>
#include <qbytearray.h>
#include <dataanalysis.h>
#include <ap.h>
#include <alglibinternal.h>
#include <alglibmisc.h>
//#include <diffequations.h>
//#include <fasttransforms.h>
//#include <integration.h>
//#include <interpolation.h>
#include <linalg.h>
#include <optimization.h>
#include <solvers.h>
#include <specialfunctions.h>
#include <statistics.h>
#include <stdafx.h>

#define PI 3.1415926535
#define MAX_CHAR   128
class mathlib
{
public:
	mathlib(void);
    ~mathlib(void);
	void whiteNoise(std::vector<float> &sequence,int n);
	void CalculatePCACoe(std::vector<std::vector<float>> &data,std::vector<std::vector<float>> &coematrix);
	void CalLDACoe(std::vector<std::vector<float>> &data,std::vector<std::vector<float>> &coematrix,int subnum,int classnum);
	void CalTRLDACoe(std::vector<std::vector<float>> &data,std::vector<std::vector<float>> &coematrix,std::vector<int> labels,int subdim,int classnum,bool is_Descend);
	bool CalCholDecomp(std::vector<std::vector<float>> &data, std::vector<std::vector<float>> &triangular,bool is_Upper);
	bool CalSEVD(std::vector<std::vector<float>> input, std::vector<float> &eigvalues, std::vector<std::vector<float>> &eigvectors);
	float CalNorm(std::vector<float> input);
	float FindtheMax(std::vector<float> inputvector,int &index);
	float FindtheMin(std::vector<float> inputvector,int &index);
	float n_hundred(int n);
	float digit_process(float t, int n);
	
	void loadDatasetFromFile(std::vector<std::vector<float>> &dataset,const QString dataSetFilename,std::vector<QString> &labels);
	void MatrixTransition(std::vector<std::vector<float>> &dataset1,std::vector<std::vector<float>> &dataset2);
	void MatrixSynchronization(std::vector<std::vector<float>> &dataset1,std::vector<std::vector<float>> &dataset2,std::vector<int> index);
	void MatrixInsertion(std::vector<std::vector<float>> &dataset,std::vector<float> input, int index);
	bool MatrixInversion(std::vector<std::vector<float>> &data,std::vector<std::vector<float>> &data_inverse);
	void MatrixScaling(std::vector<std::vector<float>> &data, float scale);
	std::vector<float> MatrixScaling(std::vector<float> data, float scale);
	void MatrixAbs(std::vector<std::vector<float>> &data);
	void MatrixAbs(std::vector<float> &data);
	std::vector<float> MatrixDiag(std::vector<std::vector<float>> input,int offset);
	float MatrixTrace(std::vector<float> diag);
	void StandardAxes(std::vector<std::vector<float>> &result, int numberdim);
	void axeScaling(std::vector<std::vector<float>> &dataset,float scalenum);
	void quicksort(std::vector<float> &v, std::vector<int> &index,int leftnum,int rightnum);
	void quicksort(std::vector<int> &v, std::vector<int> &index,int leftnum,int rightnum);
	void clearVector(std::vector<std::vector<float>> &vt);
	void vectorSort(std::vector<float>& index, std::vector<std::vector<float>>& input);
	void vectorSort(std::vector<int> &index, std::vector<float> &input,bool is_Descend);

	std::vector<int> diskmeans(std::vector<std::vector<float>> &data,std::vector<std::vector<float>> &coematrix,int classnum);
	std::vector<float> linspace(int a,int b,int size);
	std::vector<float> CalMean(std::vector<std::vector<float>> input);
	std::vector<float> randperm(int numbers,int selects);
	std::vector<int> kmeans(std::vector<std::vector<float>> input,int k);
	std::vector<float> CalUnique(std::vector<float> input);
	std::vector<int> CalUnique(std::vector<int> &input);
	std::vector<std::vector<float>> gempairs(std::vector<std::vector<float>> &data,std::vector<int> label, float gamma, float threshold, int classnum);
	std::vector<std::vector<float>> MatrixMultiply(std::vector<std::vector<float>> data_1,std::vector<std::vector<float>> data_2);
	std::vector<std::vector<float>> MatrixSum(std::vector<std::vector<float>> data_1, std::vector<std::vector<float>> data_2,bool is_Add);
	std::vector<std::vector<float>> MatrixEye(int rows,int column);
	std::vector<std::vector<float>> MatrixDeletion(std::vector<std::vector<float>> input, int index);
	std::vector<std::vector<float>> MatrixSqrt(std::vector<std::vector<float>> input);
	std::vector<std::vector<float>> MatrixMax(std::vector<std::vector<float>> data_1,std::vector<std::vector<float>> data_2);
	std::vector<std::vector<float>> MatrixDivide(std::vector<std::vector<float>> data_1,std::vector<std::vector<float>> data_2);
	std::vector<float> MatrixSum(std::vector<float> data_1,std::vector<float> data_2,bool is_Add);
	std::vector<std::vector<float>> vector_selection(std::vector<std::vector<float>> A, std::vector<std::vector<float>> B, std::vector<std::vector<float>> eigvector, int sub_dim,bool isMax);
	std::vector<std::vector<float>> vectorMultiply(std::vector<float> data_1,std::vector<float> data_2);
	std::vector<std::vector<float>> TraceRatioGeneral(std::vector<std::vector<float>> A, std::vector<std::vector<float>> B, std::vector<std::vector<float>> C, int sub_dim,bool isMax);
};

