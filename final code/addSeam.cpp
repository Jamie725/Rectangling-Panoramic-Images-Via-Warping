#include <mex.h>
#include <cstdio>
#include <cstring>
#include <climits>
#include <algorithm>
#include <omp.h>
using namespace std;

unsigned char* img;
int w, h;
int* mask;
int* outMask;
unsigned char* outImg;
int* xDispMap;
int* yDispMap;
int* outXDispMap;
int* outYDispMap;

void readPPM(const char* name, unsigned char** img, int* w, int* h) {
	FILE* file = fopen(name, "rb");
	fscanf(file, "P6\n%d %d\n255\n", w, h);
	int size = (*w) * (*h) * 3;
	*img = new unsigned char[size];
	fread(*img, sizeof(unsigned char), size, file);
	fclose(file);
}
void writePPM(const char* name, const unsigned char* img, int w, int h) {
	FILE* file = fopen(name, "wb");
	fprintf(file, "P6\n%d %d\n255\n", w, h);
	int size = w * h * 3;
	fwrite(img, sizeof(unsigned char), size, file);
	fclose(file);
}
void writePGM(const char* name, const unsigned char* img, int w, int h) {
	FILE* file = fopen(name, "wb");
	fprintf(file, "P5\n%d %d\n255\n", w, h);
	int size = w * h;
	fwrite(img, sizeof(unsigned char), size, file);
	fclose(file);
}
void energy(const unsigned char* image, int W, int H, int* energy, int* mask){
	memset(energy, 0, sizeof(energy[0])*W*H);
	//int addValue = max(W,H)*3*10;
	int addValue = 1000000;
	#pragma omp parallel for
	for(int i=0; i<H; ++i)
		for(int j=0; j<W-1; ++j)
			for(int k=0; k<3; ++k)
				energy[i*W+j] += abs((int)image[(i*W+(j+1))*3+k] - (int)image[(i*W+j)*3+k]);
	#pragma omp parallel for
	for(int i=0; i<H-1; ++i)
		for(int j=0; j<W; ++j)
			for(int k=0; k<3; ++k)
				energy[i*W+j] += abs((int)image[((i+1)*W+j)*3+k] - (int)image[(i*W+j)*3+k]);
	#pragma omp parallel for
	for(int i=0; i<H; ++i)
		for(int j=0; j<W; ++j)
			energy[i*W+j] += mask[i*W+j] * addValue;
}
///////////////////////Vertical Seam////////////////////////////
void accumulateVerEn(int * energy, int W, int H, int* verticalEnergy, int* verticalRoute){
	/////////Initialize////////////
	memset(verticalEnergy, 0, sizeof(verticalEnergy[0])*W*H);
	memset(verticalRoute, 0, sizeof(verticalRoute[0])*W*H);
	#pragma omp parallel for
	for(int i=0; i<W; ++i)
		verticalEnergy[i] = energy[i];

	////////Accumulte the min energy and find the route//////////////
	for(int i=1; i<H; ++i){
		for(int j=0; j<W; ++j){
			if(j-1<0){
				verticalEnergy[i*W+j] = energy[i*W+j] + min(verticalEnergy[(i-1)*W+j], verticalEnergy[(i-1)*W+j+1]);//e(i,j) + min( M(i-1,j) , M(i-1, j+1))
				verticalRoute[i*W+j] = (verticalEnergy[(i-1)*W+j]<=verticalEnergy[(i-1)*W+j+1]) ? 0 : 1;
			}
			else if(j+1>=W){
				verticalEnergy[i*W+j] = energy[i*W+j] + min(verticalEnergy[(i-1)*W+j-1], verticalEnergy[(i-1)*W+j]);//e(i,j) + min( M(i-1,j-1) , M(i-1, j))
				verticalRoute[i*W+j] = (verticalEnergy[(i-1)*W+j-1]<verticalEnergy[(i-1)*W+j]) ? (-1) : 0;
			}
			else{
				verticalEnergy[i*W+j] = energy[i*W+j] + min(min(verticalEnergy[(i-1)*W+j-1], verticalEnergy[(i-1)*W+j]), verticalEnergy[(i-1)*W+j+1]);//e(i,j) + min( M(i-1,j-1) , M(i-1, j), M(i-1, j+1))
				verticalRoute[i*W+j] = (verticalEnergy[(i-1)*W+j-1]<verticalEnergy[(i-1)*W+j]) ? ((verticalEnergy[(i-1)*W+j-1]<=verticalEnergy[(i-1)*W+j+1]) ? (-1) : 1) : ((verticalEnergy[(i-1)*W+j]<=verticalEnergy[(i-1)*W+j+1]) ? 0 : 1);
			}
		}
	}
}
void optimalVertSeam(int * energy, int W, int H, int* optVertSeam){
	/////////Initialize//////////
	memset(optVertSeam, 0, sizeof(optVertSeam[0])*H);
	int* verticalEnergy = new int [W*H];
	int* verticalRoute = new int [W*H];
	int minimum = INT_MAX;
	int minIndex =  INT_MAX;

	////////Find minimum seam/////////
	accumulateVerEn(energy, W, H, verticalEnergy, verticalRoute);
	for(int i=0; i<W; ++i){
		if(minimum>verticalEnergy[(H-1)*W+i]){
			minimum = verticalEnergy[(H-1)*W+i];
			minIndex = i;
		}
	}

	//////////Backtrack seam route/////////
	int lastIndex = (H-1)*W + minIndex;
	optVertSeam[H-1] = lastIndex;
	for(int i=H-2; i>=0; --i ){
		optVertSeam[i] = lastIndex - W + verticalRoute[lastIndex];
		lastIndex = optVertSeam[i];
	}
	delete[] verticalEnergy;
	delete[] verticalRoute;
}
void addVertSeam(int* optVertSeam, unsigned char * img, int* xDispMap, int* yDispMap,int* mask, int W, int H){
	#pragma omp parallel for
	for(int i=0; i<H; ++i){
		int seamPosition = 0, column = 0;
		seamPosition = optVertSeam[i]*3;
		column = optVertSeam[i] % W;
		////////////Shift//////////////
		for(int j=W-1; j>column; --j){
			//--displacementMap[i*W+j];
			xDispMap[i*W+j] = xDispMap[i*W+(j-1)];						
			yDispMap[i*W+j] = yDispMap[i*W+(j-1)];						
			mask[i*W+j] = mask[i*W+(j-1)];
			for(int k=0; k<3; ++k)
				img[(i*W+j)*3+k] = img[(i*W+(j-1))*3+k];
		}
		///////////Add seam/////////////
		mask[optVertSeam[i]] = 0;
		if(column-1<0){
			for(int j=0; j<3; ++j)
				img[seamPosition+j] = img[seamPosition+3+j];
		}
		else if(column +1 >= W){
			for(int j=0; j<3; ++j)
				img[seamPosition+j] = img[seamPosition-3+j];
		}
		else{
			for(int j=0; j<3; ++j){
				img[seamPosition+j] = ((double)img[seamPosition+3+j] + (double)img[seamPosition-3+j])/2 + 0.5;
			}
		}
	}
}
void main() {
	//////////////Initialize/////////////
	//readPPM("sakura.ppm", &img, &w, &h);
	int * energyMap = new int [w*h];
	int * optVertSeam = new int [h];

	////////////Find the energy map//////////////
	energy(img, w,h, energyMap, mask);

	////////////Find optimal vertical energy///////////
	optimalVertSeam(energyMap, w,  h, optVertSeam);
#if 1
	// Check if seam lies on background
	bool oops = false;
	for (int i = 0; i < h; ++i) {
		if (mask[optVertSeam[i]] != 0) {
			oops = true;
			break;
		}
	}
	if (oops) {
		printf("******** Warning: Seam passes background area !!! ********\n");
	}
#endif

	///////////Add vertical seam////////////////
	addVertSeam(optVertSeam, img, xDispMap, yDispMap, mask, w,  h);
int color[3];

#if 1
	color[0] = rand()%256;
	color[1] = rand()%256;
	color[2] = rand()%256;
	// Draw red on the seam
	for (int i = 0; i < h; ++i) {		
		img[optVertSeam[i]*3] = color[0];
		img[optVertSeam[i]*3+1] = color[1];
		img[optVertSeam[i]*3+2] = color[2];
	}
#endif
#if 0
	// Draw red on the seam
	for (int i = 0; i < h; ++i) {		
		img[optVertSeam[i]*3] = 255;
		img[optVertSeam[i]*3+1] = 0;
		img[optVertSeam[i]*3+2] = 0;
	}
#endif
	///////////Output image////////////////
	memcpy(outImg, img, sizeof(outImg[0])*w*h*3);
	memcpy(outMask, mask, sizeof(outMask[0])*w*h);
	memcpy(outXDispMap, xDispMap, sizeof(outXDispMap[0])*w*h);
	memcpy(outYDispMap, yDispMap, sizeof(outYDispMap[0])*w*h);
	delete[] energyMap;
	delete[] optVertSeam;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	/// column major!
	// input
	img = (unsigned char*)mxGetData(prhs[0]);
	int nDims = mxGetNumberOfDimensions(prhs[0]);
	const int* dims = mxGetDimensions(prhs[0]);
	w = dims[1];
	h = (nDims == 3)? dims[2] : 1;
	//printf("Sub image size = %dx%d\n", w, h);
	mask 	 = (int*)mxGetData(prhs[1]);
	xDispMap = (int*)mxGetData(prhs[2]);
	yDispMap = (int*)mxGetData(prhs[3]);

	// output
	/////output: outImg
	plhs[0] = mxCreateNumericArray(nDims, dims, mxUINT8_CLASS, mxREAL);
	outImg = (unsigned char*)mxGetData(plhs[0]);
	/////output: mask
	plhs[1] = mxCreateNumericMatrix(w, h, mxINT32_CLASS, mxREAL);
	outMask = (int*)mxGetData(plhs[1]);
	/////output: displacementMap
	plhs[2] = mxCreateNumericMatrix(w, h, mxINT32_CLASS, mxREAL);
	outXDispMap = (int*)mxGetData(plhs[2]);
	plhs[3] = mxCreateNumericMatrix(w, h, mxINT32_CLASS, mxREAL);
	outYDispMap = (int*)mxGetData(plhs[3]);

	// main function
	main();
}
