#include "mex.h"
#include <math.h>
#include <string.h>

void padist(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]);
mxArray* getArraysof(mxArray* A, int i, int j, int nb);
mxArray* getArrayof(mxArray* A, int y,int number);
mxArray* constructArrayof(mxArray* A, int y, mxArray* inds);
mxArray* constructArraysof(mxArray* A, int y, int z, mxArray* inds);
double cosine(mxArray *A, int m, int n);
double eucdist(mxArray *A, mwSize m, mwSize n);
double norm(double *X,int n);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  /*  check for proper number of arguments */
    if (nrhs<3) {
        mexErrMsgIdAndTxt("stats:padistmex:TooFewInputs",
                          "Two input arguments required.");
    } else if(nlhs>1) {
        mexErrMsgIdAndTxt("stats:padistmex:TooManyOutputs",
                          "Too many output arguments.");
    }
     padist(nlhs, plhs, nrhs, prhs);
}

void padist(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *X,*S,*I;
    char *form, *distance;
    mxArray *dist,*sort,*result[2];
    int numCoords,numPoints;
    int rows, cols;
    int i,j;
    
    mxArray *rhs[3];
    //mxArray *duplc = mxDuplicateArray(prhs[0]); 
    rhs[0] = mxDuplicateArray(prhs[0]);
    rhs[1] = mxDuplicateArray(prhs[1]);
    rhs[2] = mxDuplicateArray(prhs[2]);
    
    //Distance string inlezen
    distance = (char*) mxCalloc(mxGetN(rhs[2])+1, sizeof(char)); //mxCalloc is similar to malloc in C
    mxGetString(rhs[2],distance,mxGetN(rhs[2])+1);
    
    //output matrix (m x n) initialiseren
    numPoints = mxGetM(rhs[0]); //het aantal rijen
    numCoords = mxGetN(rhs[0]); //het aantal kolommen
    int numDists = (numPoints * (numPoints-1)) / 2;
    if (numDists >= (double)MWSIZE_MAX) {
        mexErrMsgIdAndTxt("stats:padistmex:OutputTooLarge",
                          "Distance matrix has more elements than the maximum allowed size in MATLAB.");
    }
    plhs[0] = mxCreateNumericMatrix(1, numDists, mxGetClassID(prhs[0]), mxREAL);
    double* d = (double*) mxGetData(plhs[0]);
    //Matrix A klaarmaken
    //X = mxGetPr(prhs[0]);
    
    cols = mxGetN(rhs[0]); //aantal kolommen
    rows = mxGetM(rhs[0]); //aantal rijen
    
    //De matrix A sorteren, met S en I als resultaat
    int trans = 2;

    mxArray *input[2];
    mxArray *tmp = mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(tmp) = trans;
    input[0] = rhs[0];
    input[1] = tmp;
    mexCallMATLAB(2,result,2,input,"sort");
    mxArray* Sdup = mxDuplicateArray(result[0]);
    mxArray* Idup = mxDuplicateArray(result[1]);
   // S = mxGetPr(Sdup);
   // I = mxGetPr(Idup);
    double* num = mxGetPr(rhs[1]);
    int antl = (int) *num;
    for (i=0; i < rows; i++){
        for (j=i+1; j <rows; j++){
            mxArray* row1 = getArrayof(Idup,i,antl);
           mxArray* row2= getArrayof(Idup,j,antl);
           mxArray *inputs[2];
           mxArray *fin;
           inputs[0] = row1;
           inputs[1] = row2;
           mexCallMATLAB(1,&fin,2,inputs,"union");
           mxArray* inds = mxDuplicateArray(fin);
                           //inds ok, gecheckt
           mxArray* rows = constructArraysof(rhs[0],i,j,inds);
           double z;          
           if (strcmp(distance,"euclidean") == 0)
           {
               z = eucdist(rows,2,mxGetN(rows));
            }

           else if(strcmp(distance,"cosine") == 0)
           {
           z = cosine(rows,2,mxGetN(rows));
           }
           else printf("no matching distance");
          mxDestroyArray(rows);
          mxDestroyArray(inds);
          mxDestroyArray(rhs);
          *(d++) = z;
        }
    }
}

double cosine(mxArray *A, int m, int n)
{
    /*
     d = 1 - sum(XI.*XJ,2);   % Cosine & Corr & RankCorr
     */

    /* This actually calculates the dot product of pairs of vectors.  It
     * assumes that the data have been properly preprocessed:  ranked for
     * Spearman's rank correlation distance, normalized to zero mean for
     * both linear and rank correlation distances, and normalized to unit
     * length for all three distances.
     */
   int k;
   double   theSum;
   double *x =mxGetPr(A);
   double norm1,norm2;
   theSum = 0;
   for (k=0;k<n;k++){
        theSum += x[2*k]*x[2*k+1];
        norm1 += pow(x[2*k],2);
        norm2 += pow(x[2*k+1],2);         
   }
        norm1 = sqrt(norm1);
        norm2 = sqrt(norm2);
        double result = 1-theSum/(norm1*norm2);
        return result;
}

double norm(double *X,int n)
{
 int k;
 double sum = 0;
    for (k=0;k < n; k++){
        sum += X[k]*X[k];
    }
    return sqrt(sum);
}


double eucdist(mxArray *A, mwSize m, mwSize n)
{
    /*
     d = sqrt(sum((XI-XJ).^2,2));            % Euclidean
     */
    int i,j,k;
    double   theSum,Y;
    double   *XI, *XJ, *XI0;
    double *x = mxGetPr(A);
    XI = x;
    XJ = XI+n;

    theSum = 0;
    for (k=0;k<n;k++){
        Y = XI[k]-XJ[k];
        theSum += Y*Y;
    }
    return sqrt(theSum);  
}
mxArray* constructArraysof(mxArray* A, int y, int z, mxArray* inds){
    double *X, *ind;
    int listsize;
    X = mxGetPr(A);
    ind = mxGetPr(inds);
     listsize = mxGetN(inds);
     mxArray *result = mxCreateDoubleMatrix(2, listsize, mxREAL);
     double *ptn = mxGetPr(result);
     int cols, rows, i, j;
     cols = mxGetN(A); //aantal kolommen
     rows = mxGetM(A); //aantal rijen
     j = 0;
     for (i = 0; i < listsize; i++){
          ptn[2*j] = X[(int)(ind[i]-1)*rows+y];
          ptn[2*j+1] = X[(int)(ind[i]-1)*rows+z];
         j++;
     }  
     return result;
}

mxArray* constructArrayof(mxArray* A, int y, mxArray* inds){
    double *X, *ind;
    int listsize;
    X = mxGetPr(A);
    ind = mxGetPr(inds);
    listsize = mxGetN(inds);
    mxArray *result = mxCreateDoubleMatrix(1, listsize, mxREAL);
    double *ptn = mxGetPr(result);
    int cols, rows, i, j;
   // cols = mxGetN(A); //aantal kolommen
    rows = mxGetM(A); //aantal rijen
    j = 0;
    for (i = 0; i < listsize; i++){
       ptn[j] = X[(int)(ind[i]-1)*rows+y];
       j++;
    }
    return result;
}

mxArray* getArrayof(mxArray* A, int y,int number){
    double *X;
    X = mxGetPr(A);
    mxArray *result = mxCreateDoubleMatrix(1, number, mxREAL);
     double *ptn = mxGetPr(result);
     int cols, rows, i, j;
     cols = mxGetN(A); //aantal kolommen
     rows = mxGetM(A); //aantal rijen
     j = 0;
     for (i = cols-number; i < cols; i++){
        ptn[j] = X[i*rows+y];
        j++;
     }   
    return result;
}


mxArray* getArraysof(mxArray* A, int y, int z,int number){
    double *X;
    X = mxGetPr(A);
    mxArray *result = mxCreateDoubleMatrix(2, number, mxREAL);
    double *ptn = mxGetPr(result);
    int cols, rows, i, j;
    cols = mxGetN(A); //aantal kolommen
    rows = mxGetM(A); //aantal rijen
    j = 0;
    for (i = cols-number; i < cols; i++){
        ptn[2*j] = X[i*rows+y];
        ptn[2*j+1] = X[i*rows+z];
       j++;
    }
    return result;
}
