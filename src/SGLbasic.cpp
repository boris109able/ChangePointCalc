#include <iostream>
#include <math.h>
#include <stdlib.h>


using namespace std;

extern "C" {

 int single = 0;

  //////////////////////////////////

void gradCalc(int *nrow, double *eta, double *y, double *ldot)
{
  for(int i = 0; i < nrow[0]; i++)
  {
    ldot[i] = 2*(eta[i] - y[i])/nrow[0];
  }
}


double squareSumDiff(int *nrow, double *eta, double *y)
{
  double squareSum = 0;
  for(int i = 0; i < nrow[0]; i++)
  {
    squareSum = squareSum + pow(eta[i] - y[i], 2); 
  }
  return squareSum/nrow[0];
}


void solveOneGroup(double *X, double *y,  int *nrow, int *ncol,  double *beta, double *lambda1, double *lambda2, int *innerIter, double *thresh, double *ldot, double *gamma, double *eta, int* betaIsZero, int& groupChange, int* isActive, int* useGroup, double *step, int *reset)
{
	int n = nrow[0];
	int p = ncol[0];
	double zeroCheck = 0;
	double check = 0;
	int count = 0;
	double t = step[0];
	double diff = 1;
	double norm = 0;
	double uOp = 0;
	double Lnew = 0;
	double Lold = 0;
	double sqNormG = 0;
	double iProd = 0;
	double *etaNew = NULL;
	etaNew = new double[n];
	double *etaNull = NULL;
	etaNull = new double[n];
	int groupStart = 0;
	//  int reset = 20;  
	//int i = 0;
	for (int i = 0; i<n; i++)
		//for(int irand = 0; irand < n; irand++)
	{
		//i = (int) floor(rand() / (RAND_MAX+0.0) * (nrow[0]+0.0));
		groupStart = i*p;
		if (useGroup[i] == 1)
		{
			//startInd = rangeGroupInd[i];    
			// Setting up null gradient calc to check if group is 0
			for (int k = 0; k < n; k++)
			{
				etaNull[k] = eta[k];
				if (k >= i)
				{
					for (int j = 0; j < p; j++)
					{
						etaNull[k] = etaNull[k] - X[k + nrow[0] * j] * beta[j + groupStart];
					}
				}
			}
			// Calculating Null Gradient
			gradCalc(nrow, etaNull, y, ldot);
			double *grad = NULL;
			grad = new double[p];
			for (int j = 0; j < p; j++)
			{
				grad[j] = 0;
				for (int k = i; k < n; k++)
				{
					grad[j] = grad[j] + X[k + n * j] * ldot[k];
				}
				if (grad[j] <= lambda1[0] && grad[j] >= -lambda1[0])
				{
					grad[j] = 0;
				}
				else if (grad[j] > lambda1[0])
				{
					grad[j] = grad[j] - lambda1[0];
				}
				else if (grad[j] < -lambda1[0])
				{
					grad[j] = grad[j] + lambda1[0];
				}
			}
			//now grad contains S()
			//check if the current group is zero
			zeroCheck = 0;
			for (int j = 0; j < p; j++)
			{
				zeroCheck = zeroCheck + pow(grad[j], 2);
			}
			if (zeroCheck <= pow(lambda2[0], 2))
			{
				if (betaIsZero[i] == 0)
				{
					for (int k = i; k < n; k++)
					{
						for (int j = 0; j < p; j++)
						{
							eta[k] = eta[k] - X[k + nrow[0] * j] * beta[j + groupStart];
						}
					}
				}
				betaIsZero[i] = 1;
				for (int j = 0; j < p; j++)
				{
					beta[j + groupStart] = 0;
				}
			}
			else
			{
				if (isActive[i] == 0)
				{
					groupChange = 1;
				}
				isActive[i] = 1;
				double *theta = new double[p];
				for (int j = 0; j < p; j++)
				{
					theta[j] = beta[j + groupStart];//store old beta values for the group
				}
				betaIsZero[i] = 0;
				double *z = NULL;
				z = new double[p];
				double *U = NULL;
				U = new double[p];
				double *G = NULL;
				G = new double[p];
				count = 0;
				check = 1000;
				//while(count <= innerIter[0] && check > thresh[0])
				while (check > thresh[0])
				{
					count++;
					gradCalc(nrow, eta, y, ldot);
					for (int j = 0; j < p; j++)
					{
						grad[j] = 0;
						for (int k = i; k < nrow[0]; k++)
						{
							grad[j] = grad[j] + X[k + n * j] * ldot[k];
						}
					}
					//t is already initialized to be step which is 1 by default
					Lold = squareSumDiff(nrow, eta, y);
					diff = -1;
					//back-tracking
					while (diff < 0)
					{
						for (int j = 0; j < p; j++)
						{
							z[j] = beta[j + groupStart] - t * grad[j];
							if (z[j] <= lambda1[0] * t && z[j] >= -lambda1[0] * t)
							{
								z[j] = 0;
							}
							else if (z[j] > lambda1[0] * t)
							{
								z[j] = z[j] - lambda1[0] * t;
							}
							else if (z[j] < -lambda1[0] * t)
							{
								z[j] = z[j] + lambda1[0] * t;
							}
						}
						norm = 0;
						for (int j = 0; j < p; j++)
						{
							norm = norm + pow(z[j], 2);
						}
						norm = sqrt(norm);
						if (norm != 0) { uOp = (1 - t*lambda2[0] / norm); }
						else { uOp = 0; }
						if (uOp < 0) { uOp = 0; }
						for (int j = 0; j < p; j++)
						{
							U[j] = uOp*z[j];
							G[j] = 1 / t *(beta[j + groupStart] - U[j]);
						}
						// Setting up betaNew and etaNew in direction of Grad for descent step
						for (int k = 0; k < n; k++)
						{
							etaNew[k] = eta[k];
							if (k >= i)
							{
								for (int j = 0; j < p; j++)
								{
									etaNew[k] = etaNew[k] - t*G[j] * X[k + n*j];
								}
							}
						}
						Lnew = squareSumDiff(nrow, etaNew, y);
						sqNormG = 0;
						iProd = 0;
						for (int j = 0; j < p; j++)
						{
							sqNormG = sqNormG + pow(G[j], 2);
							iProd = iProd + grad[j] * G[j];
						}
						diff = Lold - Lnew - t * iProd + t / 2 * sqNormG;
						t = t * gamma[0];
					}
					t = t / gamma[0];//reverse the last line effect in the while loop
					check = 0;
					for (int j = 0; j < p; j++)
					{
						check = check + fabs(theta[j] - U[j]);
						for (int k = i; k < n; k++)
						{
							eta[k] = eta[k] - X[k + n*j] * beta[j + groupStart];
						}
						//beta[j + rangeGroupInd] = U[j] + count%reset[0]/(count%reset[0]+3) * (U[j] - theta[j + rangeGroupInd[i]]);
						beta[j + groupStart] = theta[j] + count / (count + 3.0) * (U[j] - theta[j]);
						theta[j] = U[j];
						for (int k = i; k < n; k++)
						{
							eta[k] = eta[k] + X[k + n*j] * beta[j + groupStart];
						}
					}
				}
				delete[] theta;
				delete[] z;
				delete[] U;
				delete[] G;
			}
			delete[] grad;
		}
	}
	delete[] etaNew;
	delete[] etaNull;
	//delete [] theta;
}
int solveSGL(double *X, double* y, int* nrow, int* ncol, double *lambda1, double *lambda2, double *beta, int *innerIter, int *outerIter, double *thresh, double *outerThresh, double *eta, double *gamma, int *betaIsZero, double *step, int *reset)
{
  int np = nrow[0]*ncol[0];
  int n = nrow[0];
  //int p = ncol[0];
  double *ldot = NULL;
  ldot = new double[n];
  int groupChange = 1;
  int* isActive = NULL;
  isActive = new int[n];
  int* useGroup = NULL;
  useGroup = new int[n];
  int* tempIsActive = NULL;
  tempIsActive = new int[n];  
  for(int i = 0; i < n; i++)
  {
    isActive[i] = 0;
    useGroup[i] = 1;
  }
  // outer most loop creating response etc...
  int outermostCounter = 0;
  double outermostCheck = 100000;
  double* outerOldBeta = NULL;
  outerOldBeta = new double[np];
  while(groupChange == 1)
  {
    groupChange = 0;
    solveOneGroup(X, y, nrow, ncol,  beta,  lambda1, lambda2, innerIter, thresh, ldot, gamma, eta, betaIsZero, groupChange, isActive, useGroup, step, reset);
    while(outermostCounter < outerIter[0] && outermostCheck > outerThresh[0])
    //while(outermostCheck > outerThresh[0])
    {
      outermostCounter ++;
      for(int i = 0; i < np; i++)
	    {
	      outerOldBeta[i] = beta[i];
	    }
      for(int i = 0; i < n; i++)
	    {
	      tempIsActive[i] = isActive[i];
	    }
      solveOneGroup(X, y, nrow, ncol, beta, lambda1, lambda2, innerIter, thresh, ldot, gamma, eta, betaIsZero, groupChange, isActive, tempIsActive, step, reset);
	    outermostCheck = 0;
      for(int i = 0; i < np; i++)
	    {
	      outermostCheck = outermostCheck + fabs(outerOldBeta[i] - beta[i]);
	    }
    }
  }
  delete [] outerOldBeta;
  delete [] ldot;
  delete [] isActive;
  delete [] useGroup;
  delete [] tempIsActive;
  return 1;
}

//***************************************************************************************************************************************************************************
//logistic regression
 void logitGradCalc(int *nrow, double *prob, int *y, double *ldot)
  {
    for(int i = 0; i < nrow[0]; i++)
      {
	ldot[i] = (-y[i] + prob[i])/nrow[0]; /* OR MAYBE NOT? */
      }
  }
  
 void pCalc(int *nrow, double *eta, double *prob)
 {
	 for (int i = 0; i < nrow[0]; i++)
	 {
		 //prob[i] = exp(eta[i]) / (1 + exp(eta[i]));
		 prob[i] = 1 / (1 + exp(-eta[i]));
	 }
 }
  
  
  double logitNegLogLikelihoodCalc(int *nrow, double *prob, int *y)
  {
	  double logLik = 0;
	  for (int i = 0; i < nrow[0]; i++)
	  {
		  logLik = logLik + y[i] * log(prob[i]) + (1 - y[i]) * log(1 - prob[i]);
	  }

	  return -logLik / nrow[0];  /* OR MAYBE NOT? */
  }
  
  void betaZeroSolve(int *nrow, double *betaZero, double *eta, double *prob, double *thresh, int *innerIter, int *y)//numerically solve beta_0, the constant shift.
  {
	  double diff = 10;
	  double num = 0;
	  double denom = 0;
	  int count = 0;

	  while (pow(diff, 2) > pow(thresh[0], 2) && count < innerIter[0])
	  {
		  pCalc(nrow, eta, prob);
		  diff = 0;

		  for (int i = 0; i < nrow[0]; i++)
		  {
			  num += y[i] - prob[i];
			  denom += prob[i] * (1 - prob[i]);
		  }
		  diff = num / denom;
		  betaZero[0] = betaZero[0] - diff;

		  for (int i = 0; i < nrow[0]; i++)
		  {
			  eta[i] = eta[i] + diff;
		  }
	  }
  }


  void logitSolver(double *X, int *y, int *nrow, int *ncol, double *beta, double *lambda1, double *lambda2, int *innerIter, double *thresh, double *ldot, double *gamma, double *eta, int* betaIsZero, int& groupChange, int* isActive, int* useGroup, double *prob, double *betaZero, double *step)
{
	int n = nrow[0];
	int p = ncol[0];
	//int np = n*p;
	double *theta = new double[p];
	//double *thetaNew = new double[p];
	//int startInd = 0;
	double zeroCheck = 0;
	double check = 0;
	int count = 0;
	double t = step[0];
	double diff = 1;
	double norm = 0;
	double uOp = 0;
	double Lnew = 0;
	double Lold = 0;
	double sqNormG = 0;
	double iProd = 0;
	double *etaNew = NULL;
	etaNew = new double[n];
	double *etaNull = NULL;
	etaNull = new double[n];
	int groupStart = 0;
	for (int i = 0; i < n; i++)
    //int i;
    //for(int irand = 0; irand < n; irand++)		
	{
        //i = (int) floor(rand() / (RAND_MAX+0.0) * (nrow[0]+0.0));
		groupStart = i*p;
		if (useGroup[i] == 1)
		{
            //irand--;
			// Setting up null gradient calc to check if group is 0
			for (int k = 0; k < n; k++)
			{
				etaNull[k] = eta[k];
				if (k >= i) {
					for (int j = 0; j < p; j++)
					{
						etaNull[k] = etaNull[k] - X[k + n * j] * beta[j + groupStart];
					}
				}
			}
			// Calculating Null Gradient
			pCalc(nrow, etaNull, prob);
			logitGradCalc(nrow, prob, y, ldot);

			double *grad = NULL;
			grad = new double[p];

			for (int j = 0; j < p; j++)
			{
				grad[j] = 0;
				for (int k = i; k < n; k++)
				{
					grad[j] = grad[j] + X[k + n * j] * ldot[k];
				}
				if (grad[j] < lambda1[0] && grad[j] > -lambda1[0])
				{
					grad[j] = 0;
				}
				if (grad[j] > lambda1[0])
				{
					grad[j] = grad[j] - lambda1[0];
				}
				if (grad[j] < -lambda1[0])
				{
					grad[j] = grad[j] + lambda1[0];
				}
				if (pow(grad[j], 2) == pow(lambda1[0], 2))
				{
					grad[j] = 0;
				}
			}

			zeroCheck = 0;
			for (int j = 0; j < p; j++)
			{
				zeroCheck = zeroCheck + pow(grad[j], 2);
			}

			if (zeroCheck <= pow(lambda2[0], 2))   //Or not?
			{
				if (betaIsZero[i] == 0)
				{
					for (int j = 0; j < p; j++) {
						for (int k = i; k < n; k++) {
							eta[k] = eta[k] - X[k + n * j] * beta[j+groupStart];
						}
					}
				}
				betaIsZero[i] = 1;
				for (int j = 0; j < p; j++)
				{
					beta[j + groupStart] = 0;
				}
			}
			else
			{
				if (isActive[i] == 0)
				{
					groupChange = 1;
				}
				isActive[i] = 1;

				for (int k = 0; k < p; k++)
				{
					theta[k] = beta[k+groupStart];
				}

				betaIsZero[i] = 0;
				double *z = NULL;
				z = new double[p];
				double *U = NULL;
				U = new double[p];
				double *G = NULL;
				G = new double[p];
				//double *betaNew = NULL;
				//betaNew = new double[p];

				count = 0;
				check = 1000000;

				while (count <= innerIter[0] && check > thresh[0])
				{
					count++;

					pCalc(nrow, eta, prob);
					logitGradCalc(nrow, prob, y, ldot);

					for (int j = 0; j < p; j++)
					{
						grad[j] = 0;
						for (int k = i; k < n; k++)
						{
							grad[j] = grad[j] + X[k + n * j] * ldot[k];
						}

					}

					diff = -1;
					//	      t = 0.5;
					pCalc(nrow, eta, prob);
					Lold = logitNegLogLikelihoodCalc(nrow, prob, y);

					// Back-tracking

					while (diff < 0)
					{
						for (int j = 0; j < p; j++)
						{
							z[j] = beta[j + groupStart] - t * grad[j];
							if (z[j] < lambda1[0] * t && z[j] > -lambda1[0] * t)
							{
								z[j] = 0;
							}
							if (z[j] > lambda1[0] * t)
							{
								z[j] = z[j] - lambda1[0] * t;
							}
							if (z[j] < -lambda1[0] * t)
							{
								z[j] = z[j] + lambda1[0] * t;
							}
						}

						norm = 0;
						for (int j = 0; j < p; j++)
						{
							norm = norm + pow(z[j], 2);
						}
						norm = sqrt(norm);

						if (norm != 0){
							uOp = (1 - t * lambda2[0] / norm);   //Or not?
						}
						else{ uOp = 0; }

						if (uOp < 0)
						{
							uOp = 0;
						}

						for (int j = 0; j < p; j++)
						{
							U[j] = uOp*z[j];//l(U(beta_k,l,t)
							G[j] = 1 / t *(beta[j + groupStart] - U[j]);
						}

						// Setting up betaNew and etaNew in direction of Grad for descent step
						/*
						for (int j = 0; j < p; j++)
						{
							thetaNew[j] = beta[j+groupStart];
						}
						for (int j = 0; j < p]; j++)
						{
							thetaNew[j] = beta[j + groupStart] - t * G[j];
						}

						for (int j = rangeGroupInd[i] + groupLen[i]; j < ncol[0]; j++)
						{
							thetaNew[j] = beta[j];
						}*/
						for (int k = 0; k < nrow[0]; k++)
						{
							etaNew[k] = eta[k];
							if (k >= i) {
								for (int j = 0; j < p; j++)
								{
									etaNew[k] = etaNew[k] - t*G[j] * X[k + n * j];
								}
							}							
						}

						pCalc(nrow, etaNew, prob);
						Lnew = logitNegLogLikelihoodCalc(nrow, prob, y);

						sqNormG = 0;
						iProd = 0;

						for (int j = 0; j < p; j++)
						{
							sqNormG = sqNormG + pow(G[j], 2);
							iProd = iProd + grad[j] * G[j];
						}

						diff = Lold - Lnew - t * iProd + t / 2 * sqNormG;

						t = t * gamma[0];
					}
					t = t / gamma[0];

					check = 0;

					for (int j = 0; j < p; j++)
					{
						check = check + fabs(theta[j] - U[j]);
						for (int k = i; k < n; k++)
						{
							eta[k] = eta[k] - X[k + n * j] * beta[j + groupStart];
						}
						beta[j + groupStart] = U[j] + count / (count + 3) * (U[j] - theta[j]);
						theta[j] = U[j];

						for (int k = i; k < n; k++)
						{
							eta[k] = eta[k] + X[k + n * j] * beta[j + groupStart];
						}
					}
				}
				delete[] z;
				delete[] U;
				delete[] G;
				//delete[] betaNew;
			}
			delete[] grad;
		}
	}
	betaZeroSolve(nrow, betaZero, eta, prob, thresh, innerIter, y);

	delete[] etaNew;
	delete[] etaNull;
	delete[] theta;
	//delete[] thetaNew;
}
  


  int logitNest(double *X, int* y, int *nrow, int *ncol, double *lambda1, double *lambda2, double *beta, int *innerIter, int *outerIter, double *thresh, double *outerThresh, double *eta, double *gamma, int *betaIsZero, double* betaZero, double *step)
{
	int n = nrow[0];
	int p = ncol[0];
	int np = n*p;
	double oldBetaZero = betaZero[0];
	double* prob = NULL;
	prob = new double[n];
	//nullBeta never used, delete
	//double* nullBeta = NULL;
	//nullBeta = new double[ncol[0]];

	double *ldot = NULL;
	ldot = new double[n];
	int groupChange = 1;
	int* isActive = NULL;
	isActive = new int[n];
	int* useGroup = NULL;
	useGroup = new int[n];
	int* tempIsActive = NULL;
	tempIsActive = new int[n];

	for (int i = 0; i < n; i++)
	{
		isActive[i] = 0;
		useGroup[i] = 1;
	}

	// outer most loop creating response etc...
	int outermostCounter = 0;
	double outermostCheck = 1000000;
	double* outerOldBeta = NULL;
	outerOldBeta = new double[np];

	while (groupChange == 1)
	{
		groupChange = 0;
		logitSolver(X, y, nrow, ncol, beta, lambda1, lambda2, innerIter, thresh, ldot, gamma, eta, betaIsZero, groupChange, isActive, useGroup, prob, betaZero, step);
		while (outermostCounter < outerIter[0] && outermostCheck > outerThresh[0])
		{
			outermostCounter++;
			for (int i = 0; i < np; i++)
			{
				outerOldBeta[i] = beta[i];
			}
			oldBetaZero = betaZero[0];

			for (int i = 0; i < n; i++)
			{
				tempIsActive[i] = isActive[i];
			}

			logitSolver(X, y, nrow, ncol, beta, lambda1, lambda2, innerIter, thresh, ldot, gamma, eta, betaIsZero, groupChange, tempIsActive, isActive, prob, betaZero, step);
			/*for (int i = 0; i < np; i++) {
				cout << beta[i] << " ";
			}
			cout << endl;
			system("PAUSE");*/
			outermostCheck = 0;
			for (int i = 0; i < np; i++)
			{
				outermostCheck = outermostCheck + fabs(outerOldBeta[i] - beta[i]);
			}
			outermostCheck = outermostCheck + fabs(oldBetaZero - betaZero[0]);
		}
	}
	//delete [] nullBeta;
	delete[] outerOldBeta;
	delete[] ldot;
	delete[] isActive;
	delete[] useGroup;
	delete[] tempIsActive;
	delete[] prob;

	return 1;
}
}

//int logisticSGL(
/*
int main() 
{
  int* n = new int[1];
  int* p = new int[1];
  n[0] = 1000; p[0] = 20;
  double* X = new double[n[0]*p[0]];
  double* y = new double[n[0]];
  for(int k=0; k<n[0]; k++)
  {
    for(int j=0; j<p[0]; j++)
    {
      X[k+j*n[0]] = rand() / (RAND_MAX+0.0);
    }
    y[k] = X[k];
  }
  double* lambda1 = new double[1]; lambda1[0] = 0.001;
  double* lambda2 = new double[1]; lambda2[0] = 0.009;
  double* beta = new double[n[0]*p[0]];
  for(int i=0; i<n[0]*p[0]; i++)
  {
    beta[i] = 0;
  }
  double* thresh = new double[1]; thresh[0] = 0.01;
  double* outerThresh = new double[1]; outerThresh[0] = 0.01;
  double* eta = new double[n[0]];
  for(int k=0; k<n[0]; k++)
  {
    eta[k] = 0;
  }
  double* gamma = new double[1]; gamma[0] = 0.8;
  double* step = new double[1]; step[0] = 1;
  int* innerIter = new int[1]; innerIter[0] = 100;
  int* outerIter = new int[1]; outerIter[0] = 100;
  int* betaIsZero = new int[n[0]];
  for(int k=0; k<n[0]; k++)
  {
    betaIsZero[k] = 1;
  }
  int* reset = new int[1]; reset[0] = 20;
  solveSGL(X, y, n, p, lambda1, lambda2, beta, innerIter, outerIter, thresh, outerThresh, eta, gamma, betaIsZero, step, reset);
  for(int i=0; i<n[0]*p[0]; i++)
  {
    cout<<beta[i]<<" ";
    if(i>0 && i%p[0]==0)
    {
      cout<<endl;
    }
  }
  system("PAUSE");
  return 1;
}
*/
/*
int main() 
{
  int* n = new int[1];
  int* p = new int[1];
  int* s = new int[1];
  n[0] = 1000; p[0] = 200; s[0] = 10;
  double* X = new double[n[0]*p[0]];
  int* y = new int[n[0]];
  double tmp;
  for(int k=0; k<n[0]; k++)
  {
    for(int j=0; j<p[0]; j++)
    {
      X[k+j*n[0]] = rand() / (RAND_MAX+0.0);
    }
	tmp = 0;
	for (int l = 0; l < s[0]; l++) {
		tmp += X[k + l*n[0]];
	}
	if (k < n[0] / 2) {	
		y[k] = tmp>0 ? 1 : 0;
	}
	else {
		y[k] = tmp>0 ? 0 : 1;
	}
  }
  double* lambda1 = new double[1]; lambda1[0] = 0.001;
  double* lambda2 = new double[1]; lambda2[0] = 0.009;
  double* beta = new double[n[0]*p[0]];
  for(int i=0; i<n[0]*p[0]; i++)
  {
    beta[i] = 0;
  }
  double* thresh = new double[1]; thresh[0] = 0.01;
  double* outerThresh = new double[1]; outerThresh[0] = 0.01;
  double* eta = new double[n[0]];
  for(int k=0; k<n[0]; k++)
  {
    eta[k] = 0;
  }
  double* gamma = new double[1]; gamma[0] = 0.8;
  double* step = new double[1]; step[0] = 1;
  int* innerIter = new int[1]; innerIter[0] = 100;
  int* outerIter = new int[1]; outerIter[0] = 100;
  int* betaIsZero = new int[n[0]];
  for(int k=0; k<n[0]; k++)
  {
    betaIsZero[k] = 1;
  }
  int* reset = new int[1]; reset[0] = 20;
  double* betaZero = new double[1]; betaZero[0] = 0;
  logitNest(X, y, n, p, lambda1, lambda2, beta, innerIter, outerIter, thresh, outerThresh, eta, gamma, betaIsZero, betaZero, step);

  for(int i=0; i<n[0]; i++) {
	  tmp = 0;
	  for (int j = 0; j < p[0]; j++) {
		  tmp += pow(beta[i*p[0] + j], 2);
	  }
	  if (tmp > 0.00001) {
		  cout << i << " ";
	  }
  }
  cout << endl;
  for(int i=0; i<n[0]*p[0]; i++)
  {  
    if(i>0 && i%p[0]==0)
    {
      cout<<endl;
    }
    cout<<beta[i]<<" ";
  }
  system("PAUSE");
  return 1;
}*/
