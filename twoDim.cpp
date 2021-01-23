/*==================================================================================================================================
                                                     twoDim.cpp
====================================================================================================================================

 C++-code accompanying:
 
 (Signalling architectures can prevent cancer evolution).
 
 Written by:
 Leonardo Oña
 Department of Ecology, School of Biology/Chemistry.
 University of Osnabrück
 Germany
 
 Program version
 13/03/2019    :

Instructions for compiling and running the program
		
	Versions of this program were compiled and run on Windows and Mac, using Microsoft Visual C++
	2010 and XCode, respectively. The code is written in standard C++ and should be compatible with 
	other compilers. 

=================================================================================================================================*/
//
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>


#include "random.h"



using namespace std;





//PARAMETERS

const double b=0.01;                         //Benefit
const double c=0.00;                        //Cost
//const int R=30;                           //Radius
const int Rmin=1;                           //Min Radius in a loop
const int Rmax=5;                          //Max Radius in a loop
const int Rstep=1;                          //R step in a loop
const int nR= ((Rmax-Rmin)/Rstep)+1;        //Number of data per Radius
const int L=2;                              //Range of local competition L
const int iSize=10;
const int jSize=10;

const int PopSize = iSize * jSize;
//const int first1=PopSize/2;               //Position of the unique mutant signalling cell in the begining of the simulation

const int nRuns=1000;                        //Number of runs


static unsigned int call_count = 0;         //Count number of calls to the Replacement Function (i.e generation step)







//FUNCTIONS

//The function "IntegerRandom" returns a random integer from min (inclusive) to max (exclusive): [min,max). When called with one argument gives [0,max)
int IntegerRandom(const int max, const int min=0)
{return static_cast<int>(rnd::uniform()*(max-min)+min);}





int Sum(const bool *pnArray, const int nLength)
{
    int val =0;
    for (int i=0; i<nLength; i++)
        val+=pnArray[i];
    return val;
}





int Sum(const int *pnArray, const int nLength)
{
    int val =0;
    for (int i=0; i<nLength; i++)
        val+=pnArray[i];
    return val;
}



int SumArr2(const bool pnArray[][jSize], const int nLength, const int mLength)
{
    int val =0;
    for (int i=0; i<nLength; i++){
        for (int j=0; j<mLength; j++) {
            val+=pnArray[i][j];
        }
    }
    return val;
}




int MatToVet(const int ik, const int jk)
{
    int ret=0;
    int counn=0;
    for (int i=0; i<iSize; i++) {
        for (int j=0; j<jSize; j++) {
            if (i==ik && j==jk) {
                ret= counn;
            }
            counn++;
        }
    }
    return ret;
}



int VetToMat(const int valueFromVector, const int ZeroOrOne)
{
    int MapVM[PopSize][3]={0};
    
    int count2=-1;
    int count3=0;
    
    for (int i=0; i<PopSize; i++) {
        count2++;
        for (int j=0; j<3; j++) {
            if (j==0) {
                MapVM[i][j]=count2;
            } else
                if (j==1){
                    MapVM[i][j]= count3%iSize;
                    if (j==1 && count2%jSize == (jSize-1)) {
                        count3++;
                    }
                }
                else
                    if (j==2) {
                        MapVM[i][j]= count2%jSize;
                    }
        }
    }
    return MapVM[valueFromVector][ZeroOrOne+1];
}










int ChooseAnElement(const bool StateArray[iSize][jSize], const int R)
{
    //**************** Initialise Fitness Array ********************************
    double MatrixFitness[iSize][jSize]={0};
    for (int i=0; i<iSize; i++){
        for (int j=0; j<jSize; j++) {
            MatrixFitness[i][j]=1;
        }
    }
    //**************** This Part A: Cumulative Model ****************************
    /*
    for (int i=0; i<iSize; i++) {
        for (int j=0; j<jSize; j++) {
            if (StateArray[i][j]==1) {
                for (int k = -R; k <= R ; k++) {
                    for (int l = -R ; l <= R ; l++) {
                        if ((i+k)>=0 && (i+k)<iSize && (j+l)>=0 && (j+l)<jSize) {
                            MatrixFitness[i+k][j+l] += b * StateArray[i][j];
                        }
                    }
                }
            }
        }
    }
     
     */
    
    //**************** This Part B: Non-Cumulative Model *************************
    
    for (int i=0; i<iSize; i++) {
        for (int j=0; j<jSize; j++) {
            if (StateArray[i][j]==1) {
                for (int k = -R; k <= R ; k++) {
                    for (int l = -R ; l <= R ; l++) {
                        if ((i+k)>=0 && (i+k)<iSize && (j+l)>=0 && (j+l)<jSize) {
                        MatrixFitness[i+k][j+l]=b;
                        }
                    }
                }
            }
        }
    }
    
    //******************** Cost Valid for Any Model **************************
    for (int i=0; i<iSize; i++){
        for (int j=0; j<jSize; j++) {
            if (StateArray[i][j]==1) {
                MatrixFitness[i][j]+=(-c)*StateArray[i][j];
            }
        }
    }
    //**************** Transforming Matrix Fitness to Array Fitness*************
    double FitnessArray[PopSize]={0};
    int counter=0;
    for (int i=0; i<iSize; i++) {
        for (int j=0; j<jSize; j++) {
            FitnessArray [counter]= MatrixFitness[i][j];
            counter++;
        }
    }
    //******************** Cumulative vector **************************
    double FitnessCum[PopSize]={0};
    for (int i=0; i<PopSize; i++){FitnessCum[i]=FitnessArray[i];}
    for (int i=1; i<PopSize; i++){FitnessCum[i] += FitnessCum[i-1];}
    //******************** Return an Element weighted sampled **************************
    return rnd::sample_1(FitnessCum, PopSize);
}




int NeighbourWithinL(const int addresoffocal)
{
    int a[2]={0,0};
    do{
        a[0] = VetToMat(addresoffocal,0) - L + IntegerRandom(2*L+1);
        a[1] = VetToMat(addresoffocal,1) - L + IntegerRandom(2*L+1);
    } while ((a[0]==VetToMat(addresoffocal,0) && a[1]==VetToMat(addresoffocal,1)) || a[0]<0 || a[1]<0 || a[0]>(iSize-1) || a[1]>(jSize-1));
    return MatToVet(a[0], a[1]);
}


//void do_something(int (&array)[board_width][board_height]);

void Replacement(bool (&MatrixSt)[iSize][jSize], const int focal, const int replace)
{
    call_count++;

    
    if (MatrixSt[VetToMat(focal, 0)][VetToMat(focal, 1)]!=MatrixSt[VetToMat(replace, 0)][VetToMat(replace, 1)]) {
        if (MatrixSt[VetToMat(focal, 0)][VetToMat(focal, 1)]==1) {
            vector<int> VectorZeroAddress;
            vector<int> VectorNOnes;
            
            for (int i=0; i<iSize; i++) {
                for (int j=0; j<jSize; j++) {
                    if (MatrixSt[i][j]==0) {
                        int counter2=0;
                        for (int k=-1+i; k<2+i; k++) {
                            for (int l=-1+j; l<2+j; l++) {
                                if (k>=0 && k < iSize && l>=0 && l<jSize && MatrixSt[k][l]==1)
                                    counter2++;
                            }
                        }
                        if (counter2 !=0) {
                            VectorZeroAddress.push_back(MatToVet(i,j));
                            VectorNOnes.push_back(counter2);
                        }
                    }
                }
            }
            //############################
            
            vector<int> VectorAddresMaximum;
            
            for (int i=0; i<(VectorZeroAddress.size()); i++) {
                if (VectorNOnes[i]==*max_element(begin(VectorNOnes),end(VectorNOnes))) {
                    VectorAddresMaximum.push_back(VectorZeroAddress[i]);
                }
            }
            
            
            int choo=VectorAddresMaximum[IntegerRandom(static_cast<int>(VectorAddresMaximum.size()))];
            
            MatrixSt[VetToMat(choo, 0)][VetToMat(choo, 1)]= 1 ;
            
            
            VectorZeroAddress.erase (VectorZeroAddress.begin(), VectorZeroAddress.begin()+VectorZeroAddress.size());
            VectorNOnes.erase (VectorNOnes.begin(), VectorNOnes.begin()+VectorNOnes.size());
            VectorAddresMaximum.erase (VectorAddresMaximum.begin(), VectorAddresMaximum.begin()+VectorAddresMaximum.size());
        }
        else{
            vector<int> VectorZeroAddress;
            vector<int> VectorNOnes;
            
            for (int i=0; i<iSize; i++) {
                for (int j=0; j<jSize; j++) {
                    if (MatrixSt[i][j]==1) {
                        int counter2=0;
                        for (int k=-1+i; k<2+i; k++) {
                            for (int l=-1+j; l<2+j; l++) {
                                if (k>=0 && k < iSize && l>=0 && l<jSize && MatrixSt[k][l]==0)
                                    counter2++;
                            }
                        }
                        if (counter2 !=0) {
                            VectorZeroAddress.push_back(MatToVet(i,j));
                            VectorNOnes.push_back(counter2);
                        }
                    }
                }
            }
            //############################
            
            vector<int> VectorAddresMaximum;
            
            for (int i=0; i<(VectorZeroAddress.size()); i++) {
                if (VectorNOnes[i]==*max_element(begin(VectorNOnes),end(VectorNOnes))) {
                    VectorAddresMaximum.push_back(VectorZeroAddress[i]);
                }
            }
            
            
            int choo=VectorAddresMaximum[IntegerRandom(static_cast<int>(VectorAddresMaximum.size()))];
            
            MatrixSt[VetToMat(choo, 0)][VetToMat(choo, 1)]= 0 ;
            
            
            VectorZeroAddress.erase (VectorZeroAddress.begin(), VectorZeroAddress.begin()+VectorZeroAddress.size());
            VectorNOnes.erase (VectorNOnes.begin(), VectorNOnes.begin()+VectorNOnes.size());
            VectorAddresMaximum.erase (VectorAddresMaximum.begin(), VectorAddresMaximum.begin()+VectorAddresMaximum.size());
            

            
        }
    
    }
    
}







int main()
{
    
    ofstream outputfile ("output.csv");
    
    

    
    double fixProbR[nR]={0.0};
    double avFixTimeR[nR]={0.0};
    
    
    for (int R=Rmin; R <= Rmax ; R=R+Rstep) {
    
    
    
    
    int whoWon[nRuns]={0};
    int fixTime[nRuns]={0};
    
    
    for (int i=0; i<nRuns; i++){
  
        bool MatrixState [iSize][jSize]={0};
        MatrixState [iSize/2][jSize/2]=1;

    do {
        int temp=ChooseAnElement(MatrixState,R);
        Replacement(MatrixState, temp,NeighbourWithinL(temp));
    } while(SumArr2(MatrixState, iSize,jSize)!=0 && SumArr2(MatrixState, iSize,jSize)!=PopSize);
    
    
        whoWon[i]=MatrixState[1][1];
        fixTime[i]=call_count;
        
        
        call_count=0;
        
    }
    
    
        
        
        
        
        
        //****The following vector save cases where 1 reach fixation*****
        
        int OnlyCasesOneFixTime[nRuns]={0};
        
        int Counter1=0;
        
        for (int i=0; i<nRuns; i++) {
            if(whoWon[i]==1){
                Counter1++;
                OnlyCasesOneFixTime[i]=fixTime[i];
            }
        } // End of the For loop for Runs values
        //****************************************************************
        
        fixProbR[(R-Rmin)/Rstep]= double(Sum(whoWon,nRuns))/nRuns;
        avFixTimeR[(R-Rmin)/Rstep]= (double(Sum(OnlyCasesOneFixTime,nRuns))/Counter1)/PopSize;
        
    }  // End of the For loop for R values

        
        
        
    
        
    cout << endl;
    cout << endl;
    
    
    cout << "Radius:"<< endl;
    for (int i=Rmin; i<=Rmax; i=i+Rstep) {
        cout << i << " ";
    }
    
    cout << endl;
    cout << endl;
    
    
    
    cout << "Fix Prob dep R:"<< endl;
    for (int i=0; i<=(Rmax-Rmin)/Rstep; i++) {
        cout << fixProbR[i] << " ";
    }
    
    cout << endl;
    cout << endl;
    
    cout << "Av Fix Time dep R:"<< endl;
    for (int i=0; i<=(Rmax-Rmin)/Rstep; i++) {
        cout << avFixTimeR[i] << " ";
    }
    
    
    
    cout << endl;
    cout << endl;
    cout << endl;
    

    
    for (int i=0; i<=(Rmax-Rmin)/Rstep; i++)  {
        outputfile << i << "," << fixProbR[i] << "," << avFixTimeR[i] << endl;
    }
    
    
 
    
    
    outputfile.close();
    
    
    return 0;
}








