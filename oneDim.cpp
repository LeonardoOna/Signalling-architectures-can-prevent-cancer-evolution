/*==================================================================================================================================
                                                     oneDim.cpp
====================================================================================================================================

C++-code accompanying:	
		 
		(Signalling architectures can prevent cancer evolution).

Written by:
        Leonardo Oña
       	Department of Ecology, School of Biology/Chemistry.
        University of Osnabrück
        Germany

Program version
		13/03/2019	:

Instructions for compiling and running the program
		
	Versions of this program were compiled and run on Windows and Mac, using Microsoft Visual C++
	2010 and XCode, respectively. The code is written in standard C++ and should be compatible with 
	other compilers. 

=================================================================================================================================*/


#include <cmath>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>


#include "random.h"

using namespace std;



//PARAMETERS

const double b=0.1;                         //Benefit
const double c=0.000;                        //Cost
//const int R=30;                           //Radius
const int Rmin=0;                           //Min Radius in a loop
const int Rmax=50;                          //Max Radius in a loop
const int Rstep=4;                          //R step in a loop
const int nR= ((Rmax-Rmin)/Rstep)+1;        //Number of data per Radius
const int L=10;                             //Range of local competition L
const int PopSize=100;                      //Population Size
const int first1=PopSize/2;                 //Position of the unique mutant signalling cell in the beginning of the simulation

const int nRuns=1000;                        //Number of runs


static unsigned int call_count = 0;         //Count number of calls to the Replacement Function (i.e generation step)




//FUNCTIONS

//The function "IntegerRandom" returns a random integer from min (inclusive) to max (exclusive): [min,max). When called with one argument gives [0,max)
int IntegerRandom(const int max, const int min=0)
{return static_cast<int>(rnd::uniform()*(max-min)+min);}


//The function "NeighbourWithinL" returns a value within the range [-L,L] (of an address) receiving the address of the focal element
/*int NeighbourWithinL(int addresoffocal)
{
    int kk=0;
    if (addresoffocal-L<0){
        kk= IntegerRandom(addresoffocal+L+1);
    }else if (addresoffocal+L>PopSize-1) {
        kk =IntegerRandom(PopSize,addresoffocal-L);
    } else{
        kk =IntegerRandom(addresoffocal+L+1,addresoffocal-L);
    }
    return kk;
}*/




//The function "NeighbourWithinL" returns a value within the range [-L,L] (of an address) receiving the address of the focal element
/*int NeighbourWithinL(const int addresoffocal)
{
    if (addresoffocal-L<0){
        return IntegerRandom(addresoffocal+L+1);
    }else if (addresoffocal+L>PopSize-1) {
        return IntegerRandom(PopSize,addresoffocal-L);
    } else{
        return IntegerRandom(addresoffocal+L+1,addresoffocal-L);
    }
}*/






//The function "NeighbourWithinL" returns a value within the range [-L,L] (of an address) receiving the address of the focal element
int NeighbourWithinL(int addresoffocal)
{
    int kk=0;
    if (addresoffocal-L<0){
        do{kk=IntegerRandom(addresoffocal+L+1);}while(kk==addresoffocal);
    }else if (addresoffocal+L>PopSize-1) {
        do{kk =IntegerRandom(PopSize,addresoffocal-L);} while(kk==addresoffocal);
    } else{
        do{kk =IntegerRandom(addresoffocal+L+1,addresoffocal-L);}while(kk==addresoffocal);
    }
    return kk;
}












int ChooseAnElement(const bool StateArray[PopSize], const int R)
{
    //**************** Initialise Fitness Array ********************************
    double fitnessArray[PopSize]={0};
    for (int i=0; i<PopSize; i++) {fitnessArray[i]=1;}
    //**************** This Part A: Cumulative Model ***************************
    
    /*
    for (int i=0; i<PopSize; i++) {
        for (int j=-R; j<(R+1); j++) {
            if ((i+j)>=0 && (i+j)<(PopSize)) {
                fitnessArray[i] += b * StateArray[i+j];
            }
        }
        
    }
     
     */
    
    
    //**************** This Part B: Non-Cumulative Model *************************
    
     for (int i=0; i<PopSize; i++) {
        for (int j=-R; j<(R+1); j++) {
            if ((i+j)>=0 && StateArray[i]==1 && (i+j) < PopSize) {
                fitnessArray[i+j] = 1+b;
            }
        }
     }
    
    
    
    
    //******************** Cost Valid for Any Model **************************
    for (int i=0; i<PopSize; i++) {fitnessArray[i]+= (-c) * StateArray[i];}
    //******************** Cumulative vector **************************
    double fitnessCum[PopSize]={0};
    for (int i=0; i<PopSize; i++){fitnessCum[i]=fitnessArray[i];}
    for (int i=1; i<PopSize; i++){fitnessCum[i] += fitnessCum[i-1];}
    //******************** Return an Element weighted sampled **************************
    return rnd::sample_1(fitnessCum, PopSize);
}









void Replacement(bool *pnArray, const int nLength, const int focal, const int replace)
{
    call_count++;   //Count number of calls (i.e generation step)
  
    // Vector "vect" will contain all the positions in "anArray" where there are "1".
    vector<int> vect;
    int eleg=0;
    
    for (int nCount=0; nCount < PopSize; nCount++)
        
        if (pnArray[nCount]==1) {
            vect.push_back(nCount);
        }
    
    
    
    
    if (pnArray[focal] != pnArray[replace]){
      
        if (pnArray[focal]==1) {
            
            
            
            /* Situation 1 0  */
            // There are 4 situations: 1) 1's in the middle, 2) 1' in corner r, 3) 1's in cornere left, 4) 1's all over.
            
            //This is what is going to go for the situation [1 0]
            if ((vect[0]!=0 && vect[(vect.size()-1)]!=(PopSize-1))) {                       //1) 1's in the middle
                (rnd::bernoulli() == 0) ? (eleg=(vect[0]-1)) : (eleg= (vect[(vect.size()-1)]+1)) ;
                pnArray[eleg]=1;
            }
            else if ((vect[0]==0 && vect[(vect.size()-1)]!=(PopSize-1))) {                  //2) 1' in corner r
                (eleg= (vect[(vect.size()-1)]+1));
                pnArray[eleg]=1;
            }
            else if ((vect[0]!=0 && vect[(vect.size()-1)]==(PopSize-1))) {                  //3) 1's in cornere left
                (eleg= (vect[0]-1));
                pnArray[eleg]=1;
            }
            else
                cout << "There is a fixation of the 1 type !!!!!!" << endl;                 //4) 1's all over
            
            
            
        } else {
            
            
            /* Situation 0 1 */
            // There are 4 situations: 1) 1's in the middle, 2) 1' in corner r, 3) 1's in cornere left, 4) 1's all over.
            
            //This is what is going to go for the situation [1 0]
            if ((vect[0]!=0 && vect[(vect.size()-1)]!=(PopSize-1))) {                       //1) 1's in the middle
                (rnd::bernoulli() == 0) ? (eleg=(vect[0])) : eleg= vect[(vect.size()-1)] ;
                //cout << "eleg posright vect 1: " << vect[(vect.size()-1)] << endl;
                //cout << "eleg posright vect 2: " << eleg << endl;
                pnArray[eleg]=0;
            }
            else if ((vect[0]==0 && vect[(vect.size()-1)]!=(PopSize-1))) {                  //2) 1' in corner r
                eleg= vect[(vect.size()-1)];
                pnArray[eleg]=0;
            }
            else if ((vect[0]!=0 && vect[(vect.size()-1)]==(PopSize-1))) {                  //3) 1's in cornere left
                eleg= vect[0];
                pnArray[eleg]=0;
            }
            else
                cout << "There is a fixation of the 1 type !!!!!!" << endl;                 //4) 1's all over
            
        }
        
    }

    vect.erase (vect.begin(), vect.begin()+vect.size());                    // Erasing vector
    
}



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





//MAIN



int main()

{
    ofstream outputfile ("output.csv");
    
    double fixProbR[nR]={0.0};
    double avFixTimeR[nR]={0.0};

    
    for (int R=Rmin; R <= Rmax ; R=R+Rstep) {

        int whoWon[nRuns]={0};
        int fixTime[nRuns]={0};

        
        
        for (int i=0; i<nRuns; i++){
    
            bool ExampleArray[PopSize]={0};
            ExampleArray[first1]=1;
    
    
            do {
                int temp=ChooseAnElement(ExampleArray,R);
                Replacement(ExampleArray, PopSize, temp, NeighbourWithinL(temp));
            } while(Sum(ExampleArray, PopSize)!=0 && Sum(ExampleArray, PopSize)!=PopSize);
    
        
            whoWon[i]=ExampleArray[1];
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










