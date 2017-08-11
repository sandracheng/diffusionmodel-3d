/*SNO+ 3D diffusion model for water in acrylic
 
  Author:  Sandra Cheng. 10184319.  <14wlsc@queensu.ca>
  REVISION HISTORY: see GitHub. Project name is: diffusionmodel-3d
************
  Based on Thomas' initial 3D Diffusion Model code. This is for water in
  an acrylic piece, like a dogbone.       
  Acrylic dogbone, with the bottom front left corner at (0,0,0). 
    
    /|\ z-axis
     |
     | 
    ______        
   |\_____\     
   | |    |
   | |    | 
   | |    |
   | |    |        __  y-axis
   | |    |         /|
   \  \  /         /
    | |  |        /
    | |  |       /
    | |  |      /
    | |  |     /
    | |  |    /
    | |  |   /
    / /  \  /
   | |    |/
   | |    |
   | |   /|
   | |  / |
   | | /  |
    \|/___|  ______________\ x-axis
                           /         
  The code is based off 1000 cubic elements to make a 10x10x10 acrylic cube.
  So number of divisions=divs=10. 
  The outer layer has a random concentration, 1st inner layer will be the 
  water layer with maximum concentration, and the 2nd inner layer and
  beyond will be the acrylic.


  [divs=10 is not hard-coded in, but divs must be at least 5 to make sense. So for 
  divs=5 we have an outer layer of 61 elements. 2nd layer has 26, 3rd has 1. 
  The water layer is the 2nd layer, the 3rd layer is the acrylic.             

  For divs=n, outer layer has n^3-(n-1)^3 elements. 2nd layer is a cube
  with divs=n-2;1st layer is a cube of divs=n-4. 
  e.g. A 5x5x5 cube can enclose a 3x3x3 which can enclose a 1x1x1. 
  A 6x6x6 cube can enclose a 4x4x4 which can enclose a 2x2x2.]
 
                                                            
  The 10x10x10 acrylic cube may then be normalized to any size acrylic dogbone.
  Origin will still be at bottom left front corner of the cube. 
     _______________
    /990        999/|
   / |            / |
  /__|___________/  |
  |90|         99|  |
  |  |           |  |
  |  |           |  |
  |  |           |  |
  |  900---------|--/909
  | /            | /
  |/0___________9|/
  The corners will be as above.
  **if looking down on the top face as above:   990.....999
                                                 :        :
                                                90.......99
       
  **if kept flat and rotate to the left,the corners will be:   990.......90
                                                                :        :
                                                               900.......0
  **if kept flat and rotate to the right,the corners will be:  99......999
                                                                :       :
                                                                9......909
                                  
  **if looking through the front face, the back corners will be:  990....999
                                                                   :      :
                                                                  900....909
  We get 2 plots from the code, concentration for a single element vs. time
  and total mass of H2O in acrylic vs. time.   
*/

#include <stdio.h> 
#include <math.h>      
#include <stdlib.h> 

#include "TClass.h"
#include "TRandom3.h" //new PRNG based on Mersenne Twister
#include "TUUID.h"

#include "TCanvas.h"
#include "TH1.h"

void diffusionmodel3dnewarray() {                                                                    
 const int divs=5;//n.o. elements aka n.o. divisions in the big acrylic "cube"
 const int totTime=50;
 //const int timestep=1; 
 int timePassed;
          
 // measurements are taken based on whatever acrylic sample you're doing
 // these are based off Tom's 3rd dogbone and are in metres
 double xTotLen =0.01888;
 double yTotLen =0.00617;
 double zTotLen =0.20320;

 //x,y,and z lengths of each cube
 double xLen= xTotLen/divs; 
 double yLen= yTotLen/divs;
 double zLen= zTotLen/divs;

 //volume of each element 
 double elemVol= xLen*yLen*zLen;
 double maxConc= 0.0262599/(xTotLen*yTotLen*zTotLen); //max conc of water in moles/unit volume. 
 //the decimal # is given by (2% of acrylic in mass) divided by water's molecular weight
 
 double elemConc[divs*divs*divs]; //array for current conc in each element
                                  //each element is defined by its position
 double elemConcMaster[divs*divs*divs*totTime];  //updated array that
  //has all the conc for each elem at all times. 

 double DcoeffWaterRT= pow(4.22,-12); //diffusion constant for water in acrylic  
 //number is Tom's RT 'soak' D value, 3.65x10^-7  m^2/day converted into m^2/s
    
 int i,j,k; //counter for x,y and z positions respectively for the entire cube
 int x,y,z; //counter for the x,y and z  positions for the 2nd layer
 int elemCount; //counter for the total number of elements

 const int scndLayerCount= divs-2; //counter for the second outside layer
 int n,b,v,p,c,s,f,g; //just counters.
 int moleculeCount;
 const double holderConc=0.003; //this is just a random conc for the outer first layer

 //imagine having 'leader' molecules and 'follower' molecules
 //wherever the leader molecules go, their group follows...
 //fakeAvoNum is the number of leader molecules, and 10^23/fakeAvoNum = follower molecules 
 //need to implement this cause can't loop for 10^23 times.
 
 const double fakeAvoNum=pow(10,5); //representation of Avogadro's number, number of leaders
 int moleculeGrpInCube; //number of molecule groups led by 'leaders' followed by 'followers' in elem0
 const double truAvoNum=pow(10,23); //the real magnitude of Avogadro's number.
 const double followerMolecules= truAvoNum/fakeAvoNum; //number of followers

 int chosenElem; //the elem that the lead (and so the follower) molecules 'choose' to go into

 int frontFaceCounter=0;
 int moleculeFFCounter=0;

 TCanvas *c1 = new TCanvas;           
 Int_t nBins = 1000;
 Double_t lowlim = -1500;
 Double_t uplim = 1500;

 TH1* h1 = new TH1D("Legend","Title", nBins,lowlim,uplim);
 h1->SetFillColor(kOrange);

// double DcoeffAir; //diffusion constant for saturated acrylic in dry are
// double airTime; //time it has been out in air; not sure if need
   
 //need to reset the concentrations for the 1st and 2nd layers always for
 //both the current and update array. 1st layer=0, 2nd layer=maxConc.
 //all of inside starts at conc=0   
 for (timePassed=1; timePassed<=totTime; timePassed++) { 
  for (elemCount=0; elemCount<divs*divs*divs; elemCount++) { 
   //1st layer aka the cube faces
   //could have grouped front and back together but
   //need to separate front and back in order to do specific kind of testing
   if (elemCount <= divs*divs-1) {         //front face
    elemConc[elemCount]= holderConc;   
    elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs]=holderConc;
   } 

   if (elemCount >= divs*divs*(divs-1)){    //back face
    elemConc[elemCount]= holderConc;   
    elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs]=holderConc;
   } 

   for (j=0; j<divs; j++) {    
    for (k=0; k<divs*divs; k++) {
     if (k<divs) {
      elemConc[j*divs*divs+k*divs]= holderConc;  //for left face
      elemConcMaster[j*divs*divs+k*divs+(timePassed-1)*divs*divs*divs]= holderConc;

      elemConc[j*divs*divs+k*divs+divs-1]= holderConc;  //for right face
      elemConcMaster[j*divs*divs+k*divs+divs-1+(timePassed-1)*divs*divs*divs]=holderConc;

      elemConc[j*divs*divs+k]= holderConc; //for bottom face
      elemConcMaster[j*divs*divs+k+(timePassed-1)*divs*divs*divs]=holderConc;
     }
     else if (k>=divs*(divs-1)) {
      elemConc[j*divs*divs+k]= holderConc; //for top face
      elemConcMaster[j*divs*divs+k+(timePassed-1)*divs*divs*divs]=holderConc;
     }
    }
   }
  //for the 2nd layer, the water layer           
   for (z=0; z<scndLayerCount;z++) {
    for (y=0; y<scndLayerCount;y++) {
    
     if (0<z && z<scndLayerCount-3 && 0<y && y<scndLayerCount-3) { 
      elemConc[divs*(divs+z+1)+y+1]=0;   //2nd layer front face                              
      elemConc[divs*(divs*scndLayerCount+z+1)+y+1]= 0;   //2nd layer other faces                              
     }

     if (0<z && z<scndLayerCount-3 && 0<y && y<scndLayerCount-3) { 
      elemConc[divs*(divs+z+1)+y+1]=maxConc;   //2nd layer front face                              
      elemConcMaster[divs*(divs+z+1)+y+1+(timePassed-1)*divs*divs*divs]=maxConc;

      elemConc[divs*(divs*scndLayerCount+z+1)+y+1]= maxConc;   //2nd layer back face                              
      elemConcMaster[divs*(divs*scndLayerCount+z+1)+y+1+(timePassed-1)*divs*divs*divs]= maxConc;
     }

     if (z==0 && y==0) {
      for (n=0; n<scndLayerCount; n++) {
       for (b=0; b<scndLayerCount; b++) {
        int n_r = n*divs*divs;

        elemConc[divs*(divs+1)+1+n_r+divs*b]= maxConc; //2nd layer left face
        elemConcMaster[divs*(divs+1)+1+n_r+divs*b+(timePassed-1)*divs*divs*divs]=maxConc;  
       
       }
      }
      for (v=0;v<scndLayerCount;v++){
       for (p=0; p<scndLayerCount;p++){
        int p_r=p*divs*divs;

        elemConc[divs*(divs+1)+1+v+p_r]=maxConc; //2nd layer bottom face
        elemConcMaster[divs*(divs+1)+1+v+p_r+(timePassed-1)*divs*divs*divs]= maxConc; 
       }
      }
     }//end bracket for z==0 and y==0

     if (z==0 && y==(scndLayerCount-1)) {
      for (c=0; c<scndLayerCount; c++) {
       for (s=0; s<scndLayerCount; s++) {
        int c_r = c*divs*divs+divs*s;

        elemConc[divs*(divs+1)+y+1+c_r]= maxConc; //2nd layer right face
        elemConcMaster[divs*(divs+1)+y+1+c_r+(timePassed-1)*divs*divs*divs]=maxConc; 
       }    
      }
     }

     if (z==(scndLayerCount-1) && y==0) {
      for (f=0;f<z;f++) {
       for (g=0; g<scndLayerCount;g++) {
        int g_r=  g*divs*divs;

        elemConc[divs*(divs+z+1)+1+f+g_r] = maxConc; //2nd layer top face
        elemConcMaster[divs*(divs+z+1)+1+f+g_r+(timePassed-1)*divs*divs*divs]=maxConc; 
       }
      }
     } //end bracket for z==scndLayerCount-1 and y==0
    }
   }        
                                                     
//cout << "time: "<< timePassed << " elem: " << elemCount << " elem conc: " << elemConc[elemCount] << endl;
//up to here the initialization values are corrrrrrrrect yay
  }  //end bracket for first elemCount                                                         

  for (elemCount=0; elemCount<divs*divs*divs; elemCount++) {                  
//  printf("elem:%d; elem conc: %f\n", elemCount, elemConc[elemCount]); // works!  

   moleculeGrpInCube = (int)(elemConc[elemCount]*elemVol*fakeAvoNum);
   chosenElem=0; //reset to 0 every element count loop. 
 
//make random walk occur for all layers but the first; just reset the 2nd layer
//every timestep
   if (elemConcMaster[elemCount]!=holderConc) {
    for (moleculeCount=1; moleculeCount<=moleculeGrpInCube; moleculeCount++){                     
/* need to know position of the 26 surrounding cube elements. will name it 
   by the number scheme as above, with elements labelled 0 to 26.
   The central cube number is 13                   
   
   slices of a front-facing cube:
   
   1st             2nd             3rd
   -----------     ------------    -----------
   |6   7   8|     |15  16  17|    |24  25  26|
   |3   4   5|     |12  13  14|    |21  22  23|
   |0   1   2|     |9   10  11|    |18  19  20|
   -----------     ------------    ------------

               
   cube13= elemConc[elemCount];              main cube             
   cube12= elemConc[elemCount-1];            left        
   cube14= elemConc[elemCount+1];            right
   cube16= elemConc[elemCount+divs];         top    
   cube10= elemConc[elemCount-divs];         bot
   cube4= elemConc[elemCount-divs*divs];     front
   cube22= elemConc[elemCount+divs*divs];    back
   cube0= elemConc[elemCount-divs*divs-1-divs];
   cube1= elemConc[elemCount-divs*divs-divs];         
   cube2= elemConc[elemCount-divs*divs+1-divs];
   cube3= elemConc[elemCount-divs*divs-1];
   cube5= elemConc[elemCount-divs*divs+1];
   cube6= elemConc[elemCount-divs*divs-1+divs];  
   cube7= elemConc[elemCount-divs*divs+divs];           
   cube8= elemConc[elemCount-divs*divs+1+divs];
   cube9= elemConc[elemCount-1-divs];   
   cube11= elemConc[elemCount+1-divs];
   cube15= elemConc[elemCount+divs-1];
   cube17= elemConc[elemCount+divs+1];
   cube18= elemConc[elemCount-1-divs+divs*divs];
   cube19= elemConc[elemCount-divs+divs*divs];
   cube20= elemConc[elemCount-divs+divs*divs+1];
   cube21= elemConc[elemCount-1+divs*divs];
   cube23= elemConc[elemCount+divs*divs+1];      
   cube24= elemConc[elemCount+divs-1+divs*divs];
   cube25= elemConc[elemCount+divs+divs*divs];
   cube26= elemConc[elemCount+divs+1+divs*divs];                         */ 

 // cout << elemCount << " " << elemConc[elemCount] << endl;
    unsigned int StartSeed;
    gRandom->SetSeed(StartSeed);
        
    Int_t r = (gRandom->Rndm())*27; 
//    cout << timePassed << " " << elemCount << " "  << moleculeCount << " " << r << endl;
//    cout << r << " moleculeCount:" << moleculeCount << "/" << moleculeGrpInCube << endl; 
   
    //above is a PRNG to get an int from [0, 27-1]
    //for each molecule to move into one of the adjacent 26 cubes or stay.
    //all if statements match the random number to the cube it moves to

     if (r==0 && (elemCount-divs*divs-1-divs)>=0 && (elemCount-divs*divs-1-divs)<divs*divs*divs) {   
     //if the elem exists, aka is between index 0 and divs*divs*divs-1 then run
                                 
     //Ed says to scale outside the if statements. So let's just choose the element only.
     //at end of if statements, change concs of both elems and update master array for both... 
      chosenElem= elemCount-divs*divs-1-divs;
     }                                                   
     else if (r==1 && (elemCount-divs*divs-divs)>=0 && (elemCount-divs*divs-divs)<divs*divs*divs) {                               
      chosenElem= elemCount-divs*divs-1-divs;
     }                 
     else if (r==2 && (elemCount-divs*divs-divs+1)>=0 && (elemCount-divs*divs-divs+1)<divs*divs*divs) {    
      chosenElem= elemCount-divs*divs-divs+1;
     }    
     else if (r==3 && (elemCount-divs*divs-1)>=0 && (elemCount-divs*divs-1)<divs*divs*divs)  {                                    
      chosenElem= elemCount-divs*divs-1;
     }
     else if (r==4 && (elemCount-divs*divs)>=0 && (elemCount-divs*divs)<divs*divs*divs) {                                    
      chosenElem= elemCount-divs*divs;
     }
     else if (r==5 && (elemCount-divs*divs+1)>=0 && (elemCount-divs*divs+1)<divs*divs*divs) {                                    
      chosenElem= elemCount-divs*divs+1;
     }
     else if (r==6 && (elemCount-divs*divs-1+divs)>=0 && (elemCount-divs*divs-1+divs)<divs*divs*divs) {                          
      chosenElem= elemCount-divs*divs-1+divs;
     }
     else if (r==7 && (elemCount-divs*divs+divs)>=0 && (elemCount-divs*divs+divs)<divs*divs*divs) {                                     
      chosenElem= elemCount-divs*divs+divs;
     }
     else if (r==8 && (elemCount-divs*divs+1+divs)>=0 && (elemCount-divs*divs+1+divs)<divs*divs*divs) {                              
      chosenElem= elemCount-divs*divs+1+divs;
     }
     else if (r==9 && (elemCount-1-divs)>=0 && (elemCount-1-divs)<divs*divs*divs) {                                    
      chosenElem= elemCount-1-divs;
     }                                            
     else if (r==10 && (elemCount-divs)>=0 && (elemCount-divs)<divs*divs*divs) {                                    
      chosenElem= elemCount-divs;
     }
     else if (r==11 && (elemCount+1-divs)>=0 && (elemCount+1-divs)<divs*divs*divs) {                                    
      chosenElem= elemCount+1-divs;
     }
     else if (r==12 && (elemCount-1)>=0 && (elemCount-1)<divs*divs*divs) {                                    
      chosenElem= elemCount-1;
     }                 
     else if (r==14 && (elemCount+1)>=0 && (elemCount+1)<divs*divs*divs) {                                    
      chosenElem= elemCount+1;
     }
     else if (r==15 && (elemCount+divs-1)>=0 && (elemCount+divs-1)<divs*divs*divs) {                                    
      chosenElem= elemCount+divs-1;
     }
     else if (r==16 && (elemCount+divs)>=0 && (elemCount+divs)<divs*divs*divs) {                                    
      chosenElem= elemCount+divs;
     }
     else if (r==17 && (elemCount+divs+1)>=0 && (elemCount+divs+1)<divs*divs*divs) {                                    
      chosenElem= elemCount+divs+1;
     }
     else if (r==18 && (elemCount-1-divs+divs*divs)>=0 && (elemCount-1-divs+divs*divs)<divs*divs*divs) {                                    
      chosenElem= elemCount-1-divs+divs*divs;
     }
     else if (r==19 && (elemCount+divs*divs-divs)>=0 && (elemCount+divs*divs-divs)<divs*divs*divs) {                                    
      chosenElem= elemCount+divs*divs-divs;
     }
     else if (r==20 && (elemCount+1+divs*divs-divs)>=0 && (elemCount+1+divs*divs-divs)<divs*divs*divs) {                                   
      chosenElem= elemCount+1+divs*divs-divs; 
     }
     else if (r==21 && (elemCount-1+divs*divs)>=0 && (elemCount-1+divs*divs)<divs*divs*divs) {                                    
      chosenElem= elemCount-1+divs*divs;  
     }
     else if (r==22 && (elemCount+divs*divs>=0) && (elemCount+divs*divs)<divs*divs*divs) {                                    
      chosenElem= elemCount+divs*divs;
     }                               
     else if (r==23 && (elemCount+divs*divs+1)>=0 && (elemCount+divs*divs+1)<divs*divs*divs) {                                   
      chosenElem= elemCount+divs*divs+1;
     }
     else if (r==24 && (elemCount+divs+divs*divs-1)>=0 && (elemCount+divs+divs*divs-1)<divs*divs*divs) {                                    
      chosenElem= elemCount-divs+divs*divs+1;
     }
     else if (r==25 && (elemCount+divs+divs*divs)>=0 && (elemCount+divs+divs*divs)<divs*divs*divs) {      
      chosenElem= elemCount+divs+divs*divs; 
     }
     else if (r==26 && (elemCount+divs+1+divs*divs)>=0 && (elemCount+divs+1+divs*divs)<divs*divs*divs) {                              
      chosenElem= elemCount+divs+1+divs*divs; 
     }
     else if (r==13){   //for r=13, the molecule stays put                                  
      chosenElem= elemCount; 
     }
                        
   //for the chosen elem, update conc for the current and master array.   
   elemConc[chosenElem]= (elemConc[chosenElem]*elemVol*truAvoNum+followerMolecules)/(elemVol*truAvoNum);
   elemConcMaster[chosenElem+(timePassed-1)*divs*divs*divs]= elemConc[chosenElem];  
 
//   cout << timePassed << " " << elemCount << " chosen elem conc: " << elemConc[chosenElem] << endl;   
  
   //adjust elem13 accordingly 
   elemConc[elemCount]= (elemConc[elemCount]*elemVol*truAvoNum-followerMolecules)/(elemVol*truAvoNum);
   elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs]= elemConc[elemCount];  

//   cout << timePassed << " " << elemCount << " elem 0 conc: " << elemConc[elemCount] << endl;
   //the formula above should work for r==13 as well; we get a net 0.
   
   } //end bracket for the molecule count per elem 

// cout << "time: " << timePassed << " elem:" << elemCount << " conc:" << elemConc[elemCount] << " no. moleculeGrps:" << moleculeGrpInCube << endl; 

  //want to calculate the total amount of water passing into the LAB layer
  //so we sum up the amount of water in the inner cube(s)
  //just for a quick calculation, I'll do 5 layers so the inner cube
  //is just the one cube, cube number 62. 
  //I'll also only set the front face to be maxConc, else is 0.
   if (elemCount==62) { 
    cout << "time passed:" << timePassed << " element conc:" << elemConc[elemCount] <<endl;  
   }     
     

/*  //can also do it for the entire front face of the inner cube.
  if (elemCount<divs*divs+divs*divs && elemConcMaster[elemCount]==0 && timePassed==totTime) {
   frontFaceCounter+=1;
   moleculeFFCounter+=elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs]*(xLen*yLen*zLen);
   cout << moleculeFFCounter/(frontFaceCounter*xLen*yLen*zLen) << endl;
  }                
  */
    
//  h1->Fill(elemConc[elemCount]);
//  h1->Draw("hist");   
// if ( elemConcUpdate[elemCount] != 1451.7 && elemConcUpdate[elemCount] != 0 ) {
//   cout << elemCount << "     " <<  elemConcUpdate[elemCount] << endl;
// }
// h1->Draw("hist");

/*
  //decide at what time you want to see the conc in all the elems
  //sidenote: can also do it for one element for all time...
  
    if (timePassed==totTime) {
      h1->Fill(elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs]);
      h1->Draw("hist");
    }      */
 

/*    
  //can do a graph for conc vs time for a single element
 
   if (elemCount==45){
    const Int_t dataPtCount = totTime;
    Int_t currTime[dataPtCount];
    Double_t currConc[dataPtCount];

    for (Int_t counter =0; counter<dataPtCount;counter++) {
     Int_t currTime[counter] = timePassed;
     Double_t currConc[counter] = elemConc[elemCount];

    }


    TGraph *gr1 = new TGraph(dataPtCount, currTime, currConc);
    TCanvas *c2 = new TCanvas("c2","Test",200,10,600,400);  
  
    gr1->Draw("AC");
   }  */

   } //end bracket for the 2nd 0<elemCount<divs*divs*divs
   
  } //end bracket for the if statement separating 1st layer from rest
  

 } // end bracket for total time                 
// cout << "it works, yan." << endl; 

} //end bracket 4 int main(void)
  


//ignore the below - it has not been checked fully.                        

/*this function is for calculating the new conc at each elem in the dry
  acrylic area. So not the boundaries. Each elem is surrounded by 6 
  other cubic elems on the left, right, top, bottom, front, and back. 
  The diffusion process, though faciltated by conc, is still randomized.
 function name: cubeElemUpdate
 param1: current conc in the elem we are interested in; given by
  earlier part of the code in the main function
 param2: conc in the elem left to theElem
 param3: conc in the elem right to theElem
 param4: conc in the elem above theElem
 param5: conc in the elem below theElem
 param6: conc in the elem in front of theElem
 param7: conc in the elem behind theElem
 param8: side length of each cube element
 param9: the timestep for each random walk
 param10: diffusion constant for water in acrylic   
 output: we get the new concentration in theElem, in quantity per volume
double cubeElemUpdate (double theElemConc, double leftElemConc,\     
 double rightElemConc, double topElemConc, double botElemConc,\
 double frontElemConc,double backElemConc,double cSize, double timestep,\
 double DcoeffWater)                                
{                                                   
 double elemConcUpdate = 0;     */
/*   elemConcUpdate = DcoeffWater*timestep*(leftElemConc+rightElemConc+ \
  topElemConc+botElemConc+frontElemConc+backElemConc-6*theElemConc) \
  / (cSize*cSize)     
 //need to simulate with rng here; if conc on left > conc on right, then...
 if theElemConc < 
} //end fxn */              
