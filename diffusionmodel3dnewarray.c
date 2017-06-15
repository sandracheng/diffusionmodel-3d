/*SNO+ 3D diffusion model for water in acrylic
 
  Author:  Sandra Cheng. 10184319.  <14wlsc@queensu.ca>

  REVISION HISTORY: N/A

************
  Based on Thomas' initial 3D Diffusion Model code. This is for water in
  an acrylic dogbone.       

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

  The code will use 1000 cubic elements to make a 10x10x10 acrylic cube.
  So number of divisions =10. However this is not hard-coded in.    
  The 10x10x10 acrylic cube may then be normalized to any size acrylic dogbone.
  Origin is still at bottom left front corner of the cube. 
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

  The front face corners will be as above.

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
#include <time.h>   
#include <stdlib.h> 

#include "TCanvas.h"
#include "TH1.h"

//void output(FILE *fp_out, double time, double x, double y, double z, \
// double cubeElem, double avgCubeElem);
//double cubeElemUpdate (double theElemConc, double leftElemConc,\
 double rightElemConc,double topElemConc, double botElemConc, \
 double frontElemConc, double backElemConc,double width, \
 double timestep, double DcoeffWater);                

int diffusionmodel3dnewarray() {                                                                    
 const int divs=4;//n.o. elements aka n.o. divisions in the big acrylic cube
 double blockxLen; //measured x length of acrylic dogbone
 double blockyLen; //measured y length of acrylic dogbone
 double blockzLen; //measured z length of acrylic dogbone
 const int totTime=100;
 //const int timestep=1; 
 int timePassed=0;          
 double maxConc= 1451.7008; //max concentration of water in moles. Placeholder #.
 double elemConc[divs*divs*divs]; //array for current conc in each element
                                  //each element is defined by its position
 double elemConcMaster[divs*divs*divs*totTime];  //updated array that
  //has all the conc for each elem at all times.
 
 double DcoeffWater; //diffusion constant for water in acrylic     
 int i,j,k; //counter for x,y and z positions respectively
 int elemCount; //counter for the total number of elements
 double cSize= 0.5; //sidelength of cube; all cubes are same size.
 int moleculeInCube;

 TCanvas *c1 = new TCanvas;
 Int_t nBins = 1000;
 Double_t lowlim = -1500;
 Double_t uplim = 1500;

 TH1I *h1 = new TH1I("Legend","Title", nBins,lowlim,uplim);
 h1->SetFillColor(kOrange);
 //
// double DcoeffAir; //diffusion constant for saturated acrylic in dry are
// doule airTime; //time it has been out in air; not sure if need
           
 srand(time(NULL));

 //set the concentrations for all the boundary faces and the inside
 //for both the current and update array. boundary=conc, inside=0   
 for (timePassed=1; timePassed<=totTime; timePassed++) {     
  for (elemCount=0; elemCount<divs*divs*divs; elemCount++) {                    

    if (elemCount <= divs*divs-1 || elemCount >= divs*divs*(divs-1)) {
     elemConc[elemCount]= maxConc; //for front and back faces   
     elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs] = maxConc;
    } 
                     
    for (j=0; j<divs; j++) {    
     for (k=0; k<divs*divs; k++) {
      if (k<divs) {
       elemConc[j*divs*divs+k*divs]= maxConc;  //for left face
       elemConcMaster[j*divs*divs+k*divs+(timePassed-1)*divs*divs*divs]= maxConc;

       elemConc[j*divs*divs+k*divs+divs-1] = maxConc;  //for right face
       elemConcMaster[j*divs*divs+k*divs+divs-1+(timePassed-\
        1)*divs*divs*divs]=maxConc;

       elemConc[j*divs*divs+k] = maxConc; //for bottom face
       elemConcMaster[j*divs*divs+k+(timePassed-1)*divs*divs*divs] = maxConc;
      }
      else if (k >= divs*(divs-1)) {
       elemConc[j*divs*divs+k]= maxConc; //for top face
       elemConcMaster[j*divs*divs+k+(timePassed-1)*divs*divs*divs] = maxConc;
      }
     }
    }
  //cout << "elem: " << elemCount << " conc: " << elemConc[elemCount]<< endl; works.                  
  }                                                            
  
  for (elemCount=0; elemCount<divs*divs*divs; elemCount++) {                  
  
//   moleculeInCube = elemConc[elemCount]/(cSize*cSize*cSize); 
//   printf("elem:%d; elem conc: %f\n", elemCount, elemConc[elemCount]); // works!           
 
   //make random walking part only for not a boundary cube for all time
   //random walk need to occur for the boundary parts for t=1
   if (elemConcMaster[elemCount]!=maxConc) {
    moleculeInCube = elemConc[elemCount]*(cSize*cSize*cSize); 
   
    for (int moleculeCount=1; moleculeCount<= moleculeInCube; moleculeCount++) { 
 /*need to know position of the 26 surrounding cube elements. will name it 
   by the number scheme as above, with elements are labelled 0 to 26.
   the central cube number is 13                   
                  
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
   cube26= elemConc[elemCount+divs+1+divs*divs];                           */  

    int r = rand() % 27; //random number between 0 and 26; 27 possibilities 
    //for each molecule to move into one of the adjacent 26 cubes or stay.
    //all if statements match the random number to the cube it moves to

     if (r==0) {                                    
     //new conc of cube 0 relative to the elem=previous conc+(1molecule/totVol)
     //new conc of the elem=previous conc-(1molecule/totVol).
     //update master array accordingly for both.
      elemConc[elemCount-divs*divs-1-divs]= elemConc[elemCount-divs*divs-1-\
       divs]+1/(cSize*cSize*cSize); 
      elemConcMaster[elemCount-divs*divs-1-divs+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount-divs*divs-1-divs];           
      
      elemConc[elemCount]= elemConc[elemCount]-1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs]=\
       elemConc[elemCount];
     }                                                   
     else if (r==1) {                                    
      elemConc[elemCount-divs*divs-divs] = elemConc[elemCount-divs*divs-divs]\
       +1/(cSize*cSize*cSize);
      elemConcMaster[elemCount-divs*divs-divs+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount-divs*divs-divs];  
         
      elemConc[elemCount] = elemConc[elemCount]-1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount-divs*divs-divs];  
     }                 
     else if (r==2) {                                
      elemConc[elemCount-divs*divs-divs+1]=elemConc[elemCount-divs*divs-divs+1]\
       +1/(cSize*cSize*cSize);
      elemConcMaster[elemCount-divs*divs-divs+1+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount-divs*divs-divs+1];
 
      elemConc[elemCount]= elemConc[elemCount]-1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount];
     }    
     else if (r==3) {                                    
     elemConc[elemCount-divs*divs-1] = elemConc[elemCount-divs*divs-1]\
      +1/(cSize*cSize*cSize);
     elemConcMaster[elemCount-divs*divs-1+(timePassed-1)*divs*divs*divs]\
      =elemConc[elemCount-divs*divs-1];

     elemConc[elemCount] = elemConc[elemCount]-1/(cSize*cSize*cSize);
     elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount];
     }
     else if (r==4) {                                    
     elemConc[elemCount-divs*divs] = elemConc[elemCount-divs*divs]\
      +1/(cSize*cSize*cSize);
     elemConcMaster[elemCount-divs*divs+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount-divs*divs];

     elemConc[elemCount] = elemConc[elemCount]-1/(cSize*cSize*cSize);
     elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount];
     }
     else if (r==5) {                                    
      elemConc[elemCount-divs*divs+1] = elemConc[elemCount-divs*divs+1]\
       +1/(cSize*cSize*cSize); 
      elemConcMaster[elemCount-divs*divs+1+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount-divs*divs+1];

      elemConc[elemCount] = elemConc[elemCount]-1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount];
     }
     else if (r==6) {                                    
      elemConc[elemCount-divs*divs-1+divs] = elemConc[elemCount-divs*divs-1+divs]\
       +1/(cSize*cSize*cSize);
      elemConcMaster[elemCount-divs*divs+1+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount-divs*divs+1];

      elemConc[elemCount]=elemConc[elemCount]-1/(cSize*cSize*cSize);
      elemConcMaster[elemCount-divs*divs+1+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount-divs*divs+1];
     }
     else if (r==7) {                                    
      elemConc[elemCount-divs]=elemConc[elemCount-divs]+1/(cSize*cSize*cSize);
      elemConcMaster[elemCount-divs+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount-divs];

      elemConc[elemCount]=elemConc[elemCount]-1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount];
     }
     else if (r==8) {                                    
      elemConc[elemCount-divs*divs+1+divs] = elemConc[elemCount-divs*divs+1+divs]\
       +1/(cSize*cSize*cSize);
      elemConcMaster[elemCount-divs*divs+1+divs+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount-divs*divs+1+divs];

      elemConc[elemCount] = elemConc[elemCount]-1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount];
     }
     else if (r==9) {                                    
      elemConc[elemCount-1-divs] = elemConc[elemCount-1-divs]\
       +1/(cSize*cSize*cSize);
      elemConcMaster[elemCount-divs-1+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount-divs-1];

      elemConc[elemCount] = elemConc[elemCount]-1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount];
     }                                            
     else if (r==10) {                                    
      elemConc[elemCount-divs] = elemConc[elemCount-divs]\
       +1/(cSize*cSize*cSize);
      elemConcMaster[elemCount-divs+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount-divs];

      elemConc[elemCount]=elemConc[elemCount]-1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount];
     }
     else if (r==11) {                                    
      elemConc[elemCount+1-divs] = elemConc[elemCount+1-divs]\
       +1/(cSize*cSize*cSize);
      elemConcMaster[elemCount-divs+1+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount-divs+1];

      elemConc[elemCount] = elemConc[elemCount]-1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount];
     }
     else if (r==12) {                                    
      elemConc[elemCount-1] = elemConc[elemCount-1]\
       +1/(cSize*cSize*cSize);
      elemConcMaster[elemCount-1+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount-1];

      elemConc[elemCount] = elemConc[elemCount]-1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount];
     }                 
     else if (r==14) {                                    
      elemConc[elemCount+1] = elemConc[elemCount+1] \
       +1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+1+(timePassed-1)*divs*divs*divs] \
       =elemConc[elemCount+1];

      elemConc[elemCount] = elemConc[elemCount]-1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs] \
       =elemConc[elemCount];
     }
     else if (r==15) {                                    
      elemConc[elemCount+divs-1] = elemConc[elemCount+divs-1] \
       +1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+divs-1+(timePassed-1)*divs*divs*divs] \
       =elemConc[elemCount+divs-1];

      elemConc[elemCount] = elemConc[elemCount]-1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs] = \ 
       elemConc[elemCount];
     }
     else if (r==16) {                                    
      elemConc[elemCount+divs] = elemConc[elemCount+divs] \
       +1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+divs+(timePassed-1)*divs*divs*divs] \
       =elemConc[elemCount+divs];

      elemConc[elemCount] = elemConc[elemCount]-1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs] \
       =elemConc[elemCount];
     }
     else if (r==17) {                                    
      elemConc[elemCount+divs+1] = elemConc[elemCount+divs+1] \
       +1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+divs+1+(timePassed-1)*divs*divs*divs] \
       =elemConc[elemCount+divs+1];

      elemConc[elemCount] = elemConc[elemCount]-1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs] \
       =elemConc[elemCount];
     }
     else if (r==18) {                                    
      elemConc[elemCount-1-divs+divs*divs] = elemConc[elemCount-1+ \ 
       divs*divs-divs]+1/(cSize*cSize*cSize);
      elemConcMaster[elemCount-1-divs+divs*divs+(timePassed- \ 
       1)*divs*divs*divs]=elemConc[elemCount-1-divs+divs*divs];

      elemConc[elemCount] =elemConc[elemCount]-1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs] \
       =elemConc[elemCount];
     }
     else if (r==19) {                                    
      elemConc[elemCount+divs*divs-divs] = elemConc[elemCount+divs*divs-divs] \
       +1/(cSize*cSize*cSize);
      elemConcMaster[elemCount-divs+divs*divs+(timePassed-1)*divs*divs*divs] \
       =elemConc[elemCount-divs+divs*divs];

      elemConc[elemCount] =elemConc[elemCount]-1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount];
     }
     else if (r==20) {                                   
      elemConc[elemCount+1+divs*divs-divs] = elemConc[elemCount+ \ 
       1+divs*divs-divs]+1/(cSize*cSize*cSize);
      elemConcMaster[elemCount-divs+divs*divs+1+ \ 
       (timePassed-1)*divs*divs*divs]=elemConc[elemCount-divs+divs*divs+1];

      elemConc[elemCount] =elemConc[elemCount]-1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount];
     }
     else if (r==21) {                                    
      elemConc[elemCount-1+divs*divs] = elemConc[elemCount-1+divs*divs] \
       +1/(cSize*cSize*cSize);
      elemConcMaster[elemCount-1+divs*divs+(timePassed-1)*divs*divs*divs] \
       =elemConc[elemCount-1+divs*divs];

      elemConc[elemCount] = elemConc[elemCount]-1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs] \
       =elemConc[elemCount];
     }
     else if (r==22) {                                    
      elemConc[elemCount+divs*divs] = elemConc[elemCount+divs*divs] \
       +1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+divs*divs+(timePassed-1)*divs*divs*divs] \
       =elemConc[elemCount];

      elemConc[elemCount] = elemConc[elemCount]-1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs] \
       =elemConc[elemCount];
     }                               
     else if (r==23) {                                   
      elemConc[elemCount+divs*divs+1] = elemConc[elemCount+divs*divs+1] \
       +1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+divs*divs+1+(timePassed-1)*divs*divs*divs] \
       =elemConc[elemCount+divs*divs+1];

      elemConc[elemCount] = elemConc[elemCount]-1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount];
     }
     else if (r==24) {                                    
      elemConc[elemCount-divs+divs*divs+1] = elemConc[elemCount+divs*divs-divs+1]\
       +1/(cSize*cSize*cSize);
      elemConcMaster[elemCount-divs+divs*divs+1+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount-divs+divs*divs+1];

      elemConc[elemCount] = elemConc[elemCount]-1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount];
     }
     else if (r==25) {                                    
      elemConc[elemCount+divs*divs+divs] = elemConc[elemCount+divs*divs+divs]\
       +1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+divs+divs*divs+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount+divs+divs*divs];

      elemConc[elemCount] = elemConc[elemCount]-1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount];
     }
     else if (r==26) {                                    
      elemConc[elemCount+divs+1+divs*divs] = elemConc[elemCount+divs+1+divs*divs]\
       +1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+1+divs+divs*divs+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount+divs+1+divs*divs];

      elemConc[elemCount] = elemConc[elemCount]-1/(cSize*cSize*cSize);
      elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount];
     }
     else if (r==13){   //for r=13, where the molecule stays put                                  
      elemConc[elemCount]= elemConc[elemCount]; 
      elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs]\
       =elemConc[elemCount];
     }
     else {
     break;
     }
 cout << r << endl;
// cout << "elemcount:" << elemCount << " r:"  << r << " conc:" << elemConc[elemCount] << " time: " << timePassed << endl;  
//     h1->Fill(elemConc[elemCount]);
//     h1->Draw("hist");   
  
    } //end bracket for the molecule count per elem    
 //  h1->Fill(elemConc[elemCount]);
//  h1->Draw("hist");   
// if ( elemConcUpdate[elemCount] != 1451.7 && elemConcUpdate[elemCount] != 0 ) {
//   cout << elemCount << "     " <<  elemConcUpdate[elemCount] << endl;
// }
// h1->Draw("hist");
   } //end bracket for if statement to separate nonboundary&boundary  
      
   //cout << "time: " << timePassed << " " << elemConcMaster[elemCount+(timePassed-1)*divs*divs*divs] << endl; 
   
   //decide at what time you want to see the conc in all the elems
   //sidenote: can also do it for one element for all time...
   if (timePassed==totTime) {
    h1->Fill(elemConc[elemCount]);
    h1->Draw("hist");
    cout << elemCount << " " << elemConc[elemCount] << endl; 
   }  

                                                            
  } //end bracket for the 0<elemCount<divs*divs*divs                         
 

 } // end bracket for total time                 
 

} //end bracket 4 int main(void)
                





                         

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
