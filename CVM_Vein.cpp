/**********
  Shuji Ishihara  2026
  
  Simulation codes for Sugimura et al. 2026
  "Cell size control emerges from the vein-dependent coordinated divisions of distinct cell groups in Drosophila wing"
  bioRxiv doi: https://doi.org/10.64898/2025.12.26.696565

**********/

#include<stdio.h>
#include<math.h>
#include<time.h>
#include<stdlib.h>
#include<string.h>
#include<vector>
#include<span> // c++20
using namespace std;
#include <random>
#include <unistd.h>
#include <algorithm>

// Random variables
std::random_device rnd; 
int random_seed = rnd(); // seed
std::mt19937 mt(random_seed);
std::uniform_real_distribution<> unidist(0.0, 1.0);
std::exponential_distribution<> expdist(1.0);

/******************************************
          Main Variables
********************************************/ 


/***  Set outputfile and parameter values  ***/
#include "CVM_Parameters.cpp"

/***  Include libraries for simulation  ***/
#include "lib/head_cell_vertex.cpp"
#include "lib/TimeEvolution.cpp" 
#include "lib/head_cell_div.cpp"
#include "lib/gnuplot.cpp"
#include "lib/libVein.cpp"

/***  Set Scenario dependent cell division ***/
#include "CVM_CellDivScenario.cpp"
#include "CVM_DataOutput.cpp"

/*-------   main function   -------*/
int main( int argc, char *argv[] )
{
  int INIT_CELL_NUMBER;
  FILE *fp; // For output DATAOUTFILE
  FILE *gp; // For gnuplot

  /*** Set initial configuration ***/
  load_initialconfiguraton( argv[1] );   // Load Initial Configuration, defined in lib/IO.cpp  
  Set_BXY(X);           // Boundary reconfiguration, defined in lib/head_cell_vertex.cpp
  Set_Veinregion(0.2);  // Determine Vein region, defined in lib/libVein.cpp

  /*** Set Parameters ***/
  BVo = (LX*LY)/(double)CELL_NUMBER;
  Vo.assign( CELL_NUMBER, BVo ); // set Vo[CELL_NUMBER] 
  BLo = ShpIdx*sqrt(BVo); 
  Lo.assign(CELL_NUMBER, BLo); 

  /****************************
        Time evolution
  ****************************/
  InitialRelaxation( 10.0 ); // relaxation without cell division, just for initialization
  /***  Set cell division plan, initial guess ***/
  INIT_CELL_NUMBER = CELL_NUMBER; // initial cell number 
  for( size_t c=0; c<CELL_NUMBER; ++c ){ 
    cell[c].cell_div_time = Set_CellDIvision_Time( cell[c] ); // initial assignment of planned cell division time, lib/CellDivScenario.cpp
  }

  /*** For Output  ***/
  fp = fopen(DATAOUTFILE,"w");  
  Output_Parameters( random_seed, fp ); // lib/IO.cpp
  gp = getgnuplot( GNUPLOT_USE );       // lib/gnuplot.cpp

  /***   start of time evolutio ***/
  double pt = 0.0; // time 
  for( int tc=0; pt<SIMULATION_TIME; tc++ ){
    static RK2Workspace rk2wk;
    double dx[E_NUM],dy[E_NUM];

    /*** Output gnuplot ***/
    if( tc%TIME_INTERVAL==0 ){
      Set_BXY(X);
      Output_gnuplot( gp, X, Scenario, pt);
    }

    /*** T1 vertex reconnection ***/
    for( size_t i=0; i<E_NUM; ++i )dx[i] = ( *edge[i].x1-*edge[i].x2-edge[i].bx*LX );
    for( size_t i=0; i<E_NUM; ++i )dy[i] = ( *edge[i].y1-*edge[i].y2-edge[i].by*LY );
    for( size_t i=0; i<E_NUM; ++i ){
      if( dx[i]*dx[i]+dy[i]*dy[i]<1.0e-8 ){
        if( Compare_Replace(X,i) ){ // reconnect if energy decreases, defined in lib/head_cell_vertex.cpp
          // fprintf(stderr,"--%f %d (%d-%d)\n",pt,(int)i,edge[i].vtx1,edge[i].vtx2);
        }
      }
    }

    /*** Solve ODE ***/
    RK2( X, pt, rk2wk ); // Runge-Kutta 2nd order , lib/ODE.cpp    
    if( tc%128==0 )Set_BXY(X); 

    /*** Cell_division and target cell area ***/
    Set_Cell_Division_Triggered_By_Neighbours(pt);  // lib/CellDivScenario.cpp
    Cell_Division(pt, INIT_CELL_NUMBER);
    
    /*** data output to DATAOUTFILE  ***/
    if(tc%1000==0){
      DataOutput(fp, pt, INIT_CELL_NUMBER);           
    }
  }/*-----  End of Time evolution  -----*/


  /*** Finalize ***/
  fclose(fp);
  Output_gnuplot( gp, X, Scenario, pt);

  if(gp!=nullptr){
    fflush(gp); // バッファに格納されているデータを吐き出す（必須）
    pclose(gp);
  }

  return 0;
}/*-----  End of main function -----*/













