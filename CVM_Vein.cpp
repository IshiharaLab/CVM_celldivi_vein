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
#include "CellDivScenario.cpp"


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
    Cell_Division_Triggered_By_Neighbours(pt);  // lib/CellDivScenario.cpp

    for( size_t dc=0; dc<INIT_CELL_NUMBER; ++dc ){ // cell division
      switch(cell[dc].cell_type){
        case 'v': continue;
        case 'i': 
          if( cell[dc].cell_div_time < pt ){             
            div_by_nucleation++;
            cell_div_counter++;                       
            cell[dc].cell_type = 'm'; // transition to mitotic phase
          }
          continue;
        case 'm':
           Vo[dc] = 1.5*BVo;
           Lo[dc] = ShpIdx*sqrt(1.5*BVo); 
          if( cell[dc].cell_div_time+MITOTIC_TIME < pt ){ // cell division
            PRE_CELL_DIVISION( X, Vo, dc, pt );       
            fprintf( stderr, " cell div. %f %d %d\n", pt,cell_div_counter,INIT_CELL_NUMBER );
          }
          continue;
        case 'd': continue;
        default:  ERROR("Initial cell_type not defined, This line should not be reached\n");
      }
    }
    
    /*** data output to DATAOUTFILE  ***/
    if(tc%1000==0){
      // time and energy, # of accumulated number of cell divisions
      double energy = all_energy( X, edge, cell );
      fprintf( fp, "%.2f %f %d %d %d   %d  ", pt, energy, cell_div_counter,div_by_nucleation,div_by_triggered,div_by_nucleation+div_by_triggered );

      //number of cells in type
      int i_cell = 0, m_cell = 0, d_cell = 0;
      for(size_t c=0;c<CELL_NUMBER;++c ){
        switch( cell[c].cell_type ){
          case 'v': break;
          case 'i': i_cell++; break;
          case 'm': m_cell++; break;
          case 'd': d_cell++; break;
        }
      }
      fprintf( fp, "%d %d %d ",i_cell,m_cell,d_cell);

      // Mechanical state, statistic
      if( i_cell==0 ){
        fprintf( fp, " - - " );
      }else{
        double icarea_mean = 0.0, icarea_std = 0.0;
        for(size_t c=0;c<CELL_NUMBER;++c){
          if(cell[c].cell_type=='i'){
            double area = cell_area(cell[c],x,y);
            icarea_mean += area;
            icarea_std += area*area;
          }
        }
        icarea_mean /= (double)i_cell;
        icarea_std /= (double)i_cell;
        icarea_std -= icarea_mean *icarea_mean;
        fprintf(fp, " %f %f ",icarea_mean, sqrt(icarea_std));
      }

      if( m_cell==0 ){
        fprintf(fp, " - - ");
      }else{
        double mcarea_mean = 0.0, mcarea_std = 0.0;
        for( size_t c=0; c<CELL_NUMBER; ++c ){
          if(cell[c].cell_type=='m'){
            double area = cell_area(cell[c],x,y);
            mcarea_mean += area;
            mcarea_std += area*area;
          }
        }
        mcarea_mean /= (double)m_cell;
        mcarea_std /= (double)m_cell;
        mcarea_std -= mcarea_mean *mcarea_mean;
        fprintf(fp," %f %f ",mcarea_mean, sqrt(mcarea_std));
      }
  
      if(i_cell==0){ fprintf(fp, " - "); }
      else{
        double i_volenergy = 0.0;
        for(int i=0;i<INIT_CELL_NUMBER;++i){
          if( cell[i].cell_type=='i' ){
            double area = cell_area( cell[i], x, y );
            i_volenergy += 0.5*KV*(area-Vo[i])*(area-Vo[i]);
          }
        }      
        //i_volenergy /=(double)i_cell;
        fprintf( fp, " %f %f ", i_volenergy, i_volenergy/(double)i_cell);
      }
      fprintf(fp, "\n");
      fflush(fp);
      //printf("%f\n",pt);
     
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













