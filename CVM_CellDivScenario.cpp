/*************************************************************
  Control cell division rates depending on different Scenarios 
*************************************************************/

extern const char Scenario; //  defined in CVM_Parameters.cpp
// extern int cell_div_counter = 0, div_by_nucleation = 0, div_by_triggered = 0;
extern const double MITOTIC_TIME;
inline double celldiv_space(const double cy){ return 1.0-fabs(1-2*(cy/LY)); }

double Set_CellDIvision_Time( CELL dcell ) // Planned cell div time: called once as initial setting  
{
  CELL_MOMENT cmt; // cell moment, ./lib/geometry.cpp
  cmt.set_Cell_Centers(dcell); // centroid position
  cmt.cy = fmod(cmt.cy+LY,LY); // 

  switch( dcell.cell_type ){
    case 'v': return INFINITY; // never divided
    case 'i':
      switch(Scenario){
        case 'B': return 800*celldiv_space(cmt.cy)*expdist(mt); // distance deoendent, stochastic Sc B
        case 'C': return 150*celldiv_space(cmt.cy) - 40.0; // distance deoendent, deterministic  Sc C
        case 'A': return INFINITY;// do nothing
      }
    case 'm': ERROR("Initial cell_type = 'm', This line should not be reached\n");
    case 'd': ERROR("Initial cell_type = 'd', This line should not be reached\n");
    default:  ERROR("Initial cell_type not defined, This line should not be reached\n");  
  }
}

void Set_Cell_Division_Triggered_By_Neighbours( double const pt )
{    
  if( Scenario == 'C' )return;
  double VCELLDIVRATE = 0.0, NCELLDIVRATE = 0.0;
  if( Scenario == 'A' ){ VCELLDIVRATE=0.2*Dt;  NCELLDIVRATE = 0.01*Dt; }
  else if( Scenario == 'B' ){ VCELLDIVRATE = 0.0*Dt; NCELLDIVRATE = 0.01*Dt; }
 
  vector<char> nst( CELL_NUMBER,'z' );
  for( size_t c=0; c<CELL_NUMBER; ++c ){    
    if(cell[c].cell_type == 'i'){
      int *nc = cell[c].nc;
      for( size_t j=0;j<cell[c].jnum;++j ){
        if( cell[nc[j]].cell_type == 'v' && unidist(mt)<VCELLDIVRATE ){ // neighbour of vein cell  
          nst[c] = 'm';
          cell[c].cell_div_time = pt;
          div_by_nucleation++;
          cell_div_counter++;
          break;
        }else if( (cell[nc[j]].cell_type == 'm' || cell[nc[j]].cell_type == 'd') && unidist(mt)<NCELLDIVRATE ){          
          nst[c] = 'm';
          cell[c].cell_div_time = pt;
          div_by_triggered++;
          cell_div_counter++;
          break;
        }
      }
    }
  }

  for( size_t c=0; c<CELL_NUMBER; ++c ){
    if( nst[c] == 'm' ){ cell[c].cell_type = 'm'; }
  }
}



// perform cell division
void Cell_Division(const double pt, const int INIT_CELL_NUMBER)
{
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
}
