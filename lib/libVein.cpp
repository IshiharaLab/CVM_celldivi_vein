extern int cell_div_counter, div_by_nucleation, div_by_triggered;
extern const double MitoticTime; // mitotic phase duration

double Set_CellDIvision_Time(CELL dcell); // Planned cell div time: called once as initial setting  
void Cell_Division_Triggered_By_Neighbours(double const pt);

void Set_Veinregion(double margin)
{
  for( size_t c=0; c<CELL_NUMBER; ++c ){
    CELL_MOMENT cmt; // CELL_MOMENT geometry.cpp
    cmt.set_Cell_Moment( cell[c] );        
    if( cmt.cy < margin*LY || cmt.cy> (1-margin)*LY ){ cell[c].cell_type = 'v'; }
    else{ cell[c].cell_type = 'i'; }  
  }
} 
