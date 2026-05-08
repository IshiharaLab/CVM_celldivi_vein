
void Output_Parameters(int random_seed, FILE *fp)
{
  fprintf(fp,"# Scenario = %c\n", Scenario); //
  fprintf(fp,"# seed = %d\n", random_seed); //  乱数生成のseed, 出力しとく
  fprintf(fp,"## LX*LY = %f x %f \n",LX,LY);
  fprintf(fp,"## KV = %e \n",KV);
  fprintf(fp,"## Vo = %e \n",BVo);
  fprintf(fp,"## KL = %e \n",KL);
  fprintf(fp,"## ShpIde = %e \n",ShpIdx);  
}



void DataOutput(FILE *fp, const double pt, const int INIT_CELL_NUMBER)
{
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
  
}