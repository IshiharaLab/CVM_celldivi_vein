#include <optional>
#include <string>

FILE* open_gnuplot(double LX, double LY) // for initialization 
{
    FILE *gp = popen("gnuplot -persist", "w");
    if (!gp) return NULL;
    fprintf(gp, "unset key\n");
    fprintf(gp, "set xrange [0:%f]\n", LX);
    fprintf(gp, "set yrange [%f:%f]\n", 0.1*LY, 0.9*LY);
    fprintf(gp, "set size ratio 0.8\n");
    fprintf(gp, "unset tics\n");
    fprintf(gp, "set term qt font \"Arial 20\"\n" );
    return gp;
}

FILE* open_gnuplot(double LX, double LY, const char* OUTPUT_GIFFILE)
{
    FILE *gp = open_gnuplot(LX, LY);
    if (!gp) return nullptr;
    fprintf(gp, "set term gif animate delay 5\n");
    fprintf(gp, "set output '%s'\n", OUTPUT_GIFFILE);
    return gp;
}

FILE* getgnuplot(int gpcontrol) 
{
  if ( gpcontrol == 1 ) { return open_gnuplot(LX, LY); }
  else if( gpcontrol == 2){ return open_gnuplot(LX, LY, OUTGIFFILE); }
  else{ return nullptr; }
}

void boundary_check(const CELL tcell, int &bx, int &by)
{
  bx = 0;  by = 0;
  for(int j=0; j<tcell.jnum; ++j){    
    if(tcell.bx[j]!=0){ bx = tcell.bx[j]; }
    if(tcell.by[j]!=0){ by = tcell.by[j]; }
  }
}




void print_cell_polygon(FILE* gp, const CELL& c, const double *x,const double *y, double LX, double LY)
{
  // original image
  for (int j = 0; j < c.jnum; ++j) {
    fprintf(gp, "%f %f\n",x[c.vtx[j]] + LX * c.bx[j],y[c.vtx[j]] + LY * c.by[j]);
  }
  fprintf(gp,"\n");

  // periodic image
  int bx, by;
  boundary_check(c, bx, by);
  for (int j = 0; j < c.jnum; ++j) {
    fprintf(gp, "%f %f\n",x[c.vtx[j]] + LX * (c.bx[j] - bx),y[c.vtx[j]] + LY * (c.by[j] - by));
  }
  fprintf(gp,"\n");
}

auto color_of = [](char type) -> const char* {
    switch (type) {
        case 'v': return "blue";
        case 'i': return "gray";
        case 'm': return "green";
        case 'd': return "yellow";
        default:  return nullptr;
    }
};


/* output */
void Output_gnuplot( FILE *gp, vector<double> &X, const char Sc, const double ptime )
{
  if( gp == nullptr) return;
  double *x = X.data(), *y = X.data()+V_NUM;
  Set_BXY(X);
  
  fprintf(gp, "set title 'Scenario %c : time = %4.1f  min.' offset 0, -1.3 font 'Arial, 20' \n",Sc, 2*ptime);
  fprintf(gp, "set style fill\n");
  //fprintf(gp, "set style fill solid 0.2\n");

  fprintf(gp, "plot ");
  for( int i = 0; i < CELL_NUMBER; ++i) {
    fprintf(gp, "'-' with filledcurves closed fc rgb \"%s\" notitle, ", color_of(cell[i].cell_type));     
  }
  fprintf(gp,"'-' with lines lt -1 \n");

  for(size_t i=0;i<CELL_NUMBER;++i){
    print_cell_polygon( gp, cell[i], x, y, LX, LY);
    fprintf(gp, "e\n");
  }

  for( size_t i=0;i<E_NUM;++i ){
    if(edge[i].bx==0 && edge[i].by==0){
      fprintf(gp,"%f %f\n%f %f\n\n",x[edge[i].vtx1],y[edge[i].vtx1],x[edge[i].vtx2],y[edge[i].vtx2]);
    }else if(edge[i].bx!=0 &&edge[i].by==0){
      fprintf(gp,"%f %f\n%f %f\n\n",x[edge[i].vtx1],y[edge[i].vtx1],x[edge[i].vtx2]+LX*edge[i].bx,y[edge[i].vtx2]);
      fprintf(gp,"%f %f\n%f %f\n\n",x[edge[i].vtx1]-LX*edge[i].bx,y[edge[i].vtx1],x[edge[i].vtx2],y[edge[i].vtx2]);
    }else if(edge[i].bx==0 &&edge[i].by!=0){
      fprintf(gp,"%f %f\n%f %f\n\n",x[edge[i].vtx1],y[edge[i].vtx1],x[edge[i].vtx2],y[edge[i].vtx2]+LY*edge[i].by);
      fprintf(gp,"%f %f\n%f %f\n\n",x[edge[i].vtx1],y[edge[i].vtx1]-LY*edge[i].by,x[edge[i].vtx2],y[edge[i].vtx2]);
    }else{
      fprintf(gp,"%f %f\n%f %f\n\n",x[edge[i].vtx1],y[edge[i].vtx1],x[edge[i].vtx2]+LX*edge[i].bx,y[edge[i].vtx2]+LY*edge[i].by);
      fprintf(gp,"%f %f\n%f %f\n\n",x[edge[i].vtx1]-LX*edge[i].bx,y[edge[i].vtx1]-LY*edge[i].by,x[edge[i].vtx2],y[edge[i].vtx2]);
    }
  }
  fprintf(gp,"e\n");

  fflush(gp);
  usleep(10000);
}
