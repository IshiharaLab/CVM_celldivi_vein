
double cell_area( const CELL& cell, const double *x, const double *y )
{
  int pn = cell.jnum-1;
  double area=0.0;
  for(size_t j=0;j<pn;++j){
    area+=(x[cell.vtx[j]]+cell.bx[j]*LX)*(y[cell.vtx[j+1]]+cell.by[j+1]*LY)-(x[cell.vtx[j+1]]+cell.bx[j+1]*LX)*(y[cell.vtx[j]]+cell.by[j]*LY);
  }
  area+=(x[cell.vtx[pn]]+cell.bx[pn]*LX)*(y[cell.vtx[0]]+cell.by[0]*LY)-(x[cell.vtx[0]]+cell.bx[0]*LX)*(y[cell.vtx[pn]]+cell.by[pn]*LY);
  return area/2;
}


double boundary_length( const CELL &tcell, const double *x, const double *y ) // not used
{
  int pn = tcell.jnum-1;
  double lgth=0.0;
  double dx[NMAX],dy[NMAX];
  for( int j=0; j<pn; ++j ){ dx[j] = (x[tcell.vtx[j+1]]+tcell.bx[j+1]*LX)-(x[tcell.vtx[j]]+tcell.bx[j]*LX); }
  dx[pn] = (x[tcell.vtx[0]]+tcell.bx[0]*LX) - (x[tcell.vtx[pn]]+tcell.bx[pn]*LX);
  for( int j=0; j<pn; ++j ){ dy[j] = (y[tcell.vtx[j+1]]+tcell.by[j+1]*LY)-(y[tcell.vtx[j]]+tcell.by[j]*LY); }
  dy[pn] = (y[tcell.vtx[0]]+tcell.by[0]*LY) - (y[tcell.vtx[pn]]+tcell.by[pn]*LY);
  for(size_t i=0;i<tcell.jnum;++i){ lgth += sqrt(dx[i]*dx[i]+dy[i]*dy[i]);  }
  return lgth;
}

struct CELL_MOMENT{
  double area,tarea[NMAX];
  double cx,cy;
  double mxx,mxy,myy;
  double Nxx,Nxy,Nyy; // Nxx not normalized moment; not stubstaracted by cx*cx;
  double area_dx[NMAX],area_dy[NMAX];
  double ncx_dx[NMAX],ncx_dy[NMAX],ncy_dx[NMAX],ncy_dy[NMAX];
  double nxx_dx[NMAX],nxx_dy[NMAX],nxy_dx[NMAX],nxy_dy[NMAX],nyy_dx[NMAX],nyy_dy[NMAX];
  double l1,l2;

  CELL_MOMENT(){} // default;
  CELL_MOMENT(const CELL& e_cell){  set_Cell_Moment(e_cell); }

  void set_Cell_Centers(const CELL &e_cell){
    double *x=X.data(),*y=X.data()+V_NUM;    
    area = 0.0; cx=0.0; cy=0.0;

    for(size_t j=0;j<e_cell.jnum;++j){
      int nj = (j==e_cell.jnum-1 ? 0: j+1 );
      double x1,x2,y1,y2,tarea;
      x1 = x[e_cell.vtx[j] ]+e_cell.bx[j ]*LX;
      x2 = x[e_cell.vtx[nj]]+e_cell.bx[nj]*LX;
      y1 = y[e_cell.vtx[j] ]+e_cell.by[j ]*LY;
      y2 = y[e_cell.vtx[nj]]+e_cell.by[nj]*LY;
      tarea = (x1*y2-x2*y1);
      area += tarea/2;
      cx += tarea*(x1+x2);
      cy += tarea*(y1+y2);    
    }
    cx/=(6*area);
    cy/=(6*area);
  }

  void set_Cell_Moment( const CELL &e_cell ){
    double *x=X.data(),*y=X.data()+V_NUM;    
    area = 0.0; cx=0.0; cy=0.0; mxx=0.0; mxy=0.0; myy=0.0;

    //  moment計算は　"On the Calculation of Arbitrary Moments of Polygons" Carsten Steger　参照
    for(size_t j=0;j<e_cell.jnum;++j){
      int nj = (j!=e_cell.jnum-1 ? j+1 : 0 );
      double x1,x2,y1,y2,tarea;
      x1 = x[e_cell.vtx[j] ]+e_cell.bx[j ]*LX;
      x2 = x[e_cell.vtx[nj]]+e_cell.bx[nj]*LX;
      y1 = y[e_cell.vtx[j] ]+e_cell.by[j ]*LY;
      y2 = y[e_cell.vtx[nj]]+e_cell.by[nj]*LY;
      tarea = (x1*y2-x2*y1);
      area += tarea/2;
      cx += tarea*(x1+x2);
      cy += tarea*(y1+y2);
      mxx += tarea*(x1*x1+x1*x2+x2*x2);
      mxy += tarea*(2*x1*y1+x1*y2+x2*y1+2*x2*y2);
      myy += tarea*(y1*y1+y1*y2+y2*y2);
    }
    cx/=(6*area);
    cy/=(6*area);
    mxx/=(12*area); mxx-=cx*cx;
    mxy/=(24*area); mxy-=cx*cy;
    myy/=(12*area); myy-=cy*cy;    
  }

}; 

//vector<CELL_MOMENT> cell_moment; // donot use //



