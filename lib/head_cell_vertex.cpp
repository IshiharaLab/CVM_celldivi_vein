
/***************************
        MACRO
***************************/
#define ERROR(error_text)  { fprintf(stderr,"%s\n",error_text);  exit(10); }
#define SWAP(type,a,b) { type temp = a; a = b; b = temp; }



/***************************
      MAIN STRUCTUREs
***************************/
struct VERTEX{
  double *x,*y;
  int nv1,nv2,nv3; /* neighbour vertoices */
  int e1,e2,e3; /* neighbour edge */
  int c1,c2,c3;	/* neighbour cellxf */
  
  void reset() { // for initialization, called by void set_vertex()
    nv1 = nv2 = nv3 = -1;
    e1  = e2  = e3  = -1;
    c1  = c2  = c3  = -1;    
  }
};

struct EDGE{
  int vtx1,vtx2;
  double *x1,*y1,*x2,*y2;
  char bx,by;  /** (-1,0,1) **/
  int c1,c2;

  void set_ptrs(char tbx, char tby){
    bx = tbx;
    by = tby;
    x1 = X.data()+vtx1;
    y1 = X.data()+V_NUM+vtx1;
    x2 = X.data()+vtx2;
    y2 = X.data()+V_NUM+vtx2;
  }
};


#define NMAX 32
struct CELL{
  int jnum;
  int vtx[NMAX]; // vectorを使った方が良い？　: なんか遅くなる。。。  
  char bx[NMAX],by[NMAX]; // for boundary condition  
  int nc[NMAX]; // neighboring cells  !! Order is not aligend !! 
  // Variables for cell type and divisions process
  char cell_type; // 'v' vein  'i' intervein/initial  'm' mitotic phase  'd' divided intervein
  int cell_div_count; // counting division events
  double cell_div_time; // cell division time
  CELL(): cell_div_count(0), cell_div_time(0.0) {}  // 初期値
  //CELL_CONTROL ctl;
};

typedef struct VERTEX VERTEX;
typedef struct EDGE EDGE;
typedef struct CELL CELL;
typedef struct CELL_MOMENT CELL_MOMENT; // 使っていない

vector<VERTEX> vtx;
vector<EDGE> edge;
vector<CELL> cell;

 void Assign_Variables()
 {
  N = 2*V_NUM;
  X.resize(N);   
  Vo.resize(CELL_NUMBER);
  Lo.resize(CELL_NUMBER);
  vtx.resize(V_NUM);
  edge.resize(E_NUM);
  cell.resize(CELL_NUMBER);
  x = X.data();
  y = X.data()+V_NUM;
 }



/***************************
        FUNCTION  
***************************/
// defined in this file
void Set_Vertex( const vector<double> &X );
void Set_BXY( vector<double> &X );
bool Compare_Replace( vector<double> &X, const int ex );
int Replace( const int ex,vector<double> &X, vector<VERTEX> &e_vtx, vector<EDGE> &e_edge, vector<CELL> &e_cell );

// defined in IO.cpp
void Set_initialconfiguraton(char *inputfile);
void output(vector<double> &X,FILE *out_file);
void out_geometry(const vector<double> &X,FILE *out_file);
void OutputTensionPressure(const vector<double> &X,FILE *out_file);

// defined in geometry.cpp
double cell_area( const CELL& cell, const double *x, const double *y );
double boundary_length( const CELL &tcell, const double *x, const double *y );
void set_Cell_Moment( CELL_MOMENT &tcell_moment, const CELL &e_cell );

// definied in ODE.cpp
extern void onestep(vector<double> &X, vector<double> &kX, const double pt);
extern double all_energy(const vector<double> &X, const vector<EDGE> &e_edge, const vector<CELL> &e_cell);

// for vein 
void Set_NeighbourCells_RandomOrder( int cid, CELL &e_cell ); // set cell.nc[]


#include "IO.cpp"
#include "geometry.cpp"

/************************************************************************/

void Set_Vertex(const vector<double> &X)
{
  for(size_t i=0;i<V_NUM;++i){
    vtx[i].x=x+i;
    vtx[i].y=y+i;
    vtx[i].reset();
  }
 
  for(size_t i=0;i<E_NUM;++i){
    if(vtx[edge[i].vtx1].nv3==-1){
      vtx[edge[i].vtx1].nv3=edge[i].vtx2;
      vtx[edge[i].vtx1].e3=i;
    }else if(vtx[edge[i].vtx1].nv2==-1){
      vtx[edge[i].vtx1].nv2=edge[i].vtx2;
      vtx[edge[i].vtx1].e2=i;
    }else{
      vtx[edge[i].vtx1].nv1=edge[i].vtx2;
      vtx[edge[i].vtx1].e1=i;
    }
    if(vtx[edge[i].vtx2].nv3==-1){
      vtx[edge[i].vtx2].nv3=edge[i].vtx1;
      vtx[edge[i].vtx2].e3=i;
    }else if(vtx[edge[i].vtx2].nv2==-1){
      vtx[edge[i].vtx2].nv2=edge[i].vtx1;
      vtx[edge[i].vtx2].e2=i;
    }else{
      vtx[edge[i].vtx2].nv1=edge[i].vtx1;
      vtx[edge[i].vtx2].e1=i;
    }
  }

  for(size_t i=0;i<CELL_NUMBER;++i){
    for(size_t j=0;j<cell[i].jnum;++j){
      if(vtx[cell[i].vtx[j]].c3==-1){vtx[cell[i].vtx[j]].c3=i;}
      else if(vtx[cell[i].vtx[j]].c2==-1){vtx[cell[i].vtx[j]].c2=i;}
      else{vtx[cell[i].vtx[j]].c1=i;}
    }
  }

  for(size_t i=0; i<E_NUM; ++i ){
    edge[i].c1=-1;
    edge[i].c2=-1;
  }

  for(size_t i=0;i<CELL_NUMBER;++i){
    for(size_t j=0;j<cell[i].jnum;++j){
      int j1,j2,k;
      j1 = cell[i].vtx[j];
      j2 = ( j==cell[i].jnum-1 ? cell[i].vtx[0] : cell[i].vtx[j+1] );
      for( k=0; k<E_NUM; k++ ){
        if( (edge[k].vtx1==j1 && edge[k].vtx2==j2) || (edge[k].vtx1==j2 && edge[k].vtx2==j1) ){
          if( edge[k].c1==-1){ edge[k].c1=i; } // cell[i].ne[j] = k; }
          else if( edge[k].c2==-1){ edge[k].c2=i; } // cell[i].ne[j]=k; }
          break;
        }
      }
    }
  }

  for(size_t i=0;i<CELL_NUMBER;++i){
    Set_NeighbourCells_RandomOrder( i, cell[i] );
  }

}


void Set_NeighbourCells_RandomOrder( int cid, CELL &e_cell )
{
  int jj = 0;  
  for( int i=0; i<NMAX; ++i)e_cell.nc[i] = -1;
  for(int e=0; e<E_NUM;++e){
    if( edge[e].c1 == cid ){ e_cell.nc[jj] = edge[e].c2; ++jj; }
    else if( edge[e].c2 == cid){ e_cell.nc[jj] = edge[e].c1; ++jj; }
  }
  if(jj != e_cell.jnum)ERROR("Mismatch in Set_NeigbourCells_RandomOrder\n");
}


//void Set_BXY(double x[],double y[])
void Set_BXY(vector<double> &X)
{
  double *x= X.data(), *y = X.data()+V_NUM;
  double dx[E_NUM],dy[E_NUM];  
  
  for( size_t i = 0; i < V_NUM; ++i ){
    x[i] -= LX * std::floor(x[i] / LX);
    y[i] -= LY * std::floor(y[i] / LY);
  }

  for( size_t i=0;i<E_NUM;++i ){
    dx[i]=(*edge[i].x1-*edge[i].x2)/LX;
    dy[i]=(*edge[i].y1-*edge[i].y2)/LY;
  }
  
  for(size_t i=0;i<E_NUM;++i){
    edge[i].bx = ( dx[i] > 0.5 ) - ( dx[i] < -0.5 );
    edge[i].by = ( dy[i] > 0.5 ) - ( dy[i] < -0.5 );
  }
  
  for( size_t i=0; i<CELL_NUMBER;++i ){
    cell[i].bx[0] = 0;
    cell[i].by[0] = 0; 
    for(size_t j=1; j<cell[i].jnum; ++j ){
      double dx_ij = (x[cell[i].vtx[j]] - x[cell[i].vtx[0]])/LX;
      double dy_ij = (y[cell[i].vtx[j]] - y[cell[i].vtx[0]])/LY;
      cell[i].bx[j] = -( dx_ij > 0.5 ) + (dx_ij < -0.5 );
      cell[i].by[j] = -( dy_ij > 0.5 ) + (dy_ij < -0.5 );          
    }
  }

  // for check
  double asum = 0.0;
  for( size_t i=0; i<CELL_NUMBER;++i ){ asum += cell_area(cell[i],x,y); }

  if( fabs(asum-LX*LY) > 1.0e-10 ){
    fprintf(stderr,"%e\n",asum-LX*LY);
    fprintf(stderr,"Set_BXY error\n");
    exit(10);
  }

}

bool Compare_Replace(vector<double> &X, const int ex)
{
  vector<VERTEX> t_vtx = vtx;
  vector<EDGE> t_edge = edge;
  vector<CELL> t_cell = cell; 
  double energy,t_energy;

  Replace(ex,X,t_vtx,t_edge,t_cell);

  energy=all_energy(X,edge,cell);
  t_energy=all_energy(X,t_edge,t_cell);

  if(t_energy<energy){    
    vtx = t_vtx;
    edge = t_edge;
    cell = t_cell;
    for(size_t i=0;i<CELL_NUMBER;++i){
      Set_NeighbourCells_RandomOrder( i, cell[i] );
    }
    return true;
  }
  return false;

}


int Replace(const int ex, vector<double> &X, vector<VERTEX> &e_vtx, vector<EDGE> &e_edge, vector<CELL> &e_cell)
{
  double *x=X.data(),*y=X.data()+V_NUM;
  int i,j;
  int p1,p2,p3,p4,tp;
  int e1=0,e2=0,e3=0,e4=0;
  int et1,et2;
  int nc1,nc2,nc3,nc4;

  et1 = e_edge[ex].vtx1;
  et2 = e_edge[ex].vtx2;

  if(e_vtx[et1].c1!=e_vtx[et2].c1 && e_vtx[et1].c1!=e_vtx[et2].c2 && e_vtx[et1].c1!=e_vtx[et2].c3){
    nc3 = e_vtx[et1].c1;
    nc1 = e_vtx[et1].c2;
    nc2 = e_vtx[et1].c3;
  }else if(e_vtx[et1].c2!=e_vtx[et2].c1 && e_vtx[et1].c2!=e_vtx[et2].c2 && e_vtx[et1].c2!=e_vtx[et2].c3){
    nc3 = e_vtx[et1].c2;
    nc1 = e_vtx[et1].c1;
    nc2 = e_vtx[et1].c3;
  }else if(e_vtx[et1].c3!=e_vtx[et2].c1 && e_vtx[et1].c3!=e_vtx[et2].c2 && e_vtx[et1].c3!=e_vtx[et2].c3){
    nc3 = e_vtx[et1].c3;
    nc1 = e_vtx[et1].c1;
    nc2 = e_vtx[et1].c2;
  }else{
    return 0;
  }

  if(cell[nc1].jnum==3 || cell[nc2].jnum==3)return 0;

  if(e_vtx[et2].c1!=nc1 && e_vtx[et2].c1!=nc2){
    nc4 = e_vtx[et2].c1;
  }else if(e_vtx[et2].c2!=nc1 && e_vtx[et2].c2!=nc2){
    nc4 = e_vtx[et2].c2;
  }else if(e_vtx[et2].c3!=nc1  && e_vtx[et2].c3!=nc2){
    nc4 = e_vtx[et2].c3;
  }

  for(i=0;i<e_cell[nc3].jnum;++i){
    if(e_cell[nc3].vtx[i]==et1){
      p1 = ( i!=e_cell[nc3].jnum-1 ? e_cell[nc3].vtx[i+1] : e_cell[nc3].vtx[0] );
      p2 = ( i!=0 ? e_cell[nc3].vtx[i-1] : e_cell[nc3].vtx[e_cell[nc3].jnum-1] );
      break;
    }
  }

  for(i=0;i<e_cell[nc4].jnum;++i){
    if(e_cell[nc4].vtx[i]==et2){
      p4 = ( i!=e_cell[nc4].jnum-1 ? e_cell[nc4].vtx[i+1] : e_cell[nc4].vtx[0] );
      p3 = ( i!=0 ? e_cell[nc4].vtx[i-1] : e_cell[nc4].vtx[e_cell[nc4].jnum-1] );
      break;
    }
  }

  for(i=0;e_cell[nc1].vtx[i]!=et1;++i);
  tp = ( i!=e_cell[nc1].jnum-1 ? e_cell[nc1].vtx[i+1] : e_cell[nc1].vtx[0] );
  if(tp!=et2){	SWAP(int,nc1,nc2);  }

  if( e_edge[e_vtx[et1].e1].vtx1==p1 || e_edge[e_vtx[et1].e1].vtx2==p1 )e1=e_vtx[et1].e1;
  else if(e_edge[e_vtx[et1].e2].vtx1==p1 || e_edge[e_vtx[et1].e2].vtx2==p1 )e1=e_vtx[et1].e2;
  else if(e_edge[e_vtx[et1].e3].vtx1==p1 || e_edge[e_vtx[et1].e3].vtx2==p1 )e1=e_vtx[et1].e3;

  if( e_edge[e_vtx[et1].e1].vtx1==p2 || e_edge[e_vtx[et1].e1].vtx2==p2)e2=e_vtx[et1].e1;
  else if( e_edge[e_vtx[et1].e2].vtx1==p2 || e_edge[e_vtx[et1].e2].vtx2==p2 )e2=e_vtx[et1].e2;
  else if( e_edge[e_vtx[et1].e3].vtx1==p2 || e_edge[e_vtx[et1].e3].vtx2==p2 )e2=e_vtx[et1].e3;

  if(e_edge[e_vtx[et2].e1].vtx1==p3 || e_edge[e_vtx[et2].e1].vtx2==p3 )e3=e_vtx[et2].e1;
  else if(e_edge[e_vtx[et2].e2].vtx1==p3 || e_edge[e_vtx[et2].e2].vtx2==p3 )e3=e_vtx[et2].e2;
  else if(e_edge[e_vtx[et2].e3].vtx1==p3 || e_edge[e_vtx[et2].e3].vtx2==p3 )e3=e_vtx[et2].e3;

  if(e_edge[e_vtx[et2].e1].vtx1==p4 || e_edge[e_vtx[et2].e1].vtx2==p4 )e4=e_vtx[et2].e1;
  else if(e_edge[e_vtx[et2].e2].vtx1==p4 || e_edge[e_vtx[et2].e2].vtx2==p4 )e4=e_vtx[et2].e2;
  else if(e_edge[e_vtx[et2].e3].vtx1==p4 || e_edge[e_vtx[et2].e3].vtx2==p4 )e4=e_vtx[et2].e3;

  /*** not yet changed so far:  ***/
  {
    int ne2bx,ne2by,ne3bx,ne3by;

    /** e_edge bx,by**/
    if(e_edge[e2].vtx1==et1){
      ne2bx=-e_edge[ex].bx+e_edge[e2].bx;
      ne2by=-e_edge[ex].by+e_edge[e2].by;
    }else if(e_edge[e2].vtx2==et1){
      ne2bx=-e_edge[ex].bx-e_edge[e2].bx;
      ne2by=-e_edge[ex].by-e_edge[e2].by;
    }else{
      fprintf(stderr," Wrong -- e2 assigned error : in head.h, Replace\n ");
      exit(10);
    }

    if(e_edge[e3].vtx1==et2){
      ne3bx=e_edge[ex].bx+e_edge[e3].bx;
      ne3by=e_edge[ex].by+e_edge[e3].by;
    }else if(e_edge[e3].vtx2==et2){
      ne3bx=e_edge[ex].bx-e_edge[e3].bx;
      ne3by=e_edge[ex].by-e_edge[e3].by;
    }else{
      fprintf(stderr,"Wrong -- e3 assigned error : in head.h, Replace\n");
      exit(10);
    }
    /* 	if(fabs(ne2bx)>1 || fabs(ne2by)>1 || fabs(ne3bx)>1 || fabs(ne3by)>1)return 0; */

    e_edge[e2].bx=(char)ne2bx;
    e_edge[e2].by=(char)ne2by;
    e_edge[e3].bx=(char)ne3bx;
    e_edge[e3].by=(char)ne3by;

    /*** e_edge ***/
    e_edge[e2].vtx1=et2;
    e_edge[e2].vtx2=p2;
    e_edge[e2].x1=x+et2;
    e_edge[e2].y1=y+et2;
    e_edge[e2].x2=x+p2;
    e_edge[e2].y2=y+p2;
    e_edge[e3].vtx1=et1;
    e_edge[e3].vtx2=p3;
    e_edge[e3].x1=x+et1;
    e_edge[e3].y1=y+et1;
    e_edge[e3].x2=x+p3;
    e_edge[e3].y2=y+p3;
    e_edge[ex].c1 = nc3;
    e_edge[ex].c2 = nc4;

    /*** e_cell ***/
    if(e_cell[nc3].jnum==NMAX-1 || e_cell[nc4].jnum==NMAX-1){
      output(X,stdout);
      fprintf(stderr,"Too many cell.jnum --  in head.h, Replace\n");
      exit(10);
    }

    j=0;
    for(i=0;i<e_cell[nc1].jnum;++i){
      if(e_cell[nc1].vtx[i]==et2)++j;
      e_cell[nc1].vtx[i]=e_cell[nc1].vtx[j];
      e_cell[nc1].bx[i]=e_cell[nc1].bx[j];
      e_cell[nc1].by[i]=e_cell[nc1].by[j];
      ++j;
    }
    e_cell[nc1].jnum--;

    j=0;
    for(i=0;i<e_cell[nc2].jnum;++i){
      if(e_cell[nc2].vtx[i]==et1)++j;
      e_cell[nc2].vtx[i]=e_cell[nc2].vtx[j];
      e_cell[nc2].bx[i]=e_cell[nc2].bx[j];
      e_cell[nc2].by[i]=e_cell[nc2].by[j];
      ++j;
    }
    e_cell[nc2].jnum--;

    j=e_cell[nc3].jnum-1;
    for(i=e_cell[nc3].jnum;i>=0;i--){
      if(i!=e_cell[nc3].jnum && e_cell[nc3].vtx[i]==et1){
        e_cell[nc3].vtx[i]=et2;
        e_cell[nc3].bx[i]=e_cell[nc3].bx[i]+e_edge[ex].bx;
        e_cell[nc3].by[i]=e_cell[nc3].by[i]+e_edge[ex].by;
        i--;
      }
      e_cell[nc3].vtx[i]=e_cell[nc3].vtx[j];
      e_cell[nc3].bx[i]=e_cell[nc3].bx[j];
      e_cell[nc3].by[i]=e_cell[nc3].by[j];
      j--;
    }
    e_cell[nc3].jnum++;

    j=e_cell[nc4].jnum-1;
    for(i=e_cell[nc4].jnum;i>=0;i--){
      if(i!=e_cell[nc4].jnum && e_cell[nc4].vtx[i]==et2){
        e_cell[nc4].vtx[i]=et1;
        e_cell[nc4].bx[i]=e_cell[nc4].bx[i]-e_edge[ex].bx;
        e_cell[nc4].by[i]=e_cell[nc4].by[i]-e_edge[ex].by;
        i--;
      }
      e_cell[nc4].vtx[i]=e_cell[nc4].vtx[j];
      e_cell[nc4].bx[i]=e_cell[nc4].bx[j];
      e_cell[nc4].by[i]=e_cell[nc4].by[j];
      j--;
    }
    e_cell[nc4].jnum++;

    /**** vtx beighbour ****/
    e_vtx[et1].nv1=p1;
    e_vtx[et1].nv2=et2;
    e_vtx[et1].nv3=p3;
    e_vtx[et1].e1=e1;
    e_vtx[et1].e2=ex;
    e_vtx[et1].e3=e3;
    e_vtx[et1].c1=nc4;
    e_vtx[et1].c2=nc1;
    e_vtx[et1].c3=nc3;

    e_vtx[et2].nv1=p4;
    e_vtx[et2].nv2=et1;
    e_vtx[et2].nv3=p2;
    e_vtx[et2].e1=e4;
    e_vtx[et2].e2=ex;
    e_vtx[et2].e3=e2;
    e_vtx[et2].c1=nc3;
    e_vtx[et2].c2=nc2;
    e_vtx[et2].c3=nc4;

    if(e_vtx[p2].nv1==et1){e_vtx[p2].nv1=et2;}
    else if(e_vtx[p2].nv2==et1){e_vtx[p2].nv2=et2;}
    else if(e_vtx[p2].nv3==et1){e_vtx[p2].nv3=et2;}

    if(e_vtx[p3].nv1==et2){e_vtx[p3].nv1=et1;}
    else if(e_vtx[p3].nv2==et2){e_vtx[p3].nv2=et1;}
    else if(e_vtx[p3].nv3==et2){e_vtx[p3].nv3=et1;}

    //Set_NeigbourCells_RandomOrder( nc1, e_cell[nc1] );
    //Set_NeigbourCells_RandomOrder( nc2, e_cell[nc2] );
    //Set_NeigbourCells_RandomOrder( nc3, e_cell[nc3] );
    //Set_NeigbourCells_RandomOrder( nc4, e_cell[nc4] );
  }

  {
    for(i=0;i<e_cell[nc1].jnum;++i){
      if(fabs(e_cell[nc1].bx[i])>1 ||fabs(e_cell[nc1].by[i])>1 )return 2;
    }
    for(i=0;i<e_cell[nc2].jnum;++i){
      if(fabs(e_cell[nc2].bx[i])>1 ||fabs(e_cell[nc2].by[i])>1 )return 2;
    }
    for(i=0;i<e_cell[nc3].jnum;++i){
      if(fabs(e_cell[nc3].bx[i])>1 ||fabs(e_cell[nc3].by[i])>1 )return 2;
    }
    for(i=0;i<e_cell[nc4].jnum;++i){
      if(fabs(e_cell[nc4].bx[i])>1 ||fabs(e_cell[nc4].by[i])>1 )return 2;
    }
    return 1;
  }

}
