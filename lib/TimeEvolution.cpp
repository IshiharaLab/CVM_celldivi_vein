//defined in ODE.cpp
#define Dt 0.001
const double HDt=0.5*Dt;
const double dt=Dt/6.0;
//void adaptiveCKrungekutta( double X[],double Xerr[],double kX[],double Xout[],double h,double pt );
//void adaptive(double *X,double *h,double *pt,double *dh);
#include "ODE.cpp"
void RK2( vector<double>& X, double& pt, RK2Workspace& rk2wk);




/******************************************
      Energy and time evolution equations
********************************************/
inline double all_energy(const vector<double> &X, const vector<EDGE> &e_edge, const vector<CELL> &e_cell)
{
  const double *x=X.data(),*y=X.data()+V_NUM;
  double v_energy=0.0, r_energy=0.0;
  //double dx[E_NUM],dy[E_NUM],dist[E_NUM]; // ,theta[E_NUM];
  double area[CELL_NUMBER],length[CELL_NUMBER];
  double energy=0.0;

  //for(size_t i=0;i<E_NUM;++i){
  //  dist[i] = hypot( (*e_edge[i].x1-*e_edge[i].x2-e_edge[i].bx*LX), (*e_edge[i].y1-*e_edge[i].y2-e_edge[i].by*LY) );     
  //  //theta[i] = atan2(dy[i],dx[i]);
  //}
  //for(size_t i=0;i<E_NUM;++i)l_energy+=KCr[i]*dist[i];

  for(size_t i=0;i<CELL_NUMBER;++i){
    area[i] = cell_area( e_cell[i],x,y );
    length[i] = boundary_length( e_cell[i], x, y );
    v_energy += (area[i]-Vo[i])*(area[i]-Vo[i]);
    r_energy += (length[i]-Lo[i])*(length[i]-Lo[i]);
  }
  energy = 0.5*KV*v_energy + 0.5*KL*r_energy;

  return energy;
}


void onestep( vector<double> &X, vector<double> &kX, const double pt )
{
  const size_t n = X.size(); // = 2*V_NUM, 拡張を考えてこうしている
  const double* x = X.data();
  const double* y = X.data() + V_NUM;
  kX.assign(n, 0.0);  // サイズを n に揃えてゼロ初期化
  double *kx = kX.data(), *ky = kX.data() + V_NUM; 

  vector<double> clength(CELL_NUMBER,0.0);
  
  // 各細胞の面積と周囲長
  for(size_t i=0;i<CELL_NUMBER;++i){
    clength[i] = boundary_length( cell[i], x, y );
  }

  // 張力由来の圧力項
  for ( size_t i = 0; i < (size_t)E_NUM; ++i ) {
    const int v1 = edge[i].vtx1;
    const int v2 = edge[i].vtx2;
    const int c1 = edge[i].c1;
    const int c2 = edge[i].c2;

    const double dx = (x[v1] - x[v2] - edge[i].bx * LX);
    const double dy = (y[v1] - y[v2] - edge[i].by * LY);
    double L = hypot(dx,dy);
    if (L < 1e-15) L = 1e-15;               // 0割れ防止（極小閾値）
    const double invlen = 1.0 / L;
     
    const double tx = dx * invlen;     // 単位接線ベクトル
    const double ty = dy * invlen;
    const double tension = KL * ( clength[c1]-Lo[c1] + clength[c2]-Lo[c2]);

    kx[v1] -= tension * tx;
    ky[v1] -= tension * ty;
    kx[v2] += tension * tx;
    ky[v2] += tension * ty;
  }

  // 2. 体積（面積）由来の圧力項
  for ( size_t ci = 0; ci < (size_t)CELL_NUMBER; ++ci) {
    const int pn = cell[ci].jnum - 1;
    if (pn < 1) continue;

    int*  __restrict ij  = cell[ci].vtx;
    char* __restrict cbx = cell[ci].bx;
    char* __restrict cby = cell[ci].by;

    const double areaVo = KV * (cell_area(cell[ci], x, y) - Vo[ci]);
  
    kx[ij[0]] -= areaVo * ( y[ij[1]] + cby[1]*LY - (y[ij[pn]] + cby[pn]*LY) );
    ky[ij[0]] -= areaVo * ( x[ij[pn]] + cbx[pn]*LX - (x[ij[1]] + cbx[1]*LX) );  
    for (int j = 1; j < pn; ++j) {
      kx[ij[j]] -= areaVo * ( y[ij[j+1]] + cby[j+1]*LY - (y[ij[j-1]] + cby[j-1]*LY) );
      ky[ij[j]] -= areaVo * ( x[ij[j-1]] + cbx[j-1]*LX - (x[ij[j+1]] + cbx[j+1]*LX) );
    }
    kx[ij[pn]] -= areaVo * ( y[ij[0]] + cby[0]*LY - (y[ij[pn-1]] + cby[pn-1]*LY) );
    ky[ij[pn]] -= areaVo * ( x[ij[pn-1]] + cbx[pn-1]*LX - (x[ij[0]] + cbx[0]*LX) );
  }
  // (void)pt; // ここでは未使用
}


void InitialRelaxation(double RELAX_TIME)
{
  double pt = 0.0;
  for( int tc=0; pt<RELAX_TIME; tc++ ){
    static RK2Workspace rk2wk;    
    double dx[E_NUM],dy[E_NUM];
    for(size_t i=0;i<E_NUM;++i){
      dx[i]=(*edge[i].x1-*edge[i].x2-edge[i].bx*LX);
      dy[i]=(*edge[i].y1-*edge[i].y2-edge[i].by*LY);    
      if( dx[i]*dx[i]+dy[i]*dy[i]<1.0e-8 ){
        if( Compare_Replace(X,i) ){ // 頂点入れ替え(エネルギーが減れば), lib/head_cell_vertex.cpp参照
          fprintf(stderr,"--%f %d (%d-%d)\n",pt,(int)i,edge[i].vtx1,edge[i].vtx2);
        }
      }
    }

    //adaptive(X,&h,&pt,&dh);
    RK2( X, pt, rk2wk ); // Runge-Kutta 2nd order , lib/ODE.cpp参照
    //double NG = 1.0, NS = 0.1;
    //for( size_t i=0; i<E_NUM; ++i ){ KCr[i]=KCr[i]+NG*Dt*(KC-KCr[i])+NS*KC*sqrt(2.0*NG*Dt)*normdist(mt); }
    if( tc%128==0 )Set_BXY(X);
  }
}