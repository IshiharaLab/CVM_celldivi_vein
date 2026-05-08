struct DIV_CELL_DATA{
  int dc;
  int e1,e2;
  double mx1,my1,mx2,my2;
};
typedef struct DIV_CELL_DATA DIV_CELL_DATA;

int PRE_CELL_DIVISION( vector<double> &X, vector<double> &Vo, int divc);
void CELL_DIVISION( vector<double> &X,vector<double> &Vo,DIV_CELL_DATA dcst ); // called in PRE_CELL_DIVISION

void cross_point(double &crsx,double &crsy, double A0,double B0,double X0,double Y0, double A1,double B1,double X1,double Y1)
{
  const double a0 = A0*X0-B0*Y0;
  const double a1 = A1*X1-B1*Y1;
  const double det = -A0*B1+A1*B0;
  crsx = (-B1*a0+B0*a1)/det;
  crsy = (-A1*a0+A0*a1)/det;
}

struct cell_div_workingspace{
  int ja1,ja2,jb1,jb2; 
} cdwk;

/** Crossing point of "tedge" and A(x-X0)-B(y-Y0)=0 **/
DIV_CELL_DATA Div_Centerd_Angle( int dc, const CELL &dcell, vector<double> &X, double angle )
{
  double *x = X.data(), *y = X.data()+V_NUM;
  DIV_CELL_DATA tdcst;
  int cf;
  double dvcx[NMAX],dvcy[NMAX];
  double mcx = 0.0, mcy =0.0;
  double tmx,tmy;
  double dx1,dy1;

  tdcst.dc = dc; // dividing cell
  for(size_t j=0;j<dcell.jnum;++j){
    const double x1 = x[dcell.vtx[j]]+dcell.bx[j]*LX;    
    const double y1 = y[dcell.vtx[j]]+dcell.by[j]*LY;    
    dvcx[j] = x1;
    dvcy[j] = y1;
    mcx += x1;
    mcy += y1;
  }
  mcx /= dcell.jnum;
  mcy /= dcell.jnum;

  //fprintf(stderr,"# div_at= %f %f angle= %f \n",mcx,mcy,angle); 
  cf=0;
  for( int i=0; i<dcell.jnum; i++ ){
  	int ii = ( i!=dcell.jnum-1 ? i+1 : 0 );
    dx1 = dvcx[ii] - dvcx[i];
    dy1 = dvcy[ii] - dvcy[i];   

	  cross_point( tmx, tmy, sin(angle), cos(angle), mcx, mcy, dy1, dx1, dvcx[i], dvcy[i] );

	  if( (dvcx[i]<tmx && dvcx[ii]>tmx) || (dvcx[i]>tmx && dvcx[ii]<tmx) ){
  	  if(cf==0){
	      cdwk.ja1 = dcell.vtx[i];
    		cdwk.ja2 = dcell.vtx[ii];
  	    tdcst.mx1 = tmx;
      	tdcst.my1 = tmy;
	      cf++;
    	}else{
		    cdwk.jb1 = dcell.vtx[i];
        cdwk.jb2 = dcell.vtx[ii];
      	tdcst.mx2 = tmx;
    	  tdcst.my2 = tmy;
        cf++;
	    }
	  }
  }
  if(cf!=2){
    ERROR(" Error at cf in head_cell_div: abnormal edge-crossing error by cell division can casuse this error. Try again. \n" );
  }

  cf=0;
  for( int i=0; cf<2; i++ ){
    if( edge[i].vtx1==cdwk.ja1 && edge[i].vtx2==cdwk.ja2){ tdcst.e1=i; ++cf;}
    if( edge[i].vtx1==cdwk.ja2 && edge[i].vtx2==cdwk.ja1){ tdcst.e1=i; ++cf;}
    if( edge[i].vtx1==cdwk.jb1 && edge[i].vtx2==cdwk.jb2){ tdcst.e2=i; ++cf;}
    if( edge[i].vtx1==cdwk.jb2 && edge[i].vtx2==cdwk.jb1){ tdcst.e2=i; ++cf;}
  }  
  return tdcst;
}




int PRE_CELL_DIVISION( vector<double> &X, vector<double> &Vo, int divc, const double pt)
{
  //fprintf(stderr,"# cell %d is dividing at time %f\n", divc, pt );

  /* determine the cell divsion points */
  DIV_CELL_DATA dcst = Div_Centerd_Angle( divc, cell[divc], X, 2*M_PI*unidist(mt) );
  
  /** Division **/
  CELL_DIVISION( X, Vo, dcst);  // CELL_NUMBER is updated inside 

  /* Parameter changes after cell division*/
  {
    int c1 = divc, c2 = CELL_NUMBER-1;
    Vo.push_back(Vo[c1]); // 
    Lo.push_back(Lo[c1]);
    Vo[c1] = BVo/2.0; 
    Vo[c2] = Vo[divc]; 
    Lo[c1] = ShpIdx*sqrt(Vo[c1]); 
    //Lo[c1] = BLo/sqrt(2.0);
    Lo[c2] = Lo[c1]; 
    
    cell[c1].cell_type = 'd'; 
    cell[c2].cell_type = 'd';
  	cell[c1].cell_div_count +=1;     /* cell division occurs only once in this model  */
  	cell[c2].cell_div_count += 1;
    cell[c1].cell_div_time = pt;
    cell[c2].cell_div_time = pt;  
  }  

  return divc;    
}


void CELL_DIVISION( vector<double> &X, vector<double> &Vo, DIV_CELL_DATA dcst )
{ 
  int nc1,nc2; // neighbor cells of div,   
  vector<double> py(V_NUM); 
  // const int ja1 = cdwk.ja1; // unused
  const int ja2 = cdwk.ja2; 
  // const int jb1 = cdwk.jb1; // unsed
  const int jb2 = cdwk.jb2;

  for( int i=0; i<V_NUM; i++ )py[i]=y[i];

  CELL_NUMBER++; // +1
  V_NUM=2*CELL_NUMBER; // +2 
  E_NUM=3*CELL_NUMBER; // +3
  
  X.insert( X.end(), 4, 0.0 ); //  +4
  x = X.data(); 
  y = X.data()+V_NUM;
  
  for(int i=0; i<V_NUM-2; ++i)y[i] = py[i];
  const int jn1 = V_NUM-2;
  const int jn2 = V_NUM-1;
  x[jn1] = dcst.mx1;
  y[jn1] = dcst.my1;
  x[jn2] = dcst.mx2;
  y[jn2] = dcst.my2;
  
  vtx.resize(V_NUM);   // +2
  edge.resize(E_NUM);  // +3
  cell.resize(CELL_NUMBER); // +1
 
  /** vertex **/
  for(size_t i=0;i<V_NUM;++i){
    vtx[i].x = x + i;
    vtx[i].y = y + i;
  }
  
  /** edge **/
  const int en1 = E_NUM-3;  
  const int en2 = E_NUM-2;  
  const int en3 = E_NUM-1;
  edge[en3].vtx1=jn1;  edge[en3].vtx2=jn2;
  edge[en1].vtx1=jn1;  edge[en1].vtx2=edge[dcst.e1].vtx2;  edge[dcst.e1].vtx2=jn1;
  edge[en2].vtx1=jn2;  edge[en2].vtx2=edge[dcst.e2].vtx2;  edge[dcst.e2].vtx2=jn2;

  for(size_t i=0;i<E_NUM;++i){ // All edges are reassubged (bacause X is changed)
    edge[i].x1 = x + edge[i].vtx1;
    edge[i].y1 = y + edge[i].vtx1;
    edge[i].x2 = x + edge[i].vtx2;
    edge[i].y2 = y + edge[i].vtx2;
  }
      
  double dx[E_NUM],dy[E_NUM]; 
  for(size_t i=0;i<E_NUM;++i){
    dx[i]=(*edge[i].x1-*edge[i].x2)/LX;
    dy[i]=(*edge[i].y1-*edge[i].y2)/LY;
    edge[i].bx = ( dx[i] > 0.5 ) - (dx[i] < -0.5 );
    edge[i].by = ( dy[i] > 0.5 ) - (dy[i] < -0.5 );
  }  
  
  // set neibouring cell nc1, nc2
  if( edge[dcst.e1].c1 == dcst.dc ) nc1 = edge[dcst.e1].c2;
  else if( edge[dcst.e1].c2 == dcst.dc ) nc1 = edge[dcst.e1].c1;
  else ERROR("Some Errors!?");

  //printf("%d %d\n",edge[dcst.e2].c1, edge[dcst.e2].c2);
  if(edge[dcst.e2].c1 == dcst.dc) nc2 = edge[dcst.e2].c2;
  else if(edge[dcst.e2].c2 == dcst.dc) nc2 = edge[dcst.e2].c1;
  else ERROR("Some Errors!?");
  
  /** for check
  CELL *cp = &cell[dcst.dc];
  for(int i=0; i<cp->jnum; ++i){  
    printf("%f %f\n",x[cp->vtx[i]],y[cp->vtx[i]] );
  }
  printf("%f %f\n",x[cp->vtx[0]],y[cp->vtx[0]] );
  printf("\n"); **/

  /** cell **/
  const int nc = CELL_NUMBER-1;
  {
    CELL tcell[2],scell;
    int cc;

    tcell[0] = cell[dcst.dc];
    tcell[1] = cell[dcst.dc];
	    
    tcell[0].jnum=0;
    tcell[1].jnum=0;
    cc=0;

    for(int i=0;i<cell[dcst.dc].jnum;i++){
      if(cell[dcst.dc].vtx[i]==ja2){
        tcell[cc].vtx[tcell[cc].jnum]=jn1;
        tcell[cc].jnum++;
        cc=1-cc;
        tcell[cc].vtx[tcell[cc].jnum]=jn1;
        tcell[cc].jnum++;
      }
      if(cell[dcst.dc].vtx[i]==jb2){
        tcell[cc].vtx[tcell[cc].jnum]=jn2;
        tcell[cc].jnum++;
        cc=1-cc;
        tcell[cc].vtx[tcell[cc].jnum]=jn2;
        tcell[cc].jnum++;
      }
      tcell[cc].vtx[tcell[cc].jnum] = cell[dcst.dc].vtx[i];
      tcell[cc].jnum++;
    }

    /*** change neighbor cells sc1, sc2***/
    // nc1 
    scell = cell[nc1];
    scell.jnum=0;
    for(size_t i=0;i<cell[nc1].jnum;i++){
      scell.vtx[scell.jnum]=cell[nc1].vtx[i];
      scell.jnum++;
      if(cell[nc1].vtx[i]==ja2){
        scell.vtx[scell.jnum]=jn1;
        scell.jnum++;        
      }      
    }
    cell[nc1]=scell;    
  
    // nc2
    scell = cell[nc2];
    scell.jnum=0;
    for(size_t i=0;i<cell[nc2].jnum;i++){
      scell.vtx[scell.jnum]=cell[nc2].vtx[i];
      scell.jnum++;
      if(cell[nc2].vtx[i]==jb2){
        scell.vtx[scell.jnum]=jn2;
        scell.jnum++;
      }
    }       
    cell[nc2]=scell;

    cell[dcst.dc] = tcell[0];
    cell[nc] = tcell[1];
    
  }

  /** for check
  cp = &cell[dcst.dc];
  for(int i=0; i<cp->jnum; ++i){  
    printf("- - %f %f\n",x[cp->vtx[i]],y[cp->vtx[i]] );
  }
  printf("- - %f %f\n",x[cp->vtx[0]],y[cp->vtx[0]] );
  printf("\n");

  cp = &cell[nc];
  for(int i=0; i<cp->jnum; ++i){  
    printf("- - %f %f\n",x[cp->vtx[i]],y[cp->vtx[i]] );
  }
  printf("- - %f %f\n",x[cp->vtx[0]],y[cp->vtx[0]] );
  printf("\n");

  cp = &cell[nc1];
  for(int i=0; i<cp->jnum; ++i){  
    printf("- - - - %f %f\n",x[cp->vtx[i]],y[cp->vtx[i]] );
  }
  printf("- - - - %f %f\n",x[cp->vtx[0]],y[cp->vtx[0]] );
  printf("\n");

  cp = &cell[nc2];
  for(int i=0; i<cp->jnum; ++i){  
    printf("- - - - %f %f\n",x[cp->vtx[i]],y[cp->vtx[i]] );
  }
  printf("- - - - %f %f\n",x[cp->vtx[0]],y[cp->vtx[0]] );
  printf("\n");
  */

  /** Re-Set Topology**/  
  Set_BXY(X);
  Set_Vertex(X);  /** Calculate All **/

  /** Check **/
  for(int i=0;i<V_NUM;i++){
    if((edge[vtx[i].e1].vtx1==i && edge[vtx[i].e1].vtx2==vtx[i].nv1) || (edge[vtx[i].e1].vtx2==i && edge[vtx[i].e1].vtx1==vtx[i].nv1));else ERROR("E ERROR1");
    if((edge[vtx[i].e2].vtx1==i && edge[vtx[i].e2].vtx2==vtx[i].nv2) || (edge[vtx[i].e2].vtx2==i && edge[vtx[i].e2].vtx1==vtx[i].nv2));else ERROR("E ERROR2");
    if((edge[vtx[i].e3].vtx1==i && edge[vtx[i].e3].vtx2==vtx[i].nv3) || (edge[vtx[i].e3].vtx2==i && edge[vtx[i].e3].vtx1==vtx[i].nv3));else ERROR("E ERROR3");
  }    
  
}
