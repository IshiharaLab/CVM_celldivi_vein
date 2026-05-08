
/***********************************/
void load_initialconfiguraton(char *inputfile)
{  
  FILE *file;
  char *csl,line[2000];
  int cn=0,en=0,jn=0;
  int tbx,tby;
  bool with_celltype = false; // cell type information flag

  file = fopen(inputfile,"r");
  while((csl=fgets(line,sizeof(line),file))!=NULL){
    if(strcmp(line,"\n")==0){
      ;
    }else if(strstr(line,"#")!=NULL){
      //fprintf(stderr,"%s",line);
      if(strstr(line,"V_NUM")!=NULL){
        csl=strtok(line,"=");  csl=strtok(NULL,"=");
        sscanf(csl,"%d",&jn);
        V_NUM = jn;
        fprintf(stderr,"# CELL_NUMBER=%d ,,, OK\n",CELL_NUMBER);
        jn=0;
      }else if(strstr(line,"E_NUM")!=NULL){
        csl=strtok(line,"=");  csl=strtok(NULL,"=");
        sscanf(csl,"%d",&en);
        E_NUM = en;
        en=0;
      }else if(strstr(line,"C_NUM")!=NULL){
        csl=strtok(line,"=");  csl=strtok(NULL,"=");
        sscanf(csl,"%d",&cn);
        CELL_NUMBER = cn;
        cn=0;
      }else if(strstr(line,"LX")!=NULL){
        csl=strtok(line," ");  csl=strtok(NULL," ");  csl=strtok(NULL," ");
        sscanf(csl,"%lf",&LX);
      }else if(strstr(line,"LY")!=NULL){
        csl=strtok(line," ");  csl=strtok(NULL," ");  csl=strtok(NULL," ");
        sscanf(csl,"%lf",&LY);
      }else if(strstr(line,"CELL_TYPE")!=NULL){
        with_celltype = true;
        printf("%s",line);      
      }else if (strstr(line,"------")!=NULL){
        // check
        if( CELL_NUMBER == 0){ fprintf(stderr,"!!! CELL_NUMBER = 0\n");  exit(1); }
        if( V_NUM != 2*CELL_NUMBER){ fprintf(stderr,"!!! V_NUM != 2*CELL_NUMBER\n");  exit(2); }
        if( E_NUM != 3*CELL_NUMBER){ fprintf(stderr,"!!! E_NUM != 3*CELL_NUMBER\n");  exit(2); }

        fprintf(stderr,"CELL_NUMBER = %d\n",CELL_NUMBER);
        fprintf(stderr,"V_NUM = %d\n",V_NUM);
        fprintf(stderr,"E_NUM = %d\n",E_NUM);
        fprintf(stderr,"LX = %f, LY = %f\n",LX,LY);
        fprintf(stderr,"\n");

        Assign_Variables();  
      }
    }else{
      if( strstr(line,"J[")!=NULL ){
        while(*csl != ' ')csl++;
        sscanf(csl,"%lf%lf",x+jn,y+jn);
        jn++;
      }else if(strstr(line,"E[")!=NULL){
        while(*csl != ' ')csl++;
        sscanf(csl,"%d%d%d%d",&edge[en].vtx1,&edge[en].vtx2,&tbx,&tby);
        edge[en].set_ptrs((char)tbx,(char)tby);
        en++;
      }else if(strstr(line,"C[")!=NULL){      
        while(*csl != ' ')csl++;
        sscanf(csl,"%d",&cell[cn].jnum);
        while(*csl != ':')csl++;
        csl+=2;
        for(int i=0;i<cell[cn].jnum;++i){
          sscanf(csl,"%d",&cell[cn].vtx[i]);
          while(*csl != ' ')csl++;
          csl++;
        }
        while(*csl != ':')csl++;
        csl+=2;
        for(int i=0;i<cell[cn].jnum;++i){
          sscanf(csl,"%d",&tbx);
          cell[cn].bx[i]=(char)tbx;
          while(*csl != ' ')csl++;
          csl++;
        }
        while(*csl != ':')csl++;
        csl+=2;
        for(int i=0;i<cell[cn].jnum;++i){
          sscanf(csl,"%d",&tby);
          cell[cn].by[i]=(char)tby;
          while(*csl != ' ')csl++;
          csl++;
        }        
        if( with_celltype ){ // if there is cell type information ; with_celltype = true          
          while(*csl != ':')csl++;
          csl+=2;        
          sscanf(csl,"%c",&cell[cn].cell_type); // read cell type          
        }else{
          cell[cn].cell_type = 'i'; // default intervein
        }               
        cn++;
      }
    }
  }
  fclose(file);
  //exit(2);

  Set_Vertex(X);  /* put at last */  

}



void output(vector<double> &X,FILE *out_file)
{
  Set_BXY(X); 
  double *x=X.data(),*y=X.data()+V_NUM;

  fprintf(out_file,"# # FOR PERIODIC BOUNDARY CONDITION\n");
  fprintf(out_file,"# V_NUM= %d \n",V_NUM);
  fprintf(out_file,"# E_NUM= %d \n",E_NUM);
  fprintf(out_file,"# C_NUM= %d \n",CELL_NUMBER);
  fprintf(out_file,"# LX= %e \n",LX);
  fprintf(out_file,"# LY= %e \n",LY);
  fprintf(out_file,"# CELL_TYPE= 1 \n"); // cell_type あり、
  fprintf(out_file,"# ----------------------- \n");
  for(int i=0;i<V_NUM;++i){
    fprintf(out_file,"J[%d] %e %e\n",i,x[i],y[i]);
  }
  fprintf(out_file,"\n");

  for(int i=0;i<E_NUM;++i){
    fprintf(out_file,"E[%d] %d %d %d %d\n",i,edge[i].vtx1,edge[i].vtx2,edge[i].bx,edge[i].by);
  }

  fprintf(out_file,"\n");
  for(int i=0;i<CELL_NUMBER;++i){
    fprintf(out_file,"C[%d] %d : ",i,cell[i].jnum);
    for(size_t j=0;j<cell[i].jnum;++j)fprintf(out_file,"%d ",cell[i].vtx[j]);
    fprintf(out_file," : ");
    for(size_t j=0;j<cell[i].jnum;++j)fprintf(out_file,"%d ",cell[i].bx[j]);
    fprintf(out_file," : ");
    for(size_t j=0;j<cell[i].jnum;++j)fprintf(out_file,"%d ",cell[i].by[j]);
    fprintf(out_file," : %c ",cell[i].cell_type); // this line is added to output cell type
    fprintf(out_file,"\n");
  }
  fflush(out_file);

}



void OutputTensionPressure(const vector<double> &X,FILE *out_file)
{
  const double *x=X.data(),*y=X.data()+V_NUM;
  double dx[E_NUM],dy[E_NUM],dist[E_NUM],tension[E_NUM];
  double area[CELL_NUMBER],press[CELL_NUMBER],length[CELL_NUMBER];
  double tenav=0.0,preav=0.0;

  for(size_t i=0;i<E_NUM;++i)dx[i]=(*edge[i].x1-*edge[i].x2-edge[i].bx*LX);
  for(size_t i=0;i<E_NUM;++i)dy[i]=(*edge[i].y1-*edge[i].y2-edge[i].by*LY);
  for(size_t i=0;i<E_NUM;++i)dist[i]=sqrt(dx[i]*dx[i]+dy[i]*dy[i]);
  for(size_t i=0;i<CELL_NUMBER;++i)area[i]=cell_area(cell[i],x,y);
  for(size_t i=0;i<CELL_NUMBER;++i)length[i]=boundary_length( cell[i], x, y );

  for(size_t i=0;i<CELL_NUMBER;++i){
    press[i]=-KV*(area[i]-Vo[i]);
    preav += press[i];
  }
  preav/=(double)CELL_NUMBER;
  for(size_t i=0;i<E_NUM;++i){
    tension[i]=KL*(length[edge[i].c1]-Lo[edge[i].c1]+length[edge[i].c2]-Lo[edge[i].c2]);
    tenav+=tension[i];
  }
  tenav/=E_NUM;
  for(int i=0;i<E_NUM;++i)fprintf(out_file,"%d %f %f %f - - - -\n",i,dist[i],tension[i],tension[i]/tenav );
  fprintf(out_file,"\n");
  for(int i=0;i<CELL_NUMBER;++i)fprintf(out_file,"- - - - %d %f %f %f\n",i,area[i],press[i],(press[i]-preav)/tenav);
}


void out_geometry(const vector<double> &X,FILE *out_file)
{
  const double *x=X.data(),*y=X.data()+V_NUM;

  for(size_t i=0;i<E_NUM;++i){
    if( edge[i].bx==0 && edge[i].by==0 ){
      fprintf(out_file,"%f %f\n%f %f\n\n",x[edge[i].vtx1],y[edge[i].vtx1],x[edge[i].vtx2],y[edge[i].vtx2]);
    }else if(edge[i].bx!=0 &&edge[i].by==0){
      fprintf(out_file,"%f %f\n%f %f\n\n",x[edge[i].vtx1],y[edge[i].vtx1],x[edge[i].vtx2]+LX*edge[i].bx,y[edge[i].vtx2]);
      fprintf(out_file,"%f %f\n%f %f\n\n",x[edge[i].vtx1]-LX*edge[i].bx,y[edge[i].vtx1],x[edge[i].vtx2],y[edge[i].vtx2]);
    }else if(edge[i].bx==0 &&edge[i].by!=0){
      fprintf(out_file,"%f %f\n%f %f\n\n",x[edge[i].vtx1],y[edge[i].vtx1],x[edge[i].vtx2],y[edge[i].vtx2]+LY*edge[i].by);
      fprintf(out_file,"%f %f\n%f %f\n\n",x[edge[i].vtx1],y[edge[i].vtx1]-LY*edge[i].by,x[edge[i].vtx2],y[edge[i].vtx2]);
    }else{
      fprintf(out_file,"%f %f\n%f %f\n\n",x[edge[i].vtx1],y[edge[i].vtx1],x[edge[i].vtx2]+LX*edge[i].bx,y[edge[i].vtx2]+LY*edge[i].by);
      fprintf(out_file,"%f %f\n%f %f\n\n",x[edge[i].vtx1]-LX*edge[i].bx,y[edge[i].vtx1]-LY*edge[i].by,x[edge[i].vtx2],y[edge[i].vtx2]);
    }
  }
}

/***********************************/
