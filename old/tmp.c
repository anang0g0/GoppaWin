OP gcd(OP a,OP b){
  OP r={0},s={0};

  s=a;
  
  printpol(o2v(a));
  printf(" a==============\n");
  printpol(o2v(b));
  printf(" b==============\n");
  //  exit(1);
  /* ユークリッドの互除法 */
  r = omod(a , b);
  if(deg(o2v(r))>0){
  printpol(o2v(r));
  printf(" r==============\n");
  }
  // exit(1);
  while(deg(o2v(r))>=0){
    printpol(o2v(r));
    printf(" r2==============\n");

    a = b;
    b = r;
     if(deg(o2v(b))>0)
      r = omod(a ,b);
      //  } 
    printpol(o2v(a));
    printf(" aa==============\n");
    if(deg(o2v(r))==0 && r.t[0].a==1)
      return r;
    if(deg(o2v(r))==0 && r.t[0].a==0)
      return a;
  }
  printpol(o2v(b));
  printf(" be=================\n");
  //  return b;
  //    exit(1);
}




unsigned int ipow(unsigned int q,unsigned int u){
unsigned int i,m=1;

  for(i=0;i<u;i++)
    m*=q;

  printf("in ipow====%d\n",m);

  return m;
}





OP benor(int ww,int nn){
  OP w={0},ff={0},f={0},tt={0};
  int flg=0,i,j,k;
  vec v={0};

  while(1){
    for(i=0;i<DEG;i++)
      v.x[i]=0;
    
    v.x[0]=1;
    v.x[nn]=1;
    k=0;
    
    while(k<ww){
      j=xor128()%10;
      if(j!=0 && j!=10 && v.x[j]==0){
	v.x[j]=1;
	k++;
      }
    }
  printpol(v);
  printf(" v=============\n");
  tt=v2o(v);
  //  exit(1);
  j=0;
  //    tt=irr(1,10);
  for(i=1;i<deg(v)/2+1;i++){
    f.t[0].a=1;
    f.t[0].n=1;
    f.t[1].a=1;
    f.t[1].n=ipow(2,i);
    
    ff=gcd(tt,f);
    if(deg(o2v(ff))>0){
      printpol(o2v(ff));
      printf(" lcm==============\n");
    }
    if(deg(o2v(ff))==0)
      j++;
    printpol(o2v(ff));
    printf(" ff=============\n");
    printf("j=%d\n",j);
    printpol(v);
    printf(" vv===================\n");
  }
  if(j==deg(v)/2){
    printpol(o2v(tt));
    printf(" irr?============\n");
    return tt;
  }
  j=0;
  //exit(1);
    }

  
}



OP keyfinder(void){
  int i,j,k,l,c;
  unsigned long a,x,count=1;
  //  unsigned short cc[K]={0};
  unsigned short m[K],mm[T]={0},dd[K*D]={0};
  time_t timer;
  FILE *fp,*fq;

  unsigned short g2[7]={1,0,9,0,0,6,4};
  //  unsigned short s[K]={0}; //{4,12,7,8,11,13};
  unsigned short jj[T*2]={0};
  unsigned short ee[10]={1,2,3,4,5,6,7,8,9,10};
  short zz[D]={0};
  //  unsigned short zz[T]={10,97,114,105,97,98,108,101,32,80,111,108,121,110,111,109};
  int y,flg,o1=0;
  OP f={0},h={0},r={0},w={0},aa[8]={0},tt={0},ff={0},g1={0};
  EX hh={0};
  vec v;
  unsigned short d=0;
  time_t t;
  unsigned short gg[K+1]={0};
  oterm rr={0};
  OP r1={0},r2={0},t1={0},t2={0},a1={0},b1={0},a2={0},b2={0};


  
 label:
  for(i=0;i<K+1;i++)
    g[i]=0;
  ginit();
  
  
    w=setpol(g,K+1);
    oprintpol(w);

    #pragma omp parallel for  
  for(i=0;i<D;i++){
    a=trace(w,i);
    if(a==0){
      printf("trace 0 @ %d\n",i);
      goto label;
      //  exit(1);
    }
  }
  printf("@");
  //keygen(g);
  key2(g);

  fq=fopen("H.key","rb");
  
  fread(dd,2,K*D,fq);
  #pragma omp parallel for
  for(i=0;i<D;i++){
  for(j=0;j<K;j++)
    mat[j][i]=dd[K*i+j];
    }
  
  
  for(j=0;j<D;j++){
    flg=0;
    for(i=0;i<K;i++){
      //printf("%d,",mat[i][j]);
      if(mat[i][j]>0)
	flg=1;
      //      printf("\n");
    }
    if(flg==0)
      printf("0 is %d\n",j);
  }
  
  //  exit(1);   
  for(i=0;i<D;i++)
    zz[i]=0;
  
  
  j=0;
  while(j<T*2){
    l=xor128()%D;
    printf("l=%d\n",l);
    if(0==zz[l]){
      zz[l]=1;
      j++;
    }
  }
  

  
  printf("zz=");
  for(i=0;i<D;i++)
    printf("%d,",zz[i]);
  printf("\n");
  //    exit(1);
  //  
  for(i=0;i<K;i++){
    syn[i]=0;
    //#pragma omp parallel for
    for(j=0;j<D;j++){
      //   printf("%u,",zz[jj[j]]);
      syn[i]^=gf[mlt(fg[zz[j]],fg[mat[i][j]])];
    }
       printf("syn%d,",syn[i]);
  }
  printf("\n");
  //    exit(1);  
  for(i=0;i<K;i++)
    printf("mat[%d][1]=%d\n",i,mat[i][1]);
  printf("\n");
  //    exit(1);
  
    
  f=setpol(syn,K);
  printpol(o2v(f));
  printf(" syn=============\n");
  //   exit(1);

  tt.t[0].n=1;
  tt.t[0].a=1;

  
  ff=inv(f,w);
  printpol(o2v(ff));
  printf("locater==========\n");
  //exit(1);
  r2=oadd(ff,tt);
  printpol(o2v(r2));
  printf(" h+x==============\n");
  //  exit(1);
  g1=osqrt(r2,w);
    printpol(o2v(g1));
  printf(" g1!=========\n");
  if(LT(g1).n==0 && LT(g1).a==0){
    printpol(o2v(w));
    printf(" badkey=========\n");
    goto label;
  }

  
  return w;
}


OP ogcd(OP f,OP g){
  OP h={0},ww={0};
  oterm a,b;
  int i=0;


  //oprintpol((f));
  //oprintpol((g));
  //  exit(1);
  
  for(i=0;i<T;i++){// && deg(o2v(f))>=deg(o2v(g))){
  //while(deg(o2v(g))<T){
    //    memset(ss.t,0,DEG);
    h=omod(f,g);
    ww=odiv(f,g);
  f=g;
  g=h;
  //  i++;
  //   if(deg(o2v(g))<T)
  // break;
  }
  // exit(1);

    
  return h;
}



unsigned char chk(OP f){
  int i,j,flg=0;
  vec x={0};

  x=o2v(f);
  for(i=0;i<DEG;i++){
    if(x.x[i]>0){
      flg=1;
      return 1;
    }
  }
  if(flg==0)
    return 0;
  
}
