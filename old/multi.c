//date      :  20160310,20191218,20191220,20191221,20191223,20191224,20191225,20191229
//auther    : the queer who thinking about cryptographic future
//code name : OVP - One Variable Polynomial library with OpenMP friendly
//status    : now in debugging (ver 0.8)
// gcdの停止条件を修正した。vxの停止条件を修正した。
//date      :  20160310,20191218,20191220,20191221,20191223,20191224,20191225,20191229,20191230
//auther    : the queer who thinking about cryptographic future
//code name : OVP - One Variable Polynomial library with OpenMP friendly
//status    : now in debugging (ver 0.8)
// 0ベクトルが出ないように生成多項式のトレースチェックを入れた。
   
//date      :  20160310
//auther    : the queer who thinking about cryptographic future
//code name : OVP - One Variable Polynomial library with OpenMP friendly
//status    : now in debugging (ver 0.1)

#include "chash.cpp"
#include "lu.c"



#define K 256
#define DEG K*2
#define T K/2
#define E 13
#define D 6688

unsigned char tmp[E*K][D]={0};
unsigned char pub[E*K][D]={0};
unsigned char BH[E*K][D]={0};
unsigned short c[K+1]={0};

static unsigned short g[K+1]={0};
  /*
  {1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,
		       0,0,0,0,0,0,0,0,0,0,
		       		       0,0,0,0,0,0,0,0,0,0,
		       		       0,0,0,0,0,0,0,0,0,0,
		       		       0,0,0,0,0,0,0,0,0,0,
		       0,1};
  */
  //  unsigned short g[K+1]={1,1,0,1,1,0,0,1,1,0,1};
//unsigned short g[K+1]={1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};

//,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
//unsigned short g[K+1]={1,0,0,0,1,0,1};

unsigned char A[D][D]={0};
//={1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1}; //={1,5,0,1,7,3,15}; //={1,2,9,4,0,6,4}; // //
unsigned short P[D]={0};
unsigned short inv_P[D]={0};
unsigned short uu;


typedef struct {
  unsigned short n;
  unsigned short a;
} oterm;

typedef struct{
  oterm t[DEG];
} OP;

typedef struct {
  unsigned int x[DEG];
} vec;

typedef union {
  unsigned long long int u[K/4];
  unsigned short s[K];
} SU;

OP ss={0};
OP hh={0};
OP qq={0};
OP p0={0};
SU mat[D]={0};
SU m2[D]={0};
SU syn={0};

static void ginit(void){

  //  g[170]=g[166]=g[162]=g[154]=g[144]=g[142]=g[140]=g[130]=g[128]=g[126]=g[118]=g[114]=g[108]=g[106]=g[104]=g[98]=g[96]=g[94]=g[88]=g[86]=g[82]=g[80]=g[78]=g[72]=g[70]=g[64]=g[62]=g[60]=g[58]=g[54]=g[48]=g[44]=g[42]=g[40]=g[38]=g[36]=g[34]=g[32]=g[24]=g[22]=g[20]=g[12]=g[8]=g[2]=g[0]=1;

  //    g[192]=g[188]=g[186]=g[182]=g[181]=g[180]=g[179]=g[177]=g[176]=g[175]=g[173]=g[170]=g[169]=g[166]=g[164]=g[162]=g[160]=g[159]=g[157]=g[153]=g[149]=g[148]=g[147]=g[146]=g[145]=g[143]=g[140]=g[135]=g[132]=g[131]=g[123]=g[122]=g[120]=g[119]=g[118]=g[117]=g[116]=g[115]=g[114]=g[113]=g[112]=g[111]=g[106]=g[101]=g[99]=g[98]=g[95]=g[94]=g[93]=g[91]=g[90]=g[87]=g[81]=g[80]=g[78]=g[77]=g[76]=g[73]=g[70]=g[61]=g[57]=g[56]=g[54]=g[52]=g[51]=g[50]=g[47]=g[46]=g[43]=g[39]=g[37]=g[33]=g[31]=g[29]=g[27]=g[26]=g[23]=g[17]=g[16]=g[12]=g[10]=g[9]=g[7]=g[5]=g[4]=g[2]=g[0]=1;

  //      g[176]=g[174]=g[171]=g[169]=g[165]=g[161]=g[158]=g[156]=g[154]=g[153]=g[152]=g[149]=g[147]=g[146]=g[145]=g[144]=g[143]=g[142]=g[141]=g[140]=g[139]=g[136]=g[135]=g[134]=g[132]=g[131]=g[128]=g[127]=g[124]=g[122]=g[119]=g[112]=g[111]=g[109]=g[108]=g[100]=g[93]=g[92]=g[90]=g[89]=g[88]=g[87]=g[85]=g[80]=g[77]=g[76]=g[75]=g[73]=g[68]=g[66]=g[64]=g[62]=g[61]=g[60]=g[59]=g[58]=g[57]=g[56]=g[55]=g[54]=g[50]=g[48]=g[46]=g[43]=g[42]=g[39]=g[36]=g[34]=g[27]=g[26]=g[24]=g[19]=g[14]=g[10]=g[9]=g[8]=g[0]=1;
  
  g[256]=g[244]=g[242]=g[240]=g[232]=g[228]=g[224]=g[214]=g[212]=g[210]=g[208]=g[206]=g[204]=g[198]=g[192]=g[186]=g[182]=g[180]=g[174]=g[170]=g[164]=g[162]=g[158]=g[154]=g[150]=g[146]=g[142]=g[138]=g[136]=g[134]=g[130]=g[126]=g[120]=g[118]=g[114]=g[110]=g[104]=g[98]=g[94]=g[90]=g[74]=g[70]=g[68]=g[62]=g[60]=g[50]=g[42]=g[34]=g[30]=g[26]=g[14]=g[12]=g[10]=g[8]=g[6]=g[2]=g[0]=1;
}




void random_permutation(unsigned short* a){
	int i,j,x;
	for(i = 0; i < D; i++){
		a[i] = i;
	}
	for(i = 0; i < D - 2; i++){
		j = (rand() % (D-1-i)) + i + 1;

		x = a[j];
		a[j] = a[i];
		a[i] = x;
	}
	if(a[D-1] == D-1){
		a[D-1] = a[D-2];
		a[D-2] = D - 1;
	}


}


void P2Mat(unsigned short P[M]){
int i,j;
	
  for(i=0;i<D;i++)
      A[i][P[i]]=1;
}


unsigned short b2B(unsigned short b[E]){
  int i,j,k;
  unsigned short a=0;

  for(i=0;i<E;i++)
    a^=(b[i]<<i);

  return a;
}


unsigned short oinv(unsigned short a){
  int i;

  for(i=0;i<M;i++){
    if(gf[mlt(fg[a],i)]==1)
      return (unsigned short)i;
  }

}


unsigned short equ(unsigned short a,unsigned short b){
  int i;

  for(i=0;i<M;i++){
    if(gf[mlt(fg[a],fg[i])]==b)
      break;
  }
  return i;
}


int deg(vec a){
  int i,n=0;

  for(i=0;i<DEG;i++){
    if(a.x[i]>0)
      n=i;
  }
  
  if(n>0){
    return n;
  }else{
    return 0;
  }

}


vec o2v(OP f){
  vec a={0};
  int i,count=0;

  //#pragma omp parallel for
  for(i=0;i<DEG;i++){
    if(f.t[i].a>0)
      a.x[f.t[i].n]=f.t[i].a;
  }

  return a;
}


OP v2o(vec a){
  int i,count=0;
  OP f={0};
  
  //#pragma omp parallel for       
  for(i=0;i<DEG;i++){
    if(a.x[i]>0){
      f.t[i].n=i;
      f.t[i].a=a.x[i];
    }
  }

  return f;
}


vec i2v(unsigned int n){
  vec v;
  int i;
  
  for(i=0;i<32;i++){
  v.x[i]=n%2;
  n=(n>>1);
  }
	
  return v;
}


unsigned int v2i(vec v){
  unsigned int d=0,i;

  for(i=0;i<32;i++){
    d=(d<<1);
    d^=v.x[i];
  }

  return d;
}


OP init_op(OP f){
  int i;

  for(i=0;i<deg(o2v(f))+1;i++){
    f.t[i].a=0;
    f.t[i].n=0;
}
    
  return f;
}


vec init_vec(vec a){
int i;

for(i=0;i<deg(a)+1;i++)
  a.x[i]=0;

return a;
}


vec Setvec(int n){
int i,a,b;
vec v={0};


 for(i=0;i<n;i++){
  v.x[n-1-i]=c[i];
 }

 
  return v;
}


void printpol(vec a){
  int i,n;

  n=deg(a);
  if(n<0){
    printf("baka\n");
    exit(1);
    }
  

  for(i=n;i>-1;i--){
    if(a.x[i]>0){
      printf("%u",a.x[i]);
      if(i>0)
	printf("x^%d",i);
    if(i>0)
      printf("+");
    }
  }
  //  printf("\n");

  return;
}


OP oadd(OP f,OP g){
  vec a={0},b={0},c={0};
  int i,k;
  OP h={0};

  a=o2v(f);
  b=o2v(g);
  if(deg(a)>=deg(b)){
    k=deg(a)+1;
  }else{

    k=deg(b)+1;

}
  //#pragma omp parallel for       
  for(i=0;i<k;i++)
    c.x[i]=a.x[i]^b.x[i];
  h=v2o(c);

  return h;
}


OP oterml(OP f,oterm t){
  int i;
  OP h={0};
  vec test;
  unsigned short n;

  //#pragma omp parallel for
  for(i=0;i<deg(o2v(f))+1;i++){
    
    h.t[i].n=f.t[i].n+t.n;
    h.t[i].a=gf[mlt(fg[f.t[i].a],fg[t.a])];
    
  }

  return h;
}


OP omul(OP f,OP g){
  int i,count=0;
  oterm t;
  OP h={0},e={0},r={0};


  for(i=0;i<deg(o2v(g))+1;i++){
    t=g.t[i];
    e=oterml(f,t);
    h=oadd(h,e);
  }

  return h;
}


oterm LT(OP f){
  int i;
  oterm t;

  for(i=0;i<deg(o2v(f))+1;i++){
    if(f.t[i].a>0){
      t.n=f.t[i].n;
      t.a=f.t[i].a;
    }
  }

  return t;
}


oterm LTdiv(OP f,oterm t){
  oterm tt,s;

  tt=LT(f);
  if(tt.n<t.n){
    s.n=0;
    s.a=0;
  } else if(tt.n==t.n){
    s.n=0;
    s.a=equ(t.a,tt.a);
  }else if(tt.n>t.n){
    s.n=tt.n-t.n;
    s.a=equ(t.a,tt.a);
    //printf("%u\n",s.a);
}

  return s;    
}


OP omod(OP f,OP g){
  int i=0,j,n;
  OP h={0},e={0};
  oterm a,b={0},c;

 

  if(deg(o2v(f))==0 || deg(o2v(g))==0){
    printf("baka^\n");
    exit(1);
  }
  //  printf("in omod\n");
  //exit(1);  

  //  a=LT(f);
  b=LT(g);
  //printpol(o2v(g));
  while(deg(o2v(f))>=deg(o2v(g))){

    c=LTdiv(f,b);
    //    printf("%dx^%d\n",c.a,c.n);
    //    exit(1);
    ss.t[i]=c;
    i++;
    h=oterml(g,c);
      f=oadd(f,h);

      }
  //exit(1);

  return f;
}


OP opow(OP f,int n){
  int i;
  OP g={0};

  g=f;
  //memcpy(g.t,f.t,sizeof(f.t));
  
  for(i=1;i<n;i++)
    g=omul(g,f);

  return g;
}


OP opowmod(OP f,OP mod,int n){
  OP g;
  int i;
  
  g=omod(opow(f,n),mod);
  
  return g;
}


unsigned short trace(OP f,unsigned short x){
  int i;
    unsigned short u=0;

  for(i=0;i<deg(o2v(f))+1;i++){
    u^=gf[mlt(fg[f.t[i].a],mltn(f.t[i].n,fg[x]))];
  }
  
  return u;
}


// invert of polynomial
OP inv(OP a,OP I){
  OP d,x={0},s={0},q={0},r={0},t={0};
  int i;

  d=I;
  //memcpy(d.t,I.t,sizeof(d.t));
 
  //  d = I;
  //  x = {0};
  s.t[0].a = 1;
  while (deg(o2v(r)) > 0){
    memset(ss.t,0,DEG);
    r=omod(d , a);
    //printpol(o2v(ss));
    //    exit(1);
    q=t;
    //memcpy(q.t,ss.t,sizeof(q.t));

    // r = d % a;
    d=a;
    //memcpy(d.t , a.t,sizeof(d.t));
    a=r;
    //memcpy(a.t , r.t,sizeof(a.t));
     //a = r;
    t = oadd(x,omul(q , s));
    x=s;
    //memcpy(x.t,s.t,sizeof(x.t));
    //    x = s;
    s=t;
    //memcpy(s.t,t.t,sizeof(s.t));
    //    s = t;
  }
  //  gcd = d;  // $\gcd(a, n)$
  memset(ss.t,0,DEG);
  omod(I,d);
  
    memset(ss.t,0,DEG);
 return omod(oadd(x , I) , ss);
}


OP ToHorner(OP f){
  OP h;
  int i;
  vec v;

  v=o2v(f);


  return h;
}


OP vx(OP f,OP g){
  OP h={0};
  OP v[K]={0},vv={0};
  oterm a,b;
  int i,j;

  v[0].t[0].a=0;
  v[0].t[1].n=0;
  v[1].t[0].a=1;
  v[1].t[1].n=0;

  a=LT(f);
  b=LT(g);
  //printf("in vx\n");
  //  exit(1);
  i=2;
       for(i=2;i<T+2;i++){
    memset(ss.t,0,DEG);
    if(deg(o2v(f))>=deg(o2v(g))){
  h=omod(f,g);
    }
    
  printf("ss=");
  v[i]=oadd(v[i-2],omul(ss,v[i-1]));
  printf("-------");
  memset(f.t,0,sizeof(f.t));
  f=g;
  //memcpy(f.t,g.t,sizeof(g.t));
  memset(g.t,0,sizeof(g.t));
  g=h;
  
  if(deg(o2v(v[i]))<=T){
    vv=v[i];
  } else {
    
    printf("-------");
    break;
  }
   }

  
  return vv;
}



vec genrandompol(int n){
  vec x={0};
  int i,j=0,k;

  x=init_vec(x);
  for(i=0;i<K+1;i++)
    x.x[i]=xor128()%2;
  
  return x;
}


OP ogcd(OP f,OP g){
  OP h;
  oterm a,b;
  int i=0;

  a=LT(f);
  b=LT(g);
  //printpol(o2v(f));
  //printpol(o2v(g));
  //  exit(1);
  
  while(i<T && deg(o2v(f))>=deg(o2v(g))){
    memset(ss.t,0,DEG);
    h=omod(f,g);
  //printpol(o2v(h));
  //  exit(1);
  f=g;
  g=h;
  //  memcpy(f.t,g.t,sizeof(g.t));
  //memcpy(g.t,h.t,sizeof(h.t));
  /*
  a=LT(f);
  printf("%dx%d\n",a.a,a.n);
  b=LT(g);
  printf("%dx%d\n",b.a,b.n);
  */
  i++;
   if(deg(o2v(g))<=T-1)
     break;
  }
  // exit(1);

    
  return h;
}


OP bibun(vec a){
  OP w[T]={0};
  OP l={0},t={0};
 int i,j,k,n;
 vec tmp={0};
 
 
 n=deg(a);
 printf("n=%d\n",n);
 if(n==0){
   printf("baka8\n");
  exit(1);
 }
 
 for(i=0;i<T;i++){
   w[i].t[0].a=a.x[i];
   w[i].t[0].n=0;
   w[i].t[1].a=1;
   w[i].t[1].n=1;
   //printpol(o2v(w[i]));
 }
 //  exit(1);

 for(i=0;i<T;i++){
   tmp.x[0]=1;
   t=v2o(tmp);
   for(j=0;j<T;j++){
     if(i!=j)
       t=omul(t,w[j]);
   }
   //printpol(o2v(t));
   if(deg(o2v(t))==0){
     printf("baka9\n");
     exit(1);
   }
   l=oadd(l,t);
 }

 
  return l;
}


vec chen(OP f){
  vec e={0};
  int i,count=0,n,x=0;
  unsigned short y[256]={0},z;
  oterm d={0};
  OP g={0};

  
//  e=o2v(f);
n=deg(o2v(f));
  for(x=0;x<D;x++){
    z=0;
    for(i=0;i<n+1;i++){
       if(f.t[i].a>0)
      z^=gf[mlt(mltn(f.t[i].n,fg[x]),fg[f.t[i].a])];
    }
    if(z==0){
      e.x[count]=x;
      count++;
    }
    //    printf("%d\n",x);
  }
  //printpol(e);
  // exit(1);
 

return e;
}

OP decode(OP f,OP s){
int i,j,k;
 OP r={0},h={0},w={0},e={0},l={0};
 oterm t1,t2,d1,a,b;
 vec x={0};
 unsigned short d=0;
 
 
 printf("in decode\n");
 r=vx(f,s);
 // exit(1);
 
 if(deg(o2v(r))==0){
   printf("baka12\n");
    exit(1);
 }
 k=0;
 x=chen(r);

 for(i=0;i<T;i++){
  printf("x[%d]=1\n",x.x[i]);
  if(x.x[i]==0)
    k++;
  if(k>1){
    printf("baka0\n");
    exit(1);
  }
 }
 //  printf("\n");
 
  printf("あっ、でる！\n");  
  exit(1);
 
  if(deg(o2v(r))<T){
    printpol(o2v(r));
   printf("baka5 deg(r)<T\n");
   exit(1);
 }
 
 
 // printpol(x);
 //   exit(1);

 w=bibun(x);
 //  w=oterml(w,d1);
 //printpol(o2v(w));
printf("@@@@@@@@@\n");
//exit(1);

 h=ogcd(f,s);
 //printpol(o2v(h));
 // exit(1);
t1=LT(r);
t2.a=t1.a;
t2.n=0;
/*
 if(deg(o2v(w))==0){
   printf("baka3");
   exit(1);
 }
*/
 if(deg(o2v(w))==0){
   printpol(o2v(w));
 }
l=oterml(w,t2);
/*
 if(deg(o2v(l))==0){
   printf("baka5\n");
   exit(1);
 }
*/
// printf("\n");
 
 for(i=0;i<T;i++){
   //printf("%d ",trace(l,x.x[i]));
   if(trace(l,x.x[i])==0 || trace(h,x.x[i])==0){
     printpol(o2v(l));
     printpol(o2v(h));
     if(trace(l,x.x[i])==0 && trace(h,x.x[i])>0){
     printf("l=%d is 0 deg(x)=%d\n",i,deg(x));
       }else if(trace(l,x.x[i])>0 && trace(h,x.x[i])==0){
	 printf("h=0\n");
     }else if(trace(l,x.x[i])==0 && trace(h,x.x[i])==0){
       printf("l=0 & h=0\n");
     }
     if(x.x[i]==1)
       printf("x=1\n");
     exit(1);
   }
 }
 // printf("\n\n");
 for(i=0;i<T;i++){
   //printf("%d ",trace(h,x.x[i]));
   if(trace(h,x.x[i])==0){
     printpol(o2v(l));
     printpol(o2v(h));
     if(deg(x)<T-1){
	 printf("h=%d %d is 0 too few error deg(x)=%d\n",i,x.x[i],deg(x));
	 exit(1);
     }
   }
 }

 printf("%d\n",deg(x)+1);
 
 //    exit(1);
for(i=0;i<deg(x)+1;i++){
  //  if(x.x[i]>0){
  e.t[i].a=gf[mlt(fg[trace(h,x.x[i])],oinv(trace(l,x.x[i])))];
e.t[i].n=x.x[i];

// }
 }

   for(i=0;i<T;i++)
    if(gf[trace(h,x.x[i])]==0)
      printf("h=0");
   //printf("\n");
  for(i=0;i<T;i++)
    if(gf[oinv(trace(l,x.x[i]))]==0)
      printf("l=0\n");
  //  printf("\n");

  // exit(1);

 /*
   if(x.x[i]==0){
     printf("baka!");
     exit(1);
   }
 */

  // printf("\n");

 
return e;
}




OP setpol(unsigned short *f,int n){
OP g;
vec a;
 int i;

 
 memset(c,0,sizeof(c));
   // for(i=0;i<K+1;i++)
   // c[i]=0;
 memcpy(c,f,sizeof(c));
 //  for(i=0;i<n;i++)
 //  c[i]=f[i];
  a=Setvec(n);
  g=v2o(a);


return g;
}


void det(){
  OP f,h={0},w,u;
  unsigned short cc[K+1]={0},d[2]={0};
  SU HH[D]={0};
  int i,j,a,b;
  oterm t={0};
  vec e;
  
  memcpy(cc,g,sizeof(g));
  k=cc[K];
  w=setpol(g,K+1);

  //#pragma omp parallel for       
  for(i=0;i<D;i++){
  
  a=trace(w,i);
  cc[K]=k;
  
  cc[K]^=a;
  f=setpol(cc,K+1);
  
  h.t[0].a=i;
  h.t[0].n=0;
  h.t[1].a=1;
  h.t[1].n=1;
  
  memset(ss.t,0,sizeof(ss.t));
  omod(f,h);
  
  b=oinv(a);
  t.a=gf[b];
  t.n=0;
  
  u=oterml(ss,t);
  e=o2v(u);
  

#pragma omp parallel for 
  for(j=0;j<K;j++){
      HH[i].s[j]=e.x[K-1-j];  
    }

  }

  
  
  //#pragma omp parallel for
  
  for(j=0;j<D;j++){
#pragma omp parallel for
    for(i=0;i<K/4;i++){
	m2[j].u[i]=mat[j].u[i]=HH[j].u[i];
      }
    }
  
  
}



void bdet(){
  int i,j,k,l,ii;
  unsigned char dd[E*K]={0};
  FILE *ff;
  

  ff=fopen("Hb.key","wb");
  

    for(i=0;i<D;i++){
      for(j=0;j<K;j++){
	l=mat[i].s[j];
#pragma omp parallel for 
	for(k=0;k<E;k++){
	  BH[j*E+k][i]=l%2;
	  l=(l>>1);
	}
      }
    }

  for(i=0;i<D;i++){
    for(j=0;j<E*K;j++){
      printf("%d,",BH[j][i]);
      dd[j]=BH[j][i];
    }
    fwrite(dd,1,E*K,ff);
    printf("\n");
  }
  fclose(ff);
}
 

void pubkeygen(){
  int i,j,k,l;
  FILE *fp;
  unsigned char dd[E*K]={0};

  
  fp=fopen("pub.key","wb");
  for(i=0;i<E*K;i++){
    for(j=0;j<D;j++){
      #pragma omp parallel for 
      for(k=0;k<E*K;k++){
	tmp[i][j]^=cl[i][k]&BH[k][j];
      }
    }
  }
  P2Mat(P);

  for(i=0;i<E*K;i++){
      for(j=0;j<D;j++){
#pragma omp parallel for 
      for(k=0;k<D;k++)
      pub[i][j]^=tmp[i][k]&A[k][j];
    }
  }
  for(i=0;i<D;i++){
#pragma omp parallel for 
    for(j=0;j<E*K;j++){
     dd[j]=pub[j][i];
    }
    fwrite(dd,1,E*K,fp);
  }
  fclose(fp);  
}


void Pgen(){
  unsigned int i,j;
  FILE *fp;
  

  fp=fopen("P.key","wb");
  random_permutation(P);
  for(i=0;i<D;i++)
    inv_P[P[i]]=i;
  fwrite(P,2,D,fp);
  fclose(fp);
  
  for(i=0;i<D;i++)
    printf("%d,",inv_P[i]);
  printf("\n");
  

  fp=fopen("inv_P.key","wb");
  fwrite(inv_P,2,D,fp);
  fclose(fp);

}


void keygen(){
  
  int i;
  FILE *fp;
  

  det();
  makeS();
  bdet();
  Pgen();
  pubkeygen();  
 
}

void key2(){
  FILE *fp;
  unsigned short dd[K]={0};
  int i,j,k,l;

  
  printf("鍵を生成中です。４分程かかります。\n");
  fp=fopen("H.key","wb");
  det();
  
  for(i=0;i<D;i++){
    for(j=0;j<K;j++){
  	dd[j]=mat[i].s[j];
      }
    fwrite(dd,2,K,fp);
    }
  

  
  fclose(fp);
  
}


void encrypt(){


}


void decrypt(){



}


int main(int argc,char **argv){
  int i,j,k,l,c;
  unsigned long a,x,count=1;
  //  unsigned short cc[K]={0};
  unsigned short m[K],mm[T]={0},dd[K*D]={0};
  time_t timer;
  FILE *fp,*fq;

  unsigned short g2[7]={1,0,9,0,0,6,4};
  //  unsigned short s[K]={0}; //{4,12,7,8,11,13};
  unsigned short jj[T]={0};
  unsigned short ee[10]={1,2,3,4,5,6,7,8,9,10};
  unsigned short zz[D]={0};
  //  unsigned short zz[T]={10,97,114,105,97,98,108,101,32,80,111,108,121,110,111,109};
  int y,flg,o1=0;
  OP f={0},h={0},r={0},w={0},aa[T]={0},tt={0};
  vec v;
  unsigned short d=0;
  time_t t;
  unsigned short D2[K][D]={0};  

  //  unsigned short syn[K]={4,12,7,8,11,13};
  //  unsigned short g[K+1]={1,0,0,0,1,0,1};

  //  makegf(D);
  //  makefg(D);
  srand(clock()+time(&t));
  ginit();
  w=setpol(g,K+1);
#pragma omp parallel for  
  for(i=0;i<D;i++){
    a=trace(w,i);
    if(a==0){
      printf("trace 0\n");
      exit(1);
    }
  }
  keygen();
  //key2();
  for(i=0;i<K;i++){
    for(j=0;j<D;j++){
      printf("%d,",mat[j].s[i]);
      printf("\n");
    }
  }
  //  exit(1);
  
  

  /*
  do{
  flg=0;
  v=genrandompol(K+1);
  for(i=0;i<K+1;i++)
    g[i]=(unsigned short)v.x[i];
    for(i=0;i<K+1;i++){
      w=setpol(g,K+1);
      a=trace(w,i);
    if(a==0)
      flg=1;
      }
  } while(deg(v)<K || flg==1);
  */

  w=setpol(g,K+1);
  printf("すげ、オレもうイキそ・・・\n");

  uu=0;
  for(i=0;i<D;i++){
  a=trace(w,i);
  if(a==0){
    printf("trace 0\n");
        exit(1);
  }
  }
  //printpol(o2v(w));
  //    exit(1);

  /*
  fp=fopen(argv[1],"rb");
  fq=fopen(argv[2],"wb");
  
  
  while((c=fread(zz,1,T,fp))>0){
  
  for(i=0;i<D;i++){
    d=trace(w,(unsigned short)i);
    printf("%d,",d);
    if(d==0){
      printf("%i bad trace 0\n",i);
      exit(1);
    }
  }
  printf("\n");
  //exit(1);
  */
  
  j=0;
  while(j<T){
    flg=0;
    if(jj[j]==0){
      l=xor128()%D;
   	for(k=0;k<T;k++){
	  if(l==jj[k])
	    flg=1;
	}
	if(flg==0){
	  jj[j]=l;
	  j++;
	}
    }
  }
  for(i=0;i<T;i++)
    zz[jj[i]]=1;
  
  
  // det(g);


  fq=fopen("H.key","rb");

  fread(dd,2,K*D,fq);
#pragma omp parallel for
  for(i=0;i<D;i++){
    for(j=0;j<K;j++){
      mat[i].s[j]=dd[K*i+j];
    }
  }
  

  /*
  for(i=0;i<K;i++){
    for(j=0;j<D;j++){
      printf("%d,",mat[j].s[i]);
      printf("\n");
    }
  }
  */
  /*
  for(j=0;j<D;j++){
    flg=0;
    for(i=0;i<K;i++){
	//printf("%d,",mat[i][j]);
      if(mat[j].s[i]>0)
	    flg=1;
	  //      printf("\n");
    }
  }
  */

    if(flg==0)
      printf("0 is %d\n",j);


  //  exit(1);   

  
  printf("zz=");
  //
  
  for(i=0;i<K;i++){
    syn.s[i]=0;
    //#pragma omp parallel for
    for(j=0;j<T;j++){
      printf("%u,",zz[j]);	  
      syn.s[i]^=gf[mlt(fg[zz[jj[j]]],fg[mat[jj[j]].s[i]])];
    }
  }

  
    //    printf("%d,",syn[i]);

printf("\n");
  //    exit(1);  
  for(i=0;i<K;i++){
    for(j=0;j<16;j++)
      printf("mat[%d][1]=%d\n",i,mat[1].s[i]);
  printf("\n");
  //  exit(1);
  }
  
  f=setpol(syn.s,K);

  r=decode(w,f);

  
  
  for(i=0;i<T;i++){
    mm[i]=r.t[i].a;
    if(i==0){
     printf("e=%d %d %s\n",r.t[i].a,r.t[i].n,"う");
    }else{     printf("e=%d %d %s\n",r.t[i].a,r.t[i].n,"お");
}
    if(r.t[i].a==0 || r.t[i].n==0){
      printf("------------------\n");
      
       printf("err=%d i=%d\n",o1,i);
      exit(1);
    }
  }
  printf("         っ！！\n");
    for(i=0;i<D;i++){
    if(zz[i]>0)
      o1++;
  }
  printf("err=%d\n",o1);


  return 0;
  
}
