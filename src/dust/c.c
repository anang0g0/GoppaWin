//date 20200331:xgcdの終了条件が２つになってしまう。ogcdとxgcdで使い分ける。
//date : 20200326 鍵生成が det と deta で高い確率で一致する。detaは並列処理。
//date 20200229 : pattarson algorithm implementation ver 1.0
// xgcd & osqrtを追加した
//date      :  20160310,20191218,20191220,20191221,20191223,20191224,20191225,20191229,20191230
//auther    : the queer who thinking about cryptographic future
//code name :  一変数多項式演算ライブラリのつもり
//status    : now in debugging (ver 0.8)
// 0ベクトルが出ないように生成多項式のトレースチェックを入れた。
//date      :  20160310
//auther    : the queer who thinking about cryptographic future
//code name : OVP - One Variable Polynomial library with OpenMP friendly
//status    : now in debugging (ver 0.1)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <sys/mman.h>

#include <omp.h>

#include "8192.h"
#include "global.h"
#include "struct.h"



#include "chash.c"
#include "lu.c"
#include "sha3.c"


extern unsigned long xor128(void);
extern int mlt(int x,int y);
extern int mltn(int n,int a);
extern void makeS();

//#pragma omp threadprivate(mat)
//シンドロームのコピー
unsigned short sy[K]={0};
//Goppa多項式
static unsigned short g[K+1]={0};
unsigned short zz[N]={0};


/*
static unsigned short g[K+1]={1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,
		       0,0,0,0,0,0,0,0,0,0,
		       		       0,0,0,0,0,0,0,0,0,0,
		       		       0,0,0,0,0,0,0,0,0,0,
		       		       0,0,0,0,0,0,0,0,0,0,
		       0,1};
*/

//ランダム多項式の生成
static void ginit (void){
  int j, count = 0;
  unsigned short gg[K + 1] = { 0 };
  //unsigned short g[K+1]={0};

  printf ("in ginit\n");
  

  g[K] = 1; //xor128();
  g[0] = rand()%N;
  while (count < ((K / 2) - 1))
    {
      printf ("@\n");
      j = rand() % (K - 1);
      if (j < K && j > 0 && g[j] == 0)
	{
	  g[j] = rand()%N;
	  count++;
	}
    }

  
  for (j = 0; j < K + 1; j++)
  gg[j] = g[K - j];

  memcpy(g,gg,sizeof(g));
  

}



//ランダム置換の生成（Niederreoter 暗号における置換）
void random_permutation (unsigned short *a){
  int i, j, x;

  for (i = 0; i < N; i++)
    {
      a[i] = i;
    }
  for (i = 0; i < N - 2; i++)
    {
      j = (rand () % (N - 1 - i)) + i + 1;

      x = a[j];
      a[j] = a[i];
      a[i] = x;
    }
  if (a[N - 1] == N - 1)
    {
      a[N - 1] = a[N - 2];
      a[N - 2] = N - 1;
    }


}


//配列から置換行列への変換
void P2Mat (unsigned short P[N]){
  int i;

  for (i = 0; i < N; i++)
    A[i][P[i]] = 1;
}


unsigned short b2B (unsigned short b[E]){
  int i;
  unsigned short a = 0;

  for (i = E-1; i > -1 ; i--)
    a ^= (b[E-i-1] << i);

  return a;
}


//有限体の元の逆数
unsigned short oinv (unsigned short a){
  unsigned short  i;


  for (i = 0; i < M; i++)
    {
      if (gf[mlt (fg[a], i)] == 1)
	return (unsigned short) i;
    }
  
  printf("no return %u\n",a);
  return 0;
}


//aに何をかけたらbになるか
unsigned short equ (unsigned short a, unsigned short b){
  int i;


  for (i = 0; i < M; i++)
    {
      if (gf[mlt (fg[a], fg[i])] == b)
	break;
    }
  return i;
}


//最終の項までの距離
int distance (OP f){
  int i, j=0;

  for (i = 0; i < DEG; i++)
    {
      if (f.t[i].a > 0)
	      j = i;
    }

  return j;
}


//項の数
int terms (OP f){
  int i, count = 0;

  for (i = 0; i < DEG; i++)
    if (f.t[i].a > 0)
      count++;

  return count;
}


//多項式の次数(degのOP型)
int odeg( OP f){
  int i, j = 0, k;


  //k=terms(f);
  for (i = 0; i < 512; i++)
    {
      if (j < f.t[i].n && f.t[i].a > 0)
	j = f.t[i].n;
    }

  return j;
}


//多項式の次数(default)
int deg (vec a){
  int i, n = 0;

  //#pragma omp parallel for
  for (i = 0; i < 512; i++)
    {
      if (a.x[i] > 0)
	n = i;
    }


  return n;
}


//OP型からベクトル型への変換
vec o2v (OP f){
  vec a = { 0 };
  int i;

  //  printf("o2v");
  for (i = 0; i < 512; i++)
    {
      if (f.t[i].a > 0 && f.t[i].n < DEG)
	      a.x[f.t[i].n] = f.t[i].a;
    }


  return a;
}


//ベクトル型からOP型への変換
OP v2o (vec a){
  int i;
  OP f = { 0 };


  for (i = 0; i < DEG; i++)
    {
      if (a.x[i] > 0)
	{
	  f.t[i].n = i;
	  f.t[i].a = a.x[i];
	}
    }


  return f;
}


//整数からベクトル型への変換
vec i2v (unsigned int n){
  vec v;
  int i;


  for (i = 0; i < 32; i++)
    {
      v.x[i] = n % 2;
      n = (n >> 1);
    }


  return v;
}


//ベクトル型から整数への変換
unsigned int v2i (vec v){
  unsigned int d = 0, i;


  for (i = 0; i < 32; i++)
    {
      d = (d << 1);
      d ^= v.x[i];
    }


  return d;
}


//配列からベクトル表現の多項式へ変換する
vec Setvec (int n){
  int i;
  vec v = { 0 };


  for (i = 0; i < n; i++)
    {
      v.x[n - 1 - i] = c[i];
    }


  return v;
}


void printvec(vec v){
  int i,j;
  
  for(i=0;i<deg(v)+1;i++){
    if(v.x[i]>0)
    printf("g[%d]=%d;\n",i,v.x[i]);
  }

}


//多項式を表示する（OP型）
void oprintpol (OP f){
  int i, n;

  n = distance (f);
  printf ("n=%d\n", n);
  printf ("terms=%d\n", terms (f));
  printf ("deg=%d\n", odeg(f));

//exit(1);

  for (i = n; i > -1; i--)
    {
      if (f.t[i].a > 0)
	printf ("%ux^%u+", f.t[i].a, f.t[i].n);
    }
//printf("\n");
//exit(1);
  return;
}


//多項式を表示する(default)
void printpol (vec a){
  int i, n;

  n = deg (a);
  if (n < 0)
    {
      printf ("baka\n");
      exit (1);
    }


  for (i = n; i > -1; i--)
    {
      if (a.x[i] > 0)
	{
	  printf ("%u", a.x[i]);
	  //if (i > 0)
	    printf ("x^%d", i);
	    //if (i > 0)
	    printf ("+");
	}
    }
  //  printf("\n");

  return;
}

OP vadd(vec f,vec g){
  int i,j;

  for(i=0;i<512;i++)
    f.x[i]^=g.x[i];

  return v2o(f);
}


//多項式の足し算
OP oadd (OP f, OP g){
  vec a = { 0 }  , b =  { 0 }  , c =  { 0 };
  int i, j, k,l=0;
  OP h = { 0 };

  a = o2v (f);
  //exit(1);
  b = o2v (g);
  //  oprintpol((g));
  //  exit(1);
  if (odeg (f) >= odeg (g))
    {
      k = odeg (f) + 1;
    }
  else
    {

      k = odeg (g) + 1;

    }
  for (i = 0; i < k; i++){
    c.x[i]=a.x[i]^b.x[i];
  }
  // c.x[i] = a.x[i] ^ b.x[i];
  h = v2o (c);

  return h;
}


//項の順序を降順に揃える
OP sort (OP f){
  oterm o = { 0 };
  int i, j, k;


  k = terms (f);
  for (i = 0; i < k + 1; ++i)
    {
      for (j = i + 1; j < k + 1; ++j)
	{
	  if (f.t[i].n > f.t[j].n)
	    {
	      o = f.t[i];
	      f.t[i] = f.t[j];
	      f.t[j] = o;
	    }
	}
    }

  return f;
}


//リーディングタームを抽出(OP型）
oterm oLT (OP f){
  int i, k,j;
  oterm  s =  { 0 };


  k = terms (f);
  s = f.t[0];
  for (i = 0; i < k + 1; i++)
    {
      //printf("a=%d %d\n",f.t[i].a,f.t[i].n);
      if (f.t[i].a > 0)
	{
	  printf ("in LT=%d %d\n", s.a, s.n);
	  for (j = i; j < k + 1; j++)
	    {
	      if (s.n < f.t[j].n)
		{
		  s.n = f.t[j].n;
		  s.a = f.t[j].a;
		}

	       //  else{
	        // t=s;
	        // }

	    }
	}
    }
  //  exit(1);


  return s;
}



//多項式を足し算する（OP型）
OP add (OP f, OP g){
//  vec a={0},b={0},c={0};
  unsigned long long int i, j, n1 = 0, n2 = 0, m1 = 0, count = 0;
  OP h = { 0 };
  oterm o1 = { 0 }, o2 = { 0 };


  n1 = terms (f);
  printf ("n1=%d\n", n1);
  n2 = terms (g);
  printf ("n2=%d\n", n2);
  if (n1 > n2)
    {

    }

  oprintpol (f);
  printf (" fff==============\n");
  oprintpol (g);
  printf (" ggg==============\n");
  o1 = oLT (f);
  o2 = oLT (g);
  printf ("LTadd==%d %d\n", o1.n, o2.n);
  m1 = n1 + n2;
  printf ("m1=%d\n", m1);
  // exit(1);

  for (i = 0; i < n1 + 1; i++)
    {
      for (j = 0; j < n2 + 1; j++)
	{
	  if (f.t[i].n == g.t[j].n && g.t[j].a > 0 && f.t[i].a > 0)
	    {
	      o1 = oLT (f);
	      o2 = oLT (g);
	      printf ("LT==%d %d\n", o1.n, o2.n);
	      printf ("f.n==%d %d %d %d\n", f.t[i].n, g.t[j].n, i, j);
	      f.t[i].a = 0;
	      g.t[j].a = 0;
	    }
	}
    }
  for (i = 0; i < n2 + 1; i++)
    {
      if (g.t[i].a > 0)
	{
	  h.t[count++] = g.t[i];
	  g.t[i].a = 0;
	}
    }
  for (i = 0; i < n1 + 1; i++)
    {
      if (f.t[i].a > 0)
	{
	  h.t[count++] = f.t[i];
	  f.t[i].a = 0;
	}

    }

  h = sort (h);
  /*
     for (i=0; i<count; ++i) {
     for (j=i+1; j<count; ++j) {
     if (h.t[i].n > h.t[j].n) {
     oo =  h.t[i];
     h.t[i] = h.t[j];
     h.t[j] = oo;
     }
     }
     }
   */
  if (odeg( h) > 0)
    oprintpol (h);
  printf (" addh==============\n");
  //   exit(1);

  return h;
}


//多項式を項ずつ掛ける
OP oterml (OP f, oterm t){
  int i, k;
  OP h = { 0 };
  vec test;
  unsigned short n;

  k = distance (f);
  for (i = 0; i < k + 1; i++)
    {
      h.t[i].n = f.t[i].n + t.n;
      h.t[i].a = gf[mlt (fg[f.t[i].a], fg[t.a])];
    }

  
  return h;
}


//多項式の掛け算
OP omul (OP f, OP g){
  int i, count = 0, k;
  oterm t = { 0 };
  OP h = { 0 }, e = { 0 }, r = { 0 };

  //printf("omul\n");
  if (deg( o2v(f)) > deg( o2v(g)))
    {
      //printf("q");
      k = deg( o2v(f));
    }
  else
    {
      //printf("q");      
      k = deg( o2v(g));
    }
  //printf("a");
  for (i = 0; i < k + 1; i++)
    {
      
      t = g.t[i];
      e = oterml (f, t);
      h = oadd (h, e);
    }

  //  printpol(o2v(h));
  //printf(" debug======\n");
  //exit(1);

  return h;
}


//リーディグタームを抽出(default)
oterm LT (OP f){
  int i, k;
  oterm t = { 0 };


  k = odeg( (f));
  for (i = 0; i < k + 1; i++)
    {
      //printf("a=%d %d\n",f.t[i].a,f.t[i].n);
      if (f.t[i].a > 0)
	{
	  t.n = f.t[i].n;
	  t.a = f.t[i].a;
	}
    }

  return t;
}


//多項式の最後の項を抽出
oterm LT2 (OP f){
  int i, k;
  oterm t = { 0 };

  t.n = f.t[0].n;
  t.a = f.t[0].a;

  return t;
}


//多項式を単行式で割る
oterm LTdiv (OP f, oterm t){
  oterm tt = { 0 } , s = { 0 };

  tt = LT (f);
  if (tt.n < t.n)
    {
      s.n = 0;
      s.a = 0;
    }
  else if (tt.n == t.n)
    {
      s.n = 0;
      s.a = equ (t.a, tt.a);
    }
  else if (tt.n > t.n)
    {
      s.n = tt.n - t.n;
      s.a = equ (t.a, tt.a);
      //printf("%u\n",s.a);
    }
  else if (t.n == 0 && t.a > 0)
    {
      s.a = gf[mlt (fg[tt.a], oinv (t.a))];
      s.n = tt.n;
    }

  return s;
}


//モニック多項式にする
OP coeff (OP f, unsigned short d){
  int i, j, k;
  vec a, b;


  k=odeg((f))+1;
  for (i = 0; i < k; i++)
    f.t[i].a = gf[mlt (fg[f.t[i].a], oinv (d))];


  return f;

}


//多項式の剰余を取る
OP omod (OP f, OP g){
  int i = 0, j, n, k;
  OP h = { 0 }, e = { 0 };
  oterm a, b = { 0 }, c = { 0 };


  n = LT (g).n;
  if (LT(f).n == 0)
    {
      printf ("baka^\n");
      //return f;
      exit (1);
    }
  if (LT(g).n == 0)
    {
      printf ("baka('A`)\n");
      //return g;
      exit (1);
    }
  if (LT(f).n < LT(g).n)
    {
      //    exit(1);
      return f;
    }

  printf("in omod\n");
  //exit(1);

  k = LT (g).n;
  b = LT (g);



  printf ("b=========%dx^%d\n", b.a, b.n);
  printpol (o2v (g));
  if (b.a == 0 && b.n == 0)
    {
      printf ("double baka\n");
      exit (1);
    }

  printf ("\nin omod g=============%d\n", odeg( (g)));
  while (LT(f).n > 0 && LT(g).n > 0)
    {
      printf ("in!\n");
      //    exit(1);

      c = LTdiv (f, b);
      //   printf("c========%dx^%d\n",c.a,c.n);
      //    exit(1);

      printpol (o2v (g));
      printf ("\ng=================\n");

      h = oterml (g, c);
      printpol (o2v (h));
      printf ("\n");
      printf ("modh===================\n");

      printpol (o2v (f));
      printf ("\nmodf===================\n");
      //     exit(1);

      f = oadd (f, h);
      if (odeg( (f)) > 0)
	printpol (o2v (f));
      printf ("\nff1=====================\n");
      if (odeg((f)) == 0 || odeg((g)) == 0)
	{
	  printf ("blake1\n");
	  break;
	}

      if (c.n == 0 || b.n==0)
	break;
    }


  return f;
}


//多項式の商を取る
OP odiv (OP f, OP g){
  int i = 0, j, n, k;
  OP h = { 0 }, e = { 0 }, tt = { 0 }, o = { 0 };
  oterm a, b = { 0 }, c = { 0 };


  o = f;

  if (LT(f).n == 0 && LT(g).a == 0)
    {
      printf ("baka^\n");
      //return f;
      exit (1);
    }
  if (LT (g).a == 0)
    {
      printf ("baka--\n");
      exit (1);
    }
  if (LT(g).n == 0 && LT(g).a > 1)
    return coeff (f, LT(g).a);

  printf ("in odiv\n");
  //printpol(o2v(g));

  //exit(1);
  k = odeg(g);
  b = LT (g);
  printpol (o2v (g));
  printf ("in odiv1 g===========%d %d\n", b.a, b.n);
  if (b.a == 1 && b.n == 0)
    return f;
  if (b.a == 0 && b.n == 0)
    {
      printf ("baka in odiv\n");
      exit (1);
    }
  if (odeg( (f)) < odeg( (g)))
    {
      return f;
      //  a=LT(f);
    }
  printf ("odiv in b=========%dx^%d\n", b.a, b.n);
  printpol (o2v (g));
  printf ("\nf===================\n");
  // e=omul(g,g);
  //printpol(o2v(e));
  //printf("\ng^2================\n");

  printf ("\nin odiv2 g=============%d\n", odeg( (g)));

  i = 0;
  while (LT(f).n > 0 && LT(g).n > 0){
      //  printf("in!\n");
      //    exit(1);

      c = LTdiv (f, b);
      //    if(c.a>0){
      // printf("in odev c========%dx^%d\n",c.a,c.n);
      //    exit(1);
      if(c.n<DEG)
      tt.t[c.n] = c;

      //printpol(o2v(g));
      //printf("\ng=================\n");

      h = oterml (g, c);
      //printpol(o2v(h));
      //printf("\n");
      //printf("h===================\n");
      //printpol(o2v(f));
      //printf("\nf===================\n");
      //     exit(1);

      f = oadd (f, h);

      //printpol(o2v(f));
      //printf("\nff=====================\n");
      if (odeg((f)) == 0 || odeg((g))==0)
	{
	  printf ("blake2\n");
	  break;
	}

      if (c.n == 0)
	break;
    }
  //exit(1);



  return tt;
}


//多項式のべき乗
OP opow (OP f, int n){
  int i;
  OP g = { 0 };


  g = f;
  //memcpy(g.t,f.t,sizeof(f.t));

  for (i = 1; i < n; i++)
    g = omul (g, f);


  return g;
}


//多項式のべき乗余
OP opowmod (OP f, OP mod, int n){
  OP g;
  int i;

  g = omod (opow (f, n), mod);


  return g;
}


//多項式の代入値
unsigned short trace (OP f, unsigned short x){
  int i,d;
  unsigned short u = 0;


  d = odeg((f));

  for (i = 0; i < d + 1; i++)
    {
      u ^= gf[mlt (fg[f.t[i].a], mltn (f.t[i].n, fg[x]))];
    }


  return u;
}


//停止コマンド
void wait(void)
{
    int n;                    // 読み込む変数はローカルに取るべし
    printf(" (enter number and hit return) "); // 何か表示させたほうが良いだろう
    fflush(stdout);           // just in case
    scanf("%d", &n);          // fgets(line, LINESIZE, stdin); という手も
}


// invert of polynomial
OP inv (OP a, OP n){
  OP d = { 0 }, x = { 0 } , s = { 0 } , q = { 0 } , r = { 0 } , t = { 0 } , u = { 0 } , v = { 0 } , w = { 0 } , tt = { 0 } , gcd = { 0 };
  oterm b = { 0 };
  vec vv = { 0 } , xx = { 0 };


  if (odeg( (a)) > odeg( (n)))
    {
      printf ("baka_i\n");
      exit (1);
    }
  if (LT (a).a == 0)
    {
      printf (" a ga 0\n");
      exit (1);
    }


  tt=n;
  printf ("n=============\n");
  //printpol (o2v (n));
  printf ("\n");
  printf ("a=============\n");
  //printpol (o2v (a));
  printf ("\n");
  //  exit(1);

  d = n;
  x.t[0].a = 0;
  x.t[0].n = 0;
  s.t[0].a = 1;
  s.t[0].n = 0;
  while (odeg((a)) > 0)
    {
      if(odeg((a))>0)
      r = omod (d, a);
      if(LT(a).a==0)
	break;
      if(LT(a).a>0)
      q = odiv (d, a);

      d = a;
      a = r;
      t = oadd (x, omul (q, s));
      //printpol (o2v (a));
      printf ("\nin roop a==================%d\n", odeg( (a)));
      printf ("\n");

      x = s;
      s = t;
    }
  // exit(1);
  //  if(LT(a).a>0){
  d = a;
  a = r;
  //printpol (o2v (a));
  printf ("\nin roop a|==================%d\n", odeg( (a)));
  printf ("\n");

  x = s;
  s = t;

  //printpol (o2v (d));
  printf ("\nout1================\n");
  gcd = d;			// $\gcd(a, n)$
  //printpol (o2v (gcd));
  printf ("\n");
  //printpol (o2v (n));
  printf ("\n");
  printf ("out2===============\n");

  printf ("before odiv\n");
  //w=tt;

  b = LT (w);
  //printpol (o2v (w));
  printf ("\nw=======%d %d\n", b.a, b.n);
  //w=tt;
  v = oadd (x, n);
  //printpol (o2v (v));
  printf ("\n");
  if (LT (v).a == 0)
    {
      printf ("v=============0\n");
    }

  printf ("d==============\n");
  //  } //end of a>0
  w = tt;
  //printpol (o2v (v));
  printf ("\n");
  printf ("ss==============\n");
  //       exit(1);
  // if(odeg((w))>0)
  if (LT (v).n > 0 && LT(w).n>0)
    {
      u = omod (v, w);
    }
  else
    {
      printpol(o2v(v));
      printf(" v===========\n");
      printpol (o2v (x));
      printf (" x==0?\n");
      printpol (o2v (n));
      printf (" n==0?\n");

      exit (1);
    }
  //caution !!
  if (LT (u).a > 0 && LT(d).a > 0)
    {
      u = odiv (u, d);
    }

  if (LT (u).a == 0 || LT (d).a == 0)
    {
      printf ("inv div u or d==0\n");
      // exit(1);
    }
  //u=coeff(u,d.t[0].a);
  //printpol (o2v (u));
  printf ("\nu==================\n");
  if (LT (u).a == 0)
    {
      printf ("no return at u==0\n");
      exit (1);
    }


  return u;
}


//error locater for decode
OP vx(OP f,OP g){
  OP h={0},ww={0};
  OP v[K]={0},vv={0};
  oterm a,b;
  int i,j;

  v[0].t[0].a=0;
  v[0].t[1].n=0;
  v[1].t[0].a=1;
  v[1].t[1].n=0;

  //printf("in vx\n");
  //  exit(1);
  i=0;

  while(1){
    if(odeg((g))==0)
      break;
    // if(odeg((g))>0)
    h=omod(f,g);
    printpol(o2v(h));
    printf(" modh vx==============\n");
    if(LT(g).a==0)
      break;
    //if(LT(g).a>0)
    ww=odiv(f,g);
    printf("ww======= ");
    printpol(o2v(ww));
    printf("\n");
    v[i+2]=oadd(v[i],omul(ww,v[i+1]));
    printf("-------");
    f=g;
    g=h;

    vv=v[i+2];
    printf("vv==");
    printpol(o2v(vv));
    printf(" ll========\n");

    if(odeg((vv))==T)
      break;
    i++;
  }
  printpol(o2v(vv));
  printf(" vv============\n");
  //exit(1);

  return vv;
}




//整数のべき乗
unsigned int ipow (unsigned int q, unsigned int u){
  unsigned int i, m = 1;

  for (i = 0; i < u; i++)
    m *= q;

  printf ("in ipow====%d\n", m);

  return m;
}


//０多項式かどうかのチェック
unsigned char chk (OP f){
  int i, flg = 0;
  vec x = { 0 };

  x = o2v (f);
  for (i = 0; i < DEG; i++)
    {
      if (x.x[i] > 0)
	{
	  flg = 1;
	  return 1;
	}
    }
  if (flg == 0)
    return 0;

exit(1);
}


//decode用の多項式の最大公約数
OP ogcd (OP f, OP g){
  OP h;
  //oterm a, b;
  int i = 0;


  //oprintpol((f));
  //oprintpol((g));
  //  exit(1);

  for(i=0;i<T;i++)
    {
      if(odeg((g))==0)
	break;
      h = omod (f, g);
      if(odeg((h))==T-1){
	printpol(o2v(h));
	printf(" in ogcd=============\n");
	//wait();
	//break;
	return h;
      }
      f = g;
      g = h;
    }
  // exit(1);


  return h;
}


//拡張ユークリッドアルゴリズム
EX xgcd (OP f, OP g){
  OP h = { 0 } , ww = { 0 } , *v, *u;
  oterm a, b;
  int i = 0, j, k,flg=0;
  EX e = { 0 };


  v = (OP*)malloc (sizeof (OP) * (DEG));
  u = (OP*)malloc (sizeof (OP) * (DEG));
  memset (v, 0, sizeof (OP)*DEG);
  memset (u, 0, sizeof (OP)*DEG);


  u[0].t[0].a = 1;
  u[0].t[0].n = 0;
  u[1].t[0].a = 0;
  u[1].t[0].n = 0;
  u[2].t[0].a = 1;
  u[2].t[0].n = 0;

  v[0].t[0].a = 0;
  v[0].t[0].n = 0;
  v[1].t[0].a = 1;
  v[1].t[0].n = 0;


  printpol(o2v(f));
  printf(" f===============\n");
  printpol(o2v(g));
  printf(" s===============\n");
  // exit(1);


  k = 0;
  i=0;
  while(1){
     if (LT (g).n == 0){
       printf("v[%d]=%d skipped\n",i,odeg((v[i])));
       //exit(1);
       break;
     }
    
     if(odeg((g))>0)
      h = omod (f, g);
    
     if(LT(g).a>0)
      ww = odiv (f, g);

      v[i+2] = oadd (v[i], omul (ww, v[i+1]));
      u[i+2] = oadd (u[i], omul (ww, u[i+1]));
      printf ("i+1=%d %d %d g=%d\n", i + 1,odeg((v[i])),T-1,odeg((g)));
      f = g;
      g = h;

      //if(
	 if(odeg((f))==T-1 || odeg((v[i]))==T-1){
	break;
      }
      i++;
   }

  //v[i]=odiv(v[i],h);
  //u[i]=odiv(u[i],h);
  // h.t[0].a=1;
  //h.t[0].n=0;
  printf ("i=%d\n", i);
  printpol (o2v (v[i]));
  printf (" v=============\n");
  printpol (o2v (u[i]));
  printf (" u=============\n");
  printpol (o2v (f));
  printf (" f=============\n");
  //exit(1);

  
  e.d = f;
  e.v=v[i];
  e.u=u[i];

  free(v);
  free(u);


  return e;
}


//拡張ユークリッドアルゴリズム(Tで止まらない)
EX gcd (OP f, OP g){
  OP h = { 0 } , ww = { 0 } , v[3]={0}, u[3]={0};
  oterm a, b;
  int i = 0, j, k;
  EX e = { 0 };

  /*
  v = malloc (sizeof (OP) * (DEG));
  u = malloc (sizeof (OP) * (DEG));
  memset (v, 0, sizeof (OP)*DEG);
  memset (u, 0, sizeof (OP)*DEG);
  */

  u[0].t[0].a = 1;
  u[0].t[0].n = 0;
  u[1].t[0].a = 0;
  u[1].t[0].n = 0;
  u[2].t[0].a = 1;
  u[2].t[0].n = 0;

  v[0].t[0].a = 0;
  v[0].t[0].n = 0;
  v[1].t[0].a = 1;
  v[1].t[0].n = 0;


  //printpol(o2v(f));
  //printpol(o2v(g));
  //  exit(1);


  k = 0;
  //i=1;
  while(odeg((g))>0)
    {
      if (LT (g).a == 0)
	break;
      if(odeg((g))>0)
      h = omod (f, g);
      if(LT(g).a==0)
	break;
      if(LT(g).a>0)
      ww = odiv (f, g);

      v[2] = oadd (v[0], omul (ww, v[1]));
      u[2] = oadd (u[0], omul (ww, u[1]));
      printf ("i+1=%d\n", i + 1);
      f = g;
      g = h;

    }

  //v[i]=odiv(v[i],h);
  //u[i]=odiv(u[i],h);
  // h.t[0].a=1;
  //h.t[0].n=0;
  printf ("i=%d\n", i);
  printpol (o2v (v[i]));
  printf (" v=============\n");
  printpol (o2v (u[i]));
  printf (" u=============\n");
  printpol (o2v (h));
  printf (" h=============\n");
  //exit(1);

  e.d = h;
  e.v = v[1];
  e.u = u[1];

  //free(u);
  //free(v);

  return e;
}

OP init_pol (OP f)
{
  int i;

  for (i = 0; i < DEG; i++)
    {
      f.t[i].a = 0;
      f.t[i].n = 0;
    }

  return f;
}



//多項式の形式的微分
OP bibun (vec a){
  OP w[T * 2] = { 0 };
  OP l = { 0 } , t = { 0 };
  int i, j,  n,id;
  vec tmp = { 0 };



  n = deg (a);
  printf ("n=%d\n", n);
  if (n == 0)
    {
      printf ("baka8\n");
      //  exit(1);
    }

  //
  for (i = 0; i < T; i++)
    {
      w[i].t[0].a = a.x[i];
      w[i].t[0].n = 0;
      w[i].t[1].a = 1;
      w[i].t[1].n = 1;
      //printpol(o2v(w[i]));
    }
  //  exit(1);

  tmp.x[0] = 1;  
  //

  //#pragma omp parallel for private(i,j)
  for (i = 0; i < T; i++)
    {
      t = v2o (tmp);
      //
      for (j = 0; j < T; j++)
	{
	  // #pragma omp critical
	  if (i != j){
	    t = omul (t, w[j]);
	}
	}
	
      //printpol(o2v(t));

      if (odeg( (t)) == 0)
	{
	  printf ("baka9\n");
	  // exit(1);
	}
      l = oadd (l, t);
    }
  

  return l;
}


//chen探索
vec chen (OP f){
  vec e = { 0 };
  int i, count = 0, n, x = 0;
  unsigned short z;


  n = odeg( (f));
//exit(1);
//#pragma omp parallel for private(i)
  for (x = 0; x < N; x++)
    {
      z = 0;
      #pragma omp parallel for reduction (^:z)
      for (i = 0; i < n + 1; i++)
	{
	  if (f.t[i].a > 0)
	    z ^= gf[mlt (mltn (f.t[i].n, fg[x]), fg[f.t[i].a])];
	}
      if (z == 0)
	{
	  e.x[count] = x;
	  count++;
	}
      //    printf("%d\n",x);
    }


  return e;
}


//ユークリッドアルゴリズムによる復号関数
OP decode (OP f, OP s){
  int i, j, k,count=0;
  OP r = { 0 }, w ={ 0 }, e = { 0 }, l = { 0 };
  oterm t1, t2, d1, a, b;
  vec x = { 0 };
  unsigned short d = 0;
  OP h = { 0 };
  EX hh = { 0 };


  printf ("in decode\n");
  printpol (o2v (s));
  printf ("\nsyn===========\n");
  r = vx (f, s);
  //h=ogcd(f,s);

  if (odeg( (r)) == 0)
    {
      printf ("baka12\n");
      exit (1);
    }
  k = 0;
  // exit(1);
  x = chen (r);
  // exit(1);



  for (i = 0; i < T; i++)
    {
      printf ("x[%d]=1\n", x.x[i]);
      if (x.x[i] == 0)
	k++;
      if (k > 1)
	{
	  printf ("baka0\n");
	  printvec(o2v(f));
	  for(i=0;i<N;i++)
	    printf("%d,",zz[i]);
	  exit (1);
	  //return f;
	}
    }
  //exit(1);

  //  printf("\n");

  printf ("あっ、でる！\n");
  //  exit(1);

  if (odeg( (r)) < T)
    {
      printpol (o2v (r));
      printf ("baka5 deg(r)<T\n");
      exit (1);
    }


  w = bibun (x);
  //exit(1);
  //  w=oterml(w,d1);
  printpol(o2v(w));
  printf ("@@@@@@@@@\n");
  //exit(1);

  hh = xgcd (f, s);
  printpol (o2v (hh.d));
  //wait();

  
  //  exit(1);
  t1 = LT (r);

  t2.a = t1.a;
  t2.n = 0;

  if (odeg( (w)) == 0)
    {
      printpol (o2v (w));
    }
  l = oterml (w, t2);

  j=deg(x)+1;
  printf ("%d\n", j);

  //    exit(1);

  for (i = 0; i < j; i++)
    {
      if(x.x[i]>0){
      e.t[i].a = gf[mlt (fg[trace (hh.d, x.x[i])], oinv (trace (l, x.x[i])))];
      //e.t[i].a = gf[mlt (fg[trace (h, x.x[i])], oinv (trace (l, x.x[i])))];
      e.t[i].n = x.x[i];
      }
    }
  printpol(o2v(f));
  printf(" f============\n");
  printpol(o2v(l));
  printf(" l============\n");
  //  exit(1);

  
  for (i = 0; i < T; i++)
    if (gf[trace (h, x.x[i])] == 0)
      printf ("h=0");
  //printf("\n");
  for (i = 0; i < T; i++)
    if (gf[oinv (trace (l, x.x[i]))] == 0)
      printf ("l=0\n");
  //  printf("\n");


  return e;
}



//配列の値を係数として多項式に設定する
OP setpol (unsigned short f[], int n){
  OP g;
  vec a;
  int i;


  memset (c, 0, sizeof (c));
  memcpy (c, f, 2*n);
  a = Setvec (n);

  g = v2o (a);


  return g;
}



unsigned short tr[N]={0};
unsigned short ta[N]={0};

void det2(int i,unsigned short g[]){
  OP f[16]={0}, h[16] = { 0 }, w,u[16]={0};
  unsigned short cc[K + 1] = { 0 }, d[2] = {0};
  int  j, a, b ,k,t1,l=0,flg=0,id;
  oterm t[16] = { 0 };
  vec e[16]={0};
  OP ww[16] = { 0 };


  memcpy(cc,g,sizeof(cc));
  
  k = cc[K];
  w = setpol (g, K + 1);
  
  //  omp_set_num_threads(8);
  id=omp_get_thread_num();
  

  h[id].t[0].n = 0;
  h[id].t[1].a = 1;
  h[id].t[1].n = 1;
  t[id].n = 0;
 
    
  f[id]= setpol (cc, K + 1);

  //a = trace (w, i);
  // cc[K] = k;
  
  cc[K] = k^ta[i];
  //tr[i];
  f[id]= setpol (cc, K + 1);
  
  //f.t[0].a=k^ta[i]; //cc[K];
  h[id].t[0].a = i;
  
  ww[id] = odiv (f[id], h[id]);
  
  //b = oinv (a);
  t[id].a = gf[tr[i]];
  u[id] = oterml (ww[id], t[id]);
  e[id] = o2v (u[id]);

  //#pragma omp critical
  //for (j = 0; j < K; j++){
  //mat[i][j]= e[id].x[j];
  //}
  
  //#pragma omp critical(mat)
  memcpy(mat[i],e[id].x,sizeof(mat[i]));

  
}





void f1(unsigned short g[]){
  int i;
  
  for (i = 0; i < 836; i++)
    {
      det2(i,g);
    }


}
void f2(unsigned short g[]){
  int i;
  
  for (i = 836; i < 1672; i++)
    {
      det2(i,g);
    }


}
void f3(unsigned short g[]){
  int i;
  
  for (i = 1672; i < 2508; i++)
    {
      det2(i,g);
    }


}
void f4(unsigned short g[]){
  int i;
  
  for (i = 2508; i < 3344; i++)
    {
      det2(i,g);
    }


}
void f5(unsigned short g[]){
  int i;
  
  for (i = 3344; i <4180 ; i++)
    {
      det2(i,g);
    }


}
void f6(unsigned short g[]){
  int i;
  
  for (i = 4180; i <5016 ; i++)
    {
      det2(i,g);
    }


}
void f7(unsigned short g[]){
  int i;
  
  for (i = 5016; i <5852 ; i++)
    {
      det2(i,g);
    }


}
void f8(unsigned short g[]){
  int i;
  
  for (i = 5852; i < 6688; i++)
    {
      det2(i,g);
    }


}


int detb(unsigned short g[]){
  int i,j,flg=0;
  /*
  for(i=0;i<N;i++){
    for(j=0;j<K;j++)
      mat[i][j]=0;
  }
  */
#pragma omp parallel num_threads(8)
  {
    #pragma omp sections
    {
      //if(omp_get_thread_num() == 0){
#pragma omp section
      f1(g);
      
      //if(omp_get_thread_num() == 1){
#pragma omp section
      f2(g);
      //}
      //if(omp_get_thread_num() == 2){
#pragma omp section
      f3(g);
      //}
      //if(omp_get_thread_num() == 3){
#pragma omp section
      f4(g);
      //}
      //if(omp_get_thread_num() == 4){
#pragma omp section
      f5(g);
      //}
      //if(omp_get_thread_num() == 5){
#pragma omp section
      f6(g);
      //}
      //if(omp_get_thread_num() == 6){
#pragma omp section
      f7(g);
      //}
      //if(omp_get_thread_num() == 7){
#pragma omp section
      f8(g);
      //}
    }
  }
  printf("enf of detb\n");  
  for(j=0;j<N;j++){
    flg=0;
    for(i=0;i<K;i++){
      //printf("%d,",mat[i][j]);
      if(mat[j][i]>0)
	flg=1;
      //      printf("\n");
    }
    if(flg==0){
      printf("0 is %d\n",j);
      //exit(1);
      return -1;
    }
  }

  return 0;
}



//パリティチェック行列を生成する
int deta (unsigned short g[]){
  int i, j, a, b ,k,t1,l=0,flg=0,id;


  //
  //
  #pragma omp parallel num_threads(8)
    {
#pragma omp for schedule(static)
    for (i = 0; i < N; i++)
      {
	det2(i,g);
      }
    }
  for(j=0;j<N;j++){
    flg=0;
    for(i=0;i<K;i++){
      //printf("%d,",mat[i][j]);
      if(mat[j][i]>0)
	flg=1;
      //      printf("\n");
    }
    if(flg==0){
      printf("0 is %d\n",j);
      //exit(1);
      return -1;
    }
  }
  printf("end2\n");
  // exit(1);
  return 0;
}



unsigned short *base;

//パリティチェック行列を生成する
void det (unsigned short g[]){
  OP f, h = { 0 }, w,u;
  unsigned short cc[K + 1] = { 0 }, d[2] = {0},pp[20][K]={0};
  int i, j, a, b ,k,t1,l=0,flg=0,count=0;
  oterm t = { 0 };
  vec e;
  

  memcpy(cc,g,sizeof(cc));
  /*
  for (i = 0; i < K + 1; i++)
    {
      cc[i] = g[i];
      printf ("%d,", g[i]);
    }
    */
  //printf ("\n");
  //  exit(1);
  //    cc[i]=g[i];
  k = cc[K];
  w = setpol (g, K + 1);

  OP ww = { 0 };

  h.t[0].n = 0;
  h.t[1].a = 1;
  h.t[1].n = 1;
  t.n = 0;
  t1=2*T;
  // #pragma omp parallel for
  /*
  unsigned short tr[N];
  unsigned short ta[N];
  for(i=0;i<N;i++){
    ta[i] = trace (w, i);
    if(ta[i]==0){
  printf("%d %d\n",i,ta[i]);
  exit(1);
}   
    tr[i] = oinv (ta[i]);    
}
  */
    
  //
  f= setpol (cc, K + 1);

  for (i = 0; i < N; i++)
    {

      a = trace (w, i);
      // cc[K] = k;

      cc[K] = k^a;
      //tr[i];
      f= setpol (cc, K + 1);

      //f.t[0].a=k^ta[i]; //cc[K];
      h.t[0].a = i;
      
      ww = odiv (f, h);

      b = oinv (a);
      t.a = gf[b];


      u = oterml (ww, t);
      e = o2v (u);

      // #pragma omp parallel for
      //for (j = 0; j < K; j++)
      //mat[i][j]= e.x[j];
      memcpy(mat[i],e.x,sizeof(mat[i]));
    }
    
    
  for(j=0;j<N;j++){
    flg=0;
    for(i=0;i<K;i++){
      //printf("%d,",mat[i][j]);
      if(mat[j][i]>0)
	flg=1;
      //      printf("\n");
    }
    if(flg==0){
      printf("0 is %d\n",j);
      //exit(1);
    }
  }
  //exit(1);

}


//バイナリ型パリティチェック行列を生成する
void bdet (){
  int i, j, k, l;
  unsigned char dd[E * K] = { 0 };
  FILE *ff;


  ff = fopen ("Hb.key", "wb");


  for (i = 0; i < N; i++)
    {
      for (j = 0; j < K; j++)
	{
	  l = mat[i][j];
	  #pragma omp parallel for
	  for (k = 0; k < E; k++)
	    {
	      BH[j * E + k][i] = l % 2;
	      l = (l >> 1);
	    }
	}
    }
  for (i = 0; i < N; i++)
    {
      #pragma omp parallel for
      for (j = 0; j < E * K; j++)
	{
	  //  printf("%d,",BH[j][i]);
	  dd[j] = BH[j][i];
	}
      fwrite (dd, 1, E * K, ff);
      //printf("\n");
    }
  fclose (ff);

}




//Niederreiter暗号の公開鍵を作る
void pubkeygen (){
  int i, j, k, l;
  FILE *fp;
  unsigned char dd[E * K] = { 0 };


  fp = fopen ("pub.key", "wb");
#pragma omp parallel for private(j,k)
  for (i = 0; i < E * K; i++)
    {
      for (j = 0; j < N; j++)
	{
	  //tmp[i][j]=0;
	  l=0;
	  //#pragma omp parallel for reduction (^:l)
	  for (k = 0; k < E * K; k++)
	    {      
	      tmp[i][j]^= cl[i][k] & BH[k][j];
	      //l^= cl[i][k] & BH[k][j];
	    }
	  //tmp[i][j]=l;
	}
      //}
    }
  P2Mat (P);
  
#pragma omp parallel for
  for (i = 0; i < E * K; i++)
    {
      //  for(j=0;j<N;j++){
      #pragma omp parallel for
      for (k = 0; k < N; k++)
	pub[i][k] = tmp[i][P[k]];	//&A[k][j];
      //    }
    }

#pragma omp parallel for private(j)
  for (i = 0; i < N; i++)
    {
      for (j = 0; j < E * K; j++)
	{
	  dd[j] = pub[j][i];
	}
      fwrite (dd, 1, E * K, fp);
    }
  fclose (fp);

}


//秘密置換を生成する
void Pgen (){
  unsigned int i, j;
  FILE *fp;


  fp = fopen ("P.key", "wb");
  random_permutation (P);
  for (i = 0; i < N; i++)
    inv_P[P[i]] = i;
  fwrite (P, 2, N, fp);
  fclose (fp);

  for (i = 0; i < N; i++)
    printf ("%d,", inv_P[i]);
  printf ("\n");


  fp = fopen ("inv_P.key", "wb");
  fwrite (inv_P, 2, N, fp);
  fclose (fp);

}


//鍵生成
void key2 (unsigned short g[]){
  FILE *fp;
  unsigned short dd[K] = { 0 };
  int i,j,k;

  printf ("鍵を生成中です。４分程かかります。\n");
    fp = fopen ("H.key", "wb");
    i=0;
    do{
      i=deta (g);
    }while(i==-1);
      
  //exit(1);
  for (i = 0; i < N; i++)
    {
      for (j = 0; j < K; j++)
	dd[j] = mat[i][j];
      fwrite (dd, 2, K, fp);

    }
  fclose (fp);
  fp=fopen("sk.key","wb");
  fwrite(g,2,K+1,fp);
  fclose(fp);
  
}


//すべての鍵を生成する
void keygen (unsigned short *g){
  int i;
  FILE *fp;


  key2 (g);
  printf ("end of ky2\n");
  makeS ();
  printf ("end of S\n");
  bdet ();
  printf ("end of bdet\n");
  Pgen ();
  printf ("end of Pgen\n");
  pubkeygen ();

}


//ハッシュ１６進表示
static void byte_to_hex(uint8_t b, char s[23]) {
    unsigned i=1;
    s[0] = s[1] = '0';
    s[2] = '\0';
    while(b) {
        unsigned t = b & 0x0f;
        if( t < 10 ) {
            s[i] = '0' + t;
        } else {
            s[i] = 'a' + t - 10;
        }
        i--;
        b >>= 4;
    }
}




//有限体の元の平方を計算する
int isqrt (unsigned short u){
  int i, j, k;


  for (i = 0; i < M; i++)
    {
      if (gf[mlt (i, i)] == u)
	return i;
    }

printf("来ちゃいけないところに来ました\n");
exit(1);
}


//多項式の平方を計算する
OP osqrt(OP f,OP w){
  int i,j,k,jj,n;
  OP even={0},odd={0},h={0},r={0},ww={0},s={0},tmp={0},t={0};
  oterm o={0};
  vec v={0};

  j=0;
  jj=0;
  k=distance(f);
  for(i=0;i<k+1;i++){
    if(f.t[i].n%2==0){
      even.t[j].n=f.t[i].n/2;
      even.t[j++].a=gf[isqrt(f.t[i].a)];
      printf("a=%d %d\n",f.t[i].a,i);
    }
    if(f.t[i].n%2==1){
      odd.t[jj].n=(f.t[i].n-1)/2;
      odd.t[jj++].a=gf[isqrt(f.t[i].a)];
      printf(" odd %d\n",i);
    }
  }
  printpol(o2v(even));
  printf(" T0============\n");
  printpol(o2v(odd));
  printf(" T1============\n");
  // exit(1);


  k=odeg((w));
  printf("%d\n",k);
  //exit(1);
  j=0;
  jj=0;
  for(i=0;i<k+1;i++){
    if(w.t[i].n%2==0){
      h.t[j].a=gf[isqrt(w.t[i].a)];
      h.t[j++].n=w.t[i].n/2;
      printf("wi==%d %d\n",(w.t[i].n/2),i);
    }
    if(w.t[i].n%2==1){
      r.t[jj].a=gf[isqrt(w.t[i].a)];
      r.t[jj++].n=(w.t[i].n-1)/2;

    }
  }
  //  printpol(o2v(r));
  printf(" sqrt(g1)=======\n");

  //  exit(1);
  if(LT(r).n>0){
    s=inv(r,w);
  }else if(LT(r).n==0){
    printpol(o2v(r));
    printf(" r======0\n");
    wait();
    exit(1);
  }

  tmp=omod(omul(s,r),w);
  if(odeg((tmp))>0){
    printpol(o2v(tmp));
    printf(" r is not inv==========\n");
    wait();
     exit(1);
  }
  if(LT(h).n>0){
    ww=omod(omul(h,s),w);
  }
  if(LT(h).n==0){
    printf("h=========0\n");
    exit(1);
  }
  /*
  if(LT(ww).n==0 && LT(ww).a==0){
    printpol(o2v(s));
    printf(" s===========\n");
    printpol(o2v(w));
    printf(" w==============\n");
    printpol(o2v(r));
    printf(" r===========\n");
    printpol(o2v(h));
    printf(" h============\n");
    printpol (o2v (ww));
    printf (" ww==============\n");
    printf(" wwが０になりました。error\n");
    wait();
    return ww;;
    // exit(1);
  }
  */
  tmp=omod(omul(ww,ww),w);
  if(LT(tmp).n==1){
    printpol(o2v(ww));
    printf(" ww succsess!===========\n");
  }else{
    printpol(o2v(tmp));
    printf(" mod w^2==========\n");
    printpol(o2v(ww));
    printf(" ww^2 failed!========\n");
    printpol(o2v(s));
    printf(" g1^-1==============\n");
    printpol(o2v(w));
    printf(" w==============\n");
    printpol(o2v(h));
    printf(" g0===========\n");
    printpol(o2v(r));
    printf(" r===========\n");
    printf("この鍵では逆元が計算できません。error");
    wait();
    return ww;
    // exit(1);
  }


    //    exit(1);
  //  printpol(o2v(s));
  printf(" g1^-1=========\n");
  // printpol(o2v(h));
  printf(" g0=========\n");
  //exit(1);
  // printpol(o2v(ww));
  printf(" ww==========\n");
  //  exit(1);
   h=ww;
  if(odeg(omod(omul(h,ww),w))==1){
    ww=h;
    h=omod(oadd(even,omul(ww,odd)),w);
    return h;
  }else if(LT(ww).a==0){
    printf("vaka\n");
    exit(1);
  }

  // printpol(o2v(ww));
  printf(" 来ちゃだめなところに来ました\n");

  exit(1);
}


//パターソンアルゴリズムでバイナリGoppa符号を復号する
vec pattarson (OP w, OP f){
  OP g1 = { 0 }, ll = { 0 };
  int i, j, k, l, c,n,count=0;
  int flg, o1 = 0;
  OP tt = { 0 }, ff = { 0 },h={0};
  EX hh = { 0 };
  vec v;
  oterm rr = { 0 };
  OP  r2 =  { 0 },  b2 = { 0 };
  //unsigned short g[K+1]={2,2,12,1,2,8,4,13,5,10,8,2,15,10,7,3,5};


  tt.t[0].n = 1;
  tt.t[0].a = 1;


  ff = inv (f, w);
//  printpol (o2v (ff));
  b2=omod(omul(ff,f),w);
  if(odeg((b2))>0){
    printf("逆元が計算できません。error\n");
    wait();
    exit(1);
  }
  printf ("locater==========\n");
  //exit(1);
  r2 = oadd (ff, tt);
//  printpol (o2v (r2));
  printf (" h+x==============\n");
  //  exit(1);
  g1 = osqrt (r2, w);
  printpol (o2v (g1));
  b2=omod(omul(g1,g1),w);
  if(LT2(b2).a!=LT2(r2).a){
    printpol(o2v(w));
    printf(" w============\n");
    printpol(o2v(r2));
    printf(" r2============\n");
    printpol(o2v(g1));
    printf(" g1============\n");
    printf(" g1は平方ではありません。error");
    wait();
    exit(1);
  }
  printf (" g1!=========\n");
  if (LT (g1).n == 0 && LT (g1).a == 0)
    {
      printpol (o2v (w));
      printf (" badkey=========\n");
      printpol(o2v(g1));
      printf(" g1============\n");
      printf("平方が０になりました。 error\n");
      wait();
      exit(1);
    }
  //exit(1);
  hh = xgcd (w, g1);
  flg=0;
 aa:
  ff = omod (omul (hh.v, g1), w);
 printpol (o2v (ff));
  printf (" beta!=========\n");
  if (odeg( (ff)) != K / 2)
    {
      flg=1;
      printpol(o2v(w));
      printf (" locater function failed!! error\n");
	printf("cannot correct(bad key) error============\n");
	wait();
	//return -1;
    }


  //printpol (o2v (hh.v));
  printf (" alpha!=========\n");
  //exit(1);
  if(odeg((ff))==K/2){
    ll = oadd (omul (ff, ff), omul (tt, omul (hh.v, hh.v)));
  }else if(odeg((ff))==1){
    ll=oadd (omul (ff, ff), omul (tt, omul (hh.v, hh.v)));
  }else{
    printf("locate degree is !=K/2 %d\n",odeg((ff)));
    exit(1);
  }
  if (odeg( (ll)) == 0)
    {
      printf (" locater degree is 0\n");
      exit (1);
    }
  printf ("あっ、でる・・・！\n");
  count = 0;
  printpol(o2v(ll));
  printf(" ll=================\n");
  printpol(o2v(f));
  printf(" syn=================\n");
  v = chen (ll);
  if(count<2*T-1){
    printf(" decode failed error===========\n");
    for(i=0;i<N;i++){
      if(trace(ll,i)==0)
	printf("x=%d\n",i);
    }
  }


  return v;
}


//512bitの秘密鍵を暗号化
void encrypt (char buf[],unsigned char sk[64])
{
  const uint8_t *hash={0};
  sha3_context c={0};
  int image_size=512,i;
  FILE *fp;
//  unsigned short dd=0;



  printf("plain text=");
  for(i=0;i<64;i++)
    printf("%u,",sk[i]);
  printf("\n");
  //  puts(buf);
  //printf("\n");
  //exit(1);

  //scanf("%s",buf);
  sha3_Init256(&c);
  sha3_Update(&c, (char *)buf, strlen(buf));
  hash = sha3_Finalize(&c);

  //j=0;

  for(i=0; i<64; i++) {
    printf("%d", hash[i]);
      //char s[3];
      //byte_to_hex(hash[i],s);

      sk[i]^=hash[i];

  }
  printf("\nencrypt sk=");
  for(i=0;i<64;i++)
    printf("%d,",sk[i]);
  printf("\n");

  fp=fopen("enc.sk","wb");
  fwrite(sy,2,K,fp);
  fwrite(sk,1,64,fp);
  fclose(fp);

}


void decrypt (OP w)
{
  FILE *fp;
  int i,j;
  unsigned char sk[64]={0},err[N]={0};
  unsigned short buf[K]={0},tmp[K]={0};
  OP f={0};
  vec v={0};
  const uint8_t *hash={0};
  sha3_context c={0};
  int image_size=512;


  j=0;
  fp=fopen("enc.sk","rb");

  fread(tmp,2,K,fp);
  fread(sk,1,64,fp);
  fclose(fp);

  for(i=0;i<K;i++)
    buf[i]=tmp[K-i-1];

  printf("in decrypt\n");
  f=setpol(buf,K);
  v=pattarson(w,f);
  //v=o2v(h);
  //  exit(1);


  j=0;
  if(v.x[1]>0 && v.x[0]==0){
    err[0]=1;
    j++;
  }

  printf("j=%d\n",j);
  printf("after j\n");
  for(i=j;i<2*T;i++){
    if(v.x[i]>0 && v.x[i]<N){
      err[v.x[i]]=1;
    }
  }

  char buf0[8192]={0},buf1[10]={0};

  //#pragma omp parallel for
  for(i=0;i<N;i++){
    snprintf(buf1, 10, "%d",err[i] );
    strcat(buf0,buf1);
  }
  puts(buf0);
  printf("vector=%d\n",strlen(buf0));
  //exit(1);
  printf("cipher sk2=");
  for(i=0;i<64;i++)
    printf("%u,",sk[i]);
  printf("\n");

  sha3_Init256(&c);
  sha3_Update(&c, (char *)buf0, strlen(buf0));
  hash = sha3_Finalize(&c);

  j=0;
  printf("hash=");
  for(i=0; i<64; i++) {
    printf("%d", hash[i]);
      //char s[3];
      //byte_to_hex(hash[i],s);

      sk[i]^=hash[i];
  }
  printf("\ndecript sk=");
  for(i=0;i<64;i++)
    printf("%u,",sk[i]);
  printf("\n");
  //  exit(1);

  return;
}



OP synd(unsigned short zz[]){
  unsigned short syn[K]={0},s=0;
  int i,j,t1;
  OP f={0};


  printf("in synd\n");

#pragma omp parallel for //num_threads(8)
  for(i=0;i<K;i++){
    syn[i]=0;
    s=0;
    //#pragma omp parallel num_threads(8)
    for(j=0;j<N;j++){
      syn[i]^=gf[mlt(fg[zz[j]],fg[mat[j][i]])];
    }
    sy[i]=syn[i];//=s;
       printf("syn%d,",syn[i]);
  }
  printf("\n");
  
   for(int j = 0; j < K/2; j++){
    t1 = syn[j];
    syn[j] = syn[K - j - 1];
    syn[K - j - 1] = t1;
  }
  

  f=setpol(syn,K);
  printpol(o2v(f));
  printf(" syn=============\n");
  //  exit(1);

  return f;
}


//64バイト秘密鍵の暗号化と復号のテスト
void test(OP w,unsigned short zz[]){
  int i;
  vec v={0};
  const uint8_t *hash;
  sha3_context c;
  //int image_size=512;
  OP f={0};
  FILE *fp;


  fp=fopen("aes.key","rb");
  /*
      static char base64[] = {
	'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
	'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
	'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X',
	'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f',
	'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
	'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
	'w', 'x', 'y', 'z', '0', '1', '2', '3',
	'4', '5', '6', '7', '8', '9', '+', '/',
      };
  */
      char buf[8192]={0},buf1[10]={0};
      unsigned char sk[64]={0};
      // unsigned short s[K]={0};
      //fread(sk,1,32,fp);
      for(i=0;i<64;i++)
      sk[i]=i+1;

      for(i=0;i<N;i++){
	snprintf(buf1, 10, "%d",zz[i] );
	strcat(buf,buf1);
      }
      puts(buf);
      printf("vector=%u\n",strlen(buf));
      //exit(1);


      printf("sk0=");
      for(i=0;i<64;i++)
	printf("%u,",sk[i]);
      printf("\n");
      //exit(1);

      f=synd(zz);
      v=o2v(f);
      //printf("v=");
      for(i=0;i<K;i++){
	sy[i]=v.x[i];
	printf("%d,",sy[i]);
      }
      printf("\n");

      encrypt(buf,sk);
      decrypt(w);


      sha3_Init256(&c);
      sha3_Update(&c, (char *)buf, strlen(buf));
      hash = sha3_Finalize(&c);


}


//言わずもがな
int main (void){
  int i, j, k, l;
  int count = 0;
  FILE *fp, *fq;
  //unsigned short zz[N] = { 0 };
  int flg, o1 = 0;
  OP f = { 0 }, r = { 0 }, w = { 0 },ff={0},tt={0},h={0};
  EX hh = { 0 };
  vec v;
  //  unsigned short dd[K*N] = {0},gg[K+1]={0};
  time_t t;
  OP r1 = { 0 }, r2 = { 0 };
  OP g1 = { 0 },tmp={0};
  unsigned short s[6]={1,2,3,4,5,6};
  vec oo[3]={0};
  vec u={0};

//公開鍵matはグローバル変数でスタック領域に取る。
//ヒープ領域は使うときはここを有効にしてください。
/*
  mat = malloc (N * sizeof (unsigned short *));
  base=malloc(sizeof(unsigned short)*K*N);
  #pragma omp parallel for
  for(i=0;i<N;i++){
  //連続したメモリ空間にデータを配置
  mat[i]=base+i*K;
  memset(mat[i],0,K);
  }
*/
  
  srand(clock () + time (&t));
  printf ("@");
  //getkey();
  //  exit(1);

label:

  //makeS();
  //exit(1);
  for(i=0;i<3;i++){
    for(j=0;j<512;j++)
    v.x[j]=0;
    }

  for(i=0;i<3;i++){
    u.x[i]=rand()%16;
    oo[i].x[1]=1;
    oo[i].x[0]=rand()%16;
    printpol(oo[i]);
    printf("\n");
  }
  //exit(1);
  
  w=setpol(s,6);
  f=v2o(oo[0]);

  for(i=1;i<3;i++){
    h=v2o(oo[i]);
    f=omul(h,f);
    printf("2\n");
  }
  // exit(1);
  
  h=v2o(u);

  w=omul(f,w);
  printpol(o2v(w));
  printf("2\n");
  w=odiv(w,h);
  h=omod(w,h);
  r=bibun(o2v(f));
  for(i=0;i<8192;i++){
    if(o2v(f).x[i]>0)
    printf("%d\n", gf[mlt (fg[trace (w, o2v(f).x[i])], oinv (trace (r, (o2v(f)).x[i])))]);
  }


  exit(1);
  
  for (i = 0; i < K + 1; i++)
   g[i] = 0;
  ginit ();

  /*
  fp=fopen("sk.key","rb");
  fread(g,2,K+1,fp);
  fclose(fp);
  */

  w = setpol (g, K + 1);
  oprintpol (w);
  //exit(1);

  
  for (i = 0; i < N; i++)
    {
      ta[i] = trace (w, i);
      if (ta[i] == 0)
	{
	  printf ("trace 0 @ %d\n", i);
	  goto label;
	  //exit(1);
	}
    }
  
#pragma omp parallel for
  for(i=0;i<N;i++)
    tr[i]=oinv(ta[i]);
  printf ("@");

  memset(mat,0,sizeof(mat));


  //keygen(g);
  //鍵をファイルに書き込むためにはkey2を有効にしてください。
  //どうしても早くしたい人はdeta()にすること。defaultはdet()
  det(g);
  //exit(1);

  //gccの場合、並列化すると鍵生成が不完全になる。その場合、デバッグデータを出力して終了。
  //再現性がないので余り役に立たない。その場合detを使うこと。
  /*
  i=0;
  do{
    memset(mat,0,sizeof(mat));
    i=deta(g);
  }while(i==-1);
  */
    
  lab:
  //exit(1);
  //key2 (g);
  

  //固定した鍵を使いたい場合はファイルから読み込むようにしてください。  
  /*  
  fq = fopen ("H.key", "rb");
  fread (dd, 2, K * N, fq);
  //#pragma omp parallel for
  for (i = 0; i < N; i++)
    {
      for (j = 0; j < K; j++)
	mat[i][j] = dd[K * i + j];
    }
  fclose (fq);
  */
  
  //#pragma omp parallel for 
  for (j = 0; j < N; j++)
    {
      flg = 0;
      for (i = 0; i < K; i++)
	{
	  //printf("%d,",mat[i][j]);
	  if (mat[j][i] > 0)
	    flg = 1;
	  //      printf("\n");
	}
      if (flg == 0){
	printf ("0 is %d\n", j);
	exit(1);
      }
    }
  


  printf ("すげ、オレもうイキそ・・・\n");


//decode bigin
//#pragma omp parallel for
  for(j=0;j<N;j++){
    flg=0;
    for(i=0;i<K;i++){
      //printf("%d,",mat[i][j]);
      if(mat[j][i]>0)
	flg=1;
      //      printf("\n");
    }
    if(flg==0){
      printf("0 is %d\n",j);
      exit(1);
    }
  }
  //exit(1);

    k=0;
    while(1){
   
      o1=0;

    count=0;

    //  exit(1);

    
    //  for(i=0;i<N;i++)
    memset(zz,0,sizeof(zz));

         
  j=0;
  while(j<T){
    l=xor128()%N;
    //printf("l=%d\n",l);
    if(0==zz[l] && l>0){
      zz[l]=l;
      j++;
    }
  }
    
    
    
  f=synd(zz);
  printpol(o2v(f));
  printf("\n");
  //wait();
  //  exit(1);
  count=0;
  for(i=0;i<N;i++){
    if(zz[i]>0)
      count++;
  }
  printf("%d\n",count);
  //exit(1);

    
  r=decode(w,f);
  //  exit(1);

  for(i=0;i<T;i++){
    if(i==0){
     printf("e=%d %d %s\n",r.t[i].a,r.t[i].n,"う");
    }else if(r.t[i].a==r.t[i].n){
      printf("e=%d %d %s\n",r.t[i].a,r.t[i].n,"お");
    }else if(r.t[i].a!=r.t[i].n){
	printpol(o2v(w));
	printf(" goppa polynomial==============\n");
	for(j=0;j<N;j++)
	  printf("%d,",zz[j]);
	printf("\n");
	for(j=0;j<T;j++)
	  printf("e=%d %d %s\n",r.t[j].a,r.t[j].n,"お");
	exit(1);
	for(l=0;l<N;l++){
	  // printf("%d,",zz[l]);
	  if(zz[l]>0 && zz[l]==l)
	    count++;
	}
	if(count<T){
	  printf("error pattarn too few\n");
	  exit(1);
	}

}
    if(r.t[i].a==0){
      printf("------------------\n");
       printf("err=%d i=%d\n",o1,i);
      printpol(o2v(w));
      printf(" w==============\n");
      printpol(o2v(f));
      printf(" w==============\n");
      for(l=0;l<N;l++)
	printf("%d,",zz[l]);
      printf("\n");
       exit(1);
    }
  }
  o1=0;
  //  #pragma omp parallel for
    for(i=0;i<N;i++){
      if(zz[i]>0)
      o1++;
  }
  printf("err=%dっ！！\n",o1);

  //goto label;
   
  
  
  printf("パターソンアルゴリズムを実行します。何か数字を入れてください。\n");
  //exit(1);
  //wait();
  
  
    //flg=0;
  //  while(1){
  
  
  //for (i = 0; i < N; i++)
  memset(zz,0,sizeof(zz));
  //	zz[i] = 0;
  
  
  j = 0;
  while (j < T * 2)
    {
      l = xor128 () % N;
      printf ("l=%d\n", l);
      if (0 == zz[l]){
	zz[l] = 1;
	j++;
      }
    }
  
  test(w,zz);
  //wait();

  //exit(1);

  for(i=0;i<N;i++)
    printf("%d,",zz[i]);
  printf("\n");
  
  f=synd(zz);
  

  printpol (o2v (f));
  printf (" syn=============\n");
  //exit(1);
  
  
  //バグトラップのためのコード（冗長）
  
	  hh=gcd(w,f);
	  if(odeg((hh.d))>0){
	  printf(" s,wは互いに素じゃありません。\n");
	  wait();
	  goto label;
	  }
	  

      tt.t[0].n = 1;
      tt.t[0].a = 1;


      ff = inv (f, w);
      tmp=omod(omul(ff,f),w);
      if(odeg((tmp))>0){
	printpol(o2v(tmp));
	printf(" inv(h+x)============\n");
	printpol(o2v(w));
	printf(" w============\n");
	printf("この多項式では逆元計算ができません。");
	printf("count=%d\n",k);
	wait();
	goto label;
      }
      
      
      r2 = oadd (ff, tt);
      printpol (o2v (r2));
      printf (" h+x==============\n");
      //  exit(1);
      g1 = osqrt (r2, w);
      printpol (o2v (g1));
      printf (" g1!=========\n");
      r1 = omod(omul(g1,g1),w);
      printpol(o2v(r1));
      printf(" g1^2 mod w===========\n");
      printpol(o2v(r2));
      printf(" r2===========same?\n");
      //scanf("%d",&n);
      if(odeg((r1))!=odeg((r2)) && odeg((g1))>0){
	printpol(o2v(w));
	printf(" w===========\n");
	printpol(o2v(r2));
	printf(" r2===========\n");
	printf("平方根の計算に失敗しました。\n");
	printf("count=%d\n",k);
	wait();
	goto label;
	//exit(1);
      }
      if (odeg((g1)) == 0)
	{
	  printpol(o2v(g1));
	  printf(" sqrt(h+x)==============\n");
	  printpol (o2v (w));
	  printf (" badkey=========\n");
	  printf ("平方根が０になりました。\n");
	  printf("count=%d\n",k);
	  wait ();
	  //exit(1);
	  goto label;
	}

      
      hh = xgcd (w, g1);
      ff = omod (omul (hh.v, g1), w);
      printpol (o2v (ff));
      printf (" beta!=========\n");
      if (odeg( (ff)) != K / 2)
	{
	  printf("deg(l)!=K/2=%d %d %d\n",odeg((ff)),K,k);
	  //exit(1);
	  goto label;
	}
      
      //バグトラップ（ここまで）

      count=0;
      //復号化の本体
      v=pattarson (w, f);
      //エラー表示
      for (i = 0; i < T * 2; i++)
	{
	  if(i==0 && v.x[i]>0){
	    printf ("%d う\n",v.x[i]);
	    count++;
	  }
	  if (v.x[i] > 0 && i>0){
	    ++count;
	  printf ("%d お\n", v.x[i]);
	  }
	}
      printf ("err=%dっ!! \n", count);
      if(count<255){
	printpol(o2v(w));
	for(i=0;i<N;i++){
	  if(zz[i]>0)
	    printf("%d,",zz[i]);
	}
	printf("\n");
      }
      
      //goto lab;
      //wait();
      
      break;
  }

    //free(base);
    //free(mat);
    
    
    return 0;
}


