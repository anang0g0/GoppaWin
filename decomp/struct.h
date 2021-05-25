//符号のパラーメータの指定。通常[N,K,T]として、
//Nは符号の長さ、Kが符号の次元、Tは訂正エラー数
//を表す。ここではDは符号長にしている。
#define N 6688 // (^^;)?
#define M 8192 //有限体の元の数
#define K (128*2) //符号の次元
#define DEG (2*K) //多項式の次数
#define T (K/2) //エラーの数
#define E (13) //拡大体のビット数
#define D (6688) //符号長（短縮符号）
#define F 256*13 //2040

//monomial
typedef struct{
  unsigned short n; //単項式の次数
  unsigned short a; //単項式の係数
} oterm;

//polynomial
typedef struct{
  oterm t[DEG]; //単項式の配列として多項式を表現する
} OP;

typedef struct{
  unsigned int x[DEG]; //配列の添字を次数に、配列の値を係数に持つ多項式の表現
} vec;

//extra gcd
typedef struct{
  OP u; //inverse of polynomial?
  OP v; //error locater
  OP d; //gcd
} EX;


typedef union{ //test(SIMN)
  unsigned long long int u[4];
  unsigned short s[16];
} SU;

