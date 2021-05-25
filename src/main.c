#include "oplib.c"

int main(void)
{
  int i, j = 0, count = 0;
  unsigned short zz[N] = {0}, z1[N] = {0};
  OP f = {0}, r = {0}, w = {0};
  vec v, x = {0};
  MTX R = {0}, O = {0}, Q = {0};
  unsigned short s[K + 1] = {0};
  arrayul ary={0};
  //srand(clock());
/*
//while(1){
chash(ary.d);
printf("%d\n",ary.t[1]);
//}
exit(1);
*/
  Pgen();
  do
  {
    memset(S.x, 0, sizeof(S.x));
    memset(inv_S.x, 0, sizeof(inv_S.x));
    for (i = 0; i < K * E; i++)
    {
      for (j = 0; j < K * E; j++)
        S.x[i][j] = xor128() % 2;
    }
  } while (is_reg(S, inv_S.x) == -1);

  //exit(1);
  if (K > N)
    printf("configuration error! K is bigger than N\n");

  // （謎）
  memset(mat, 0, sizeof(mat));

  // 公開鍵を生成する(Berlekamp-Massey用)
  R = pk_gen();
  while(1){
  // エラーベクトルの初期化
  memset(zz, 0, sizeof(zz));
  //重みTのエラーベクトルを生成する
  mkerr(zz, T);
  // 暗号文の生成(s=eH)
  x = sin2(zz, R);
  // 復号化１(m'=sS^{-1})
  r = dec(x.x);
  v = o2v(r);
  for (i = 0; i < K; i++)
    s[i + 1] = v.x[i];

  // Berlekamp-Massey Algorithm
  f = bma(s, K);
  x = chen(f);
  // 平文の表示(m=m'P^{-1})
  ero2(x);
  //wait();
  }

  

  return 0;
}
