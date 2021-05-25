K=GF(2^13,'a')
BIN.<x> = K[]
#GF(2)[]


while 1:
      poly = BIN.random_element(256);
      if poly.is_irreducible():
      	 break;
print(poly)
 
