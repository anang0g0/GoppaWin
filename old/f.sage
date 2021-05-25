#G=GF(256,'a')
#G=GF(2)
G.<a> = GF(2^4)
R.<x> = PolynomialRing(G)
#K.<a>=GF(2**4)
#R.<x> = PolynomialRing(G,implementation='NTL')







g=x^2+x+G('a')
g0=x^2+G('a')
f=sqrt(g0)
print("f=",f)
g1=G('1')


R=1+f*g1*G('a^7')
w=(f*g1)%g
print("w^2%g=",(w^2)%g)
print("w=",w)
r=1+w*G('a^7')
print("R=",r)


es=G('a^10')*x+G('a^13')
v=G('a^14')*x+1
print("s^-1+x===",v)
print("r^2=",r*r%g)

