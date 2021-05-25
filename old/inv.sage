#G=GF(256,'a')
#G=GF(2)
G.<a> = GF(2^13)
R.<x> = PolynomialRing(G)
#K.<a>=GF(2**4)
#R.<x> = PolynomialRing(G,implementation='NTL')







gg=x^6+x^4+x^3+x+1
g0=x^6+x^4+1
g1=x^2+1
f=sqrt(g0)
ff=sqrt(g1)
h=inverse_mod(ff,gg)
print("g1^-1=",h)
w=h*f%gg
print("w*w%gg=",(w^2)%gg)
print("w==",w)
es=G('a^12')*x^5+G('a^6')*x^4+G('a')*x^3+G('a^7')*x^2+G('a^9')*x+G('a')
hh=inverse_mod(es,gg)
v=hh+x
print("s^-1+x===",v)

T0=G('a^6')*x^4+G('a^7')*x^2+G('a')
T1=G('a^12')*x^4+G('a')*x^2+G('a^9')
#print("gcd=",gcd(T0,x^32+x))
f=sqrt(T0)
h=sqrt(T1)
print("T0=",f)
print("T1=",h)
r=(f+w*h)%gg
print("r=",r)
#print("r^2%gg=",r*r%gg)
#print("h+x==",v)



