F.<x>=GF(2)[]
P.<x>=GF(2^13)[]

def GetGoppaPolynomial(polynomial_ring, polynomial_degree):
	while 1:
		irr_poly = polynomial_ring.random_element(polynomial_degree);
		irr_poly_list = irr_poly.list(); irr_poly_list[-1] = 1;
		irr_poly = polynomial_ring(irr_poly_list);
		if irr_poly.degree() != polynomial_degree:
			continue;
		elif irr_poly.is_irreducible():
			break;				
		else :
			continue;
	
	return irr_poly;


def BinRepr(poly):
	try: 
		poly_ls = list(poly);
	except TypeError:
		bin_repr = bin(poly.integer_representation())[2:];
		bin_repr = bin_repr[::-1];
		print(bin_repr);
		return;
		
	else: 
		for _ in poly_ls:
			bin_repr = bin(_.integer_representation())[2:];
			bin_repr = bin_repr[::-1];
			print(bin_repr);
			

g=GetGoppaPolynomial(P,256)
#f=BinRepr(g)
print('g(x)=',g,'is irreducible');

