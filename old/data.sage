#DataCollectionScript.py
#!/opt/sage-5.6/sage -python
import time
from sage.all import *
############################################
#Functions & Classes to Import
############################################
load("./SAGE Classes/GoppaCode.sage")
load("./SAGE Classes/McElieceCryptosystem.sage")
load("./SAGE Classes/SternsAlgorithm.sage")
def GetRandomMessage(message_length):
message = matrix(GF(2), 1, message_length);
for i in range(message_length):
message[0,i] = randint(0,1);
return message;
def SternsAttack(crypto, encrypted_message, t, p, l):
#Attacker has knowledge of PK and y (and presumably t), and is looking for a codeword of weight w
PK = copy(crypto.public_key());
y = encrypted_message;
#Calculate a parity check matrix for the code G + {0,y}, where y = mG + e
H = (PK.stack(y)).right_kernel().basis_matrix();
w = t;
#Find a weight w codeword
weight_w_codeword = SternsAlgorithm(H, w, p, l);
#Decrypt the message using the codeword found via Stern’s Algorithm
decrypted_message = PK.solve_left((y-weight_w_codeword));
return decrypted_message;
############################################
#Beginning of Script
############################################
#The goal of the script will be to write the parameters and results to a .csv file format
#Table Headers
print "\nn,m,t,Patterson’s Algorithm Wall Time,Patterson’s Algorithm CPU Time,p,l,Stern’s Attack Wall
Time,Stern’s Attack CPU Time,Polynomial";
for k in range(3):
m = k+6;
n = 2**m;
t = floor((2+(2**m-1)/m)/2);
F_2m = GF(n,’Z’);
PR_F_2m = PolynomialRing(F_2m,’X’);
for _ in range(10):
while 1:
irr_poly = PR_F_2m.random_element(t);
if irr_poly.is_irreducible():
break;
crypto = McElieceCryptosystem(n,m,irr_poly);
for i in range(5):
#Get a random message and encrypt it.
message = GetRandomMessage(crypto.goppa_code().generator_matrix().nrows());
encrypted_message = crypto.Encrypt(message);
48
#Patterson’s Algorithm Decoding
patterson_start_wall_time = time.time();
patterson_start_cpu_time = time.clock();
crypto.Decrypt(encrypted_message);
patterson_end_cpu_time = time.clock();
patterson_end_wall_time = time.time()
patterson_wall_time = patterson_end_wall_time - patterson_start_wall_time;
patterson_cpu_time = patterson_end_cpu_time - patterson_start_cpu_time;
H = (copy(crypto.public_key()).stack(encrypted_message)).right_kernel().basis_matrix();
#Stern’s Attack
for j in range(2):
p = j+1;
k_2 = floor((H.ncols()-H.nrows())/2);
l_min = floor(log(k_2,2))-1;
for k in range(7):
l = l_min + k;
if (H.nrows() < l):
l = H.nrows();
stern_start_wall_time = time.time();
stern_start_cpu_time = time.clock();
SternsAttack(crypto,encrypted_message,t,p,l);
stern_end_cpu_time = time.clock();
stern_end_wall_time = time.time()
stern_wall_time = stern_end_wall_time - stern_start_wall_time;
stern_cpu_time = stern_end_cpu_time - stern_start_cpu_time;
#n,m,t,Iteration,Patterson’s Algorithm,p,l,Stern’s Attack,Polynomial";
print "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" % (n,m,t, \
patterson_wall_time,patterson_cpu_time,\
p,l,stern_wall_time,stern_cpu_time,irr_poly);