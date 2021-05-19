n=200
eta=669
alpha=80
nh=alpha
rho=300
theta=2
gam=n*eta


def genpseudoprime(eta,etamin=100):
  if eta<2*etamin:
    return random_prime(2^eta,False,511*2^eta/512)
  else:
    return random_prime(2^etamin,False,511*2^etamin/512)*genpseudoprime(eta-etamin,etamin)

def test():
  p=[genpseudoprime(eta) for i in range(n)]
  x0=prod(p)
  gam=x0.nbits()
  g=[random_prime(2^alpha,False,2^(alpha-1)) for i in range(n)]
  invg=[inverse_mod(g[i],p[i]) for i in range(n)]
  h=[ZZ.random_element(2^(nh-1),2^nh) for i in range(n)]
  
  print "n=",n
  print "gamma=",gam
  #print "x0=",x0
  #print "p=",p
  #print "g=",g

  nu=eta-nh-rho
  rhoR=gam-nu
  
  for theta in [1,2,60,80,100,120,140,180][-1:]:
    print "\ntheta=",theta
    ell=1
    d=ell
    while theta*alpha*(1/d+1/ell)+0.03*(d+ell)>nu-10:
      ell+=1
      d=ell
    print "ell=",ell
    pg=prod([g[i] for i in range(theta)])

    # extended attack with multiple vectors
    m=Matrix([[ZZ.random_element(2^alpha) for i in range(theta)]+[0 for i in range(n-theta)]
       for j in range(ell)])
    r=Matrix([[ZZ.random_element(2^rho) for i in range(n)] for j in range(ell)])

    x=Matrix([[ZZ.random_element(2^alpha) for i in range(theta)] for k in range(d)])
    w=Matrix([[(sum([h[i]*m[j,i]*x[k,i]*invg[i]*x0/p[i] for i in range(theta)])+
             sum([h[i]*r[j,i]*x0/p[i] for i in range(n)])) % x0 
            for j in range(ell)] for k in range(d)])

    print "lattice dimension=",ell+d
    M=Matrix(ZZ,ell+d,ell+d)
    for i in range(ell):
      M[i,i]=1
      for k in range(d):
        M[i,ell+k]=w[k,i] >> rhoR
    for i in range(ell,ell+d):
      M[i,i]=x0 >> rhoR

    t=cputime(subprocesses=True)  
    ML=M.LLL()
    MLg=Matrix(ZZ,ML[:ell,:ell])
    print "With multiple vectors: we should have theta*alpha*(1/d+1/ell)<nu"
    print "theta*alpha*(1/d+1/ell)+0.03*(d+ell)=",N(theta*alpha*(1/d+1/ell)+0.03*(d+ell)),"nu=",nu
    rpg=abs(MLg.det())
    print "rec pg=",rpg==pg
    print "time=",cputime(subprocesses=True)-t