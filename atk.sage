n=4
eta=60
alpha=20
nh=10
rho=20
theta=2
gam=n*eta
ell=4
d=3

def test():
  p=[random_prime(2^eta,False,2^(eta-1)) for i in range(n)]
  x0=prod(p)
  g=[random_prime(2^alpha,False,2^(alpha-1)) for i in range(n)]
  invg=[inverse_mod(g[i],p[i]) for i in range(n)]
  h=[ZZ.random_element(2^(nh-1),2^nh) for i in range(n)]
  
  m=[ZZ.random_element(2^alpha) for i in range(theta)]+[0 for i in range(n-theta)]
  r=[ZZ.random_element(2^rho) for i in range(n)]
  w=(sum([h[i]*m[i]*invg[i]*x0/p[i] for i in range(theta)])+
     sum([h[i]*r[i]*x0/p[i] for i in range(n)])) % x0
  
  print "x0=",x0
  print "p=",p
  print "g=",g
  print "m=",m
  print "w=",w
  pg=prod([g[i] for i in range(theta)])
  print "pg=",pg

  print "\nBasic attack: we should have nh+rho+2*theta*alpha<eta"
  print "nh+rho+2*theta*alpha=%d, eta=%d" % (nh+rho+2*theta*alpha,eta)
 
  rhoR=gam-eta+nh+rho
  B=2^rhoR
  M=Matrix([[B,w],
            [0,x0]])
  ML=M.LLL()
  print "rec pg=",abs(ML[0,0]/B),abs(ML[0,0]/B)==pg

  # extended attack
  m=Matrix([[ZZ.random_element(2^alpha) for i in range(theta)]+[0 for i in range(n-theta)]
     for j in range(ell)])
  r=Matrix([[ZZ.random_element(2^rho) for i in range(n)] for j in range(ell)])

  w=[(sum([h[i]*m[j,i]*invg[i]*x0/p[i] for i in range(theta)])+
             sum([h[i]*r[j,i]*x0/p[i] for i in range(n)])) % x0 
            for j in range(ell)]

  M=Matrix(ZZ,ell+1,ell+1)
  for i in range(ell):
    M[i,i]=B
    M[i,-1]=w[i]
  M[-1,-1]=x0
  ML=M.LLL()
  MLg=ML[:ell,:ell]/B
  print "\nExtended attack: we should have theta*alpha*(1+1/ell)+nh+rho<eta"
  print "theta*alpha*(1+1/ell)+nh+rho=",N(theta*alpha*(1+1/ell)+nh+rho),"eta=",eta
  print "rec pg=",abs(MLg.det()),abs(MLg.det())==pg

  # with multiple vectors
  x=Matrix([[ZZ.random_element(2^alpha) for i in range(theta)] for k in range(d)])
  w=Matrix([[(sum([h[i]*m[j,i]*x[k,i]*invg[i]*x0/p[i] for i in range(theta)])+
             sum([h[i]*r[j,i]*x0/p[i] for i in range(n)])) % x0 
            for j in range(ell)] for k in range(d)])
  M=Matrix(ZZ,ell+d,ell+d)
  for i in range(ell):
    M[i,i]=B
    for k in range(d):
      M[i,ell+k]=w[k,i]
  for i in range(ell,ell+d):
    M[i,i]=x0

  ML=M.LLL()
  MLg=Matrix(ZZ,ML[:ell,:ell]/B)
  print "\nWith multiple vectors: we should have theta*alpha*(1/d+1/ell)+nh+rho<eta"
  print "theta*alpha*(1/d+1/ell)+nh+rho=",N(theta*alpha*(1/d+1/ell)+nh+rho),"eta=",eta
  rpg=abs(MLg.det())
  print "rec pg=",rpg,rpg==pg
  
  print "With factoring:"
  print "  Normalized messages:"
  for i in range(theta):
    mg=Matrix(Integers(g[i]),m[:,i]).T
    print " ",mg/mg[0,0]

  print "  Recovered messages:"
  for i in range(theta):
    MLgi=MLg.change_ring(Integers(g[i]))
    print " ",MLgi.right_kernel().matrix()

  print "Without factoring:"
  print "  Normalized message:"
  v=Matrix(Integers(pg),[[crt([m[j,i] for i in range(theta)],[g[i] for i in range(theta)]) for j in range(ell)]])
  print " ",v[0]/v[0,0]

  print "  Recovered message:"
  MLext=Matrix(ZZ,ell,2*ell)
  MLext[:ell,:ell]=MLg
  for i in range(ell):
    MLext[i,ell+i]=pg
  print " ",MLext.right_kernel().matrix()[0][:ell]
