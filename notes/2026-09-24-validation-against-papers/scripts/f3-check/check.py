import numpy as np, itertools
rng = np.random.default_rng(1)
g = np.diag([1.,-1,-1,-1])
def dot(a,b): return a@g@b
# Levi-Civita with eps_{0123}=+1 (lower indices), as in the code's explicit T3 formula
eps = np.zeros((4,4,4,4))
for p in itertools.permutations(range(4)):
    eps[p] = np.linalg.det(np.eye(4)[list(p)])
def E_lower(P,q):   # i eps_{mu nu rho sigma} P^rho q^sigma  (true lower components)
    return 1j*np.einsum('mnrs,r,s->mn', eps, P, q)
def code_E(P,q):    # the code's explicit component formula
    T=np.zeros((4,4),complex)
    T[0,1]=P[2]*q[3]-P[3]*q[2]; T[0,2]=-P[1]*q[3]+P[3]*q[1]; T[0,3]=P[1]*q[2]-P[2]*q[1]
    T[1,2]=-P[3]*q[0]+P[0]*q[3]; T[1,3]=P[2]*q[0]-P[0]*q[2]; T[2,3]=-P[1]*q[0]+P[0]*q[1]
    T=T-T.T; return 1j*T
def kin():
    E=6500.; P1=np.array([E,0,0,E]); P2=np.array([E,0,0,-E])
    while True:
        x1,x2=rng.uniform(0.05,0.6,2); p1=x1*P1; p2=x2*P2
        def jet():
            pt=rng.uniform(20,200); y=rng.uniform(-4,4); ph=rng.uniform(0,2*np.pi)
            return np.array([pt*np.cosh(y),pt*np.cos(ph),pt*np.sin(ph),pt*np.sinh(y)])
        j1,j2=jet(),jet(); K=p1+p2-j1-j2; mH=125.
        if K[0]>0 and dot(K,K)>4*mH**2: break
    # split K into two Higgs
    M=np.sqrt(dot(K,K)); pst=np.sqrt(M**2/4-mH**2); c=rng.uniform(-1,1); s=np.sqrt(1-c*c); ph=rng.uniform(0,2*np.pi)
    kr=np.array([M/2, pst*s*np.cos(ph), pst*s*np.sin(ph), pst*c])
    b=K[1:]/K[0]; gam=K[0]/M
    def boost(v):
        bp=b@v[1:]; b2=b@b
        return np.concatenate(([gam*(v[0]+bp)], v[1:]+((gam-1)*bp/b2+gam*v[0])*b))
    k1=boost(kr); k2=K-k1
    q1=j1-p1; q2=j2-p2
    return q1,q2,P1,P2,k1,k2
def Mtensor(k1,k2,q1,A,B,C,lowerit):
    L=(lambda v: g@v) if lowerit else (lambda v: v)
    a=L(2*k1+q1); bb=L(k2-k1-q1); c=L(2*k2+q1); d=L(k1-k2-q1)
    Mlow=A*g+B*np.outer(a,bb)+C*np.outer(c,d)
    return g@Mlow@g           # raise both indices
def S(W1,W2,M):  # W1_{mu nu} M^{mu rho} M*^{nu sigma} W2_{rho sigma}
    return np.einsum('mn,mr,ns,rs->',W1,M,M.conj(),W2)
q1,q2,P1,P2,k1,k2=kin()
print("momentum conservation:", np.abs(k1+k2+q1+q2).max(), " k1^2-mH^2:", dot(k1,k1)-125**2)
print("code eps formula == true lower eps:", np.abs(code_E(P1,q1)-E_lower(P1,q1)).max())
A,B,C=(rng.normal(size=3)+1j*rng.normal(size=3))
for lab,(A_,B_,C_) in {"AA only":(A,0,0),"BB only":(0,B,0),"A+B":(A,B,0),"A+B+C":(A,B,C)}.items():
    E1=E_lower(P1,q1)/(2*dot(P1,q1)); E2=E_lower(P2,q2)/(2*dot(P2,q2))
    corr=S(E1,E2,Mtensor(k1,k2,q1,A_,B_,C_,True))
    code=S(E1,E2,Mtensor(k1,k2,q1,A_,B_,C_,False))
    print(f"F3F3 {lab:8s}: correct {corr.real: .6e}  code-convention {code.real: .6e}  ratio {code.real/corr.real if corr.real else float('nan'): .6f}")
# F1F1 sanity: must agree
for lab,(A_,B_,C_) in {"A+B+C":(A,B,C)}.items():
    G1=lambda P,q,lowerit: ( (np.outer(g@q,g@q) if lowerit else np.outer(q,q))/dot(q,q) - g )
    corr=S(G1(P1,q1,True),G1(P2,q2,True),Mtensor(k1,k2,q1,A_,B_,C_,True))
    code=S(G1(P1,q1,False),G1(P2,q2,False),Mtensor(k1,k2,q1,A_,B_,C_,False))
    print(f"F1F1 {lab}: correct {corr.real: .6e} code {code.real: .6e}")
