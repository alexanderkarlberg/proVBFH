import numpy as np
exec(open('coeffs.py').read().split("for (k,l,name)")[0].replace("int(sys.argv[1]) if len(sys.argv)>1 else 1","3"))
MW,GW,MH,GH=80.379,2.085,125.,4.03e-3
q1k1=dot(q1+k1,q1+k1); q1k2=dot(q1+k2,q1+k2); k12=dot(k1+k2,k1+k2)
def amps(GW):
    A=2*(2*MW**2/(q1k1-MW**2+1j*MW*GW)+2*MW**2/(q1k2-MW**2+1j*MW*GW)+3*MH**2/(k12-MH**2+1j*MH*GH)+1)
    B=1/(q1k1-MW**2+1j*MW*GW)*MW**2/(MW**2-1j*MW*GW)
    C=1/(q1k2-MW**2+1j*MW*GW)*MW**2/(MW**2-1j*MW*GW)
    return A,B,C
for GWv in [GW, 0.0]:
    A,B,C=amps(GWv)
    M=Mtensor(k1,k2,q1,A,B,C,True)
    tot={}
    for k in (1,2,3):
        for l in (1,2,3):
            tot[(k,l)]=S(G(P1,q1,k,True),G(P2,q2,l,True),M)
    print(f"Gamma_W={GWv}:  " + "  ".join(f"F{k}F{l}={tot[(k,l)].real: .3e}" for (k,l) in tot))
    print("   max |Im| of individual terms:", max(abs(v.imag) for v in tot.values()))
