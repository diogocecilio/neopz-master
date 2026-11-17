import numpy as np, math
rng=np.random.default_rng(1)
R=np.array([[1/math.sqrt(3)]*3,[math.sqrt(2/3),-1/math.sqrt(6),-1/math.sqrt(6)],[0,1/math.sqrt(2),-1/math.sqrt(2)]])
for k in range(5):
    E=10**rng.uniform(3,6); nu=rng.uniform(0,.45); G=E/(2*(1+nu)); K=E/(3*(1-2*nu))
    phi=math.radians(rng.uniform(10,40)); c=rng.uniform(1,50); sn=math.sin(phi); cs=math.cos(phi)
    s=np.sort(rng.normal(0,50,3))[::-1]; rho=rng.uniform(1,80); beta=rng.uniform(0,math.pi/3)
    xi=(12*c*cs-math.sqrt(6)*rho*(3+sn)*math.cos(beta)+3*math.sqrt(2)*rho*(-1+sn)*math.sin(beta))/(4*math.sqrt(3)*sn)
    ss=np.array([xi,rho*math.cos(beta),rho*math.sin(beta)])
    # check that sigma* lies on Phi1 = 0 (Eq. 51-53)
    sig=R.T@ss; phi1=(sig[0]-sig[2])+(sig[0]+sig[2])*sn-2*c*cs
    d=R@s-ss; Drhw=np.diag([3*K,2*G,2*G]); d2=d@np.linalg.solve(Drhw,d)
    s1,s2,s3=s
    e54=(1/(sn**2*K))*(3*math.sqrt(2)*rho*math.cos(beta)*(sn+3)+math.sqrt(3)*(-3*math.sqrt(2)*rho*math.sin(beta)*(sn-1)-12*c*cs+4*(s1+s2+s3)*sn))**2 \
        +54/G*(math.sqrt(2)*(s2-s3)-2*rho*math.sin(beta))**2+6/G*(math.sqrt(6)*(2*s1-s2-s3)-6*rho*math.cos(beta))**2
    print(f"Phi1(sigma*)={phi1:.2e}  d2={d2:.6e}  Eq54={e54:.6e}  Eq54/d2={e54/d2:.6f}")
