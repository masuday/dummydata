using LinearAlgebra
using Distributions

X0=[
 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
]'
X1=[
 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0
 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1
]'
X2=[
 1 1 1 0 0 1 1 0 0 0 1 1 0 0 0
 0 0 0 1 1 0 0 1 1 1 0 0 1 1 1
]'
X12 = [
 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0
 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0
 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0
 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0
 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0
 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1
]'
X3 = [
 350 400 360 350 340 390 340 410 430 390 400 320 330 390 420
]'
y = vec([
 970 1000  980  980  970  990  950  980  990  980  990  940  930 1000 1000
]')
n = length(y)

# full model
Xfull = [X0 X1 X2 X12 X3]
bfull = (Xfull'*Xfull) \ Xfull'*y

# constraints as SAS
X1m = [(X1[:,1]-X1[:,3]) (X1[:,2]-X1[:,3]) ]
X2m = X2[:,1] - X2[:,2]
X12m = [(X12[:,1]-X12[:,2]-X12[:,5]+X12[:,6]) (X12[:,3]-X12[:,4]-X12[:,5]+X12[:,6])]
X = [X0 X1m X2m X12m X3]
b = (X'*X) \ X'*y

# TypeIII SS
SST = y'*y
SSR = b'*X'*X*b
SSM = y'*y - y'*ones(n,1)*inv(ones(n,1)'*ones(n,1))*ones(n,1)'*y
SSE = SST - SSR
Q1 = X*inv(X'*X)*X'    # SSR = y'*Q1*y
dfe = n - rank(X)

# mu
Xr = [X1m X2m X12m X3]
br = (Xr'*Xr) \ Xr'*y
Q2 = Xr*inv(Xr'*Xr)*Xr'
dfr = rank(Q1) - rank(Q2)
R = y'*(Q1 - Q2)*y
Fval = (R/dfr)/(SSE/dfe)
pval = 1 - cdf(FDist(dfr,dfe),Fval)

# effect A
Xr = [X0 X2m X12m X3]
br = (Xr'*Xr) \ Xr'*y
Q2 = Xr*inv(Xr'*Xr)*Xr'
dfr = rank(Q1) - rank(Q2)
R = y'*(Q1 - Q2)*y
Fval = (R/dfr)/(SSE/dfe)
pval = 1 - cdf(FDist(dfr,dfe),Fval)

# effect B
Xr = [X0 X1m X12m X3]
br = (Xr'*Xr) \ Xr'*y
Q2 = Xr*inv(Xr'*Xr)*Xr'
dfr = rank(Q1) - rank(Q2)
R = y'*(Q1 - Q2)*y
Fval = (R/dfr)/(SSE/dfe)
pval = 1 - cdf(FDist(dfr,dfe),Fval)

# effect A*B
Xr = [X0 X1m X2m X3]
br = (Xr'*Xr) \ Xr'*y
Q2 = Xr*inv(Xr'*Xr)*Xr'
dfr = rank(Q1) - rank(Q2)
R = y'*(Q1 - Q2)*y
Fval = (R/dfr)/(SSE/dfe)
pval = 1 - cdf(FDist(dfr,dfe),Fval)

# effect cov
Xr = [X0 X1m X2m X12m]
br = (Xr'*Xr) \ Xr'*y
Q2 = Xr*inv(Xr'*Xr)*Xr'
dfr = rank(Q1) - rank(Q2)
R = y'*(Q1 - Q2)*y
Fval = (R/dfr)/(SSE/dfe)
pval = 1 - cdf(FDist(dfr,dfe),Fval)
