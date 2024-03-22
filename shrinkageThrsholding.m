function m=shrinkageThrsholding(m,iter,iterMax,a)
m=m;
m0=m;
mMax=max(abs(m(:)));
T=a*(iterMax-iter)/iterMax;
m1=(abs(m)-T*mMax)>0;
m=m1.*m;
end
