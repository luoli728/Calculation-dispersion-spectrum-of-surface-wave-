function [m,M_fp] = inverse_Radon_ista(d,dt,h,q,N,flow,fhigh,mu)
% get taup domain results
% input:
% d -- The seismic data
% dt --Time sampling interval
% h -- Seismic trace spacing interval
% q -- densities of the model
% N -- number of modes
% flow -- The lowest frequency
% fhigh -- The highest frequency
% output:
% m -- taup domain results basa on ISTA
% M_fp -- f-p domain results

[nt,nh] = size(d);
nq = max(size(q));
if N==2; h=h/max(abs(h));end
nfft = 2*(2^nextpow2(nt));
%1. The Fourier transform
D = fft(d,nfft,1);
%2. determine the frequency range: ifreq/dt * nt) for the frequency of the real value
M = zeros(nfft,nq);
i = sqrt(-1);
ilow  = floor(flow*dt*nfft)+1;
if ilow < 2; ilow=2; end
ihigh = floor(fhigh*dt*nfft)+1;
if ihigh > floor(nfft/2)+1; ihigh=floor(nfft/2)+1; end
ilow = max(ilow,2);
a = 0.9;
iterMax =2;
%3. Calculation of each frequency domain values 
for ifreq=ilow:ihigh
    lpcg = mu*eye(nq);
    f = 2.*pi*(ifreq-1)/nfft/dt;
    L = exp(i*f*(h.^N)'*q);
    d = D(ifreq,:)';
    r0 = norm(d);
    r = d;
    [~,Nq]=size(L);
    m =zeros(Nq,1);
    dn = zeros(size(d));
    for iter = 1:iterMax
        resi(iter) =  dot(r,r);
        res_r = norm(r-L*m)./r0;
        if(res_r>0.01)
            dd = r-L*m;
            g = lsqr(L,dd,1e-6,1,lpcg);
            m = shrinkageThrsholding(g,iter,100,a); 
            %%(the iterative threshold shrinkage algorithm)
            dn = dn + L*m;
            Mean = mean(g)/nq;
            resfact = sum((g-Mean).*conj(g-Mean))/(nq-1);
            for j = 1:nq
                lpcg(j,j) = 1.0/(resfact+g(j)*conj(g(j)));
            end
        else
            continue;
        end
    end
    gp = lsqr(L,dn,1e-6,1,lpcg);
    r = dn;
    fp(ifreq,:) = abs(gp/max(abs(gp)));
    M(ifreq,:) = gp;
    M(nfft+2-ifreq,:) = conj(gp)';

end
M_fp = fp';
M(nfft/2+1,:) = zeros(1,nq);
m = real(ifft(M,[],1));
m = m(1:nt,:);
return
end