function [m,M_fp] = inverse_normal_radon(d,dt,h,q,N,flow,fhigh,mu);
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
% m -- taup domain results
% M_fp -- f-p domain results

[nt,nh] = size(d);
nq = max(size(q));

if N==2; h=h/max(abs(h));end;
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
Q = eye(nq);
ilow = max(ilow,2);

%3. Calculation of each frequency domain values 
for ifreq=ilow:ihigh
    f = 2.*pi*(ifreq-1)/nfft/dt;
    L = exp(i*f*(h.^N)'*q);
    y = D(ifreq,:)';
    x =  L'*y;
    M(ifreq,:) = x;
    fp(ifreq,:) = abs(x'/max(abs(x')));
    M(nfft+2-ifreq,:) = conj(x)';
end
M_fp = fp';
M(nfft/2+1,:) = zeros(1,nq);
m = real(ifft(M,[],1));
m = m(1:nt,:);
return