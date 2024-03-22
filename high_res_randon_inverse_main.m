%This program is used to calculate the high resolution 
%dispersion energy spectrum of surface wave based on Ista
% This is a demo for reproducing figure 1-3 in the paper
% Date: July 7th, 2003.

%%
clear ;
close all;

% 1.The layered model parameters
%vs_initial:vs,vp_initial:vp,h_initial:h,freq:frequency,rho:density 
freq = [4:1:50];
vs_initial = [390,540,660];    
vp_initial = [800,1100,1300];
h_initial  = [25,25];
rho = [1950,1950,2000];

%2.load the seismic data
%generate synthetic surface wave seismograms via  
%finite difference solution to the elastic wave equation

load 'shot_m2c.mat'

%%A total trace of 120, the shot point position 12, 
%%the number of sampling points nt=200,dt=0.003,dx=2,
Data=shot_m2c(:,20:120);
seis=Data;
[nt,nx]=size(seis);
dt=0.003;
dx=2;
figure;
set(gcf,'Position', 1.0e+03 *[ 0.9600    0.5741    1.1339    0.5669]);
y=[1:1:nt]*dt;x=[20:1:120]*dx;
wiggle1(y,x,seis,'2');
set(gca,'Fontsize',18);
xlabel('Offset (m)');ylabel('Time (s)');


%3.compute the dispersion spectrum.
%Energy dispersive imaging: x-t Transform to f-v（taup）
t=[1:nt]*dt;
h=[2:2:2*nx];
v_min=50;v_max=1049;dv=1;
v=[v_min:dv:v_max];
p=1./v;
df=1/dt/nt;
%%Normal randon transform 
[taup,fp] = inverse_normal_radon(seis,dt,h,p,1,1,80,2);

%%Randon transform based on ista

[taup1,fp1] = inverse_Radon_ista(seis,dt,h,p,1,1,80,1);

%4.The theory of dispersion curve calculation
%The algorithms are based on:
% Xiaofei Chen(1993) 'A systematic and efficient method of computing
% normal modes for multilayered half-space'.
%
% Hisada,Y,(1994).'An Efficient Method for Computing Green's Functions for
% a Layered Half-Space with Sources and Receivers at Close Depths'.
%
% Lai,C.G.(1998)'Simultaneous Inversion of Rayleigh Phase Velocity and
% Attenuation for Near-Surface Site Characterization'.
bos = getbos(vs_initial,vp_initial);
DV = 1;   
moden=1;
vr1= mat_disperse(h_initial , rho, vp_initial, vs_initial, freq);


%%The dispersion energy mapping (normal)
%%Tuap--
figure;
imagesc(p,t,taup)
set(gca,'Fontsize',16);
colormap('gray');
xlabel('dip(m/s)^-^1 ');ylabel('Time(s)');
%%F-V--
figure;
resolution=100;
Aplot = abs(fp);
[~,nf]=size(Aplot);
[nt,~] = size(seis);
nfft = 2*(2^nextpow2(nt));
m_fv = zeros(size(Aplot));
fplot = m_fv;
cplot = m_fv;

for idx = 1:nf    
     f_temp = (idx-1)/nfft/dt;
     fplot(:,idx) = f_temp;
     cplot(:,idx) = v;
end

[~,ch] = contourf(fplot,cplot,Aplot,resolution);
set(ch,'LineStyle','none');
set(ch,'edgecolor','none');
colormap(jet)
shading flat
set(gca,'XTick',0:5:50);
set(gca,'FontSize',16);
xlim([0 50])
xlabel('Frequency(Hz)','FontSize',16,'Fontweight','normal','color','k');
ylabel('Velocity (m/s)','FontSize',16,'Fontweight','normal','color','k');
hold on;
plot(freq,vr1,'g.','LineWidth',4);

%%The dispersion energy mapping(ISTA)
%%Tuap--
figure;
imagesc(p,t,taup1)
set(gca,'Fontsize',16);
colormap('gray');
xlabel('dip(m/s)^-^1 ');ylabel('Time(s)');
%%F-V--


figure;
resolution=100;
Aplot = abs(fp1);
[nv,nf]=size(Aplot);
[nt,nh] = size(seis);
nfft = 2*(2^nextpow2(nt));
m_fv = zeros(size(Aplot));
fplot = m_fv;
cplot = m_fv;

for idx = 1:nf    
     f_temp = (idx-1)/nfft/dt;
     fplot(:,idx) = f_temp;
     cplot(:,idx) = v;
end

[~,ch] = contourf(fplot,cplot,Aplot,resolution);
set(ch,'LineStyle','none');
set(ch,'edgecolor','none');
colormap(jet)
shading flat
set(gca,'XTick',0:5:50);
set(gca,'FontSize',16);
xlim([0 50])
xlabel('Frequency(Hz)','FontSize',16,'Fontweight','normal','color','k');
ylabel('Velocity (m/s)','FontSize',16,'Fontweight','normal','color','k');
hold on;
plot(freq,vr1,'g.','LineWidth',4);

