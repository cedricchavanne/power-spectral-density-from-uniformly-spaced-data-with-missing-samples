% import data from Kang et al. (2003) downloaded from
% https://engineering.jhu.edu/meneveau/hot-wire-data.html
% Cédric Chavanne, 2026-01-27.
% cedric.chavanne@ensta.org

close all
clear all
datadir = '../../data/Kang2003/M20H1/';
m = 1200000; % data length for one file
n = 30; % number of files
u = repmat(nan,n*m,1);
dt = 1/40e3;

for j=1:n
    datafile = [datadir,'m20h1-',num2str(j,'%02d'),'.dat'];
    [u1,u2,u3,u4,v1,v2,v3,v4] = textread(datafile,'%n%n%n%n%n%n%n%n');
    u(1+(j-1)*m:j*m) = u1;
end

figure(1)
clf
plot(u)

README = char('Decaying turbulence data from Kang et al. (2003)',...
    'downloaded from https://engineering.jhu.edu/meneveau/hot-wire-data.html',...
    'dt: sampling interval (s)',...
    'u: streamwise velocity for probe 1 (m/s)');

save ../../data/Kang2003.mat dt u README
