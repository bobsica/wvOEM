function [R,yf,J]  = makeR(Q, R, x, iter)
% makeR - required for OEM to call the forward model and calculate/update the
% Jacobians
%
% -Usage-
%	[R,yf,J]  = makeR(Q, R, x, iter)
%
% -Inputs-
%	Q - retrieval input data structure
%	R - structure for the Jacobians of model parameters
%   x - retrieved parameters (see below)
%   iter - required and calculated in oem, the iteration numbers of the solution
%
% -Outputs-
%	R - old retrieval Jacobians in, new out
%	yf - forward model calculated on data grid
%   J - retrieval Jacobians
%
% retrieving 2*m + 10 parameters 
% x(1:m) is q 
% x(m+1:2*m) is OD 
% x(end-9) is analog channel WV lidar constant
% x(end-8) is analog channel N2 lidar constant
% x(end-7) is digital channel N2 lidar constant
% x(end-6) is the Angstrom exponent 
% x(end-5) is dead time for SH 
% x(end-4) is dead time for SN
% x(end-3) is H analog background 
% x(end-2) is N analog background 
% x(end-1) is H digital background 
% x(end) is N digital background

m = length(Q.zRET);
[SHA,SNA,SH,SN] = forwardModelWV(Q,x);

if ~isempty(find(isnan(x)) == 1)
    'after FM: Nans in retrieval vector (FMwv(n).m)'
    iter
    stop
end

% data structure
yf = [SHA; SNA; SH; SN];
mdata = length(yf);
m1 = length(SHA);
m2 = 2*m1;
m3 = length(SH);
m4 = 2*m3;

n = 2*m + 10; %3;
Kernel = zeros(mdata,n); 

% Jacobians for q
for j = 1:2*m %m
    [dSHdxA,dSNdxA,dSHdx,dSNdx] = derivSHSN2(Q,x,j,@forwardModelWV);
    Kernel(1:m1,j) = dSHdxA;
    Kernel(m1+1:m2,j) = dSNdxA;
    Kernel(m2+1:m2+m3,j) = dSHdx;
    Kernel(m2+m3+1:mdata,j) = dSNdx;
end

% CN / CHA
[d1,d2,d3,d4] = derivBN2(Q,x,n-9,@forwardModelWV);
Kernel(1:m1,n-9) = d1;
Kernel(m1+1:m2,n-9) = d2;
[d1,d2,d3,d4] = derivBN2(Q,x,n-8,@forwardModelWV);
Kernel(1:m1,n-8) = d1;
Kernel(m1+1:m2,n-8) = d2;
[d1,d2,d3,d4] = derivBN2(Q,x,n-7,@forwardModelWV);
Kernel(m2+1:m2+m3,n-7) = d3;
Kernel(m2+m3+1:mdata,n-7) = d4;

% angstrom exponent
[d1,d2,d3,d4] = derivBN2(Q,x,n-6,@forwardModelWV);
Kernel(1:m1,n-6) = d1;
Kernel(m1+1:m2,n-6) = d2;
Kernel(m2+1:m2+m3,n-6) = d3;
Kernel(m2+m3+1:mdata,n-6) = d4;

% dead times
[dSHdxA,dSNdxA,dSHdx,dSNdx] = derivSHSN2(Q,x,n-5,@forwardModelWV);
Kernel(m2+1:m2+m3,n-5) = dSHdx;
[dSHdxA,dSNdxA,dSHdx,dSNdx] = derivSHSN2(Q,x,n-4,@forwardModelWV);
Kernel(m2+m3+1:mdata,n-4) = dSNdx;

% backgrounds
[d1,d2,d3,d4] = derivBN2(Q,x,n-3,@forwardModelWV);
Kernel(1:m1,n-3) = d1;
[d1,d2,d3,d4] = derivBN2(Q,x,n-2,@forwardModelWV);
Kernel(m1+1:m2,n-2) = d2;
[d1,d2,d3,d4] = derivBN2(Q,x,n-1,@forwardModelWV);
Kernel(m2+1:m2+m3,n-1) = d3;
[d1,d2,d3,d4] = derivBN2(Q,x,n,@forwardModelWV);
Kernel(m2+m3+1:mdata,n) = d4;


% b parameters Air density
Kair = zeros(mdata,n);
Kolap = zeros(mdata,n);
KairA = zeros(mdata,n);
KolapA = zeros(mdata,n);
dSHdz = -(SH - Q.backH) ./ Q.zDATAn;
dSNdz = -(SN - Q.backN) ./ Q.zDATAn;
dSHdzA = -(SHA - Q.backHA) ./ Q.zDATAnA;
dSNdzA = -(SNA - Q.backNA) ./ Q.zDATAnA;
dOlap = Q.olapD;
dOlapA = Q.olapDA ;
dzdo = 1./dOlap;
ff = find(isinf(dzdo));
dzdo(ff) = 0;
dzdoA = 1./dOlapA;
ff = find(isinf(dzdoA));
dzdoA(ff) = 0;
dSHolapi = dSHdz .* dzdo;
dSNolapi = dSNdz .* dzdo;
dSHolapiA = dSHdzA .* dzdoA;
dSNolapiA = dSNdzA .* dzdoA;

nAir = Q.nN ./ Q.N2rat;
nAirA = Q.nNA ./ Q.N2rat;
warning off
p = polyfit(Q.zDATAn,log(nAir),5); % take derivative of polynomial fit to radiosonde n
pA = polyfit(Q.zDATAnA,log(nAirA),5);
warning on
pd = polyder(p);
pdA = polyder(pA);
dndz = exp(polyval(p,Q.zDATAn)) .* polyval(pd,Q.zDATAn);
dndzA = exp(polyval(p,Q.zDATAnA)) .* polyval(pdA,Q.zDATAnA);
% find dF_i/dN on data grid
term1 = (SH-Q.backH) ./ nAir;
term2 = -((SH-Q.backH) .* nAir ./ dndz) .* (Q.sigmaR + Q.sigmaH);
dSHdn = term1 + term2;
term1 = (SN-Q.backN) ./ nAir .* Q.N2rat;
term2 = -((SN-Q.backN) .* nAir ./ dndz) .* (Q.sigmaR + Q.sigmaN);
dSNdn = term1 + term2;
term1A = (SHA-Q.backHA) ./ nAirA;
term2A = -((SHA-Q.backHA) .* nAirA ./ dndzA) .* (Q.sigmaR + Q.sigmaH);
dSHdnA = term1A + term2A;
term1A = (SNA-Q.backNA) ./ nAirA .* Q.N2rat;
term2A = -((SNA-Q.backNA) .* nAirA ./ dndzA) .* (Q.sigmaR + Q.sigmaN);
dSNdnA = term1A + term2A;
% dn_i / dn_j correction to dF/dn_i to go from data to retrieval grid note
% for zRET < or > zDATA you can't interpolate, so am using ifac for linear
% extrapolation
dzDATA = Q.zDATAn(2) - Q.zDATAn(1);
for j = 1:m
    for i = 1:m1-1
        if (Q.zRET(j) >= Q.zDATAn(i)) && (Q.zRET(j) < Q.zDATAn(i+1))
            ifac = dzDATA ./ (Q.zDATAn(i+1) - Q.zRET(j));
        elseif Q.zRET(j) < Q.zDATAn(1)
            ifac = dzDATA ./ (Q.zDATAn(i+1) - Q.zRET(j));
            % extrapolation: denominator is larger than interpolation since
            % it is from the i+1 point
        elseif Q.zRET(j) > Q.zDATAn(end)
            ifac = -dzDATA ./ (Q.zDATAn(i) - Q.zRET(j));
        else
            ifac = 0; 
        end
        if (Q.zRET(j) >= Q.zDATAnA(i)) && (Q.zRET(j) < Q.zDATAnA(i+1))
            ifacA = dzDATA ./ (Q.zDATAnA(i+1) - Q.zRET(j));
        elseif Q.zRET(j) < Q.zDATAnA(1)
            ifacA = dzDATA ./ (Q.zDATAnA(i+1) - Q.zRET(j));
            % extrapolation: denominator is larger than interpolation since
            % it is from the i+1 point
        elseif Q.zRET(j) > Q.zDATAnA(end)
            ifacA = -dzDATA ./ (Q.zDATAnA(i) - Q.zRET(j));
        else
            ifacA = 0; 
        end
        k = i + m1;
        kk = i + m2;
        kkk = i + m2+m3;
        KairA(i,j) = dSHdnA(i) .* ifacA;
        KolapA(i,j) = dSHolapiA(i) .* ifacA;
        KairA(k,j) = dSNdnA(i) .* ifacA;
        KolapA(k,j) = dSNolapiA(i) .* ifacA;
        Kair(kk,j) = dSHdn(i) .* ifac;
        Kolap(kk,j) = dSHolapi(i) .* ifac;
        Kair(kkk,j) = dSNdn(i) .* ifac;
        Kolap(kkk,j) = dSNolapi(i) .* ifac;
    end
end
KairA(m1,1:m) = Kair(m1-1,1:m);
KolapA(m1,1:m) = KolapA(m1-1,1:m);
KairA(2*m1,1:m) = KairA(2*m1-1,1:m);
KolapA(2*m1,1:m) = KolapA(2*m1-1,1:m);
Kair(m2+m3,1:m) = Kair(m2+m3-1,1:m);
Kolap(m2+m3,1:m) = Kolap(m2+m3-1,1:m);
Kair(mdata,1:m) = Kair(mdata-1,1:m);
Kolap(mdata,1:m) = Kolap(mdata-1,1:m);

Kslope = (SH-Q.backH) ./ Q.slope;
KslopeA = (SHA-Q.backHA) ./ Q.slopeA;

dSHdSigmaR = -(SH-Q.backH) .* Q.tauR ./ Q.sigmaR; 
% breaking Rayleigh into 2 parts (for SH and SN)
dSNdSigmaR = -(SN-Q.backN) .* Q.tauR ./ Q.sigmaR;
dSHdSigmaH = - (SH-Q.backH) .* Q.tauH ./ Q.sigmaH; % SH ./ Q.sigmaH 
dSNdSigmaN = - (SN-Q.backN) .* Q.tauN ./ Q.sigmaN; % SN ./ Q.sigmaN
dSHdSigmaRA = -(SHA-Q.backHA) .* Q.tauRA ./ Q.sigmaR; 
% breaking Rayleigh into 2 parts (for SH and SN)
dSNdSigmaRA = -(SNA-Q.backNA) .* Q.tauRA ./ Q.sigmaR;
dSHdSigmaHA = - (SHA-Q.backHA) .* Q.tauHA ./ Q.sigmaH; % SH ./ Q.sigmaH 
dSNdSigmaNA = - (SNA-Q.backNA) .* Q.tauNA ./ Q.sigmaN; % SN ./ Q.sigmaN

R.KsigmaSHR = dSHdSigmaR;
R.KsigmaSNR = dSHdSigmaR;
R.KsigmaN = dSNdSigmaN;
R.KsigmaH = dSHdSigmaH;
R.Kair = Kair;
R.Kslope = Kslope;
R.Kolap = Kolap;
R.KsigmaSHRA = dSHdSigmaRA;
R.KsigmaSNRA = dSHdSigmaRA;
R.KsigmaNA = dSNdSigmaNA;
R.KsigmaHA = dSHdSigmaHA;
R.KairA = KairA;
R.KslopeA = KslopeA;
R.KolapA = KolapA;

J = Kernel;

if ~isempty(find(isnan(J)) == 1)
    'after FM: Nans in kernel (FMwv(n).m)'
    iter
    dbstop in makeR at 226
end

return
