function [R,yf,J]  = makeJ(Q, R, x, iter)
% makeR - required for OEM to call the forward model and calculate/update the
% Jacobians
%
% -Usage-
%	[R,yf,J]  = makeJ(Q, R, x, iter)
%
% -Inputs-
%	Q - retrieval input data structure
%	R - structure for the Jacobians of model parameters - computed now after
%	retrieval in function R  = makeParameterJacobians(Q,x)
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

J = Kernel;

return

'quack quack'
ssttoppp

% b parameters Air density
KairD = zeros(mdata,n);
KolapD = zeros(mdata,n);
KairA = zeros(mdata,n);
KolapA = zeros(mdata,n);

dSHolapi = (SH - exp(x(end-1))) ./ Q.olap;
dSNolapi = (SN - exp(x(end))) ./ Q.olap;
dSHolapiA = (SHA - x(end-3)) ./ Q.olapA;
dSNolapiA = (SNA - x(end-2)) ./ Q.olapA;

% dSHdz = -(SH - exp(x(end-1))) ./ Q.zDATAn;
% dSNdz = -(SN - exp(x(end))) ./ Q.zDATAn;
% dSHdzA = -(SHA - x(end-3)) ./ Q.zDATAnA;
% dSNdzA = -(SNA - x(end-2)) ./ Q.zDATAnA;
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
%term2 = -((SH-exp(x(end-1))) .* nAir ./ dndz) .* (Q.sigmaR + Q.sigmaH);
term1 = (SH-exp(x(end-1))) ./ nAir;
term2 = -(SH-exp(x(end-1))) .* (Q.sigmaR + Q.sigmaH) .* nAir ./ dndz;
dSHdn = term1 + term2; %zeros(size(nAir)); 
term1 = (SN-exp(x(end))) ./ Q.nN;
term2 = -(SN-exp(x(end))) .* (Q.sigmaR + Q.sigmaN) .* nAir ./ dndz;
dSNdn = term1 + term2; %zeros(size(nAir)); 
term1 = (SHA-x(end-3)) ./ nAirA;
term2 = -(SHA-x(end-3)) .* (Q.sigmaR + Q.sigmaH) .* nAirA ./ dndzA;
dSHdnA = term1 + term2; %zeros(size(nAirA));
term1 = (SNA-x(end-2)) ./ Q.nNA;
term2 = -(SNA-x(end-2)) .* (Q.sigmaR + Q.sigmaN) .* nAirA ./ dndzA;
dSNdnA = term1 + term2;

% dSHdn = (SH-exp(x(end-1))) .* (Q.sigmaR.*nAir);
% term1 = (SN-exp(x(end))) ./ nAir .* Q.N2rat;
% term2 = -((SN-exp(x(end))) .* nAir ./ dndz) .* (Q.sigmaR + Q.sigmaN);
% dSNdn = term1 + term2;
% term1A = (SHA-x(end-3)) ./ nAirA;
% term2A = -((SHA-x(end-3)) .* nAirA ./ dndzA) .* (Q.sigmaR + Q.sigmaH);
% dSHdnA = term1A + term2A;
% term1A = (SNA-x(end-2)) ./ nAirA .* Q.N2rat;
% term2A = -((SNA-x(end-2)) .* nAirA ./ dndzA) .* (Q.sigmaR + Q.sigmaN);
% dSNdnA = term1A + term2A;
% dn_i / dn_j correction to dF/dn_i to go from data to retrieval grid note
% for zRET < or > zDATA you can't interpolate, so am using ifac for linear
% extrapolation
dzDATA = Q.zDATAn(2) - Q.zDATAn(1);
ll = 0;
for l = 1:2
    for jj = ll+1:m+ll
        j = jj - ll;
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
            KairA(i,jj) = dSHdnA(i) .* ifacA;
            KolapA(i,jj) = dSHolapiA(i) .* ifacA;
            KairA(k,jj) = dSNdnA(i) .* ifacA;
            KolapA(k,jj) = dSNolapiA(i) .* ifacA;
            KairD(kk,jj) = dSHdn(i) .* ifac;
            KolapD(kk,jj) = dSHolapi(i) .* ifac;
            KairD(kkk,jj) = dSNdn(i) .* ifac;
            KolapD(kkk,jj) = dSNolapi(i) .* ifac;
        end
    end
    ll = ll + m;
end
R.Kair = KairD + KairA;
R.Kolap = KolapD + KolapA;

% KairA(m1,1:m) = Kair(m1-1,1:m);
% KolapA(m1,1:m) = KolapA(m1-1,1:m);
% KairA(2*m1,1:m) = KairA(2*m1-1,1:m);
% KolapA(2*m1,1:m) = KolapA(2*m1-1,1:m);
% Kair(m2+m3,1:m) = Kair(m2+m3-1,1:m);
% Kolap(m2+m3,1:m) = Kolap(m2+m3-1,1:m);
% Kair(mdata,1:m) = Kair(mdata-1,1:m);
% Kolap(mdata,1:m) = Kolap(mdata-1,1:m);

Kolap = zeros(mdata,n);
Kair = zeros(mdata,n);
for i = 1:mdata
    for jj = 1:2*m
        j = jj;
        if j > m
            j = j - m;
        end
        QQ = Q;
        dn = 1e-4 .* Q.olapRET(j);
        QQ.olapRET(j) = Q.olapRET(j) + dn;
        [SHAj, SNAj, SHj, SNj] = forwardModelWV(QQ,x);
        dSHdx = (SHj - SH) ./ dn;
        dSNdx = (SNj - SN) ./ dn;
        dSHdxA = (SHAj - SHA) ./ dn;
        dSNdxA = (SNAj - SNA) ./ dn;
        Kolap(1:m1,jj) = dSHdxA;
        Kolap(m1+1:m2,jj) = dSNdxA;
        Kolap(m2+1:m2+m3,jj) = dSHdx;
        Kolap(m2+m3+1:mdata,jj) = dSNdx;
        QQ = Q;
        dn = 1e-4 .* Q.nNret(j);
        QQ.nNret(j) = Q.nNret(j) + dn;
        [SHAj, SNAj, SHj, SNj] = forwardModelWV(QQ,x);
        dSHdx = (SHj - SH) ./ dn;
        dSNdx = (SNj - SN) ./ dn;
        dSHdxA = (SHAj - SHA) ./ dn;
        dSNdxA = (SNAj - SNA) ./ dn;
        Kair(1:m1,jj) = dSHdxA;
        Kair(m1+1:m2,jj) = dSNdxA;
        Kair(m2+1:m2+m3,jj) = dSHdx;
        Kair(m2+m3+1:mdata,jj) = dSNdx;       
    end
end
% R.Kair = Kair;
% R.Kolap = Kolap;

% Rayleigh cross section error
if Q.logAlpha
    ODj = exp(interp1(Q.zRET,x(m+1:2*m),Q.zDATAn,'linear'));
    ODjA = exp(interp1(Q.zRET,x(m+1:2*m),Q.zDATAnA,'linear'));
else
    ODj = interp1(Q.zRET,x(m+1:2*m),Q.zDATAn,'linear');
    ODjA = interp1(Q.zRET,x(m+1:2*m),Q.zDATAnA,'linear');
end

dSHAsigmaRay = -((SHA-x(end-3)) .* (ODjA./Q.sigmaR) +...
    (SHA-x(end-3)).*((ODjA.*(Q.lambdaH./Q.lambda).^-x(end-6))./Q.sigmaH));
dSNAsigmaRay = -((SNA-x(end-2)) .* (ODjA./Q.sigmaR) +...
    (SNA-x(end-2)).*((ODjA.*(Q.lambdaN./Q.lambda).^-x(end-6))./Q.sigmaN));
dSHsigmaRay = -((SH-exp(x(end-1))) .* (ODj./Q.sigmaR) +...
    (SH-exp(x(end-1))).*((ODj.*(Q.lambdaH./Q.lambda).^-x(end-6))./Q.sigmaH));
dSNsigmaRay = -((SN-exp(x(end))) .* (ODj./Q.sigmaR) +...
    (SN-exp(x(end))).*((ODj.*(Q.lambdaN./Q.lambda).^-x(end-6))./Q.sigmaN));
R.KsigmaRay = [dSHAsigmaRay; dSNAsigmaRay; dSHsigmaRay; dSNsigmaRay]; 

% dSHdSigmaR = -(SH-exp(x(end-1))) .* Q.tauR ./ Q.sigmaR; 
% % breaking Rayleigh into 2 parts (for SH and SN)
% dSNdSigmaR = -(SN-exp(x(end))) .* Q.tauR ./ Q.sigmaR;
% dSHdSigmaH = - (SH-exp(x(end-1))) .* Q.tauH ./ Q.sigmaH; % SH ./ Q.sigmaH 
% dSNdSigmaN = - (SN-exp(x(end))) .* Q.tauN ./ Q.sigmaN; % SN ./ Q.sigmaN
% dSHdSigmaRA = -(SHA-x(end-3)) .* Q.tauRA ./ Q.sigmaR; 
% % breaking Rayleigh into 2 parts (for SH and SN)
% dSNdSigmaRA = -(SNA-x(end-2)) .* Q.tauRA ./ Q.sigmaR;
% dSHdSigmaHA = - (SHA-x(end-3)) .* Q.tauHA ./ Q.sigmaH; % SH ./ Q.sigmaH 
% dSNdSigmaNA = - (SNA-x(end-2)) .* Q.tauNA ./ Q.sigmaN; % SN ./ Q.sigmaN

% slope Jacobian
KslopeD = (SH-exp(x(end-1))) ./ Q.slope;
KslopeA = (SHA-x(end-3)) ./ Q.slopeA;
KslopeNA = zeros(size(KslopeA));
KslopeN = zeros(size(KslopeD));
R.Kslope = [KslopeA; KslopeNA; KslopeD; KslopeN];

% R.KsigmaSNR = dSHdSigmaR;
% R.KsigmaN = dSNdSigmaN;
% R.KsigmaH = dSHdSigmaH;
% R.KsigmaSHRA = dSHdSigmaRA;
% R.KsigmaSNRA = dSHdSigmaRA;
% R.KsigmaNA = dSNdSigmaNA;
% R.KsigmaHA = dSHdSigmaHA;
%R.KairA = KairA;
%R.KslopeA = KslopeA;
%R.KolapA = KolapA;


if ~isempty(find(isnan(J)) == 1)
    'after FM: Nans in kernel (FMwv(n).m)'
    iter
    dbstop in makeR at 226
end

return
