function R  = makeParameterJacobians(Q,x)
% makeParameterJacobians
%   pulled this out of makeR.m to speed up code, only do b parameter jacobians
%   at end
%	current b parameters:
%       air density
%       overlap
%       Rayleigh scatter cross section
%       slope
%
% -Usage-
%	R = makeRarameterJacobians(Q,x)
%
% -Inputs-
%	Q retrieval a priori information
%	x retrieved parameters
%
% -Outputs-
%	R Jacobians for b parameters

m = length(Q.zRET);
[SHA,SNA,SH,SN] = forwardModelWV(Q,x);
% data structure
yf = [SHA; SNA; SH; SN];
mdata = length(yf);
m1 = length(SHA);
m2 = 2*m1;
m3 = length(SH);
m4 = 2*m3;
n = 2*m + 10; %3;

% Jacobian for air density
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
Kair = [dSHdnA; dSNdnA; dSHdn; dSNdn];
R.Kair = diag(Kair);

% overlap Jacobian
dSHolapi = (SH - exp(x(end-1))) ./ Q.olap;
dSNolapi = (SN - exp(x(end))) ./ Q.olap;
dSHolapiA = (SHA - x(end-3)) ./ Q.olapA;
dSNolapiA = (SNA - x(end-2)) ./ Q.olapA;
kOlap = [dSHolapiA; dSNolapiA; dSHolapi; dSNolapi];
R.Kolap = diag(kOlap);

% Rayleigh cross section Jacobian
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

% slope Jacobian
KslopeD = (SH-exp(x(end-1))) ./ Q.slope;
KslopeA = (SHA-x(end-3)) ./ Q.slopeA;
KslopeNA = zeros(size(KslopeA));
KslopeN = zeros(size(KslopeD));
R.Kslope = [KslopeA; KslopeNA; KslopeD; KslopeN];

return
