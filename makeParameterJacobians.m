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

% KairD = zeros(mdata,n);
% KolapD = zeros(mdata,n);
% KairA = zeros(mdata,n);
% KolapA = zeros(mdata,n);
% 
% dSHolapi = (SH - exp(x(end-1))) ./ Q.olap;
% dSNolapi = (SN - exp(x(end))) ./ Q.olap;
% dSHolapiA = (SHA - x(end-3)) ./ Q.olapA;
% dSNolapiA = (SNA - x(end-2)) ./ Q.olapA;

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
Kair = [dSHdnA; dSNdnA; dSHdn; dSNdn];
R.Kair = diag(Kair);

% dzDATA = Q.zDATAn(2) - Q.zDATAn(1);
% ll = 0;
% for l = 1:2
%     for jj = ll+1:m+ll
%         j = jj - ll;
%         for i = 1:m1-1
%             if (Q.zRET(j) >= Q.zDATAn(i)) && (Q.zRET(j) < Q.zDATAn(i+1))
%                 ifac = dzDATA ./ (Q.zDATAn(i+1) - Q.zRET(j));
%             elseif Q.zRET(j) < Q.zDATAn(1)
%                 ifac = dzDATA ./ (Q.zDATAn(i+1) - Q.zRET(j));
%                 % extrapolation: denominator is larger than interpolation since
%                 % it is from the i+1 point
%             elseif Q.zRET(j) > Q.zDATAn(end)
%                 ifac = -dzDATA ./ (Q.zDATAn(i) - Q.zRET(j));
%             else
%                 ifac = 0; 
%             end
%             if (Q.zRET(j) >= Q.zDATAnA(i)) && (Q.zRET(j) < Q.zDATAnA(i+1))
%                 ifacA = dzDATA ./ (Q.zDATAnA(i+1) - Q.zRET(j));
%             elseif Q.zRET(j) < Q.zDATAnA(1)
%                 ifacA = dzDATA ./ (Q.zDATAnA(i+1) - Q.zRET(j));
%                 % extrapolation: denominator is larger than interpolation since
%                 % it is from the i+1 point
%             elseif Q.zRET(j) > Q.zDATAnA(end)
%                 ifacA = -dzDATA ./ (Q.zDATAnA(i) - Q.zRET(j));
%             else
%                 ifacA = 0; 
%             end
%             k = i + m1;
%             kk = i + m2;
%             kkk = i + m2+m3;
%             KairA(i,jj) = dSHdnA(i) .* ifacA;
%             KolapA(i,jj) = dSHolapiA(i) .* ifacA;
%             KairA(k,jj) = dSNdnA(i) .* ifacA;
%             KolapA(k,jj) = dSNolapiA(i) .* ifacA;
%             KairD(kk,jj) = dSHdn(i) .* ifac;
%             KolapD(kk,jj) = dSHolapi(i) .* ifac;
%             KairD(kkk,jj) = dSNdn(i) .* ifac;
%             KolapD(kkk,jj) = dSNolapi(i) .* ifac;
%         end
%     end
%     ll = ll + m;
% end
% R.Kair = KairD + KairA;
% R.Kolap = KolapD + KolapA;

% Kolap = zeros(mdata,n);
% Kair = zeros(mdata,n);
% for jj = 1:2*m
%     j = jj;
%     if j > m
%         j = jj - m;
%     end
%     QQ = Q;
%     dn = 1e-12 .* Q.olapRET(j);
%     QQ.olapRET(j) = Q.olapRET(j) + dn;
%     QQ.olap = interp1(Q.zRET,Q.olapRET,Q.zDATAn,'linear');
%     QQ.olapA = interp1(Q.zRET,Q.olapRET,Q.zDATAnA,'linear');
%     [SHAj, SNAj, SHj, SNj] = forwardModelWV(QQ,x);
%     dSHdx = (SHj - SH) ./ dn;
%     dSNdx = (SNj - SN) ./ dn;
%     dSHdxA = (SHAj - SHA) ./ dn;
%     dSNdxA = (SNAj - SNA) ./ dn;
%     Kolap(1:m1,jj) = dSHdxA;
%     Kolap(m1+1:m2,jj) = dSNdxA;
%     Kolap(m2+1:m2+m3,jj) = dSHdx;
%     Kolap(m2+m3+1:mdata,jj) = dSNdx;
%     QQ = Q;
%     dn = 1e-12 .* Q.nNret(j);
%     QQ.nNret(j) = Q.nNret(j) + dn;
%     QQ.nN = interp1(Q.zRET,Q.nNret,Q.zDATAn,'linear');
%     QQ.nNA = interp1(Q.zRET,Q.nNret,Q.zDATAnA,'linear');
%     [SHAj, SNAj, SHj, SNj] = forwardModelWV(QQ,x);
%     dSHdx = (SHj - SH) ./ dn;
%     dSNdx = (SNj - SN) ./ dn;
%     dSHdxA = (SHAj - SHA) ./ dn;
%     dSNdxA = (SNAj - SNA) ./ dn;
%     Kair(1:m1,jj) = dSHdxA;
%     Kair(m1+1:m2,jj) = dSNdxA;
%     Kair(m2+1:m2+m3,jj) = dSHdx;
%     Kair(m2+m3+1:mdata,jj) = dSNdx;       
% end
%  R.Kair = Kair;
% R.Kolap = Kolap;

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
