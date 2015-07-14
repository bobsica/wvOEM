function [dSHdxA,dSNdxA,dSHdx,dSNdx] = derivSHSN2(Q,x,j,theFM)
% derivative of forward model with respect to vmr

[SHA, SNA, SH, SN] = theFM(Q,x);

if ~isempty(find(isnan(x)) == 1)
    'after FM: Nans in retrieval vector derivSHSN2'
    stop
end

dn = 1.e-4 .* x(j); % 1e-5
xpert = x;
if x(j) == 0 % trap for tau's where tau(1) = 0
    dn = 1.e-4 .* x(j+1);
end
xpert(j) = x(j) + dn;

[SHAj, SNAj, SHj, SNj] = theFM(Q,xpert);
%dFdRat = (ratj - rat) ./ dn;
dSHdx = (SHj - SH) ./ dn;
dSNdx = (SNj - SN) ./ dn;
dSHdxA = (SHAj - SHA) ./ dn;
dSNdxA = (SNAj - SNA) ./ dn;
%dFdRat = dSHdx./SH - dSNdx./SN;

return

