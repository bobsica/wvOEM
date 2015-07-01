function [dFdSHA,dFdSNA,dFdSH,dFdSN] = derivBN2(Q,x,j,theFM)
% derivative of SN with respect to n_tot

[SHA,SNA,SH,SN] = theFM(Q,x);

dn = 1e-5 .* x(j);
xpert = x;
xpert(j) = x(j) + dn;

[SHAj,SNAj,SHj,SNj] = theFM(Q,xpert);

dFdSHA = (SHAj - SHA) ./ dn;
dFdSNA = (SNAj - SNA) ./ dn;
dFdSH = (SHj - SH) ./ dn;
dFdSN = (SNj - SN) ./ dn;

return

