function [SHA, SNA, SH , SN] = forwardModelWV2d0(Q,x)
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

% retrieved parameters must be on data grid
m = length(Q.zRET);
n = 2*m + 10; % q, tau + 2 back + Angstrom + CN
lambda = 0.3547; lambdaH = 0.40749; lambdaN = 0.38669;

if Q.logWV
    xj = exp(interp1(Q.zRET,x(1:m),Q.zDATAn,'linear'));
    xjA = exp(interp1(Q.zRET,x(1:m),Q.zDATAnA,'linear'));
else
    xj = interp1(Q.zRET,x(1:m),Q.zDATAn,'linear');
    xjA = interp1(Q.zRET,x(1:m),Q.zDATAnA,'linear');
end
if Q.logAlpha
    ODj = exp(interp1(Q.zRET,x(m+1:2*m),Q.zDATAn,'linear'));
    ODjA = exp(interp1(Q.zRET,x(m+1:2*m),Q.zDATAnA,'linear'));
else
    ODj = interp1(Q.zRET,x(m+1:2*m),Q.zDATAn,'linear');
    ODjA = interp1(Q.zRET,x(m+1:2*m),Q.zDATAnA,'linear');
end

tauR = Q.tauRno .* exp(-ODj);
tauH = Q.tauHno .* exp(-(ODj.*(lambdaH./lambda).^-x(end-6)));
tauN = Q.tauNno .* exp(-(ODj.*(lambdaN./lambda).^-x(end-6)));
Q.tauR4J = tauR;
Q.tauH4J = tauH;
Q.tauN4J = tauN;
tauRA = Q.tauRnoA .* exp(-ODjA);
tauHA = Q.tauHnoA .* exp(-(ODjA.*(lambdaH./lambda).^-x(end-6)));
tauNA = Q.tauNnoA .* exp(-(ODjA.*(lambdaN./lambda).^-x(end-6)));

% counts (analog)
CNpA = exp(x(end-8));
CHpA = exp(x(end-9));
SHA = Q.olapA.*CHpA./Q.zDATAnA.^2 .* ((xjA.*(Q.nNA./Q.N2rat))).*tauRA.*tauHA + x(end-3);
SNA = Q.olapA.*CNpA./Q.zDATAnA.^2 .* (Q.nNA.*tauRA.*tauNA) + x(end-2);

% counts (digital)
CNp = exp(x(end-7));
CHp = Q.slope .* CNp; % note slope is in vmr units
SHtrue = Q.olap.*CHp./Q.zDATAn.^2 .* ((xj.*(Q.nN./Q.N2rat))).*tauR.*tauH; % + x(end-1);
SNtrue = Q.olap.*CNp./Q.zDATAn.^2 .* (Q.nN.*tauR.*tauN); % + x(end);

% dead time correction
SHHz = Q.y2Hz .* SHtrue;
SH = SHtrue.*exp(-SHHz.*x(end-5).*1e-9) + x(end-1); 
SNHz = Q.y2Hz .* SNtrue;
SN = SNtrue .*exp(-SNHz.*x(end-4).*1e-9) + x(end);

return
