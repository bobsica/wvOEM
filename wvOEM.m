% wvOEM
% Optimal Estimation Method applied to water vapour lidar
% currently specifically for RALMO
% R. Sica
% this version started on github 2 July 2015, with v1.0.0 as of 3 July 2015
%
% v1.0.1: flag to control how background is calculated.
% v1.1.1: backgrounds are calculated from uncorrected data; also measurement
% variance is allowed to be < background variance. New flag for which background
% variance to use for analog as well as digital channel.

VERSION = '1-1-0'
%0305 12: 3000, 1200, 20, 80, 1750, true, true
%0305 00: 5000, 1300, 20, 80, 1600, true, true
date = 20090905; %20150305;
nb = '12';
dextsp = [nb '30'];
oemStop = 6000; %14000; % 0308 8000/11000; 0305 2500, 5000; 200906 - 14000 200905 6000
oemStopA = 3000; %5000; % 0308 5000; 0305 1200, 1300; 200906 - 5000, 200905 3000 day
in.LRfree = 50; % was 20 on 0305, 0308 50, 200905-6 50
in.LRpbl = 80; 
in.LRtranHeight = 1500; %2500; % m, height the above 2 hand off to each other
% 1800 0308, 900/800 0305; 200906 - 1500, 200905 1800
varAV = true;
varAVA = false; % was f
% true use variance of average (night) or variance of measurements (day)

dataPath = '/Users/BobSica/Dropbox/matlab/matlabWork/fromMCH/ralmodata/';
outPath = '/Users/BobSica/Dropbox/matlab/matlabWork/fromMCH/ralmoOEMwvOutput/';
diaryFile = [outPath 'diary' int2str(date) dextsp 'LT-v' VERSION '.markdown'];
if exist(diaryFile) ~= 0
    delete(diaryFile);
end 
diary(diaryFile)
dext = [nb '30chan2.mat']; % extension for data file
dextout = [nb '30chan2-v' VERSION '.mat']; % extension for output file with version
dexts = [nb '00.mat'];
fext = [nb '30chan2.fig'];
fextout = [nb '30chan2-v' VERSION '.fig'];
gext = [nb '30combined.mat'];
oemGo = 50; %50; %50;
oemGoA = 100; %300;
pieceWise = true;
deadTimeH = 4; %4.0; % ns
dfacDeadH = 0.01; % 0308 0.01, 0305 0.001
dfacDeadN = dfacDeadH; %0.01;
deadTimeN = 4;
normStop = 25000;
aposteriori = false; % if true use a posteriori analog variance
% varMask = 0; % variance increase
% maskLow = 100000; %500; % height of low mask
% maskHigh = 0; % height of high mask
asl = 490; %really 491, but some files use 490; % Payerne

reality = true;
if reality
    dataPath = '/Users/BobSica/Dropbox/matlab/matlabWork/fromMCH/ralmodata/';
    ralmoFile = [dataPath 'ralmo' int2str(date) gext];
    load(ralmoFile);
end
savefigs = true
savedat = true
logAlpha = false;
logWV = true;
cf = 3; %tent function on covariance
mAir = 28.966;
mWV = 18.02;
corrLh = 84; %100;
corrLalpha = 84; %100 ; % 500 100 2000
oemPath = './'; % for saving plots

% inputs for ralmo data
in.pieceWise = pieceWise;
in.date = date;
in.slope = 37.88;
in.slopeA = in.slope ./ 3; % 3 is nominal, not accurate 2.75; 
% ad hoc factor 3.1096; %./ 2.569 (night); ./ 3.1096 (day)
% .* (2.76/2*0.9); % ad hoc correction until I get the correct number;
% units are g/kg, slope is corrected to vmr in makeRealWVlog
in.coRET = 2; % 2, was 3, this coadds retrieval grid
in.coAddData = 6; % 6
in.zgo = oemGo;% 59; 0 means start where overlap > 0.1
%60; % 2500 %1000; %75; % 240; 100; %60;
in.zgoA = oemGoA;
in.zstop = normStop; % this is where you cut data for normalization %11e3;
in.zOEM = oemStop; % 8000; this is where you cut at end of makeQ for retrieval
in.zOEMA = oemStopA;
zCNnorm = 3000;
in.zCNnorm = min(zCNnorm,oemStopA);
in.beamDiv = 0.09e-3; % mrad
in.logAlpha = logAlpha; % required to set switch in FM in makeQ9
in.aposteriori = aposteriori;
% in.varMask = varMask;
% in.maskLow = maskLow;
% in.maskHigh = maskHigh;
in.Aoffset = 0; %1.1133e4;
in.asl = asl;
in.varAV = varAV;
in.varAVA = varAVA;
in.dexts = dexts;
in.dextsp = dextsp;
in.dext = dext

% initialize R, retrieval structure, for Q.pack, though we don't use this
R = []; R.jq = {}; R.ji = {}; iter = 0;

% O, input structure
O = makeO;

% Q structure
[Q,y,yvar] = makeQ(in);
Q.logAlpha = logAlpha;
Q.logWV = logWV;
Q.DeadTimeH = deadTimeH;
Q.DeadTimeN = deadTimeN;

m = length(Q.zRET);
%n = 2*m + 10; % retrieving 5 overlap parameters and 2 backgrounds;
n = 2*m + 10; 
% backgrounds (4) + angstrom exponent + C_SN + C_SNA + dead time (N2 & WV)

mchanA = length(Q.zDATAnA);
mchanD = length(Q.zDATAn);
mdata = length(y);

if logAlpha
    if logWV
        x1 = log(Q.qvTrueRET);
        x2 = log(Q.odRret);
    else
        x1 = Q.qvTrueRET;
        x2 = log(Q.odRret);
    end
else
    if logWV
        x1 = log(Q.qvTrueRET);
        x2 = Q.odRret;
    else
        x1 = Q.qvTrueRET;
        x2 = Q.odRret;
    end
end

x_a = [x1; x2; log(Q.CHpA); log(Q.CNpA); log(Q.CNp); Q.Ang; Q.DeadTimeH; Q.DeadTimeN;...
    Q.backHA; Q.backNA; log(Q.backH); log(Q.backN)];

% Covariances
% yvar(1:2*mchan) = 1000*yvar(1:2*mchan);
% 'ad hoc 1000x analog variance increase'
Se = diag(yvar);
dfacq = 0.5; % percent error in vmr
%dfacCN = 0.01;
dfacAir = 0.01;
dfacSigmaN = 0.003; % ISSI recommend
dfacSigmaH = 0.003; % ISSI recommend
dfacSigmaR = 0.003; % ISSI recommend
dfacSigRamH = 0.1; 
% from Inaba and Kobayasi, but these would be the slope variance if known
dfacSigRamN = 0.1;
dfacDiv = 0.1; %0.1;
%dfacCoefs = 0.1; % 05;
dfacTheta = 0.1; %0.1;
dfacSlope = 0.05;
dfacAlpha = 0.1; %0.5;
dfacOD = 0.25; %0.20;
dfacAng = 0.01; % 0.05
dfacCNp = 0.1; % was 0.25;
dfacCHp = 0.5; 
dfacOlap = 0.01;

if logWV
    varlogq = (dfacq .* ones(size(Q.zRET))).^2;
    % since log(vmr) you don't multiple by Q.nHretvmr
else
    varlogq = (dfacq .* Q.qvTrueRET).^2;
end
if logAlpha
    varAlpha = (dfacAlpha .* ones(size(Q.zRET))).^2;
    varOD = (dfacOD .* ones(size(Q.zRET))).^2;
else
    varAlpha1 = (dfacAlpha .* Q.alphaRret).^2;
%    'varAlpha2 set to 0 in makeQ11.m'
    varAlpha2 = (Q.asrCorErrRet).^2;
    varAlpha = varAlpha1 + varAlpha2;
    % 2nd term is uncertainty of baseline for "0" aerosols 
    varOD = (dfacOD .* Q.odRret).^2;
    varOD(1) = varOD(2);
    % can't have variance of 1st value be 0, but Q.odRret(1) = 0
end 

varOla = (dfacOlap.*Q.olap).^2;
varOlap = [varOla; varOla]; % for 2*m retrieval parameters for Solap
varCNp = dfacCNp.^2; %(dfacCNp .* Q.CNp).^2;
varCHp = dfacCHp.^2; %(dfacCNp .* Q.CNp).^2;
varAng = (dfacAng .* Q.Ang).^2;
varDiv = (dfacDiv*Q.beamDiv).^2;
%varcoefs = (dfacCoefs.*Q.coefs).^2;
varTheta = (dfacTheta.*Q.Theta).^2;
varDTH = (dfacDeadH.*Q.DeadTimeH).^2;
varDTN = (dfacDeadN.*Q.DeadTimeN).^2;
varHdig = Q.backVarH./(Q.backH.^2);
varNdig = Q.backVarN./(Q.backN.^2);

%'reducing background variance by 0.1/0.1' vars2 = [varlogq; varAlpha;
%varDT; Q.backVarH; Q.backVarN]; vars2 = [varlogq; varAlpha; Q.backVarH;
%Q.backVarN];
vars2 = [varlogq; varOD; varCHp; varCNp; varCNp; varAng; varDTH; varDTN;...
     Q.backVarHA; Q.backVarNA; varHdig; varNdig];

dzRET = Q.zRET(2) - Q.zRET(1);
lc = (round(corrLh ./ dzRET) .* dzRET) .* ones(size(Q.zRET));
lcalpha = (round(corrLalpha ./ dzRET) .* dzRET) .* ones(size(Q.zRET));

Airvar = (dfacAir .* Q.nNret).^2; % only perturbing N2, not air

S_a = zeros(n,n);
Sair = zeros(n,n);
Solap = zeros(n,n);
SolapA = zeros(n,n);

SsigmaN = (dfacSigmaN.*Q.sigmaN).^2;
SsigmaH = (dfacSigmaH.*Q.sigmaH).^2;
SsigmaR = (dfacSigmaR.*Q.sigmaR).^2;
%SsigRamH = (dfacSigRamH.*Q.sigRamH).^2; SsigRamN =
%(dfacSigRamN.*Q.sigRamN).^2;
Sslope = (dfacSlope.*Q.slope).^2;
%S_DT = varDT;

% mvarL = find(Q.zDATAn < in.maskLow);
% mvarH = find(Q.zDATAn > in.maskHigh);
% mvarLR = find(Q.zRET < in.maskLow);
% mvarHR = find(Q.zRET > in.maskHigh);
for i=1:m
    for j=1:m
        sigprod = sqrt(vars2(i) .* vars2(j));
        sigAirvar = sqrt(Airvar(i) .* Airvar(j));
        sigOlapvar = sqrt(varOlap(i) .* varOlap(j));
        diffz = Q.zRET(i) - Q.zRET(j);
        sumlc = lc(i) + lc(j);
        shape(1) = exp(-(2.*diffz./sumlc).^2); %Gaussian correlation function
        shape(2) = exp(-2.*abs(diffz)./sumlc); %Exponential correlation function
        shape(3) = (1-(1-exp(-1)).*2.*abs(diffz)./sumlc); % Tent correlation func
        if shape(3) < 0
            shape(3) = 0;
        end
        S_a(i,j) = sigprod .* shape(cf); %Gaussian correlation function
        Sair(i,j) = sigAirvar .* shape(cf);
        Solap(i,j) = sigOlapvar .* shape(cf);
%         SolapA(i,j) = Solap(i,j);
%         if Q.zRET(j) < Q.zRET(mvarLR(end)+1) || Q.zRET(j) > Q.zRET(mvarHR(1))
%             SolapA(i,j) = 0;
%         end
    end
%     if Q.zRET(i) < Q.zRET(mvarLR(end)+1) || Q.zRET(i) > Q.zRET(mvarHR(1))
%         SolapA(i,j) = 0;
%     end
end
% using different correlation length for alpha
for i=m+1:2*m
    for j=m+1:2*m
        sigprod = sqrt(vars2(i) .* vars2(j));
        sigOlapvar = sqrt(varOlap(i) .* varOlap(j));
        diffz = Q.zRET(i-m) - Q.zRET(j-m);
        sumlc = lcalpha(i-m) + lcalpha(j-m);
        shape(1) = exp(-(2.*diffz./sumlc).^2); %Gaussian correlation function
        shape(2) = exp(-2.*abs(diffz)./sumlc); %Exponential correlation function
        shape(3) = (1-(1-exp(-1)).*2.*abs(diffz)./sumlc); % Tent correlation func
        if shape(3) < 0
            shape(3) = 0;
        end
        S_a(i,j) = sigprod .* shape(cf); %Gaussian correlation function
        Solap(i,j) = sigOlapvar .* shape(cf);
%         SolapA(i,j) = Solap(i,j);
%         if Q.zRET(j-m) < Q.zRET(mvarLR(end)+1) || Q.zRET(j-m) > Q.zRET(mvarHR(1))
%             SolapA(i,j) = 0;
%         end
    end
%     if Q.zRET(i-m) < Q.zRET(mvarLR(end)+1) || Q.zRET(i-m) > Q.zRET(mvarHR(1))
%         SolapA(i,j) = 0;
%     end
end

S_a(n-9,n-9) = vars2(n-9);
S_a(n-8,n-8) = vars2(n-8);
S_a(n-7,n-7) = vars2(n-7);
S_a(n-6,n-6) = vars2(n-6);
S_a(n-5,n-5) = vars2(n-5);
S_a(n-4,n-4) = vars2(n-4);
S_a(n-3,n-3) = vars2(n-3);
S_a(n-2,n-2) = vars2(n-2);
S_a(n-1,n-1) = vars2(n-1);
S_a(n,n) = vars2(n);

S_ainv = [];
Seinv = [];

[X,R] = oem(O,Q,R,@makeR,S_a,Se,S_ainv,Seinv,x_a,y);
if logWV
    X.vmr = exp(X.x(1:m));
else
    X.vmr = X.x(1:m);
end
%X.alpha = exp(X.x(m+1:2*m));
'X.cost'
X.cost

degF = trace(X.A(1:m,1:m));
fini = round(degF);

if logAlpha
    xaAlpha = exp(x_a(m+1:2*m) + Q.odnormR);
    xAlpha = exp(X.x(m+1:2*m) + Q.odnormR);
else
    xaAlpha = x_a(m+1:2*m); % done in makeQ + Q.odnormR (from ground to start);
    xAlpha = X.x(m+1:2*m); % + Q.odnormR;
end

disp(' ')
outB =[Q.backHA X.x(end-3) (X.x(end-3)-Q.backHA)/Q.backHA*100];
str = ['Change in H2O analog background term from ', num2str(outB(1),'%0.5g'),...
    ' to ',num2str(outB(2),'%0.5g'),' +/- ',num2str(X.e(end-3),'%0.5g')];
disp(str)

disp(' ')
outB =[Q.backNA X.x(end-2) (X.x(end-2)-Q.backNA)/Q.backNA*100];
str = ['Change in N2 analog background term from ', num2str(outB(1),'%0.5g'),...
    ' to ',num2str(outB(2),'%0.5g'),' +/- ',num2str(X.e(end-2),'%0.5g')];
disp(str)

disp(' ')
outB =[Q.backH exp(X.x(end-1)) (exp(X.x(end-1))-Q.backH)/Q.backH*100];
str = ['Change in H2O digital background term from ', num2str(outB(1),'%0.5g'),...
    ' to ',num2str(outB(2),'%0.5g'),' +/- ',num2str(X.e(end-1)*100,'%0.5g'), ' %'];
disp(str)

disp(' ')
outB =[Q.backN exp(X.x(end)) (exp(X.x(end))-Q.backN)/Q.backN*100];
str = ['Change in N2 digital background term from ', num2str(outB(1),'%0.5g'),...
    ' to ',num2str(outB(2),'%0.5g'),' +/- ',num2str(X.e(end)*100,'%0.5g'),' %'];
disp(str)

disp(' ')
outB =[Q.CNpA exp(X.x(end-8)) (exp(X.x(end-8))-Q.CNpA)/Q.CNpA*100];
str = ['Change in N2 analog channel lidar constant from ', num2str(outB(1),...
    '%0.5g'), ' to ',num2str(outB(2),'%0.5g'),' +/- ',...
    num2str(X.e(end-8).*exp(X.x(end-8)),'%0.5g')];
disp(str)

disp(' ')
outB =[Q.CHpA exp(X.x(end-9)) (exp(X.x(end-9))-Q.CHpA)/Q.CHpA*100];
str = ['Change in WV analog channel lidar constant from ', num2str(outB(1),...
    '%0.5g'), ' to ',num2str(outB(2),'%0.5g'),' +/- ',...
    num2str(X.e(end-9).*exp(X.x(end-9)),'%0.5g')];
disp(str)

disp(' ')
outB =[Q.CNp exp(X.x(end-7)) (exp(X.x(end-3))-Q.CNp)/Q.CNp*100];
str = ['Change in N2 lidar digital channel constant from ',...
    num2str(outB(1),'%0.5g'), ' to ',num2str(outB(2),'%0.5g'),' +/- ',...
    num2str(X.e(end-7).*exp(X.x(end-7)),'%0.5g')];
disp(str)

outSlopeA = exp(X.x(end-9)) ./ exp(X.x(end-8));

disp(' ')
outB =[Q.Ang X.x(end-6) (X.x(end-6)-Q.Ang)/Q.Ang*100];
str = ['Change in Angstrom Exponent from ', num2str(outB(1),'%0.5g'), ' to ',...
    num2str(outB(2),'%0.5g'),' +/- ',num2str(X.e(end-6),'%0.5g')];
disp(str)

disp(' ')
outB =[Q.DeadTimeH X.x(end-5) (X.x(end-5)-Q.DeadTimeH)/Q.DeadTimeH*100];
str = ['Change in H2O dead time from ', num2str(outB(1),'%0.5g'), ' to ',...
    num2str(outB(2),'%0.5g'),' +/- ',num2str(X.e(end-5),'%0.5g'),' ns'];
disp(str)

disp(' ')
outB =[Q.DeadTimeN X.x(end-4) (X.x(end-4)-Q.DeadTimeN)/Q.DeadTimeN*100];
str = ['Change in N2 dead time from ', num2str(outB(1),'%0.5g'), ' to ',...
    num2str(outB(2),'%0.5g'),' +/- ',num2str(X.e(end-4),'%0.5g'),' ns'];
disp(str)

% remove retrieval points outside of data grid
mbLo = find(Q.zRET < Q.zDATAn(1));
mbHi = find(Q.zRET > Q.zDATAn(end));
mbH = mbHi(1) - 1;
mbL = mbLo(end) + 1;
mbH2 = m - mbH;

% plot Jacobians and averaging kernels
handfig(1) = figure;
subplot(2,2,1)
plot(X.J(1:mchanA,mbL:mbH),Q.zDATAnA./1000)
if logWV
    title('ln(q) - Analog')
else
    title('q - Analog')
end
xlabel('Jacobian')
ylabel('Altitude (km)')

subplot(2,2,2)
%plot(X.J(mchan+1:mdata,(m+mbL:2*m-mbH2)),Q.zDATAn./1000)
plot(X.J(mchanA+1:2*mchanA,(m+3:2*m-mbH2)),Q.zDATAnA./1000)
if logAlpha
    title('ln(opticalDepth) - Analog')
else
    title('opticalDepth - Analog')
end

xlabel('Jacobian')
ylabel('Altitude (km)')

subplot(2,2,3)
plot(X.J(2*mchanA+1:2*mchanA+mchanD,mbL:mbH),Q.zDATAn./1000)
if logWV
    title('ln(q) - Digital')
else
    title('q - Digital')
end
xlabel('Jacobian')
ylabel('Altitude (km)')

subplot(2,2,4)
%plot(X.J(mchan+1:mdata,(m+mbL:2*m-mbH2)),Q.zDATAn./1000)
plot(X.J(2*mchanA+mchanD+1:mdata,(m+3:2*m-mbH2)),Q.zDATAn./1000)
if logAlpha
    title('ln(opticalDepth) - Digital')
else
    title('opticalDepth - Digital')
end

xlabel('Jacobian')
ylabel('Altitude (km)')

handfig(2) = figure;
subplot(1,2,1)
plot(X.A(mbL:mbH,mbL:mbH),Q.zRET(mbL:mbH)./1000,':')
hold on
unit = ones(size(Q.zRET(mbL:mbH)));
response = X.A(mbL:mbH,mbL:mbH)*unit;
plot(response,Q.zRET(mbL:mbH)./1000,'r:')
if logWV
    title('ln(q)')
else
    title('q')
end
xlabel('Averaging Kernels')
ylabel('Altitude (km)')
%axis([-0.1 1.1 Q.zRET(1)./1000 Q.zRET(end)./1000])
xlim([-0.1 1.1])
pltx = get(gca,'XLim');
plot(pltx,[Q.zRET(fini) Q.zRET(fini)]./1000,'k--')

subplot(1,2,2)
plot(X.A((m+mbL+2:2*m-mbH2),(m+mbL+2:2*m-mbH2)),Q.zRET(mbL+2:mbH)./1000,':')
hold on
unit = ones(size(Q.zRET(mbL+2:mbH)));
response = X.A((m+mbL+2:2*m-mbH2),(m+mbL+2:2*m-mbH2))*unit;
plot(response,Q.zRET(mbL+2:mbH)./1000,'r:')
if logAlpha
    title('ln(opticalDepth)')
else
    title('opticalDepth')
end
xlabel('Averaging Kernels')
ylabel('Altitude (km)')
%axis([-0.1 1.1 Q.zRET(1)./1000 Q.zRET(end)./1000])
xlim([-0.1 1.1])
pltx = get(gca,'XLim');
plot(pltx,[Q.zRET(fini) Q.zRET(fini)]./1000,'k--')

% calculate and plot vertical resolution; skip 1st 2 and last kernels
wwidth = zeros(size(Q.zRET));
k = 0;
for j = 1:m % 3:m-1
    k = k + 1;
    wwidth(k) = fwhmquiet(Q.zRET,X.A(1:m,j));
    if isnan(wwidth(k))
        wwidth(k) = 0;
    elseif wwidth(k) == -999
        wwidth(k) = 0;
    end
end
width = wwidth(mbL:mbH);
handfig(3) = figure;
plot(width,Q.zRET(mbL:mbH)./1000)
hold on
xlabel('Vertical Resolution (m)')
ylabel('Alitude (km)')
%axis([0 300 Q.zRET(1)./1000 Q.zRET(end)./1000])
pltx = get(gca,'XLim');
plot(pltx,[Q.zRET(fini) Q.zRET(fini)]./1000,'k--')

% plot count residuals
handfig(4) = figure;
subplot(2,2,1)
plot((X.yf(1:mchanA)-y(1:mchanA))./y(1:mchanA).*100,Q.zDATAnA./1000,'b')
hold on
plot(sqrt(yvar(1:mchanA))./y(1:mchanA)*100,Q.zDATAnA./1000,'r')
plot(-sqrt(yvar(1:mchanA))./y(1:mchanA)*100,Q.zDATAnA./1000,'r')
% pltx = get(gca,'XLim');
% plot(pltx,[maskLow maskLow]./1000,'k--')
% plot(pltx,[maskHigh maskHigh]./1000,'k--')
xlabel('Percent difference of Analog SH')
ylabel('Altitude (km)')
%xlim([-30 30])
%ylim([0 1.05.*oemStop./1000])

subplot(2,2,2)
plot((X.yf(mchanA+1:2*mchanA)-y(mchanA+1:2*mchanA))./y(mchanA+1:2*mchanA).*100,...
    Q.zDATAnA./1000,'b')
hold on
plot(sqrt(yvar(mchanA+1:2*mchanA))./y(mchanA+1:2*mchanA)*100,Q.zDATAnA./1000,'r')
plot(-sqrt(yvar(mchanA+1:2*mchanA))./y(mchanA+1:2*mchanA)*100,Q.zDATAnA./1000,'r')
% pltx = get(gca,'XLim');
% plot(pltx,[maskLow maskLow]./1000,'k--')
% plot(pltx,[maskHigh maskHigh]./1000,'k--')
xlabel('Percent difference of Analog SN')
ylabel('Altitude (km)')
%ylim([0 1.05.*oemStop./1000])

subplot(2,2,3)
plot((X.yf(2*mchanA+1:2*mchanA+mchanD)-y(2*mchanA+1:2*mchanA+mchanD))...
    ./y(2*mchanA+1:2*mchanA+mchanD).*100,Q.zDATAn./1000,'b')
hold on
plot(sqrt(Q.yTrueH)./Q.yTrueH*100,Q.zDATAn./1000,'r')
%legend('Residuals','Poisson Noise','Location','Best')
plot(-sqrt(Q.yTrueH)./Q.yTrueH*100,Q.zDATAn./1000,'r')
plot(sqrt(yvar(2*mchanA+1:2*mchanA+mchanD))./Q.yTrueH*100,Q.zDATAn./1000,'gx')
plot(-sqrt(yvar(2*mchanA+1:2*mchanA+mchanD))./Q.yTrueH*100,Q.zDATAn./1000,'gx')
xlabel('Percent difference of SH')
ylabel('Altitude (km)')
%xlim([-100 100])
ylim([0 1.05.*oemStop./1000])

subplot(2,2,4)
%fd1 = find(Q.zDATAn > oemGoA);
%fdh = find(Q.zDATAn > oemStopA);
plot((X.yf(2*mchanA+mchanD+1:mdata)-y(2*mchanA+mchanD+1:mdata))./...
    y(2*mchanA+mchanD+1:mdata).*100,Q.zDATAn./1000,'b')
% plot((X.yf(2*mchanA+mchanD+1:2*mchanA+mchanD+fdh(1)-1).*...
%(in.slope./in.slopeA)-y(2*mchanA+mchanD+1:2*mchanA+mchanD+fdh(1)-1))./...
%     y(2*mchanA+mchanD+1:2*mchanA+mchanD+fdh(1)-1).*100,Q.zDATAn(1:fdh(1)-1)./1000,'b')
% plot((X.yf(2*mchanA+mchanD+1+fdh(1):end)-y(2*mchanA+mchanD+1+fdh(1):end))./...
%     y(2*mchanA+mchanD+1+fdh(1):end).*100,Q.zDATAn(fdh(1)+1:end)./1000,'b')
hold on
plot(sqrt(Q.yTrueN)./Q.yTrueN*100,Q.zDATAn./1000,'r')
%legend('Residuals','Poisson Noise','Location','Best')
plot(-sqrt(Q.yTrueN)./Q.yTrueN*100,Q.zDATAn./1000,'r')
plot(sqrt(yvar(2*mchanA+mchanD+1:mdata))./y(2*mchanA+mchanD+1:mdata)...
    *100,Q.zDATAn./1000,'gx')
plot(-sqrt(yvar(2*mchanA+mchanD+1:mdata))./y(2*mchanA+mchanD+1:mdata)...
    *100,Q.zDATAn./1000,'gx')
xlabel('Percent difference of SN')
ylabel('Altitude (km)')
xlim([-20 20])
ylim([0 1.05.*oemStop./1000])

% plot retrievals of n and q
handfig(5) = figure;
if reality
%    lambda = 0.3547; lambdaH = 0.40749; lambdaN = 0.38669; tH = xAlpha .*
%    (lambda./lambdaH).^(-X.x(end-2)); tN = xAlpha .*
%    (lambda./lambdaN).^(-X.x(end-2)); TrH =
%    interp1(Q.zRET,exp(-tH),Q.zDATAn,'linear'); TrN =
%    interp1(Q.zRET,exp(-tN),Q.zDATAn,'linear'); wvTret = in.slope .*
%    (TrH./TrN) .* (Q.tauHno./Q.tauNno) .* ((Q.SHcoadd -
%    Q.backH)./(Q.SNcoadd - Q.backN));
%    semilogx(out.q,(out.z-490)./1000,'bx') % RALMO standard plot
    nRalmo = length(ralmo.t) ./ ralmo.records;
    itst = 1:1:length(ralmo.t);
    mR = mod(itst,nRalmo);
    fmR = find(mR ==0);
    rTime = ralmo.t(fmR);
    fRalmo = find(rTime > Q.ralmoTimeEnd);
    iRalmo = fRalmo(1);
    'RALMO time ='
    datevec(rTime(iRalmo))
    'Analysis time ='
    datevec(Q.ralmoTimeEnd)
    gg = 1 + (iRalmo-1)*nRalmo;
    semilogx(ralmo.q(gg:nRalmo+gg-1),(ralmo.z(gg:nRalmo+gg-1)-490)./1000,'bx')
    hold on
else
%    semilogx(Q.qvTrueRET.*1000.*(mWV./mAir),Q.zRET./1000); % was ./(mAir
%    ./ mWV)
    hold on
end
fd1 = find(Q.zRET > oemGoA);
fdh = find(Q.zRET > oemStopA);
% % % comment lines are for plots with step changes
% % % this one swapped in.slopeA & in.slope
% % semilogx(X.vmr(mbL:fd1(1)-1).*1000.*(mAir./mWV).*(in.slope./in.slopeA)...
% %     ,Q.zRET(mbL:fd1(1)-1)./1000,'r','LineWidth',2)
% % % was this at night? (scaling reversed)
% % % semilogx(X.vmr(fd1(1):fdh(1)-1).*1000.*(mAir./mWV).*(in.slopeA./in.slope)...
% % %     ,Q.zRET(fd1(1):fdh(1)-1)./1000,'r','LineWidth',2)
% % semilogx(X.vmr(fd1(1):fdh(1)-1).*1000.*(mAir./mWV).*(in.slope./in.slopeA)...
% %     ,Q.zRET(fd1(1):fdh(1)-1)./1000,'r','LineWidth',2)
% % semilogx(X.vmr(fdh(1)+2:mbH).*1000.*(mAir./mWV).*(in.slope./in.slopeA)...
% %     ,Q.zRET(fdh(1)+2:mbH)./1000,'r','LineWidth',2)
% % %semilogx(X.vmr(mbL:mbH).*1000.*(mAir./mWV).*(in.slope./in.slopeA)...
% % %    ,Q.zRET(mbL:mbH)./1000,'r','LineWidth',2)
semilogx(X.vmr(mbL:mbH).*1000.*(mWV./mAir).*(outSlopeA./in.slope),...
    Q.zRET(mbL:mbH)./1000,'r','LineWidth',2)
if logWV
%    semilogx(exp(x_a(mbL:mbH)).*1000.*(mWV./mAir),Q.zRET(mbL:mbH)./1000);
%    % was ./(mAir ./ mWV)
else
%    semilogx(x_a(mbL:mbH).*1000.*(mWV./mAir),Q.zRET(mbL:mbH)./1000); % was
%    ./(mAir ./ mWV)
end
semilogx(Q.wvTradNoA,Q.zDATAnA./1000,'yo:');
semilogx(Q.wvTradNo,Q.zDATAn./1000,'cx:');
semilogx(Q.mmrSnd,Q.zsnd./1000,'g');
semilogx(Q.qvTrueRET(mbL:mbH).*1000.*(mWV./mAir),Q.zRET(mbL:mbH)./1000,'k-.');

warning off
ylabel('Altitude (km)')
xlabel('WVmmr (g/kg)')
warning on

% if reality
%     legend('RALMO','Retrieval','Ratio Analog','Ratio Digital','Radiosonde',...
%         'a priori','Location','Best')
% else
%     legend('True','Retrieval','a priori','Location','Best')
% end
ylim([Q.zRET(1) Q.zRET(end)]./1000)
%pltx = get(gca,'XLim'); plty = get(gca,'YLim'); axis([0 pltx(2) 2
%plty(2)]) pltx = get(gca,'XLim'); plot(pltx,[Q.zRET(fini)
%Q.zRET(fini)]./1000,'k--')
if logWV
    pltx(1) = min(exp(x_a(mbL:mbH)).*1000.*(mWV./mAir));%./5;
    pltx(2) = max(exp(x_a(mbL:mbH)).*1000.*(mWV./mAir));%.*5;
else
    pltx(1) = min(x_a(mbL:mbH).*1000.*(mWV./mAir));%./5;
    pltx(2) = max(x_a(mbL:mbH).*1000.*(mWV./mAir));%.*5;
end
plty = get(gca,'YLim');
pltx = [.001 10];
axis([.001 10 0 in.zOEM./1000]);
semilogx(pltx,[round(Q.zRET(fini)) floor(Q.zRET(fini))]./1000,'k--')

handfig(6) = figure;
subplot(1,2,1)
plot(xaAlpha(mbL:mbH),Q.zRET(mbL:mbH)./1000)
hold on
plot(xAlpha(mbL:mbH),Q.zRET(mbL:mbH)./1000)
%xlabel('Aerosol Extinction (10^6 m^{-1})')
xlabel('Optical Depth at 355 nm')
ylabel('Altitude (km)')
legend('a priori','Retrieval')
hold on
pltx = get(gca,'XLim');
plot(pltx,[Q.zRET(fini) Q.zRET(fini)]./1000,'k--')

subplot(1,2,2)
xaAlphaD = interp1(Q.zRET(mbL:mbH),xaAlpha(mbL:mbH),Q.zDATAn,'linear');
plot(Q.tauRno.*exp(-xaAlphaD),Q.zDATAn./1000)
hold on
xAlphaD = interp1(Q.zRET(mbL:mbH),xAlpha(mbL:mbH),Q.zDATAn,'linear');
plot(Q.tauRno.*exp(-xAlphaD),Q.zDATAn./1000)
%xlabel('Aerosol Extinction (10^6 m^{-1})')
xlabel('Transmission at 355 nm')
ylabel('Altitude (km)')
legend('a priori','Retrieval')
hold on
pltx = get(gca,'XLim');
plot(pltx,[Q.zRET(fini) Q.zRET(fini)]./1000,'k--')

handfig(7) = figure;
alp = derivative(xAlphaD(3:end-2))./derivative(Q.zDATAn(3:end-2));
alpR = derivative(-log(Q.tauRno))./derivative(Q.zDATAn);
salp = smooth(alp,10);
plot(alpR*1e6,Q.zDATAn./1000)
hold on
plot(salp*1e6,Q.zDATAn(3:end-2)./1000)
xlabel 'Extinction (10^{6} m^{-1})'
ylabel 'Altitude(km)'
title 'Extinction'
legend('Molecular Extinction','Aerosol Extinction')
plot([0 0],[0 round(Q.zDATAn(end)./1000)],'k:')

handfig(8) = figure;
hold on
if logAlpha
    plot(X.e(m+1:2*m).*100,Q.zRET./1000)
    plot(X.es(m+1:2*m).*100,Q.zRET./1000)
else
    plot(abs(X.e(m+1:2*m)./X.x(m+1:2*m).*100),Q.zRET./1000)
    plot(abs(X.es(m+1:2*m)./X.x(m+1:2*m).*100),Q.zRET./1000)
end
xlabel 'Uncertainty (%)'
ylabel 'Altitude (km)'
legend('Optical Depth Smoothing','Optical Depth Total')

SolapA = Solap;

digWVgo = 2*mchanA + 1;
digN2go = 2*mchanA + mchanD +1;
anN2go = mchanA + 1;

SxsigmaN = X.G(1:m,mchanA+1:2*mchanA)*R.KsigmaNA*SsigmaN*R.KsigmaNA'...
    *X.G(1:m,mchanA+1:2*mchanA)' + X.G(1:m,digN2go:mdata)*R.KsigmaN...
    *SsigmaN*R.KsigmaN'...
    *X.G(1:m,digN2go:mdata)';
SxsigmaH = X.G(1:m,1:mchanA)*R.KsigmaHA*SsigmaH*R.KsigmaHA'*X.G(1:m,1:mchanA)' + ...
    X.G(1:m,digWVgo:digN2go-1)*R.KsigmaH*SsigmaH*R.KsigmaH'...
    *X.G(1:m,digWVgo:digN2go-1)';
SxsigmaSHR = X.G(1:m,1:mchanA)*R.KsigmaSHRA*SsigmaR*R.KsigmaSHRA'*X.G(1:m,1:mchanA)' ...
    + X.G(1:m,digWVgo:digN2go-1)*R.KsigmaSHR*SsigmaR*R.KsigmaSHR'*...
    X.G(1:m,digWVgo:digN2go-1)';
SxsigmaSNR = X.G(1:m,anN2go:2*mchanA)*R.KsigmaSNRA*SsigmaR*R.KsigmaSNRA'*...
    X.G(1:m,anN2go:2*mchanA)' + X.G(1:m,digN2go:mdata)*R.KsigmaSNR*SsigmaR...
    *R.KsigmaSNR'*X.G(1:m,digN2go:mdata)';
SxAirH = X.G(:,1:mchanA)*R.KairA(1:mchanA,:)*Sair*R.KairA(1:mchanA,:)'...
    *X.G(:,1:mchanA)' + X.G(:,digWVgo:digN2go-1)...
    *R.Kair(digWVgo:digN2go-1,:)*Sair*R.Kair(digWVgo:digN2go-1,:)'...
    *X.G(:,digWVgo:digN2go-1)';
SxAirN = X.G(:,mchanA+1:2*mchanA)*R.KairA(mchanA+1:2*mchanA,:)...
    *Sair*R.KairA(mchanA+1:2*mchanA,:)'*X.G(:,mchanA+1:2*mchanA)'...
    + X.G(:,digN2go:mdata)*R.Kair(digN2go:mdata,:)...
    *Sair*R.Kair(digN2go:mdata,:)'*X.G(:,digN2go:mdata)';
%SxSlopeA = X.G(:,1:mchanA)*R.KslopeA*Sslope*R.KslopeA'*X.G(:,1:mchanA)'...
%    + X.G(:,anN2go:digWVgo-1)*R.KslopeA*Sslope*R.KslopeA'*X.G(:,anN2go:digWVgo-1)';
SxSlope = X.G(:,digN2go:mdata)*R.Kslope*Sslope*R.Kslope'*X.G(:,digN2go:mdata)'...
    + X.G(:,digWVgo:digN2go-1)*R.Kslope*Sslope*R.Kslope'*X.G(:,digWVgo:digN2go-1)';
%SxSlope = SxSlopeA.*(in.slopeA./in.slope).^2 + SxSlopeD;
SxOlapH = X.G(:,1:mchanA)*R.KolapA(1:mchanA,:)...
    *Solap*R.KolapA(1:mchanA,:)'*X.G(:,1:mchanA)'...
    + X.G(:,digWVgo:digN2go-1)*R.Kolap(digWVgo:digN2go-1,:)...
    *Solap*R.Kolap(digWVgo:digN2go-1,:)'*X.G(:,digWVgo:digN2go-1)';
SxOlapN = X.G(:,mchanA+1:2*mchanA)*R.KolapA(mchanA+1:2*mchanA,:)...
    *SolapA*R.KolapA(mchanA+1:2*mchanA,:)'*X.G(:,mchanA+1:2*mchanA)'...
    + X.G(:,digN2go:mdata)*R.Kolap(digN2go:mdata,:)...
    *Solap*R.Kolap(digN2go:mdata,:)'*X.G(:,digN2go:mdata)';
SxOlapN = X.G(:,digN2go:mdata)*R.Kolap(digN2go:mdata,:)...
    *Solap*R.Kolap(digN2go:mdata,:)'*X.G(:,digN2go:mdata)';

sigmaNerr = sqrt(diag(SxsigmaN(1:m,1:m)));
sigmaHerr = sqrt(diag(SxsigmaH(1:m,1:m)));
sigmaRerrH = sqrt(diag(SxsigmaSHR(1:m,1:m)));
sigmaRerrN = sqrt(diag(SxsigmaSNR(1:m,1:m)));

AirerrH = sqrt(diag(SxAirH(1:m,1:m)));
AirerrN = sqrt(diag(SxAirN(1:m,1:m)));
SlopeErr = sqrt(diag(SxSlope(1:m,1:m)));

OlaperrH = sqrt(diag(SxOlapH(1:m,1:m)));
OlaperrN = sqrt(diag(abs(SxOlapN(1:m,1:m)))); % last element negative?

handfig(9) = figure;
title('wvmr')
'errors only correct for log(vmr) retrieval'
% note since X.x is the log(vmr), sigma_X.x = sigma_vmr / vmr
semilogx(X.eo(mbL:mbH).*100,Q.zRET(mbL:mbH)./1000)
hold on
semilogx(X.es(mbL:mbH).*100,Q.zRET(mbL:mbH)./1000)
semilogx(sigmaRerrH(mbL:mbH).*100,Q.zRET(mbL:mbH)./1000,'--')
semilogx(sigmaRerrN(mbL:mbH).*100,Q.zRET(mbL:mbH)./1000,'--')
semilogx(sigmaHerr(mbL:mbH).*100,Q.zRET(mbL:mbH)./1000,'--')
semilogx(sigmaNerr(mbL:mbH).*100,Q.zRET(mbL:mbH)./1000,'--')
semilogx(AirerrH(mbL:mbH).*100,Q.zRET(mbL:mbH)./1000,'--')
semilogx(AirerrN(mbL:mbH).*100,Q.zRET(mbL:mbH)./1000,'--')
semilogx(SlopeErr(mbL:mbH).*100,Q.zRET(mbL:mbH)./1000,'--')
semilogx(OlaperrH(mbL:mbH).*100,Q.zRET(mbL:mbH)./1000,'-.')
semilogx(OlaperrN(mbL:mbH).*100,Q.zRET(mbL:mbH)./1000,'-.')
%plot(DTerr(mbL:mbH).*100,Q.zRET(mbL:mbH)./1000,'--')
totH = sqrt((X.eo(mbL:mbH)).^2+(X.es(mbL:mbH)).^2+(sigmaNerr(mbL:mbH)).^2+...
    (sigmaHerr(mbL:mbH)).^2+(sigmaRerrH(mbL:mbH)).^2+(sigmaRerrN(mbL:mbH)).^2+...
    (AirerrH(mbL:mbH)).^2+(AirerrN(mbL:mbH)).^2+(SlopeErr(mbL:mbH)).^2+...
    (OlaperrH(mbL:mbH)).^2+(OlaperrN(mbL:mbH)).^2); %+(DTerr(mbL:mbH)).^2);
semilogx(totH.*100,Q.zRET(mbL:mbH)./1000,'k')
hleg = legend('Statistical','Smoothing','\Gamma_{mol}, \sigma_{R_H}',...
    '\Gamma_{mol}, \sigma_{R_N}','\Gamma_{mol}, \sigma_{N}',...
    '\Gamma_{mol}, \sigma_{H}','\Gamma_{mol}, n_{air_H}','\Gamma_{mol}, n_{air_N}',...
    'Calibration Const','Overlap H', 'Overlap N','Total',...
    'Location','BestOutside');
set(hleg,'FontSize',12);
xlabel('Fractional Uncertainty (%)')
ylabel('Altitude (km)')
pltx = get(gca,'XLim');
semilogx(pltx,[Q.zRET(fini) Q.zRET(fini)]./1000,'k--')
%xlim([1e-4 50])

if savedat
    Qwv = Q;
    Xwv = X;
    Rwv = R;
    save([outPath 'wvOEM' int2str(date) dextout], 'Qwv', 'Xwv', 'Rwv')
end

% note this has to redone for optical density save figs
if savefigs
    savefig(handfig,[outPath 'wvOEM' int2str(date) fextout])
end
diary off
