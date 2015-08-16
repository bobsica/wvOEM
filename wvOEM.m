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

%v2-0-0
%090600: 17500, 12000, 50, 80, 2300, 787.5, 787.5, var:false/false, 52.5*6/52.5, 3/5, OD
%0905 12: 17500, 12000, 50, 80, 1500, 787.5, 787.5, var:false/false, 52.5*6/52.5, 3/5, OD
%200305 12: 17500, 12000, 20, 20, 1000, 787.5, 787.5, var:false/false, 52.5/52.5, 3/5, OD
%300305 12: 17500, 12000, 50, 80, 1200, 787.5, 787.5, var:false/false, 52.5/52.5, 3/5, OD
% cutoff: 20-0305: 3, 0905: 2, 0906: 2, 30-0305: 3

VERSION = '2-0-1'
date = 20090906; %20090906; %20150305
nb = '00';
oemStop = 17500; 
oemStopA = 12000;
coRET = 3;
coAddData = 5;
oemGo = 6*(3.5.*coRET.*coAddData);
oemGoAreal = 1*(3.5.*coRET.*coAddData); %20090906-6* for no spike
in.LRfree = 50; % was 20 on 0305, 0308 50, 200905-6 50
in.LRpbl = 80; % 50 on 0305; was 80 on otherwise
in.LRtranHeight = 2300;
corrLh = 787.5; % 90 (night) 360 (1/6grid, use 360); %360; %90 %100;
corrLalpha = 787.5; %90 0906-2000
varAV = false; % true - use variance of average, false - variance of measurements
varAVA = false; 
uFakvec = 2; % 0 use max, 1 is 0.8 cut, 2 is sum 0.8, 3 is Dof
% logOutput turns on diary, savefigs for figure file, savedat for .mat file
logOutput = false
savefigs = false
savedat = false

dextsp = [nb '30'];
dataPath = '/Users/BobSica/Dropbox/matlab/matlabWork/fromMCH/ralmodata/';
outPath = '/Users/BobSica/Dropbox/matlab/matlabWork/fromMCH/ralmoOEMwvOutput/';
if logOutput
    diaryFile = [outPath 'diary' int2str(date) dextsp 'LT-v' VERSION '.markdown'];
    if exist(diaryFile) ~= 0
        delete(diaryFile);
    end
    diary(diaryFile)
end
dext = [nb '30chan2.mat']; % extension for data file
dextout = [nb '30chan2-v' VERSION '.mat']; % extension for output file with version
dexts = [nb '00.mat'];
fext = [nb '30chan2.fig'];
fextout = [nb '30chan2-v' VERSION '.fig'];
gext = [nb '30combined.mat'];
zAoffset = 10; % 10 bins, i.e. 10*3.75=37.5 m, make dig and anal heights agree
oemGoA = oemGoAreal; % + zAoffset;
pieceWise = true;
deadTimeH = 4; %4.0; % ns
dfacDeadH = 0.01; % 0308 0.01, 0305 0.001
dfacDeadN = dfacDeadH; %0.01;
deadTimeN = 4;
normStop = 25000;
aposteriori = false; % if true use a posteriori analog variance
oStretch = 1.2; % 1.2; % stretch or shrink overlap
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
logAlpha = false; 
logWV = true;
cf = 3; %tent function on covariance
mAir = 28.966;
mWV = 18.02;
oemPath = './'; % for saving plots

% inputs for ralmo data
in.go = 3; % plus/minus in.go points in pieceWise
in.pieceWise = pieceWise;
in.date = date;
in.slope = 34; %30.14; %35; % 2015 37.88; 34 is adhoc, (30+38)/2
in.slopeA = in.slope ./ 3; % 3 is nominal, not accurate 2.75; 
% ad hoc factor 3.1096; %./ 2.569 (night); ./ 3.1096 (day)
% .* (2.76/2*0.9); % ad hoc correction until I get the correct number;
% units are g/kg, slope is corrected to vmr in makeRealWVlog
in.coRET = coRET; % 2, was 3, this coadds retrieval grid
in.coAddData = coAddData; % 6
in.zgo = oemGo;% 59; 0 means start where overlap > 0.1
%60; % 2500 %1000; %75; % 240; 100; %60;
in.zgoA = oemGoA;
in.zstop = normStop; % this is where you cut data for normalization %11e3;
in.zOEM = oemStop; % 8000; this is where you cut at end of makeQ for retrieval
in.zOEMA = oemStopA;
zCNnorm = 3000;
in.zCNnorm = min(zCNnorm,oemStopA);
%in.beamDiv = 0.09e-3; % mrad
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
in.oStretch = oStretch;
in.zAoffset = zAoffset;
in.corrLh = corrLh;
in.corrLalpha = corrLalpha;
in.uFakvec = uFakvec;
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

% 'BAD'
% Q.odRret = Q.odRret*2.5;
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
dfacSigmaN = 0.003; % ISSI recommend
dfacSigmaH = 0.003; % ISSI recommend
dfacSigmaR = 0.003; % ISSI recommend
dfacSigRamH = 0.1; 
% from Inaba and Kobayasi, but these would be the slope variance if known
dfacSigRamN = 0.1;
%dfacDiv = 0.1; %0.1;
%dfacCoefs = 0.1; % 05;
%dfacTheta = 0.1; %0.1;
dfacSlope = 0.05;
dfacAlpha = 0.25; %0.5;
dfacOD = 0.5; % 0.25, 0.20;
dfacAng = 0.1; % 0.01 0.05
dfacCNp = 0.1; % was 0.25;
dfacCHp = 0.5; 

dfacOlap0 = 0.1; %0.1; %0.25; %0.01;
zogo = 2500;
if Q.zDATAnA(end) > zogo && Q.zDATAn(1) < zogo
    fzogoA = find(Q.zDATAnA < zogo);
    zzogoA = Q.zDATAnA(fzogoA(end)+1);
    fzogoD = find(Q.zDATAn < zogo);
    zzogoD = Q.zDATAn(fzogoD(end)+1);
    %varOlapA = (dfacOlap0.*exp(-(Q.zDATAnA-zzogoA)./500) .* Q.olapA).^2;
    varOlapA = zeros(size(Q.zDATAnA));
    varOlapA(1:fzogoA(end)) = (dfacOlap0 .* Q.olapA(1:fzogoA(end))).^2;
    varOlapA(fzogoA(end)+1:end) = 1e-4 .* Q.olapA(fzogoA(end)+1:end);
    %varOlapD = (dfacOlap0*exp(-(Q.zDATAn-zzogoD)./500) .* Q.olap).^2;
    varOlapD = zeros(size(Q.zDATAn));
    varOlapD(1:fzogoD(end)) = (dfacOlap0 .* Q.olap(1:fzogoD(end))).^2;
    varOlapD(fzogoD(end)+1:end) = 1e-4 .* Q.olap(fzogoD(end)+1:end);
else
    varOlapA = (dfacOlap0 .* Q.olapA).^2;
    varOlapD = (dfacOlap0 .* Q.olap).^2;
end
Solap = diag([varOlapA; varOlapA; varOlapD; varOlapD]);

dfacAir = 0.01;
varAirA = (dfacAir .* Q.nNA./Q.N2rat).^2;
varAirD = (dfacAir .* Q.nN./Q.N2rat).^2;
Sair = diag([varAirA; varAirA; varAirD; varAirD]);

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
%    varAlpha2 = (Q.asrCorErrRet).^2;
    varAlpha = varAlpha1; % + varAlpha2;
    % 2nd term is uncertainty of baseline for "0" aerosols 
    varOD = (dfacOD .* Q.odRret).^2;
    varOD(1) = varOD(2);
    % can't have variance of 1st value be 0, but Q.odRret(1) = 0
end 

% varOla = (dfacOlap.*Q.olapRET).^2;
% varOlap = [varOla; varOla]; % for 2*m retrieval parameters for Solap
varCNp = dfacCNp.^2; %(dfacCNp .* Q.CNp).^2;
varCHp = dfacCHp.^2; %(dfacCNp .* Q.CNp).^2;
varAng = (dfacAng .* Q.Ang).^2;
%varDiv = (dfacDiv*Q.beamDiv).^2;
%varcoefs = (dfacCoefs.*Q.coefs).^2;
%varTheta = (dfacTheta.*Q.Theta).^2;
varDTH = (dfacDeadH.*Q.DeadTimeH).^2;
varDTN = (dfacDeadN.*Q.DeadTimeN).^2;
% varHdig = Q.backVarH./(Q.backH.^2);
% varNdig = Q.backVarN./(Q.backN.^2);

%'reducing background variance by 0.1/0.1' vars2 = [varlogq; varAlpha;
%varDT; Q.backVarH; Q.backVarN]; vars2 = [varlogq; varAlpha; Q.backVarH;
%Q.backVarN];
vars2 = [varlogq; varOD; varCHp; varCNp; varCNp; varAng; varDTH; varDTN;...
     Q.backVarHA; Q.backVarNA; Q.backVarH; Q.backVarN];

dzRET = Q.zRET(2) - Q.zRET(1);
lc = (round(corrLh ./ dzRET) .* dzRET) .* ones(size(Q.zRET));
lcalpha = (round(corrLalpha ./ dzRET) .* dzRET) .* ones(size(Q.zRET));

SsigmaN = (dfacSigmaN.*Q.sigmaN).^2;
SsigmaH = (dfacSigmaH.*Q.sigmaH).^2;
SsigmaR = (dfacSigmaR.*Q.sigmaR).^2;
Sslope = (dfacSlope.*Q.slope).^2;

S_a = zeros(n,n);
for i=1:m
    for j=1:m
        sigprod = sqrt(vars2(i) .* vars2(j));
        diffz = Q.zRET(i) - Q.zRET(j);
        sumlc = lc(i) + lc(j);
        shape(1) = exp(-(2.*diffz./sumlc).^2); %Gaussian correlation function
        shape(2) = exp(-2.*abs(diffz)./sumlc); %Exponential correlation function
        shape(3) = (1-(1-exp(-1)).*2.*abs(diffz)./sumlc); % Tent correlation func
        if shape(3) < 0
            shape(3) = 0;
        end
        S_a(i,j) = sigprod .* shape(cf); %Gaussian correlation function
    end
end
% option for using different correlation length for alpha
for i=m+1:2*m
    for j=m+1:2*m
        sigprod = sqrt(vars2(i) .* vars2(j));
        diffz = Q.zRET(i-m) - Q.zRET(j-m);
        sumlc = lcalpha(i-m) + lcalpha(j-m);
        shape(1) = exp(-(2.*diffz./sumlc).^2); %Gaussian correlation function
        shape(2) = exp(-2.*abs(diffz)./sumlc); %Exponential correlation function
        shape(3) = (1-(1-exp(-1)).*2.*abs(diffz)./sumlc); % Tent correlation func
        if shape(3) < 0
            shape(3) = 0;
        end
        S_a(i,j) = sigprod .* shape(cf); %.* shape(cf); %Gaussian correlation function
    end
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

[X,R] = oem(O,Q,R,@makeJ,S_a,Se,S_ainv,Seinv,x_a,y);
R  = makeParameterJacobians(Q,X.x);
if logWV
    X.vmr = exp(X.x(1:m));
else
    X.vmr = X.x(1:m);
end
if logAlpha
    xaAlpha = exp(x_a(m+1:2*m)); % + Q.odnormR);
    xAlpha = exp(X.x(m+1:2*m)); % + Q.odnormR);
    X.xAlpha = xAlpha; % only used for save file
else
    xaAlpha = x_a(m+1:2*m); % done in makeQ + Q.odnormR (from ground to start);
    xAlpha = X.x(m+1:2*m); % + Q.odnormR;
    X.xAlpha = xAlpha; % only used for save file
end

'X.cost'
X.cost

degF1 = trace(X.A(1:m,1:m));
degF2 = trace(X.A(m+1:2*m,m+1:2*m));

finiDoF = round(min(degF1,degF2));

% remove retrieval points outside of data grid
mbLo = find(Q.zRET < min(Q.zDATAnA(1),Q.zDATAn(1)));
mbHi = find(Q.zRET > max(Q.zDATAnA(end),Q.zDATAn(end)));
mbH = mbHi(1) - 1;
mbL = mbLo(end) + 1;
mbH2 = m - mbH;

% calculate error matrices
Sxsigma = X.G*R.KsigmaRay*SsigmaN*R.KsigmaRay'*X.G';
SxAir = X.G*R.Kair*Sair*R.Kair'*X.G';
SxSlope = X.G*R.Kslope*Sslope*R.Kslope'*X.G';
SxOlap = X.G*R.Kolap*Solap*R.Kolap'*X.G';

% for q
if logWV
    sigmaRayErrq = sqrt(diag(Sxsigma(1:m,1:m)));
    AirErrq = sqrt(diag(SxAir(1:m,1:m)));
    SlopeErrq = sqrt(diag(SxSlope(1:m,1:m)));
    OlapErrq = sqrt(diag(SxOlap(1:m,1:m)));
    obsErrq = X.eo(1:m);
else
    sigmaRayErrq = sqrt(diag(Sxsigma(1:m,1:m))) ./ X.vmr;
    AirErrq = sqrt(diag(SxAir(1:m,1:m))) ./ X.vmr;
    SlopeErrq = sqrt(diag(SxSlope(1:m,1:m))) ./ X.vmr;
    OlapErrq = sqrt(diag(SxOlap(1:m,1:m))) ./ X.vmr;
    obsErrq = X.eo(1:m) ./ X.vmr;
end
totErrq = sqrt(obsErrq(mbL:mbH).^2 + sigmaRayErrq(mbL:mbH).^2....
+ AirErrq(mbL:mbH).^2 + SlopeErrq(mbL:mbH).^2 + OlapErrq(mbL:mbH).^2);

% for OD
if logAlpha
    sigmaRayErro = sqrt(diag(Sxsigma(m+1:2*m,m+1:2*m)));
    AirErro = sqrt(diag(SxAir(m+1:2*m,m+1:2*m)));
    SlopeErro = sqrt(diag(SxSlope(m+1:2*m,m+1:2*m)));
    OlapErro = sqrt(diag(SxOlap(m+1:2*m,m+1:2*m)));
    obsErro = X.eo(m+1:2*m);
else
    axAlpha = abs(xAlpha); % alpha can be negative
    sigmaRayErro = sqrt(diag(Sxsigma(m+1:2*m,m+1:2*m))) ./ axAlpha;
    AirErro = sqrt(diag(SxAir(m+1:2*m,m+1:2*m))) ./ axAlpha;
    SlopeErro = sqrt(diag(SxSlope(m+1:2*m,m+1:2*m))) ./ axAlpha;
    OlapErro = sqrt(diag(SxOlap(m+1:2*m,m+1:2*m))) ./ axAlpha;
    obsErro = X.eo(m+1:2*m) ./ axAlpha;
end
totErro = sqrt(obsErro(mbL:mbH).^2 + sigmaRayErro(mbL:mbH).^2....
+ AirErro(mbL:mbH).^2 + SlopeErro(mbL:mbH).^2 + OlapErro(mbL:mbH).^2);

% output to screen
disp(' ')
outB =[round(degF1) round(degF2) length(Q.zRET)];
str = ['Degrees of Freedom (wv/OD): ', num2str(outB(1),'%0.5g'),...
    ' , ',num2str(outB(2),'%0.5g'), ' of ' num2str(outB(3)), ' degrees'];
disp(str)

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

% plot raw data
handfig(1) = figure;
subplot(1,2,1)
plot(Q.y2Hz./1e6*y(1:mchanA),Q.zDATAnA./1000,'b')
hold on
plot(Q.y2Hz./1e6*y(1+mchanA:2*mchanA),Q.zDATAnA./1000,'r')
xlabel('ADC Count Rate (MHz)')
ylabel('Altitude (km)')
legend('Analog WV','Analog N2')
subplot(1,2,2)
semilogx(Q.y2Hz./1e6*y(1+2*mchanA:2*mchanA+mchanD),Q.zDATAn./1000,'b')
hold on
semilogx(Q.y2Hz./1e6*y(1+2*mchanA+mchanD:mdata),Q.zDATAn./1000,'red')
xlabel('Photocount Rate (MHz)')
ylabel('Altitude (km)')
legend('Digital Water Vapour','Digital N2')

% plot Jacobians and averaging kernels
handfig(2) = figure;
subplot(2,2,1)
plot(X.J(1:mchanA,mbL:mbH),Q.zDATAnA./1000)
xlabel('Jacobian')
ylabel('Altitude (km)')

subplot(2,2,2)
%plot(X.J(mchan+1:mdata,(m+mbL:2*m-mbH2)),Q.zDATAn./1000)
plot(X.J(mchanA+1:2*mchanA,(m+3:2*m-mbH2)),Q.zDATAnA./1000)
xlabel('Jacobian')
ylabel('Altitude (km)')

subplot(2,2,3)
plot(X.J(2*mchanA+1:2*mchanA+mchanD,mbL:mbH)./1e4,Q.zDATAn./1000)
xlabel('Jacobian')
ylabel('Altitude (km)')

subplot(2,2,4)
%plot(X.J(mchan+1:mdata,(m+mbL:2*m-mbH2)),Q.zDATAn./1000)
plot(X.J(2*mchanA+mchanD+1:mdata,(m+3:2*m-mbH2))./1e5,Q.zDATAn./1000)
xlabel('Jacobian')
ylabel('Altitude (km)')

handfig(3) = figure;
subplot(1,2,1)
hk = plot(X.A(mbL:mbH,mbL:mbH),Q.zRET(mbL:mbH)./1000);
hold on
unit = ones(size(Q.zRET(mbL:mbH)));
response = X.A(mbL:mbH,mbL:mbH)*unit;
fak = find(diag(X.A(mbL:mbH,mbL:mbH)) >= 0.8);
fak2 = find(response >= 0.8);
if isempty(fak2)
    fak2 = 1;
end
fakvec = [fak(end); fak2(end); finiDoF]'
if uFakvec == 0
    fini = finiDoF;
else
    fini = fakvec(uFakvec);
end
%plot(response,Q.zRET(mbL:mbH)./1000,'r:')
set(hk,'LineWidth',0.5)
xlabel('WV Averaging Kernel')
ylabel('Altitude (km)')
%axis([-0.1 1.1 Q.zRET(1)./1000 Q.zRET(end)./1000])
xlim([-0.1 1.1])
pltx = get(gca,'XLim');
plot(pltx,[Q.zRET(fini) Q.zRET(fini)]./1000,'k--')

subplot(1,2,2)
hk = plot(X.A((m+mbL+2:2*m-mbH2),(m+mbL+2:2*m-mbH2)),Q.zRET(mbL+2:mbH)./1000);
hold on
unit = ones(size(Q.zRET(mbL+2:mbH)));
response = X.A((m+mbL+2:2*m-mbH2),(m+mbL+2:2*m-mbH2))*unit;
%plot(response,Q.zRET(mbL+2:mbH)./1000,'r:')
set(hk,'LineWidth',0.5)
xlabel('OD Averaging Kernel')
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
zw = Q.zRET(mbL:mbH);
fpltw = find(width ~= 0);
handfig(4) = figure;
plot(width(fpltw),zw(fpltw)./1000)
hold on
xlabel('Vertical Resolution (m)')
ylabel('Alitude (km)')
%axis([0 300 Q.zRET(1)./1000 Q.zRET(end)./1000])
pltx = get(gca,'XLim');
plot([0 pltx(2)],[Q.zRET(fini) Q.zRET(fini)]./1000,'k--')

% plot count residuals
handfig(5) = figure;
subplot(2,2,1)
plot((X.yf(1:mchanA)-y(1:mchanA))./y(1:mchanA).*100,Q.zDATAnA./1000,'b')
hold on
plot(sqrt(yvar(1:mchanA))./y(1:mchanA)*100,Q.zDATAnA./1000,'r')
plot(-sqrt(yvar(1:mchanA))./y(1:mchanA)*100,Q.zDATAnA./1000,'r')
% pltx = get(gca,'XLim');
% plot(pltx,[maskLow maskLow]./1000,'k--')
% plot(pltx,[maskHigh maskHigh]./1000,'k--')
xlabel('WV Residuals - Analog (%)')
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
xlabel('N2 Residuals - Analog (%)')
ylabel('Altitude (km)')
%ylim([0 1.05.*oemStop./1000])

subplot(2,2,3)
plot((X.yf(2*mchanA+1:2*mchanA+mchanD)-y(2*mchanA+1:2*mchanA+mchanD))...
    ./y(2*mchanA+1:2*mchanA+mchanD).*100,Q.zDATAn./1000,'b')
hold on
plot(sqrt(Q.yTrueH)./Q.yTrueH*100,Q.zDATAn./1000,'r')
%legend('Residuals','Poisson Noise','Location','Best')
plot(-sqrt(Q.yTrueH)./Q.yTrueH*100,Q.zDATAn./1000,'r')
% plot(sqrt(yvar(2*mchanA+1:2*mchanA+mchanD))./Q.yTrueH*100,Q.zDATAn./1000,'gx')
% plot(-sqrt(yvar(2*mchanA+1:2*mchanA+mchanD))./Q.yTrueH*100,Q.zDATAn./1000,'gx')
xlabel('WV Residuals - Digital (%)')
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
% plot(sqrt(yvar(2*mchanA+mchanD+1:mdata))./y(2*mchanA+mchanD+1:mdata)...
%     *100,Q.zDATAn./1000,'gx')
% plot(-sqrt(yvar(2*mchanA+mchanD+1:mdata))./y(2*mchanA+mchanD+1:mdata)...
%     *100,Q.zDATAn./1000,'gx')
xlabel('N2 Residuals - Digital (%)')
ylabel('Altitude (km)')
%xlim([-20 20])
ylim([0 1.05.*oemStop./1000])

% plot retrievals of n and q
handfig(6) = figure;
if reality
    lralmot = length(ralmo.t);
    if exist('ralmo.records','var') == 0
        ralmo.records = 1;
        lralmot = length(ralmo.z);
    end
    nRalmo = lralmot ./ ralmo.records;
    itst = 1:1:lralmot;
    mR = mod(itst,nRalmo);
    fmR = find(mR ==0);
    if fmR > length(ralmo.t)
        rTime = ralmo.t(1);
    else
        rTime = ralmo.t(fmR);
    end
    fRalmo = find(rTime > Q.ralmoTimeEnd);
    iRalmo = fRalmo(1);
    'RALMO time ='
    datevec(rTime(iRalmo))
    'Analysis time ='
    datevec(Q.ralmoTimeEnd)
    gg = 1 + (iRalmo-1)*nRalmo;
    hwv1=semilogx(ralmo.q(gg:nRalmo+gg-1),(ralmo.z(gg:nRalmo+gg-1)-490)./1000,'b');
    hold on
else
%    semilogx(Q.qvTrueRET.*1000.*(mWV./mAir),Q.zRET./1000); % was ./(mAir
%    ./ mWV)
    hold on
end
fd1 = find(Q.zRET > oemGoA);
fdh = find(Q.zRET > oemStopA);
hwv2 = semilogx(X.vmr(mbL:mbH).*1000.*(mWV./mAir).*(outSlopeA./in.slope),...
    Q.zRET(mbL:mbH)./1000,'r','LineWidth',2);
hwv3 = semilogx(Q.wvTradNoA,Q.zDATAnA./1000,'yo:');
hwv4 = semilogx(Q.wvTradNo,Q.zDATAn./1000,'cx:');
hwv5 = semilogx(Q.mmrSnd,Q.zsnd./1000,'g');
hwv6=semilogx(Q.qvTrueRET(mbL:mbH).*1000.*(mWV./mAir),Q.zRET(mbL:mbH)./1000,'k-.');

warning off
ylabel('Altitude (km)')
xlabel('WV Mixing Ratio (g/kg)')
warning on

hleg = legend([hwv1,hwv2,hwv5,hwv6],'Traditional','OEM','Sonde','a priori',...
    'Location','SouthWest');
set(hleg,'FontSize',8,'Box','off');

ylim([Q.zRET(1) Q.zRET(end)]./1000)
% if logWV
%     pltx(1) = min(exp(x_a(mbL:mbH)).*1000.*(mWV./mAir));%./5;
%     pltx(2) = max(exp(x_a(mbL:mbH)).*1000.*(mWV./mAir));%.*5;
% else
%     pltx(1) = min(x_a(mbL:mbH).*1000.*(mWV./mAir));%./5;
%     pltx(2) = max(x_a(mbL:mbH).*1000.*(mWV./mAir));%.*5;
% end
plty = get(gca,'YLim');
pltx = get(gca,'XLim'); %[.0001 10];
% axis([.0001 10 0 in.zOEM./1000]);
semilogx(pltx,[round(Q.zRET(fini)) floor(Q.zRET(fini))]./1000,'k--')

handfig(7) = figure;
fbad = find(isnan(Q.mmrSnd) == 1);
fgd = fbad(1)-1;
fhg = find(Q.zRET < Q.zsnd(fgd));
ralRET = interp1((ralmo.z(gg:nRalmo+gg-1)-490),ralmo.q(gg:nRalmo+gg-1),...
    Q.zRET(mbL:fhg(end)),'linear');
sndRET = interp1(Q.zsnd(1:fgd),Q.mmrSnd(1:fgd),Q.zRET(mbL:fhg(end)),'linear');
xvmr = X.vmr(mbL:fhg(end)).*1000.*(mWV./mAir).*(outSlopeA./in.slope);
plot((xvmr-sndRET)./sndRET*100,Q.zRET(mbL:fhg(end))./1000,'ro:')
hold on
plot((ralRET-sndRET)./sndRET*100,Q.zRET(mbL:fhg(end))./1000,'bx:')
ylim([0 round(Q.zRET(fini)./1000)+1])
plot([0 0],[0 round(Q.zRET(fini)./1000)+1],'k:')
%pltx = get(gca,'XLim');
%plot(pltx,[Q.zRET(fini) Q.zRET(fini)]./1000,'k--')
xlabel('Percent Difference Versus Sonde')
legend('OEM','RALMO')
ylabel('Altitude (km)')

handfig(8) = figure;
plot(exp(-xaAlpha(mbL:mbH)),Q.zRET(mbL:mbH)./1000)
hold on
plot(exp(-xAlpha(mbL:mbH)),Q.zRET(mbL:mbH)./1000,'r')
plot(Q.tauRnoA,Q.zDATAnA./1000,'g')
xlabel('Transmission (355 nm)')
ylabel('Altitude (km)')
hleg = legend('a prioiri aerosol','retrieved aerosol','molecular',...
    'Location','East');
set(hleg,'FontSize',8,'Box','off');
plot(Q.tauRno,Q.zDATAn./1000,'g')
if logAlpha
    Tup = exp(-xAlpha(mbL:mbH)).*sqrt((X.eo(mbL+m:mbH+m))); %exp(-xAlpha(mbL:mbH))+
%    Tdn = exp(-xAlpha(mbL:mbH))-exp(-xAlpha(mbL:mbH)).*sqrt((X.eo(mbL+m:mbH+m)));
else
    Tup = sqrt((X.eo(mbL+m:mbH+m))); % exp(-xAlpha(mbL:mbH)) +
%    Tdn = exp(-xAlpha(mbL:mbH)) - sqrt((X.eo(mbL+m:mbH+m)));
end

horiz_errorbar(Q.zRET(mbL:5:mbH)'./1000,exp(-xAlpha(mbL:5:mbH)),Tup(1:5:end),'r*')
%jbfilly(Q.zRET(mbL:mbH)'./1000,Tup',Tdn','r','r',0,.25);
xlim([0 2.25])
ylim([0 Q.zRET(end)./1000])
pltx = get(gca,'XLim');
hold on
plot(pltx,[Q.zRET(fini) Q.zRET(fini)]./1000,'k--')

handfig(9) = figure;
%subplot(1,2,1)
plot(xaAlpha(mbL:mbH),Q.zRET(mbL:mbH)./1000)
hold on
plot(xAlpha(mbL:mbH),Q.zRET(mbL:mbH)./1000)
xlabel('Optical Depth at 355 nm')
ylabel('Altitude (km)')
legend('a priori','Retrieval','Location','NorthWest')
if logAlpha
    Tup = xAlpha(mbL:mbH)+exp(-xAlpha(mbL:mbH)).*sqrt((X.eo(mbL+m:mbH+m)));
    Tdn = xAlpha(mbL:mbH)-exp(-xAlpha(mbL:mbH)).*sqrt((X.eo(mbL+m:mbH+m)));
else
    Tup = xAlpha(mbL:mbH) + sqrt((X.eo(mbL+m:mbH+m)));
    Tdn = xAlpha(mbL:mbH) - sqrt((X.eo(mbL+m:mbH+m)));
end
jbfilly(Q.zRET(mbL:mbH)'./1000,Tup',Tdn','r','r',0,.25);
pltx = get(gca,'XLim');
plot(pltx,[Q.zRET(fini) Q.zRET(fini)]./1000,'k--')

handfig(10) = figure;
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

% plot errors
handfig(11) = figure;
% note since X.x is the log(vmr), sigma_X.x = sigma_vmr / vmr
plot(obsErrq(mbL:mbH).*100,Q.zRET(mbL:mbH)./1000)
hold on
plot(sigmaRayErrq(mbL:mbH).*100,Q.zRET(mbL:mbH)./1000,'--')
plot(AirErrq(mbL:mbH).*100,Q.zRET(mbL:mbH)./1000,'--')
plot(SlopeErrq(mbL:mbH).*100,Q.zRET(mbL:mbH)./1000,'--')
plot(OlapErrq(mbL:mbH).*100,Q.zRET(mbL:mbH)./1000,'-.')
plot(totErrq.*100,Q.zRET(mbL:mbH)./1000,'k')
hleg = legend('Statistical','\sigma_{Rayleigh}', 'Air Density',...
    'Calibration','Overlap','Total','Location','Best');
set(hleg,'FontSize',8,'Box','off');
xlabel('Mixing Ratio Uncertainty (%)')
ylabel('Altitude (km)')
pltx = get(gca,'XLim');
plot(pltx,[Q.zRET(fini) Q.zRET(fini)]./1000,'k--')
%xlim([0 50])

handfig(12) = figure;
%subplot(1,2,2)
% optical depth and transmission have the same percentage error
% below is only correct for alpha retrieval, not log alpha
plot(obsErro(mbL:mbH).*100,Q.zRET(mbL:mbH)./1000)
hold on
plot(sigmaRayErro(mbL:mbH).*100,Q.zRET(mbL:mbH)./1000,'--')
plot(AirErro(mbL:mbH).*100,Q.zRET(mbL:mbH)./1000,'--')
plot(SlopeErro(mbL:mbH).*100,Q.zRET(mbL:mbH)./1000,'--')
plot(OlapErro(mbL:mbH).*100,Q.zRET(mbL:mbH)./1000,'-.')
plot(totErro.*100,Q.zRET(mbL:mbH)./1000,'k')
hleg = legend('Statistical','\sigma_{Rayleigh}', 'Air Density',...
    'Calibration','Overlap','Total','Location','Best');
set(hleg,'FontSize',8,'Box','off');
%xlabel('Optical Depth Uncertainty (%)')
xlabel('Transmission Uncertainty (%)')
ylabel('Altitude (km)')
pltx = get(gca,'XLim');
plot(pltx,[Q.zRET(fini) Q.zRET(fini)]./1000,'k--')
%xlim([0 2])

handfig(13) = figure;
plot(Q.asrDATA(20:end),Q.zDATAn(20:end)./1000,'b')
hold on
plot(Q.asrDATAA,Q.zDATAnA./1000,'b')
xlabel('Backscatter Ratio')
ylabel('Altitude (km)')

handfig(14) = figure;   % extinction coefficient plots
subplot(1,2,1)
Npoly = 4;                 % Order of polynomial fit
Fsg = 35;                 % Window length
[bsg,gsg] = sgolay(Npoly,Fsg);   % Calculate S-G coefficients
xflip = flip(xAlpha);
sdflip = flip(sqrt(X.eo(1+m:2*m))); %variance of xAlpha
HalfWin  = ((Fsg+1)/2) -1;
for n = (Fsg+1)/2:length(xAlpha)-(Fsg+1)/2,
  SG1(n) = dot(gsg(:,2),xflip(n - HalfWin:n + HalfWin));
  if ~logAlpha
    SG1N(n) = dot(gsg(:,2),sdflip(n - HalfWin:n + HalfWin)); 
  else
      'Errors for extinction wrong for logAlpha true'
      stooopppp
  end
end
SG1 = SG1 ./ (Q.zRET(1)-Q.zRET(2)); % Turn differential into derivative
%varReduce = 1./sum(gsg(:,2).^2'); % (1) - (2) because array is flipped!
SG1N = SG1N ./ (Q.zRET(1)-Q.zRET(2));
%sqrt(SG1N./varReduce)./(abs(Q.zRET(1)-Q.zRET(2)));
xaAlphaD = interp1(Q.zRET(mbL:mbH),xaAlpha(mbL:mbH),Q.zDATAn,'linear');
xAlphaD = interp1(Q.zRET(mbL:mbH),xAlpha(mbL:mbH),Q.zDATAn,'linear');
alp = derivative(xAlphaD(3:end-2))./derivative(Q.zDATAn(3:end-2));
alpRET = derivative(xAlpha(mbL:mbH))./derivative(Q.zRET(mbL:mbH));
alpR = derivative(-log(Q.tauRno))./derivative(Q.zDATAn);
salp = smooth(alpRET,19);
plot(alpR*1e6,Q.zDATAn./1000)
hold on
% plot(salp*1e6,Q.zRET(mbL:mbH)./1000)
plot(SG1*1e6,flip(Q.zRET(1:end-(Fsg+1)/2)./1000))
%plot(SG1N*1e6,flip(Q.zRET(1:end-(Fsg+1)/2)./1000))
xlabel 'Extinction (10^{6} m^{-1})'
ylabel 'Altitude(km)'
%title 'Extinction'
hleg = legend('Molecular','Aerosol');
set(hleg,'FontSize',8,'Box','off');
%plot([0 0],[0 round(Q.zDATAn(end)./1000)],'k:')
pltx = get(gca,'XLim');
plot(pltx,[Q.zRET(fini) Q.zRET(fini)]./1000,'k--')
subplot(1,2,2)
plot(alpR*1e6,Q.zDATAn./1000)
hold on
%plot(salp*1e6,Q.zRET(mbL:mbH)./1000)
xlabel 'Extinction (10^{6} m^{-1})'
ylabel 'Altitude(km)'
plty = get(gca,'YLim');
ylim([2.5 plty(2)])
horiz_errorbar(flip(Q.zRET(1:10:end-(Fsg+1)/2)./1000),SG1(1:10:end)*1e6...
    ,SG1N(1:10:end)*1e6)
hold on
pltx = get(gca,'XLim');
plot(pltx,[Q.zRET(fini) Q.zRET(fini)]./1000,'k--')

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
if logOutput
    diary off
end

