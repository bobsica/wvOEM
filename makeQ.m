function [Q,y,yvar] = makeQ(in)
% makeQ(in)
% -Usage-
%	[Q,y,yvar] = makeQ(in)
%
% -Inputs-
%	in - data structure of parameters passed from oemWV
%
% -Outputs-
%	Q - retrieval input structure
%   y - n channels of data
%   yvar - the data's variance
%
% currently 4 channels, retrieving 2*m + 10 parameters 
%
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

mAir = 28.966;
mWV = 18.02;
N2rat = 0.781; % (ISSI, amount of N2 in air)
slope = (in.slope.*mAir) ./ (N2rat.*mWV);
slopeA = (in.slopeA.*mAir) ./ (N2rat.*mWV);
% whoops above were inverted until 23Jun15
Rate = 30;
%believeOlap = 5e-3;

% Raman cross sections
nu0 = 1./(3371e-8); % known at 337.1 nm
sig0N = 3.5e-34; % m2/sr
sig0H = 7.8e-34;
nuN = 1./(3866.9e-8);
nuH = 1./(4074.9e-8);
Hshift = 3657; % cm-1
Nshift = 2330.7;
N20 = (nu0 - Nshift).^4;
H2O0 = (nu0 - Hshift).^4;
sigRamN = sig0N .* (nuN - Nshift).^4 ./ N20; % m2/sr
sigRamH = sig0H .* (nuH - Hshift).^4 ./ H2O0;

% Rayleigh cross sections
boltz = 1.3806488e-23;
clight = 299792458; %ISSI value
lambda = 0.3547;
lambdaH = 0.40749;
lambdaN = 0.38669;
A=4.02*10^(-28);B=-0.3228;C=0.389;D=0.09426;
expR = 4+B+C*lambda+D/lambda;
sigmaR = A / lambda^(expR).*1e-4; % m^2
expN = 4+B+C*lambdaN+D/lambdaN;
sigmaN = A / lambdaN^(expN).*1e-4;
expH = 4+B+C*lambdaH+D/lambdaH;
sigmaH = A / lambdaH^(expH).*1e-4;

% for get_RALMO files (in MHz)
dataPath = '/Users/BobSica/Dropbox/matlab/matlabWork/fromMCH/ralmodata/';
ralmoFile = [dataPath 'S0S3' int2str(in.date) in.dext];
load(ralmoFile)
load './ralmoFixedOverlap.mat'

dzRaw =  zCounts(2) - zCounts(1);
y2HzRaw = clight ./ (2.*(deltaTime.*Rate).*dzRaw);

S.Sn2 = N2counts; % already in counts ./ (y2HzRaw./1e6);
S.Swv = WVcounts; % ./ (y2HzRaw./1e6);
S.Sn2A = N2countsA;
S.SwvA = WVcountsA;
S.z = zCounts;
tmpwv = y2HzRaw .* WVcounts;
S.cSwv = (tmpwv./(1-tmpwv.*4e-9)) ./ y2HzRaw;
tmpn2 = y2HzRaw .* N2counts;
S.cSn2 = (tmpn2./(1-tmpn2.*4e-9)) ./ y2HzRaw;

% coadd
if in.coAddData == 1
    'this needs to be fixed/check'
    ssttoopp2
    if in.zgo == 0
        oolap = ones(size(S.z));
        fltt = find(S.z <= ralmoO.zoverlap(end));
        oolap(fltt) = interp1(ralmoO.zoverlap,ralmoO.overlap,S.z(fltt),'linear');
        op1 = find(oolap > believeOlap);
        dgo = find(S.z >= S.z(oolap(1)));
    else
        dgo = find(S.z >= in.zgo);
    end
    dgo = find(S.z >= in.zgo);
    dstp = find(S.z <= in.zstop); %1.1.*hiNorm);
    zN = S.z(dgo(1):dstp(end));
    SHcoadd = S.Swv(dgo(1):dstp(end));
    SNcoadd = S.Sn2(dgo(1):dstp(end));
    SHcoaddA = S.SwvA(dgo(1):dstp(end));
    SNcoaddA = S.Sn2A(dgo(1):dstp(end));
%     SHerrA = S.SwvvarA(dgo(1):dstp(end));
%     SNerrA = S.Sn2varA(dgo(1):dstp(end));
else
    [tmpaddN,zzN] = coadd(S.Sn2,S.z,in.coAddData);
    [tmpaddH,zzH] = coadd(S.Swv,S.z,in.coAddData);
    [tmpaddNA,zzN] = coadd(S.Sn2A,S.z,in.coAddData);
    [tmpaddHA,zzH] = coadd(S.SwvA,S.z,in.coAddData);
    if in.zgo == 0
        'check this option'
        ssssttttoooopppp
        oolap = ones(size(zzN));
        fltt = find(zzN <= ralmoO.zoverlap(end));
        oolap(fltt) = interp1(ralmoO.zoverlap,ralmoO.overlap,zzN(fltt),'linear');
        op1 = find(oolap > believeOlap);
        dgo = find(zzN >= zzN(op1(1)));
    else
        dgo = find(zzN >= in.zgo);
        dgoA = find(zzN >= in.zgoA);
    end
    dstp = find(zzN <= in.zstop);
    SHcoadd = in.coAddData .* tmpaddH;
    SNcoadd = in.coAddData .* tmpaddN;
    SHcoaddA = in.coAddData .* tmpaddHA;
    SNcoaddA = in.coAddData .* tmpaddNA;
%     SHerrA = tmpaddHAvar(dgo(1):dstp(end)) .* SHcoaddA(dgo(1):dstp(end));
%     SNerrA = tmpaddNAvar(dgo(1):dstp(end)) .* SNcoaddA(dgo(1):dstp(end));
    
%     backH = backH .* in.coAddData; backVarH = in.coAddData.^2 .* backVarH;
%     backN = backN .* in.coAddData; backVarN = in.coAddData.^2 .* backVarN;

% background and background variance
    zHback = 50e3; zNback = 50e3;
    findBH = find(zzN > zHback);
    backH = mean(SHcoadd(findBH(1):end));
    findBN = find(zzN > zNback);
    backN = mean(SNcoadd(findBN(1):end));
    % variance of average or measurement
    if in.varAV
        backVarH = (std(SHcoadd(findBH(1):end)) ...
            ./ sqrt(length(SHcoadd(findBH(1):end)))).^2;
        backVarN = (std(SNcoadd(findBH(1):end)) ...
            ./ sqrt(length(SNcoadd(findBN(1):end)))).^2;
    else
        backVarH = (std(SHcoadd(findBH(1):end))).^2; % variance of measurement
        backVarN = (std(SNcoadd(findBN(1):end))).^2;
    end
    findBH = find(zzN > zHback);
    backHA = mean(SHcoaddA(findBH(1):end)-in.Aoffset);
    findBN = find(zzN > zNback);
    backNA = mean(SNcoaddA(findBN(1):end)-in.Aoffset);
    if in.varAVA
        backVarHA = (std(SHcoaddA(findBH(1):end)-in.Aoffset)...
         ./ sqrt(length(SHcoaddA(findBH(1):end)))).^2;
        backVarNA = (std(SNcoaddA(findBN(1):end)-in.Aoffset)...
         ./ sqrt(length(SNcoaddA(findBN(1):end)))).^2;
    else
        backVarHA = (std(SHcoaddA(findBH(1):end)-in.Aoffset)).^2;
        backVarNA = (std(SNcoaddA(findBN(1):end)-in.Aoffset)).^2;
    end

    SHcoadd = SHcoadd(dgo(1):dstp(end));
    SNcoadd = SNcoadd(dgo(1):dstp(end));
    SHcoaddA = SHcoaddA(dgoA(1):dstp(end));
    SNcoaddA = SNcoaddA(dgoA(1):dstp(end));
    zN = zzN(dgo(1):dstp(end));
    zNA = zzN(dgoA(1):dstp(end)) - in.zAoffset;
end    
% logic above is coadding analog decrease by 1/sqrt(bins), and the
% uncertainty is the standard deviation relative to the analog
% backscattered photons, so you take the background out for a percentage to
% multiple the counts by.

dzDATA = zN(2) - zN(1);
y2Hz = clight ./ (2.*(deltaTime.*Rate).*dzDATA);

sndChk = datevec(timeEnd);
%if sndChk(4) < 10
    sndDate = in.date;
%else
%    sndDate = in.date + 1;
%end
sndFile = [dataPath 'snd' int2str(sndDate) in.dexts];
%sndFile = [dataPath 'RS92v' int2str(in.date+1) '000000.mat'];
load(sndFile);
%load('./ralmodata/fromAH20120717/snd.mat')
% if RS92v uncomment next 2 lines
%snd.ptu = RS92; snd.ptu.gph = snd.ptu.z;

%snd.ptu = out.snd;
itop = length(snd.ptu.gph);
fzb = find(snd.ptu.gph(1:itop) < in.asl);
if isempty(fzb)
    fzb = 0;
end
sndgph = snd.ptu.gph(fzb(end)+1:itop);
sndT = snd.ptu.T(fzb(end)+1:itop);
sndp = snd.ptu.p(fzb(end)+1:itop);
sndrh = snd.ptu.rh(fzb(end)+1:itop);

[zSndA,isnd] = sort(sndgph);
[zSnd,ia,ic] = unique(zSndA);
zsndASL = zSnd - in.asl;
Tsndi = sndT(isnd(ia));
psndi = sndp(isnd(ia));
nsndi = psndi.*100 ./ (boltz .* Tsndi); % psnd is in mb
RHsndi = sndrh(isnd(ia));
mmri = rh2mr(Tsndi,psndi,RHsndi);
mmriVol = mmri .* (mAir./mWV);
nden = exp(interp1(zsndASL,log(nsndi),zN,'linear'));
ndenA = exp(interp1(zsndASL,log(nsndi),zNA,'linear'));

% fill in density data below first sample and above last if required for
% last
ScaleHeight = 8.771e+3;
N0 = 2.504e25; % m^-3
i = 1;
while zN(i) < zsndASL(1)
    nden(i) = nsndi(1) .* exp(-(zN(i)-zsndASL(1))./ScaleHeight);
    i = i + 1;
end 
tden = find(isnan(nden) == 1);
if isempty(tden)
    nNz = N2rat .* nden;
else % if you need air density above sonde burst use barometric law
    normDen = nden(tden(1)-1);
    nBaro = normDen .* exp(-(zN-zN(tden(1)-1))./ScaleHeight); % Molecular profile
    nNz = zeros(size(zN));
    nNz(1:tden(1)-1) = N2rat .* nden(1:tden(1)-1); % Leblanc et al
    nNz(tden(1):end) = N2rat .* nBaro(tden(1):end);
end

i = 1;
while zNA(i) < zsndASL(1)
    ndenA(i) = nsndi(1) .* exp(-(zNA(i)-zsndASL(1))./ScaleHeight);
    i = i + 1;
end 
tden = find(isnan(nden) == 1);
if isempty(tden)
    nNzA = N2rat .* ndenA;
else % if you need air density above sonde burst use barometric law
    normDen = ndenA(tden(1)-1);
    nBaroA = normDen .* exp(-(zNA-zNA(tden(1)-1))./ScaleHeight); % Molecular profile
    nNzA = zeros(size(zN));
    nNzA(1:tden(1)-1) = N2rat .* ndenA(1:tden(1)-1); % Leblanc et al
    nNzA(tden(1):end) = N2rat .* nBaroA(tden(1):end);
end

% wv mixing ratio; use US standard for shape, sonde at zson m for
% normalization
zsong = 500; %zN(10); % picked zN(10) to avoid zN(1) less than first sonde point
zsons = 1000;
load('./USstandardWV')
usmmr = interp1(wvUSstandard.z,wvUSstandard.mmr,zsndASL,'linear');
fmmrSndg = find(zsndASL > zsong);
fmmrSnds = find(zsndASL > zsons);
mmrMod = trapz(zsndASL(fmmrSndg(1):fmmrSnds(1)-1),usmmr(fmmrSndg(1):fmmrSnds(1)-1));
mmrSnd = trapz(zsndASL(fmmrSndg(1):fmmrSnds(1)-1),mmri(fmmrSndg(1):fmmrSnds(1)-1));
usmmrNorm = (mmrSnd./mmrMod) .* usmmr;


%usmmrNorm = (mmrSnd./mmrMod) .* usmmr;
qtrue = interp1(zsndASL,usmmrNorm,zN,'linear'); % mass mixing ratio
qvtrue = (qtrue./1000) .* (mAir./mWV);
nH = qvtrue .* nden; 

molePath = ['/Users/BobSica/Dropbox/matlab/matlabWork/fromMCH/ralmodata/alphaBeta'...
    int2str(in.date) '.mat'];
load(molePath);
beta_mol_DATA = interp1(alphaBeta.z-in.asl,alphaBeta.beta_mol,zN,'linear');
alpha_mol_DATA = interp1(alphaBeta.z-in.asl,alphaBeta.alpha_mol,zN,'linear');

'below asr max use asr max value (3000 m max)'
fcut = find(asr.z < 3000);
[amx,imx] = max(asr.asr(fcut));
fasrMasr = find(asr.z < asr.z(imx));
asr.asr(fasrMasr) = amx; 

asrDATA = interp1(asr.z,asr.asr,zN,'linear');
beta_mol_DATAA = interp1(alphaBeta.z-in.asl,alphaBeta.beta_mol,zNA,'linear');
alpha_mol_DATAA = interp1(alphaBeta.z-in.asl,alphaBeta.alpha_mol,zNA,'linear');
asrDATAA = interp1(asr.z,asr.asr,zNA,'linear');
%zzASR = asr.z; zASR = asr.z;
flow = find(zN < asr.z(1));
if ~isempty(flow)
 asrDATA(flow) = asrDATA(flow(end)+1);
end
fhi = find(zN > asr.z(end));
if ~isempty(fhi)
  asrDATA(fhi) = 1;
end
% fchk = find(isnan(asrDATA));
% fchk2 = find(isnan(beta_mol_DATA));
% fchk3 = find(isnan(alpha_mol_DATA));
% if ~isempty(fchk | fchk2 | fchk3)
%    'asrDATA or alpha/beta from model has NaNs in makeQ9, stop'
%    stop
% end
flowA = find(zNA < asr.z(1));
if ~isempty(flowA)
 asrDATAA(flowA) = asrDATAA(flowA(end)+1);
end
fhiA = find(zNA > asr.z(end));
if ~isempty(fhiA)
  asrDATAA(fhiA) = 1;
end

LR = in.LRfree * ones(size(asrDATA));
fff = find(zN < in.LRtranHeight);
LR(fff) = in.LRpbl;

LRA = in.LRfree * ones(size(asrDATAA));
fffA = find(zNA < in.LRtranHeight);
LRA(fffA) = in.LRpbl;

'unsmoothed a prior ASR'
%asrDATAs = asrDATA;
asrDATAs = asrDATA; %smooth(asrDATA,90); %was 45
asrDATAsA = asrDATAA; %smooth(asrDATAA,90);
%'smooth a priori for asr' fend = find(asrDATA == 1); plin =
%polyfit(zN(1:fend(1)),asrDATA(1:fend(1)),1); asrDATAs = polyval(plin,zN);
fneg = find(asrDATAs < 1);
asrDATAs(fneg) = 1;
fnegA = find(asrDATAsA < 1);
asrDATAsA(fnegA) = 1;
%basrZ = find(zN > 8000); asrDATAs(basrZ) = 0; basr = find(asrDATAs < 0);
%asrDATAs(basr) = 0;
alphaAer = LR .* (beta_mol_DATA .* (asrDATAs-1));
znoAer = find(zN > 12000); % was 3000 for 20130122
alphaAer(znoAer) = 0;
'asr set to 0 > 12000'
%if in.logAlpha
    fl0 = find(alphaAer <= 0);
    alphaAer(fl0) = 1e-12;
%end
alphaAerA = LRA .* (beta_mol_DATAA .* (asrDATAsA-1));
fl0A = find(alphaAerA <= 0);
alphaAerA(fl0A) = 1e-12;

%'setting extinction = 0' alphaAer = zeros(size(zN)); 'ad hoc extinction'
%alphaAer = 1e-5.*ones(size(zN)); alphaAer = alphaAer./10; alphaCorErr =
%((sdasrCor./asrCor).*alpha_mol_DATA); % estimate of background SD
alphaCorErr = 0;
z0 = 0:.1:zN(1);
alpha0 = alphaAer(1) .* ones(size(z0));
odnorm = trapz(z0,alpha0);
odAer = cumtrapz(zN,alphaAer) + odnorm; 
odAerH = cumtrapz(zN,alphaAer.*lambda./lambdaH) + odnorm;
odAerN = cumtrapz(zN,alphaAer.*lambda./lambdaN) + odnorm;
% note in v11 we don't include the bit from the ground up in odAer,
% it is in the retrieved "C*tau"
z0A = 0:.1:zNA(1);
alpha0A = alphaAerA(1) .* ones(size(z0A));
odnormA = trapz(z0A,alpha0A);
odAerA = cumtrapz(zNA,alphaAerA) + odnormA; 
odAerHA = cumtrapz(zNA,alphaAerA.*lambda./lambdaH) + odnormA;
odAerNA = cumtrapz(zNA,alphaAerA.*lambda./lambdaN) + odnormA;

figure
subplot(2,1,1)
plot(asrDATAs,zN./1000)
hold on
plot(asrDATAsA,zNA./1000)
xlabel 'ASR (\beta_{tot}/\beta_{mol})'
ylabel 'Altitude (km)'
%ylim([0 12.5])
subplot(2,1,2)
semilogx(alphaAer*1e6,zN./1000);
hold on
semilogx(alphaAerA*1e6,zNA./1000);
%hold on plot(asrCorErr*1e6,zN./1000,'--')
xlabel 'Aerosol Extinction (10^6 m^{-1})'
ylabel 'Altitude (km)'
%ylim([0 12.5])

% find transmission molecular
intnBaro = ScaleHeight .* N0 .* (1 - exp(-(zN(1)./ScaleHeight))); 
tauR0no = exp(-sigmaR.*intnBaro); 
tauH0no = exp(-sigmaH.*intnBaro); 
tauN0no = exp(-sigmaN.*intnBaro); 
tauRno = tauR0no .* exp(-cumtrapz(sigmaR.*nNz./N2rat).*dzDATA);
tauHno = tauH0no .* exp(-cumtrapz(sigmaH.*nNz./N2rat).*dzDATA); 
tauNno = tauN0no .* exp(-cumtrapz(sigmaN.*nNz./N2rat).*dzDATA); 
tauR = tauRno .* exp(-odAer);
tauH = tauHno .* exp(-odAerH);
tauN = tauNno .* exp(-odAerN);

intnBaroA = ScaleHeight .* N0 .* (1 - exp(-(zNA(1)./ScaleHeight))); 
tauR0noA = exp(-sigmaR.*intnBaroA); 
tauH0noA = exp(-sigmaH.*intnBaroA); 
tauN0noA = exp(-sigmaN.*intnBaroA); 
tauRnoA = tauR0noA .* exp(-cumtrapz(sigmaR.*nNzA./N2rat).*dzDATA);
tauHnoA = tauH0noA .* exp(-cumtrapz(sigmaH.*nNzA./N2rat).*dzDATA); 
tauNnoA = tauN0noA .* exp(-cumtrapz(sigmaN.*nNzA./N2rat).*dzDATA); 
tauRA = tauRnoA .* exp(-odAerA);
tauHA = tauHnoA .* exp(-odAerHA);
tauNA = tauNnoA .* exp(-odAerNA);

% overlap, olapD is height derivative of the overlap for makeR.m
olap = ones(size(zN));
olapD = zeros(size(zN));
olapA = ones(size(zNA));
olapDA = zeros(size(zNA));
% overlap stretch factor
ralmoO.zoverlap = in.oStretch .* ralmoO.zoverlap;
flt = find(zN <= ralmoO.zoverlap(end));
olap(flt) = interp1(ralmoO.zoverlap,ralmoO.overlap,zN(flt),'linear');
olapD(flt) = interp1(ralmoO.zoverlap,ralmoO.overlapD,zN(flt),'linear');
fltA = find(zNA <= ralmoO.zoverlap(end));
olapA(fltA) = interp1(ralmoO.zoverlap,ralmoO.overlap,zNA(fltA),'linear');
olapDA(fltA) = interp1(ralmoO.zoverlap,ralmoO.overlapD,zNA(fltA),'linear');

%'background out of variance digital calculation'
yObsN = y2Hz .* (SNcoadd); %-backN);
yObsH = y2Hz .* (SHcoadd); %-backH);
backNHz = y2Hz .* backN;
backHHz = y2Hz .* backH;
yTrueN = (yObsN./(1-yObsN*4e-9)) ./ y2Hz;
backTN = (backNHz./(1-backNHz*4e-9)) ./ y2Hz;
yTrueH = (yObsH./(1-yObsH*4e-9)) ./ y2Hz;
backTH = (backHHz./(1-backHHz*4e-9)) ./ y2Hz;

% cut export to what you specified in the beginning
dendN = find(zN <= in.zOEM);
Q.zDATAn = zN(1:dendN(end));
Q.SNcoadd = SNcoadd(1:dendN(end));
Q.SHcoadd = SHcoadd(1:dendN(end));
Q.yTrueN = yTrueN(1:dendN(end));
Q.yTrueH = yTrueH(1:dendN(end));
Q.backTH = backTH;
Q.backTN = backTN;

%dgoNA = find(zNA >= in.zgoA);
dendNA = find(zNA <= in.zOEMA);
Q.zDATAnA = zNA(1:dendNA(end));
Q.SNcoaddA = SNcoaddA(1:dendNA(end));
Q.SHcoaddA = SHcoaddA(1:dendNA(end));
y = [Q.SHcoaddA; Q.SNcoaddA; Q.SHcoadd; Q.SNcoadd];

% fix digital error at high rates
if in.pieceWise
    lzD = length(zN);
    lzA = length(zNA);
    go = 6; %12; %24
    stop = go-1;
    j = 0;
    for i = go:lzA-stop
        j = j + 1;
        [pp,spp,ppregress] = fitlinenp(zNA(i-stop:i+stop),SHcoaddA(i-stop:i+stop));
        tmp = pp(1).*zNA(i-stop:i+stop) + pp(2);
        varWVA(i) = (std(SHcoaddA(i-stop:i+stop) - tmp)).^2;
        [pp,spp,ppregress] = fitlinenp(zNA(i-stop:i+stop),SNcoaddA(i-stop:i+stop));
        tmp = pp(1).*zNA(i-stop:i+stop) + pp(2);
        varN2A(i) = (std(SNcoaddA(i-stop:i+stop) - tmp)).^2;
    end
    j = 0;
    for i = go:lzD-stop
        j = j + 1;
        [pp,spp,ppregress] = fitlinenp(zN(i-stop:i+stop),yTrueH(i-stop:i+stop));
        tmp = pp(1).*zN(i-stop:i+stop) + pp(2);
        varWV(i) = (std(yTrueH(i-stop:i+stop) - tmp)).^2;
        [pp,spp,ppregress] = fitlinenp(zN(i-stop:i+stop),yTrueN(i-stop:i+stop));
        tmp = pp(1).*zN(i-stop:i+stop) + pp(2);
        varN2(i) = (std(yTrueN(i-stop:i+stop) - tmp)).^2;
    end
    
    WVvarA = zeros(size(SHcoaddA));
    N2varA = zeros(size(SNcoaddA));
    WVvarA(go:lzA) = varWVA;
    WVvarA(1:go-1) = WVvarA(go);
    WVvarA(lzA+1:end) = WVvarA(lzA);
    N2varA(go:lzA) = varN2A;
    N2varA(1:go-1) = N2varA(go);
    N2varA(lzA+1:end) = N2varA(lzA);
    WVvar = zeros(size(SHcoadd));
    N2var = zeros(size(SNcoadd));
    WVvar(go:lzD) = varWV;
    WVvar(1:go-1) = WVvar(go);
    WVvar(lzD+1:end) = WVvar(lzD);
    N2var(go:lzD) = varN2;
    N2var(1:go-1) = N2var(go);
    N2var(lzD+1:end) = N2var(lzD);

%     fSH = find(WVvar < backVarH);
%     fSN = find(N2var < backVarN);
%     fSHA = find(WVvarA < backVarHA);
%     fSNA = find(N2varA < backVarNA);
%' let variance be < background'
   
else
% ad hoc DT error estimate

'dead code!'
stoooopp

%     yTrueN1 = (yObsN./(1-yObsN*4.04e-9)) ./ y2Hz;
%     yTrueH1 = (yObsH./(1-yObsH*4.04e-9)) ./ y2Hz;
% 
%     trueErrN = (abs(yTrueN1-yTrueN)).^2;
%     trueErrH = (abs(yTrueH1-yTrueH)).^2;
% 
%     warning off
%     figure
%     semilogx(Q.yTrueN-backN,Q.zDATAn)
%     hold on
%     semilogx(trueErrN,Q.zDATAn)
%     legend('Signal Variance','DT Variance'); title 'N2 Digital'
%     figure
%     semilogx(Q.yTrueH-backH,Q.zDATAn)
%     hold on
%     semilogx(trueErrH,Q.zDATAn)
%     legend('Signal Variance','DT Variance'); title 'WV Digital'
%     warning on
% 
%     fdtN = find(trueErrN > (Q.SNcoadd + backVarN));
%     SNerr(fdtN) = trueErrN(fdtN);
%     fdtH = find(trueErrH > (Q.SHcoadd + backVarH));
%     SHerr(fdtH) = trueErrH(fdtH);
% 
%     SHerr = Q.SHcoadd;
%     SNerr = Q.SNcoadd;
%     fSH = find(SHerr < sqrt(backVarH));
%     fSN = find(SNerr < sqrt(backVarN));
%     fSHA = find(SHerrA < sqrt(backVarHA));
%     fSNA = find(SNerrA < sqrt(backVarNA));
%     SHerr(fSH) = backVarH;
%     SNerr(fSN) = backVarN;
%     SHerrA(fSHA) = backVarHA;
%     SNerrA(fSNA) = backVarNA;
end

% use a posteriori analog variance, but either way apply covariance mask
if in.aposteriori
    dataPath = '/Users/BobSica/Dropbox/matlab/matlabWork/fromMCH/ralmodata/';
    ralmoFile = [dataPath 'S0S3' int2str(in.date) in.dextsp 'aPostDA.mat'];
    load(ralmoFile)
    yvar = [WVvarA(1:dendNA(end)); N2varA(1:dendNA(end)); WVvar(1:dendN(end));...
       N2var(1:dendN(end))];   
else
% can't have a 0 variance or variance < background
    f00 = find(WVvarA < backVarHA);
    WVvarA(f00) = backVarHA;
    f00 = find(N2varA < backVarNA);
    N2varA(f00) = backVarNA;
    f00 = find(WVvar < backVarH);
    WVvar(f00) = backVarH;
    f00 = find(N2var < backVarN);
    N2var(f00) = backVarN;
    WVvarT = yTrueH;
    f00 = find(WVvarT <= backVarH);
    WVvarT(f00) = backVarH;
    N2varT = yTrueN;
    f00 = find(N2varT <= backVarH);
    N2varT(f00) = backVarN;
    
    if in.pieceWise
%        yvar = [WVvarA(1:dendNA(end)); N2varA(1:dendNA(end));...
%               WVvar(1:dendN(end)); N2var(1:dendN(end))];
       yvar = [WVvarA(1:dendNA(end)); N2varA(1:dendNA(end));...
           WVvarT(1:dendN(end)); N2varT(1:dendN(end))];
%        yvar = [smooth(WVvarA(1:dendNA(end)),24); 
%           smooth(N2varA(1:dendNA(end)),36);...
%            WVvarT; N2varT];     
    else
        'dead code II'
        stooop
       yvar = [SHerrA(1:dendNA(end)); SNerrA(1:dendNA(end)); SHerr; SNerr];
    end
end

figure
semilogx(WVvarA,zNA)
hold on
semilogx(WVvarT,zN)
semilogx(yTrueH,zN)
semilogx(N2varA,zNA)
semilogx(N2var,zN)
semilogx(yTrueN,zN)
legend('pWVA','WVvar','SHtrue','pN2A','N2var','SNtrue')
title('Piecewise Variance (A/D) + Poisson')

% Find a priori CN' at height zCNnorm
fmmrz = find(zN > in.zCNnorm);
fmmrzA = find(zNA > in.zCNnorm);
CNp = ((SNcoadd(fmmrz(1)) - backN) .* zN(fmmrz(1)).^2) ./ (nNz(fmmrz(1))...
    .* tauR(fmmrz(1)) .* tauN(fmmrz(1)) .* olap(fmmrz(1))); 
CHp = slope .* CNp; % note slope is in vmr units
CNpA = ((SNcoaddA(fmmrzA(1)) - backNA) .* zNA(fmmrzA(1)).^2)...
    ./ (nNzA(fmmrzA(1)) .* tauRA(fmmrzA(1)) .* tauNA(fmmrzA(1))...
    .* olapA(fmmrzA(1)));
%CHpA = ((SHcoaddA(fmmrz(1)) - backHA) .* zN(fmmrz(1)).^2) ./ (nNz(fmmrz(1))...
%    .* tauR(fmmrz(1)) .* tauH(fmmrz(1)) .* olap(fmmrz(1)));
CHpA = slopeA .* CNpA; % note slope is in vmr units

% 'wrong A variance'
% yvar = [Q.SHcoaddA; Q.SNcoaddA; SHerr; SNerr];
Q.CNpA = CNpA;
Q.CHpA = CHpA; %slope .* CNpA;
Q.CNp = CNp;
Q.slope = slope;
Q.slopeA = slopeA;
Q.CHp = slope .* CNp;
% note this slope is fixed in makeRealWVlog like: 
%CHp = (0.781.*x(end-2).*1000.*Q.mWV) ./ (Q.slope(mmr).*Q.mAir);

Q.olap = olap(1:dendN(end));
Q.olapD = olapD(1:dendN(end));
Q.olapA = olapA(1:dendNA(end));
Q.olapDA = olapDA(1:dendNA(end));
Q.tauR = tauR(1:dendN(end));
Q.tauH = tauH(1:dendN(end));
Q.tauN = tauN(1:dendN(end));
Q.tauRA = tauRA(1:dendNA(end));
Q.tauHA = tauHA(1:dendNA(end));
Q.tauNA = tauNA(1:dendNA(end));
Q.tauRno = tauRno(1:dendN(end));
Q.tauHno = tauHno(1:dendN(end));
Q.tauNno = tauNno(1:dendN(end));
Q.tauRnoA = tauRnoA(1:dendNA(end));
Q.tauHnoA = tauHnoA(1:dendNA(end));
Q.tauNnoA = tauNnoA(1:dendNA(end));

Q.Ang = 1; % a priori Angstrom exponent
%Q.tauR0 = tauR0; Q.tauH0 = tauH0; Q.tauN0 = tauN0;
Q.odnormR = odnorm;
Q.odnormRA = odnormA;
%Q.odnormH = odnormH; Q.odnormN = odnormN;
Q.backH = backH;
Q.backN = backN;
Q.backHA = backHA;
Q.backNA = backNA;
Q.nN = nNz(1:dendN(end));
Q.nNA = nNzA(1:dendNA(end));
Q.sigmaH = sigmaH;
Q.sigmaN = sigmaN;
Q.sigmaR = sigmaR;
Q.sigRamN = sigRamN;
Q.sigRamH = sigRamH;
Q.backVarN = backVarN;
Q.backVarH = backVarH;
Q.backVarNA = backVarNA;
Q.backVarHA = backVarHA;
Q.y2Hz = y2Hz;
Q.N2rat = N2rat;
Q.mAir = mAir;
Q.mWV = mWV;
%'fix Q.alphaRdata'
Q.alphaRdata = alphaAer;
%Q.tauRdata = tauAer; Q.alphaRdata = alpha_mol_DATA(1:dendN(end))./10;
%Q.asrCorErr = asrCorErr;
Q.alphaCorErr = alphaCorErr;
Q.mzsnd = zsndASL;
%Q.nsndi = nsndi; Q.LR = LR(1:dendN(end)); Q.beta_mol_DATA =
%beta_mol_DATA(1:dendN(end)); Q.alpha_mol_DATA =
%alpha_mol_DATA(1:dendN(end));
Q.nsnd = nsndi;
Q.zsnd = zsndASL;
Q.Tsnd = Tsndi;
Q.RHsnd = RHsndi;
Q.vmrSnd = mmriVol;
Q.mmrSnd = mmri;
Q.ralmoTimeEnd = timeEnd;
Q.qvTrue = qvtrue;

% make data from forward model; x(1:m) is log(q) on retrieval grid;
% x(m+1:2*m) is extinction Retrieval grid
% digital
% dzRET = in.coRET .* dzDATA;
% altRET = (Q.zDATAn(1)-dzDATA:dzRET:in.zOEM)';
% if altRET(end) <= Q.zDATAn(end)
%     zRET = [altRET; altRET(end)+dzRET];
% else
%     zRET = altRET;
% end
% % analog
% altRETA = (Q.zDATAnA(1)-dzDATA:dzRET:in.zOEMA)';
% if altRETA(end) <= Q.zDATAnA(end)
%     zRETA = [altRETA; altRETA(end)+dzRET];
% else
%     zRETA = altRETA;
% end

dzRET = in.coRET .* dzDATA;
altRET = (min(Q.zDATAn(1),Q.zDATAnA(1))-dzDATA:dzRET:max(in.zOEM,in.zOEMA))';
if altRET(end) <= max(Q.zDATAn(end),Q.zDATAnA(end))
    zRET = [altRET; altRET(end)+dzRET];
else
    zRET = altRET;
end

Q.zRET = zRET;
nzRET = N2rat .* exp(interp1(zsndASL,log(nsndi),zRET,'linear'));
% nden = exp(interp1(zsndASL,log(nsndi),zN,'linear'));
% nNz = N2rat .* nden; % Leblanc et al
Q.nNret = nzRET; % Leblanc et al
Q.alphaRret = interp1(zN,alphaAer,Q.zRET,'linear','extrap');
Q.alphaRdata = alphaAer(1:dendN(end));
Q.asrCorErrRet = alphaCorErr; 
Q.odRret = interp1(zN,odAer,Q.zRET,'linear', .1.*odAer(1)); %'extrap');
Q.odRdata = odAer(1:dendN(end));

qtrueRET = interp1(zsndASL,usmmrNorm,zRET,'linear'); % mass mixing ratio
qvtrueRET = (qtrueRET./1000) .* (mAir./mWV);
Q.qvTrueRET = qvtrueRET;

qtrue = interp1(zsndASL,usmmrNorm,zN,'linear'); % mass mixing ratio
qvtrue = (qtrue./1000) .* (mAir./mWV);

Q.tauRret = interp1(zN,tauR,Q.zRET,'linear','extrap');
Q.tauHret = interp1(zN,tauH,Q.zRET,'linear','extrap');
Q.tauNret = interp1(zN,tauN,Q.zRET,'linear','extrap');

olapDret = zeros(size(Q.zRET));
fltr = find(Q.zRET <= ralmoO.zoverlap(end));

figure
subplot(1,2,1)
plot(Q.tauR,Q.zDATAn./1000,'r')
hold on
%plot(Q.tauRret,Q.zRET./1000,'rx')
plot(Q.tauH,Q.zDATAn./1000,'b')
%plot(Q.tauHret,Q.zRET./1000,'bx')
plot(Q.tauN,Q.zDATAn./1000,'g')
%plot(Q.tauNret,Q.zRET./1000,'gx')
plot(tauRno(1:dendN(end)),Q.zDATAn./1000,'r--')
plot(tauHno(1:dendN(end)),Q.zDATAn./1000,'b--')
plot(tauNno(1:dendN(end)),Q.zDATAn./1000,'g--')
plot(exp(-odAer(1:dendN(end))),Q.zDATAn./1000,'r:')
plot(exp(-odAerH(1:dendN(end))),Q.zDATAn./1000,'b:')
plot(exp(-odAerN(1:dendN(end))),Q.zDATAn./1000,'g:')
xlabel 'Transmission'
ylabel 'Altitude(km)'
legend('\tau_R','\tau_H','\tau_N','\tau_{Rmol}','\tau_{Hmol}',...
    '\tau_{Nmol}','\tau_{Ra}','\tau_{Ha}','\tau_{Na}',...
    'Location','BestOutside')
title 'Transmission (Molecular/Aerosol/Total; Digital)'
subplot(1,2,2)
plot(Q.tauRA,Q.zDATAnA./1000,'r')
hold on
%plot(Q.tauRret,Q.zRET./1000,'rx')
plot(Q.tauHA,Q.zDATAnA./1000,'b')
%plot(Q.tauHret,Q.zRET./1000,'bx')
plot(Q.tauNA,Q.zDATAnA./1000,'g')
%plot(Q.tauNret,Q.zRET./1000,'gx')
plot(tauRnoA(1:dendNA(end)),Q.zDATAnA./1000,'r--')
plot(tauHnoA(1:dendNA(end)),Q.zDATAnA./1000,'b--')
plot(tauNnoA(1:dendNA(end)),Q.zDATAnA./1000,'g--')
plot(exp(-odAerA(1:dendNA(end))),Q.zDATAnA./1000,'r:')
plot(exp(-odAerHA(1:dendNA(end))),Q.zDATAnA./1000,'b:')
plot(exp(-odAerNA(1:dendNA(end))),Q.zDATAnA./1000,'g:')
xlabel 'Transmission'
ylabel 'Altitude(km)'
%legend('\tau_R','\tau_H','\tau_N','\tau_{Rmol}','\tau_{Hmol}',...
%     '\tau_{Nmol}','\tau_{Ra}','\tau_{Ha}','\tau_{Na}',...
%     'Location','SouthWest')
title 'Transmission (Molecular/Aerosol/Total; Analog)'

%figure
tauHR = tauH./tauR;
tauNR = tauN./tauR;
tauHRno = tauHno./tauRno;
tauNRno = tauNno./tauRno;

% lbackNA = polyval(pN,Q.zDATAnA); %Q.backNA; %
% lbackHA = polyval(pH,Q.zDATAnA); %Q.backHA; %
Q.wvTrad = in.slope .* (Q.tauH./Q.tauN) .* ((Q.yTrueH - Q.backTH)...
    ./(Q.yTrueN - Q.backTN));
Q.wvTradNo = in.slope .* (Q.tauHno(1:dendN(end))...
    ./ Q.tauNno(1:dendN(end)))...
    .* ((Q.yTrueH - Q.backTH)./(Q.yTrueN - Q.backTN));
Q.wvTradNoA = in.slopeA .* (Q.tauHno(1:dendNA(end))...
    ./ Q.tauNno(1:dendNA(end)))...
    .* ((Q.SHcoaddA - Q.backHA)./(Q.SNcoaddA - Q.backNA));
% lwvTradNoA = in.slopeA .* (Q.tauHno(1:dendNA(end))...
%     ./ Q.tauNno(1:dendNA(end)))...
%     .* ((Q.SHcoaddA - lbackHA)./(Q.SNcoaddA - lbackNA));

% figure
% plot(Q.wvTradNo,Q.zDATAn)
% hold on
% plot(Q.wvTradNo,Q.zDATAn,'o')
% plot(Q.wvTradNoA,Q.zDATAnA)
% plot(lwvTradNoA,Q.zDATAnA,'x')
% xlabel('WV (g/kg)')
% legend('Digital WV','Dig WV-molecular only','Analog Constant','Analog linear')
% ylim([0 6000])
% 
% figure
% wvint = interp1(Q.zDATAn,Q.wvTradNo,Q.zDATAnA,'linear');
% ratwv = Q.wvTradNoA ./ wvint;
% lr = find(Q.zDATAnA > 800); % 0
% hr = find(Q.zDATAnA > 1600); % 800
% [pRat,sigmapRat,regressRat] = fitlinenp(Q.zDATAnA(1:hr(1)),ratwv(1:hr(1)));
% rlin = polyval(pRat,Q.zDATAnA);
% plot(ratwv,Q.zDATAnA)
% hold on
% plot(ratwv-rlin+1,Q.zDATAnA)

return

