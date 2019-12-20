% script for calibrating muogenic 10Be and 26Al production parameters based on depth profile data
% based on calbhcore.m from CRONUScalc: https://bitbucket.org/cronusearth/cronus-calc
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License, version 2, as published by the Free Software Foundation (www.fsf.org).
% Jakob Heyman - 2016-2019 (jakob.heyman@gu.se)

clear all;
close all;

% What version is this?
ver = '201912';

% Choose site. Options: Beacon, LeymonHigh, LeymonLow, LaCiotat
site = 'Beacon';

% run and load expage constants
make_consts_expage;
load consts_expage;

% display site
fprintf(1,'Site: %s\n',site);

% load data. All data is derived from Balco (2017):
% http://hess.ess.washington.edu/repository/muons2016/
mucal = mucalib_data(site);

% convert 10Be concentrations according to standards
[testi,stdi] = ismember(mucal.std10,consts.std10); % find index of standard conversion factors
mult10 = consts.std10_cf(stdi); % pick out conversion factor
mucal.N10 = mucal.N10 .* mult10;
mucal.N10unc = mucal.N10unc .* mult10;

% convert 26Al concentrations according to standards
[testi,stdi] = ismember(mucal.std26,consts.std26); % find index of standard conversion factors
mult26 = consts.std26_cf(stdi); % pick out conversion factor
mucal.N26 = mucal.N26 .* mult26;
mucal.N26unc = mucal.N26unc .* mult26;

% set uncertainty to minimum 2.9% (10Be) or 4.9% (26Al)
unchigh10 = mucal.N10unc./mucal.N10 > 0.029;
unclow10 = (unchigh10 < 1);
mucal.N10unc = mucal.N10unc.*unchigh10 + mucal.N10.*unclow10.*0.029;
unchigh26 = mucal.N26unc./mucal.N26 > 0.049;
unclow26 = (unchigh26 < 1);
mucal.N26unc = mucal.N26unc.*unchigh26 + mucal.N26.*unclow26.*0.049;

% initial guess at the parameters
pinit = mucal.pinit;
% Parameters in pinit are:
%   p(1)    erosion rate ([g/cm2]/ka)
%   p(2)    attenuation length for 10-Be (g/cm2)
%   p(3)    attenuation length for 26-Al (g/cm2)
%   p(4)    fstar10 (scaled by 1.0e-3)
%   p(5)    sigma010 (scaled by 1.0e-30)
%   p(6)    star26 (scaled by 1.0e-3)
%   p(7)    sigma026 (scaled by 1.0e-30)

% pre-calculated muon P at depth (without fstar and sigma0)
mucalpre = mucalib_depthcalc_precalc(site);
mucal.dz = mucalpre.dz;
mucal.Pfast10d = mucalpre.Pfast10d;
mucal.Pneg10d = mucalpre.Pneg10d;
mucal.Pfast26d = mucalpre.Pfast26d;
mucal.Pneg26d = mucalpre.Pneg26d;

% Get the number of samples.
nsamples10 = length(mucal.N10);
nsamples26 = length(mucal.N26);
nresids = nsamples10+nsamples26;
mucal.nresids = nresids;
% We have 7 parameters to fit.
npars=7;
mucal.npars = npars;

% Decay constants
l10 = consts.l10;
l26 = consts.l26;

% Prefs
Pref10 = consts.Pref10;
Pref26 = consts.Pref26;

% calculate atmospheric pressure
atm = ERA40atm(mucal.lat,mucal.lon,mucal.elv);

% Age Relative to t0=2010 - LSD tv from LSDfix
% tv = [0:10:50 60:100:50060 51060:1000:2000060 logspace(log10(2001060),7,200)];

% Fix w,Rc,SPhi, for sp and nu prod rate scaling 10 Ma back in time
LSDfix = LSD_fix(mucal.lat,mucal.lon,1E7,-1,mucal.samplingyr,consts);

mucal.tv = LSDfix.tv;

% decay factors
mucal.dcf10 = exp(-mucal.tv.*l10);
mucal.dcf26 = exp(-mucal.tv.*l26);

% surface spallation production
Psp = P_sp_expage(atm,LSDfix.Rc,LSDfix.SPhi,LSDfix.w,consts,1,1);
mucal.Psp10 = Psp.sp10 .* Pref10;
mucal.Psp26 = Psp.sp26 .* Pref26;

% muon production (without sigma0 and fstar!)
% This is pre-computed and the data is included in mucalib_depthcalc_precalc.m to save time
%
%Pmud = mucalib_Pmu(mucal.dz,atm,LSDfix.RcEst,consts.SPhiInf,1,1,consts);
%mucal.Pfast10d = Pmud.P_fast10';
%mucal.Pneg10d = Pmud.P_neg10';
%mucal.Pfast26d = Pmud.P_fast26';
%mucal.Pneg26d = Pmud.P_neg26';

% Use LM to find optimal values of erosion rate and attenuation length.
[pstar,iter]=mucalib_lm('mucalib_fun','mucalib_jac',pinit,1.0e-5,100,mucal);

% Compute the residual and J at the optimal parameters.
rstar=mucalib_fun(pstar,mucal);
Jstar=mucalib_jac(pstar,mucal);

% Compute Chi2 and pvalue.
chi2=norm(rstar,2)^2;
pvalue=1-chi2cdf(chi2,nresids-npars);

% Compute the covariance matrix for the fitted parameters.
covp=inv(Jstar'*Jstar);
sigmapstar=sqrt(diag(covp));

% Compute the correlation matrix for the fitted parameters.
for i=1:npars 
    for j=1:npars
        corp(i,j)=covp(i,j)/(sqrt(covp(i,i))*sqrt(covp(j,j)));
    end
end
disp('Correlations between fitted parameters');
corp

% Print out the fitted parameters.
fprintf(1,'Chi2 = %f   p-value = %f\n',[chi2; pvalue]);
fprintf(1,'Erosion Rate = %f ± %f ((g/cm2)/ka)\n',pstar(1),sigmapstar(1));
fprintf(1,'10Be Lsp = %f ± %f (g/cm2)\n',pstar(2),sigmapstar(2));
fprintf(1,'26Al Lsp = %f ± %f (g/cm2)\n',pstar(3),sigmapstar(3));
fprintf(1,'fstar10 = %f ± %f E-3\n',pstar(4),sigmapstar(4));
fprintf(1,'sigma010 = %f ± %f E-30\n',pstar(5),sigmapstar(5));
fprintf(1,'fstar26 = %f ± %f E-3\n',pstar(6),sigmapstar(6));
fprintf(1,'sigma026 = %f ± %f E-30\n',pstar(7),sigmapstar(7));


% plotting =========================================================================================
plotz = (0:10:7000);
tvzv = pstar(1).*1E-3 .* mucal.tv; % depth vector for tv (g/cm2)
tvzm = bsxfun(@plus,tvzv',plotz); % depth matrix for tv (one col per depth sample)

% sp production
Psp10rm = exp(-tvzm./pstar(2)); % spallation 10 depth prod ratio
Psp26rm = exp(-tvzm./pstar(3)); % spallation 26 depth prod ratio
Psp10m = bsxfun(@times,Psp10rm,mucal.Psp10'); % spallation 10Be production matrix
Psp26m = bsxfun(@times,Psp26rm,mucal.Psp26'); % spallation 26Al production matrix

% 10Be muon prod
Pfast10m = interp1(mucal.dz,mucal.Pfast10d,tvzm,'pchip') .* pstar(5).*1E-30; % Pfast10 surface prod
Pneg10m = interp1(mucal.dz,mucal.Pneg10d,tvzm,'pchip') .* pstar(4).*1E-3; % Pneg10 surface prod

% 26Al muon prod
Pfast26m = interp1(mucal.dz,mucal.Pfast26d,tvzm,'pchip') .* pstar(7).*1E-30; % Pfast26 surface prod
Pneg26m = interp1(mucal.dz,mucal.Pneg26d,tvzm,'pchip') .* pstar(6).*1E-3; % Pneg26 surface prod

% full production plus decay
Pfull10m = Psp10m + Pfast10m + Pneg10m; % full P10 matrix
Pfull26m = Psp26m + Pfast26m + Pneg26m; % full P26 matrix
Pfull10lm = bsxfun(@times,Pfull10m,mucal.dcf10'); % full P10 including decay matrix
Pfull26lm = bsxfun(@times,Pfull26m,mucal.dcf26'); % full P26 including decay matrix

% calculate full N by integration
N10v = trapz(mucal.tv',Pfull10lm); % calculated N10
N26v = trapz(mucal.tv',Pfull26lm); % calculated N26

% calculate sp production and N
Psp10lm = bsxfun(@times,Psp10m,mucal.dcf10'); % Psp10 including decay matrix
Psp26lm = bsxfun(@times,Psp26m,mucal.dcf26'); % Psp26 including decay matrix
Nsp10 = trapz(mucal.tv',Psp10lm); % calculated N10
Nsp26 = trapz(mucal.tv',Psp26lm); % calculated N26

% calculate mu production and N
Pfast10lm = bsxfun(@times,Pfast10m,mucal.dcf10'); % Pfast10 including decay matrix
Pfast26lm = bsxfun(@times,Pfast26m,mucal.dcf26'); % Pfast26 including decay matrix
Pneg10lm = bsxfun(@times,Pneg10m,mucal.dcf10'); % Pneg10 including decay matrix
Pneg26lm = bsxfun(@times,Pneg26m,mucal.dcf26'); % Pneg26 including decay matrix
Nfast10 = trapz(mucal.tv',Pfast10lm); % calculated N10
Nfast26 = trapz(mucal.tv',Pfast26lm); % calculated N26
Nneg10 = trapz(mucal.tv',Pneg10lm); % calculated N10
Nneg26 = trapz(mucal.tv',Pneg26lm); % calculated N26

% plot 10Be figure
figure; hold on; box on;
semilogx(Nsp10,plotz,'color','red');
semilogx(Nneg10,plotz,'color','green');
semilogx(Nfast10,plotz,'color','blue');
semilogx(N10v,plotz,'color','black');
semilogx(mucal.N10,mucal.Nz10,'.','color','black','markersize',15);
set(gca,'xaxislocation','top');
set(gca,'XScale','log'); % fix for matlab
xlabel('^{10}Be (atoms/g)');
ylabel('Shielding depth (g/cm^{2})');
legend('spallation','slow muons','fast muons','location','southeast');
axis([1E3 1E8 0 7000],'ij');
hold off;

% plot 26Al figure
figure; hold on; box on;
semilogx(Nsp26,plotz,'color','red');
semilogx(Nneg26,plotz,'color','green');
semilogx(Nfast26,plotz,'color','blue');
semilogx(N26v,plotz,'color','black');
semilogx(mucal.N26,mucal.Nz26,'.','color','black','markersize',15);
set(gca,'xaxislocation','top');
set(gca,'XScale','log'); % fix for matlab
xlabel('^{26}Al (atoms/g)');
ylabel('Shielding depth (g/cm^{2})');
legend('spallation','slow muons','fast muons','location','southeast');
axis([1E4 5E8 0 7000],'ij');
hold off;
% end plotting =====================================================================================
