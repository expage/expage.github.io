function depthcalcET()

% Depth profile 10Be/26Al exposure age calculator calculating the best fit erosion rate E and
% exposure duration T within the range Emin-Emax/Tmin-Tmax (one period of exposure with constant
% erosion rate). The calculation is done using the expage calculator production rates (nuclide-
% specific LSD scaling) and assuming one period of exposure. At least two concentrations must be
% given for the calculator to work.
% This is free software: you can use/copy/modify/distribute as long as you keep it free/open.
% Jakob Heyman - 2015-2018 (jakob.heyman@gu.se)

clear;
close all;
tic();

% age and erosion ranges to test ===================================================================
Tmin = 1; % yr
Tmax = 500000; % yr
Emin = 0; % mm/ka
Emax = 50; % mm/ka
% ==================================================================================================

% number of E and T points for the test grid
ETgrid = 100;

% What version is this?
ver = '201902';

% read input file
% NOTE! sample thickness in standard input is here changed to sample mid-point depth (cm)
% erosion is included to have the same input as for depthcalc.m but it is not used here
[sample_name,lat,long,elv,aa,depth,rho,othercorr,erosion,N10,delN10,be_stds,N26,delN26,al_stds,...
    samplingyr] = textread('input.txt','%s %n %n %n %s %n %n %n %n %n %n %s %n %n %s %n',...
    'commentstyle','matlab');

% run and load expage constants
make_consts_expage;
load consts_expage;

% constants
Pref10 = consts.P10_ref_nu; delPref10 = consts.delP10_ref_nu;
Pref26 = consts.P26_ref_nu; delPref26 = consts.delP26_ref_nu;
% Decay constant
l10 = consts.l10; dell10 = consts.dell10;
l26 = consts.l26; dell26 = consts.dell26;

% convert 10Be concentrations according to standards
for i = 1:numel(N10);
    mult10(i,1) = consts.be_stds_cfs(strcmp(be_stds(i),consts.be_stds_names));
end;
N10 = N10 .* mult10;
delN10 = delN10 .* mult10;

% convert 26Al concentrations according to standards
for i = 1:numel(N26);
    mult26(i,1) = consts.al_stds_cfs(strcmp(al_stds(i),consts.al_stds_names));
end;
N26 = N26 .* mult26;
delN26 = delN26 .* mult26;

% fix longitude values
long(find(long < 0)) = long(find(long < 0)) + 360;

% fix sample pressure
std_v = strcmp(aa,'std');
ant_v = strcmp(aa,'ant');
pre_v = strcmp(aa,'pre');
pressure(std_v) = ERA40atm(lat(std_v),long(std_v),elv(std_v));
pressure(ant_v) = antatm(elv(ant_v));
pressure(pre_v) = elv(pre_v);

% check and fix inputs =============================================================================
diffstr = '';

if sum(lat ~= lat(1)) > 0;
    diffstr = [diffstr,' latitude'];
end;
if sum(long ~= long(1)) > 0;
    diffstr = [diffstr,' longitude'];
end;
if sum(elv ~= elv(1)) > 0;
    if (strcmp(aa(1),'pre'));
        diffstr = [diffstr,' atm-pressure'];
    else;
        diffstr = [diffstr,' elevation'];
    end;
end;
if sum(rho ~= rho(1)) > 0;
    diffstr = [diffstr,' density'];
end;
if sum(othercorr ~= othercorr(1)) > 0;
    diffstr = [diffstr,' shielding'];
end;

% if differences...
if numel(diffstr) > 0;
    fprintf(1,'The input samples have differences in the following parameter');
    if numel(diffstr) > 14; fprintf(1,'s'); end;
    fprintf(1,':%s\n',diffstr);
end;

lat = mean(lat);
long = mean(long);
elv = mean(elv);
shield = mean(othercorr);
rho = mean(rho);
samplingyr = mean(samplingyr);
% ==================================================================================================

% use mean atmospheic pressure
atm = mean(pressure);

% Set nucl to 0 for both 10/26
nucl10 = 0; nucl26 = 0;

% Check and count if 10Be and 26Al is measured
n10 = sum(N10+delN10 > 0);
n26 = sum(N26+delN26 > 0);
if n10>1; nucl10 = 1; end;
if n26>1; nucl26 = 1; end;

% mesh for plotting
xx = linspace(Emin,Emax,5e3);
yy = linspace(Tmin,Tmax,5e3)';

% Vectors for T and E
Tvect = linspace(Tmin,Tmax,ETgrid);
Evect = linspace(Emin,Emax,ETgrid) .* 1e-4; % fix cm/yr

% E and T limits for plotting
ETlims = [Emin Emax Tmin Tmax];

% for loop counters
j10 = 0;
j26 = 0;

% Age Relative to t0=2010 - LSD tv from LSDfix
% tv = [0:10:50 60:100:50060 51060:1000:2000060 logspace(log10(2001060),7,200)];

% Fix w,Rc,SPhi, for sp and mu prod rate scaling
LSDfix = LSD_fix(lat,long,Tmax+2010-samplingyr,-1,consts);

% time vector tv1
tv1 = LSDfix.tv;

% adjust tv, Rc, and SPhi to sampling year
if samplingyr <= 2010;
    clipidx = min(find(tv1 > 2010-samplingyr));
    tv = [2010-samplingyr tv1(clipidx:end)];
    Rc = interp1(tv1,LSDfix.Rc,tv);
    SPhi = interp1(tv1,LSDfix.SPhi,tv);
    tv = tv - 2010 + samplingyr;
else; % assume 2010 value for all years >2010
    Rc = [LSDfix.Rc(1) LSDfix.Rc];
    SPhi = [LSDfix.SPhi(1) LSDfix.SPhi];
    tv = [0 (tv1 + samplingyr - 2010)];
end;

% calculate min and max depth (cm) and make shielding depth vector
mind = min(depth);
maxd = max(depth) + Tmax.*Emax.*1e-4;
dz = linspace(mind,maxd,100);

% muon production
fprintf(1,'calculating muon P...');
P_mu = P_mu_expage(dz.*rho,atm,LSDfix.RcEst,consts.SPhiInf,nucl10,nucl26,consts,'no');
fprintf(1,' done!\n');

% pick out Pmu if data exists
if nucl10 == 1; Pmu10 = P_mu.mu10 .* shield; end;
if nucl26 == 1; Pmu26 = P_mu.mu26 .* shield; end;

% spallation surface production scaling
Psp0 = LSDspal(atm,Rc,SPhi,LSDfix.w,nucl10,nucl26,consts);

% pick out Psp if data exists
if nucl10 == 1; Psp010 = Psp0.sp10 .* shield; end;
if nucl26 == 1; Psp026 = Psp0.sp26 .* shield; end;

% interpolate Lsp using CRONUScalc method (Sato 2008; Marrero et al. 2016)
Lsp = rawattenuationlength(atm,Rc);

% display iteration
fprintf(1,'iterating up to %.0f (number of E and T points):\n',ETgrid);

% loop for Tvect
for i = 1:numel(Tvect);
    % display i
    if i/10 == round(i/10);
        fprintf(1,'%.0f ',i);
    end;
    
    % fix tv2 and Lsp2
    clipidx = max(find(tv < Tvect(i)));
    tv2 = [tv(1:clipidx) Tvect(i)];
    Lsp2 = interp1(tv,Lsp,tv2);
    
    if nucl10 == 1;
        % interpolate Psp0 over time
        Psp010T = interp1(tv,Psp010,tv2);
    end;
    if nucl26 == 1;
        % interpolate Psp0 over time
        Psp026T = interp1(tv,Psp026,tv2);
    end;
    
    % fix depth matrix for tv2 and all E in Evect
    depthmT = bsxfun(@times,tv2,Evect');
    
    % loop for individual samples
    for j = 1:numel(depth)
        % add depth of sample
        depthm = depthmT + depth(j);
        
        if N10(j)+delN10(j) > 0;
            % add 1 to counter
            j10 = j10+1;
            
            % interpolate Pmu over time for all E
            Pmu10m = interp1(dz',Pmu10',depthm','pchip');
            
            % calculate Psp over time for all E
            Psp10m = repmat(Psp010T.*Pref10,numel(Evect),1)' .* ...
                exp(-rho.*depthm./repmat(Lsp2,numel(Evect),1))';
            
            % decay matrix
            dcf10m = repmat(exp(-tv2.*l10),numel(Evect),1)';
            
            % calculate N for all E
            N10T = trapz(tv2',Psp10m.*dcf10m + Pmu10m.*dcf10m);
            
            % calculate N diff for all E
            N10diff(j10,:) = ((N10(j)-N10T) ./ delN10(j)).^2;
        end;
        
        if N26(j)+delN26(j) > 0;
            % add 1 to counter
            j26 = j26+1;
            
            % interpolate Pmu over time for all E
            Pmu26m = interp1(dz',Pmu26',depthm','pchip');
            
            % calculate Psp over time for all E
            Psp26m = repmat(Psp026T.*Pref26,numel(Evect),1)' .* ...
                exp(-rho.*depthm./repmat(Lsp2,numel(Evect),1))';
            
            % decay matrix
            dcf26m = repmat(exp(-tv2.*l26),numel(Evect),1)';
            
            % calculate N for all E
            N26T = trapz(tv2',Psp26m.*dcf26m + Pmu26m.*dcf26m);
            
            % calculate N diff for all E
            N26diff(j26,:) = (N26(j)-N26T).^2 ./ delN26(j).^2;
        end;
    end;
    
    % calculate Rchi2 for all E and clear Ndiff
    if nucl10 == 1;
        Rchi2m10(i,:) = sum(N10diff) ./ (n10-1);
        clear N10diff;
    end;
    if nucl26 == 1;
        Rchi2m26(i,:) = sum(N26diff) ./ (n26-1);
        clear N26diff;
    end;
end;

% fix mm/ka
Evect = Evect .* 1e4;

fprintf(1,'\nT-range: %.0f - %.0f yr\n',Tmin,Tmax)
fprintf(1,'E-range: %.1f - %.1f mm/ka\n',Emin,Emax)

% contour line reduced chi-square values to add to minimum reduced chi-square value
cl = [0.05 0.5 1 2];

% if 10Be exists
if nucl10 == 1;
    % fix parameters
    Rchi2m = Rchi2m10;
    n = n10;
    nstr1 = '10Be';
    nstr2 = '^{10}Be';
    
    % do plotting
    plot_ET(Rchi2m,n,Evect,Tvect,xx,yy,cl,ETlims,nstr1,nstr2);
end;

% if 26Al exists
if nucl26 == 1;
    % fix parameters
    Rchi2m = Rchi2m26;
    n = n26;
    nstr1 = '26Al';
    nstr2 = '^{26}Al';
    
    % do plotting
    plot_ET(Rchi2m,n,Evect,Tvect,xx,yy,cl,ETlims,nstr1,nstr2);
end;

toc()


% subfunction plot_ET ==============================================================================
function plot_ET(Rchi2m,n,Evect,Tvect,xx,yy,cl,ETlims,nstr1,nstr2);
    % interpolate to finer mesh with spline
    Rchi2mi = interp2(Evect,Tvect,Rchi2m,xx,yy,'spline');
    
    % find index and minimum chi square
    [idxT,idxE] = find(Rchi2mi == min(min(Rchi2mi)));
    Tbest = yy(idxT);
    Ebest = xx(idxE);
    Rchisq = Rchi2mi(idxT,idxE);
    Pvalue = 1 - chi2cdf(Rchisq.*(n-1),n-1);
    
    % display output
    fprintf(1,'%s   E = %.1f mm/ka   T = %.0f yr   Rchi2 = %.2f   P-value = %.2f\n',...
        nstr1,Ebest,Tbest,Rchisq,Pvalue);
    
    % contour lines - display lines for Rchisq plus values in vector cl
    clp = cl + Rchisq;
    
    % plotting ============================================
    figure;
    hold on;
    box on;
    plot(Ebest,Tbest,'.','color','black','markersize',15);
    contour(xx,yy,Rchi2mi,clp);
    axis(ETlims,'square');
    xlabel([nstr2,' erosion (mm/ka)']);
    ylabel([nstr2,' exposure duration (yr)']);
    hold off;
% end subfunction plot_ET ==========================================================================
