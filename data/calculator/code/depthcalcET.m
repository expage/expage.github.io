function depthcalcET()

% Depth profile 10Be/26Al exposure age calculator calculating the best fit erosion rate E and
% exposure duration T within the range Emin-Emax/Tmin-Tmax (one period of exposure with constant
% erosion rate). The calculation is done using the expage calculator production rates (nuclide-
% specific LSD scaling) and assuming one period of exposure. At least two concentrations must be
% given for the calculator to work.
% This is free software: you can use/copy/modify/distribute as long as you keep it free/open.
% Jakob Heyman - 2015-2023 (jakob.heyman@gu.se)

clear;
close all;
tic();

% age and erosion ranges to test ===================================================================
Tmin = 200000; % yr (must be larger than 0!)
Tmax = 500000; % yr
Emin = 0; % mm/ka
Emax = 3; % mm/ka
% ==================================================================================================

% number of E and T points for the test grid
ETgrid = 100;

% What version is this?
ver = '202306';

% fix input ========================================================================================
% variable names for input with variable names in first line
varnames = {'sample','Pflag','std10','std26','isostP','lat','long','elv','depth','dens','shield',...
    'erosion','N10','N10unc','N26','N26unc','samplingyr','pressure'};
vartypes = {'%s','%s','%s','%s','%s','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n',...
    '%n'};
% read input file
fid = fopen('input.txt');
varsin = strsplit(fgetl(fid)); % read first line
if(ismember(varsin,varnames)); % if first line contain only variable names
    [testi,vari] = ismember(varsin,varnames); % find index of varnames
    typestr = vartypes{vari(1)}; % fix type string
    for i = 2:numel(vari); % fix type string
        typestr = [typestr ' ' vartypes{vari(i)}];
    end;
elseif numel(varsin) == 16; % if no variable names in first line
    frewind(fid); % read from first line
    varsin = {'sample','lat','long','elv','Pflag','depth','dens','shield','erosion','N10',...
        'N10unc','std10','N26','N26unc','std26','samplingyr'};
    typestr = '%s %n %n %n %s %n %n %n %n %n %n %s %n %n %s %n';
else;
    fprintf(1,'ERROR! Something is wrong with the input\n');
    return;
end;
indata = textscan(fid,typestr,'CommentStyle','%','TreatAsEmpty','-'); % scan data
for i = 1:numel(varsin); % fix variables
    samplein.(varsin{i}) = indata{i};
end;
fclose(fid);
% ==================================================================================================

% input fix
sample = samplein.sample;
lat = samplein.lat;
long = samplein.long;
depth = samplein.depth;
dens = samplein.dens;
shield = samplein.shield;
samplingyr = samplein.samplingyr;
if isfield(samplein,'elv');
    elv = samplein.elv;
    Pflag = samplein.Pflag;
end;

% run and load expage constants
make_consts_expage;
load consts_expage;

% Prefs
Pref10 = consts.Pref10; Pref10unc = consts.Pref10unc;
Pref26 = consts.Pref26; Pref26unc = consts.Pref26unc;
% decay constant
l10 = consts.l10; l10unc = consts.l10unc;
l26 = consts.l26; l26unc = consts.l26unc;

% if there is no N10 or N26 in input: fill with 0
if isfield(samplein,'N10'); N10 = samplein.N10; else; N10(1:numel(sample),1) = 0; end;
if isfield(samplein,'N10unc'); N10unc = samplein.N10unc; else; N10unc(1:numel(sample),1) = 0; end;
if isfield(samplein,'std10'); std10 = samplein.std10; else; std10(1:numel(sample),1) = {'0'}; end;
if isfield(samplein,'N26'); N26 = samplein.N26; else; N26(1:numel(sample),1) = 0; end;
if isfield(samplein,'N26unc'); N26unc = samplein.N26unc; else; N26unc(1:numel(sample),1) = 0; end;
if isfield(samplein,'std26'); std26 = samplein.std26; else; std26(1:numel(sample),1) = {'0'}; end;

% if there is NaN in N10 and N26: replace with 0
samplein.N10(isnan(samplein.N10)) = 0;
samplein.N10unc(isnan(samplein.N10unc)) = 0;
samplein.N26(isnan(samplein.N26)) = 0;
samplein.N26unc(isnan(samplein.N26unc)) = 0;

% convert 10Be concentrations according to standards
[testi,stdi] = ismember(std10,consts.std10); % find index of standard conversion factors
mult10 = consts.std10_cf(stdi); % pick out conversion factor
N10 = N10 .* mult10;
N10unc = N10unc .* mult10;

% convert 26Al concentrations according to standards
[testi,stdi] = ismember(std26,consts.std26); % find index of standard conversion factors
mult26 = consts.std26_cf(stdi); % pick out conversion factor
N26 = N26 .* mult26;
N26unc = N26unc .* mult26;

% fix longitude values
long(find(long < 0)) = long(find(long < 0)) + 360;

% fix sample pressure
if isfield(sample,'pressure') == 0;
    % if there is no pressure flag: use std
    if isfield(sample,'Pflag') == 0; Pflag(1:numel(sample),1) = {'std'}; end;
    stdv = strcmp(Pflag,'std');
    antv = strcmp(Pflag,'ant');
    prev = strcmp(Pflag,'pre');
    pressure(stdv) = ERA40atm(lat(stdv),long(stdv),elv(stdv));
    pressure(antv) = antatm(elv(antv));
    pressure(prev) = elv(prev);
else;
    pressure = samplein.pressure;
end;

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
if sum(dens ~= dens(1)) > 0;
    diffstr = [diffstr,' density'];
end;
if sum(shield ~= shield(1)) > 0;
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
shield = mean(shield);
dens = mean(dens);
samplingyr = mean(samplingyr);
% ==================================================================================================

% use mean atmospheic pressure
atm = mean(pressure);

% Set nucl to 0 for both 10/26
nucl10 = 0; nucl26 = 0;

% Check and count if 10Be and 26Al is measured
n10 = sum(N10+N10unc > 0);
n26 = sum(N26+N26unc > 0);
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
% tv = [0:10:50 60:100:2960 3060:200:74860 75060:1000:799060 800060:2000:2000060 1E7];

% Fix tv, Rc, RcEst, SPhi, and w for sp and mu prod rate scaling
LSDfix = LSD_fix(lat,long,Tmax,-1,samplingyr,consts);

% fix variables
tv = LSDfix.tv;
Rc = LSDfix.Rc;
SPhi = LSDfix.SPhi;

% calculate min and max depth (cm) and make shielding depth vector
mind = min(depth);
maxd = max(depth) + Tmax.*Emax.*1e-4;
dz = linspace(mind,maxd,100);

% muon production
fprintf(1,'calculating muon P...');
Pmu = P_mu_expage(dz.*dens,atm,LSDfix.RcEst,consts.SPhiInf,nucl10,nucl26,0,consts,'no');
fprintf(1,' done!\n');

% pick out Pmu if data exists
if nucl10 == 1; Pmu10 = Pmu.mu10 .* shield; end;
if nucl26 == 1; Pmu26 = Pmu.mu26 .* shield; end;

% spallation surface production scaling
Psp0 = P_sp_expage(atm,Rc,SPhi,LSDfix.w,consts,nucl10,nucl26,0);

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
        
        if N10(j)+N10unc(j) > 0;
            % add 1 to counter
            j10 = j10+1;
            
            % interpolate Pmu over time for all E
            Pmu10m = interp1(dz',Pmu10',depthm','pchip');
            
            % calculate Psp over time for all E
            Psp10m = repmat(Psp010T.*Pref10,numel(Evect),1)' .* ...
                exp(-dens.*depthm./repmat(Lsp2,numel(Evect),1))';
            
            % decay matrix
            dcf10m = repmat(exp(-tv2.*l10),numel(Evect),1)';
            
            % calculate N for all E
            N10T = trapz(tv2',Psp10m.*dcf10m + Pmu10m.*dcf10m);
            
            % calculate N diff for all E
            N10diff(j10,:) = ((N10(j)-N10T) ./ N10unc(j)).^2;
        end;
        
        if N26(j)+N26unc(j) > 0;
            % add 1 to counter
            j26 = j26+1;
            
            % interpolate Pmu over time for all E
            Pmu26m = interp1(dz',Pmu26',depthm','pchip');
            
            % calculate Psp over time for all E
            Psp26m = repmat(Psp026T.*Pref26,numel(Evect),1)' .* ...
                exp(-dens.*depthm./repmat(Lsp2,numel(Evect),1))';
            
            % decay matrix
            dcf26m = repmat(exp(-tv2.*l26),numel(Evect),1)';
            
            % calculate N for all E
            N26T = trapz(tv2',Psp26m.*dcf26m + Pmu26m.*dcf26m);
            
            % calculate N diff for all E
            N26diff(j26,:) = (N26(j)-N26T).^2 ./ N26unc(j).^2;
        end;
    end;
    
    % calculate Rchi2 for all E and clear Ndiff
    if nucl10 == 1;
        Rchi2m10(i,:) = sum(N10diff) ./ (n10-1);
        clear N10diff; j10 = 0;
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
