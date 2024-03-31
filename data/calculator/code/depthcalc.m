function depthcalc()

% Depth profile 10Be/26Al exposure age calculator calculating the best-fit depth profile exposure
% age calculated using relative concentration uncertainties (absolute unertainties can be used - see
% lines 18-19 below). The calculation is done using the expage calculator production rates (nuclide-
% specific LSD scaling). At least two concentrations must be given for the calculator to work.
% This is free software: you can use/copy/modify/distribute as long as you keep it free/open.
% Jakob Heyman - 2015-2024 (jakob.heyman@gu.se)

clear;
close all;
tic();

% What version is this?
ver = '202403';

% make choices here ================================================================================
% use absolute concentration uncertainties for depth profile matching (1 = yes)
absunc = 0;

% plotting? (1 = yes)
plotting = 1; % disabling plotting speeds up the function
% ==================================================================================================

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

% run and load expage constants
make_consts_expage;
load consts_expage;

% Prefs
Pref10 = consts.Pref10; Pref10unc = consts.Pref10unc;
Pref26 = consts.Pref26; Pref26unc = consts.Pref26unc;
% decay constant
l10 = consts.l10; l10unc = consts.l10unc;
l26 = consts.l26; l26unc = consts.l26unc;

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

% if there is no erosion in input: assume zero erosion
if isfield(samplein,'erosion'); erosion = samplein.erosion; else; erosion = 0; end;

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

% fix erosion rate unit (mm/ka -> cm/yr)
erosion = erosion .* 1E-4;

% check and fix inputs =============================================================================
diffstr = '';

if sum(lat ~= lat(1)) > 0;
    diffstr = [diffstr,' latitude'];
end;
if sum(long ~= long(1)) > 0;
    diffstr = [diffstr,' longitude'];
end;
if sum(elv ~= elv(1)) > 0;
    if (strcmp(Pflag(1),'pre'));
        diffstr = [diffstr,' atm-pressure'];
    else;
        diffstr = [diffstr,' elevation'];
    end;
end;
if sum(dens ~= dens(1)) > 0;
    diffstr = [diffstr,' density'];
end;
if sum(erosion ~= erosion(1)) > 0;
    diffstr = [diffstr,' erosion'];
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
shield = mean(shield);
dens = mean(dens);
erosion = mean(erosion);
samplingyr = mean(samplingyr);
% ==================================================================================================

% use mean atmospheic pressure
atm = mean(pressure);

% Set nucl to 0 for both 10/26
nucl10 = 0; nucl26 = 0;

% Check if 10Be and 26Al is measured
n10test = find(N10+N10unc > 0);
n26test = find(N26+N26unc > 0);
N10 = N10(n10test);
N26 = N26(n26test);
N10unc = N10unc(n10test);
N26unc = N26unc(n26test);
dv10 = depth(n10test);
dv26 = depth(n26test);
sample10 = sample(n10test);
sample26 = sample(n26test);
if numel(N10)>1; nucl10 = 1; end;
if numel(N26)>1; nucl26 = 1; end;

% Age Relative to t0=2010 - LSD tv from LSDfix
% tv = [0:10:50 60:100:2960 3060:200:74860 75060:1000:799060 800060:2000:2000060 1E7];

% Fix tv, Rc, RcEst, SPhi, and w for sp and mu prod rate scaling
LSDfix = LSD_fix(lat,long,1E7,-1,samplingyr,consts);

% fix variables
tv = LSDfix.tv;
Rc = LSDfix.Rc;
SPhi = LSDfix.SPhi;

% muon production
if erosion > 0;
    % depth vector for Pmu calculation
    derosion = linspace(min(depth),1e7.*erosion + max(depth) + 1,50);
    
    % calculate Pmu for depth vector
    Pmu = P_mu_expage(derosion.*dens,atm,LSDfix.RcEst,consts.SPhiInf,nucl10,nucl26,0,consts,'no');
    
    % interpolate Pmu for individual samples
    if nucl10 == 1;
        for i = 1:numel(dv10);
            Pmu10(i,:) = interp1(derosion,Pmu.mu10,dv10(i)+tv.*erosion,'pchip') .* shield;
        end;
    end;
    if nucl26 == 1;
        for i = 1:numel(dv26);
            Pmu26(i,:) = interp1(derosion,Pmu.mu26,dv26(i)+tv.*erosion,'pchip') .* shield;
        end;
    end;
else; % no erosion
    Pmu = P_mu_expage(depth.*dens,atm,LSDfix.RcEst,consts.SPhiInf,nucl10,nucl26,0,consts,'no');
    
    % pick out Pmu if data exists
    if nucl10 == 1; Pmu10 = Pmu.mu10(n10test)' .* shield; end;
    if nucl26 == 1; Pmu26 = Pmu.mu26(n26test)' .* shield; end;
end;

% spallation production scaling
Psp0 = P_sp_expage(atm,Rc,SPhi,LSDfix.w,consts,nucl10,nucl26,0);

% interpolate Lsp using CRONUScalc method (Sato 2008; Marrero et al. 2016)
Lsp = rawattenuationlength(atm,Rc);

% pick out Psp if data exists
if nucl10 == 1;
    for i = 1:numel(N10);
        Psp10(i,:) = Psp0.sp10 .* Pref10 .* shield .* exp(-dens.*(tv.*erosion+dv10(i))./Lsp);
    end;
end;
if nucl26 == 1;
    for i = 1:numel(N26);
        Psp26(i,:) = Psp0.sp26 .* Pref26 .* shield .* exp(-dens.*(tv.*erosion+dv26(i))./Lsp);
    end;
end;

% age calculation
if nucl10 == 1; % if 10Be exists
    % fix parameters for age calculation
    l = l10;
    N = N10;
    Nunc = N10unc;
    Psp = Psp10;
    Pmu = Pmu10;
    sample = sample10;
    Pref = Pref10;
    Prefunc = Pref10unc;
    nstr = '10Be';
    
    % calculate best dpeth profile age
    [age10,ageunc_int10,ageunc_ext10] = ...
        get_depthage(tv,l,N,Nunc,Psp,Pmu,Lsp,sample,Pref,Prefunc,1E7,absunc,nstr);
end;

if nucl26 == 1; % if 26Al exists
    % fix parameters for age calculation
    l = l26;
    N = N26;
    Nunc = N26unc;
    Psp = Psp26;
    Pmu = Pmu26;
    sample = sample26;
    Pref = Pref26;
    Prefunc = Pref26unc;
    nstr = '26Al';
    
    % calculate best dpeth profile age
    [age26,ageunc_int26,ageunc_ext26] = ...
        get_depthage(tv,l,N,Nunc,Psp,Pmu,Lsp,sample,Pref,Prefunc,6E6,absunc,nstr);
end;


% plotting =========================================================================================
if plotting == 1 && nucl10+nucl26 >= 1;
    % make depth vect for plotting
    dplot = linspace(0,max(depth).*1.25,50);
    
    if erosion > 0;
        % pick out oldest age of 10/26
        agetest = [];
        if nucl10 == 1; agetest(end+1) = age10 + ageunc_ext10; end;
        if nucl26 == 1; agetest(end+1) = age26 + ageunc_ext26; end;
        agemax = max(agetest);
        
        %make depth vect for mu
        dplotmu = linspace(0,max(depth).*1.25 + agemax.*erosion,50);
    else; % if no erosion
        % make depth vect for mu
        dplotmu = dplot;
    end;
    
    % muon production for depth profile points
    fprintf('calculating depth profile P from muons...');
    Pmuplot = P_mu_expage(dplotmu.*dens,atm,LSDfix.RcEst,consts.SPhiInf,nucl10,nucl26,0,consts,'no');
    fprintf(' done!\n');
    
    if nucl10 == 1; % if 10Be exists and plotting = 1
        % fix parameters for plotting
        N = N10;
        Nunc = N10unc;
        l = l10;
        Psp0n = Psp0.sp10;
        Pmun = Pmuplot.mu10;
        dv = dv10;
        age = age10;
        ageunc_int = ageunc_int10;
        ageunc_ext = ageunc_ext10;
        Pref = Pref10;
        nstr = '^{10}Be';
        
        % do plotting
        depthplot(tv,l,Psp0n,Pmun,Lsp,dens,shield,erosion,dplot,dplotmu,N,Nunc,dv,...
            age,ageunc_int,ageunc_ext,Pref,nstr);
    end;
    
    if nucl26 == 1; % if 26Al exists and plotting = 1
        % fix parameters for plotting
        N = N26;
        Nunc = N26unc;
        l = l26;
        Psp0n = Psp0.sp26;
        Pmun = Pmuplot.mu26;
        dv = dv26;
        age = age26;
        ageunc_int = ageunc_int26;
        ageunc_ext = ageunc_ext26;
        Pref = Pref26;
        nstr = '^{26}Al';
        
        % do plotting
        depthplot(tv,l,Psp0n,Pmun,Lsp,dens,shield,erosion,dplot,dplotmu,N,Nunc,dv,...
            age,ageunc_int,ageunc_ext,Pref,nstr);
    end;
end; % end plotting
% ==================================================================================================

toc()


% subfunction get_depthage =========================================================================
function [age,ageunc_int,ageunc_ext] = ...
    get_depthage(tv,l,N,Nunc,Psp,Pmu,Lsp,sample,Pref,Prefunc,maxtt,absunc,nstr);
    dcf = exp(-tv.*l); % decay factor;
    for i = 1:numel(N);
        Ntv(i,:) = cumtrapz(tv,(Psp(i,:).*dcf + Pmu(i,:).*dcf)); % potential N back in time
        % test saturation
        if N(i) <= max(Ntv(i,:)); % if not saturated
            tt(i) = interp1(Ntv(i,:),tv,N(i)); % individual sample exposure age
            
            % uncertainty estimate
            % A with integrated average Lsp
            if tt(i) > 0;
                Lsp_avg = interp1(tv,cumtrapz(tv,Lsp.*exp(-l.*tv)),min(tt(i),max(tv)))/...
                    interp1(tv,cumtrapz(tv,exp(-l.*tv)),min(tt(i),max(tv)));
                A = l./Lsp_avg;
                FP = (N(i).*A)./(1 - exp(-A.*tt(i)));
                dtdN = 1./(FP - N(i).*A);
                deltt(i) = sqrt(dtdN.^2 * Nunc(i).^2);
            else; % t = 0, estimate uncertainty based on conc + unc
                deltt(i) = interp1(Ntv(i,:),tv,N(i)+Nunc(i));
            end;
        else;
            fprintf(1,'Sample %s appears to be saturated in %s!\n',sample{i},nstr);
            tt(i) = maxtt;
            deltt(i) = maxtt;
        end;
    end;
    
    % find min and max value for interpolation and make new time vector
    minage = round(min(tt) - (max(tt)-min(tt)));
    maxage = round(max(tt) + (max(tt)-min(tt)));
    if minage < 0; minage = 0; end;
    if maxage > 1e7; maxage = 1e7; end;
    tvint = (minage:1:maxage);
    
    % interpret N vectors
    for i = 1:numel(N);
        Ntvm(i,:) = interp1(tv,Ntv(i,:),tvint);
    end;
    
    % calculate chi square vector
    Nm = repmat(N,1,numel(tvint));
    Nuncm = repmat(Nunc,1,numel(tvint));
    chi2v = sum((Ntvm-Nm).^2 ./ Nuncm.^2);
    
    % calculate chi square vector for relative N uncertainties
    Nuncm_rel = repmat(Nunc./N,1,numel(tvint));
    chi2v_rel = sum((Ntvm-Nm).^2 ./ Nuncm_rel.^2);
    
    % find minimum chi square and index
    [chimin_rel,idx] = min(chi2v_rel);
    chimin = chi2v(idx);
    
    if absunc == 1;
        % find minimum chi square and index
        [chimin,idx] = min(chi2v);
    end;
    
    % find best fit age (depth profile age)
    age = tvint(idx);
    
    % estimate depth profile age uncertainty
    ageunc_int = sqrt(1./sum((deltt./tt).^-2) .* 1./(numel(N)-1) .* ...
        sum((age-tt).^2./(deltt./tt).^2));
    if absunc == 1;
        ageunc_int = sqrt(1./sum(deltt.^-2) .* 1./(numel(N)-1) .* sum((age-tt).^2./deltt.^2));
    end;
    ageunc_ext = sqrt(ageunc_int.^2 + (Prefunc./Pref.*age).^2);
    
    Rchi2 = chimin./(numel(N)-1);
    Pvalue = 1 - chi2cdf(chimin,numel(N)-1);
    
    fprintf(1,'%s age = %.0f ± %.0f (%.0f) yr    R-chisq = %.3f    P-value = %.3f\n',...
        nstr,age,ageunc_ext,ageunc_int,Rchi2,Pvalue);
% end subfunction get_depthage =====================================================================


% subfunction depthplot ============================================================================
function depthplot(tv,l,Psp0n,Pmun,Lsp,dens,shield,erosion,dplot,dplotmu,N,Nunc,dv,...
    age,ageunc_int,ageunc_ext,Pref,nstr);
    % clip Psp0 and Lsp for depth profile age
    % 1. actual age
    clipidx_tt = max(find(tv < age));
    tvplot_tt = [tv(1:clipidx_tt) age];
    Psp0_tt = interp1(tv,Psp0n,tvplot_tt);
    Lsp_tt = interp1(tv,Lsp,tvplot_tt);
    % 2. min age internal unc
    clipidx_intmin = max(find(tv < age-ageunc_int));
    tvplot_intmin = [tv(1:clipidx_intmin) age-ageunc_int];
    Psp0_intmin = interp1(tv,Psp0n,tvplot_intmin);
    Lsp_intmin = interp1(tv,Lsp,tvplot_intmin);
    % 3. max age internal unc
    clipidx_intmax = max(find(tv < age+ageunc_int));
    tvplot_intmax = [tv(1:clipidx_intmax) age+ageunc_int];
    Psp0_intmax = interp1(tv,Psp0n,tvplot_intmax);
    Lsp_intmax = interp1(tv,Lsp,tvplot_intmax);
    % 4. min age external unc
    clipidx_extmin = max(find(tv < age-ageunc_ext));
    tvplot_extmin = [tv(1:clipidx_extmin) age-ageunc_ext];
    Psp0_extmin = interp1(tv,Psp0n,tvplot_extmin);
    Lsp_extmin = interp1(tv,Lsp,tvplot_extmin);
    % 5. max age external unc
    clipidx_extmax = max(find(tv < age+ageunc_ext));
    tvplot_extmax = [tv(1:clipidx_extmax) age+ageunc_ext];
    Psp0_extmax = interp1(tv,Psp0n,tvplot_extmax);
    Lsp_extmax = interp1(tv,Lsp,tvplot_extmax);
    
    % calculate decay factors
    dcf_tt = exp(-tvplot_tt.*l);
    dcf_intmin = exp(-tvplot_intmin.*l);
    dcf_intmax = exp(-tvplot_intmax.*l);
    dcf_extmin = exp(-tvplot_extmin.*l);
    dcf_extmax = exp(-tvplot_extmax.*l);
    
    % pick out Pmu
    for i = 1:numel(dplot);
        % muon production for depth profile points
        Pmu_tt(i,:) = interp1(dplotmu,Pmun.*shield,dplot(i)+tvplot_tt.*erosion);
        Pmu_intmin(i,:) = interp1(dplotmu,Pmun.*shield,dplot(i)+tvplot_intmin.*erosion);
        Pmu_intmax(i,:) = interp1(dplotmu,Pmun.*shield,dplot(i)+tvplot_intmax.*erosion);
        Pmu_extmin(i,:) = interp1(dplotmu,Pmun.*shield,dplot(i)+tvplot_extmin.*erosion);
        Pmu_extmax(i,:) = interp1(dplotmu,Pmun.*shield,dplot(i)+tvplot_extmax.*erosion);
    end;
    
    % calculate Psp for depth profile points
    for i = 1:numel(dplot);
        Psp_tt(i,:) = Psp0_tt .* Pref .* shield .* ...
            exp(-dens.*(tvplot_tt.*erosion+dplot(i))./Lsp_tt);
        Psp_intmin(i,:) = Psp0_intmin .* Pref .* shield .* ...
            exp(-dens.*(tvplot_intmin.*erosion+dplot(i))./Lsp_intmin);
        Psp_intmax(i,:) = Psp0_intmax .* Pref .* shield .* ...
            exp(-dens.*(tvplot_intmax.*erosion+dplot(i))./Lsp_intmax);
        Psp_extmin(i,:) = Psp0_extmin .* Pref .* shield .* ...
            exp(-dens.*(tvplot_extmin.*erosion+dplot(i))./Lsp_extmin);
        Psp_extmax(i,:) = Psp0_extmax .* Pref .* shield .* ...
            exp(-dens.*(tvplot_extmax.*erosion+dplot(i))./Lsp_extmax);
    end;
    
    % calculate N after age yr
    for i = 1:numel(dplot);
        Nplot_tt(i,:) = trapz(tvplot_tt,(Psp_tt(i,:).*dcf_tt + Pmu_tt(i).*dcf_tt));
        Nplot_intmin(i,:) = trapz(tvplot_intmin,(Psp_intmin(i,:).*dcf_intmin ...
            + Pmu_intmin(i).*dcf_intmin));
        Nplot_intmax(i,:) = trapz(tvplot_intmax,(Psp_intmax(i,:).*dcf_intmax ...
            + Pmu_intmax(i).*dcf_intmax));
        Nplot_extmin(i,:) = trapz(tvplot_extmin,(Psp_extmin(i,:).*dcf_extmin ...
            + Pmu_extmin(i).*dcf_extmin));
        Nplot_extmax(i,:) = trapz(tvplot_extmax,(Psp_extmax(i,:).*dcf_extmax ...
            + Pmu_extmax(i).*dcf_extmax));
    end;
    
    % fix uncertainty areas
    Nintunc_x = [Nplot_intmin' flip(Nplot_intmax')];
    Nextunc_x = [Nplot_extmin' flip(Nplot_extmax')];
    N_y = [dplot flip(dplot)];
    
    figure;
    hold on;
    box on;
    
    % plot uncertainty areas
    patch(Nextunc_x,N_y,[0.9 0.9 0.9],'EdgeColor','none');
    patch(Nintunc_x,N_y,[0.75 0.75 0.75],'EdgeColor','none');
    
    % plot standard depth profile
    plot(Nplot_tt,dplot,'color','black');
    plot(N,dv,'.','color','red','markersize',15);
    for i = 1:numel(N);
        plot([N(i)-Nunc(i),N(i)+Nunc(i)],[dv(i),dv(i)],'color','red');
    end;
    
    axis('ij');
    set(gca,'xaxislocation','top');
    ylim([0 max(dplot)])
    xlabel([nstr,' (atoms/g)']);
    ylabel('Depth (cm)');
    set(gca,'layer','top'); % plot axis on top
    hold off;
% end subfunction depthplot ========================================================================
