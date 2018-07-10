function depthcalc()

% Depth profile 10Be/26Al exposure age calculator calculating the best fit depth profile exposure
% age. The calculation is done using the expage calculator production rates (nuclide-specific LSD
% scaling). At least two concentrations must be given for the calculator to work.
% This is free software: you can use/copy/modify/distribute as long as you keep it free/open.
% Jakob Heyman - 2015-2018 (jakob.heyman@gu.se)

clear;
close all;
tic();

% What version is this?
ver = '201806';

% plotting? (1 = yes) ==============================================================================
plotting = 1; % disabling plotting speeds up the function
% ==================================================================================================

% read input file
% NOTE! sample thickness in standard input is here changed to sample mid-point depth (cm)
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
    if (strcmp(aa(1),'pre'));
        diffstr = [diffstr,' atm-pressure'];
    else;
        diffstr = [diffstr,' elevation'];
    end;
end;
if sum(rho ~= rho(1)) > 0;
    diffstr = [diffstr,' density'];
end;
if sum(erosion ~= erosion(1)) > 0;
    diffstr = [diffstr,' erosion'];
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
erosion = mean(erosion);
samplingyr = mean(samplingyr);
% ==================================================================================================

% use mean atmospheic pressure
atm = mean(pressure);

% Set nucl to 0 for both 10/26
nucl10 = 0; nucl26 = 0;

% Check if 10Be and 26Al is measured
n10test = find(N10+delN10 > 0);
n26test = find(N26+delN26 > 0);
N10 = N10(n10test);
N26 = N26(n26test);
delN10 = delN10(n10test);
delN26 = delN26(n26test);
dv10 = depth(n10test);
dv26 = depth(n26test);
sample10 = sample_name(n10test);
sample26 = sample_name(n26test);
if numel(N10)>1; nucl10 = 1; end;
if numel(N26)>1; nucl26 = 1; end;

% Age Relative to t0=2010 - LSD tv from LSDfix
% tv = [0:10:50 60:100:50060 51060:1000:2000060 logspace(log10(2001060),7,200)];

% Fix w,Rc,SPhi, for sp and mu prod rate scaling
LSDfix = LSD_fix(lat,long,1e7,-1,consts);

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

% muon production
if erosion > 0;
    % depth vector for Pmu calculation
    derosion = linspace(min(depth),1e7.*erosion + max(depth),50);
    
    % calculate Pmu for depth vector
    P_mu = P_mu_LSD(derosion.*rho,atm,LSDfix.RcEst,consts.SPhiInf,nucl10,nucl26,consts,'no');
    
    % interpolate Pmu for individual samples
    if nucl10 == 1;
        for i = 1:numel(dv10);
            Pmu10(i,:) = interp1(derosion,P_mu.Be,dv10(i)+tv.*erosion,'pchip') .* shield;
        end;
    end;
    if nucl26 == 1;
        for i = 1:numel(dv26);
            Pmu26(i,:) = interp1(derosion,P_mu.Al,dv26(i)+tv.*erosion,'pchip') .* shield;
        end;
    end;
else; % no erosion
    P_mu = P_mu_LSD(depth.*rho,atm,LSDfix.RcEst,consts.SPhiInf,nucl10,nucl26,consts,'no');
    
    % pick out Pmu if data exists
    if nucl10 == 1; Pmu10 = P_mu.Be'(n10test) .* shield; end;
    if nucl26 == 1; Pmu26 = P_mu.Al'(n26test) .* shield; end;
end;

% spallation surface production scaling
Psp0 = LSDspal(atm,Rc,SPhi,LSDfix.w,nucl10,nucl26,consts);

% interpolate Lsp using CRONUScalc method (Sato 2008; Marrero et al. 2016)
Lsp = rawattenuationlength(atm,Rc);

% pick out Psp if data exists
if nucl10 == 1;
    for i = 1:numel(N10);
        Psp10(i,:) = Psp0.Be .* Pref10 .* shield .* exp(-rho.*(tv.*erosion+dv10(i))./Lsp);
    end;
end;
if nucl26 == 1;
    for i = 1:numel(N26);
        Psp26(i,:) = Psp0.Al .* Pref26 .* shield .* exp(-rho.*(tv.*erosion+dv26(i))./Lsp);
    end;
end;

% age calculation
if nucl10 == 1; % if 10Be exists
    % fix parameters for age calculation
    l = l10;
    N = N10;
    delN = delN10;
    Psp = Psp10;
    Pmu = Pmu10;
    sample = sample10;
    Pref = Pref10;
    delPref = delPref10;
    nstr = '10Be';
    
    % calculate best dpeth profile age
    [age10,ageunc_int10,ageunc_ext10] = ...
        get_depthage(tv,l,N,delN,Psp,Pmu,Lsp,sample,Pref,delPref,nstr);
end;

if nucl26 == 1; % if 26Al exists
    % fix parameters for age calculation
    l = l26;
    N = N26;
    delN = delN26;
    Psp = Psp26;
    Pmu = Pmu26;
    sample = sample26;
    Pref = Pref26;
    delPref = delPref26;
    nstr = '26Al';
    
    % calculate best dpeth profile age
    [age26,ageunc_int26,ageunc_ext26] = ...
        get_depthage(tv,l,N,delN,Psp,Pmu,Lsp,sample,Pref,delPref,nstr);
end;



% plotting =========================================================================================
if plotting == 1 && nucl10+nucl26 >= 1;
    % make depth vect for plotting
    dplot = linspace(0,max(depth).*1.25,50);
    
    if erosion > 0;
        % pick out oldest age of 10/26
        agetest = [];
        if nucl10 == 1; agetest(end+1) = age10 + age10unc_ext; end;
        if nucl26 == 1; agetest(end+1) = age26 + age26unc_ext; end;
        agemax = max(agetest);
        
        %make depth vect for mu
        dplotmu = linspace(0,max(depth).*1.25 + agemax.*erosion,50);
    else; % if no erosion
        % make depth vect for mu
        dplotmu = dplot;
    end;
    
    % muon production for depth profile points
    fprintf('calculating depth profile P from muons...');
    Pmuplot = P_mu_LSD(dplotmu.*rho,atm,LSDfix.RcEst,consts.SPhiInf,nucl10,nucl26,consts,'no');
    fprintf(' done!\n');
    
    if nucl10 == 1; % if 10Be exists and plotting = 1
        % fix parameters for plotting
        N = N10;
        delN = delN10;
        l = l10;
        Psp0n = Psp0.Be;
        Pmun = Pmuplot.Be;
        dv = dv10;
        age = age10;
        ageunc_int = ageunc_int10;
        ageunc_ext = ageunc_ext10;
        Pref = Pref10;
        nstr = '^{10}Be';
        
        % do plotting
        depthplot(tv,l,Psp0n,Pmun,Lsp,rho,shield,erosion,dplot,dplotmu,N,delN,dv,...
            age,ageunc_int,ageunc_ext,Pref,nstr);
    end;
    
    if nucl26 == 1; % if 26Al exists and plotting = 1
        % fix parameters for plotting
        N = N26;
        delN = delN26;
        l = l26;
        Psp0n = Psp0.Al;
        Pmun = Pmuplot.Al;
        dv = dv26;
        age = age26;
        ageunc_int = ageunc_int26;
        ageunc_ext = ageunc_ext26;
        Pref = Pref26;
        nstr = '^{26}Al';
        
        % do plotting
        depthplot(tv,l,Psp0n,Pmun,Lsp,rho,shield,erosion,dplot,dplotmu,N,delN,dv,...
            age,ageunc_int,ageunc_ext,Pref,nstr);
    end;
end; % end plotting
% ==================================================================================================

toc()


% subfunction get_depthage =========================================================================
function [age,ageunc_int,ageunc_ext] = ...
    get_depthage(tv,l,N,delN,Psp,Pmu,Lsp,sample,Pref,delPref,nstr);
    dcf = exp(-tv.*l); % decay factor;
    for i = 1:numel(N);
        Ntv(i,:) = cumtrapz(tv,(Psp(i,:).*dcf + Pmu(i,:).*dcf)); % potential N back in time
        % test saturation
        if N(i) <= max(Ntv(i,:)); % if not saturated
            tt(i) = interp1(Ntv(i,:),tv,N(i)); % individual sample exposure age
        else;
            fprintf(1,'Sample %s appears to be saturated in %s!\n',sample{i},nstr);
            tt(i) = 1E7;
        end;
        
        % uncertainty estimate
        % A with integrated average Lsp
        if tt(i) > 0;
            Lsp_avg = interp1(tv,cumtrapz(tv,Lsp.*exp(-l.*tv)),min(tt(i),max(tv)))/...
                interp1(tv,cumtrapz(tv,exp(-l.*tv)),min(tt(i),max(tv)));
            A = l./Lsp_avg;
            FP = (N(i).*A)./(1 - exp(-A.*tt(i)));
            dtdN = 1./(FP - N(i).*A);
            deltt(i) = sqrt(dtdN.^2 * delN(i).^2);
        else; % t = 0, estimate uncertainty based on conc + unc
            deltt(i) = interp1(Ntv(i,:),tv,N(i)+delN(i));
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
    delNm = repmat(delN,1,numel(tvint));
    chi2v = sum((Ntvm-Nm).^2 ./ delNm.^2);
    
    % find minimum chi square and index
    [chimin,idx] = min(chi2v);
    
    % find best fit age (depth profile age)
    age = tvint(idx);
    
    % estimate depth profile age uncertainty
    ageunc_int = sqrt(1./sum(deltt.^-2) .* 1./(numel(N)-1) .* sum((age-tt).^2./deltt.^2));
    ageunc_ext = sqrt(ageunc_int.^2 + (delPref./Pref.*age).^2);
    
    Rchi2 = chimin./(numel(N)-1);
    Pvalue = 1 - chi2cdf(chimin,numel(N)-1);
    
    fprintf(1,'%s age = %.0f ± %.0f (%.0f) yr    R-chisq = %.3f    P-value = %.3f\n',...
        nstr,age,ageunc_ext,ageunc_int,Rchi2,Pvalue);
% end subfunction get_depthage =====================================================================


% subfunction depthplot ============================================================================
function depthplot(tv,l,Psp0n,Pmun,Lsp,rho,shield,erosion,dplot,dplotmu,N,delN,dv,...
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
            exp(-rho.*(tvplot_tt.*erosion+dplot(i))./Lsp_tt);
        Psp_intmin(i,:) = Psp0_intmin .* Pref .* shield .* ...
            exp(-rho.*(tvplot_intmin.*erosion+dplot(i))./Lsp_intmin);
        Psp_intmax(i,:) = Psp0_intmax .* Pref .* shield .* ...
            exp(-rho.*(tvplot_intmax.*erosion+dplot(i))./Lsp_intmax);
        Psp_extmin(i,:) = Psp0_extmin .* Pref .* shield .* ...
            exp(-rho.*(tvplot_extmin.*erosion+dplot(i))./Lsp_extmin);
        Psp_extmax(i,:) = Psp0_extmax .* Pref .* shield .* ...
            exp(-rho.*(tvplot_extmax.*erosion+dplot(i))./Lsp_extmax);
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
    patch(Nextunc_x,N_y,'facecolor',[0.9 0.9 0.9],'EdgeColor','none');
    patch(Nintunc_x,N_y,'facecolor',[0.75 0.75 0.75],'EdgeColor','none');
    
    % plot standard depth profile
    plot(Nplot_tt,dplot,'color','black');
    plot(N,dv,'.','color','red','markersize',15);
    for i = 1:numel(N);
        plot([N(i)-delN(i),N(i)+delN(i)],[dv(i),dv(i)],'r');
    end;
    
    axis('ij');
    set(gca,'xaxislocation','top');
    ylim([0 max(dplot)])
    xlabel([nstr,' (atoms/g)']);
    ylabel('Depth (cm)');
    set(gca,'layer','top'); % plot axis on top
    hold off;
% end subfunction depthplot ========================================================================