function expage()

% expage 10Be and 26Al exposure age calculator.
% Read and fix input, use get_1026_age to calculate exposure rates, and fix and write output.
%
% Based on code written by Greg Balco for the CRONUS calculator v. 2.
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License, version 2, as published by the Free Software Foundation (www.fsf.org).
%
% Jakob Heyman - 2015-2018 (jakob.heyman@gu.se)

clear;
tic();

% What version is this?
ver = '201806';

% plotting? (1 = yes) ==============================================================================
plotpointages = 0; % plots exposure ages as points
plotprobdens = 0; % plots exposure ages as probability density curves (for single groups of ages)
% ==================================================================================================

% read input file
[samplein.sample_name,samplein.lat,samplein.long,samplein.elv,samplein.aa,samplein.thick,...
    samplein.rho,samplein.othercorr,samplein.E,samplein.N10,samplein.delN10,samplein.be_stds,...
    samplein.N26,samplein.delN26,samplein.al_stds,samplein.samplingyr] = ...
    textread('input.txt','%s %n %n %n %s %n %n %n %n %n %n %s %n %n %s %n','commentstyle','matlab');

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
for i = 1:numel(samplein.N10);
    be_mult(i,1) = consts.be_stds_cfs(strcmp(samplein.be_stds(i),consts.be_stds_names));
end;
samplein.N10 = samplein.N10 .* be_mult;
samplein.delN10 = samplein.delN10 .* be_mult;

% convert 26Al concentrations according to standards
for i = 1:numel(samplein.N26);
    al_mult(i,1) = consts.al_stds_cfs(strcmp(samplein.al_stds(i),consts.al_stds_names));
end;
samplein.N26 = samplein.N26 .* al_mult;
samplein.delN26 = samplein.delN26 .* al_mult;

% fix longitude values
samplein.long(find(samplein.long < 0)) = samplein.long(find(samplein.long < 0)) + 360;

% fix sample pressure
std_v = strcmp(samplein.aa,'std');
ant_v = strcmp(samplein.aa,'ant');
pre_v = strcmp(samplein.aa,'pre');
samplein.pressure(std_v) = ERA40atm(samplein.lat(std_v),samplein.long(std_v),samplein.elv(std_v));
samplein.pressure(ant_v) = antatm(samplein.elv(ant_v));
samplein.pressure(pre_v) = samplein.elv(pre_v);

% fix erosion rate unit (mm/ka -> cm/yr)
samplein.E = samplein.E .* 1E-4;

% initial Lsp for simple age calculation
Lsp1 = consts.Lsp;

% fix simple thickness scaling factor
thick_v = find(samplein.thick > 0);
samplein.thickSF1(1:numel(samplein.thick),1) = 1;
samplein.thickSF1(thick_v) = (Lsp1./(samplein.rho(thick_v).*samplein.thick(thick_v))).*...
    (1 - exp(((-1.*samplein.rho(thick_v).*samplein.thick(thick_v))./Lsp1)));

% fix for output
output(1,1) = {'sample'};
outn(1) = 1;
if sum(samplein.N10+samplein.delN10) > 0;
    output(1,end+1:end+3) = {'10age(yr)','10uncext(yr)','10uncint(yr)'};
    outn(1) = max(outn)+1;
    outn(2) = max(outn)+2;
end;
if sum(samplein.N26+samplein.delN26) > 0;
    output(1,end+1:end+3) = {'26age(yr)','26uncext(yr)','26uncint(yr)'};
    outn(3) = max(outn)+1;
    outn(4) = max(outn)+2;
end;

% pick out samples one by one
for i = 1:numel(samplein.lat);
    sample.sample_name = samplein.sample_name(i);
    sample.lat = samplein.lat(i);
    sample.long = samplein.long(i);
    sample.thick = samplein.thick(i);
    sample.rho = samplein.rho(i);
    sample.othercorr = samplein.othercorr(i);
    sample.E = samplein.E(i);
    sample.N10 = samplein.N10(i);
    sample.delN10 = samplein.delN10(i);
    sample.be_stds = samplein.be_stds(i);
    sample.N26 = samplein.N26(i);
    sample.delN26 = samplein.delN26(i);
    sample.al_stds = samplein.al_stds(i);
    sample.samplingyr = samplein.samplingyr(i);
    sample.pressure = samplein.pressure(i);
    sample.thickSF1 = samplein.thickSF1(i);
    
    % write sample name to output
    output(i+1,1) = sample.sample_name;
    
    % Set nucl and mt to 0 for both 10/26 and check if there is N10/N26
    nucl10 = 0; nucl26 = 0; mt10 = 0; mt26 = 0;
    if (sample.N10 + sample.delN10) > 0; nucl10 = 1; end;
    if (sample.N26 + sample.delN26) > 0; nucl26 = 1; end;
    
    if nucl10 + nucl26 == 0;
        continue;
    end;
    
    % display sample name
    fprintf(1,'%.0f. %s',i,sample.sample_name{1});
    
    % Find P scaling factor according to Stone/Lal
    P_St_SF = stone2000(sample.lat,sample.pressure,1) * sample.thickSF1 * sample.othercorr;
    
    % if 10Be measured: calculate max time
    if nucl10 == 1;
        [tsimple10,mt10] = get_mt(sample,Pref10,P_St_SF,l10,Lsp1,sample.N10);
    end;
    
    % if 26Al measured: calculate max time
    if nucl26 == 1;
        [tsimple26,mt26] = get_mt(sample,Pref26,P_St_SF,l26,Lsp1,sample.N26);
    end;
    
    % pick largest of mt10 and mt26 as max time
    mt = max(mt10,mt26);
    
    % Age Relative to t0=2010 - LSD tv from LSDfix
    % tv = [0:10:50 60:100:50060 51060:1000:2000060 logspace(log10(2001060),7,200)];
    
    % Fix w,Rc,SPhi, for sp and mu prod rate scaling
    LSDfix = LSD_fix(sample.lat,sample.long,mt,-1,consts);
    
    % time vector tv1
    tv1 = LSDfix.tv;
    
    % adjust tv, Rc, and SPhi to sampling year
    if sample.samplingyr <= 2010;
        clipidx = min(find(tv1 > 2010-sample.samplingyr));
        tv = [2010-sample.samplingyr tv1(clipidx:end)];
        Rc = interp1(tv1,LSDfix.Rc,tv);
        SPhi = interp1(tv1,LSDfix.SPhi,tv);
        tv = tv - 2010 + sample.samplingyr;
    else; % assume 2010 value for all years >2010
        Rc = [LSDfix.Rc(1) LSDfix.Rc];
        SPhi = [LSDfix.SPhi(1) LSDfix.SPhi];
        tv = [0 (tv1 + sample.samplingyr - 2010)];
    end;
    
    % Production from muons
    if sample.E <= 0;
        P_mu = P_mu_LSD(sample.thick.*sample.rho./2,sample.pressure,LSDfix.RcEst,...
            consts.SPhiInf,nucl10,nucl26,consts,'no');
        if nucl10 == 1; sample.mu10 = P_mu.Be .* sample.othercorr; end;
        if nucl26 == 1; sample.mu26 = P_mu.Al .* sample.othercorr; end;
    else;
        tv_z = (tv.*sample.E + sample.thick./2) .* sample.rho; % time - depth vect (g/cm^2)
        if nucl10 == 1;
            sample.mu10 = get_PmuE(sample,tv_z,tsimple10,LSDfix.RcEst,consts,1,0);
        end;
        if nucl26 == 1;
            sample.mu26 = get_PmuE(sample,tv_z,tsimple26,LSDfix.RcEst,consts,0,1);
        end;
    end;
    
    % spallation production scaling
    LSDnu = LSDspal(sample.pressure,Rc,SPhi,LSDfix.w,nucl10,nucl26,consts);
    
    % interpolate Lsp using CRONUScalc method (Sato 2008; Marrero et al. 2016)
    sample.Lsp = rawattenuationlength(sample.pressure,Rc);
    
    % Thickness scaling factor.
    if sample.thick > 0;
        thickSF = (sample.Lsp./(sample.rho.*sample.thick)).*...
            (1 - exp(((-1.*sample.rho.*sample.thick)./sample.Lsp)));
    else;
        thickSF = 1;
    end;
    
    % include in sample
    sample.tv = tv;
    sample.dpfs = exp(-tv.*sample.E.*sample.rho./sample.Lsp); % spallation depth dependence
    
    if nucl10 == 1;
        % sample surface spallation production rate over time
        sample.sp = LSDnu.Be.*Pref10.*thickSF.*sample.othercorr;
        
        % sample muon P
        sample.mu = sample.mu10;
        
        % various parameters
        sample.N = sample.N10; sample.delN = sample.delN10;
        sample.Pref = Pref10; sample.delPref = delPref10;
        sample.l = l10;
        
        % get age
        results = get_1026_age(sample,1E7);
        
        % fill output
        output(i+1,outn(1):outn(2)) = results.outstr;
        
        % fill plot matrix
        plotm10(i,1) = i;
        plotm10(i,2:4) = results.num;
        
        % display age
        fprintf(1,strcat(' \t10Be',results.outdisp{1}));
    end;
    
    if nucl26 == 1;
        % sample surface spallation production rate over time
        sample.sp = LSDnu.Al.*Pref26.*thickSF.*sample.othercorr;
        
        % sample muon P
        sample.mu = sample.mu26;
        
        % various parameters
        sample.N = sample.N26; sample.delN = sample.delN26;
        sample.Pref = Pref26; sample.delPref = delPref26;
        sample.l = l26;
        
        % get age
        results = get_1026_age(sample,6E6);
        
        % fill output
        output(i+1,outn(3):outn(4)) = results.outstr;
        
        % fill plot matrix
        plotm26(i,1) = i;
        plotm26(i,2:4) = results.num;
        
        % display age
        fprintf(1,strcat(' \t26Al',results.outdisp{1}));
    end;
    
    fprintf(1,'\n')
    clear sample;
end;

% fix and save output ==============================================================================
if sum(samplein.N10 + samplein.delN10 + samplein.N26 + samplein.delN26) > 0;
    % fix output string
    outstr = '%s';
    for j = 1:size(output,2)-1;
        outstr = strcat(outstr,'\t%s');
    end;
    outstr = strcat(outstr,'\n');

    % fill empty cells with '-'
    nullidx = cellfun(@isempty,output);
    output(nullidx) = {'-'};
    
    % write out-expage.txt
    out = fopen('out-expage.txt','w');
    for i = 1:size(output,1);
        fprintf(out,outstr,output{i,:});
    end;
    fclose(out);
end;
% ==================================================================================================

% plotting ========================================
if plotpointages == 1;
    if exist('plotm10');
        plot_points(plotm10,10);
    end;
    if exist('plotm26');
        plot_points(plotm26,26);
    end;
end;
if plotprobdens == 1;
    if exist('plotm10');
        plot_probdens(plotm10,10,Pref10,delPref10);
    end;
    if exist('plotm26');
        plot_probdens(plotm26,26,Pref26,delPref26);
    end;
end;
% =================================================

toc()
clear;
% end expage function ==============================================================================


% subfunction get_mt ===============================================================================
function [tsimple,mt] = get_mt(sample,Pref,P_St_SF,l,Lsp1,N);
    % Find P according to Stone/Lal - no muon production here!
    % Pref used as ref prod rate
    P_St = Pref * P_St_SF;
    A = l + sample.rho * sample.E / Lsp1;
    if (N < (P_St./A)); % if not saturated: do calculation
        tsimple = (-1/A)*log(1-(N * A / P_St));
        mt = tsimple .* 2;
        if mt < 12060;
            mt = 12060; % don't allow unreasonably short times
        elseif mt > 1e7;
            mt = 1e7;
        end;
    else; % if saturated: use full tv
        tsimple = 5e6;
        mt = 1e7;
    end;
% end subfunction get_mt ===========================================================================


% subfunction get_PmuE =============================================================================
function out = get_PmuE(sample,tv_z,tsimple,RcEst,consts,nucl10,nucl26);
    aged = sample.E .* tsimple; % depth at t simple
    mu_z = [(linspace(0,aged,9) + sample.thick./2).*sample.rho max(tv_z)];
    P_mu_d = P_mu_LSD(mu_z,sample.pressure,RcEst,consts.SPhiInf,nucl10,nucl26,consts,'no');
    if nucl10 == 1; Pmud = P_mu_d.Be; elseif nucl26 == 1; Pmud = P_mu_d.Al; end;
    out = interp1(mu_z,Pmud,tv_z,'pchip') .* sample.othercorr; % P_mu
% end subfunction get_PmuE =========================================================================


% subfunction get_1026_age =========================================================================
function results = get_1026_age(sample,maxage);
    % Calculate N(t) including decay and erosion
    dcf = exp(-sample.tv.*sample.l); % decay factor;
    N_nu = cumtrapz(sample.tv,(sample.sp.*dcf.*sample.dpfs + sample.mu.*dcf)); % pot N back in time
    
    % Look for saturation
    if sample.N <= max(N_nu); % if not saturated
        t_nu = min(interp1(N_nu,sample.tv,sample.N),maxage); % get age by reverse-interpolation
    else;
        t_nu = maxage;
    end;
    
    % uncertainty propagation based on CRONUS v. 2 =============================
    if t_nu > 0;
        % A with integrated average Lsp
        Lsp_avg = interp1(sample.tv,cumtrapz(sample.tv,sample.Lsp.*exp(-sample.l.*sample.tv)),...
            min(t_nu,max(sample.tv)))/interp1(sample.tv,cumtrapz(sample.tv,...
            exp(-sample.l.*sample.tv)),min(t_nu,max(sample.tv)));
        A = sample.l + sample.rho * sample.E ./Lsp_avg;
        % do most of computation
        FP = (sample.N.*A)./(1 - exp(-A.*t_nu));
        delFP = (sample.delPref/sample.Pref) * FP;
        dtdN = 1./(FP - sample.N.*A);
        dtdP = -sample.N./(FP.*FP - sample.N.*A.*FP);
        % make respective delt's
        delt_ext_nu = sqrt(dtdN.^2 * sample.delN.^2 + dtdP.^2 * delFP.^2);
        delt_int_nu = sqrt(dtdN.^2 * sample.delN.^2);
        FP_nu10 = FP;
    else; % t = 0, estimate uncertainty based on conc + unc
        delt_int_nu = interp1(N_nu,sample.tv,sample.N+sample.delN);
        delt_ext_nu = delt_int_nu * (1 + sample.delPref/sample.Pref);
    end; % end uncertainty block ================================================
    
    % fix results
    results.num(1) = round(t_nu);
    results.num(2) = round(delt_ext_nu);
    results.num(3) = round(delt_int_nu);
    results.outstr(1) = {num2str(t_nu,'%.0f')};
    results.outstr(2) = {num2str(delt_ext_nu,'%.0f')};
    results.outstr(3) = {num2str(delt_int_nu,'%.0f')};
    results.outdisp = strcat({' = '},results.outstr{1},{' ± '},results.outstr{2},{' yr'});
    if t_nu+delt_int_nu > maxage; % notification for saturated samples
        results.outdisp = strcat(results.outdisp{1},{' (saturated!)'});
    end;
% end subfunction get_1026_age =====================================================================
    

% subfunction plot_points ==========================================================================
function plot_points(plotm,nucl);
    % This subfunction plots exposure ages as points with external uncertainties
    figure;
    hold on;

    % remove samples without N
    rmidx = find(plotm(:,3) == 0);
    plotm(rmidx,:) = [];

    % fix y legend and max y
    if nucl == 10;
        nucl_legend = '^{10}Be exposure age (yr)';
        ymaxmax = 1E7;
    elseif nucl == 26;
        nucl_legend = '^{26}Al exposure age (yr)';
        ymaxmax = 5E6;
    end;
        
    % plot uncertainty lines
    ymax = [];
    for i = 1:size(plotm,1);
        x = plotm(i,1);
        if plotm(i,3) > 0;
            y1 = plotm(i,2) + plotm(i,3);
            y2 = plotm(i,2) - plotm(i,3);
            ymax(end+1) = y1;
            plot([x x],[y1 y2],'color','black');
        elseif plotm(i,3) == -1;
            plot([x x],[ymaxmax plotm(i,2)],'color','black');
            ymax(end+1) = ymaxmax;
        end;
    end;

    % plot points
    plot(plotm(:,1),plotm(:,2),'.','markersize',15,'color','black');

    % fix plot
    axis([0 size(plotm,1)+1 0 max(ymax)]);
    set(gca,'xtick',[]);
    set(gca,'xticklabel',[]);
    xlabel('Samples');
    ylabel(nucl_legend);
    hold off;
% end subfunction plot_points ======================================================================


% subfunction plot_probdens ========================================================================
function plot_probdens(plotm,nucl,Pref,delPref);
    % This subfunction plots probability density curves
    % red curves: single age probability density curves using internal uncertainty
    % black curve: summed probability density curve
    % black vertical line: weighted mean age
    % grey area: weighted uncertainty with propagated production rate uncertainty added
    figure;
    hold on;

    % fix x legend and max x
    if nucl == 10;
        nucl_legend = '^{10}Be exposure age (yr)';
        xmaxmax = 1E7;
    elseif nucl == 26;
        nucl_legend = '^{26}Al exposure age (yr)';
        xmaxmax = 5E6;
    end;

    % find number of saturated samples
    nsat = numel(find(plotm(:,3) < 0));

    % remove saturated samples and samples without n10
    rmidx = find(plotm(:,3) <= 0);
    plotm(rmidx,:) = [];

    % vectors and matrices for probability estimation
    ages = plotm(:,2);
    uncs = plotm(:,4);
    agemin = min(ages - 4.*uncs);
    agemax = max(ages + 4.*uncs);
    timev = linspace(agemin,agemax,500);
    timem = repmat(timev,numel(ages),1);
    agesm = repmat(ages,1,numel(timev));
    uncsm = repmat(uncs,1,numel(timev));

    % estimate probability distr and sum
    probdensmatr = normpdf(timem,agesm,uncsm);
    probsum = sum(probdensmatr);

    if numel(ages)>1;
        % calculate weigthed mean age
        wage = sum(ages./uncs.^2)/sum(1./uncs.^2);
        wunc = (sum(1./uncs.^2.*(ages-wage).^2)/sum(1./uncs.^2)*(numel(ages)-1)/numel(ages))^0.5;
        wunc = sqrt(wunc^2 + (wage*delPref/Pref)^2); % add prodrate uncertainty
        
        % plot grey age uncertainty region
        uncage = linspace(max(wage-wunc,agemin),min(wage+wunc,agemax),100);
        uncprob = interp1(timev,probsum,uncage,'pchip');
        uncage = [uncage,wage+wunc,wage-wunc];
        uncprob = [uncprob 0 0];
        patch(uncage,uncprob,[0.85 0.85 0.85],'EdgeColor','none');
        
        % plot Pref line
        wageprob = interp1(timev,probsum,wage,'pchip');
        plot([wage wage],[0 wageprob],'color','black');
    end;

    % plot individual sample prob dens curves
    for i = 1:size(probdensmatr,1);
        plot(timev,probdensmatr(i,:),'color','red');
    end

    if numel(ages)>1;
        % plot summed prob dens curve
        plot(timev,sum(probdensmatr),'color','black');
    end;

    % fix plot
    xlim([max(agemin,0) min(agemax,xmaxmax)]);
    set(gca,'ytick',[]);
    set(gca,'yticklabel',[]);
    xlabel(nucl_legend);
    ylabel('Relative probability');
    set(gca,'layer','top'); % plot axis on top
    if nsat > 0; % if saturated samples: display number of excluded samples
        if nsat == 1; satstr = 'sample'; else satstr = 'samples'; end;
        xlims = xlim;
        ylims = ylim;
        xmin = xlims(1) + 0.75*(xlims(2)-xlims(1));
        text(xmin,0.9*ylims(2),{[num2str(nsat,'%.0f'),' saturated'],[satstr,' excluded']});
    end;
    hold off;
% end subfunction plot_probdens ====================================================================
