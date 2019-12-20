function expage()

% expage 10Be and 26Al exposure age calculator.
% Read and fix input, use get_1026_age to calculate exposure rates, and fix and write output.
%
% Based on code written by Greg Balco for the CRONUS calculator v. 2.
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License, version 2, as published by the Free Software Foundation (www.fsf.org).
%
% Jakob Heyman - 2015-2019 (jakob.heyman@gu.se)

clear;
tic();

% What version is this?
ver = '201912';

% plotting? (1 = yes) ==============================================================================
plotpointages = 0; % plots exposure ages as points
plotprobdens = 0; % plots exposure ages as probability density curves (for single groups of ages)
% ==================================================================================================

% fix input ========================================================================================
% variable names for input with variable names in first line
varnames = {'sample','Pflag','std10','std26','isostP','lat','long','elv','thick','dens','shield',...
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
    varsin = {'sample','lat','long','elv','Pflag','thick','dens','shield','erosion','N10',...
        'N10unc','std10','N26','N26unc','std26','samplingyr'};
    typestr = '%s %n %n %n %s %n %n %n %n %n %n %s %n %n %s %n';
else;
    fprintf(1,'ERROR! Something is wrong with the input\n');
    return;
end;
indata = textscan(fid,typestr,'CommentStyle','%'); % scan data
for i = 1:numel(varsin); % fix variables
    samplein.(varsin{i}) = indata{i};
end;
fclose(fid);
% ==================================================================================================

% run and load expage constants
make_consts_expage;
load consts_expage;

% decay constant
l10 = consts.l10; l10unc = consts.l10unc;
l26 = consts.l26; l26unc = consts.l26unc;

% if there is no erosion in input: assume zero erosion
if isfield(samplein,'erosion') == 0; samplein.erosion(1:numel(samplein.sample),1) = 0; end;

% if there is no N10 or N26 in input: fill with 0
if isfield(samplein,'N10') == 0; samplein.N10(1:numel(samplein.sample),1) = 0; end;
if isfield(samplein,'N10unc') == 0; samplein.N10unc(1:numel(samplein.sample),1) = 0; end;
if isfield(samplein,'std10') == 0; samplein.std10(1:numel(samplein.sample),1) = {'0'}; end;
if isfield(samplein,'N26') == 0; samplein.N26(1:numel(samplein.sample),1) = 0; end;
if isfield(samplein,'N26unc') == 0; samplein.N26unc(1:numel(samplein.sample),1) = 0; end;
if isfield(samplein,'std26') == 0; samplein.std26(1:numel(samplein.sample),1) = {'0'}; end;

% convert 10Be concentrations according to standards
[testi,stdi] = ismember(samplein.std10,consts.std10); % find index of standard conversion factors
mult10 = consts.std10_cf(stdi); % pick out conversion factor
samplein.N10 = samplein.N10 .* mult10;
samplein.N10unc = samplein.N10unc .* mult10;

% convert 26Al concentrations according to standards
[testi,stdi] = ismember(samplein.std26,consts.std26); % find index of standard conversion factors
mult26 = consts.std26_cf(stdi); % pick out conversion factor
samplein.N26 = samplein.N26 .* mult26;
samplein.N26unc = samplein.N26unc .* mult26;

% fix longitude values
samplein.long(find(samplein.long < 0)) = samplein.long(find(samplein.long < 0)) + 360;

% fix sample pressure
if isfield(samplein,'pressure') == 0;
    % if there is no pressure flag: use std
    if isfield(samplein,'Pflag') == 0; samplein.Pflag(1:numel(samplein.sample),1) = {'std'}; end;
    stdv = strcmp(samplein.Pflag,'std');
    antv = strcmp(samplein.Pflag,'ant');
    prev = strcmp(samplein.Pflag,'pre');
    samplein.pressure(stdv) = ERA40atm(samplein.lat(stdv),samplein.long(stdv),samplein.elv(stdv));
    samplein.pressure(antv) = antatm(samplein.elv(antv));
    samplein.pressure(prev) = samplein.elv(prev);
end;

% fix erosion rate unit (mm/ka -> cm/yr)
samplein.erosion = samplein.erosion .* 1E-4;

% fix isostatic adjustment input
if isfield(samplein,'isostP');
    samplein.isostP = isost_data_load(samplein.isostP);
end;
            
% initial Lsp for simple age calculation
Lsp1 = consts.Lsp;

% fix simple thickness scaling factor
thick_v = find(samplein.thick > 0);
samplein.thickSF1(1:numel(samplein.thick),1) = 1;
samplein.thickSF1(thick_v) = (Lsp1./(samplein.dens(thick_v).*samplein.thick(thick_v))).*...
    (1 - exp(((-1.*samplein.dens(thick_v).*samplein.thick(thick_v))./Lsp1)));

% fix for output
output(1,1) = {'sample'};
outn(1) = 1;
if sum(samplein.N10+samplein.N10unc) > 0;
    output(1,end+1:end+3) = {'10age(yr)','10uncext(yr)','10uncint(yr)'};
    outn(1) = max(outn)+1;
    outn(2) = max(outn)+2;
end;
if sum(samplein.N26+samplein.N26unc) > 0;
    output(1,end+1:end+3) = {'26age(yr)','26uncext(yr)','26uncint(yr)'};
    outn(3) = max(outn)+1;
    outn(4) = max(outn)+2;
end;

% pick out samples one by one
for i = 1:numel(samplein.lat);
    sample.sample = samplein.sample(i);
    sample.lat = samplein.lat(i);
    sample.long = samplein.long(i);
    sample.thick = samplein.thick(i);
    sample.dens = samplein.dens(i);
    sample.shield = samplein.shield(i);
    sample.erosion = samplein.erosion(i);
    sample.N10 = samplein.N10(i);
    sample.N10unc = samplein.N10unc(i);
    sample.N26 = samplein.N26(i);
    sample.N26unc = samplein.N26unc(i);
    sample.samplingyr = samplein.samplingyr(i);
    sample.pressure = samplein.pressure(i);
    sample.thickSF1 = samplein.thickSF1(i);
    if isfield(samplein,'isostP');
        sample.isostP = samplein.isostP{i};
        sample.elv = samplein.elv(i);
        sample.Pflag = samplein.Pflag(i);
        if strcmp(sample.isostP,'-') || strcmp(sample.Pflag,'pre');
            sample = rmfield(sample,'isostP');
        end;
    end;
    
    % write sample name to output
    output(i+1,1) = sample.sample;
    
    % Set nucl and mt to 0 for both 10/26 and check if there is N10/N26
    nucl10 = 0; nucl26 = 0; mt10 = 0; mt26 = 0;
    if (sample.N10 + sample.N10unc) > 0; nucl10 = 1; end;
    if (sample.N26 + sample.N26unc) > 0; nucl26 = 1; end;
    
    % if no 10Be or 26Al: move on
    if nucl10 + nucl26 == 0;
        continue;
    end;
    
    % Prefs
    Pref10 = consts.Pref10; Pref10unc = consts.Pref10unc;
    Pref26 = consts.Pref26; Pref26unc = consts.Pref26unc;
    
    % display sample name
    fprintf(1,'%.0f. %s',i,sample.sample{1});
    
    % Find P scaling factor according to Stone/Lal
    P_St_SF = stone2000(sample.lat,sample.pressure,1) * sample.thickSF1 * sample.shield;
    
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
    
    % Fix tv, Rc, RcEst, SPhi, and w for sp and mu prod rate scaling
    LSDfix = LSD_fix(sample.lat,sample.long,mt,-1,sample.samplingyr,consts);
    
    % include tv in sample
    sample.tv = LSDfix.tv;
    
    % Production from muons
    if sample.erosion <= 0;
        Pmu = P_mu_expage(sample.thick.*sample.dens./2,sample.pressure,LSDfix.RcEst,...
            consts.SPhiInf,nucl10,nucl26,consts,'no');
        if nucl10 == 1; sample.mu10 = Pmu.mu10 .* sample.shield; end;
        if nucl26 == 1; sample.mu26 = Pmu.mu26 .* sample.shield; end;
    else;
        tv_z = (sample.tv.*sample.erosion + sample.thick./2) .* sample.dens; % T-d vect (g/cm^2)
        if nucl10 == 1;
            sample.mu10 = get_PmuE(sample,tv_z,tsimple10,LSDfix.RcEst,consts,1,0);
        end;
        if nucl26 == 1;
            sample.mu26 = get_PmuE(sample,tv_z,tsimple26,LSDfix.RcEst,consts,0,1);
        end;
    end;
    
    % fix atmospheric pressure if using isostatic adjustment
    if isfield(sample,'isostP');
        % calculate elevation development
        sample.elvv = isost_elv(sample.isostP,sample);
        % calculate atmospheric pressure
        if strcmp(sample.Pflag,'std');
            sample.pressure = ERA40atm(sample.lat,sample.long,sample.elvv);
        elseif strcmp(sample.Pflag,'ant');
            sample.pressure = antatm(sample.elvv);
        end;
        % change Pref to isostatic calibration values
        Pref10 = consts.Pref10iso; Pref10unc = consts.Pref10isounc;
    end;
    
    % spallation production scaling
    Psp = P_sp_expage(sample.pressure,LSDfix.Rc,LSDfix.SPhi,LSDfix.w,consts,nucl10,nucl26);
    
    % interpolate Lsp using CRONUScalc method (Sato 2008; Marrero et al. 2016)
    sample.Lsp = rawattenuationlength(sample.pressure,LSDfix.Rc);
    
    % Thickness scaling factor.
    if sample.thick > 0;
        thickSF = (sample.Lsp./(sample.dens.*sample.thick)).*...
            (1 - exp(((-1.*sample.dens.*sample.thick)./sample.Lsp)));
    else;
        thickSF = 1;
    end;
    
    % spallation depth dependence
    sample.dpfs = exp(-sample.tv.*sample.erosion.*sample.dens./sample.Lsp);
    
    if nucl10 == 1;
        % sample surface spallation production rate over time
        sample.sp = Psp.sp10.*Pref10.*thickSF.*sample.shield;
        
        % sample muon P
        sample.mu = sample.mu10;
        
        % various parameters
        sample.N = sample.N10; sample.Nunc = sample.N10unc;
        sample.Pref = Pref10; sample.Prefunc = Pref10unc;
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
        sample.sp = Psp.sp26.*Pref26.*thickSF.*sample.shield;
        
        % sample muon P
        sample.mu = sample.mu26;
        
        % various parameters
        sample.N = sample.N26; sample.Nunc = sample.N26unc;
        sample.Pref = Pref26; sample.Prefunc = Pref26unc;
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
if sum(samplein.N10 + samplein.N10unc + samplein.N26 + samplein.N26unc) > 0;
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
        plot_probdens(plotm10,10,Pref10,Pref10unc);
    end;
    if exist('plotm26');
        plot_probdens(plotm26,26,Pref26,Pref26unc);
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
    A = l + sample.dens * sample.erosion / Lsp1;
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
    agez = sample.erosion .* tsimple .* sample.dens; % shielding depth at t simple
    agezv = [0 25 50 75 125 250]; % shielding depth break values for number of mu_z points (5-10)
    idx = find(agez>agezv,1,'last');
    mu_z = [(linspace(0,agez,idx+3) + sample.thick.*sample.dens./2) max(tv_z)];
    Pmu_d = P_mu_expage(mu_z,sample.pressure,RcEst,consts.SPhiInf,nucl10,nucl26,consts,'no');
    if nucl10 == 1; Pmud = Pmu_d.mu10; elseif nucl26 == 1; Pmud = Pmu_d.mu26; end;
    out = interp1(mu_z,Pmud,tv_z,'pchip') .* sample.shield; % P_mu
% end subfunction get_PmuE =========================================================================


% subfunction get_1026_age =========================================================================
function results = get_1026_age(sample,maxage);
    % Calculate N(t) including decay and erosion
    dcf = exp(-sample.tv.*sample.l); % decay factor;
    N_nu = cumtrapz(sample.tv,(sample.sp.*dcf.*sample.dpfs + sample.mu.*dcf)); % pot N back in time
    
    % Look for saturation
    if sample.N <= max(N_nu); % if not saturated
        t_nu = min(interp1(N_nu,sample.tv,sample.N),maxage); % get age by reverse-interpolation
        
        % uncertainty propagation based on CRONUS v. 2 ====================================
        if t_nu > 0;
            % A with integrated average Lsp
            Lsp_avg = interp1(sample.tv,cumtrapz(sample.tv,sample.Lsp.*...
                exp(-sample.l.*sample.tv)),min(t_nu,max(sample.tv)))/interp1(sample.tv,...
                cumtrapz(sample.tv,exp(-sample.l.*sample.tv)),min(t_nu,max(sample.tv)));
            A = sample.l + sample.dens * sample.erosion ./Lsp_avg;
            % do most of computation
            FP = (sample.N.*A)./(1 - exp(-A.*t_nu));
            delFP = (sample.Prefunc/sample.Pref) * FP;
            dtdN = 1./(FP - sample.N.*A);
            dtdP = -sample.N./(FP.*FP - sample.N.*A.*FP);
            % make respective delt's
            delt_ext_nu = sqrt(dtdN.^2 * sample.Nunc.^2 + dtdP.^2 * delFP.^2);
            delt_int_nu = sqrt(dtdN.^2 * sample.Nunc.^2);
            FP_nu10 = FP;
        else; % t = 0, estimate uncertainty based on conc + unc
            delt_int_nu = interp1(N_nu,sample.tv,sample.N+sample.Nunc);
            delt_ext_nu = delt_int_nu * (1 + sample.Prefunc/sample.Pref);
        end; % end uncertainty block ======================================================
    else; % if saturated: use maxage for age and uncertainties
        t_nu = maxage;
        delt_int_nu = maxage;
        delt_ext_nu = maxage;
    end;
    
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
function plot_probdens(plotm,nucl,Pref,Prefunc);
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

    % remove saturated samples and samples without N
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
        [wage,wunc] = evm(ages',uncs');
        wunc = sqrt(wunc^2 + (wage*Prefunc/Pref)^2); % add prodrate uncertainty
        
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
    plot(timev',probdensmatr','color','red');

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


% subfunction evm ==================================================================================
function [wage,wunc] = evm(agev,uncv);
    % this subfunction calculates weighted average and uncertainty using the expected value method
    % (Birch and Singh 2014)
    % fix data (make matrices)
    agem = repmat(agev,numel(agev),1);
    uncm = repmat(uncv,numel(uncv),1);
    xm = repmat(agev',1,numel(agev));
    
    % calculate probability for each sample age
    Mui = sum(sqrt(2./(pi.*(2.*uncm).^2)).*exp(-(xm-agem).^2./(2.*uncm.^2)))./numel(agev);
    
    % calculate summed probability for all samples
    Muj = sum(Mui);
    
    % calculate sample weights
    wi = Mui./Muj;
    
    % calculate weighted mean age
    wage = sum(wi.*agev);
    
    % uncertainty estimation
    uncint = sqrt(sum(wi.^2.*uncv.^2));
    uncext = sqrt(sum(wi.*(agev-wage).^2));
    wunc = max(uncint,uncext);
% end subfunction evm ==============================================================================
