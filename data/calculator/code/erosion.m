function erosion();

% expage 10Be and 26Al erosion rate calculator.
% Read and fix input, use get_1026_erosion to calculate erosion rates, and fix and write output.
%
% Based on code written by Greg Balco for the CRONUS calculator v. 2.
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License, version 2, as published by the Free Software Foundation (www.fsf.org).
%
% Jakob Heyman - 2015-2019 (jakob.heyman@gu.se)

tic();

% What version is this?
ver = '201912';

% fix input ========================================================================================
% variable names for input with variable names in first line
varnames = {'sample','Pflag','std10','std26','lat','long','elv','thick','dens','shield','N10',...
    'N10unc','N26','N26unc','samplingyr','pressure'};
vartypes = {'%s','%s','%s','%s','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n','%n'};
% read input file
fid = fopen('input.txt');
varsin = strsplit(fgetl(fid)); % read first line
if(ismember(varsin,varnames)); % if first line contain only variable names
    [testi,vari] = ismember(varsin,varnames); % find index of varnames
    typestr = vartypes{vari(1)}; % fix type string
    for i = 2:numel(vari); % fix type string
        typestr = [typestr ' ' vartypes{vari(i)}];
    end;
elseif numel(varsin) == 15; % if no variable names in first line
    frewind(fid); % read from first line
    varsin = {'sample','lat','long','elv','Pflag','thick','dens','shield','N10','N10unc','std10',...
        'N26','N26unc','std26','samplingyr'};
    typestr = '%s %n %n %n %s %n %n %n %n %n %s %n %n %s %n';
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

% Prefs
Pref10 = consts.Pref10; Pref10unc = consts.Pref10unc;
Pref26 = consts.Pref26; Pref26unc = consts.Pref26unc;
% decay constant
l10 = consts.l10; l10unc = consts.l10unc;
l26 = consts.l26; l26unc = consts.l26unc;

% muon parameters 10
sigma010 = consts.sigma0_10nu; delsigma010 = consts.delsigma0_10nu;
k_negpartial10 = consts.k_negpartial10; delk_negpartial10 = consts.delk_negpartial10;
fstar10 = consts.fstar10nu; delfstar10 = consts.delfstar10nu;
% muon parameters 26
sigma026 = consts.sigma0_26nu; delsigma026 = consts.delsigma0_26nu;
k_negpartial26 = consts.k_negpartial26; delk_negpartial26 = consts.delk_negpartial26;
fstar26 = consts.fstar26nu; delfstar26 = consts.delfstar26nu;

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

% fix for output
output(1,1) = {'sample'};
outn(1) = 1;
if sum(samplein.N10) > 0;
    output(1,end+1:end+3) = {'10E(mm/ka)','10uncext(mm/ka)','10uncint(mm/ka)'};
    outn(1) = max(outn)+1;
    outn(2) = max(outn)+2;
end;
if sum(samplein.N26) > 0;
    output(1,end+1:end+3) = {'26E(mm/ka)','26uncext(mm/ka)','26uncint(mm/ka)'};
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
    sample.N10 = samplein.N10(i);
    sample.N10unc = samplein.N10unc(i);
    sample.N26 = samplein.N26(i);
    sample.N26unc = samplein.N26unc(i);
    sample.samplingyr = samplein.samplingyr(i);
    sample.pressure = samplein.pressure(i);
    
    % write sample name to output
    output(i+1,1) = sample.sample;
    
    % check if there is any N10 or N26
    nucl10 = 0; nucl26 = 0;
    if sample.N10 > 0; nucl10 = 1; end;
    if sample.N26 > 0; nucl26 = 1; end;
    
    if nucl10 + nucl26 == 0;
        continue;
    end;
    
    % display sample name
    fprintf(1,'%.0f. %s',i,sample.sample{1});

    % Fix tv, Rc, SPhi, and w for sp and mu prod rate scaling 10 Ma back in time
    LSDfix = LSD_fix(sample.lat,sample.long,1E7,-1,sample.samplingyr,consts);
    
    % Age Relative to t0=2010 - LSD tv from LSD_fix
    % tv = [0:10:50 60:100:50060 51060:1000:2000060 logspace(log10(2001060),7,200)];
    
    % include tv in sample
    sample.tv = LSDfix.tv;
    
    % interpolate Lsp (Sato 2008; Marrero et al. 2016)
    sample.Lsp = rawattenuationlength(sample.pressure,LSDfix.Rc);
    sample.LspAv = trapz(sample.tv,sample.Lsp)./sample.tv(end); % pick out average

    % thickness scaling factor.
    if sample.thick > 0;
        sample.thickSF = (sample.Lsp./(sample.dens.*sample.thick)).*...
            (1 - exp(((-1.*sample.dens.*sample.thick)./sample.Lsp)));
    else;
        sample.thickSF = 1;
        if numel(sample.tv) > 1; sample.thickSF(1:numel(sample.tv)) = 1; end;
    end;
    
    % spallation production scaling
    Psp = P_sp_expage(sample.pressure,LSDfix.Rc,LSDfix.SPhi,LSDfix.w,consts,nucl10,nucl26);
    
    % muon production
    P_mu_full = P_mu_expage(0,sample.pressure,LSDfix.RcEst,consts.SPhiInf,nucl10,nucl26,consts,...
        'yes');
    
    % Precompute P_mu(z) to ~200,000 g/cm2
    % log-spacing setup for the step size with a relative accuracy near 1E-3 at 1000 m/Ma E-rate.
    % z_mu from CRONUS calculator changed in the upper 30 g/cm2 to reduce steps from 101 to 80. 
    % start at the mid-depth of the sample.
    sample.z_mu = [(0:3:27) logspace(1.48,5.3,70)]+(sample.thick.*sample.dens./2);
    P_mu_z = P_mu_expage(sample.z_mu,sample.pressure,LSDfix.RcEst,consts.SPhiInf,nucl10,nucl26,...
        consts,'no');
    
    % if there is N10: do erosion rate calculation and report results
    if nucl10 == 1;
        % sample- and nuclide-specific parameters
        Nspec.N = sample.N10; Nspec.delN = sample.N10unc;
        Nspec.Psps = Psp.sp10;
        Nspec.P_fast = P_mu_full.P_fast10;
        Nspec.P_neg = P_mu_full.P_neg10;
        Nspec.Pmu_z = P_mu_z.mu10 .* sample.shield;
        
        % nuclide-specific constants
        Nspec.Pref = Pref10; Nspec.Prefunc = Pref10unc;
        Nspec.l = l10; Nspec.lunc = l10unc;
        Nspec.sigma0 = sigma010; Nspec.delsigma0 = delsigma010;
        Nspec.k_negpartial = k_negpartial10; Nspec.delk_negpartial = delk_negpartial10;
        Nspec.fstar = fstar10; Nspec.delfstar = delfstar10;
        
        % get erosion and uncertainty
        warning('off'); % turn off warnings
        results = get_1026_erosion(sample,Nspec);
        warning('on');
        
        % fill output
        output(i+1,outn(1):outn(2)) = results.outstr;
        
        % display results
        fprintf(1,' \t10Be = %s ± %s mm/ka',output{i+1,outn(1)},output{i+1,outn(1)+1});
        if results.EmMyr == 0;
            fprintf(1,' (saturated)');
        end;
    end;
    
    % if there is N26: do erosion rate calculation and report results
    if nucl26 == 1;
        % sample- and nuclide-specific parameters
        Nspec.N = sample.N26; Nspec.delN = sample.N26unc;
        Nspec.Psps = Psp.sp26;
        Nspec.P_fast = P_mu_full.P_fast26;
        Nspec.P_neg = P_mu_full.P_neg26;
        Nspec.Pmu_z = P_mu_z.mu26 .* sample.shield;
        
        % nuclide-specific constants
        Nspec.Pref = Pref26; Nspec.Prefunc = Pref26unc;
        Nspec.l = l26; Nspec.lunc = l26unc;
        Nspec.sigma0 = sigma026; Nspec.delsigma0 = delsigma026;
        Nspec.k_negpartial = k_negpartial26; Nspec.delk_negpartial = delk_negpartial26;
        Nspec.fstar = fstar26; Nspec.delfstar = delfstar26;
        
        % get erosion and uncertainty
        warning('off'); % turn off warnings
        results = get_1026_erosion(sample,Nspec);
        warning('on');
        
        % fill output
        output(i+1,outn(3):outn(4)) = results.outstr;
        
        % display results
        fprintf(1,' \t26Al = %s ± %s mm/ka',output{i+1,outn(3)},output{i+1,outn(3)+1});
        if results.EmMyr == 0;
            fprintf(1,' (saturated)');
        end;
    end;
    
    fprintf(1,'\n');
    clear sample;
end;

% fix and save output =======================
if sum(samplein.N10 + samplein.N26) > 0;
    % fix output string
    outstr = '%s';
    for j = 1:size(output,2)-1;
        outstr = strcat(outstr,'\t%s');
    end;
    outstr = strcat(outstr,'\n');

    % fill empty cells with '-'
    nullidx = cellfun(@isempty,output);
    output(nullidx) = {'-'};
    
    % write out-erosion.txt
    out = fopen('out-erosion.txt','w');
    for i = 1:size(output,1);
        fprintf(out,outstr,output{i,:});
    end;
    fclose(out);
end;
% ===========================================

toc()
clear;
% end erosion function =============================================================================


% subfunction get_1026_erosion - calculates erosion and uncertainty ================================
function results = get_1026_erosion(sample,Nspec);
P_nu = Nspec.Psps.*Nspec.Pref.*sample.shield;
P_mu = (Nspec.P_fast + Nspec.P_neg) .* sample.shield;
P_fast = Nspec.P_fast .* sample.shield;
P_neg = Nspec.P_neg .* sample.shield;
P_mu_z = Nspec.Pmu_z .* sample.shield;

% initial guess
% average P_nu
P_nuAv = trapz(sample.tv,P_nu)./sample.tv(end);
P_temp = P_nuAv + P_mu;
E_lal = sample.LspAv.*(P_temp./Nspec.N - Nspec.l);
x0 = E_lal;

% constants block for subfunction ET_objective
c3.tv = sample.tv;
c3.z_mu = sample.z_mu-(sample.thick.*sample.dens./2); % take 1/2 depth away so t will match P
c3.P_mu_z = P_mu_z;
c3.l = Nspec.l;
c3.tsf = sample.thickSF;
c3.L = sample.Lsp;

% test saturation and ask fzero for the erosion rates.
P_mud = interp1(sample.z_mu,P_mu_z,sample.thick.*sample.dens./2,'pchip'); % P_mu at 1/2 sample d
% zero erosion N
Ntest = trapz(sample.tv,(P_nu.*exp(-sample.tv.*Nspec.l) + P_mud.*exp(-sample.tv.*Nspec.l)));

% fix parameters for ET_objective
opts = optimset('fzero');
opts = optimset(opts,'tolx',1e-8,'display','off');
c3.P_sp_t = P_nu;
        
% calculate erosion rate
[x_nu,fval_nu,exitflag_nu,output] = fzero(@(x) ET_objective(x,c3,Nspec.N),x0,opts);
% diagnostics
diag_nu = ET_objective(x_nu,c3,Nspec.N,'yes');

% uncertainty propagation ===================================================================
% Common information
Pmu0 = P_mu_z(1);
delPfast = P_fast .* (Nspec.delsigma0 ./ Nspec.sigma0);
delPneg = P_neg .* sqrt((Nspec.delk_negpartial./Nspec.k_negpartial)^2 + ...
    (Nspec.delfstar./Nspec.fstar)^2);
delPmu0 = sqrt(delPfast.^2 + delPneg.^2);
L = sample.LspAv;

% SF - varying variable assignments
diag = diag_nu;
rel_delP0 = Nspec.Prefunc./Nspec.Pref;

% Conditional on having a low enough N+Ndel (saturation check)
if Nspec.N+Nspec.delN <= Ntest;
    % find what Lmu ought to be 
    Lmu = x_nu ./ ((Pmu0./diag.Nmu) - Nspec.l);
    
    % Find what P0 ought to be
    Psp0 = diag.Nsp.*(Nspec.l + (x_nu./L));
    delPsp0 = Psp0 .* rel_delP0;
    
    % Find the derivatives with respect to the uncertain parameters.  
    % Here we're calculating centered derivatives using fzero and
    % the subfunction E_simple.
    dEdN = (1e4./(sample.dens.*2.*Nspec.delN)) .* ...
        ( (fzero(@(y) E_simple(y,Psp0,Pmu0,L,Lmu,Nspec.l,(Nspec.N+Nspec.delN)),x_nu)) - ...
        (fzero(@(y) E_simple(y,Psp0,Pmu0,L,Lmu,Nspec.l,(Nspec.N-Nspec.delN)),x_nu)) );
    
    dEdPsp0 = (1e4./(sample.dens.*2.*delPsp0)) .* ...
        ( (fzero(@(y) E_simple(y,(Psp0+delPsp0),Pmu0,L,Lmu,Nspec.l,Nspec.N),x_nu)) - ...
        (fzero(@(y) E_simple(y,(Psp0-delPsp0),Pmu0,L,Lmu,Nspec.l,Nspec.N),x_nu)) );
    
    dEdPmu0 = (1e4./(sample.dens.*2.*delPmu0)) .* ...
        ( (fzero(@(y) E_simple(y,Psp0,(Pmu0+delPmu0),L,Lmu,Nspec.l,Nspec.N),x_nu)) - ...
        (fzero(@(y) E_simple(y,Psp0,(Pmu0-delPmu0),L,Lmu,Nspec.l,Nspec.N),x_nu)) );
    
    % Add in quadrature to get the uncertainties.
    delE_ext = sqrt((dEdPsp0.*delPsp0).^2+(dEdPmu0.*delPmu0).^2+(dEdN.*Nspec.delN).^2);
    delE_int = abs(dEdN .* Nspec.delN);
else;
    x_nu = 0;
    delEintN = Nspec.N-Nspec.delN;
    [delE_int,fval_nu,exitflag_nu,output] = fzero(@(x) ET_objective(x,c3,delEintN),x0,opts);
    delEextN = Nspec.N-sqrt(Nspec.delN^2 + (Nspec.N.*Nspec.Prefunc./Nspec.Pref)^2);
    [delE_ext,fval_nu,exitflag_nu,output] = fzero(@(x) ET_objective(x,c3,delEextN),x0,opts);
    delE_int =  max(0,delE_int*1e4/sample.dens); % convert to m/Ma
    delE_ext =  max(0,delE_ext*1e4/sample.dens); % convert to m/Ma
end;
% end of uncertainty block ==================================================================

% report
results.Pmu0 = Pmu0; % Surface production rate due to muons;
results.Egcm2yr = [x_nu]; % erosion rate in gcm2/yr
results.EmMyr = results.Egcm2yr * 1e4 / sample.dens;    % erosion rate in m/Ma
results.delE_int = [delE_int]; % internal error in m/Ma
results.delE_ext = [delE_ext]; % external error in m/Ma
% diagnostics
results.fzero_status = [exitflag_nu]; % exit status of fzero
results.fval = [fval_nu]; % objective function value from fzero
% outstring
if results.EmMyr == 0;
    results.outstr(1) = {'0'};
else;
    results.outstr(1) = {num2str(results.EmMyr,'%.3f')};
end;
if results.delE_int == 0;
    results.outstr(2) = {'0'};
    results.outstr(3) = {'0'};
else;
    results.outstr(2) = {num2str(results.delE_ext,'%.3f')};
    results.outstr(3) = {num2str(results.delE_int,'%.3f')};
end;
% end subfunction get_1026_erosion =================================================================


% subfunction ET_objective =========================================================================
function miss = ET_objective(E,cs,target,dFlag);
% Forward calculator for N under a particular erosion rate. Calculates expected N for a particular
% erosion rate; returns difference between predicted and measured N. Used to implicitly solve for E.
%
% needs:
%   E erosion rate - what is being solved for (g/cm2/yr)
%   target - measured number of atoms (atoms/g)
%   structure cs containing time and depth information
%       cs.tv - time vector for spallation (yr)
%       cs.P_sp_t - P_sp to match tv (or a scalar) (atoms/g/yr)
%           This is surface, not thickness-averaged, P
%       cs.z_mu - vector of depths (g/cm2) 
%       cs.P_mu_z - P(z) for muons - matches z_mu (atoms/g/yr)
%       cs.l - decay constant (1/yr)
%       cs.tsf - thickness scaling factor (nondimensional)
%       cs.L - Effective attenuation length for spallation (g/cm2)
%
% Input argument dFlag is optional - 'yes' yields a structure containing diagnostic information. 
%
% Output argument miss is difference between predicted and measured nuclide conc (atoms/g). 
%

% Checks
if nargin < 4;
    dFlag = 'no';
end;

if size(cs.tv) ~= size(cs.P_sp_t);
    error('Mismatched tv and P_sp');
end;

if size(cs.z_mu) ~= size(cs.P_mu_z);
    error('Mismatched z_mu and P_mu_z');
end;

% 1. Forward integration for muons....by trapezoidal integration
% Find the t's corresponding to z_mu at E
t_mu = cs.z_mu./E;
% If fzero asked for a negative erosion rate,
% set Nmu to the saturation concentration and carry on. 
% Should yield a negative solution for E, thus flagging saturation.
if E <= 0;
    Nmu = cs.P_mu_z(1)./cs.l;
else;
    % Actually do it
    % muons use linear average for thickness...this is dealt with upstream
    Nmu = trapz(t_mu,(cs.P_mu_z.*exp(-cs.l.*t_mu)));
end;
% The accuracy of this is set by the spacing of z_mu upstream. See the 
% hard-copy docs for details. 

% Forward integration for spallation...using integral formula
% with average P(t) in the time step...
if length(cs.P_sp_t) > 1; % time-dependent P - average in timesteps
    P1 = cs.P_sp_t(1:end-1); P2 = cs.P_sp_t(2:end); Pav = (P1+P2)./2;
    L1 = cs.L(1:end-1); L2 = cs.L(2:end); Lav = (L1+L2)./2;
    tsf1 = cs.tsf(1:end-1); tsf2 = cs.tsf(2:end); tsfav = (tsf1+tsf2)./2;
else; % non-time-dependent P, scalar
    Pav = cs.P_sp_t;
    Lav = cs.L;
    tsfav = cs.tsf;
end;
% Do integration in each time step, using integral formula in depth
t1 = cs.tv(1:end-1); t2 = cs.tv(2:end);
A = (cs.l + E./Lav);
if length(cs.P_sp_t) == 1;
    Nsp = cs.tsf.*Pav./A; % use zero-to-infinity analytical formula
else % analytical formula by pieces
    stepN = (Pav.*tsfav./A).*( exp(-A.*t1) - exp(-A.*t2) );
    Nsp = sum(stepN);
    % go to infinity at end of calculation
    finalt1 = max(cs.tv);
    finalP = cs.P_sp_t(end);
    finalN = (finalP./A(end)).*(exp(-A(end).*finalt1));
    Nsp = Nsp + finalN;
end;

% Diagnostic return option
if strcmp(dFlag,'yes');
    miss.ver = ver;
    miss.Nmu = Nmu;
    miss.Nsp = Nsp;
    if length(cs.P_sp_t) == 1;
        miss.expected_Nsp = cs.tsf.*cs.P_sp_t./A;
        miss.PP = cs.P_sp_t;
    end;
    miss.A = A;
    miss.tsf = cs.tsf;
    miss.PP = cs.P_sp_t;
else;
    N = Nmu + Nsp;
    miss = target - N;
end;
% end subfunction ET_objective =====================================================================


% subfunction E_simple =============================================================================
function out = E_simple(x,Psp,Pmu,Lsp,Lmu,l,target);
% calculates miss between N computed with simple erosion rate expression and measured N
% used in simplified error-propagation scheme
% x is erosion rate (g/cm2/yr)

N = (Psp./(l + x./Lsp)) + (Pmu./(l + x./Lmu));
out = N - target;
% end subfunction E_simple =========================================================================
