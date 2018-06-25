function erosion();

% expage 10Be and 26Al erosion rate calculator.
% Read and fix input, use get_1026_erosion to calculate erosion rates, and fix and write output.
%
% Based on code written by Greg Balco for the CRONUS calculator v. 2.
%
% This program is free software; you can redistribute it and/or modify it under the terms of the GNU
% General Public License, version 2, as published by the Free Software Foundation (www.fsf.org).
%
% Jakob Heyman - 2015-2018 (jakob.heyman@gu.se)

tic();

% What version is this?
ver = '201806';

% count number of input columns in line 1 of input
inid = fopen('input.txt','r');
line1 = fgets(inid);
line1 = regexprep(line1,' +',' ');
numcols = numel(strfind(line1,' ')) + numel(strfind(line1,sprintf('\t',''))) + 1;
fclose(inid);

% read input file
if numcols == 15; % if no erosion rate in input
    [samplein.sample_name,samplein.lat,samplein.long,samplein.elv,samplein.aa,samplein.thick,...
        samplein.rho,samplein.othercorr,samplein.N10,samplein.delN10,samplein.be_stds,...
        samplein.N26,samplein.delN26,samplein.al_stds,samplein.samplingyr] = ...
        textread('input.txt','%s %n %n %n %s %n %n %n %n %n %s %n %n %s %n','commentstyle',...
        'matlab');
else; % has to be 16 columns!
    [samplein.sample_name,samplein.lat,samplein.long,samplein.elv,samplein.aa,samplein.thick,...
        samplein.rho,samplein.othercorr,samplein.E,samplein.N10,samplein.delN10,samplein.be_stds,...
        samplein.N26,samplein.delN26,samplein.al_stds,samplein.samplingyr] = ...
        textread('input.txt','%s %n %n %n %s %n %n %n %n %n %n %s %n %n %s %n','commentstyle',...
        'matlab');
end;

% run and load expage constants
make_consts_expage;
load consts_expage;

% constants
Pref10 = consts.P10_ref_nu; delPref10 = consts.delP10_ref_nu;
Pref26 = consts.P26_ref_nu; delPref26 = consts.delP26_ref_nu;
% Decay constant
l10 = consts.l10; dell10 = consts.dell10;
l26 = consts.l26; dell26 = consts.dell26;
% muon parameters 10
sigma010 = consts.sigma0_10nu; delsigma010 = consts.delsigma0_10nu;
k_negpartial10 = consts.k_negpartial10; delk_negpartial10 = consts.delk_negpartial10;
fstar10 = consts.fstar10nu; delfstar10 = consts.delfstar10nu;
% muon parameters 26
sigma026 = consts.sigma0_26nu; delsigma026 = consts.delsigma0_26nu;
k_negpartial26 = consts.k_negpartial26; delk_negpartial26 = consts.delk_negpartial26;
fstar26 = consts.fstar26nu; delfstar26 = consts.delfstar26nu;

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
    sample.sample_name = samplein.sample_name(i);
    sample.lat = samplein.lat(i);
    sample.long = samplein.long(i);
    sample.elv = samplein.elv(i);
    sample.aa = samplein.aa(i);
    sample.thick = samplein.thick(i);
    sample.rho = samplein.rho(i);
    sample.othercorr = samplein.othercorr(i);
    sample.N10 = samplein.N10(i);
    sample.delN10 = samplein.delN10(i);
    sample.be_stds = samplein.be_stds(i);
    sample.N26 = samplein.N26(i);
    sample.delN26 = samplein.delN26(i);
    sample.al_stds = samplein.al_stds(i);
    sample.samplingyr = samplein.samplingyr(i);
    sample.pressure = samplein.pressure(i);
    
    % write sample name to output
    output(i+1,1) = sample.sample_name;
    
    % check if there is any N10 or N26
    nucl10 = 0; nucl26 = 0;
    if sample.N10 > 0; nucl10 = 1; end;
    if sample.N26 > 0; nucl26 = 1; end;
    
    if nucl10 + nucl26 == 0;
        continue;
    end;
    
    % display sample name
    fprintf(1,'%.0f. %s',i,sample.sample_name{1});

    % Fix w,Rc,SPhi, for sp and nu prod rate scaling 10 Ma back in time
    LSDfix = LSD_fix(sample.lat,sample.long,1E7,-1,consts);
    
    % Age Relative to t0=2010 - LSD tv from LSD_fix
    % tv = [0:10:50 60:100:50060 51060:1000:2000060 logspace(log10(2001060),7,200)];
    
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
    sample.tv = tv;
    
    % interpolate Lsp (Sato 2008; Marrero et al. 2016)
    sample.Lsp = rawattenuationlength(sample.pressure,Rc);
    sample.LspAv = trapz(tv,sample.Lsp)./tv(end); % pick out average

    % thickness scaling factor.
    if sample.thick > 0;
        sample.thickSF = (sample.Lsp./(sample.rho.*sample.thick)).*...
            (1 - exp(((-1.*sample.rho.*sample.thick)./sample.Lsp)));
    else;
        sample.thickSF = 1;
        if numel(tv) > 1; sample.thickSF(1:numel(tv)) = 1; end;
    end;

    % spallation production scaling
    LSDnu = LSDspal(sample.pressure,Rc,SPhi,LSDfix.w,nucl10,nucl26,consts);

    % muon production
    P_mu_full = P_mu_LSD(0,sample.pressure,LSDfix.RcEst,consts.SPhiInf,nucl10,nucl26,consts,'yes');
    
    % Precompute P_mu(z) to ~200,000 g/cm2
    % This log-spacing setup for the step size has relative accuracy near 
    % 1e-3 at 1000 m/Myr erosion rate. 
    % start at the mid-depth of the sample.
    sample.z_mu = [0 logspace(0,5.3,100)]+(sample.thick.*sample.rho./2);
    P_mu_z = P_mu_LSD(sample.z_mu,sample.pressure,LSDfix.RcEst,consts.SPhiInf,nucl10,nucl26,...
        consts,'no');
    
    % if there is N10: do erosion rate calculation and report results
    if nucl10 == 1;
        % sample- and nuclide-specific parameters
        Nspec.N = sample.N10; Nspec.delN = sample.delN10;
        Nspec.Psps = LSDnu.Be;
        Nspec.P_fast = P_mu_full.P_fast10;
        Nspec.P_neg = P_mu_full.P_neg10;
        Nspec.Pmu_z = P_mu_z.Be .* sample.othercorr;
        
        % nuclide-specific constants
        Nspec.Pref = Pref10; Nspec.delPref = delPref10;
        Nspec.l = l10; Nspec.dell = dell10;
        Nspec.sigma0 = sigma010; Nspec.delsigma0 = delsigma010;
        Nspec.k_negpartial = k_negpartial10; Nspec.delk_negpartial = delk_negpartial10;
        Nspec.fstar = fstar10; Nspec.delfstar = delfstar10;
        
        % get erosion and uncertainty
        results = get_1026_erosion(sample,Nspec);
        
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
        Nspec.N = sample.N26; Nspec.delN = sample.delN26;
        Nspec.Psps = LSDnu.Al;
        Nspec.P_fast = P_mu_full.P_fast26;
        Nspec.P_neg = P_mu_full.P_neg26;
        Nspec.Pmu_z = P_mu_z.Al .* sample.othercorr;
        
        % nuclide-specific constants
        Nspec.Pref = Pref26; Nspec.delPref = delPref26;
        Nspec.l = l26; Nspec.dell = dell26;
        Nspec.sigma0 = sigma026; Nspec.delsigma0 = delsigma026;
        Nspec.k_negpartial = k_negpartial26; Nspec.delk_negpartial = delk_negpartial26;
        Nspec.fstar = fstar26; Nspec.delfstar = delfstar26;
        
        % get erosion and uncertainty
        results = get_1026_erosion(sample,Nspec);
        
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
P_nu = Nspec.Psps.*Nspec.Pref.*sample.othercorr;
P_mu = (Nspec.P_fast + Nspec.P_neg) .* sample.othercorr;
P_fast = Nspec.P_fast .* sample.othercorr;
P_neg = Nspec.P_neg .* sample.othercorr;
P_mu_z = Nspec.Pmu_z .* sample.othercorr;

% initial guess
% average P_nu
P_nuAv = trapz(sample.tv,P_nu)./sample.tv(end);
P_temp = P_nuAv + P_mu;
E_lal = sample.LspAv.*(P_temp./Nspec.N - Nspec.l);
x0 = E_lal;

% constants block for subfunction ET_objective
c3.tv = sample.tv;
c3.z_mu = sample.z_mu-(sample.thick.*sample.rho./2); % take 1/2 depth away so t will match P
c3.P_mu_z = P_mu_z;
c3.l = Nspec.l;
c3.tsf = sample.thickSF;
c3.L = sample.Lsp;

% test saturation and ask fzero for the erosion rates.
P_mud = interp1(sample.z_mu,P_mu_z,sample.thick.*sample.rho./2,'pchip'); % P_mu at 1/2 sample d
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
rel_delP0 = Nspec.delPref./Nspec.Pref;

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
    dEdN = (1e4./(sample.rho.*2.*Nspec.delN)) .* ...
        ( (fzero(@(y) E_simple(y,Psp0,Pmu0,L,Lmu,Nspec.l,(Nspec.N+Nspec.delN)),x_nu)) - ...
        (fzero(@(y) E_simple(y,Psp0,Pmu0,L,Lmu,Nspec.l,(Nspec.N-Nspec.delN)),x_nu)) );
    
    dEdPsp0 = (1e4./(sample.rho.*2.*delPsp0)) .* ...
        ( (fzero(@(y) E_simple(y,(Psp0+delPsp0),Pmu0,L,Lmu,Nspec.l,Nspec.N),x_nu)) - ...
        (fzero(@(y) E_simple(y,(Psp0-delPsp0),Pmu0,L,Lmu,Nspec.l,Nspec.N),x_nu)) );
    
    dEdPmu0 = (1e4./(sample.rho.*2.*delPmu0)) .* ...
        ( (fzero(@(y) E_simple(y,Psp0,(Pmu0+delPmu0),L,Lmu,Nspec.l,Nspec.N),x_nu)) - ...
        (fzero(@(y) E_simple(y,Psp0,(Pmu0-delPmu0),L,Lmu,Nspec.l,Nspec.N),x_nu)) );
    
    % Add in quadrature to get the uncertainties.
    delE_ext = sqrt((dEdPsp0.*delPsp0).^2+(dEdPmu0.*delPmu0).^2+(dEdN.*Nspec.delN).^2);
    delE_int = abs(dEdN .* Nspec.delN);
else;
    x_nu = 0;
    delEintN = Nspec.N-Nspec.delN;
    [delE_int,fval_nu,exitflag_nu,output] = fzero(@(x) ET_objective(x,c3,delEintN),x0,opts);
    delEextN = Nspec.N-sqrt(Nspec.delN^2 + (Nspec.N.*Nspec.delPref./Nspec.Pref)^2);
    [delE_ext,fval_nu,exitflag_nu,output] = fzero(@(x) ET_objective(x,c3,delEextN),x0,opts);
    delE_int =  max(0,delE_int*1e4/sample.rho); % convert to m/Ma
    delE_ext =  max(0,delE_ext*1e4/sample.rho); % convert to m/Ma
end;
% end of uncertainty block ==================================================================

% report
results.Pmu0 = Pmu0; % Surface production rate due to muons;
results.Egcm2yr = [x_nu]; % erosion rate in gcm2/yr
results.EmMyr = results.Egcm2yr * 1e4 / sample.rho;    % erosion rate in m/Ma
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
