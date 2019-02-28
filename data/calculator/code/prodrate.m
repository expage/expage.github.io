function prodrate()

% Function for 10Be and 26Al reference production rate calibration.
% Calibrate site reference 10Be and/or 26Al production rate using chi-square minimization.
% This is free software: you can use/copy/modify/distribute as long as you keep it free/open.
% Jakob Heyman - 2015-2018 (jakob.heyman@gu.se)

clear;
close all;
tic();

% What version is this?
ver = '201902';

% do choices here ==================================================================================
% min and max Pref (atoms/g/yr)
Pmin10 = 3;
Pmax10 = 6;
Pmin26 = 20;
Pmax26 = 50;

% plotting / cluster test?
plotting = 1; % plot Pref probability curves (1 = yes)
Pcluster = 1; % exclude outliers to try to achieve a well-clustered group Pref (1 = yes)

% parameters for cluster test / outlier rejection
chiprob = 0.1; % lower p-value limit for chi-square probability test
mingroupn = 3; % minimum number of samples in well-clustered group
maxoutratio = 1/3; % maximum outlier ratio
% ==================================================================================================

% read and fix input file
[samplein.sample_name,samplein.lat,samplein.long,samplein.elv,samplein.aa,samplein.thick,...
    samplein.rho,samplein.othercorr,samplein.E,samplein.N10,samplein.delN10,samplein.be_stds,...
    samplein.N26,samplein.delN26,samplein.al_stds,samplein.samplingyr,samplein.calage,...
    samplein.calageunc] = textread...
    ('input.txt','%s %n %n %n %s %n %n %n %n %n %n %s %n %n %s %n %n %n','commentstyle','matlab');

% full number of samples
fulln = numel(samplein.sample_name);

% run and load expage constants
make_consts_expage;
load consts_expage;

% convert 10Be concentrations according to standards
for i = 1:fulln;
    be_mult(i,1) = consts.be_stds_cfs(strcmp(samplein.be_stds(i),consts.be_stds_names));
end;
samplein.N10 = samplein.N10 .* be_mult;
samplein.delN10 = samplein.delN10 .* be_mult;

% convert 26Al concentrations according to standards
for i = 1:fulln;
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

% Pref for simple age estimates
P_ref_10 = consts.P10_ref_nu; delP_ref_10 = consts.delP10_ref_nu;
P_ref_26 = consts.P26_ref_nu; delP_ref_26 = consts.delP26_ref_nu;

% Be-10 and Al-26 decay constants and spallation attenuation approximation
l10 = consts.l10;
l26 = consts.l26;
Lsp1 = consts.Lsp;

% fix simple thickness scaling factor
thick_v = find(samplein.thick > 0);
samplein.thickSF1(1:fulln,1) = 1;
samplein.thickSF1(thick_v) = (Lsp1./(samplein.rho(thick_v).*samplein.thick(thick_v))).*...
    (1 - exp(((-1.*samplein.rho(thick_v).*samplein.thick(thick_v))./Lsp1)));

% make Pref vectors
Pvect10 = [Pmin10:0.01:Pmax10];
Pvect26 = [Pmin26:0.02:Pmax26];

% declare age and unc matrixes
age10 = []; unc10 = []; age26 = []; unc26 = [];
calage10 = []; calunc10 = []; calage26 = []; calunc26 = [];

% number of 10Be and 26Al samples
num10 = 0;
num26 = 0;
OKsample = []; OK10 = []; OK26 = []; % variables for outlier removal

% fix for output
output(1,1) = {'sample'};
outn(1) = 1;
if sum(samplein.N10+samplein.delN10) > 0;
    output(1,end+1:end+2) = {'P10(at/g/yr)','P10unc(at/g/yr)'};
    outn(1) = max(outn)+1;
    outn(2) = max(outn)+1;
    if Pcluster == 1;
        output(1,end+1) = {'P10-included?'};
        outn(3) = max(outn)+1;
    end;
end;
if sum(samplein.N26+samplein.delN26) > 0;
    output(1,end+1:end+2) = {'P26(at/g/yr)','P26unc(at/g/yr)'};
    outn(4) = max(outn)+1;
    outn(5) = max(outn)+1;
    if Pcluster == 1;
        output(1,end+1) = {'P26-included?'};
        outn(6) = max(outn)+1;
    end;
end;

% pick out samples one by one
for i = 1:fulln;    
    sample.sample_name = samplein.sample_name(i,:);
    sample.lat = samplein.lat(i);
    sample.long = samplein.long(i);
    sample.thick = samplein.thick(i);
    sample.rho = samplein.rho(i);
    sample.othercorr = samplein.othercorr(i);
    sample.E = samplein.E(i);
    sample.N10 = samplein.N10(i);
    sample.delN10 = samplein.delN10(i);
    sample.be_stds = samplein.be_stds(i,:);
    sample.N26 = samplein.N26(i);
    sample.delN26 = samplein.delN26(i);
    sample.al_stds = samplein.al_stds(i,:);
    sample.samplingyr = samplein.samplingyr(i,:);
    sample.calage = samplein.calage(i);
    sample.calageunc = samplein.calageunc(i);
    sample.pressure = samplein.pressure(i);
    sample.thickSF1 = samplein.thickSF1(i);
    
    % write sample name to output
    output(i+1,1) = sample.sample_name;
    
    % Set nucl and mt to 0 for both 10/26 and check if there is N10/N26
    nucl10 = 0; nucl26 = 0; mt10 = 0; mt26 = 0;
    if (sample.N10 + sample.delN10) > 0; nucl10 = 1; num10 = num10 + 1; end;
    if (sample.N26 + sample.delN26) > 0; nucl26 = 1; num26 = num26 + 1; end;
    
    if nucl10 + nucl26 == 0;
        continue;
    end;
    
    % display sample name
    fprintf(1,'%.0f. %s',i,sample.sample_name{1});
    
    % Atoms/g measurement
    N10 = sample.N10; delN10 = sample.delN10;
    N26 = sample.N26; delN26 = sample.delN26;
    
    % Find P scaling factor according to Stone/Lal
    P_St_SF = stone2000(sample.lat,sample.pressure,1) * sample.thickSF1 * sample.othercorr;
    
    % if 10Be measured: calculate max time
    if nucl10 == 1;
        [tsimple10,mt10] = get_mt(sample,P_ref_10,P_St_SF,l10,Lsp1,sample.N10);
    end;
    
    % if 26Al measured: calculate max time
    if nucl26 == 1;
        [tsimple26,mt26] = get_mt(sample,P_ref_26,P_St_SF,l26,Lsp1,sample.N26);
    end;
    
    % pick largest of mt10 and mt26 as max time
    mt = max(mt10,mt26);
    
    % Age Relative to t0=2010 - LSD tv from LSDfix
    % tv = [0:10:50 60:100:50060 51060:1000:2000060 logspace(log10(2001060),7,200)];
    
    % Fix w,Rc,SPhi, for sp and nu prod rate scaling
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
        P_mu = P_mu_expage(sample.thick.*sample.rho./2,sample.pressure,LSDfix.RcEst,...
            consts.SPhiInf,nucl10,nucl26,consts,'no');
        if nucl10 == 1; sample.mu10 = P_mu.mu10 .* sample.othercorr; end;
        if nucl26 == 1; sample.mu26 = P_mu.mu26 .* sample.othercorr; end;
    else;
        tv_z = (tv.*sample.E + sample.thick./2) .* sample.rho; % time - depth vect (g/cm^2)
        if nucl10 == 1;
            aged10 = sample.E .* tsimple10; % depth at t simple
            mu_z10 = [(linspace(0,aged10,9) + sample.thick./2).*sample.rho max(tv_z)];
            P_mu_d10 = P_mu_expage(mu_z10,sample.pressure,LSDfix.RcEst,consts.SPhiInf,1,0,...
                consts,'no'); % Pmu at d
            sample.mu10 = interp1(mu_z10,P_mu_d10.mu10,tv_z,'pchip') .* sample.othercorr; % P_mu
        end;
        if nucl26 == 1;
            aged26 = sample.E .* tsimple26; % depth at t simple
            mu_z26 = [(linspace(0,aged26,9) + sample.thick./2).*sample.rho max(tv_z)];
            P_mu_d26 = P_mu_expage(mu_z26,sample.pressure,LSDfix.RcEst,consts.SPhiInf,0,1,...
                consts,'no'); % Pmu at d
            sample.mu26 = interp1(mu_z26,P_mu_d26.mu26,tv_z,'pchip') .* sample.othercorr; % P_mu
        end;
    end;
    
    % spallation production scaling
    Psp = LSDspal(sample.pressure,Rc,SPhi,LSDfix.w,nucl10,nucl26,consts);
    
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
    dpfs = exp(-tv.*sample.E.*sample.rho./sample.Lsp); % spallation depth dependence
    
    if nucl10 == 1;
        % for decay
        dcf10 = exp(-tv.*l10); % decay factor;
        
        % sample surface spallation production rate over time including decay and erosion
        sample.sp = bsxfun(@times,Psp.sp10'.*dcf10'.*dpfs'.*thickSF'.*sample.othercorr,Pvect10);
        
        % sample muon P
        sample.mu = repmat(sample.mu10'.*dcf10',1,numel(Pvect10));
        
        % full production including decay and erosion
        sample.Pfull = sample.sp + sample.mu;
        
        % various parameters
        sample.N = sample.N10; sample.delN = sample.delN10;
        sample.Pvect = Pvect10; sample.l = l10;
        
        % get ages and sample Pref
        [age10(num10,:),unc10(num10,:),Prefi,Prefiunc] = get_ages(sample,'10');
        
        % stop if no Pref
        if isnan(Prefi); return; end;
        
        % fill calage and calunc matrix
        calage10(num10,1:numel(Pvect10)) = sample.calage;
        calunc10(num10,1:numel(Pvect10)) = sample.calageunc;
        
        % fill output
        output(i+1,outn(1):outn(2)) = {num2str(Prefi,'%.3f'),num2str(Prefiunc,'%.3f')};
        
        % display Pref
        fprintf(1,' \tP10 = %s ± %s at/g/yr',output{i+1,outn(1)},output{i+1,outn(2)});
        
        % fill Pref vector for uncertainty estimation, plotting, and cluster analysis
        Pref10v(num10) = Prefi;
        Punc10v(num10) = sqrt(Prefiunc^2 - (sample.calageunc/sample.calage*Prefi)^2);
        Prow10v(num10) = i; % input row number
    end;
    
    if nucl26 == 1;
        % for decay
        dcf26 = exp(-tv.*l26); % decay factor;
        
        % sample surface spallation production rate over time including decay and erosion
        sample.sp = bsxfun(@times,Psp.sp26'.*dcf26'.*dpfs'.*thickSF'.*sample.othercorr,Pvect26);
        
        % sample muon P
        sample.mu = repmat(sample.mu26'.*dcf26',1,numel(Pvect26));
        
        % full production including decay and erosion
        sample.Pfull = sample.sp + sample.mu;
        
        % various parameters
        sample.N = sample.N26; sample.delN = sample.delN26;
        sample.Pvect = Pvect26; sample.l = l26;
        
        % get ages and sample Pref
        [age26(num26,:),unc26(num26,:),Prefi,Prefiunc] = get_ages(sample,'26');
        
        % stop if no Pref
        if isnan(Prefi); return; end;
        
        % fill calage and calunc matrix
        calage26(num26,1:numel(Pvect26)) = sample.calage;
        calunc26(num26,1:numel(Pvect26)) = sample.calageunc;
        
        % fill output
        output(i+1,outn(4):outn(5)) = {num2str(Prefi,'%.3f'),num2str(Prefiunc,'%.3f')};
        
        % display Pref
        fprintf(1,' \tP26 = %s ± %s at/g/yr',output{i+1,outn(4)},output{i+1,outn(5)});
        
        % fill Pref vector for uncertainty estimation, plotting, and cluster analysis
        Pref26v(num26) = Prefi;
        Punc26v(num26) = sqrt(Prefiunc^2 - (sample.calageunc/sample.calage*Prefi)^2);
        Prow26v(num26) = i; % input row number
    end;
    
    fprintf(1,'\n');
    clear sample;
end;

% fix for output
if num10 > 1 || num26 > 1;
    output(fulln+2,1) = {'Pref-group'};
    output(fulln+3,1) = {'Rchi2 | Pvalue'};
end;

% calculate group Pref10
if num10 > 1;
    [Pref,Prefunc,rchisq,Pvalue] = get_Pref(age10,unc10,calage10,calunc10,Pvect10,Pref10v,Punc10v);
    
    % if doing cluster analysis
    if Pcluster == 1;
        % fix variable cl for function get_cluster
        cl.Pref=Pref; cl.Prefunc=Prefunc; cl.rchisq=rchisq; cl.Pvalue=Pvalue;
        cl.age=age10; cl.unc=unc10; cl.calage=calage10; cl.calunc=calunc10; cl.Pvect=Pvect10;
        cl.Prefv=Pref10v; cl.Puncv=Punc10v; cl.Prowv=Prow10v; cl.num=num10;
        cl.maxoutratio=maxoutratio; cl.mingroupn=mingroupn; cl.chiprob=chiprob;
        
        % find clustered Pref by removing outliers
        [Pref,Prefunc,rchisq,Pvalue,OKsample,OK10] = get_cluster(cl);
        
        % mark OK samples in output
        output(OKsample+1,outn(3)) = {'X'};
        
        % mark OK Pref in output
        if Pvalue>=chiprob && numel(OKsample)>=mingroupn;
            output(fulln+2:fulln+3,outn(3)) = {'X'};
        end;
    end;
    
    % pick out full dataset P-ref, unc, and chisquare
    output(fulln+2,outn(1)) = {num2str(Pref,'%.3f')};
    output(fulln+2,outn(2)) = {num2str(Prefunc,'%.3f')};
    output(fulln+3,outn(1)) = {num2str(rchisq,'%.3f')};
    output(fulln+3,outn(2)) = {num2str(Pvalue,'%.3f')};
    
    % display group Pref
    fprintf(1,'Pref10 = %s ± %s at/g/yr    ',output{fulln+2,outn(1)},output{fulln+2,outn(2)});
    fprintf(1,'R-chi2 = %s    P-value = %s',output{fulln+3,outn(1)},output{fulln+3,outn(2)});
    if Pcluster==1 && numel(OKsample)<num10;
        fprintf(1,'    (%.0f outlier',num10-numel(OKsample));
        if num10-numel(OKsample) > 1; fprintf(1,'s)'); else; fprintf(1,')'); end;
    end;
    fprintf(1,'\n');
    
    % do plotting here...
    if plotting == 1;
        plot_Pref(Pmin10-0.5,Pmax10+0.5,0.01,Pref10v',Punc10v',Pref,Prefunc,'^{10}Be',OK10);
    end;
end;

% calculate group Pref26
if num26 > 1;
    [Pref,Prefunc,rchisq,Pvalue] = get_Pref(age26,unc26,calage26,calunc26,Pvect26,Pref26v,Punc26v);
    
    % if doing cluster analysis
    if Pcluster == 1;
        % fix variable cl for function get_cluster
        cl.Pref=Pref; cl.Prefunc=Prefunc; cl.rchisq=rchisq; cl.Pvalue=Pvalue;
        cl.age=age26; cl.unc=unc26; cl.calage=calage26; cl.calunc=calunc26; cl.Pvect=Pvect26;
        cl.Prefv=Pref26v; cl.Puncv=Punc26v; cl.Prowv=Prow26v; cl.num=num26;
        cl.maxoutratio=maxoutratio; cl.mingroupn=mingroupn; cl.chiprob=chiprob;
        
        % find clustered Pref by removing outliers
        [Pref,Prefunc,rchisq,Pvalue,OKsample,OK26] = get_cluster(cl);
        
        % mark OK samples in output
        output(OKsample+1,outn(6)) = {'X'};
        
        % mark OK Pref in output
        if Pvalue>=chiprob && numel(OKsample)>=mingroupn;
            output(fulln+2:fulln+3,outn(6)) = {'X'};
        end;
    end;
    
    % pick out full dataset P-ref, unc, and chisquare
    output(fulln+2,outn(4)) = {num2str(Pref,'%.3f')};
    output(fulln+2,outn(5)) = {num2str(Prefunc,'%.3f')};
    output(fulln+3,outn(4)) = {num2str(rchisq,'%.3f')};
    output(fulln+3,outn(5)) = {num2str(Pvalue,'%.3f')};
    
    % display group Pref
    fprintf(1,'Pref26 = %s ± %s at/g/yr    ',output{fulln+2,outn(4)},output{fulln+2,outn(5)});
    fprintf(1,'R-chi2 = %s    P-value = %s',output{fulln+3,outn(4)},output{fulln+3,outn(5)});
    if Pcluster==1 && numel(OKsample)<num26;
        fprintf(1,'    (%.0f outlier',num26-numel(OKsample));
        if num26-numel(OKsample) > 1; fprintf(1,'s)'); else; fprintf(1,')'); end;
    end;
    fprintf(1,'\n');
    
    % do plotting here...
    if plotting == 1;
        plot_Pref(Pmin26-5,Pmax26+5,0.02,Pref26v',Punc26v',Pref,Prefunc,'^{26}Al',OK26);
    end;
end;

% fix and save output ============================
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
    
    % write out-glacialE.txt
    out = fopen('out-prodrate.txt','w');
    for i = 1:size(output,1);
        fprintf(out,outstr,output{i,:});
    end;
    fclose(out);
end;
% ================================================

toc()
clear;
% end prodrate function ============================================================================


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


% subfunction get_ages =============================================================================
function [age,ageunc,Pref,Prefunc] = get_ages(sample,nuclstr);
    % Calculate N(t) including decay and erosion
    N_m = cumtrapz(sample.tv',sample.Pfull); % potential N back in time
    
    idx = 1; % start index for age
    
    % test for saturation
    if sample.N > min(max(N_m));
        idx = min(find(sample.N <= N_m(end,:)));
        age(1:idx-1) = 1E7; % set age to 10 Ma for saturation
        fprintf(1,'Sample saturated for Pref%s < %.2f at/g/yr\n',nuclstr,sample.Pvect(idx));
    end;
    
    % calculate ages
    for j = idx:numel(sample.Pvect);
        age(j) = interp1(N_m(:,j),sample.tv',sample.N);
    end;
    
    % test for too low Pmax
    if max(age) < sample.calage;
        fprintf(1,' \tToo high Pmin%s (change in prodrate.m and re-run!)\n',nuclstr);
    end;
    if min(age) > sample.calage;
        fprintf(1,' \tToo low Pmax%s (change in prodrate.m and re-run!)\n',nuclstr);
    end;
    
    % Error propagation from Balco et al. (2008) CRONUS calculator
    % A with decay-weighted average Lsp
    Lsp_avg = interp1(sample.tv,cumtrapz(sample.tv,sample.Lsp.*exp(-sample.l.*sample.tv)),age)/...
        interp1(sample.tv,cumtrapz(sample.tv,exp(-sample.l.*sample.tv)),age);
    A = sample.l + sample.rho .* sample.E ./Lsp_avg;
    FP = (sample.N.*A)./(1 - exp(-A.*age));
    dtdN = 1./(FP - sample.N.*A);
    ageunc = sqrt(dtdN.^2 .* sample.delN.^2);
    
    % find individual sample Pref and Pref uncertainty including calage uncertainty
    Pref = interp1(age,sample.Pvect,sample.calage,'pchip');
    Prefunc = sqrt((interp1(sample.Pvect,ageunc,Pref)/sample.calage*Pref)^2 + ...
        (sample.calageunc/sample.calage*Pref)^2);
% end subfunction get_ages =========================================================================


% subfunction get_Pref =============================================================================
function [Pref,Prefunc,rchisq,Pvalue] = get_Pref(age,unc,calage,calunc,Pvect,Prefv,Puncv);
    % calculate weighted ages for ages minus calage
    agetest = sum((age - calage)./unc.^2)./sum(1./unc.^2);
    
    % pick out Pref
    Pref = interp1(agetest,Pvect,0,'pchip');
    
    % find individual sample age for Pref
    ageP_sample = interp1(Pvect',age',Pref,'pchip');
    
    % find individual sample age unc for Pref
    uncP_sample = interp1(Pvect',unc',Pref,'pchip');
    
    % pick out exact chi square
    rchisq = 1/(size(age,1)-1) .* sum(((ageP_sample - calage(:,1)')./uncP_sample).^2);
    
    % calculate P value
    Pvalue = 1 - chi2cdf(rchisq.*(size(age,1)-1),size(age,1)-1);
    
    % calculate ref prod rate unc (not including calib age unc) - unbiased estimator
    Puncint = sqrt(sum(1./Puncv.^2.*(Prefv-Pref).^2)/...
        (sum(1./Puncv.^2)-(sum((1./Puncv.^2).^2)/sum(1./Puncv.^2))));
    
    % calculate calib age unc part
    calage_unc = mean(calunc(:,1)./calage(:,1));
    
    % calculate total prod rate uncertainty
    Prefunc = sqrt(Puncint^2 + (calage_unc * Pref)^2);
% end subfunction get_Pref =========================================================================


% subfunction get_cluster ==========================================================================
function [Pref,Prefunc,rchisq,Pvalue,OKsample,OKnucl] = get_cluster(cl);
    % in variables
    inPref = cl.Pref; inPrefunc = cl.Prefunc; inrchisq = cl.rchisq; inPvalue = cl.Pvalue;
    OKnucl = (1:1:numel(cl.Prowv));
    
    % define maximum number of outliers
    remove = floor(cl.num*cl.maxoutratio);
    r = 1; % outlier counter
    
    % loop for outlier removal
    while cl.Pvalue<cl.chiprob && remove>=r && cl.num-r>=cl.mingroupn;
        % calculate deviation for all individual samples
        chidev = ((cl.Prefv-cl.Pref)./cl.Puncv).^2;
        
        % find index of sample with largest dev
        [value,rmv_idx] = max(chidev);
        
        % remove outlier
        cl.age(rmv_idx,:) = []; cl.unc(rmv_idx,:) = [];
        cl.calage(rmv_idx,:) = []; cl.calunc(rmv_idx,:) = [];
        cl.Prefv(rmv_idx) = []; cl.Puncv(rmv_idx) = [];
        cl.Prowv(rmv_idx) = []; OKnucl(rmv_idx) = [];
        
        % calculate Pref
        [cl.Pref,cl.Prefunc,cl.rchisq,cl.Pvalue] = ...
            get_Pref(cl.age,cl.unc,cl.calage,cl.calunc,cl.Pvect,cl.Prefv,cl.Puncv);
        
        r = r + 1; % add 1 to outlier
    end;
    
    % if not well-clustered: use input variables and remove row numbers
    if cl.Pvalue < cl.chiprob || numel(cl.Prowv) < cl.mingroupn;
        cl.Prowv = []; OKnucl = [];
        % use input variables
        cl.Pref = inPref; cl.Prefunc = inPrefunc; cl.rchisq = inrchisq; cl.Pvalue = inPvalue;
    end;
    
    % fix output
    Pref = cl.Pref;
    Prefunc = cl.Prefunc;
    rchisq = cl.rchisq;
    Pvalue = cl.Pvalue;
    OKsample = cl.Prowv;
% end subfunction get_cluster ======================================================================


% subfunction plot_Pref ============================================================================
function plot_Pref(plmin,plmax,step,Prefv,Puncv,Pref,Prefunc,nuclstr,OKsample);
    figure;
    hold on;
    box on;
    
    % split Pref vectors in clustered and outliers
    cPrefv = Prefv(OKsample);
    cPuncv = Puncv(OKsample);
    oPrefv = Prefv; oPrefv(OKsample) = [];
    oPuncv = Puncv; oPuncv(OKsample) = [];
    
    % vectors and matrices for probability estimation
    xv = linspace(plmin,plmax,(plmax-plmin)/step+1);
    cxm = repmat(xv,numel(cPrefv),1); if numel(OKsample) == 0; cxm = []; end;
    cPrefm = repmat(cPrefv,1,numel(xv));
    cPuncm = repmat(cPuncv,1,numel(xv));
    oxm = repmat(xv,numel(oPrefv),1);
    oPrefm = repmat(oPrefv,1,numel(xv));
    oPuncm = repmat(oPuncv,1,numel(xv));
    
    % estimate probability distribution
    cprobdensmatr = normpdf(cxm,cPrefm,cPuncm);
    oprobdensmatr = normpdf(oxm,oPrefm,oPuncm);
    
    % pick out summed probability distribution
    if numel(OKsample) > 0; probsum = sum(cprobdensmatr); else; probsum = sum(oprobdensmatr); end;
    
    % plot grey Pref uncertainty region
    uncP = linspace(Pref-Prefunc,Pref+Prefunc,100);
    uncprob = interp1(xv,probsum,uncP,'pchip');
    uncP = [uncP,Pref+Prefunc,Pref-Prefunc];
    uncprob = [uncprob 0 0];
    patch(uncP,uncprob,[0.85 0.85 0.85],'EdgeColor','none');
    
    % plot Pref line
    Prefprob = interp1(xv,probsum,Pref,'pchip');
    plot([Pref Pref],[0 Prefprob],'color','black');
    
    % plot individual sample prob dens curves (outliers)
    for j = 1:size(oprobdensmatr,1);
        plot(xv',oprobdensmatr(j,:)','color','blue');
    end
    
    % plot individual sample prob dens curves (clustered samples)
    for j = 1:size(cprobdensmatr,1);
        plot(xv',cprobdensmatr(j,:)','color','red');
    end
    
    % plot summed prob dens curve
    plot(xv',probsum','color','black');
    
    xlabel(strcat({'Ref '},nuclstr,{' prodrate (atoms/g/yr)'}));
    set(gca,'ytick',[]);
    ylabel('Relative probability');
    set(gca,'layer','top'); % plot axis on top
    hold off;
% end subfunction plot_Pref ========================================================================
