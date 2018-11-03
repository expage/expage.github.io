% script calculating depth profile muon production without sigma0 and fstar
% This is free software: you can use/copy/modify/distribute as long as you keep it free/open.
% Jakob Heyman - 2017-2018 (jakob.heyman@gu.se)

clear all;
tic();

% Choose site. Options: Beacon, LeymonHigh, LeymonLow, LaCiotat
site = 'Beacon';

% run and load expage constants
make_consts_expage;
load consts_expage;

% display site
fprintf(1,'Site: %s\n',site);

% load site data from mucalib_data
mucal = mucalib_data(site);

% calculate atmospheric pressure
atm = ERA40atm(mucal.lat,mucal.lon,mucal.elv);

% Fix RcEst for muon calculation
LSDfix = LSD_fix(mucal.lat,mucal.lon,1E7,-1,consts);

% muon production (without sigma0 and fstar!)
fprintf(1,'calculating muon P (up to %.0f)\n',mucal.dz(end));
Pmud = mucalib_Pmu(mucal.dz,atm,LSDfix.RcEst,consts.SPhiInf,1,1,consts);
fprintf(1,'\ndone!\n');

fprintf(1,'writing mucalib_depthcalc_out.txt...');
out = fopen('mucalib_depthcalc_out.txt','w');
fprintf(out,'%s\t%s\t%s\t%s\t%s\n','dz(g/cm2)','Pfast10','Pneg10','Pfast26','Pneg26');
for i = 1:numel(mucal.dz);
    dz = num2str(mucal.dz(i),'%.0f');
    Pfast10 = num2str(Pmud.P_fast10(i),'%.8e');
    Pneg10 = num2str(Pmud.P_neg10(i),'%.8f');
    Pfast26 = num2str(Pmud.P_fast26(i),'%.8e');
    Pneg26 = num2str(Pmud.P_neg26(i),'%.8f');
    fprintf(out,'%s\t%s\t%s\t%s\t%s\n',dz,Pfast10,Pneg10,Pfast26,Pneg26);
end;
fclose(out);
fprintf(1,' done!\n');

toc();
