function out = mucalib_data(site);

% function fixing depth profile data for mucalib calibration

% Load sample site data
% ==================================================================================================
if strcmp(site,'Beacon');
    % Beacon Heights site, Antarctica
    out.lat = -77.85350472;
    out.lon = 160.7731844;
    out.elv = 2183;
    out.samplingyr = 2009;
    out.dz = (0:1:8000); % g/cm2
elseif strcmp(site,'LeymonHigh');
    % Leymon High site, NW Spain
    out.lat = 42.065;
    out.lon = -7.014;
    out.elv = 1277;
    out.samplingyr = 2010;
    out.dz = [(0:1:6000) (6005:5:10000) (10100:100:20000) (21000:1000:100000) ...
        (110000:10000:200000)]; % g/cm2
elseif strcmp(site,'LeymonLow');
    % Leymon Low site, NW Spain
    out.lat = 42.064;
    out.lon = -7.013;
    out.elv = 1246;
    out.samplingyr = 2010;
    out.dz = [(0:1:6000) (6005:5:10000) (10100:100:20000) (21000:1000:100000) ...
        (110000:10000:200000)]; % g/cm2
elseif strcmp(site,'LaCiotat');
    % La Ciotat site, SE France
    out.lat = 43.17926;
    out.lon = 5.57567;
    out.elv = 310;
    out.samplingyr = 2009;
    out.dz = [(0:1:6000) (6005:5:10000) (10100:100:20000) (21000:1000:100000) ...
        (110000:10000:200000)]; % g/cm2
end;

% Load measured depth profile 10Be and 26Al data (plus initial guesses for mucalib.m)
% ==================================================================================================
if strcmp(site,'Beacon');
    % An initial guess at the calibration parameters.
    % [Erate; Lambda10; Lambda26; fstar10; sigma010; fstar26; sigma026]
    out.pinit=[0.016; 142.0; 141.0; 1.9; 0.23; 12.7; 3.8];
    
    % 10Be and 26Al standards
    out.std10 = '07KNSTD';
    out.std26 = 'KNSTD';
    
    % depth profile 10Be and 26Al concentrations and depths
    % data include primary data (Phillips et al. 2016) from UW/LLNL
    % data source: supplementary material from Balco (2017)
    % samples marked with * have error-weighted mean 10/26 concentrations from double measurements
    % {sample_name  dz10(g/cm2)  N10(at/g)  N10unc(at/g)  dz26(g/cm2)  N26(at/g)  N26unc(at/g)}
    depth_profile = ...
    {'09-DV-001-BCO*'	1.6	59000224	162102	1.6	228938540	829349
    '0.02-0.03*'	5.8	57243780	383158	5.8	230764098	370523
    '0.03-0.04*'	8.2	56842257	358250	8.2	230491668	1095654
    '0.04-0.05'	10.5	56104631	323726	10.5	226829087	5740997
    '0.05-0.06'	12.8	55390344	501066	12.8	224687026	4349888
    '0.06-0.07'	15.2	54703114	295369	15.2	219755424	4226653
    '0.08-0.09*'	19.8	53072593	1205150	19.8	211165053	4616081
    '0.10-0.11'	24.5	51289913	365656	24.5	205480206	3341999
    '0.15-0.16'	36.1	48621799	301543	36.1	190179765	6081633
    '0.18-0.19*'	43.1	46129538	981550	43.1	178103427	8222479
    '0.29-0.30'	68.8	38990141	432431	68.8	148930602	5756102
    '0.33-0.34'	78.2	35366824	331275	78.2	136742494	4542771
    '0.49-0.50'	115.6	28560472	263649	115.6	111196399	5154986
    '0.53-0.54*'	124.9	27105125	251020	124.9	104122686	1399763
    '0.79-0.80'	184.8	17197202	140354	184.8	68306527	2337417
    '0.98-0.99'	228.9	12775873	95881	228.9	47823514	1519463
    '1.02-1.03'	238.2	12517042	307281	238.2	45976064	833322
    '1.51-1.52'	351.1	5464487	134272	351.1	21993257	478493
    '2.01-2.02'	463.8	2502696	29356	463.8	10172458	189576
    '2.50-2.51'	572.9	1205906	11988	572.9	4937509	111629
    '3.36-3.38'	768.2	381101	5579	768.2	1664479	54706
    '4.49-4.51'	1024	162691	3052	1024	741992	25039
    '5.02-5.04'	1144.5	132850	3407	1144.5	581568	55262
    '6.09-6.11'	1394	91265	1722	1394	473647	32644
    '7.14-7.16'	1644	78145	2237	1644	385827	23740
    '8.60-8.62'	2010	59447	1215	2010	328328	20067
    '9.98-10.00'	2313	53066	1074	2313	233303	13353
    '13.00-13.02'	3013	37886	946	3013	173483	19082
    '14.98-15.03'	3508	33814	876	3508	172025	25619
    '17.84-17.86'	4139	27094	1002	4139	122107	9153
    '20.10-20.13'	4655	25367	1940	4655	114194	11654
    '22.17-22.20'	5127	22948	1302	5127	89518	9133
    '25.20-25.22'	5832	18602	576	5832	88158	9084
    '27.46-27.50'	6360	16876	682	6360	72327	17161};
elseif strcmp(site,'LeymonHigh');
    % An initial guess at the calibration parameters.
    % [Erate; Lambda10; Lambda26; fstar10; sigma010; fstar26; sigma026]
    out.pinit=[5; 160; 160; 1.9; 0.23; 12.7; 3.8];
    
    % 10Be and 26Al standards
    out.std10 = 'NIST_27900';
    out.std26 = 'SMAL11';
    
    % depth profile 10Be and 26Al concentrations and depths
    % {sample_name  dz10(g/cm2)  N10(at/g)  N10unc(at/g)  dz26(g/cm2)  N26(at/g)  N26unc(at/g)}
    depth_profile = ...
    {'L-H-1'	2.64	339400	8700	2.64	2149900	110000
    'L-H-33'	85.14	224600	6500	85.14	1513300	77300
    'L-H-66'	168.3	147600	4900	168.3	939500	60800
    'L-H-98'	263.42	94700	3300	263.42	714200	49300
    'L-H-125'	334.63	71600	2900	334.63	444800	31900
    'L-H-171'	454.86	52300	2700	454.86	326900	25900
    'L-H-198'	530.64	44500	2300	530.64	268100	33800
    'L-H-231'	612.15	34200	2400	612.15	212500	21500
    'L-H-264'	696.96	29300	2400	696.96	186900	20600
    'L-H-297'	772.2	29300	1600	772.2	173400	17700
    'L-H-330'	858	29700	2200	858	139100	15500
    'L-H-363'	947.43	27400	1800	947.43	167800	24600
    'L-H-396'	1021.68	23900	1500	1021.68	116900	17200
    'L-H-429'	1098.24	24300	1500	1098.24	132900	34900
    'L-H-459'	1156.68	27600	2200	1156.68	162500	21800
    'L-H-495'	1277.1	27800	1700	1277.1	134900	17000
    'L-H-599'	1557.4	21000	2800	1557.4	102700	25300
    'L-H-698'	1814.8	26200	4200	1814.8	105200	26300
    'L-H-797'	2072.2	18000	2100	2072.2	109100	24600
    'L-H-899'	2337.4	19400	1600	2337.4	95100	22700
    'L-H-995'	2587	17800	1300	2587	111900	22800
    'L-H-15.45'	4017	13700	900	4017	76500	12800
    'L-H-20'	5200	12600	1700	5200	49200	11900};
elseif strcmp(site,'LeymonLow');
    % An initial guess at the calibration parameters.
    % [Erate; Lambda10; Lambda26; fstar10; sigma010; fstar26; sigma026]
    out.pinit=[3; 170; 130; 2; 0.26; 14; 3.7];
    
    % 10Be and 26Al standards
    out.std10 = 'NIST_27900';
    out.std26 = 'SMAL11';
    
    % depth profile 10Be and 26Al concentrations and depths
    % {sample_name  dz10(g/cm2)  N10(at/g)  N10unc(at/g)  dz26(g/cm2)  N26(at/g)  N26unc(at/g)}
    depth_profile = ...
    {'L-L-1'	2.63	667700	21200	2.63	3362600	182000
    'L-L-32'	82.56	448000	13000	82.56	1522100	89400
    'L-L-69'	181.47	247900	7100	181.47	1213000	75300
    'L-L-101'	265.63	180000	6000	265.63	1085600	92100
    'L-L-125'	328.75	128500	4400	328.75	826600	72500
    'L-L-165'	433.95	84400	3600	433.95	431500	38300
    'L-L-198'	520.74	65400	4600	520.74	211800	21600
    'L-L-231'	607.53	54100	2800	607.53	319900	30000
    'L-L-264'	694.32	43300	2100	694.32	210200	31200
    'L-L-297'	781.11	39700	2700	781.11	283700	38100
    'L-L-330'	851.4	36100	2200	851.4	249300	41200
    'L-L-363'	925.65	34900	2000	925.65	175000	35900
    'L-L-394'	985	35400	2900	985	218800	24400
    'L-L-429'	1102.53	30900	2000	1102.53	122000	19700
    'L-L-462'	1205.82	30100	1900	1205.82	197100	30100
    'L-L-498'	1318.7	29100	1400	1318.7	140900	22700
    'L-L-665'	1748.95	23300	1900	1748.95	189400	34000
    'L-L-698'	1835.74	24300	2100	1835.74	154600	22300
    'L-L-863'	2269.69	25900	2500	2269.69	117500	29400
    'L-L-961'	2527.43	22900	2100	2527.43	102200	17500
    'L-L-1015'	2669.45	18100	1200	2669.45	123700	20000
    'L-L-15'	3945	13500	400	3945	53800	11900
    'L-L-20'	5260	11700	2900	5260	48300	10600
    'L-L-25'	6575	8200	500	6575	52200	13600};
elseif strcmp(site,'LaCiotat');
    % An initial guess at the calibration parameters.
    % [Erate; Lambda10; Lambda26; fstar10; sigma010; fstar26; sigma026]
    out.pinit=[6; 140; 140; 1.5; 0.22; 24; 0.6];
    
    % 10Be and 26Al standards
    out.std10 = 'NIST_27900';
    out.std26 = 'SMAL11';
    
    % depth profile 10Be concentrations and depths
    % {sample_name  dz10(g/cm2)  N10(at/g)  N10unc(at/g)}
    depth_profile10 = ...
    {'CIOT-01'	17.46	116700	17000
    'CIOT-02'	124.37	67900	3300
    'CIOT-03'	233.66	36300	2900
    'CIOT-04'	358.1	28100	2600
    'CIOT-05'	466.79	22500	1700
    'CIOT-06B'	559.74	17600	900
    'CIOT-06S'	559.74	16500	900
    'CIOT-07'	700.89	15500	800
    'CIOT-08'	768.7	12000	800
    'CIOT-09'	879.54	17000	1000
    'CIOT-10'	1093.39	12500	900
    'CIOT-15'	2238.01	9100	900
    'CIOT-16'	2445.24	9100	1200
    'CIOT-17'	2706.84	8600	900};
    
    % depth profile 26Al concentrations and depths
    % {sample_name  dz26(g/cm2)  N26(at/g)  N26unc(at/g)}
    depth_profile26 = ...
    {'CIOT-02'	124.37	512500	36600
    'CIOT-03'	233.66	273000	25300
    'CIOT-04'	358.1	182400	22400
    'CIOT-05'	466.79	178700	19800
    'CIOT-07'	700.89	106700	13600
    'CIOT-08'	768.7	92500	13100
    'CIOT-09'	879.54	97000	16000
    'CIOT-10'	1093.39	107600	18500
    'CIOT-15'	2238.01	48000	11900
    'CIOT-16'	2445.24	43800	9300
    'CIOT-17'	2706.84	34400	7200};
end;

if exist('depth_profile');
    out.sample_name = depth_profile(:,1);       % sample name
    out.Nz10 = cell2mat(depth_profile(:,2));    % sample midpoint depth (g/cm2)
    out.N10 = cell2mat(depth_profile(:,3));     % 10Be conc (atoms/g)
    out.N10unc = cell2mat(depth_profile(:,4));  % 10Be conc uncertainty (atoms/g)
    out.Nz26 = cell2mat(depth_profile(:,5));    % sample midpoint depth (g/cm2)
    out.N26 = cell2mat(depth_profile(:,6));     % 26Al conc (atoms/g)
    out.N26unc = cell2mat(depth_profile(:,7));  % 26Al conc uncertainty (atoms/g)
else;
    out.sample_name10 = depth_profile10(:,1);     % sample name
    out.Nz10 = cell2mat(depth_profile10(:,2));    % sample midpoint depth (g/cm2)
    out.N10 = cell2mat(depth_profile10(:,3));     % 10Be conc (atoms/g)
    out.N10unc = cell2mat(depth_profile10(:,4));  % 10Be conc uncertainty (atoms/g)
    out.sample_name26 = depth_profile26(:,1);     % sample name
    out.Nz26 = cell2mat(depth_profile26(:,2));    % sample midpoint depth (g/cm2)
    out.N26 = cell2mat(depth_profile26(:,3));     % 26Al conc (atoms/g)
    out.N26unc = cell2mat(depth_profile26(:,4));  % 26Al conc uncertainty (atoms/g)
end;
