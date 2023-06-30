function isost_plot(lat,long,tv,tv0,data);

% function for plotting elevation change or relative sea level back in time
% input data:
% lat - latitude (WGS84), can be vector
% long - longitude (WGS84), must be same size as lat
% tv - time vector (yr)
% tv0 - reference year of tv year 0
% data - optional string or cell array input to plot one or more isostatic adjustments that can be
%        ICE7G, ICE6G, ANU2017, PD2015, Gowan2021, or the name of a file with input data

% plot relative sea level instead of elevation change? (1 = yes)
rsl = 1;

% fix longitude values
long(long<0) = long(long<0) + 360;

% check if samples are within data domains for ANU2017 and PD2015
ANU2017test = (lat>=36).*(lat<=84.5).*(long>=220).*(long<=300);
[X,Y] = sweref(lat,long);
PD2015test = (X>=-181820).*(X<=1808180).*(Y>=5841012).*(Y<=8011012);

% load data
if nargin == 5;
    if isa(data,'char'); data = {data}; end;
    for i = 1:numel(data);
        if strcmp(data{i},'ICE7G'); ICE7G = isost_data_ICE7G();
        elseif strcmp(data{i},'ICE6G'); ICE6G = isost_data_ICE6G();
        elseif strcmp(data{i},'ANU2017'); ANU2017 = isost_data_ANU2017();
        elseif strcmp(data{i},'Gowan2021'); Gowan2021 = isost_data_Gowan2021();
        elseif strcmp(data{i},'PD2015'); PD2015 = isost_data_PD2015();
        elseif strfind(data{i},'PD2015_');
            PD2015mod = isost_data_PD2015();
            parstr = {'A','T','Bv','Be'};
            for j = 1:numel(parstr);
                PD2015mod = modify_PD2015(PD2015mod,data{i},parstr{j});
            end;
        elseif isfile(data{i});
            FILE = file_data(data{i});
            filename = data{i};
        else;
            fprintf(1,'\n%s is not an allowed isostatic adjustment input!\n',data{i});
            fprintf(1,'Allowed inputs are: ICE7G, ICE6G, ANU2017, PD2015, Gowan2021, ');
            fprintf(1,'or the file name of a valid isostatic adjustment file\n');
        end;
    end;
else;
    % load ICE7G data
    ICE7G = isost_data_ICE7G();
    % load ICE6G data
    ICE6G = isost_data_ICE6G();
    % load Gowan2021 data
    Gowan2021 = isost_data_Gowan2021();
    % load ANU2017 data if within domain
    if sum(ANU2017test) > 0; ANU2017 = isost_data_ANU2017(); end;
    % load PD2015 data if within domain
    if sum(PD2015test) > 0; PD2015 = isost_data_PD2015(); end;
end;

% fix sample structure
sample.elv = 0; % use 0 as reference elevation
sample.tv = tv;
sample.samplingyr = tv0;

% fix matrices
ICE7Gtvm = []; ICE6Gtvm = []; ANU2017tvm = []; PD2015tvm = []; PD2015modtvm = []; Gowan2021tvm = [];
ICE7Gplm = []; ICE6Gplm = []; ANU2017plm = []; PD2015plm = []; PD2015modplm = []; Gowan2021plm = [];

% loop for coordinates
for i = 1:numel(lat);
    sample.lat = lat(i);
    sample.long = long(i);
    
    % fill plot matrices
    tvm(1:numel(tv),i) = tv(:);
    if exist('ICE7G');
        ICE7Gtvm(1:numel(tv),end+1) = tv(:);
        ICE7Gplm(1:numel(tv),end+1) = isost_elv(ICE7G,sample)';
    end;
    
    if exist('ICE6G');
        ICE6Gtvm(1:numel(tv),end+1) = tv(:);
        ICE6Gplm(1:numel(tv),end+1) = isost_elv(ICE6G,sample)';
    end;
    if exist('ANU2017') && ANU2017test(i)==1;
        ANU2017tvm(1:numel(tv),end+1) = tv(:);
        ANU2017plm(1:numel(tv),end+1) = isost_elv(ANU2017,sample)';
    end;
    if exist('PD2015') && PD2015test(i)==1;
        PD2015tvm(1:numel(tv),end+1) = tv(:);
        PD2015plm(1:numel(tv),i) = isost_elv(PD2015,sample)';
    end;
    if exist('PD2015mod') && PD2015test(i)==1;
        PD2015modtvm(1:numel(tv),end+1) = tv(:);
        PD2015modplm(1:numel(tv),i) = isost_elv(PD2015mod,sample)';
    end;
    if exist('Gowan2021');
        Gowan2021tvm(1:numel(tv),end+1) = tv(:);
        Gowan2021plm(1:numel(tv),end+1) = isost_elv(Gowan2021,sample)';
    end;
end;

% fix and plot
leg = [];
legin = {};
figure;
hold on; box on;
if exist('ICE7G');
    if rsl==1; ICE7Gplm = -ICE7Gplm; end;
    legtemp = plot(ICE7Gtvm,ICE7Gplm,'color',[0.7 0.7 0.7]); legin(end+1) = {'ICE7G'};
    leg(end+1) = legtemp(1);
end;
if exist('ICE6G');
    if rsl==1; ICE6Gplm = -ICE6Gplm; end;
    legtemp = plot(ICE6Gtvm,ICE6Gplm,'color','black'); legin(end+1) = {'ICE6G'};
    leg(end+1) = legtemp(1);
end;
if exist('ANU2017');
    if rsl==1; ANU2017plm = -ANU2017plm; end;
    legtemp = plot(ANU2017tvm,ANU2017plm,'color','blue'); legin(end+1) = {'ANU2017'};
    leg(end+1) = legtemp(1);
end;
if exist('PD2015');
    if rsl==1; PD2015plm = -PD2015plm; end;
    legtemp = plot(PD2015tvm,PD2015plm,'color','red'); legin(end+1) = {'PD2015'};
    leg(end+1) = legtemp(1);
end;
if exist('PD2015mod');
    if rsl==1; PD2015modplm = -PD2015modplm; end;
    legtemp = plot(PD2015modtvm,PD2015modplm,'color',[0.6 0 0]); legin(end+1) = {'PD2015mod'};
    leg(end+1) = legtemp(1);
end;
if exist('Gowan2021');
    if rsl==1; Gowan2021plm = -Gowan2021plm; end;
    legtemp = plot(Gowan2021tvm,Gowan2021plm,'color','green'); legin(end+1) = {'Gowan2021'};
    leg(end+1) = legtemp(1);
end;
if exist('FILE');
    if rsl==1;
        if isfield(FILE,'shoreline'); plfile = FILE.shoreline;
        else; plfile = -FILE.delv; end;
    else;
        if isfield(FILE,'delv'); plfile = FILE.delv;
        else; plfile = -FILE.shoreline; end;
    end;
    legtemp = plot(FILE.tv,plfile,'color','magenta'); legin(end+1) = {filename};
    leg(end+1) = legtemp(1);
end;
    
set(gca(),'xdir','reverse');
xlabel('Time (yr)');
if rsl==1; ylab = 'Shoreline (m a.s.l.)'; else; ylab = 'Elevation change (m)'; end;
ylabel(ylab);
legend(leg,legin,'location','northwest');
set(gca,'layer','top'); % plot axis on top
hold off;


% subfunction sweref ===============================================================================
function [X,Y] = sweref(lat,long)
    %Converts coordinates in lat/long WGS84 to SWEREF99TM
    lat = lat .* pi/180;
    long = long .* pi/180;
    %Constants
    a = 6378137; f = 1/298.257222101; e2 = f*(2-f); n = f/(2-f); at = a/(1+n)*(1+n^2/4+n^4/64);
    A = e2; B = 1/6*(5*e2^2-e2^3); C = 1/120*(104*e2^3 -45*e2^4); D = 1/1260*(1237*e2^4);
    long_av = pi/180*15.0; k = 0.9996; FN = 0; FE = 500000;
    b1 = 1/2*n - 2/3*n^2 + 5/16*n^3 + 41/180*n^4;
    b2 = 13/48*n^2 - 3/5*n^3 + 557/1440*n^4;
    b3 = 61/240*n^3 - 103/140*n^4;
    b4 = 49561/161280*n^4;
    % converstion from lat,long to X,Y in GRS80 (SWEREF)
    % the small difference between WGS84 and GRS80 is neglected here
    dlong = long-long_av;
    % calculating 'conformal latitude'
    lat1 = lat - sin(lat).*cos(lat).*(A + B.*sin(lat).^2 + C.*sin(lat).^4 + D.*sin(lat).^6);
    es=atan(tan(lat1)./cos(dlong));
    ns=atanh(cos(lat1).*sin(dlong));
    Y = k.*at.*(es + b1.*(sin(2.*es).*cosh(2.*ns)) + b2.*(sin(4.*es).*cosh(4.*ns)) + ...
        b3.*(sin(6.*es).*cosh(6.*ns)) + b4.*(sin(8.*es).*cosh(8.*ns))) + FN;
    X = k.*at.*(ns + b1.*(cos(2.*es).*sinh(2.*ns)) + b2.*(cos(4.*es).*sinh(4.*ns)) + ...
        b3.*(cos(6.*es).*sinh(6.*ns)) + b4.*(cos(8.*es).*sinh(8.*ns))) + FE;
% end subfunction sweref ===========================================================================


% subfunction file_data ============================================================================
function out = file_data(name);
    % read file
    fid = fopen(name);
    indata = textscan(fid,'%s','CommentStyle','%');
    fclose(fid);
    indata = indata{1};
    infields = {'tv','delv','shoreline','yr0'}; % potential input fields
    varnum = []; rmidx = [];
    % find which input variables are in the input file
    for j = 1:numel(infields);
        if sum(strcmpi(indata,infields(j))) == 1;
            varnum(end+1) = find(strcmpi(indata,infields(j)));
        else;
            rmidx(end+1) = j;
        end;
    end;
    % remove variables that are not in the input file
    infields(rmidx) = []; 
    % sort infields and varnum
    [varnum,varidx] = sort(varnum);
    infields = infields(varidx);
    % add number for last line of file
    varnum(end+1) = numel(indata)+1;
    % fill fields with data from indata
    for j = 1:numel(infields);
        out.(infields{j}) = str2num(char(indata(varnum(j)+1:varnum(j+1)-1)));
    end;
% end subfunction file_data ========================================================================


% subfunction modify_PD2015 ========================================================================
% subfunction for modification of PD2015 parameters based on name
function PD2015 = modify_PD2015(PD2015,PD2015name,parstr);
    PD2015name_cell = strsplit(PD2015name,'_');
    if cell2mat(strfind(PD2015name_cell,parstr));
        idx = find(not(cellfun('isempty',strfind(PD2015name_cell,parstr))));
        mult = str2num(erase(PD2015name_cell{idx},parstr));
        PD2015.(parstr) = PD2015.(parstr) .* mult;
    end;
% end subfunction modify_PD2015 ====================================================================
