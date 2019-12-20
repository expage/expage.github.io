function elv = isost_elv(isost,sample);

% Function for calculating isostatic elevation development.
% This is free software: you can use/copy/modify/distribute as long as you keep it free/open.
% Jakob Heyman (jakob.heyman@gu.se) 2019

% check method to use (ICE6G, ANU2017, PD2015, or data from file) and fix elv
if isfield(isost,'delv_500'); % ICE6G or ANU2017
    elv = isost_interp(isost,sample.lat,sample.long,sample.elv,sample.tv,sample.samplingyr);
elseif isfield(isost,'Bv'); % Påsse and Daniels (2015)
    elv = isost_PD2015(isost,sample.lat,sample.long,sample.elv,sample.tv,sample.samplingyr);
elseif isfield(isost,'tv'); % data from file
    elv = isost_file(isost,sample.tv,sample.elv,sample.samplingyr);
else;
    fprintf(1,'\nIsostatic adjustment input ERROR!\n');
    elv = sample.elv;
end;

% make row vector of elv
elv = elv(:)';


% subfunction isost_interp ==========================================================================
function elv = isost_interp(indata,lat,long,elv0,tv,samplingyr);
    % interpolate uplift from the input data
    % fill elv vector
    tv = tv(:);
    elv(1:numel(tv),1) = elv0;
    % uplift time vector
    indata.tv = indata.tv + samplingyr - indata.tv0;
    % field names for loops below
    Tnames = indata.Tnames;
    % fix for coordinates
    if strcmp(indata.name,'ICE6G');
        [indata,lat,long] = latlongfix(indata,Tnames,lat,long);
    elseif lat<min(indata.lat) || lat>max(indata.lat) || ...
        long<min(indata.lon) || long>max(indata.lon);
        fprintf(1,'\nthe sample is located outside the %s domain!\n',indata.name);
        elv = elv0; return;
    end;
    % interpret uplift
    delv = [0];
    for i = 1:numel(Tnames);
        delv(end+1) = interp2(indata.lon,indata.lat,indata.(Tnames{i}),long,lat,'pchip');
    end;
    % find index of tv within indata.tv (0-26 ka BP) and add interpolated elevation difference
    tvidx = find(tv>=indata.tv(1) & tv<=indata.tv(end));
    elv(tvidx) = elv(tvidx) + interp1(indata.tv,delv,tv(tvidx),'pchip');
% end subfunction isost_interp ======================================================================


% subfunction isost_PD2015 =========================================================================
function elv = isost_PD2015(PD2015,lat,long,elv0,tv,samplingyr);
    % interpolate uplift for Fennoscandia based on Påsse and Daniels (2015)
    % convert lat long to Sweref99 TM coordinates X and Y and check that they match the domain
    [X Y] = sweref(lat,long);
    if X<min(PD2015.Xv) || X>max(PD2015.Xv) || Y<min(PD2015.Yv) || Y>max(PD2015.Yv);
        fprintf(1,'\nthe sample is located outside the PD2015 domain!\n');
        elv = elv0; return;
    end;
    % interpolate A,T,Bv,Be
    A = interp2(PD2015.Xv,PD2015.Yv,PD2015.A,X,Y,'pchip');
    T = interp2(PD2015.Xv,PD2015.Yv,PD2015.T,X,Y,'pchip');
    Bv = interp2(PD2015.Xv,PD2015.Yv,PD2015.Bv,X,Y,'pchip');
    Be = interp2(PD2015.Xv,PD2015.Yv,PD2015.Be,X,Y,'pchip');
    % relate tv to 1950
    tvpd =  tv - samplingyr + 1950;
    % calculate U and E
    U = 2./pi.*A.*(atan(T/Bv)-atan((T-tvpd)./Bv)) + 2./pi.*A.*(atan(T/Be)-atan((T-tvpd)./Be));
    E = 2./pi.*29.*(atan(9200/1300)-atan((9200-tvpd)./1300)) + ...
        2./pi.*36.*(atan(13700/1350)-atan((13700-tvpd)./1350));
    elv = U-E;
    elv = elv0 - (elv-elv(1));
% end subfunction isost_PD2015 =====================================================================


% subfunction sweref ===============================================================================
function [X,Y] = sweref(lat,long)
    %Converts coordinates in lat/long WGS84 to SWEREF99TM
    lat = lat * pi/180;
    long = long * pi/180;
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
    lat1 = lat - sin(lat)*cos(lat)*(A + B*sin(lat)^2 + C*sin(lat)^4 + D*sin(lat)^6);
    es=atan(tan(lat1)/cos(dlong));
    ns=atanh(cos(lat1)*sin(dlong));
    Y = k*at*(es + b1*(sin(2*es)*cosh(2*ns)) + b2*(sin(4*es)*cosh(4*ns)) + ...
        b3*(sin(6*es)*cosh(6*ns)) + b4*(sin(8*es)*cosh(8*ns))) + FN;
    X = k*at*(ns + b1*(cos(2*es)*sinh(2*ns)) + b2*(cos(4*es)*sinh(4*ns)) + ...
        b3*(cos(6*es)*sinh(6*ns)) + b4*(cos(8*es)*sinh(8*ns))) + FE;
% end subfunction sweref ===========================================================================


% subfunction isost_file ===========================================================================
function elv = isost_file(file,tv,elv0,samplingyr);
    % subfunction for fixing elevation based on specific file
    tv = tv(:);
    elv(1:numel(tv),1) = elv0;
    if isfield(file,'delv')==0 && isfield(file,'shoreline');
        file.delv = -file.shoreline;
    end;
    if isfield(file,'delv');
        % assume that the isostatic adjustment tv starts at sampling year if not given in yr0
        if isfield(file,'yr0')==0;
            file.yr0 = samplingyr;
        end;
        % fix file.tv
        file.tv = file.tv + samplingyr - file.yr0;
        % find index of tv within file.tv and add interpolated elevation difference
        tvidx = find(tv>=file.tv(1) & tv<=file.tv(end));
        elv(tvidx) = elv(tvidx) + interp1(file.tv,file.delv,tv(tvidx),'pchip');
    else;
        fprintf(1,'\nIsostatic adjustment input ERROR! (missing delv or shoreline)\n');
    end;
% end subfunction isost_file =======================================================================


% subfunction latlongfix ===========================================================================
function [indata,lat,long] = latlongfix(indata,Tnames,lat,long);
    % fix for lat > max indata.lat
    if lat > max(indata.lat);
        % use max values for 90 deg N.
        indata.lat(2:end+1) = indata.lat; indata.lat(1) = 90;
        % add first line at top of all delv matrices
        for i = 1:numel(Tnames);
            latfix = indata.(Tnames{i});
            latfix(2:end+1,:) = latfix;
            indata.(Tnames{i}) = latfix;
        end;
    end;
    % fix for lat < min indata.lat
    if lat < min(indata.lat);
        % use min values for 90 deg S.
        indata.lat(end+1) = -90;
        % add last line at bottom of all delv matrices
        for i = 1:numel(Tnames);
            latfix = indata.(Tnames{i});
            latfix(end+1,:) = latfix(end,:);
            indata.(Tnames{i}) = latfix;
        end;
    end;
    % fix for long < min indata.lon or long > max indata.lon
    if long<min(indata.lon) || long>max(indata.lon)
        % add 360 degrees to long < min indata.long
        if long < min(indata.long); long = 360 + long; end;
        % add one step to long vector
        indata.lon(end+1) = indata.lon(end)*2-indata.lon(end-1);
        % add first column to end of all delv matrices
        for i = 1:numel(Tnames);
            longfix = indata.(Tnames{i});
            longfix(:,end+1) = longfix(:,1);
            indata.(Tnames{i}) = longfix;
        end;
    end;
% end subfunction latlongfix =======================================================================
