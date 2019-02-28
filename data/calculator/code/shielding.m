function shielding(demfile,samplefile);

% Function calculating topographic shielding based on a DEM
%
% Input:
% demfile = name of dem file
% samplefile = name of sample file
% Example: shielding('dem.asc','samples.txt')
%
% The demfile should have a projected coordinate system with unit meters for x, y, and elevation.
% It should be an ascci file with 6 or 7 header lines (from QGIS or ESRI):
%   1. ncols (or NCOLS)
%   2. nrows (or NROWS)
%   3. xllcorner (or XLLCORNER / XLLCENTER / xllcenter)
%   4. yllcorner (or YLLCORNER / YLLCENTER / yllcenter)
%   5. dx (or CELLSIZE)
%   6. dy (or NODATA_VALUE)
%   7. NODATA_value (or nothing if only 6 header lines)
%
% The samplefile should contain only data (no header) in three columns:
%   1. sample name
%   2. x coordinate (unit: m - same projection as the demfile)
%   3. y coordinate (unit: m - same projection as the demfile)
%
% Output:
% Topographic shielding values are displayed and saved in out-shielding.txt
%
% This is free software: you can use/copy/modify/distribute as long as you keep it free/open.
% Jakob Heyman - 2019 (jakob.heyman@gu.se)


% Fixed Earth radius (m) used for curvature calculation
Rearth = 6371E3;

% read and fix dem input
fprintf(1,'reading DEM file...');
fid = fopen(demfile);
hdr = textscan(fid,'%s %f',6);
hdrvars = hdr{1,1};
hdrnum = hdr{1,2};
ncols = hdrnum(1);
nrows = hdrnum(2);
xll = hdrnum(3);
yll = hdrnum(4);
% fix for cell center
if strcmp(hdrvars{3},'XLLCENTER') || strcmp(hdrvars{3},'xllcenter');
    xll = xll-dx/2; yll = yll-dy/2;
end;
% fix for CELLSIZE
if strcmp(hdrvars{5},'CELLSIZE') || strcmp(hdrvars{5},'cellsize');
    dx = hdrnum(5); dy = hdrnum(5);
else;
    dx = hdrnum(5); dy = hdrnum(6);
    nodata = textscan(fid,'%s %f',1);
end;
fnum = [repmat(['%f '],1,ncols-1),'%f'];
indem = textscan(fid,fnum);
dem = cell2mat(indem);
fclose(fid);
fprintf(1,' done!\n');

% fix X and Y vectors
Xv = (xll+dx/2:dx:xll+dx*ncols);
Yv = (yll+dy*(nrows-0.5):-dy:yll);

% read and fix samplefile
fid = fopen(samplefile);
samplein = textscan(fid,'%s %f %f');
sample = samplein{1,1};
samplex = samplein{1,2};
sampley = samplein{1,3};
fclose(fid);

% find index of samplex and sampley (x1,y1)
x1 = interp1(Xv,(1:1:numel(Xv)),samplex,'nearest');
y1 = interp1(Yv,(1:1:numel(Yv)),sampley,'nearest');

% angles for horizon (radians)
angles = [0:pi/180:(2*pi)];

% header for output
output(1,1:2) = {'sample','shielding'};

% loop for samples
for i = 1:numel(sample);
    % if sample is not within the dem: skip this round
    if isnan(x1(i)) || isnan(y1(i)) || x1(i)==1 || x1(i)==numel(Xv) || y1(i)==1 || y1(i)==numel(Yv);
        continue;
    end;
    
    % display sample name and fill output
    fprintf(1,'%s',sample{i});
    output(end+1,1) = sample{i};
    
    % find index of edge coordinates (x2,y2)
    dyv([1:91,272:361]) = Yv(1)-Yv(y1(i));
    dyv(92:271) = Yv(end)-Yv(y1(i));
    dx2 = tan(angles).*dyv;
    x2 = interp1(Xv,(1:1:numel(Xv)),Xv(x1(i))+dx2,'nearest');
    dxv(1:181) = Xv(end)-Xv(x1(i));
    dxv(182:361) = Xv(1)-Xv(x1(i));
    dy2 = dxv./tan(angles);
    y2 = interp1(Yv,(1:1:numel(Yv)),Yv(y1(i))+dy2,'nearest');
    x2(isnan(x2(1:181))) = numel(Xv); x2(isnan(x2)) = 1;
    y2(isnan(y2(1:91))) = 1; y2(isnan(y2(1:271))) = numel(Yv); y2(isnan(y2)) = 1;
    
    % sample elevation
    sampleelv = dem(y1(i),x1(i));
    
    % loop for azimuth angles
    for j = 1:numel(angles);
        % pick out index of azimuth lines
        [xl yl] = bresenham(x1(i),y1(i),x2(j),y2(j));
        linind = sub2ind(size(dem),yl,xl);
        % pick out elevations, calculate Earth curvature effect, and find max angle
        elv = dem(linind);
        distance = sqrt(((xl-x1(i)).*dx).^2 + ((yl-y1(i)).*dy).^2);
        curveminus = Rearth./cos(distance/Rearth)-Rearth;
        elvrad(j) = max(atan((elv-curveminus-sampleelv)./distance));
    end;
    
    % don't allow negative angles
    elvrad(elvrad<0) = 0;
    
    % calculate shielding
    shield = 1 - sum(sin(elvrad).^3.3./360);
    
    % display shielding and fill output
    fprintf(1,' \t%.5f\n',shield);
    output(end,2) = {num2str(shield,'%.5f')};
end;

% write out-shielding.txt
out = fopen('out-shielding.txt','w');
for i = 1:size(output,1);
    fprintf(out,'%s\t%s\n',output{i,:});
end;
fclose(out);


% subfunction bresenham ============================================================================
function [x y] = bresenham(x1,y1,x2,y2);
    % Bresenham line algorithm
    % (x1,y1): start position
    % (x2,y2): end position
    % [x y]: line coordinates from (x1,y1) to (x2,y2)
    x1=round(x1); x2=round(x2);
    y1=round(y1); y2=round(y2);
    dx=abs(x2-x1);
    dy=abs(y2-y1);
    steep=abs(dy)>abs(dx);
    if steep; t=dx; dx=dy; dy=t; end;
    % the main algorithm goes here.
    if dy==0;
        q=zeros(dx+1,1);
    else;
        q=[0;diff(mod([floor(dx/2):-dy:-dy*dx+floor(dx/2)]',dx))>=0];
    end;
    % and ends here.
    if steep;
        if y1<=y2; y=[y1:y2]'; else; y=[y1:-1:y2]'; end;
        if x1<=x2; x=x1+cumsum(q); else; x=x1-cumsum(q); end;
    else;
        if x1<=x2; x=[x1:x2]'; else; x=[x1:-1:x2]'; end;
        if y1<=y2; y=y1+cumsum(q); else; y=y1-cumsum(q); end;
    end;
% end subfunction bresenham ========================================================================
