function isost = isost_data_load(isost);

% function for loading data for isostatic elevation ajustments
% This is free software: you can use/copy/modify/distribute as long as you keep it free/open.
% Jakob Heyman - 2015-2018 (jakob.heyman@gu.se)

% find all unique input isost
single_isost = unique(isost);

% run through and load/fix data
for i = 1:numel(single_isost);
    if strcmp(single_isost(i),'-'); % '-' marks no isostatic adjustment
        continue;
    elseif strcmpi(single_isost(i),'ICE6G'); % if using the ICE6G adjustment
        ICE6G = isost_data_ICE6G();
        ICE6Gv = strcmpi(isost,'ICE6G');
        isost(ICE6Gv) = ICE6G;
    elseif strcmpi(single_isost(i),'ANU2017'); % if using the ANU2017 adjustment
        ANU2017 = isost_data_ANU2017();
        ANU2017v = strcmpi(isost,'ANU2017');
        isost(ANU2017v) = ANU2017;
    elseif strcmp(single_isost(i),'PD2015'); % if using Påsse and Daniels (2015) model
        PD2015 = isost_data_PD2015();
        PD2015v = strcmpi(isost,'PD2015');
        isost(PD2015v) = PD2015;
    else; % if using a self-provided file
        % read file
        fid = fopen(single_isost{i});
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
            outstruct.(infields{j}) = str2num(char(indata(varnum(j)+1:varnum(j+1)-1)));
        end;
        filev = strcmpi(isost,single_isost(i));
        isost(filev) = outstruct;
        % clear structure for any potential further run
        clear outstruct;
    end;
end;
