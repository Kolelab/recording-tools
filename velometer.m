% % % % % complement to nFit
% % % % % (c) Charles CH Cohen, 2014-present
% % % % % this software is released to the public under a Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 
% % % % % International license (CC BY-NC-ND 4.0, in English).
% % % % % for any questions, please email c.cohen@gmx.com

disp('Please import voltage response of interest')
uiwait(msgbox('Please import voltage response of interest', 'Import Data', 'modal'))

% % % % % startpath
mpath = mfilename('fullpath');
mfile = mfilename();
mpath = erase(mpath, mfile);

% % % % % create default path and recsetup files
if (~exist(strcat(mpath, 'path.txt'), 'file'))

    fileID = fopen(strcat(mpath, 'path.txt'), 'wt');
    fclose(fileID);
end

if (~exist(strcat(mpath, 'velosetup.txt'), 'file'))
    
    temp = 0;
    dlmwrite(strcat(mpath, 'velosetup.txt'), temp, 'delimiter', '\t');
    
end

% % % % read any delimited text file with numeric data
% % % % start with current path
fileID = fopen(strcat(mpath, 'path.txt'), 'rt');
path = fread(fileID, '*char')';
fclose(fileID);

if (isempty(path))

    [name, path] = uigetfile(strcat(mpath, '*.*'), 'Select file...');

else
    
    [name, path] = uigetfile(strcat(path, '*.*'), 'Select file...');
    
end

% % % % % write current path to path file
fileID = fopen(strcat(mpath, 'path.txt'), 'wt');
fwrite(fileID, path);
fclose(fileID);

% % % % % get pertinent velocity setup details
% % % % % start by checking if velosetup already defined
velosetup = dlmread(strcat(mpath, 'velosetup.txt'), '\t');

% % % % % if velosetup = 0, being initialized for the first time...
if (~sum(velosetup))

    velosetupqstr = {'Please enter desired dV/dt threshold (V/s).'};
    velosetupans = inputdlg(velosetupqstr, 'Velocity Setup');
    
    % % % % % create numeric velosetup
    velosetup = zeros(size(velosetupans, 1), 1);
    for i = 1:size(velosetupans, 1)        
        if (~isempty(velosetupans{i, 1}))
            velosetup(i, 1) = str2double(velosetupans{i, 1});
        else
            velosetup(i, 1) = 0;
        end
    end
    
    % % % % % label velosetup entries
    thresh = velosetup(1, 1);
    
    % % % % % write velosetup details to file
    dlmwrite(strcat(mpath, 'velosetup.txt'), velosetup, 'delimiter', '\t');

else

    % % % % thresh keeps prev value
    % % % % % label velosetup entries
    thresh = velosetup(1, 1);
    
    % % % % % user decides whether to keep thresh value
    keepveloqstr = {strcat('Threshold dV/dt (V/s)= ', num2str(thresh))};
    keepveloans = questdlg(keepveloqstr, 'Keep velocity setup?', 'Yes', 'No', 'Yes');
    
    switch keepveloans
        
        case 'No'
            
            velosetupqstr = {'Please enter desired dV/dt threshold (V/s).'};
            velosetupans = inputdlg(velosetupqstr, 'Velocity Setup');
            
            % % % % % create numeric velosetup
            velosetup = zeros(size(velosetupans, 1), 1);
            for i = 1:size(velosetupans, 1)
                if (~isempty(velosetupans{i, 1}))
                    velosetup(i, 1) = str2double(velosetupans{i, 1});
                else
                    velosetup(i, 1) = 0;
                end
            end
            
            % % % % % label velosetup entries
            thresh = velosetup(1, 1);
            
            % % % % % write velosetup details to file
            dlmwrite(strcat(mpath, 'velosetup.txt'), velosetup, 'delimiter', '\t');
    end
end


s = uiimport(fullfile(path,name));

t = s.data(:, 1);

vimp = s.data(:, 2:end);

% % % % create dv/dt
dvdt = diff(vimp)./diff(t);
td = t(2:end);

% % % create cell array for intersections
cint = cell(1, size(dvdt, 2));

% % % % % create thresh array matching vimp and t
thresha = zeros(size(td, 1), 1);

for i = 1:size(thresha, 1)
    thresha(i, 1) = thresh;
end

% % % % % store all intersections in sint
for i = 1:size(dvdt, 2)
    cint{1, i} = intersections(td, thresha, td, dvdt(:, i));
end

% % % % % store results in another cell array
cresults = cell(2, 3);
cresults{1, 1} = 'Loc';
cresults{1, 2} = 'Nseg';
cresults{1, 3} = 'Thresh';

minta = zeros(1, size(cint, 2));
for i = 1:size(cint, 2)
    minta(1,i) = cint{1,i}(1,1);
end

for i = 1:size(minta, 2)
    mint = min(minta);
end

for i = 1:size(minta, 2)
    if abs(mint-minta(1,i)) < 1e-7
        loc = i;
    end
end
        
cresults{2, 1} = loc;
cresults{2, 2} = size(minta, 2);
cresults{2, 3} = minta(1,loc);
