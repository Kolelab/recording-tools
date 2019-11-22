% % % % % complement to nFit
% % % % % (c) Charles CH Cohen, 2014-present
% % % % % this software is released to the public under a Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 
% % % % % International license (CC BY-NC-ND 4.0, in English).
% % % % % for any questions, please email c.cohen@gmx.com

disp('Please import voltage responses from each recording location (first col s, remaining V)')
uiwait(msgbox('Please import voltage responses from each recording location (first col s, remaining V)', 'Import Data', 'modal'))


% % % % % startpath
mpath = mfilename('fullpath');
mfile = mfilename();
mpath = erase(mpath, mfile);


% % % % % create default path and recsetup files
if (~exist(strcat(mpath, 'path.txt'), 'file'))

    fileID = fopen(strcat(mpath, 'path.txt'), 'wt');
    fclose(fileID);
end

if (~exist(strcat(mpath, 'recsetup.txt'), 'file'))
    
    temp = 0;
    dlmwrite(strcat(mpath, 'recsetup.txt'), temp, 'delimiter', '\t');
    
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

% % % % % get pertinent recording setup details
% % % % % start by checking if recsetup already defined
recsetup = dlmread(strcat(mpath, 'recsetup.txt'), '\t');

% % % % % if recsetup = 0, being initialized for the first time...
if (~sum(recsetup))

    recsetupqstr = {'Please enter the amlifier filter frequency applied at recording (kHz), if applicable. Otherwise leave empty.', 'Please enter the delay (ms) to current injection', 'Please enter the number of injection protocols'};
    recsetupans = inputdlg(recsetupqstr, 'Recording Setup');
    
    % % % % % create numeric recsetup
    recsetup = zeros(size(recsetupans, 1), 1);
    for i = 1:size(recsetupans, 1)        
        if (~isempty(recsetupans{i, 1}))
            recsetup(i, 1) = str2double(recsetupans{i, 1});
        else
            recsetup(i, 1) = 0;
        end
    end
    
    % % % % % label recsetup entries
    ff = recsetup(1, 1);
    delay = recsetup(2, 1);
    ninj = recsetup(3, 1);
    
    % % % % % write recsetup details to file
    dlmwrite(strcat(mpath, 'recsetup.txt'), recsetup, 'delimiter', '\t');
    
else

    % % % % ff, delay and ninj keep prev values
    % % % % % label recsetup entries
    ff = recsetup(1, 1);
    delay = recsetup(2, 1);
    ninj = recsetup(3, 1);
    
    % % % % % user decides whether to keep rec values
    keeprecqstr = {strcat('Amplifier filter frequency (kHz)= ', num2str(ff)); strcat('Delay to injection (ms)= ', num2str(delay)); strcat('Injection protocols number= ', num2str(ninj))};
    keeprecans = questdlg(keeprecqstr, 'Keep recording setup?', 'Yes', 'No', 'Yes');
    
    switch keeprecans
        
        case 'No'
            
            recsetupqstr = {'Please enter the amlifier filter frequency applied at recording (kHz), if applicable. Otherwise leave blank.', 'Please enter the delay (ms) to current injection', 'Please enter the number of injection protocols'};
            recsetupans = inputdlg(recsetupqstr, 'Recording Setup');
            
            % % % % % create numeric recsetup
            recsetup = zeros(size(recsetupans, 1), 1);
            for i = 1:size(recsetupans, 1)
                if (~isempty(recsetupans{i, 1}))
                    recsetup(i, 1) = str2double(recsetupans{i, 1});
                else
                    recsetup(i, 1) = 0;
                end
            end
            
            % % % % % label recsetup entries
            ff = recsetup(1, 1);
            delay = recsetup(2, 1);
            ninj = recsetup(3, 1);
            
            % % % % % write recsetup details to file
            dlmwrite(strcat(mpath, 'recsetup.txt'), recsetup, 'delimiter', '\t');
    end
end


s = uiimport(fullfile(path,name));

t = s.data(:, 1);

vimp = s.data(:, 2:end);


% % % % % remove leftover traces (below protocol number), if any
startcol = size(vimp, 2);
endcol = startcol;
while (mod(endcol, ninj))
    endcol = endcol - 1;
end
vimp = vimp(:, 1:endcol);


% % % % % look at an average based on surpassing applied recording filter
% % % % % frequency for reducing false positive signal the period in seconds
T = t(2,1) - t(1,1);

% % % % % the sampling frequency in Hz
sf = round(1./T);

% % % % the amplifier filter frequency in Hz, if applicable
if (ff)
    
    ff = ff*1000;
    
else
    
    ff = sf;
    recsetup(1,1) = ff;
    % % % % % update recsetup.txt
    dlmwrite(strcat(mpath, 'recsetup.txt'), recsetup, 'delimiter', '\t');
end

% % % % obtain downsampling frequency ratio
temp = sf/ff;

% % % % % set minimum number of data points over which to do stats
dmin = 30;

% % % % Set dminf as function of ff
if ((temp/2 - round(temp)/2) < sf)
    
    dminf = round(temp) + 1;
    
else
    
    if ((temp - round(temp)) < 0.5)
    
        dminf = round(temp);
        
    else
        
        dminf = round(temp) + 2;
        
    end
end


% % % % % % get baseline, average of 70% of delay time to injection
ibl = int64(0.7*delay*1e-3*sf);

% % % % create t_ms from t (in s)
t_ms = 1000*t;


% % % % % % create vsort from v sorted by ninj
vsort = cell(ninj, 1);

i = 1;
ninjcol = int64(size(vimp, 2)./ninj);

while (i <= ninj)
    
    j = i;
    k = 1;
    
    while (j <= size(vimp, 2) && k <= ninjcol)
        
        vsort{i, 1}(:, k) = vimp(:, j);
        
        j = j + ninj;
        k = k + 1;
        
    end
    
    i = i + 1;
    
end


% % % % % create vout for output v
vout = zeros(size(t, 1), ninj);


% % % % % create nsigv matrix, for storing
% % % % % the precision of v over ninj
nsigv = zeros(ninj, 1);


% % % % % begin processing each amp in vsort
for inj = 1:ninj
    
v = vsort{inj, 1}(:, :);

% % % % get a bl the size of v
bl = zeros(size(v, 1), size(v, 2));
bl(1, :) = mean(v(1:ibl, :), 1);

for i = 2:size(bl, 1)
    bl(i, :) = bl(1, :);
end

% % % % % for bookeeping purposes, an avg bl is created
bl_avg = mean(mean(bl, 2));

% % % % v - bl is created here
vbl = v - bl;


% % % % % % determine range over which to test for noisy signal
if (delay > 10)

    start = delay - 10;
    istart = int64(start*1e-3*sf);

else

    start = delay;
    istart = int64(delay);

end


% % % % % % determine if signal is positive or negative
count_pos = 0;
count_neg = 0;

idelay = int64(delay*1e-3*sf);

vbl_avg = mean(vbl, 2);
vbl_bl = mean(vbl_avg(1:idelay, 1));

% % % % Search for likely signal end
j = size(vbl_avg, 1)-dmin;
iend = size(v, 1);

while (j > idelay)

    if (vbl_avg(j-dmin:j) > mean(vbl_avg(j+1:size(vbl_avg, 1))) + 2*std(vbl_avg(j+1:size(vbl_avg, 1))))

        iend = j;
        break

    else if (vbl_avg(j-dmin:j) < mean(vbl_avg(j+1:size(vbl_avg, 1))) - 2*std(vbl_avg(j+1:size(vbl_avg, 1))))

            iend = j;
            break

        else

            j = j-1;

        end
    end
end


% % % % % % count n data points above or below baseline, to determine if
% % % % % % signal is overall pos or neg
j = idelay+1;
while (j < iend-dminf)

    if (vbl_avg(j:j+dminf) > vbl_bl)

        count_pos = count_pos + 1;

        j = j+dminf+1;

    else if (vbl_avg(j:j+dminf) < vbl_bl)

            count_neg = count_neg + 1;

            j = j+dminf+1;

        else

            j = j+1;

        end
    end
end

if (count_pos > count_neg)

    spos = 1;
    sneg = 0;

else

    spos = 0;
    sneg = 1;

end


% % % % % remove noisy traces
vbl_sd = std(vbl, 0, 2);

noise_count = zeros(1, size(vbl, 2));
if (spos)

    i = 1;
    while (i <= size(vbl, 2))

        count = 0;

        j = 1;
        while (j <= size(vbl, 1)-dminf)

            if (vbl(j:j+dminf, i) > vbl_avg(j:j+dminf, 1) + 2*vbl_sd(j:j+dminf))

                count = count + dminf;

                j = j + dminf + 1;

            else

                j = j + 1;

            end
        end

        noise_count(1, i) = count;

        i = i + 1;

    end

else

    i = 1;
    while (i <= size(vbl, 2))

        count = 0;

        j = 1;
        while (j <= size(vbl, 1)-dminf)

            if (vbl(j:j+dminf, i) < (vbl_avg(j:j+dminf, 1) - 2*vbl_sd(j:j+dminf)))

                count = count + dminf;

                j = j + dminf + 1;

            else

                j = j + 1;

            end
        end

        noise_count(1, i) = count;

        i = i + 1;

    end
end


% % % % % Remove noisy traces; those which have more than 5% noise.
alpha = 0.05;

keep_count = zeros(1, size(noise_count, 2));

for i = 1:size(keep_count, 2)

    if (noise_count(1, i) < 0.05*size(vbl, 1))

        keep_count(1, i) = 1;

    end

end

keep_counter = sum(keep_count);


keep_vbl = zeros(size(v, 1), keep_counter);
keep_bl = zeros(1, keep_counter);

i = 1;
count = 1;
while (i <= size(vbl, 2))

    if (keep_count(1, i) > 0)

        keep_vbl(:, count) = vbl(:, i);

        keep_bl(1, count) = mean(bl(1:ibl, i));

        count = count + 1;

        i = i + 1;

    else

        i = i + 1;

    end

end

% % % baselined, kept avg v in V
keep_vbl_avg = mean(keep_vbl, 2);

keep_bl_avg = mean(keep_bl, 2);

temp = zeros(size(keep_vbl_avg, 1), 1);

temp(:, :) = keep_bl_avg;

keep_bl_avg_vec = temp;

keep_vbl_avg_RMP = keep_vbl_avg + keep_bl_avg_vec;

% % % % % % Defines the junction potential in V, included in the final trace
JP = -0.015;

temp = zeros(size(keep_vbl_avg_RMP, 1), 1);

temp(:, :) = JP;

JP_vec = temp;

keep_vbl_avg_RMP_JP = keep_vbl_avg_RMP + JP_vec;

RMP = keep_bl_avg*1000;

RMP_JP = RMP + JP*1000;

% % % % % % % Convert to ms and mV with JP
keep_vbl_avg_mv = keep_vbl_avg*1000;
keep_vbl_avg_RMP_JP_mv = 1000*keep_vbl_avg_RMP_JP;

keep_count_str = [' ' num2str(keep_counter)];
total_count_str = [' ' num2str(size(v, 2))];

% % % % get nsig of v; use first column, otherwise excessive runtime
nsig = zeros(size(v, 1), 1);
for j = 1:size(v, 1)
    nsig(j, 1) = sigdigits(v(j, 1));
end

% % % % % keep max nsig, the precision of each v, in nsigmax
nsigv(inj, 1) = floor(mean(max(nsig)));

% % % % % set nsig of keep_vbl_avg_mv to max(nsig)
keep_vbl_avg_mv_round = roundsd(keep_vbl_avg_mv, nsigv(inj, 1));
vout(:, inj) = keep_vbl_avg_mv;

end


% % % % % % Create output file at same location
nameOut = strcat(strtok(name, '.'), ' avg NTR.dat');

pathnameOut = strcat(path,nameOut);

vtout = zeros(size(t, 1), size(vout, 2) + 1);

vtout(:, 1) = t_ms;
vtout(:, 2:end) = vout;

nrows = size(vtout, 1);
ncols = size(vtout, 2);

matheader = [nrows, ncols];

dlmwrite(pathnameOut, matheader, 'delimiter', ' ');
dlmwrite(pathnameOut, vtout, '-append', 'delimiter', '\t', 'precision', floor(mean(max(nsigv))));

disp(['Voltage response for ', strtok(name, '.'), ' exported: ']);
disp(strcat(' ', pathnameOut));

clear
