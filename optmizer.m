function val = getLastMonitorValue(filename)
% GETLASTMONITORVALUE  Return the second column of the last row of a
% monitor CSV as a scalar.
%
%   val = getLastMonitorValue('Mass Flow Averaged 1 Monitor Plot.csv')
%
%   The file is expected to have a header line and rows of the form:
%       iteration, value
%   Only the tail of the file is read, so this is cheap even for large files.

    if nargin < 1 || isempty(filename)
        filename = 'Mass Flow Averaged 1 Monitor Plot.csv';
    end

    fid = fopen(filename, 'r');
    if fid < 0
        error('getLastMonitorValue:openFailed', ...
              'Could not open file: %s', filename);
    end
    cleaner = onCleanup(@() fclose(fid));   %#ok<NASGU>

    % Jump to end and read a small tail chunk
    fseek(fid, 0, 'eof');
    fileSize = ftell(fid);
    chunkSize = min(4096, fileSize);
    fseek(fid, -chunkSize, 'eof');
    chunk = fread(fid, chunkSize, '*char').';

    % Split into lines, discard empty ones
    lines = regexp(chunk, '\r?\n', 'split');
    lines = lines(~cellfun(@(s) isempty(strtrim(s)), lines));
    if isempty(lines)
        error('getLastMonitorValue:emptyFile', 'File is empty: %s', filename);
    end

    % Parse the last line
    lastLine = strtrim(lines{end});
    vals = str2double(strsplit(lastLine, ','));
    if numel(vals) < 2 || any(isnan(vals(1:2)))
        error('getLastMonitorValue:parseFailed', ...
              'Could not parse numeric values from last line: "%s"', lastLine);
    end

    val = vals(2);   % scalar
end

function solution = start_simulation(starExe, macroFile)

    if nargin < 2 || isempty(starExe) ||isempty(macroFile)
        starExe   = 'D:\Program Files\Siemens\20.06.010-R8\STAR-CCM+20.06.010-R8\star\lib\win64\clang17.0vc14.2-r8\lib\starccmw.exe';
        macroFile   = 'E:\temp\blue\difusor\optimizer test\setting_up.java';
    end



    % Number of parallel processes 
    nProcs = 16;
    cmd = sprintf('"%s" -new -np %d -macro "%s"',starExe, nProcs, macroFile);
    %cmd = sprintf('"%s" -new -macro "%s"', starExe, macroFile);

    [status, cmdout] = system(cmd);
    disp(cmdout);
    if status ~= 0
        disp("abnormal exiting");
    end
    solution = getLastMonitorValue();
end


disp(start_simulation());