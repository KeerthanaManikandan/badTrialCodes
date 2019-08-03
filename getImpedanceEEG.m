function [electrodeNames,electrodeImpedances] = getImpedanceEEG(filename)
    % Read BrainAmp impedance values from text file
    %
    % Modified from script auto-generated by MATLAB on 2018/01/03 12:29:29
    % Format string for each line of text:
    formatSpec = '%1s%3s%9s%5s%12s%2s%3s%s%[^\n\r]';

    % Open the text file.
    fileID = fopen(filename,'r');

    % Read columns of data according to format string.
    textscan(fileID, '%[^\n\r]', 0, 'WhiteSpace', '', 'ReturnOnError', false);
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'ReturnOnError', false);

    % Remove white space around all cell columns.
    dataArray{1} = strtrim(dataArray{1});
    dataArray{2} = strtrim(dataArray{2});
    dataArray{4} = strtrim(dataArray{4});
    dataArray{8} = strtrim(dataArray{8});

    % Close the text file.
    fclose(fileID);
    clear fileID

    firstColumn = dataArray{:, 1};
    startRow = find(strcmp(firstColumn,'#'),1,'first');
    endRow = find(strcmp(firstColumn,'#'),1,'last');
    electrodeNames = dataArray{:, 4}(startRow:endRow);
    electrodeImpedances = cellfun(@str2num,dataArray{:, 8}(startRow:endRow),'uniformOutput',0);
    electrodeImpedances(cellfun(@isempty,electrodeImpedances)) = {NaN};

%     electrodeImpedances = [cellfun(@str2num,dataArray{:, 2}(startRow:endRow)) cell2mat(electrodeImpedances)];
    electrodeImpedances = cell2mat(electrodeImpedances);
end