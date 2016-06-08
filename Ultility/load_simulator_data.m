%% Import inten data from text file.
% Initialize variables.
foldername = 'C:\Users\CaiJunhao\Desktop\simulator\';
filename = [foldername,'inten.txt'];
% Format string for each line of text:
formatSpec = '%10f%f%[^\n\r]';
% Open the text file.
fileID = fopen(filename,'r');
% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN, 'ReturnOnError', false);
% Close the text file.
fclose(fileID);
% Allocate imported array to column variable names
t = dataArray{:, 1};
inten = dataArray{:, 2};

%% Import trans data from text file.
filename = [foldername,'trans.txt'];
% Format string for each line of text:
formatSpec = '%10f%f%[^\n\r]';
% Open the text file.
fileID = fopen(filename,'r');
% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN, 'ReturnOnError', false);
% Close the text file.
fclose(fileID);
% Allocate imported array to column variable names
trans = dataArray{:, 2};

%% Import specf data from text file.
filename = [foldername,'specf.txt'];
% Format string for each line of text:
formatSpec = '%10f%f%[^\n\r]';
% Open the text file.
fileID = fopen(filename,'r');
% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN, 'ReturnOnError', false);
% Close the text file.
fclose(fileID);
% Allocate imported array to column variable names
c = 299792.458; 
f = dataArray{:, 1};        
lambda = c./f;
specf = dataArray{:, 2};

%% Import chirp data from text file.
filename = [foldername,'chirp.txt'];
% Format string for each line of text:
formatSpec = '%10f%f%[^\n\r]';
% Open the text file.
fileID = fopen(filename,'r');
% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN, 'ReturnOnError', false);
% Close the text file.
fclose(fileID);
% Allocate imported array to column variable names
t_chirp = dataArray{:, 1};
chirp = dataArray{:, 2};

%% Import phase data from text file.
filename = [foldername,'phase.txt'];
% Format string for each line of text:
formatSpec = '%10f%f%[^\n\r]';
% Open the text file.
fileID = fopen(filename,'r');
% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN, 'ReturnOnError', false);
% Close the text file.
fclose(fileID);
% Allocate imported array to column variable names
phase = dataArray{:, 2};

%% Clear temporary variables
clearvars filename formatSpec fileID dataArray ans foldername;