function [data] = GEM2FEM(filename)
%this function converts data from the standard GEM2 format to FEMIC
%this will also read the data if the file is already in FEMIC format
%which would be the case for previously filtered data


%file = dlmread('EMline1raw.csv')

%% Open the data file.
fileID = fopen(filename,'r'); %give file an identifier
header = fgets(fileID); %read in header as single charcter string
colnames = strread(header,'%s','delimiter',','); %split the string based on commas
ncol = length(colnames); %number of columns in the header

%find the file extension - currently .csv is used for the unfiltered and .dat is used for the filtered
ext = []; %start with empty vector
for index=length(filename):-1:1 %loop backwards through filename
  if (filename(index) == '.') %break out of loop when '.' is encountered
    break;
  end
  ext = strcat(filename(index),ext); %while in loop put letters of extension in 'ext'
end

if (ext == 'csv') %then this is a raw GEM2 file
	headerSpec = repmat('%s',1,ncol); %26 strings to be read in header line
	dataSpec = repmat('%f64',1,ncol); %26 strings to be read in header line
	dataArray= textscan(fileID, dataSpec, 'delimiter', ','); %data
	sizedata = size(cell2mat(dataArray)); %number of rows and columns in the data
	nrows = sizedata(1);
	fclose(fileID); %% Close the text file.
save('dataArray.mat','dataArray');
%find the frequency info
freq_find = strfind(header, 'Hz'); %find number of frequencies /2 for imaginary and quadrature components
nfreq = length(freq_find);

freq = zeros(1,nfreq); %holds the actual frequencies
k = 1; %counting variable

for ii = freq_find
	fr = []; %start with empty vector to hold character string
	for jj = (ii-1):-1:(ii-10) %search backwards from 'Hz' to a maximum of 10 characters to get frequency
	 if (header(jj) == '_') %break out of loop when '_' is encountered
	    break;
 	 end
	  fr = strcat(header(jj),fr); %while in loop put letters of extension in 'fr'
	end
	f(k) = str2num(fr); %convert character string to numeric
	k = k + 1;
end


f = unique(f); %get rid of repeats (since there are quadrature and imaginary components to each freq)
nfreq = length(f); %the total number of frequencies
freq = zeros(nfreq,nrows);
%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
X = repmat(dataArray{:, 3},1,nfreq)';
Y = repmat(dataArray{:, 4},1,nfreq)';
freq = repmat(f,1,nrows);
id = repmat(1:nrows,nfreq,1);

k = 12;
dt = [];
for ii=1:nfreq
	inphase = dataArray{:,k};
	k = k + 1;
	quadrature = dataArray{:,k};
	k = k + 1;
	d = complex(inphase, quadrature);
	dt = [dt d];
end

dt = dt';
data(:,1)=X(:);
data(:,2)=Y(:);
data(:,3)=zeros(length(data(:,1)),1);
data(:,4)=id(:);
data(:,5)=freq(:);
data(:,6)=ones(length(data(:,1)),1);
data(:,7)=(2.*data(:,6));
data(:,8)=dt(:);
%dlmwrite('data.dat','data');
save('dataArray.mat','dataArray');
elseif (ext == 'dat') %then this is a filtered file already in FEMIC format
	data = csvread(filename); %simple data read in
	nfreq = unique(data(:,5)); %number of frequencies
	%renumber soundings
	[uniqueid,I,I2] = unique(data(:,4));
	data(:,4) = I2; 
	%sort data by sounding ID, then frequency to make things simpler
	[sortdata,I] = sortrows(data(:,[4 5]));
	data = data(I,:);
	fclose(fileID); %% Close the text file.
else
	msgbox('The file is in an uknown format');
end
