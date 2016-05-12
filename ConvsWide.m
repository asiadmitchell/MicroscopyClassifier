function[ConfusionMatrix, Error]=ConvsWide(confiles, widefiles, K)
%
%[CONFUSIONMATRIX, ERROR]= CONVSWIDE(CONFILES, WIDEFILES, K)              
%This function reads in two file groups (CONFILES and WIDEFILES) and a K
%then outputs the CONFUSIONMATRIX and ERROR using a Nearest Neighbor 
%classifier and a K-fold cross validation.                                     
%                                                 
%CONFILES and WIDEFILES are ‘.txt’ files that contain the path of all files
%associated with their respective group. Each file must be listed on a 
%separate line.                                               
%                                                 
%Example:                                             
%[ConfusionMaxtrix, Error] = ConvsWide(‘confiles.txt’, ‘widefiles.txt’, 10)
%
%   ConfusionMatrix =                                     
%       22 0                                              
%       0 47                                              
%   Error =                                               
%       1                                             
%
%
conf = textread(confiles, '%s', 'delimiter', '\n'); %Reads in ‘confiles’
numconf = length(conf);
conmtr = zeros(numconf, 1024);
ctr=1;
for i=1:numconf
  img = double(imread(conf{i}));
  FT = fft2(img);
  FT_shifted = fftshift(FT);
  FTsumcon = sum(FT_shifted);
  conmtr(i,:) = FTsumcon;
  ctr=ctr+1;
end
widef = textread(widefiles, '%s', 'delimiter', '\n'); %Reads in ‘widefiles’
numwidef = length(widef);
widemtr = zeros(numwidef, 1024);
for i=1:numwidef
  img = double(imread(widef{i}));
  FT = fft2(img);
  FT_shifted = fftshift(FT);
  FTsumwide = sum(FT_shifted); %Sums FT_shifted values in each column
  widemtr(i,:) = FTsumwide; %Matrix containing FT_shifted sums for all files containined in ‘widefile’ variable
  ctr=ctr+1;
end
tmtr = [conmtr; widemtr]; %Concatenates ‘widemtr’ and ‘conmatr’, stores in ‘tmtr’
cgroup = zeros(numconf,1); %Following lines assign group labels 
for i=1:numconf
    cgroup(i) = 'C'; %Assigns ‘C’ label to ‘confiles’
end
wgroup = zeros (numwidef, 1);
for i=1:numwidef
    wgroup(i) = 'W'; %Assigns ‘W’ label to ‘widefiles’
end
group = [cgroup; wgroup]; %Concatennates ‘cgroup’ and ‘wgroup’, stores in ‘group’
indices = crossvalind('kfold', group, K); %Gets indices for K fold cross validation
cp = classperf(group);
for i=1:K %Performs K fold cross validation using Nearest Neighbor classifier
    test = (indices==i);
    train =~ test;
    class = knnclassify(tmtr(test,:), tmtr(train,:), group(train,:)); 
    cperf = classperf(cp, class, test);
end
ConfusionMatrix = get(cperf, 'DiagnosticTable'); %Gets ‘ConfusionMatrix’ values
Error = get(cperf, 'ErrorRate'); %Gets ‘Error’
