function[cperf, ConfusionMatrix, Error]=ProteinClass(proteinA, proteinB, K)
%[cperf, ConfusionMatrix]=ProteinClass(proteinA, proteinB, K)
proteinA_file = textread(proteinA, '%s', 'delimiter', '\n');
numproteinA = length(proteinA_file);
proteinAmtr = zeros(numproteinA, 1024);
ctr=1;
for i=1:numproteinA
  img = double(imread(proteinA_file{i}));
  FT = fft2(img);
  FT_shifted = fftshift(FT);
  FTsumproteinA = sum(FT_shifted);
  proteinAmtr(i,:) = FTsumproteinA;
  ctr=ctr+1;
end
proteinB_file = textread(proteinBfiles, '%s', 'delimiter', '\n');
numproteinBf = length(proteinB_file);
proteinBmtr = zeros(numproteinB, 1024);
for i=1:numproteinB
  img = double(imread(proteinB_file{i}));
  FT = fft2(img);
  FT_shifted = fftshift(FT);
  FTsumproteinB = sum(FT_shifted);
  proteinBmtr(i,:) = FTsumproteinB;
  ctr=ctr+1;
end
tmtr = [proteinAmtr; proteinBmtr];
Agroup = zeros(numproteinA,1);
for i=1:numproteinA
    Agroup(i) = 'A';
end
Bgroup = zeros (numproteinB, 1);
for i=1:numproteinB
    Bgroup(i) = 'B';
end
group = [Agroup; Bgroup];
indices = crossvalind('kfold', group, K);
cp = classperf(group);
for i=1:K
    test = (indices==i);
    train =~ test;
    class = knnclassify(tmtr(test,:), tmtr(train,:), group(train,:));
    cperf = classperf(cp, class, test);
end
class = char(class);
cp.ErrorRate
proteinAusionMatrix = get(cperf, 'DiagnosticTable');
Error = get(cperf, 'ErrorRate');

