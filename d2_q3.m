%question 3 of project
%parameters
alpha = 0.05;
P = 5; % The order of the VAR model used for the computation of the 
        % Granger causality index (GCI) 
GCIthresh = 0.10;
CGCIthresh = 0.10;
rng(1);

%load dataset
%load data
dataset = readtable('BlockChain_Train_csv_cleaned.csv');
%load half of data
data1 = dataset(614:end,:);

%dataset size
[n,m]=size(data1);


%log returns matrix
xmatrix = data1{2:end,2:end};
xmatrixminusone = data1{1:end-1,2:end};
xM = log(xmatrix) - log(xmatrixminusone);

labelM = data1.Properties.VariableNames(2:end);

%gci
[GCIM,pGCIM] = GCI(xM,P,1);

%% Plot the GCI-causality network
% The network of weighted connections given by GCI_{X->Y}(P)
tit1txt = sprintf('GCI_{X->Y}(%d)',P);
plotnetworktitle(GCIM,[],labelM,tit1txt,1);

adj1M = pGCIM < alpha;
tit2txt = sprintf('Adjacency p(GCI_{X->Y}(%d)) < %1.2f',P,alpha);
plotnetworktitle(adj1M,[0 1],labelM,tit2txt,2);

adjfdr1M = adjFDRmatrix(pGCIM,alpha,2);
tit3txt = sprintf('FDR (a=%1.3f) GCI_{X->Y}(%d)',alpha,P);
plotnetworktitle(adjfdr1M,[0 1],labelM,tit3txt,3);

GCIthreshM = GCIM > GCIthresh;
tit4txt = sprintf('Adjacency GCI_{X->Y}(%d) > %1.2f',P,GCIthresh);
plotnetworktitle(GCIthreshM,[0 1],labelM,tit4txt,4);

%cgci
[CGCIM,pCGCIM] = CGCI(xM,P,1);

tit1txt = sprintf('CGCI_{X->Y}(%d)',P);
plotnetworktitle(CGCIM,[],labelM,tit1txt,5);

adj1M = pCGCIM < alpha;
tit2txt = sprintf('Adjacency p(CGCI_{X->Y}(%d)) < %1.2f',P,alpha);
plotnetworktitle(adj1M,[0 1],labelM,tit2txt,6);

adjfdr1M = adjFDRmatrix(pCGCIM,alpha,2);
tit3txt = sprintf('FDR (a=%1.3f) CGCI_{X->Y}(%d)',alpha,P);
plotnetworktitle(adjfdr1M,[0 1],labelM,tit3txt,7);


CGCIthreshM = CGCIM > CGCIthresh;
tit4txt = sprintf('Adjacency CGCI_{X->Y}(%d) > %1.2f',P,CGCIthresh);
plotnetworktitle(CGCIthreshM,[0 1],labelM,tit4txt,8);
