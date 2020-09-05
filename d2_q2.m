%question 2 of project
%parameters
alpha = 0.01;
% zalpha = norminv(1-alpha/2);
zalpha = 1.96;
maxtau = 20;
maxtau2 = 10;
rthresh = 0.25;
tau = 0;
npart = 3;
rng(1);
maxwordlength = 15;

%load dataset
%load data
dataset = readtable('BlockChain_Train_csv_cleaned.csv');
%load half of data
data1 = dataset(614:end,:);

%dataset size
[n,m]=size(data1);
K = m -1;
%plot all autocorrelations
% for i = 2:39
%     figure(i);
%     clf;
%     column = data1{:,i};
%     ACR = autocorrelation(column,maxtau);
%     figure(i);
%     clf;
%     plot(ACR(:,1),ACR(:,2),'.-')
%     hold on
%     plot([0 maxtau],(zalpha/sqrt(n))*[1 1],'c--')
%     plot([0 maxtau],-(zalpha/sqrt(n))*[1 1],'c--')
%     xlabel('lag \tau')
%     ylabel('r_Y(\tau)')
%     title(['autocorrelation of ',data1.Properties.VariableNames(i)],'Interpreter','none');
% end
  
%remove autocorrelation with log returns 
% for i = 2:39
%     figure(100+i);
%     clf;
%     column = log(data1{2:end,i})-log(data1{1:end-1,i});
%     ACR = autocorrelation(column,maxtau);
%     figure(i);
%     clf;
%     plot(ACR(:,1),ACR(:,2),'.-')
%     hold on
%     plot([0 maxtau],(zalpha/sqrt(n))*[1 1],'c--')
%     plot([0 maxtau],-(zalpha/sqrt(n))*[1 1],'c--')
%     xlabel('lag \tau')
%     ylabel('r_Y(\tau)')
%     title(['autocorrelation of log returns for ',data1.Properties.VariableNames(i)],'Interpreter','none');
% end

%calculate correlations with close price of bitcoin
% for i = 3:39
%     closePrice = log(data1{2:end,2})-log(data1{1:end-1,2});
%     column = log(data1{2:end,i})-log(data1{1:end-1,i});
%     
%     cceyV = mycrosscorr(closePrice,column,maxtau2);
%     figure(200+i)
%     clf
%     plot([-maxtau2:maxtau2]',cceyV,'.-')
%     hold on
%     plot([-maxtau2 maxtau2],(zalpha/sqrt(n))*[1 1],'c--')
%     plot([-maxtau2 maxtau2],-(zalpha/sqrt(n))*[1 1],'c--')
%     xlabel('lag \tau')
%     ylabel('r_{XY}(\tau)')
%     title(['Cross Correlation of ',data1.Properties.VariableNames(2),'and ',data1.Properties.VariableNames(i),'log returns'],'Interpreter','none');
%     
% 
% end

%calculate correlation matrix 

xmatrix = data1{2:end,2:end};
xmatrixminusone = data1{1:end-1,2:end};
%log returns matrix
xM = log(xmatrix) - log(xmatrixminusone);

nameM = data1.Properties.VariableNames(2:end);

% For each pair compute the correlation matrix
ccM = NaN*ones(K,K);
p1M = zeros(K,K);
if tau==0
    % The correlation matrix is symmetric
    [ccM,p1M] = corrcoef(xM);
    p1M(1:K+1:K*K) = 0;
else
    % The correlation matrix is not symmetric
    for ik=1:K-1
        for jk=ik+1:K
            [tmpM,ptmpM] = corrcoef(xM(1:end-tau,ik),xM(1+tau:end,jk));
            ccM(ik,jk) = tmpM(1,2);
            p1M(ik,jk) = ptmpM(1,2);
            [tmpM,ptmpM] = corrcoef(xM(1:end-tau,jk),xM(1+tau:end,ik));
            ccM(jk,ik) = tmpM(1,2);
            p1M(jk,ik) = ptmpM(1,2);
        end
    end
end    

tit1txt = sprintf('R_{XY}(%d)',tau);
h1 = plotnetworktitle(ccM,[],nameM,tit1txt,1);

adj1M = p1M < alpha;
tit2txt = sprintf('Adjacency p(R_{XY}(%d)) < %1.2f',tau,alpha);
h2 = plotnetworktitle(adj1M,[0 1],nameM,tit2txt,2);

adjfdr1M = adjFDRmatrix(p1M,alpha,2);
tit3txt = sprintf('FDR (a=%1.3f) R_{XY}(%d)',alpha,tau);
h3 = plotnetworktitle(adjfdr1M,[0 1],nameM,tit3txt,3);

rthreshM = abs(ccM) > rthresh;
tit4txt = sprintf('Adjacency R_{XY}(%d) > %1.2f',tau,rthresh);
h4 = plotnetworktitle(rthreshM,[0 1],nameM,tit4txt,4);