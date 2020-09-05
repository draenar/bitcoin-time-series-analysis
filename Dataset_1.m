%load data
dataset = readtable('BlockChain_Train_csv_cleaned.csv');
%load half of data
data1 = dataset(1:size(dataset,1)/2,:);


%plot every column
% for i = 2:39
%     figure(i);
%     clf;
%     column = data1{:,i};
%     days = data1{:,1};
%     plot(days, column);
%     xlabel('date');
%     ylabel('value');
%     title(data1.Properties.VariableNames(i),'Interpreter','none');
% end

%plot log returns
% for i = 2:39
%     figure(i);
%     clf;
%     column = log(data1{2:end,i})-log(data1{1:end-1,i});
%     days = data1{2:end,1};
%     plot(days, column);
%     xlabel('date');
%     ylabel('value');
%     title(data1.Properties.VariableNames(i),'Interpreter','none');
% end

%plot autocorrelations
alpha = 0.05;
% zalpha = norminv(1-alpha/2);
zalpha = 1.96;
maxtau = 20;
[n,m]=size(data1);
column = data1{:,8};
ACR = autocorrelation(column,maxtau);
figure(7)
clf
plot(ACR(:,1),ACR(:,2),'.-')
hold on
plot([0 maxtau],(zalpha/sqrt(n))*[1 1],'c--')
plot([0 maxtau],-(zalpha/sqrt(n))*[1 1],'c--')
xlabel('lag \tau')
ylabel('r_Y(\tau)')
title(sprintf('autocorrelation of stock %s, residual from AR(%d)',name1,p))
