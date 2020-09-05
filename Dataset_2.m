%load data
dataset = readtable('BlockChain_Train_csv_cleaned.csv');
%load half of data
data2 = dataset(614:end,:);


%load data
dataset = readtable('BlockChain_Train_csv_cleaned.csv');
%load half of data
data1 = dataset(1:size(dataset,1)/2,:);

%plot bitcoin value
figure(1);
clf;
bitcoin = data1.;
days = data1.Date;
plot(days, bitcoin);
xlabel('date');
ylabel('value');
title(data1.Properties.VariableNames(2));

%plot every column
data1{:,2};