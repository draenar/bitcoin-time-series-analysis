%load dataset
%load data
dataset = readtable('BlockChain_Train_csv_cleaned.csv');
%load half of data
data1 = dataset(1:size(dataset,1)/2,:);

%plot every column
for i = 2:39
    figure(i);
    clf;
    column = data1{:,i};
    days = data1{:,1};
    plot(days, column);
    xlabel('date');
    ylabel('value');
    title(data1.Properties.VariableNames(i),'Interpreter','none');
end