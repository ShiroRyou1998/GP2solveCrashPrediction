% this program is aimed to train the BP-NN model for accident prediction
% input:accident elements originData
% output:test and train set

clc
clear
%% 获得测试数据集
load('originDataSL.mat')

[numSample,numVar]=size(originData);
selectRate=0.1;
zeroNum=find(originData(:,end)==0,1); %事故为0的位置
onesNum=zeroNum-1; %事故为1的个数
selectNum=floor(selectRate*onesNum); %1的挑选数
selectIndexOne=randperm(onesNum,selectNum);
testOne=originData(selectIndexOne,:); %1的测试集

zeroNum=numSample-onesNum; %0的个数
selectNum=floor(selectRate*zeroNum); %0的挑选数
selectIndexZero=randperm(zeroNum,selectNum);
selectIndexZero=selectIndexZero+onesNum;%记录0的位置
testZero=originData(selectIndexZero,:); %0的测试集

originData([selectIndexOne,selectIndexZero],:)=[];%originData作为训练集
testData=[testOne;testZero];%总测试