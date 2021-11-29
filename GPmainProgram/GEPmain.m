%GEP programe main v1.0.0

%this program is aimed to get regression function via GEP algorithm
%attention:this program apply basic GEP and HGA method
%main.m include 3 bodies, datainput model, processing model, conduct model
%decode, select, mutate, IS, cross, 5 parts covered

%1.body established, some bug corrected
%2.select,mutate,IS,RIS,3 cross have added
%3.plot function added
%4.HGA and elite gene library strategy added

%created by sean98ma in 2020/4/14
%last edited by sean98ma in 2021/10/28

clc
clear
close all
tic

%set debug config
dbstop if error

%% INPUT & INITIALIZE
%only paramaters in this block can be modified
%input data
load('GEPvar12.mat')
%last variable of Data is y, others are xi
%name data as 'originData', #row is #samples, #cline is #variables
%if originData is 368*4, represent 368samples of 4 varibales (3x+1y)

%pop setting
popSize=50;%even#,50-100 is recomand
iterationMax=100;%100-500
indexFinish=0.0001;
countFinish=0;

%gene setting
geneHead=9;
chromNum=4;%each chrom used '+' to connect
Func=['/' '+' '-' '*' 'l' 'q' 's'];
%Const=[1.0964  0.4390   -2.1498    0.5188   -1.3976   -0.4095   -0.8592    1.5455   -0.9755];
Const=10*(rand(1,9)-0.5);
indFit='mae';% mae* mse rmse csc auc

pmutate=0.144;pis=0.1;pris=0.1;
pcrosss=0.9;pcrossd=0.4;pcrossg=0.2;
isLength=3;
selectZeroRate=0.35; %how many 0-samples are selected to train with 1-samples(0-1)

%%%%%%%%%%%%%%%%%% ALGORITHM WORK ZONE %%%%%%%%%%%%%%%%%%%%%%%%%%%
%% process paramaters
Fnary=ones(1,length(Func));
Fnary(1:4)=2;
[dataNum,varNum]=size(originData);
TailNum=varNum-1;%variables of x, e.g. y=sin(x1)+cos(x2),tailnum=2
Tail=[];protExp=['+0'];calcuExp=['('];
for i=1:TailNum
    Tail=[Tail char(64+i)];
    protExp=[protExp '*' char(64+i) ];
    calcuExp=[calcuExp char(64+i) ','];
end
Tail=[Tail '?'];calcuExp(end)=[')'];
geneTail=(max(Fnary)-1)*geneHead+1;%genecost = genetail

try
    %load('eliteLib.mat')
catch
    eliteLib=[];
end
eliteLib=[];

% select sample randomly
onesNum=find(originData(:,end)==0,1)-1;
zerosNum=dataNum-onesNum;
selectZeroIndex=onesNum+randperm(zerosNum,floor(zerosNum*selectZeroRate));
sourceDataZero=originData(selectZeroIndex,:);
sourceDataOne=originData(1:onesNum,:);
sourceData=[sourceDataOne;sourceDataZero];


%% GEP body
%initialize parameter
maxfitness=0;
bestfitness=0;
bestchrom='shiro_ryou';
bestindividual='shiro_ryou is the best';
compareAcc=1919810;

%initialize plot function
maxfitPlot=zeros(1,iterationMax);
maxvarPlot=zeros(1,iterationMax);
popfitPlot=zeros(1,iterationMax);

%creat newpop
newpop=GEPnewpop(Func,Tail,Const,geneHead,geneTail,popSize,chromNum);

%do the loop

for i=1:iterationMax
    
    %caculate fitness
    [fitnessList,varList,maxfitness,maxMathexp,maxchrom,compareAcc]=...
        GEPfitness(newpop,geneHead,geneTail,chromNum,Func,Tail,Fnary,Const,sourceData,protExp,calcuExp,indFit);
    
    %continue?
    if bestfitness==900
        popfitPlot(i)=maxfitness;
        maxfitPlot(i)=bestfitness;
        maxvarPlot(i)=bestvariance;
        bestcompare=compareAcc;
        disp('GEP find the bestchorm(1000), so end loop ahead of schedule')
        break;
    end
    
    if abs(bestfitness-maxfitness)<indexFinish
        countFinish=countFinish+1;
        if countFinish==30
            %plot info
            popfitPlot(i)=maxfitness;
            maxfitPlot(i)=bestfitness;
            maxvarPlot(i)=bestvariance;
            bestcompare=compareAcc;
            disp('GEP reach converge, so end loop ahead of schedule')
            break;
        end
    else
        countFinish=0;
    end
    
    popfitPlot(i)=maxfitness;%plot info
    
    %update the bestindividual
    if maxfitness>bestfitness
        bestfitness=maxfitness;
        bestvariance=1000/bestfitness-1;
        bestindividual=maxMathexp;
        bestchrom=maxchrom;
        bestcompare=compareAcc;
    end
    
    %update iretation plot
    maxfitPlot(i)=bestfitness;
    maxvarPlot(i)=bestvariance;
    
    %dynamic plot test
    cla;
    hold on;
    plot(1:i,popfitPlot(1:i),'b');
    plot(1:i,bestfitness*ones(1,i),'r');
    if i==1
        title('Plot 1: maxfitness of each generation')
        xlabel('genaration')
        ylabel('fitness')
        grid on
    end
    %pause(0.01)
    drawnow
    
    %select individuals that join GEPopration
    newpop=GEPselect(fitnessList,newpop);
    
    %mutate
    newpop=GEPmutate(newpop,geneHead,geneTail,chromNum,Func,Tail,Const,pmutate);
    
    %IS & RIS
    newpop=GEPdis(newpop,geneHead,geneTail,chromNum,isLength,Func,pis,pris);
    
    %cross (single double & gene)
    newpop=GEPtcross(newpop,chromNum,pcrosss,pcrossd,pcrossg);
    
    %HGA strategy
    if ~mod(i,20)
        eliteLib=[eliteLib;bestchrom];
    end
    
    if ~mod(i,40)
        newpop=GEPelite(newpop,eliteLib);
    end

end

%% result and test plot info
Price=inline(vectorize([bestindividual protExp]));
[dataNum,varNum]=size(testData);
% test data
for i=1:varNum-1
    eval([char(64+i),'=testData','(:,',num2str(i),')',';']);
end
eval(['yTest=Price' calcuExp ';']);
figure(2)
subplot(2,2,1)
scatter(1:length(yTest),yTest,15,'r','filled')
hold on
scatter(1:length(yTest),testData(:,end),18,'b')
grid on
title('Test data prediction')
legend('risk','event')
subplot(2,2,2)
[AUCtest,~]=ROCplot(testData(:,end),yTest);
% train data

for i=1:varNum-1
    eval([char(64+i),'=sourceData','(:,',num2str(i),')',';']);
end
eval(['yTest=Price' calcuExp ';']);
subplot(2,2,3)
scatter(1:length(yTest),yTest,15,'r','filled')
hold on
scatter(1:length(yTest),sourceData(:,end),18,'b')
grid on
title('Train data prediction')
subplot(2,2,4)
[AUCtrain,~]=ROCplot(sourceData(:,end),yTest);
suptitle('Prediction & ROC curve of test & train data')

%save elite gene library
save('eliteLib.mat','eliteLib')


toc