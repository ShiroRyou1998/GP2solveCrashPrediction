%GEPfitness function
%fitness function : MSE

%aimed to get fitness of a pop in 3 steps
%step1:decode the chrom
%step2:caculate the variance
%step3:caculate the fitness

%input:pop, geneinfo, chorminfo, F T C and their features, sourcedata
%output:fitness list,variance list,maxfiness and its mathexpression

%warning:if soursedata or variable numbers changed, this function need too!

function [fitnessList,varList,maxfitness,maxMathexp,maxchrom,compareAcc]=...
    GEPfitness(pop,geneHead,geneTail,chromNum,Func,Tail,Fnary,Const,sourceData,protExp,calcuExp,indFit)
if nargin<12
    indFit='mae';
end
    
    %setting -- connection symbol
    connectSym='+';
    %sourcedata processing
    [dataNum,varNum]=size(sourceData);
    
    for i=1:varNum-1
        eval([char(64+i),'=sourceData','(:,',num2str(i),')',';']);
    end
    
    yP=sourceData(:,varNum);
    xi=sourceData(:,1:i);
    
    %body
    [popSize,chromSize]=size(pop);
    geneSize=chromSize/chromNum;
    
    %initialize
    fitnessList=zeros(popSize,1);
    varList=ones(popSize,1);
    maxfitness=-114514;
    maxMathexp='an error occurred!';
    maxchrom='shiro_ryou is best';
    
    %do the loop
    
    for i=1:popSize
        popTemp=pop(i,:);
        chromExp=[];
        
        %divide chrom into genes and transfer to mathexpression
        for j=1:chromNum
            subGene=popTemp((1+(j-1)*geneSize):j*geneSize);
            mathexp=GEPdecode(subGene,geneHead,geneTail,geneSize,Func,Fnary,Tail,Const);
            chromExp=[chromExp connectSym mathexp];
        end
        chromExp(1)=[];
        
        %calculate var and fitness
        voidChromExp=[chromExp protExp];%if var changers,ND changing
        regressFun=inline(vectorize(voidChromExp));
        
        %regressP=regressFun(A,B);%if var changers,input will change
        eval(['regressP=regressFun' calcuExp ';']);
        compareAcc=regressP;
        
        %if u want to use RMSE et al.
        switch indFit
            case 'mse'
                variance=sum( (regressP-yP).^2 );%MSE
            case 'mae'
                variance=sum( abs(regressP-yP) )/dataNum;%MAE
            case 'rmse'
                variance=sqrt(sum( (regressP-yP).^2 )/dataNum);%RMSE
            case 'csc'
                variance=GEPcostSense(yP,regressP)/dataNum; %cost-sensitive class
            case 'auc'
                [variance,~]=ROCplot(yP,regressP,0);%AUC
                variance=1./variance-1;
            otherwise
                variance=sum( abs(regressP-yP) )/dataNum;%MAE
        end
        fitness=1000/(1+variance);
        
        if isnan(fitness)
            % disp('NaN occurred');
            fitness=0;
            variance=Inf;
        end
        
        fitnessList(i)=fitness;
        varList(i)=variance;
        
        if fitness>maxfitness
            maxfitness=fitness;
            maxMathexp=chromExp;
            maxchrom=popTemp;
        end
        
        
    end
    
end