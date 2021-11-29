function [AUC,youdenIndex]=ROCplot(realSample,possSample,plotSwitch)
if nargin<3
    plotSwitch=1;% need plot
end
ROCfpr=[];
ROCtpr=[];
ROCTP=[];
ROCTN=[];
ROCFP=[];
ROCFN=[];
youdenDiff=-114514;
youdenIndex=114514;

allSample=[realSample;possSample];
for threshold=max([min(possSample),0]):10^-3:min([max(possSample),1])
    predSample=(possSample>threshold);
    [TPR,FPR,TP,TN,FP,FN]=TPcalcu(realSample,predSample);
    ROCfpr=[ROCfpr FPR];
    ROCtpr=[ROCtpr TPR];
    ROCTP=[ROCTP TP];
    ROCTN=[ROCTN TN];
    ROCFP=[ROCFP FP];
    ROCFN=[ROCFN FN];
    if TPR-FPR>youdenDiff
        youdenIndex=threshold;
        youdenDiff=TPR-FPR;
    end
end

AUC = polyarea([1 ROCfpr 0 1],[1 ROCtpr 0 0]);

if plotSwitch
    plot([1 ROCfpr 0],[1 ROCtpr 0],'b-*');
    grid on
    hold on
    plot([0,1],[0,1],'r-');
    
    ylabel({'TPR'});xlabel({'FPR'});title(['AUC=' num2str(AUC,3)]);
    %annotation('textbox',...
    %    [0.577 0.343 0.125 0.1],'String',{'AUC=' num2str(AUC)});
    
    predSample=(possSample>youdenIndex);
    [TPR,FPR,~,~,~,~]=TPcalcu(realSample,predSample);
    %annotation('textbox',...
    %    [0.577 0.18 0.225 0.15],'String',{['yd=' num2str(youdenIndex) ' FPR='...
    %    num2str(FPR) ' TPR=' num2str(TPR)]});
    %scatter1=scatter(FPR,TPR,30,200,'filled');
    %datatip(scatter1,'DataIndex',1);
end

    function [TPR,FPR,TP,TN,FP,FN]=TPcalcu(realSample,predSample)
        
        trueSheet=realSample+predSample;
        TP=sum(trueSheet==2);
        TN=sum(trueSheet==0);
        
        falseSheet=realSample-predSample;
        FP=sum(falseSheet==-1);
        FN=sum(falseSheet==1);
        
        condPositive=TP+FN;
        condNegative=FP+TN;
        
        FPR=FP/condNegative;
        TPR=TP/condPositive;
        
    end
end