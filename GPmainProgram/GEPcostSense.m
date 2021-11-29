% this func is aimed to calcu the cost of FP FN fitness
% you can change the iteration step itStep to get more or less dencity
% you can change the cost matrix to get more suitable cost

function costFit=GEPcostSense(realSample,possSample)
itStep=10^-2;
costMat=[1;4];% FP FN, FN>>FP

costFit=114514;
for threshold=max([min(possSample),0]):itStep:min([max(possSample),1])
    predSample=(possSample>threshold);
    [~,~,~,~,FP,FN]=TPcalcu(realSample,predSample);
    costTemp=[FP,FN]*costMat;
    if costTemp<costFit
        costFit=costTemp;
    end
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