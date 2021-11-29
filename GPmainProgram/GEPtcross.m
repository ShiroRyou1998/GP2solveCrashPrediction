%GEPtcross function
%aimed to conduct 3 cross oprerations--single double and allelic gene

%input:pop, chrominfo, probility of 3 methods
%output:newpop

function newpop=GEPtcross(pop,chromNum,pcrosss,pcrossd,pcrossg)

[popsize,chromSize]=size(pop);
geneSize=chromSize/chromNum;

%single
selectOrder=randperm(popsize);

for i=1:2:popsize
    needle=rand;
    if needle<pcrosss
        chromPar1=pop(selectOrder(i),:);
        chromPar2=pop(selectOrder(i+1),:);
        
        pointCut=randperm(chromSize-1,1);
        chromCut1=chromPar1((pointCut+1):chromSize);
        chromCut2=chromPar2((pointCut+1):chromSize);
        
        pop(selectOrder(i),:)=[chromPar1(1:pointCut) chromCut2];
        pop(selectOrder(i+1),:)=[chromPar2(1:pointCut) chromCut1];
    end

end

%double
selectOrder=randperm(popsize);

for i=1:2:popsize
    needle=rand;
    if needle<pcrossd
        chromPar1=pop(selectOrder(i),:);
        chromPar2=pop(selectOrder(i+1),:);
        
        pointCut=randperm(chromSize-1,2);
        pointCut1=min(pointCut);
        pointCut2=max(pointCut);
        
        chromCut1=chromPar1((pointCut1+1):pointCut2);
        chromCut2=chromPar2((pointCut1+1):pointCut2);
        
        pop(selectOrder(i),:)=[chromPar1(1:pointCut1) chromCut2 ...
            chromPar1((pointCut2+1):chromSize)];
        pop(selectOrder(i+1),:)=[chromPar2(1:pointCut1) chromCut1 ...
            chromPar2((pointCut2+1):chromSize)];
    end

end

%allelic gene
selectOrder=randperm(popsize);

for i=1:2:popsize
    needle=rand;
    if needle<pcrossg
        chromPar1=pop(selectOrder(i),:);
        chromPar2=pop(selectOrder(i+1),:);
        
        needle=randperm(chromNum,1);%first&last is single, others is double
        
        if needle==1 || needle==chromNum
            %single, pointCut the first pos of a gene
            pointCut=geneSize*(needle-1);
            chromCut1=chromPar1((pointCut+1):chromSize);
            chromCut2=chromPar2((pointCut+1):chromSize);
        
            pop(selectOrder(i),:)=[chromPar1(1:pointCut) chromCut2];
            pop(selectOrder(i+1),:)=[chromPar2(1:pointCut) chromCut1];
        else
            %double
        
            pointCut1=geneSize*(needle-1);
            pointCut2=geneSize*needle-1;
        
            chromCut1=chromPar1((pointCut1+1):pointCut2);
            chromCut2=chromPar2((pointCut1+1):pointCut2);
        
            pop(selectOrder(i),:)=[chromPar1(1:pointCut1) chromCut2 ...
            chromPar1((pointCut2+1):chromSize)];
            pop(selectOrder(i+1),:)=[chromPar2(1:pointCut1) chromCut1 ...
            chromPar2((pointCut2+1):chromSize)];

        end
    end

end

newpop=pop;
end