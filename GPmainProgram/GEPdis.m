%GEPdis function
%used to conduct IS&RIS--2 IS operations in total

%input:pop,geneinfo,chrominfo,length that IS need, F
%      possibility of IS/RIS
%output:newpop

function newpop=GEPdis(prepop,geneHead,geneTail,chromNum,isLength,F,pis,pris)

[popSize,chromSize]=size(prepop);
geneSize=chromSize/chromNum;
bodySize=geneHead+geneTail;%head & tail
newpop=[];

for i=1:popSize
    popTemp=prepop(i,:);%select a chrom
    chromTemp=[];%create an empty chrom to fill
    
    for j=1:chromNum
        subGene=popTemp((1+(j-1)*geneSize):j*geneSize);%select piece of gene
        DCfield=subGene((bodySize+1):geneSize);
        subGene=subGene(1:bodySize);
        
        needle=rand;
        if needle<pis %IS
            needle=randperm(bodySize,1);
            cutNum=randperm(min([isLength,bodySize-needle+1]),1);
            cutIs=subGene(needle:(needle+cutNum-1));%get certain string
            
            needle=randperm(geneHead-1,1)+1;%insert position
            newHead=[subGene(1:needle-1) cutIs subGene(needle:geneHead)];
            newHead=newHead(1:geneHead);
            subGene=[newHead subGene((geneHead+1):bodySize)];
            
        end

        needle=rand;
        if needle<pris %RIS
            needle=randperm(geneHead-1,1)+1;
            
            for k=needle:geneHead
                if strfind(F,subGene(k))
                    %k is the head of RIS
                    cutNum=randperm(min([isLength,bodySize-k+1]),1);
                    cutRis=subGene(k:(k+cutNum-1));%get certain string
                    newHead=[cutRis subGene];
                    newHead=newHead(1:geneHead);
                    subGene=[newHead subGene((geneHead+1):bodySize)];
                    
                    break;%get out of the loop
                end
            end
        end
        
        chromTemp=[chromTemp subGene DCfield];%fill
                
    end
    
    newpop=[newpop;chromTemp];
    
end


end