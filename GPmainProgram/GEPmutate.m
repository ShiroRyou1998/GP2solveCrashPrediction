%GEPmutate function
%aimed to conduct mutation
%every gene loci will be detected

%input:pop, chrom info, gene info, FTC & their features, mutate possibility
%output:newpop

function newpop=GEPmutate(prepop,geneHead,geneTail,chromNum,Func,Tail,Const,pmutate)

[popSize,chromSize]=size(prepop);
geneSize=chromSize/chromNum;
FT=[Func Tail];
Tlength=length(Tail);
FTlength=length(FT);
Clength=length(Const);
newpop=[];

for i=1:popSize
    popTemp=prepop(i,:);%select a chrom
    chromTemp=[];%create an empty chrom to fill
    
    for j=1:chromNum
        subGene=popTemp((1+(j-1)*geneSize):j*geneSize);%select piece of gene
        
        %do mutation for every gene loci
        for k=1:geneSize
            needle=rand;
            if needle<pmutate
                %mutate opration
                if k<=geneHead
                    %head change in F&T
                    subGene(k)=FT(randperm(FTlength,1));
                    
                else if k<=(geneHead+geneTail)
                        %tail change in T
                        subGene(k)=Tail(randperm(Tlength,1));
                    else
                        %DCfield change in C
                        subGene(k)=num2str(randperm(Clength,1));
                    end
                end
                
                
                
            end
        end
        
        chromTemp=[chromTemp subGene];%fill
    end
    
    newpop=[newpop;chromTemp];
end


end