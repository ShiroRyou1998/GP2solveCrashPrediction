%GEPnewpop function
%aimed to creat new population

%input:function set and terminate set,gene info,population info
%output:newpop

function newpop=GEPnewpop(Func,Tail,Const,geneHead,geneTail,popSize,chromNum)

newpop=[];

FT=[Func Tail];
Tlength=length(Tail);
FTlength=length(FT);
Clength=length(Const);

for i=1:popSize
    chromBody=[];
    
    for j=1:chromNum
        geneBody=[];
        for k=1:geneHead
            geneBody=[geneBody FT(randperm(FTlength,1))];
        end
        
        for k=1:geneTail
            geneBody=[geneBody Tail(randperm(Tlength,1))];
        end
        
        for k=1:geneTail
            geneBody=[geneBody int2str(randperm(Clength,1))];
        end
        
        chromBody=[chromBody geneBody];
    end
    
    newpop=[newpop;chromBody];
end

end