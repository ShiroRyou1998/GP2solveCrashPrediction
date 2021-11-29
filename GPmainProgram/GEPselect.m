%GEPselect function
%select newpop via russian_roulette

%input:fitnesslist,prepop
%output:newpop

function newpop=GEPselect(fitnessList,prepop)

[popsize,~]=size(fitnessList);

% square method
fitnessList=fitnessList.^3;
% an effient method to enlarge the diff between each chrom

newpop=[];

probability=fitnessList/(sum(fitnessList));

probability=cumsum(probability);

for i=1:popsize
    needle=rand;
    for j=1:popsize
        if needle<probability(j)
            newpop=[newpop;prepop(j,:)];
            break;
        end
    end
end


end