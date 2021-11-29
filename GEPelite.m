% GEPelite function
% only used in HGA model
% before using this function, plz confirm elite gene lib has been set

function newpop=GEPelite(pop,eliteLib)
eliteRate=0.1;

[popSize,~]=size(pop);
[libSize,~]=size(eliteLib);
addNum=ceil(eliteRate*popSize)*2;

for i=1:addNum
    selectIndex=randperm(libSize,1);
    pop=[pop;eliteLib(selectIndex,:)];
end

newpop=pop;
end