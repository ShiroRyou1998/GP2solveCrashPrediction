function PDPplot(originData,Xnum,Ynum,Price,calcuExp)
load('name')
[~,varNum]=size(originData);
for i=1:varNum-1
    eval([char(64+i),'=originData','(:,',num2str(i),')',';']);
end
X=originData(:,Xnum);
Y=originData(:,Ynum);

X=linspace(min(X),max(X),20);
Y=linspace(min(Y),max(Y),20);
[X,Y]=meshgrid(X,Y);

for j=1:varNum-1
    if j==Xnum
        eval([char(64+j),'=X',';']);
    else
        if j==Ynum
            eval([char(64+j),'=Y',';']);
        else
            eval([char(64+j),'=ones(size(X)).*mean(',char(64+j),');']);
        end
    end
end

eval(['Z=Price' calcuExp ';']);
surf(X,Y,Z);
xlabel(name(Xnum))
ylabel(name(Ynum))
zlabel('risk')
