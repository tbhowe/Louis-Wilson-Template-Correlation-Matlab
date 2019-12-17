function out=columnshufflev1(template)
% shuffles columns of a matrix

[~,n] = size(template) ;
b=zeros(size(template));
idx = randperm(n) ;

for colplate=1:n
b(:,idx(colplate)) = template(:,colplate) ;
end
out=b;