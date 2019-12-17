function out=rowshufflev1(template)
% shuffles bins of a matrix row-wise

[m,n] = size(template) ;
b=zeros(size(template));
for rowplate=1:m
idx = randperm(n) ;
b(rowplate,:) = template(rowplate,:);
b(rowplate,idx) = template(rowplate,:) ;
end
out=b;