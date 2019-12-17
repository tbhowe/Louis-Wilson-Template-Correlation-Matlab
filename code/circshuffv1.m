function out=circshuffv1(template)
% shuffles bins of a matrix row-wise
warning('off')
[m,n] = size(template) ;
range=floor(n/2);
b=zeros(size(template));
r = randi([ -range range ],1,m);
for rowplate=1:m
b(rowplate,:) = circshift(template(rowplate,:),[1 r(rowplate)]);
end
out=b;