function out=rowshuffler1(template)
% shuffles bins of a matrix row-wise

[m,~] = size(template) ;
b=zeros(size(template));
ordering = randperm(m);
B = template(ordering, :);

out=B;