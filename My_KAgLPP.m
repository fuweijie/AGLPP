function [Z,eigvector,eigvalue]=My_KAgLPP(data,anchor,d,sigma)
      options = [];
      options.Metric = 'Euclidean';
      options.NeighborMode = 'KNN';
      options.k = 5;
      options.WeightMode = 'HeatKernel';
      options.PCARatio = 1;
      options.ReducedDim=d;
      options.t=7;
[Z, D, W] = AnchorGraph(data', anchor', 3, 0, 10,sigma);
[n,m]=size(Z);
% data=eye(m,m);
% 
% if (~exist('options','var'))
%    options = [];
% end
% 
% [nSmp,nFea] = size(data);
% if size(W,1) ~= nSmp
%     error('W and data mismatch!');
% end
% 
if isfield(options,'keepMean') && options.keepMean
    ;
else
    if issparse(data)
        data = full(data);
    end
    sampleMean = mean(data);nSmp=length(data);
    data = (data - repmat(sampleMean,nSmp,1));
end
%==========================

if (~exist('options','var'))
   options = [];
end

if isfield(options,'Kernel') && options.Kernel
    K = data;
    clear data;
else
    K = double(constructKernel(anchor,[],options,sigma));
end

options.Kernel = 1;

[eigvector, eigvalue]  = My_KGE_(W, D, options, K,sigma);


eigIdx = find(eigvalue < 1e-3);
eigvalue (eigIdx) = [];
eigvector(:,eigIdx) = [];



