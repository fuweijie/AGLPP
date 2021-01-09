function [Z,eigvector,eigvalue]=My_AgLPP(data,anchor,d,sigma)
      options = [];
      options.Metric = 'Euclidean';
      options.NeighborMode = 'KNN';
      options.k = 5;
      options.WeightMode = 'HeatKernel';
      options.PCARatio = 0.99;
      options.ReducedDim=d;
      
[Z, D, W] = AnchorGraph(data', anchor', 3, 0, 10,sigma);
data=anchor;

if (~exist('options','var'))
   options = [];
end

[nSmp,nFea] = size(data);
if size(W,1) ~= nSmp
    error('W and data mismatch!');
end

if isfield(options,'keepMean') && options.keepMean
    ;
else
    if issparse(data)
        data = full(data);
    end
    sampleMean = mean(data);
    data = (data - repmat(sampleMean,nSmp,1));
end
%==========================

[eigvector, eigvalue] = LGE(W, D, options, data);


eigIdx = find(eigvalue < 1e-3);
eigvalue (eigIdx) = [];
eigvector(:,eigIdx) = [];



