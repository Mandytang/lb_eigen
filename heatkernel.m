function h = heatkernel(t,vfix,evecs,evals,N)
%
% computes heat kernel from all points to a fixed point (vfix)
% for a given time t (using only the first N smallest eigenvalues
% and eigenvectors)
%
% Output:
% h   matrix rows: all vertices, cols: times in t

%make sure evals is vector:
if size(evals,1) == size(evals,2)
    evals = diag(evals);
end

%make sure first are the smallest:
if evals(end) < evals(1)
    evals = evals(end:-1:1);
    evecs = evecs(:,end:-1:1);
end

evals = evals(1:N);
evecs = evecs(:,1:N);


%make evals col vec
if size(evals,2) ~= 1
    evals = evals';
end
%make t row vec
if (size(t,1) ~= 1)
    t = t'
end


h = evecs * ( exp(-evals * t) .* repmat(evecs(vfix,:)',1,length(t))  );
