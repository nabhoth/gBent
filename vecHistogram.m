function [ bars ] = vecHistogram( functions )
%vecHistogram given is a matrix of spectra, this scrit generates a histograms
%as would hist(functions, unique8functions)) if hist could be doing that
[rows,cols] = size(functions);
urows = unique(functions, 'rows');
[ur, uc] = size(urows);
bars = zeros(1,ur);
for i=1:ur
    l = ismember(functions, urows(i,:),'rows');
    bars(1,i) = sum(l);
end

end

