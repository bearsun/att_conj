function hlslope= slope( hl )
%SLOPE Summary of this function goes here
%   Detailed explanation goes here

hlslope = cellfun(@(a) polyfit(1:10,a,1), mat2cell(hl,ones(1,9)),'uni',0);
hlslope = cell2mat(hlslope);
hlslope = hlslope(:,1);

end

