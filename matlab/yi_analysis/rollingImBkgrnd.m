function [bkgrnd]=rollingImBkgrnd(im,neighborhood)



%default neighborhood

if ~exist('neighborhood')

    neighborhood=15;

elseif isempty(neighborhood)

    neighborhood=15;

end


se = strel('disk',neighborhood);


bkgrnd = imopen(im,se);
