% Function to help choose initial angles
function bin_class = equivalent_class(bins)
bin_class = bins;
for i = 1 : size(bins,1)
    bin = bins(i,:);
    indices = find(bin==min(bin));
    if length(indices) ~= length(bin)
        count = 0;
        while length(indices) ~= 1 && count < length(bin)
            indices = mod(indices+1, length(bin));
            indices(indices==0) = length(bin);
            indices = indices(bin(indices)==min(bin(indices)));
            count = count + 1;
        end
        if length(indices) ~= 1
            indices = indices(1);
        end
        bin_class(i,:) = circshift(bin, length(bin)-indices+count+1);
    end
end
end