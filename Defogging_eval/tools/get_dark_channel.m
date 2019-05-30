% function dark_channel = get_dark_channel(image, win_size)
% 
% [m, n, ~] = size(image);
% 
% pad_size = floor(win_size/2);
% 
% padded_image = padarray(image, [pad_size pad_size], Inf);
% 
% dark_channel = zeros(m, n); 
% 
% for j = 1 : m
%     for i = 1 : n
%         patch = padded_image(j : j + (win_size-1), i : i + (win_size-1), :);
% 
%         dark_channel(j,i) = min(patch(:));
%      end
% end
% 
% end

function dark_channel = get_dark_channel(I, neighborhood_size)
%GET_DARK_CHANNEL  Compute dark channel of input image with respect to a square
%neighborhood patch using erosion.
%
%   INPUTS:
%
%   -|I|: input color or grayscale image.
%
%   -|neighborhood_size|: the side of the square patch used for erosion, in
%   pixels.
%
%   OUTPUTS:
%
%   -|I_dark|: output grayscale image of the same type, height and width as |I|.
%
%   -|I_eroded|: intermediate eroded image of the same dimensions as |I|.

se = strel('square', neighborhood_size);

I_eroded = imerode(I, se);
dark_channel = min(I_eroded, [], 3);

end