%**************************************************************
function output = compute_gauss_smooth(input, window)
% Smoothing function: 
% output = smooth(input, window)
% "Window" is the total kernel width.
% Input array must be one-dimensional.

kerneltype=1;

input_dims = ndims(input);
input_size = size(input);
if input_dims > 2 | min(input_size) > 1,
   disp('Input array is too large.');
   return
end

if input_size(2) > input_size(1),
   input = input';
   toggle_dims = 1;
else
   toggle_dims = 0;
end

if window/2 ~= round(window/2),
   window = window + 1;
end
halfwin = window/2;

input_length = length(input);

if kerneltype==1
    %********* gauss window +/- 1 sigma
    x = -halfwin:1:halfwin;
    kernel = exp(-x.^2/(window/2)^2);
    kernel = kernel/sum(kernel);
else
    %   disp('Convolving with EPSP of GROWTH=1ms and DECAY=20ms')
    Growth = 1;
    Decay = 20;
    Half_BW=round(Decay*8);
    BinSize=(Half_BW*2)+1;
    Kernel=[0:Half_BW];
    Half_Kernel=(1-(exp(-(Kernel./Growth)))).*(exp(-(Kernel./Decay)));
    Half_Kernel=Half_Kernel./sum(Half_Kernel);
    Kernel(1:Half_BW)=0;
    kernel(Half_BW+1:BinSize)=Half_Kernel;
%    kernel=Kernel.*1000;
end


padded(halfwin+1:input_length+halfwin) = input;
padded(1:halfwin) = ones(halfwin, 1)*input(1);
padded(length(padded)+1:length(padded)+halfwin) = ones(halfwin, 1)*input(input_length);

output = conv(padded, kernel);
output = output(window:input_length+window-1);

if toggle_dims == 1,
   output = output';
end

return;




