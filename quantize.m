function [yq yqs MSEq] = quantize(signal, alpha, r)
% blockwise quantization of the blen x nb matrix 'signal'
% output is a blen x nb matrix of the quantized signal
% input 'signal' is a blen x nb matrix
% alpha is the cutoff scalar for outlier removal
% r is the # of bits used in the quantized output signal

nb = length(signal(1,:));
blen = length(signal(:,1));
yq = zeros(blen,nb);
yqs = zeros(blen*nb,1); %put yq in listenable form
L = 2^r; %number of quantization levels
MSEq = zeros(1,nb);

for block = 2:nb
    dev = std(signal(:,block));
    avg = mean(signal(:,block));
    
    ytmax = avg + alpha * dev; %upper bound for ythresh
    ytmin = avg - alpha * dev; %lower bound for ythresh

    ythresh = zeros(blen,1); %truncate block
    for ind = 1:blen
        if signal(ind,block) > ytmax
            ythresh(ind) = ytmax;
        elseif signal(ind,block) < ytmin
            ythresh(ind) = ytmin;
        else
            ythresh(ind) = signal(ind,block);
        end
    end

    %direct quantization to form yq 
    q = (max(ythresh) - min(ythresh))/L;
    yq(:,block) = round(ythresh/q)*q;
    
    offset = 1 + (block-1)*blen;
    yqs(offset:(offset+blen-1)) = yq(:,block);
    
    m = yq(:,block);
    n = signal(:,block);
    MSEq(block) = (m-n)' * (m-n); % store the MSE for each block
end

MSEq = sum(MSEq) / blen; % return the MSE for the entire signal

end

