%% Thomas Kok
% Audio Signal Compression
clc; clear;
[y, Fs] = wavread('test.wav');
y = y / max(abs(y)); %normalize the sound to [-1.0, 1.0]

% break sound file into block_len sample blocks
block_len = 160; % samples per block
nb = floor(length(y)/block_len); % number of blocks in the signal
blocks = zeros(block_len,nb); % each column contains one block of sound
for ind = 1:nb
    yind = 1 + (ind - 1) * block_len;
    blocks(:,ind) = y(yind:(yind+block_len-1));
end
y = y(1:block_len*nb); %eliminate the extra data points
clear yind ind
%% direct quantization 
r = 3:6; %number of quantization bits
alpha_range = 1:.25:5; %cutoff scalar for outlier removal
alpha_list_direct = zeros(1,length(r));
MSEdirect = zeros(1,length(r));

for n = 1:length(r)
    % find the best alpha for each value of r
    leastMSE = 1000;
    for b = 1:length(alpha_range)
        [~, ~, MSE] = quantize(blocks, alpha_range(b), r(n));
        if MSE < leastMSE
            leastMSE = MSE;
            alpha_list_direct(n) = alpha_range(b);
        end
    end
    
    [~, yqs, MSEq] = quantize(blocks, alpha_list_direct(n), r(n));
    MSEdirect(n) = MSEq;
    sound(yqs,Fs); pause;
end

clear n alpha yq MSEq leastMSE alpha_range MSE b yqs
%% solve for a(k) and e(n)
% filter coefficient calculation for each block
a = zeros(10,nb);
e = zeros(block_len,nb);
for block = 2:nb
    last = fliplr(blocks(block_len - 9:block_len, block - 1)');
    col = [blocks(block_len,block - 1) ; blocks(1:block_len - 1,block)];
    A = toeplitz(col,last);
    
    a(:,block) = A\blocks(:,block); %filter coefficients
    e(:,block) = blocks(:,block) - A * a(:,block);
end

clear A block last ind B col
%% get residuals using eqn (2), a(k) and y(n)
e2 = zeros(block_len,nb);
for block = 2:nb
    coeffs = flipud(a(:,block));
    offset = 1 + (block-1)*block_len;
    
    for n = offset:(offset+block_len-1)
        ynk = (y((n-10):(n-1))');
        e2(n-offset+1,block) = y(n) - ynk*coeffs;
    end
end

% get the MSE for each block
MSEe = zeros(nb,1);
for block = 1:nb
    m = e(:,block);
    n = e2(:,block);
    
    MSEe(block) = (m-n)' * (m-n) / block_len;
end

clear block coeffs ynk offset n m e2

%% reconstruct y from e

[yhat , ~] = construct(e, a);
% sound(yhats,Fs); %listen to reconstructed signal

% get the MSE for each block in the reconstructed y(n)
MSEy = zeros(1,nb);
for block = 2:nb
    m = blocks(:,block);
    n = yhat(:,block);
    
    MSEy(block) = (m-n)' * (m-n) / block_len;
end

clear m n block yhat
%% residual quantization
alpha_range = 1:.25:5;
alpha_list_eq = zeros(1,length(r));
MSEeq = r;

for n = 1:length(r) %cutoff scalar for outlier removal
    %get the best value of alpha
    leastMSE = 1000;
    for b = 1:length(alpha_range)
        [~, ~, MSE] = quantize(e, alpha_range(b), r(n));
        if MSE < leastMSE
            leastMSE = MSE;
            alpha_list_eq(n) = alpha_range(b);
        end
    end
    
    %quantize the residuals
    [eq, ~, ~] = quantize(e, alpha_list_eq(n), r(n));
    
    % reconstruct the signal from quantized residuals
    [~, yhats] = construct(eq, a);
    sound(yhats,Fs); pause;
    MSEeq(n) = (y-yhats)' * (y-yhats) / length(y);

end

clear n yhat alpha_range leastMSE MSE b yhats

%% residual and filter coefficient quantization
MSEeq_aq = r;
alpha_range = 1:.25:5;
alpha_list_aq = r;

for n = 1:length(r) %cutoff scalar for outlier removal
    %quantize filter coefficients
    leastMSE = 1000;
    for b = 1:length(alpha_range)
        [~, ~, MSE] = quantize(a, alpha_range(b), r(n));
        if MSE < leastMSE
            leastMSE = MSE;
            alpha_list_aq(n) = alpha_range(b);
        end
    end
    %quantize the residuals and coefficients
    [eq, ~, ~] = quantize(e, alpha_list_eq(n), r(n));
    [aq, ~, ~] = quantize(a, alpha_list_aq(n), r(n));
    
    % reconstruct the signal from quantized residuals
    [~, yhats] = construct(eq, aq);
%     sound(yhat,Fs); pause;
    MSEeq_aq(n) = (y-yhats)' * (y-yhats) / length(y);

end

clear eq eqs aq aqs n yhat yhats MSEeq MSEaq r alpha_range MSE leastMSE
%% data visualisation

% subplot(2,1,1); hold on; grid on;
figure; hold on; grid on; box on;
scatter(r,MSEeq);
scatter(r,MSEdirect,'red');
legend('Residual Quantization','Direct Quantization');
title('Comparison of MSE for Direct and Residual Quantization');
xlabel('Bits of Quantization');
ylabel('Mean Square Error');
axis([0 9 0 3]);

% subplot(2,1,2); hold on; grid on;
figure; hold on; grid on; box on;
scatter(r(4:length(r)),MSEeq(4:length(r)))
scatter(r(4:length(r)),MSEdirect(4:length(r)),'red');
title('A Closer Look');
xlabel('Bits of Quantization');
ylabel('Mean Square Error');
axis([3.5 8.5 0 .08]);











