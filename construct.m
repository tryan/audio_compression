function [yhat yhats] = construct(e, a)
% constructs a signal yhat from the input residuals e
% using equation (1) and coefficients a(k)
% yhats is yhat formatted for playback with sound(yhats,Fs)

nb = length(e(1,:));
blen = length(e(:,1));
yhats = zeros(1,blen*nb);
yhat = zeros(blen,nb);
for block = 2:nb
    coeffs = flipud(a(:,block));
    offset = 1 + (block-1)*blen;
    for n = offset:(offset+blen-1)
        ynk = (yhats((n-10):(n-1)));
        yhats(n) = ynk * coeffs + e(n-offset+1,block);
    end
    
    yhat(:,block) = yhats(offset:(offset+blen-1));
end

yhats = yhats';

end

