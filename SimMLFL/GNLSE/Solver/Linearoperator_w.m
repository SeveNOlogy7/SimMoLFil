function [ LOP ] = Linearoperator_w( alpha,betaw,w )
%Linearoperator_w

LOP = -fftshift(alpha/2);
if (length(betaw) == length(w))     % If the user manually specifies beta(w)
    LOP = LOP - 1i*betaw;
    LOP = fftshift(LOP);
else
    for ii = 0:length(betaw)-1
        LOP = LOP - 1i*betaw(ii+1)*(w).^ii/factorial(ii);
    end
end

end

