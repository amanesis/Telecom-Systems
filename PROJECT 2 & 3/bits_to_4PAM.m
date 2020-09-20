function [ Xk ] = bits_to_4PAM( b )
% Xk = bits_to_4PAM(b)
% OUTPUT
% Xk: Xk, k=0,...,(N/2)-1
% INPUT
% b: sequence of bits
for k=1:2:length(b)/2
    if(b(k)==0 && b(k+1)==0)
        Xk(k)=3;
    elseif(b(k)==0 && b(k+1)==1)
        Xk(k)=1;
    elseif(b(k)==1 && b(k+1)==1)
        Xk(k)=-1;
    elseif(b(k)==1 && b(k+1)==0)
        Xk(k)=-3;
    else
        disp('Error')
        return
    end
end