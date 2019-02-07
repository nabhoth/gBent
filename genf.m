function [bFunctions, tFunctions, aFunctions, C, wlabels, tlabels] = genf( n )
%GENF this scripts generates spectra for all functions for a given n<=4 for a
%Boolean function
% additionally it generates two groups of functions, vFunctions the true
% bent function, bFunctions pseudo bent functions calculated from the truth
% vectors
W = [1 1; 1 -1];
Wn = W;
vl = (2^(n/2))/(2^n);
for i=1:n-1
    Wn = kron(W,Wn);
end
wlabels =[];
tlabels =[];
%true bent functions identified through Walsh transformfrom Walsh
%representation
bFunctions = [];
%pseudo bent functions identified through Walsh transform from true
%representation
tFunctions = [];
%all functions
aFunctions = zeros(2^2^n,2^n+2^n+2^n+2^n);
%counter of true bent functions for the bFunctions accumulator
vCounter = 0;
wCounter = 0;
%vector of functions in the walsh basis
v = ones(2^n,1);
v = -v;
%vector of functions in truth basis
w = zeros(2^n,1);
a = 1;
b = 1;
for i = 1:2^(2^n)

    %calculate the walsh transform from the truth vector
    wc = (Wn*w);
    wec = abs(wc);
    wec = wec/(2^n);
    %calculate the walsh transform from the walsh vector
    vc = (Wn*v);
    vec = abs(vc);
    vec = vec/(2^n);


    aFunctions(b,1:end) = [w', wec', v', vec'];
    b =  b + 1;
    %Bent From H transform
    if (max (vec) == vl && min (vec) == vl)
        vCounter = vCounter +1;
        nl = [w', wec', v', vec'];
        a =  a+ 1;
        bFunctions = [bFunctions; nl];
        wlabels = [wlabels; 1];
    else
        wlabels = [wlabels; 0];        
    end
    %Bent from Truth Vector    
    if (max(wec(2:end,1)')  == min(wec(2:end,1)'))
        wCounter = wCounter +1;
        nl = [w', wec', v', vec'];
        tFunctions = [tFunctions; nl];
        tlabels = [wlabels; 1];
    else
        tlabels = [wlabels; 0];        
    end
    v(1,1) = -v(1,1);
    w(1,1) = mod(w(1,1)+1,2);
    for j=2:2^n
        if (v(j-1,1) > 0)
            break;
        else
            v(j,1) = -v(j,1);
            w(j,1) = mod(w(j,1)+1,2);
        end
    end
end

[C, iv, ib] = intersect(bFunctions(:,1:16), tFunctions(:, 1:16),'rows');
[a,b] = size(bFunctions)
nvIndex = ones(a,1);
[a,b] = size(ib)
for i=1:a
    nvIndex(ib(i,1),1) = 0;
end
end

    

