function [intervalMean, intervalVariance] = IntervalStatistics(samples)

%IntervalStatistics
%   Author: Ander Gray   
%
%   IntervalStatistics is used to determine the sample mean and sample 
%   variance of an interval data set. The mean and variance of an interval
%   data set are also generally intervals.
%
%   [a b] = IntervalStatistics(samples) returns the two vectors a and b
%   where a is the sample mean and b the variance, and where samples is a
%   2 by N matrix containing the interval samples.
%   
%   The algorithm outline here follows the description of:
%   Ferson, "Experimental uncertainty estimation and statistics for data 
%   having interval uncertainty." Sandia National Laboratories
    
    if isequal(samples(:,1),samples(:,2))
       means = mean(samples(:,1));
       vars = var(samples(:,1));
       intervalMean = [means, means];
       intervalVariance = [vars,vars];
       return
    end

    intervalMean = [mean(samples(:,1)),mean(samples(:,2))];
    intervalVariance = [0,0];
    
    samplesSorted = sort(samples);
    Ys = reshape(samplesSorted',[numel(samples),1]);

    Ys = sort(Ys);

    overlap = zeros(numel(Ys)-1,1);
    Ns      = zeros(numel(overlap),1);
    Sks     = zeros(numel(overlap),1);
    Ms      = zeros(numel(overlap),1);

    N = numel(samples)/2;

    VarNull = 0;
    for i=1:numel(overlap)
        sorted = sort([Ys(i),Ys(i+1)]);
        if sorted(1)>intervalMean(2)
        elseif sorted(2)<intervalMean(1)
        else
            overlap(i)=1;
            leftBounds = samplesSorted(:,2) <= Ys(i);
            rightBounds = samplesSorted(:,1) >= Ys(i+1);
            Ns(i) = sum(leftBounds) + sum(rightBounds);
            if Ns(i) == 0
                VarNull = 1;
                break;
            end
            leftOk   = samplesSorted(:,2) .* leftBounds;
            rightOk  = samplesSorted(:,1) .* rightBounds;

            left = leftOk(leftOk ~=0);
            right = rightOk(rightOk ~=0);
            Sks(i) = sum(left) + sum(right);
            test = Sks(i)/Ns(i);
            if sorted(1)<= test && sorted(2) >= test
                if intervalMean(1)<= test && intervalMean(2) >= test
                    Ms(i) = (sum(left.^2) + sum(right.^2))/(numel(samplesSorted)/2);
                end
            end 
        end
    end

    if VarNull ==0
        MsOks      = ones(numel(Ms),1) .* Ms;
        SksFinal0  = Sks .* MsOks;
        NsFinal0   = Ns .* MsOks;

        MsFinal1   = Ms(Ms ~= 0);
        SksFinal1  = SksFinal0(SksFinal0 ~=0)./MsFinal1;
        NsFinal1   = NsFinal0(NsFinal0 ~=0)./MsFinal1;

        vars = MsFinal1 - SksFinal1.^2 ./((N .* NsFinal1));
        intervalVariance(1) = min(vars) * (N/(N-1));
    end
    p = num2cell(samples,2);
    cart = cartprod(p{:});
    vars = var(cart');
    intervalVariance(2) = max(vars);

end