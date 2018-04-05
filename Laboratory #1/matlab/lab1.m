function lab1()
    clear all;
    
    X = readFromFile(); 
    X = [X zeros(1,150) + 13];
    
    X = sort(X);
    
    minX = X(1);
    fprintf('Mmin = %s\n', num2str(minX));

    maxX = X(end);
    fprintf('Mmax = %s\n', num2str(maxX));

    R = maxX - minX;
    fprintf('R = %s\n', num2str(R));

    mu = expectation(X);
    fprintf('mu = %s\n', num2str(mu));

    sigmaSqr = populationVariance(X);
    fprintf('sigma^2 = %s\n', num2str(sigmaSqr));

    sSqr = unbiasedSampleVariance(X);
    fprintf('S^2: %s\n', num2str(sSqr));

    m = numberOfSubintervals(length(X));
    fprintf('m = %s\n ', num2str(m));
    
    intervals(X, m);
    hold on;
    f(X, mu, sSqr);
    
    figure;
    empiricF(X);
    hold on;
    F(X, mu, sSqr, m);
        
end

function X = readFromFile()
    X = csvread('data.csv');
end

function m = numberOfSubintervals(size)
    m = floor(log2(size) + 2);
end

function mu = expectation(X)
    n   = length(X);
    sum = 0;

    for i = 1:n
        sum = sum + X(i);
    end
    
    mu = sum / n;
end

function sigmaSqr = populationVariance(X)
    n   = length(X);
    sum = 0;

    for i = 1:n
        sum = sum + (X(i))^2;
    end
    
    mu  = expectation(X);
    
    sigmaSqr = sum / n - mu^2;
end

function sSqr = unbiasedSampleVariance(X)
    sigmaSqr = populationVariance(X);
    n        = length(X); 

    sSqr = n / (n - 1) * sigmaSqr;
end

function intervals(X, m)
    count = zeros(1, m+1);  
    delta = (X(end) - X(1)) / m;
    
    J = X(1):delta:X(end);
    j = 1;
    n = length(X);
    
    for i = 1:n      
        if (j ~= m)
            if ((not (X(i) >= J(j) && X(i) < J(j+1))))
                j = j + 1;
                fprintf('[%.2f;%.2f)\t', J(j-1), J(j));
            end
        end
        count(j) = count(j) + 1;
    end
    fprintf('[%2.2f;%2.2f]\n', J(m), J(m + 1));
    
    Xbuf = count(1:m+1);
    for i = 1:m+1
        Xbuf(i) = count(i) / (n*delta); 
    end
    
    stairs(J, Xbuf), grid;
end

function f(X, MX, DX)
    Y = 1 / sqrt(2*pi*DX) * exp( -power((X - MX), 2) / (2*DX));
    plot(X, Y);
end

function F(X, MX, DX, m)
    R = X(end) - X(1);
    delta = R/m;
    
    Xn = (MX - R): delta :(MX + R);
    
    Y = 1/2 * (1 + erf((Xn - MX) / sqrt(2*DX))); 
    
    
    plot(Xn, Y, 'r');
end

function empiricF(X)
    n = length(X);
    Y = zeros(1, n);
 
    for i = 1:n
        Y(i) = i / n;
    end   
    
    [yy, xx] = ecdf(X);
    yy(length(yy)+1) = 1;
    xx(length(xx)+1) = 20;
    
    stairs(xx, yy);
end