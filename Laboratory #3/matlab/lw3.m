function lw3()
    t0Array = csvread('t.csv');
    y0Array = csvread('y.csv');
    
    y1Array = OLS(y0Array, t0Array);
    
    delta = calculateDelta(y0Array, y1Array);
    fprintf('Delta = %3.2f', delta);
    
    plot(t0Array, y0Array, '.b', t0Array, y1Array, 'r');
    legend({
        'Original sample';
        'Output model';
    });
end

function Y = OLS(yArray, tArray)
    psiMatrix = makePsiMatrix(tArray);
    
    thetaArray = (psiMatrix' * psiMatrix) \ (psiMatrix' * yArray'); 
    disp(thetaArray);
    
    n = length(yArray);
    Y = zeros(1, n);
    for i = 1:n
        Y(i) = thetaArray(1) + ...
               thetaArray(2)*tArray(i) + ...
               thetaArray(3)*tArray(i)*tArray(i);
    end
end

function Psi = makePsiMatrix(tArray)    
    n   = length(tArray);
    p   = 3;
    Psi = zeros(n, p);
    for i = 1:n
        Psi(i, 1) = 1;
        Psi(i, 2) = tArray(i);
        Psi(i, 3) = tArray(i)*tArray(i);
    end
end

function X = calculateDelta(yArray, yNewArray)
    n   = length(yArray);
    sum = 0;
    for i = 1:n
        sum = sum + power((yArray(i) - yNewArray(i)), 2);
    end
    X = sqrt(sum);
end