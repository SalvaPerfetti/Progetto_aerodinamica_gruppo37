function df = finiteDifferenceCentered(x, y)
    % Controllo input
    if length(x) ~= length(y)
        error('I vettori x e y devono avere la stessa lunghezza.');
    end
    
    n = length(x);
    df = zeros(1, n);  
    
    % passo (assumendo passo uniforme)
    h = mean(diff(x));
    
    % Differenze finite centrate (punti interni)
    for i = 2:n-1
        df(i) = (y(i+1) - y(i-1)) / (2 * h);
    end
    
    % estremi (usando differenze in avanti e indietro)
    df(1) = (y(2) - y(1)) / h;      % Differenza in avanti
    df(n) = (y(n) - y(n-1)) / h;    % Differenza indietro
end
