function xm = calcMidline(x1, x2)
    % calcMidline - Calcola la linea media tra due vettori di punti.
    %
    % Sintassi:
    %    xm = calcMidline(x1, x2)
    %
    % Input:
    %    x1 - Primo vettore di punti (n x 1 o 1 x n)
    %    x2 - Secondo vettore di punti (n x 1 o 1 x n)
    %
    % Output:
    %    xm - Vettore della linea media (n x 1 o 1 x n)
    %
    % Esempio:
    %    x1 = [1, 2, 3];
    %    x2 = [4, 5, 6];
    %    xm = calcMidline(x1, x2);
    %    disp(xm); % Output: [2.5, 3.5, 4.5]
    
    % Verifica che i due vettori abbiano la stessa dimensione
    if length(x1) ~= length(x2)
        error('I due vettori devono avere la stessa lunghezza.');
    end
    
    % Calcola la linea media
    xm = (x1 + x2) / 2;
end