function xm = calcMidline(x1, x2)    
    % Verifica che i due vettori abbiano la stessa dimensione
    if length(x1) ~= length(x2)
        error('I due vettori devono avere la stessa lunghezza.');
    end
    
    % Calcola la linea media
    xm = (x1 + x2) / 2;
end