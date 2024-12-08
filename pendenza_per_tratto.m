function dylm = pendenza_per_tratto(x,y)
    % Calcola la pendenza tra punti successivi
    % Input:
    %   x - vettore delle coordinate x
    %   y - vettore delle coordinate y
    % Output:
    %   dylm - vettore delle pendenze (dy/dx) tra punti successivi

    % Calcolo delle differenze tra punti successivi
    dx = diff(x);
    dy = diff(y);

    % Calcolo delle pendenze (dy/dx)
    dylm = dy ./ dx;

    % Aggiungere NaN o 0 alla fine per mantenere stessa lunghezza
    % Opzionale: assumiamo la pendenza dell'ultimo tratto uguale all'ultimo valore
    dylm = [dylm; NaN]; % NaN indica un'assenza di valore (l'ultimo tratto non esiste)
end
