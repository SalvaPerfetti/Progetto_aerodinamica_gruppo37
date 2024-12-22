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
end
