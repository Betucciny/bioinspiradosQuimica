clear all
format longG
% Numero de individuos
Nind = 60;
% Numero de variables
Nvar = 8;
% Arreglo de tamaño Nvar con los limites inferiores correspondientes
Li = [100 1000 1000 10 10 10 10 10];
% Arreglo de tamaño Nvar con los limites superiores correspondientes
Ls = [10000 10000 10000 1000 1000 1000 1000 1000];
%Numero de generaciones del genetico
Ngen = 20000;
%Factor de cruza
Fc = 0.7;
%Factor de mutacion
Fm = 0.7;


rng('shuffle');
poblacion = crearPob(Li, Ls, Nind, Nvar);
colFO = Nvar + 1;

FO = zeros(Nind, 1);
S = zeros(Nind, 1);
for i = 1:Nind
    FO(i) = funcionObjetivo(poblacion(i,:));
    g = restdes(poblacion(i,:));
    h = restigu(poblacion(i,:));
    S(i) = SVR(g, h);
end

% archivoFO = fopen("bin/binDEBfo(30).txt", "w");
% archivoSR = fopen("bin/binDEBsr(30).txt", "w");
% final = fopen("bin/binDEB(30).txt", "w");

for p = 1:Ngen
    p;
    u = zeros(size(poblacion));
    
    % Impresion de resultados
    bindex = indiceMejor(FO, S, Nind);
    poblacion(bindex,:);
    FO(bindex);
    S(bindex);


    %Generacion de ruido
    for i = 1:Nind
        indices = randperm(Nind, 3);
        while any(indices==i)
            indices = randperm(Nind, 3);
        end
        u(i, :) = poblacion(indices(1),:) + Fm * (poblacion(indices(3),:) ...
            - poblacion(indices(2),:));
        for elem = 1:Nvar
            u(i,elem) = ajustar(u(i,elem), Li(elem), Ls(elem));
        end
    end

    newPob = zeros(size(poblacion));

    %Creacion de trial y eleccion
    for i = 1:Nind
        jrand = randi(Nvar);
        trial = zeros(1, Nvar);
        for j = 1:Nvar
            r = rand();
            if r < Fc || j==jrand
                trial(1, j) = u(i,j);
            else
                trial(1, j) = poblacion(i,j);
            end
        end

        FOtarget = FO(i);
        Starget = S(i);

        FOtrial = funcionObjetivo(trial);
        Strial = SVR(restdes(trial), restigu(trial));

        %         Reglas de DEB
        if DEB(FOtrial, Strial, FOtarget, Starget)
            newPob(i, :) = trial;
            FO(i) = FOtrial;
            S(i) = Strial;
        else
            newPob(i, :) = poblacion(i, :);
        end
        
    end 
    poblacion = newPob;    
end


bindex = indiceMejor(FO, S, Nind);
poblacion(bindex,:)
FO(bindex)
S(bindex)


function FO = funcionObjetivo(p)
    FO = p(1) + p(2) + p(3);
end

function g = restdes(p)
    g = zeros(1,3);
    g(1) = 100*p(1) - p(1)*p(6) + 833.33252*p(4) - 83333.333;
    g(2) = p(2)*p(4) - p(2)*p(7) - 1250*p(4) + 1250*p(5);
    g(3) = p(3)*p(5) - p(3)*p(8) - 2500*p(5) + 1250000;
end

function h = restigu(p)
    h = zeros(1,3);
    h(1) = 0.0025*(p(4) + p(6)) - 1;
    h(2) = 0.0025*(-p(4) + p(5) + p(7)) - 1;
    h(3) =  0.01*(-p(5) + p(8)) - 1; 
end

function s = SVR(g, h)
    s = 0;
    for i = 1:size(g,2)
        s = s + max([0 g(i)]);
    end
    for i = 1:size(h,2)
        s = s + max([0 abs(h(i))]);
    end
end

function mejor1= DEB(FO1, SVR1, FO2, SVR2)
FOtrial = FO1;
FOtarget = FO2;
Strial = SVR1;
Starget = SVR2;
if Starget == 0 && Strial == 0
    if FOtrial < FOtarget
        mejor1 = true;
        return
    else
        mejor1 = false;
        return
    end
elseif Starget ~= 0 && Strial ~= 0
    if Starget > Strial
        mejor1 = true;
        return
    else
        mejor1 = false;
        return
    end
else
    if Strial == 0
        mejor1 = true;
        return
    else
        mejor1 = false;
        return
    end
end
end

function pob = crearPob(li, ls, Nind, Nvar)
    pob = zeros(Nind, Nvar);
    for i=1:Nvar
        pob(:,i) = li(i) + (ls(i)- li(i))*rand(Nind, 1);
    end
end

function ajustado = ajustar(valor, li, ls)
    while true
        if valor < li
            valor = 2*li - valor;
        elseif valor > ls
            valor = 2*ls - valor;
        else
            break
        end
    end
    ajustado = valor;
end