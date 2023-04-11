clear all
format longG
% Numero de individuos
Nind = 60;
% Numero de variables
Nvar = 6;
% Arreglo de tamaño Nvar con los limites inferiores correspondientes
Li = [0 0 0 0 0.00001 0.00001];
% Arreglo de tamaño Nvar con los limites superiores correspondientes
Ls = [1 1 1 1 16 16];
%Numero de generaciones del genetico
Ngen = 150000;
%Factor de cruza
Fc = 0.7;
% Fc = 0.6 + (0.3)*rand();
%Factor de mutacion
Fm = 0.7;


rng('shuffle');
poblacion = crearPob(Li, Ls, Nind, Nvar);


FO = zeros(Nind, 1);
S = zeros(Nind, 1);
for i = 1:Nind
    FO(i) = funcionObjetivo(poblacion(i,:));
    g = restdes(poblacion(i,:));
    h = restigu(poblacion(i,:));
    S(i) = SVR(g, h);
end


for p = 1:Ngen
    p;
    u = zeros(size(poblacion));
    
    % Impresion de resultados
    bindex = indiceMejor(FO, S, Nind);
    poblacion(bindex,:)
    FO(bindex)
    S(bindex)

    Fm = 0.5 + (0.3)*rand();
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


function bindex = indiceMejor(FO, S, NH)
bindex = 1;
for i = 1:NH
    if DEB(FO(i), S(i), FO(bindex), S(bindex))
        bindex = i;
    end
end
end

function FO = funcionObjetivo(p)
    FO = -p(4);
end

function g = restdes(p)
    g = zeros(1,1);
    g(1) = p(5).^0.5 + p(6).^0.5 - 4;
end

function h = restigu(p)
    k1 = 0.09755988;
    k2 = 0.99*k1;
    k3 = 0.0391908;
    k4 = 0.9*k3;
    h = zeros(1,4);
    h(1) = p(1) + k1*p(2)*p(5) - 1;
    h(2) = p(2) - p(1) + k2*p(2)*p(6);
    h(3) = p(3) + p(1) + k3*p(3)*p(5) - 1; 
    h(4) = p(4) - p(3) + p(2) - p(1) + k4*p(4)*p(6);
end

function s = SVR(g, h)
    s = 0;
    for i = 1:size(g,2)
        s = s + max([0 g(i)]);
    end
    for i = 1:size(h,2)
        s = s + max([0 abs(h(i))-0.000001]);
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