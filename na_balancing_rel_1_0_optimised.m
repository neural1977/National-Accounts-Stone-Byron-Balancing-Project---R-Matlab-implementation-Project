%%        National Accounts Balancing Procedure                                              
%                                                                                 %
%         Implementation: Matlab Language                                         %
%         Date of creation: 20/01/2014                                            %
%         Date of last change: 17/02/2014                                         %
%         Release Version: 1.0                                                    %                          
%         Release Description: Translation of the original Nicolardi's procedure  % 
%                                                                                 %
%         Authors: Francesco Pugliese                                             %
%         Email: frpuglie@istat.it                                                %                          
%                                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Note Generali
% * Usando una serie di funzioni possiamo evitare di dover cancellare
% manualmente variabili non più utilizzate
% 


%% Inizializzazione delle variabile e dell'ambiente
format long g;
clear;                                                                                                                   % azzera il workspace      

%% Data, personalize the following folders according your needs
input_dir='D:/Documenti Utente/Google Drive/National Accounts/Stone Balancing/Input Balancing/Input/Data/';
output_dir='D:/Documenti Utente/Google Drive/National Accounts/Stone Balancing/Input Balancing/Output Matlab Optimised/';
variances_dir='D:/Documenti Utente/Google Drive/National Accounts/Stone Balancing/Input Balancing/Input/Variances/';
config_dir='D:/Documenti Utente/Google Drive/National Accounts/Stone Balancing/Input Balancing/';
config_file='vincoli.deck';

%% Globals
config_path=strcat(config_dir, config_file); 

read_format=1;                                                                                                           % 0: legge da file excel, 1: leggge da file .csv, 2: legge da file .dif
write_format=0;                                                                                                          % 0: scrive su file excel, 1: scrive su file .csv, 2: scrive su file .dif

csv_delimiter=';';
variances_prefix='V'; 

%%    BLOCCO LOGICO 1 
tim_g=tic;     
% tim_g=clock;      % cattura il tempo d'inizio --- Valutare uso TIMER

f = fopen(config_path,'r');                                                                                             
if f>0
    X = textscan(f,'%s','delimiter','\n');
    fclose(f);
else
    disp(['Errore nella lettura del file di configurazione: ',config_file]);
    return;
end
    
nofeq = length(X{1});                                                                                                     % numero di macro equazioni nel file di configurazione

%%     BLOCCO LOGICO 2     

%Questo è un esempio analogo al codice originale di sinteticità del
%linguaggio MATLAB - Inoltre si fa uso di preallocazione che riduce la
%frammentazione della memoria e il tempo di calcolo

Y=X{1}; % prende la prima colonna del cell array

cY =  cellfun(@(x)length(strsplit(x)),Y); % Dimensione dei campi
nofet =max(cY);
% Preallocate 
config = cell(nofeq,nofet);
config(:) = {' '}; %Inizializzo con spazio vuoto come in esempio
for i = 1:nofeq
   row = strsplit(Y{i});
   config(i,1:length(row)) = row;
end

clear nfeq;

%%     BLOCCO LOGICO 3     


%[nofeq, nofet] = size(config); Calcolati già in anticipo
nofef=nofet/3;

%% BLOCCO LOGICO 4
equat = config(:,3:3:nofet);


%% BLOCCO LOGICO 5     

% Crea l'array dei segni, usando un appoggio
appo_signs = cell2mat(config(:,1:3:nofet));
                                                                                    
%% BLOCCO LOGICO 6
signs = zeros(nofeq, nofef);
signs(appo_signs == '+') = 1;
signs(appo_signs == '-') = -1;

clear appo_signs;

direct = config(:,2:3:nofet);

%% BLOCCO LOGICO 7
tmp2 = direct;
appo_direct = direct;

direct = cell(size(appo_direct));
direct(:) = {NaN};
                                                                                     
direct(strcmp(appo_direct,'SR'))={1}; %Usiamo direttamente i numeri non ha più senso usare le stringhe e lavorare con i codici ASCII
direct(strcmp(appo_direct,'SC'))={2};
direct(strcmp(appo_direct,'SM'))={3};
direct(strcmp(appo_direct,'MM') | strcmp(appo_direct,'VC') | strcmp(appo_direct,'VR'))={4};

direct=cell2mat(direct);  % converte direct in un array numerico   

%% BLOCCO LOGICO 8
dimeq = zeros(nofeq, 3); 
tmp2=tmp2(:,1);

for i = 1:size(equat,1)
    switch read_format
        case 0
            
        case 1
            tmp=dlmread(strcat(input_dir,equat{i,1},'.csv'),csv_delimiter);
        case 2
            
        otherwise
            error('Read Format Non Corretto')
    end

    tmpr=size(tmp,1); 
    tmpc=size(tmp,2); 
    tmpe=length(tmp);
    
    switch tmp2{i} %Swithch case è più sintetico di un IF in questo caso
        case 'MM' 
        dimeq(i,1)=tmpr; 
        dimeq(i,2)=tmpc; 
   
        case 'SR'
        dimeq(i,1)=tmpr; 
        dimeq(i,2)=1; 
 
        case 'VC' 
        dimeq(i,1)=tmpe; 
        dimeq(i,2)=1; 

        case 'SC'
        dimeq(i,1)=1; 
        dimeq(i,2)=tmpc; 

        case 'VR'
        dimeq(i,1)=1; 
        dimeq(i,2)=tmpe; 

        case 'SM'
        dimeq(i,1)=1; 
        dimeq(i,2)=1; 
    end
end
dimeq(:,3) = dimeq(:,1).* dimeq(:,2);

%% BLOCCO LOGICO 10
%  COSTRUZIONE STRUMENTI PER IL CALCOLO DEL BILANCIAMENTO   
% Uso il comando UNIQUE, codice più compatto e leggibile

equatt = equat'; % Traspongo perché usiamo logica per righe
pos = find(direct' == 4); %Salvo le posizioni originali
[mats,imats] = unique(equatt(pos),'stable'); %Estraggo in ordine i valori unici originali

[rm,cm] = ind2sub(size(equatt),pos(imats)); %Ricavo le corrispondenti righe (colonne trasposte) per indicizzare dimats
dimats = dimeq(cm,1:2); %Estraggo

nofmats=length(mats);

%% BLOCCO LOGICO 11

% Definisce un array "matstin" che contiene la posizione di mats all'interno dei vincoli
matstin = zeros(size(equat)); %Probabilmente possiamo riutilizzare parte dell'approccio sopra
for i = 1 : nofmats
    matstin(strcmp(equat, mats(i))) = i;
end

%% BLOCCO LOGICO 12,13  

% Blocco che genera gli importanti array pointeq e signseq
pointeq=zeros(nofmats, nofef); 
signseq=zeros(nofmats, nofef);

for i = 1:nofmats
    tmp = zeros(size(equat));
    tmp(strcmp(equat, mats(i))) = 1;
    tmp = sum(tmp.*signs,2);                                                                                                  % somma riga
    tmploc = find(tmp~=0);
    pointeq(i,1:length(tmploc)) = tmploc;
    signseq(i,1:length(tmploc)) = tmp(tmploc);
end;
equat(strcmp(equat, ' '))={'x0'};

%% BLOCCO LOGICO 14,15

% Uso FPRINTF perchè più performante
fprintf('\n\n\nNumero di equazioni nel sistema : %.0f\n', sum(dimeq(:,3)));
fprintf('Numero di elementi della matrice contabile : %.0f\n', sum(dimats(:,1).*dimats(:,2)));

% Potremmo parallelizzare le letture (se questo ha un senso da un punto di
% vista di performance - Al momento il tutto è comunque molto rapido)

oldcont = 1e40;
contres = 0;
dimvec = sum(dimeq(:,3));                                                                                                     % dimveq contiene il numero totale di elementi presenti nelle macroequazioni

% lettura dei dati di input
for i = 1:nofmats
    switch read_format
        case 0
        case 1
            appo = dlmread(strcat(input_dir,mats{i},'.csv'),csv_delimiter);
        case 2
    end;
    eval(strcat(mats{i},' = appo;'));
end;

% lettura dei dati delle varianze
for i = 1:nofmats
    switch read_format
        case 0
        case 1
            appov = dlmread(strcat(variances_dir,variances_prefix,mats{i},'.csv'),csv_delimiter);
        case 2
    end;
    eval(strcat(variances_prefix,mats{i},' = appov;'));
end

%% BLOCCO LOGICO 16
tim=tic;

% Preparazione dei residui per il ciclo del gradiente coniugato vero e proprio 
%res=zeros(dimvec); resc=zeros(dimvec);
res=zeros(dimvec,1); 
resc=zeros(dimvec,1);
top=0; 
bot=0;
% In questo pezzo sembra applicare i vincoli matrice e di bilanciamento riga per riga. 
% In pratica legge ogni riga e applica i vincoli letti venendo direzionato dai meta-tag del linguaggio eulero. 
% Nel caso siano vincoli di matrice (tipo somma tra matrici con MM, VC o VR) effettua semplicemente una somma tra matrici, cumulandola in tmp attraverso una somma (caso 4). 
% Ogni tmp contiene il risultato delle operazioni lineari di matrice di quella riga (residuo). Tmp e Wtmp in questo caso sono matrici.  
% Nel caso siano vincoli di bilanciamento (tipo somma di riga o colonna con SR o SC) effettua una somma prima di riga o di colonna degli elementi, e poi cumula tutto in tmp (caso 1 o 2). 
% Ogni tmp contiene il risultato delle operazioni di bilanciamento  di quella riga (residuo). Tmp e Wtmp in questo caso sono vettori.
% Nel caso siano vincoli di somma degli elementi di matrice (SM) allora effettua tutta la somma e la cumula in uno scalare tmp la somma degli elementi, e poi cumula tutto in tmp (caso 3). 
% Ogni tmp contiene il risultato delle operazioni di bilanciamento  di quella riga (residuo). Tmp e Wtmp in questo caso sono scalari.

for i = 1:nofeq
  tmp = 0; 
  wtmp = 0;

  for j = 1:nofef

    switch direct(i,j) 
    case 1
                                                                                      % Caso =1 (tag SR del linguaggio EULERO)
        tmp=tmp+sum(eval(strcat(equat{i,j},'.*signs(i,j)')),2);                       % Moltiplica la matrice alla posizione i,j (del file di configurazione equat) per il segno relativo, e ne calcola la somma delle righe,   
                                                                                     % proprio come indicato dal tag SR. Per esempio di sommarighe(MATRIX * 1)
        wtmp=wtmp+sum(eval(strcat(variances_prefix,equat{i,j})),2);                                % Somma le righe del corrispondente elemento delle tolleranze (varianze?) che è preceduto da V      
    case 2
                                                                                      % Caso =2 (tag SC del linguaggio EULERO) 
        tmp=tmp+sum(eval(strcat(equat{i,j},'.*signs(i,j)')),1);                        % Moltiplica la matrice alla posizion i,j (del file di configurazione equat) per il segno relativo, e ne calcola la somma delle colonne,
                                                                                      % proprio come indicato dal tag SC. Per esempio di sommarighe(MATRIX * 1)
        wtmp=wtmp+sum(eval(strcat(variances_prefix,equat{i,j})),1);                                % Somma le righe del corrispondente elemento delle tolleranze (varianze?) che è preceduto da V
    case 3
                                                                                      % Caso =3 (tag SM del linguaggio EULERO) 
        tmp=tmp+sum(sum(eval(strcat(equat{i,j},'.*signs(i,j)'))));                    % Effettua la somma di tutti gli elementi della matrice o vettore appena letto ed inserisce il valore in tmp (attenzione è uno scalare) 
        wtmp=wtmp+sum(eval(strcat(variances_prefix,equat{i,j})));                                  % Effettua la somma di tutti gli elementi della matrice o vettore delle tolleranze appena letto ed inserisce il valore in wtmp (attenzione è uno scalare) 
    case 4
                                                                                      % Caso =4 (tag MM, VC o VR del linguaggio EULERO) 
        tmp=tmp+eval(strcat(equat{i,j},'.*signs(i,j)'));                              % Effettua la somma elemento per elemento tra matrici o vettori della riga corrente e cumula tutto in tmp.   
        wtmp=wtmp+eval(strcat(variances_prefix,equat{i,j}));                                       % La stessa operazione viene effettuata per le tolleranze (varianze). 
    otherwise
        break;                                                                        % Caso =5 Causa l'uscita dal ciclo, con questa opzione   
    end
  end
  
 
  % Definisce un array logico (vero/falso) per verificare se ci sono varianze = 0
  aa = (tmp ~= 0) & (wtmp == 0);                                                      % Verifica che per ogni residuo, la varianza non sia uguale a zero (defaul=FALSE)
  verif = 0;
  if sum(sum(aa)) ~= 0                                                               
    verif = 1;
    posiz = aa==1;
    residd = tmp(posiz);
    varian = wtmp(posiz);
    fprintf('\n\n\nNella macro-eq. %.0f vi sono residui a varianza 0 nelle microequazioni\n', i);
    fprintf('%.0f %f %f\n', posiz, resid, varian);
  end
  
  
  
%% BLOCCO LOGICO 17
  
  w = 1./sqrt(wtmp);                                                                   % Inizializza la matrice di precondizionamento linearizzata
  w(isnan(w)) = 0;
  w(w == Inf) = 0;
  eval(strcat('w',num2str(i),'=w;'));                                                  % crea la matrice w1, w2, ecc...
  
  pi=tmp.*w;                                                                           % Inizializza il parametro phi_i del procedimento di gradiente coniugato
  pi(isnan(pi)) = 0;
  pi(pi==Inf) = 0;
  eval(strcat('pi',num2str(i),'=pi;'));                                                % crea la matrice pi1, pi2, ecc...

  ri = tmp.*w;                                                                           % Inizializza il parametro rho_i del procedimento di gradiente coniugato
  ri(isnan(ri)) = 0;
  ri(ri==Inf) = 0;
  eval(strcat('ri',num2str(i),'=ri;'));                                                % crea la matrice ri1, ri2, ecc...

  li=tmp.*0;                                                                           % Inizializza il parametro lambda_i del procedimento di gradiente coniugato
  li(isnan(li)) = 0;
  li(li==Inf) = 0; %#ok<NASGU> 
  eval(strcat('li',num2str(i),'=li;'));                                                % crea la matrice li1, li2, ecc...
  
  contres = contres + sum(sum(tmp));                                                       % cumula la somma all'interno di tmp
  bot = top + numel(tmp);                                                                  % determina estremi inferiore, top, ed estremo superiore, bot, come indici per res 
  eval(strcat('res((top+1):bot)=pi',num2str(i),';'));                                  % cumula i vari pi di ogni iterazione per i (pi1, pi2, ..) in res, usa eval e parse per poter adattare i vari pi1, pi2 alla forma espressiva 
  
  top = bot;                                                                        % il primo indice diventa l'ultimo dell'iterazione precedente   
  
end

if verif == 1 
  error('errore nelle varianze: verificare il sistema\n');                             % Errore causa l'uscita dal programma   
end

% Dealloca tutte le strutture superflue
clear tmp wtmp;                                                     

% NOTA!!! per il momento è commentato, se non ce la fa a tenerli in memoria allora bisogna rileggerli alla fine del programma
%for (i in 1:nofmats) 
%  eval(parse(text=paste("rm(",mats[i],")",sep="")))                                % Libera l'area di memoria occupata da tutti gli array elencati dentro mats: es. rm(OFFERT), rm(DOMAND), ecc.. con mats(1)=OFFERT, mats(2)=DOMAND,...

% Quando somma righe (marginale di colonna) è uguale al vettore colonna (informazione esterna) e somma colonne (marginale di riga) è uguale al vettore riga (informazione esterna) allora SISTEMA BILANCIATO
if contres == 0
    error('errore di definizione: il sistema è già bilanciato\n');                     % Errore causa l'uscita dal programma   
end

%% BLOCCO LOGICO 18

%%%%%%%%%%%%%%%%%%%   INIZIO  ITERAZIONI  PER   IL  BILANCIAMENTO   %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mancano controlli e commenti, che però potrebbero risultare inutili
iter=0; limit=0.0000001; li=zeros(dimvec,1); counter=4; % tim=proc.time()

fprintf('\n\n\n*************    INIZIO ITERAZIONI PER IL BILANCIAMENTO (GRADIENTE CONIUGATO)   *****************\n\n');
pp = @plus;

while 1==1 
  counter = counter + 1;

  % all'inizio abbiamo i pi, li e ri di ogni riga: pi1, pi2, ecc... (che sono i parametri del gradiente coniugato) 
  for i = 1:nofmats 
    tout=0;
    for j = 1:size(pointeq,2) 
      if pointeq(i,j) > 0                                                                                                       % se l'indice è diverso da zero   
       % allora ad un certo punto può trovarsi a sommare i valori del vettore riga e del vettore colonna, il risultato è un oggetto bidimensionale con: 
       % numero righe = numero di righe del vettore colonna
       % numero di colonne = numero di colonne del vettore riga
   
       % funzione speakeasy sum che prevede anche il riciclo nel caso due vettore abbiano differente dimensione
       
       tout = speakeasy_recycling_sum_optimised(tout, eval(strcat('pi',num2str(pointeq(i,j)),'.*signseq(i,j).*w',num2str(pointeq(i,j)))),pp);      % da quanto ho capito, per ogni elemento presente in mats, costruisce dei tout che contengono dei pi * w calcolati nel
                                                                                                                                    % seguente modo: considera i parametri del gradiente coniugato solo per le righe interessate dagli oggetti per esempio, 
                                                                                                                                    % oggetto OFFERT si trova alla riga 1, allora moltiplica pi1 * (il segno di offert nella riga 1) * w1. Alla fine cumula      end
      end
    end 
      
    eval(strcat('tout',num2str(i),'=tout;'));                                                                                        % valuta l'espressione, e memorizza i vari tout nelle variabili tout1, tout2, tout3, ecc...
  end  
  


%%       BLOCCO LOGICO 20

  % Ciclo che riapplica i vincoli matrice e di bilanciamento riga per riga, ma questa volta alle strutture moltiplicate, es: tout1 (pi1*w1)  
  % In pratica legge ogni riga e applica i vincoli letti venendo direzionato dai meta-tag del linguaggio eulero. 
  clear tout;  
  ai=0; rr=0;

  for i = 1:nofeq
    tmp = 0;
    rr=rr+sum(sum(eval(strcat('ri',num2str(i),'.^2'))));
   
    for j = 1:nofef
        kk=matstin(i,j); 
    
        switch direct(i,j) 
        case 1
                                                                                                                    % Caso =1 (tag SR del linguaggio EULERO)
            tmp=tmp+sum(eval(strcat('tout',num2str(kk),'.*',variances_prefix,mats{kk},'.*signs(i,j)')),2);          % Moltiplica i tout (ottenuti per ogni oggetto di input) per la relativa matrice delle varianze,
                                                                                                                    % e ne calcola la somma per riga, proprio come indicato dal tag SR. 
        case 2
                                                                                                                    % Caso =2 (tag SC del linguaggio EULERO) 
            tmp=tmp+sum(eval(strcat('tout',num2str(kk),'.*',variances_prefix,mats{kk},'.*signs(i,j)')),1);          % Moltiplica i tout (ottenuti per ogni oggetto di input) per la relativa matrice delle varianze,
                                                                                                                    % e ne calcola la somma per colonna, proprio come indicato dal tag SC. 
        case 3
                                                                                                                    % Caso =3 (tag SM del linguaggio EULERO) 
            tmp=tmp+sum(sum(eval(strcat('tout',num2str(kk),'.*',variances_prefix,mats{kk},'.*signs(i,j)'))));       % Moltiplica i tout (ottenuti per ogni oggetto di input) per la relativa matrice delle varianze, ed
                                                                                                                    % effettua la somma di tutti gli elementi della matrice o vettore appena letto ed inserisce il valore in tmp (attenzione è uno scalare)          
        case 4
                                                                                                                    % Caso =4 (tag MM, VC o VR del linguaggio EULERO) 
            tmp=tmp+eval(strcat('tout',num2str(kk),'.*',variances_prefix,mats{kk},'.*signs(i,j)'));                 % Moltiplica i tout (ottenuti per ogni oggetto di input) per la relativa matrice delle varianze, ed
                                                                                                                    % effettua la somma elemento per elemento tra matrici o vettori della riga corrente e cumula tutto in tmp. 
        otherwise
            break;                                                                                                  % Caso =5 Causa l'uscita dal ciclo, con questa opzione   
        end
    end


%%       BLOCCO LOGICO 21

   eval(strcat('api',num2str(i),'=tmp.*w',num2str(i),';'));                                                          % crea la matrice api1, api2, ecc..., da quanto ho capito api coincidono con gli alpha_i

   ai=ai+sum(sum(eval(strcat('api',num2str(i),'.*pi',num2str(i)))));                                                 % Effettua il prodotto api1 * pi1 (phi_1), api2*pi2, ecc. e cumula tutto nello scalare ai
                                                                                                                     % I vari api, dovrebbero essere gli A*phi_i del gradiente coniugato, mentre ai, rappresenta il coefficiente alpha_i
                                                                                                                     % mentre i vari tout dovrebbero rappresentare l'equivalente della matrice: phi*w
                                                                                                                     % il tmp genera la matrice tout*V... ossia dovrebbe generare in output l'equivalente A*phi_i che poi viene giustamente moltiplicato per w
                                                                                                                     % il tutto bypassando ingegnosamente la costruzione della matrice G
  end
    
  ai=rr/ai;                                                                                                          % Calcola il rapporto tra la somma degli api_i*pi_i e la somma degli ri_i (rho_i). Questo coefficiente sarebbe alpha_i
    
%%      BLOCCO LOGICO 22      %


  % Dealloca tutte le strutture superflue
  clear tmp;
  for i = 1:nofmats                                                                                      
    eval(strcat('clear tout',num2str(i),';'));                                                                       % Libera l'area di memoria occupata da tutti gli array tout_1, tout_2, ecc.
  end
  
  rr1=0;
  for i = 1:nofeq
    eval(strcat('ri',num2str(i),'=ri',num2str(i),'-api',num2str(i),'.*ai;'));                                        % Qui aggiorna i rho_i secondo la consueta formula del gradiente coniugato: rho_i=rho_i-alpha_i * A %*% phi_i                                                                                                                    % A questo punto il programma di Nicolardi comincia ad assumere la parvenza della consueta procedura di gradiente coniugato
    rr1 = rr1 + sum(sum(eval(strcat('ri',num2str(i),'.^2'))));                                                           % Calcola la somma dei quadrati dei vari rho_i 
    eval(strcat('li',num2str(i),'=li',num2str(i),'+pi',num2str(i),'.*ai;'));                                         % Qui aggiorna i lambda_i secondo la consueta formula del gradiente coniugato: lambda_i=lambda_i + alpha_i * phi_i   
  end  

  bi=rr1/rr;                                                                                                         % Calcola il rapporto tra la somma degli ri^2 del passo attuale e quelli del passo precedente. a somma degli ri_i (rho_i). Questo coefficiente sarebbe beta_i.
                                                                                                                     % che è appunto un rapporto: beta_i=((t(rho_i) %*%  A %*% phi_i) / (t(phi_i) %*% A %*% phi_i)). Dunque rr coinciderebbe a t(rho_i) %
  for i = 1:nofeq 
    eval(strcat('pi',num2str(i),'=ri',num2str(i),'+pi',num2str(i),'*bi;'));                                          % Qui aggiorna i phi_i secondo la consueta formula del gradiente coniugato:  phi_i=-rho_i+beta_i * phi_i
  end  

  if (counter==5) 
    iter=iter+1;
 
%% BLOCCO LOGICO 23

    
  % ciclo di determinazione dei parametri moltiplicati per w, elemento per elemento
  for i = 1:nofmats 
    tout=0;
    for j = 1:size(pointeq,2) 
      if pointeq(i,j) > 0                                                                                                           % se l'indice è diverso da zero   
 
        % funzione speakeasy sum che prevede anche il riciclo nel caso due vettore abbiano differente dimensione
        tout=speakeasy_recycling_sum_optimised(tout, eval(strcat('li',num2str(pointeq(i,j)),'.*signseq(i,j).*w',num2str(pointeq(i,j)))),pp);     % da quanto ho capito, per ogni elemento presente in mats, costruisce dei tout che contengono dei li * w calcolati nel
                                                                                                                             % tutto il risultato in tout1 (che corrisponde al primo elemento in mats, es OFFERT) and so on...
      end
    end 
      
      eval(strcat('tout',num2str(i),'=tout;'));                                                                                     % valuta l'espressione, e memorizza i vari tout nelle variabili tout1, tout2, tout3, ecc...
  end  
   
  top=0;                                                                                                                            % indice per il vettore dei risultati finale
  
%% BLOCCO LOGICO 24

  % Ciclo che riapplica i vincoli matrice e di bilanciamento riga per riga, ma questa volta alle strutture moltiplicate, es: tout1 (li1*w1)  
  % In pratica legge ogni riga e applica i vincoli letti venendo direzionato dai meta-tag del linguaggio eulero. 
  % Praticamente riapplica lo schema tout-tout1,(2,3,ecc.) solo a phi_i e a lambda_i, ma non a rho_i, in ogni . Bisogna capire perchè. Forse perchè nei primi due rientra in gioco la G che in qualche modo se la determina con
  % facendosi guidare dal linguaggio eulero, mentre in rho_i non entra in gioco oppure rho_i se la ricava da phi_i.
    
    clear tout;  
    for i = 1:nofeq
        tmp = 0;
   
        for j = 1:nofef
            kk=matstin(i,j);                                                                                            % determina la posizione all'interno di equat di un oggetto di mats        
    
            switch direct(i,j) 
            case 1
                                                                                                                        % Caso =1 (tag SR del linguaggio EULERO)
                tmp=tmp+sum(eval(strcat('tout',num2str(kk),'.*',variances_prefix,mats{kk},'.*signs(i,j)')),2);          % Moltiplica i tout (ottenuti per ogni oggetto di input) per la relativa matrice delle varianze,
                                                                                                                        % e ne calcola la somma per riga, proprio come indicato dal tag SR. 
            case 2
                                                                                                                        % Caso =2 (tag SC del linguaggio EULERO) 
                tmp=tmp+sum(eval(strcat('tout',num2str(kk),'.*',variances_prefix,mats{kk},'.*signs(i,j)')),1);          % Moltiplica i tout (ottenuti per ogni oggetto di input) per la relativa matrice delle varianze,
                                                                                                                        % e ne calcola la somma per colonna, proprio come indicato dal tag SC. 
            case 3
                                                                                                                        % Caso =3 (tag SM del linguaggio EULERO) 
                tmp=tmp+sum(sum(eval(strcat('tout',num2str(kk),'.*',variances_prefix,mats{kk},'.*signs(i,j)'))));       % Moltiplica i tout (ottenuti per ogni oggetto di input) per la relativa matrice delle varianze, ed
                                                                                                                        % effettua la somma di tutti gli elementi della matrice o vettore appena letto ed inserisce il valore in tmp (attenzione è uno scalare)          
            case 4
                                                                                                                        % Caso =4 (tag MM, VC o VR del linguaggio EULERO) 
                tmp=tmp+eval(strcat('tout',num2str(kk),'.*',variances_prefix,mats{kk},'.*signs(i,j)'));                 % Moltiplica i tout (ottenuti per ogni oggetto di input) per la relativa matrice delle varianze, ed
                                                                                                                        % effettua la somma elemento per elemento tra matrici o vettori della riga corrente e cumula tutto in tmp. 
            otherwise
                break;                                                                                                  % Caso =5 Causa l'uscita dal ciclo, con questa opzione   
            end
        end
    
        bot=top+numel(tmp);                                                                                             % determina estremi inferiore, top, ed estremo superiore, bot, come indici per res 
        eval(strcat('resc((top+1):bot)=tmp.*w',num2str(i),';'));                                                        % cumula i vari tmp*wi di ogni iterazione per i, in resc, usa eval e parse per poter adattare i vari tmp*w1, tmp*w2 alla forma espressiva            
        top = bot;                                                                                                      % il primo indice diventa l'ultimo dell'iterazione precedente     
    end

    % Dealloca tutte le strutture superflue
    clear tmp;
    for i = 1:nofmats                                                                                      
        eval(strcat('clear tout',num2str(i),';'));                                                                       % Libera l'area di memoria occupata da tutti gli array tout_1, tout_2, ecc.
    end

    counter=0;
    cont=sum(sum(abs(res-resc)));                                                                                      % Determina in valore assoluto la differenza tra la soluzione precedente e quella corrente
    ttt=sprintf('\nNumero di iterazione : %.0f - Errore : %f\n', iter, cont);                                                 
    if cont < limit
      fprintf('\nNumero di iterazione : %.0f - Errore : %f <--- SISTEMA BILANCIATO\n', iter, cont);                                                  
      break;
    elseif cont > oldcont                                                                                              % Il sistema è andato in stallo 
      ttt=strcat(ttt,' ...oscillazione...');  
    end  
    disp(ttt);    
    oldcont=cont;                                                                                                      % Memorizza lo scarto  
  end  
end        % fine del ciclo del gradiente coniugato

fprintf('\n\n****************    FINE BILANCIAMENTO    *****************\n\n\n');                                                  
  
%% BLOCCO LOGICO 24  

%%%%%%%%%%%%%%%%%%   SCRITTURA OGGETTI BILANCIATI E TERMINE PROGRAMMA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ciclo di determinazione dei parametri moltiplicati per w, elemento per elemento
for i = 1:nofmats 
  tout=0;
  for j = 1:size(pointeq,2) 
    if pointeq(i,j) > 0                                                                                                           % se l'indice è diverso da zero   
 
      % funzione speakeasy sum che prevede anche il riciclo nel caso due vettore abbiano differente dimensione
      tout=speakeasy_recycling_sum_optimised(tout, eval(strcat('li',num2str(pointeq(i,j)),'.*signseq(i,j).*w',num2str(pointeq(i,j)))),pp);     % da quanto ho capito, per ogni elemento presente in mats, costruisce dei tout che contengono dei li * w calcolati nel
                                                                                                                                  % seguente modo: considera i parametri del gradiente coniugato solo per le righe interessate dagli oggetti per esempio, 
                                                                                                                                  % oggetto OFFERT si trova alla riga 1, allora moltiplica li1 * (il segno di offert nella riga 1) * w1. Alla fine cumula 
                                                                                                                                  % tutto il risultato in tout1 (che corrisponde al primo elemento in mats, es OFFERT) and so on...
    end
  end 
      
    eval(strcat('tout',num2str(i),'=tout;'));                                                                                     % valuta l'espressione, e memorizza i vari tout nelle variabili tout1, tout2, tout3, ecc...
end  

% Dealloca tutte le strutture superflue
clear tmp;
for i = 1:nofmats                                                                                      
  eval(strcat('tin',num2str(i),'=',variances_prefix,mats{i},'.*tout',num2str(i),';'));                                
  eval(strcat('clear tout',num2str(i),';'));                                                                                       % Libera l'area di memoria occupata da tutti gli array tout_1, tout_2, ecc.
end

temp=toc(tim);                                                                                                                     % Tempo di esecuzione
tmin=floor(temp/60);                                                                                                               % Determina i minuti
tsec=mod(temp,60);                                                                                                                 % Determina i secondi 

tout=0;
% ciclo di determinazione delle soluzioni
for i = 1:nofmats
    tout=eval(strcat(mats{i},'-tin',num2str(i)));                                                                                  % Fa la differenza tra la matrice originale e quello determinato, forse proprio come viene fatto alla fine da Byron: y=p-V %*% t(G) %*% lambda_i

    % al momento salva in csv
    switch read_format
        case 0
        case 1
            %Approccio più laborioso ma più rapido rispetto a DLMWRITE che
            %ha una struttura più completa di verifica
            FID = fopen(fullfile(output_dir,['/Q' mats{i} '.csv']),'w+');
            fprintf(FID,[repmat('%.6f;',1,size(tout,2)-1) '%.6f\n'],tout');
            fclose(FID);
%            dlmwrite(fullfile(output_dir,['/Q' mats{i} '.csv']),tout,'delimiter', csv_delimiter,'precision', '%.6f');
        case 2
    end
end

temp_g=toc(tim_g);                                                                                                                 % Tempo di esecuzione generale
tmin_g=floor(temp_g/60);                                                                                                           % Determina i minuti
tsec_g=mod(temp_g,60);                                                                                                             % Determina i secondi 

fprintf('\nNumero di iterazioni eseguite : %.0f\n', iter)                                                  
fprintf('Tempo impiegato dal ciclo del gradiente coniugato+preprocessamenti (no tempi lettura/scrittura): %.2f minuti e %.2f secondi \n', tmin, tsec);
fprintf('Tempo generale impiegato dalla procedura : %.2f minuti e %.2f secondi \n', tmin_g, tsec_g);

save bal_up_new

  
  
  