function  c=speakeasy_recycling_sum(a,b,pp)
% Questa funzione ritorna la somma riciclata tra due array come viene effettuata in speakeasy/modeleasy. 
% In pratica considera tutti i casi di array in cui speakeasy ha un comportamento differente da R e li rettifica. 
% Caso 1: se gli array sono vettori (stessa o anche di dimensioni diverse), se il primo � vettore colonna ed il secondo � vettore riga allora produce una matrice che ha come numero di righe, quello del vettore 
% colonna e come numero di colonne quello del vettore riga. Al suo interno la matrice risultante presenta ad ogni colonna, il vettore colonna sommato ad ogni elemento del vettore riga. Se il primo � un vettore 
% riga ed il secondo � un vettore colonna, il risultato � identico. 
% Caso 2: se il primo operando � una matrice nxm e il secondo � un vettore riga 1xm, allora il risultato sar� una matrice nxm, dove ogni colonna del primo operando viene sommata al secondo operando. Questa operazione
% pu� essere commutativa. 
% Caso 3: se il primo operando � una matrice nxm e il secondo � un vettore colonna nx1, allora il risultato sar� una matrice nxm, dove ogni colonna del primo operando viene sommata al secondo operando. Questa operazione
% pu� essree commutativa.

% In pratica speakeasy � in grado di "spalmare" un vettore colonna o una matrice su di un vettore riga che abbia come numero di colonne uguale al numero di colonne del primo operando  
% Di seguito vengono contemplati tutti i casi di somma tra due array che speakeasy/modeleasy � in grado di gestire di default ed R no. 

c = bsxfun(pp,a,b);