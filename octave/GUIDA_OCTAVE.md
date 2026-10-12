# Guida passo passo: eseguire il codice originale di TreeQSM con GNU Octave (Windows)

Obiettivo: far girare gli stessi file .m di TreeQSM 2.4.1 sulle nostre nuvole e confrontare i risultati con quelli del porting in R.
Non si programma nulla: si copiano alcuni comandi nella finestra di Octave, uno alla volta.

## 1. Installare Octave (una volta)
1. Aprire https://octave.org, "Download", sezione Windows, e scaricare l'installer a 64 bit (file .exe). Eseguirlo con le impostazioni predefinite.
2. Avviare "Octave" dal menu Start (la versione con interfaccia grafica). Si apre una finestra con la **Command Window** (la finestra dei comandi) al centro.
3. Digitare `version` e premere Invio: deve comparire il numero di versione. Poi `pkg list`: deve comparire un elenco di pacchetti; cercare `statistics`.
   Se manca: `pkg install -forge statistics` (serve internet; puo' chiedere di installare altri pacchetti da cui dipende).

## 2. Preparare la cartella di lavoro
1. Creare la cartella `D:\octave_treeqsm` (o un'altra, senza accenti).
2. Decomprimere `TreeQSM-master.zip` (da GitHub) dentro: deve risultare `D:\octave_treeqsm\TreeQSM-master\src`.
3. Decomprimere il pacchetto `pacchetto_matlab_treeqsm.zip` e copiare nella stessa cartella `D:\octave_treeqsm` questi file e cartelle:
   `esegui_treeqsm_octave.m`, `rif_num_v2.csv`, `rif_num_las11.csv`, `alberi_corretti_v2`, `alberi_corretti_las11`.

## 3. Dire a Octave dove sono i file
**Si digita solo quello che sta nelle caselle grigie, una riga alla volta, premendo Invio dopo ciascuna. Il testo "Deve rispondere" descrive cosa comparira', NON va digitato.**

Digita:  `cd D:/octave_treeqsm`   (nessuna risposta: e' normale)

Digita:  `addpath(genpath('D:/octave_treeqsm/TreeQSM-master/src'))`   (nessuna risposta)

Digita:  `which treeqsm`
Deve rispondere con una frase che contiene un percorso dentro `TreeQSM-master\src`, per esempio `'treeqsm' is a function from the file D:\octave_treeqsm\TreeQSM-master\src\treeqsm.m`.

Digita:  `which pdist2`
Deve rispondere con una frase che contiene un percorso (un file). Se risponde che non lo trova, digitare `pkg load statistics` e ripetere.

Digita:  `pwd`
Deve rispondere con la cartella di lavoro, per esempio `ans = D:/octave_treeqsm`.

Digita:  `ls`
Deve elencare `esegui_treeqsm_octave.m`, `rif_num_v2.csv`, `rif_num_las11.csv` e le cartelle `alberi_corretti_v2` e `alberi_corretti_las11`.

## 4. Prova minima: un solo albero (impostazioni predefinite, serve a misurare quanto e' lento)
```
esegui_treeqsm_octave
```
Octave sembrera' "bloccato": e' normale, sta lavorando. Alla fine dell'albero compare una riga del tipo `alberi_corretti_v2: albero 16 fatto (X minuti dall'inizio)`, poi un riepilogo.
**Annotare X**: ci dice se la corsa completa si puo' fare e in quanto tempo. Per interrompere: Ctrl+C nella Command Window.

## 5. Se la prova minima funziona: tre alberi, poi la corsa completa
```
alberi_prova = [4 16 35];
esegui_treeqsm_octave
```
Corsa completa (di notte, con il computer che non va in sospensione):
```
alberi_prova = [];
nuvole_da_fare = [1 2];
esegui_treeqsm_octave
```

## 6. Cosa mandare indietro
* il testo o uno screenshot di quanto compare in Command Window (anche se c'e' un errore, **il messaggio d'errore completo**);
* i file `risultati_octave_alberi_corretti_v2.csv` (e `risultati_octave_alberi_corretti_las11.csv`), che compaiono nella cartella di lavoro.

## Nota: la funzione alphaShape
Octave non ha la funzione `alphaShape`, che TreeQSM usa solo per due misure della chioma (area e volume "alpha" della chioma). Lo script `esegui_treeqsm_octave.m` contiene un sostituto minimo che restituisce NaN per quelle due misure: non toccano diametri, cilindri, altezza e volumi del fusto. I file originali di TreeQSM non sono modificati.
Se si usa una versione dello script senza il sostituto, prima di lanciarla digitare nella Command Window (tre righe, Invio dopo ciascuna):
```
function shp = alphaShape(varargin)
  shp.area = NaN; shp.volume = NaN;
end
```

## Problemi frequenti
* `'treeqsm' undefined`: il percorso non e' stato aggiunto: ripetere `addpath(genpath(...))` del punto 3.
* `'pdist2' undefined`: manca il pacchetto statistics: `pkg install -forge statistics`, poi `pkg load statistics`.
* `csvread: unable to open file`: Octave e' nella cartella sbagliata: controllare `pwd` e rifare `cd`.
* Altri errori: mandare il messaggio completo; di solito indica la funzione che Octave non capisce.
