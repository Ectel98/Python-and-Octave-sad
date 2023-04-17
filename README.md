
Sleep Apnea Detection - Python and Octave 


Partendo da un file contenente i dati relativi ad un segnale ecg, mediante Octave vengono estratti gli intervalli R-R utilizzando un codice basato sull'algoritmo di Pan–Tompkins, viene poi generanto un file che li contente questi intervalli. 
Un altro algritmo in Python elabora questo file per determinare la presenza di fenomeni di Apnea notturna.

L'algoritmo per l'elaborazione degli intervalli R-R è basato su questa pubblicazione https://archive.physionet.org/physiotools/apdet/index.shtml
Per testare la correttezza di quest'algorimto sono stati utilizzati i file del seguente database https://physionet.org/content/apnea-ecg/1.0.0/
