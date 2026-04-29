# Elementi Finiti 2025-2026

Benvenuti nella repository ufficiale del corso di **Elementi Finiti 2025-2026**. Qui troverete tutto il materiale necessario per le esercitazioni, inclusi script, notebook e documentazione.

## 📚 Informazioni sul corso
- **Docente**: [Prof. Giancarlo Sangalli](https://www-dimat.unipv.it/sangalli/)
- **Tutor**: [Ivan Bioli](https://sites.google.com/view/ivan-bioli/homepage)
- **Pagina ufficiale del corso**: [Elementi Finiti - UNIPV](https://www-dimat.unipv.it/sangalli/elementi_finiti_mat.html)

## 📌 Contenuto della repository
Questa repository contiene:
- **Materiale didattico**: appunti, slide e riferimenti utili.
- **Esercitazioni**: file Julia con esempi e problemi da svolgere.
- **Ambiente di lavoro**: configurazioni e pacchetti necessari per eseguire il codice.

## 🚀 Installazione e Setup
### 1️⃣ Installare Julia
Scaricate e installate Julia seguendo le istruzioni ufficiali:
[Download Julia](https://julialang.org/downloads/)

### 2️⃣ Installare VSCode (Opzionale ma consigliato)
Scaricate e installate VSCode:
[Download VSCode](https://code.visualstudio.com/)

Poi installate l'estensione ufficiale per Julia tramite il Marketplace di VSCode.

### 3️⃣ Clonare la repository
Se avete `git` installato, potete clonare questa repository eseguendo il comando:
```bash
git clone https://github.com/IvanBioli/elementifiniti2526.git
```
Se non avete `git`, potete scaricare il contenuto come file ZIP e decomprimerlo.

### 4️⃣ Installare i pacchetti necessari
Aprite ed eseguite il file `src_julia/exercises/ex00/ex00_packages.jl` in VSCode, oppure da terminale:
```bash
julia src_julia/exercises/ex00/ex00_packages.jl
```

Esistono due modalità di installazione — scegliete modificando la prima riga del file:

- **Pacchetti locali (ambiente di progetto):** lasciate attiva la riga
  ```julia
  Pkg.activate("elementifinitiunipv_pkg"; shared=false)
  ```
  I pacchetti vengono installati nell'ambiente isolato `elementifinitiunipv_pkg/`. Tenete questa riga **attiva** anche in tutti i file delle esercitazioni successive.

- **Pacchetti globali _(opzione consigliata)_:** commentate la riga sopra (aggiungete `#` davanti).
  I pacchetti vengono installati nell'ambiente Julia globale e sono disponibili senza attivare alcun progetto. Assicuratevi che la stessa riga sia **commentata** anche in tutti i file delle esercitazioni successive.

## 🎓 Primi Passi in Julia
Per seguire il corso con efficacia, è necessario completare il tutorial introduttivo su JuliaAcademy:
[Introduction to Julia (for programmers)](https://juliaacademy.com/courses).
Dopo aver completato il corso, riceverete un certificato ufficiale.

## ❓ Supporto
Se avete domande, contattate il docente o il tutor del corso. Potete anche aprire una issue in questa repository.

Buon lavoro!
