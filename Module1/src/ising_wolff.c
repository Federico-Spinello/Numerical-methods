/*
 * ising_wolff.c - Implementazione simulazione Ising 2D con algoritmo di Wolff
 *
 * Questo modulo implementa una singola simulazione Monte Carlo del modello di
 * Ising 2D utilizzando l'algoritmo di Wolff (cluster updates).
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../include/ising_wolff.h"
#include "../include/geometry.h"
#include "../include/random.h"

#define DIM 2  // dimensionalita del reticolo

/*
 * Calcola magnetizzazione per sito
 */
static double magn(int const * const restrict lattice, long int volume)
  {
  long int r, sum;

  sum=0;
  for(r=0; r<volume; r++)
     {
     sum+=lattice[r];
     }

  return (double) sum / (double) volume;
  }


/*
 * Calcola energia per sito
 */
static double energy(int const * const restrict lattice,
              long int const * const restrict nnp,
              long int volume)
  {
  long int r, sum;
  int i;

  sum=0;
  for(r=0; r<volume; r++)
     {
     for(i=0; i<DIM; i++)
        {
        sum+=-lattice[r]*lattice[nnp[dirgeo(r, i, volume)]];
        }
     }

  return (double) sum / (double) volume;
  }


/*
 * Costruzione cluster non ricorsiva (evita stack overflow)
 */
static void build_cluster_norec(int const * const restrict lattice,
                         long int r,
                         int * restrict occup,
                         long int * restrict pointtoocc,
                         long int * restrict clustersize,
                         long int const * const restrict nnp,
                         long int const * const restrict nnm,
                         long int volume,
                         double prob)
  {
  (void) r; // evita warnings
  int i;
  long int index, r1;
  long int oldcs, oldcsnew;

  oldcs=0; // valore iniziale, con *clustersize=1

  // se i primi vicini hanno lo stesso orientamento e non sono occupati
  // vengono aggiunti al cluster con probabilita prob

  while(*clustersize>oldcs) // ci sono siti aggiunti recentemente da controllare
       {
       oldcsnew=*clustersize;

       for(index=oldcs; index<oldcsnew; index++)
          {
          r1=pointtoocc[index];

          for(i=0; i<DIM; i++)
             {
             // avanti
             if(occup[nnp[dirgeo(r1, i, volume)]]==0 && lattice[r1]*lattice[nnp[dirgeo(r1, i, volume)]]==1)
               {
               if(myrand()<prob)
                 {
                 occup[nnp[dirgeo(r1, i, volume)]]=1;
                 pointtoocc[*clustersize]=nnp[dirgeo(r1, i, volume)];
                 (*clustersize)++;
                 }
               }

             // indietro
             if(occup[nnm[dirgeo(r1, i, volume)]]==0 && lattice[r1]*lattice[nnm[dirgeo(r1, i, volume)]]==1)
               {
               if(myrand()<prob)
                 {
                 occup[nnm[dirgeo(r1, i, volume)]]=1;
                 pointtoocc[*clustersize]=nnm[dirgeo(r1, i, volume)];
                 (*clustersize)++;
                 }
               }
             }
          }

       oldcs=oldcsnew;
       }
  }


/*
 * Esegue una singola simulazione Ising 2D con Wolff
 *
 * Questa funzione è il cuore del modulo: inizializza il reticolo,
 * esegue termalizzazione e misurazioni, e salva i risultati.
 */
int run_ising_simulation(const SimulationParams *params, SimulationResults *results)
    {
    int  L, *lattice, *occup;
    long int i, r, volume, thermalization, measurements, iter, clustersize;
    long int *nnp, *nnm, *pointtoocc;
    double beta, T, locE, locM, prob, m2_val, m4_val;
    double sumE, sumE2, sumM, sumM2, sumM4, sumM2sq, sumM4sq;
    double meanE, meanM, meanE2, meanM2, meanM4;
    double chi, C, binder;
    double varE, varM, varM2, varM4, stdE, stdM, stdM2, stdM4;
    double errE, errM, errChi, errC, errBinder;
    double *M_samples;  // array per salvare m di ogni measurement
    FILE *fp;

    const unsigned long int seed1=(unsigned long int) time(NULL);
    const unsigned long int seed2=seed1+127;

    // Estrai parametri dalla struttura
    L = params->L;
    T = params->T;
    thermalization = params->thermalization;
    measurements = params->measurements;

    // Validazione parametri
    if(L<=0)
      {
      fprintf(stderr, "ERRORE: L deve essere positivo\n");
      return -1;
      }

    if(T<=0.0)
      {
      fprintf(stderr, "ERRORE: T deve essere positiva\n");
      return -1;
      }

    if(thermalization<0)
      {
      fprintf(stderr, "ERRORE: thermalization deve essere non negativo\n");
      return -1;
      }

    if(measurements<=0)
      {
      fprintf(stderr, "ERRORE: measurements deve essere positivo\n");
      return -1;
      }

    // inizializza generatore numeri casuali
    myrand_init(seed1, seed2);

    // calcola il volume
    volume=1;
    for(i=0; i<DIM; i++)
       {
       volume*=L;
       }

    // alloca il reticolo e i primi vicini: nnp[dirgeo(r, i, volume)]= primo vicino nella direzione positiva "i" del sito r
    
    lattice=(int *)malloc((unsigned long int)(volume)*sizeof(int));
    if(lattice == NULL)
      {
      fprintf(stderr, "ERRORE: allocazione lattice fallita\n");
      return -1;
      }
    nnp=(long int *)malloc((unsigned long int)(DIM*volume)*sizeof(long int));
    if(nnp == NULL)
      {
      fprintf(stderr, "ERRORE: allocazione nnp fallita\n");
      free(lattice);
      return -1;
      }
    nnm=(long int *)malloc((unsigned long int)(DIM*volume)*sizeof(long int));
    if(nnm == NULL)
      {
      fprintf(stderr, "ERRORE: allocazione nnm fallita\n");
      free(lattice);
      free(nnp);
      return -1;
      }

    // array per tenere traccia dei siti occupati durante costruzione cluster
    // 0 = libero, 1 = occupato
    occup=(int *)malloc((unsigned long int)(volume)*sizeof(int));
    if(occup== NULL)
      {
      fprintf(stderr, "ERRORE: allocazione occup fallita\n");
      free(lattice);
      free(nnp);
      free(nnm);
      return -1;
      }

    // le prime "clustersize" entrate puntano ai siti del cluster
    pointtoocc=(long int *)malloc((unsigned long int)(volume)*sizeof(long int));
    if(pointtoocc== NULL)
      {
      fprintf(stderr, "ERRORE: allocazione pointtoocc fallita\n");
      free(lattice);
      free(nnp);
      free(nnm);
      free(occup);
      return -1;
      }

    // Alloca array solo per magnetizzazione (serve per calcolare var(m^2) e var(m^4))
    M_samples = (double *)malloc((unsigned long int)(measurements) * sizeof(double));
    if(M_samples == NULL)
      {
      fprintf(stderr, "ERRORE: allocazione array M_samples fallita\n");
      free(lattice);
      free(nnp);
      free(nnm);
      free(occup);
      free(pointtoocc);
      return -1;
      }

    // Inizializza accumulatori per calcolo statistiche
    sumE = 0.0;
    sumE2 = 0.0;
    sumM = 0.0;
    sumM2 = 0.0;
    sumM4 = 0.0;

    // inizializza nnp e nnm
    init_neighbors(nnp, nnm, L, DIM);

    // inizializza reticolo a configurazione ordinata (cold start per T<2.3, hot altrimenti)
    if(T < 2.3)
      {
      // cold start - tutti spin allineati
      for(r=0; r<volume; r++)
         {
         lattice[r]=1;
         }
      }
    else
      {
      // hot start - configurazione casuale
      for(r=0; r<volume; r++)
         {
         if(myrand() < 0.5)
           {
           lattice[r]=1;
           }
         else
           {
           lattice[r]=-1;
           }
         }
      }

    // beta = temperatura inversa
    beta = 1.0 / T;

    // probabilita di aggiunta al cluster
    prob=1.0-exp(-2.0*beta);

    // ============================================================
    // FASE DI TERMALIZZAZIONE
    // ============================================================
    for(iter=0; iter<thermalization; iter++)
       {
       // azzera array di occupazione
       for(r=0; r<volume; r++)
          {
          occup[r]=0;
          }
       clustersize=0;

       // scegli sito random come seed
       r=(long int)((double)volume*myrand());
       occup[r]=1; // r marcato come occupato
       pointtoocc[clustersize]=r; // puntatore a r aggiunto
       clustersize++;

       // costruisci cluster
       build_cluster_norec(lattice, r, occup, pointtoocc, &clustersize, nnp, nnm, volume, prob);

       // flippa il cluster
       for(r=0; r<clustersize; r++)
          {
          lattice[pointtoocc[r]]=-lattice[pointtoocc[r]];
          }
       }

    // ============================================================
    // FASE DI MISURAZIONE
    // ============================================================

    for(iter=0; iter<measurements; iter++)
       {
       // azzera array di occupazione
       for(r=0; r<volume; r++)
          {
          occup[r]=0;
          }
       clustersize=0;

       // scegli sito random come seed
       r=(long int)((double)volume*myrand());
       occup[r]=1;
       pointtoocc[clustersize]=r;
       clustersize++;

       // costruisci cluster
       build_cluster_norec(lattice, r, occup, pointtoocc, &clustersize, nnp, nnm, volume, prob);

       // flippa il cluster
       for(r=0; r<clustersize; r++)
          {
          lattice[pointtoocc[r]]=-lattice[pointtoocc[r]];
          }

       // calcola osservabili
       locE = energy(lattice, nnp, volume);
       locM = magn(lattice, volume);

       // Salva magnetizzazione per calcoli successivi
       M_samples[iter] = locM;

       // Accumula per tutti i momenti necessari
       sumE += locE;
       sumE2 += locE * locE;
       sumM += fabs(locM);  // <|m|>
       sumM2 += locM * locM;  // <m^2> (senza valore assoluto!)
       sumM4 += locM * locM * locM * locM;  // <m^4>
       }

    // ============================================================
    // CALCOLO OSSERVABILI FINALI
    // ============================================================

    // Calcola medie
    meanE = sumE / (double)measurements;
    meanM = sumM / (double)measurements;
    meanE2 = sumE2 / (double)measurements;
    meanM2 = sumM2 / (double)measurements;
    meanM4 = sumM4 / (double)measurements;

    // Varianze usando <x^2> - <x>^2 (formula standard per momenti)
    varE = meanE2 - meanE * meanE;
    varM = meanM2 - meanM * meanM;  // var(|m|) ≈ <m^2> - <|m|>^2

    // Per errore su Binder servono anche varianze di m^2 e m^4
    // Calcoliamo <(m^2)^2> e <(m^4)^2> dall'array M_samples
    sumM2sq = 0.0;  // somma di (m^2)^2
    sumM4sq = 0.0;  // somma di (m^4)^2

    for(iter=0; iter<measurements; iter++)
       {
       locM = M_samples[iter];
       m2_val = locM * locM;
       m4_val = m2_val * m2_val;

       sumM2sq += m2_val * m2_val;  // (m^2)^2
       sumM4sq += m4_val * m4_val;  // (m^4)^2
       }

    // Varianze di m^2 e m^4: var(X) = <X^2> - <X>^2
    varM2 = (sumM2sq / (double)measurements) - meanM2 * meanM2;
    varM4 = (sumM4sq / (double)measurements) - meanM4 * meanM4;

    // Deviazioni standard
    stdE = sqrt(varE);
    stdM = sqrt(varM);
    stdM2 = sqrt(varM2);
    stdM4 = sqrt(varM4);

    // Errori standard della media (SEM)
    errE = stdE / sqrt((double)measurements);
    errM = stdM / sqrt((double)measurements);

    // suscettivita magnetica: chi = beta * V * (<M^2> - <|M|>^2)
    chi = beta * (double)volume * (meanM2 - meanM * meanM);

    // calore specifico: C = beta^2 * V * (<E^2> - <E>^2)
    C = beta * beta * (double)volume * (meanE2 - meanE * meanE);

    // cumulante di Binder: U = <M^4>/<M^2>^2
    binder = meanM4 / (meanM2 * meanM2);

    // Errore su suscettivita (propagazione da fluttuazioni di m)
    errChi = beta * (double)volume * sqrt(2.0) * stdM / sqrt((double)measurements);

    // Errore su calore specifico (propagazione da fluttuazioni di E)
    errC = beta * beta * (double)volume * sqrt(2.0) * stdE / sqrt((double)measurements);

    // Errore su Binder cumulant: U = <m^4> / <m^2>^2
    // Propagazione errori: δU/U = sqrt[ (δm4/m4)^2 + 4(δm2/m2)^2 ]
    // dove δm4 = σ(m^4)/sqrt(N) e δm2 = σ(m^2)/sqrt(N)
    double errM2 = stdM2 / sqrt((double)measurements);  // SEM di <m^2>
    double errM4 = stdM4 / sqrt((double)measurements);  // SEM di <m^4>

    if(meanM2 > 1e-10 && meanM4 > 1e-10)
      {
      // Errore relativo
      double rel_err_m2 = errM2 / meanM2;
      double rel_err_m4 = errM4 / meanM4;
      double rel_err_binder = sqrt(rel_err_m4*rel_err_m4 + 4.0*rel_err_m2*rel_err_m2);
      errBinder = binder * rel_err_binder;
      }
    else
      {
      errBinder = 0.0;  // evita divisione per zero
      }

    // ============================================================
    // SALVA RISULTATI
    // ============================================================

    fp=fopen(params->output_file, "w");
    if(fp==NULL)
      {
      fprintf(stderr, "ERRORE: apertura file %s fallita\n", params->output_file);
      free(lattice);
      free(occup);
      free(pointtoocc);
      free(nnp);
      free(nnm);
      free(M_samples);
      return -1;
      }

    // scrivi header
    fprintf(fp, "# Ising 2D - Algoritmo di Wolff\n");
    fprintf(fp, "# L = %d, T = %.6f, beta = %.6f\n", L, T, beta);
    fprintf(fp, "# thermalization = %ld, measurements = %ld\n", thermalization, measurements);
    fprintf(fp, "# Colonne: T  m  err_m  E  err_E  chi  err_chi  C  err_C  binder  err_binder\n");

    // scrivi dati (11 colonne ora)
    fprintf(fp, "%.6f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n",
            T, meanM, errM, meanE, errE, chi, errChi, C, errC, binder, errBinder);

    fclose(fp);

    // Popola struttura risultati (se fornita)
    if(results != NULL)
      {
      results->mag_mean = meanM;
      results->mag2_mean = meanM2;
      results->mag4_mean = meanM4;
      results->energy_mean = meanE;
      results->energy2_mean = meanE2;
      results->mag_err = errM;
      results->energy_err = errE;
      results->chi_err = errChi;
      results->C_err = errC;
      }

    // libera memoria
    free(lattice);
    free(occup);
    free(pointtoocc);
    free(nnp);
    free(nnm);
    free(M_samples);

    return 0;
    }
