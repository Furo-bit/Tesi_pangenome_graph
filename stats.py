import pstats

# Carica i risultati del profiling
p = pstats.Stats('output.prof')

# Ordina e stampa le statistiche (per tempo cumulativo)
p.sort_stats('cumulative').print_stats(10)
