library(HMM)
library(readr)
library(stringr)


states = c("+", "-")
symbols = c("A", "C", "G", "T")

# Ovo su vrijednosti dobivene iz C++ koda
A = matrix(c(0.998544, 0.00146056, 0.0014556, 0.998539), ncol=2, byrow=TRUE)
E = matrix(c(0.15769, 0.359534, 0.333819, 0.148957,
             0.23309, 0.254988, 0.235036, 0.276886), ncol=4, byrow=TRUE)

hmm = initHMM(states, symbols, transProbs = A, emissionProbs = E)
print(hmm)

# Ovaj file predstavlja primjer na kojem se uci
# Ovo ispod su samo gluposti specificne za R
observations = read_file("chr19_0.txt")
observations = unlist(str_extract_all(observations, boundary("character")))
observations = toupper(observations)
observations = observations[!observations == '\r\n']

# Ovaj file predstavlja primjer na kojem pokrecemo algoritam
test_observations = read_file("chr19_3.txt")
test_observations = unlist(str_extract_all(test_observations, boundary("character")))
test_observations = toupper(test_observations)
test_observations = test_observations[!test_observations == '\r\n']

# BW ucenje
bw = baumWelch(hmm, observations, 300)
print(bw$hmm)

# Viterbi predikcije
viterbi1 = viterbi(bw$hmm, test_observations)
print(viterbi1)

# Spremanje u file (ime file-a proizvoljno)
fileConn<-file("chr19_3_viterbi.txt")
writeLines(paste(viterbi1,collapse=''), fileConn)
close(fileConn)
 