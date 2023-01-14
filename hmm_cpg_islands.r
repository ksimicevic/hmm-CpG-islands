library(HMM)

numSim <- 1000
symbols <- c("A", "C", "G", "T")
states <- c("+", "-")

A = matrix(c(), byrow=TRUE)
E = matrix(c(), byrow=TRUE)

hmm = initHMM(states, symbols, transProbs = A, emissionProbs = E)
print(hmm)

sim1 = simHMM(hmm, numSim)
viterb1 = viterbi(hmm, sim1$observation)
print(viterb1)
