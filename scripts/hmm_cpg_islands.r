library(HMM)
library(readr)
library(stringr)

numSim <- 10
symbols <- c("A", "C", "G", "T")
# states <- c("+", "-")
states = c("a", "c", "t", "g", "A", "C", "T", "G")

A = matrix(c(0.305849, 0.189895, 0.279249, 0.194143, 0.01, 0.01, 0.01, 0.01,
             0.283797, 0.368467, 0.106746, 0.322126, 0.01, 0.01, 0.00462963, 0.01,
             0.231064, 0.250871, 0.372331, 0.223427, 0.01, 0.01, 0.01, 0.01,
             0.179291, 0.189024, 0.237404, 0.260304, 0.01, 0.01, 0.01, 0.01,
             0.01, 0.00087108, 0.000853971, 0.01, 0.18, 0.131953, 0.178571, 0.0964467,
             0.01, 0.00087108, 0.01, 0.01, 0.352308, 0.36934, 0.327381, 0.423012,
             0.01, 0.01, 0.000853971, 0.01, 0.392308, 0.338939, 0.36045, 0.316413,
             0.01, 0.01, 0.00256191, 0.01, 0.0753846, 0.159767, 0.128968, 0.164129), ncol=8, byrow=TRUE)
E = matrix(c(1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1, 0,
             0, 0, 0, 1,
             1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1, 0,
             0, 0, 0, 1), ncol=4, byrow=TRUE)

hmm = initHMM(states, symbols, transProbs = A, emissionProbs = E)
print(hmm)

obs_true = read_file("D:\\FER\\BIOINF2\\data\\three-even-islands\\sequences\\chr19_0.txt")
obs_true = unlist(str_extract_all(obs_true, boundary("character")))
obs_true = toupper(obs_true)
obs_true = obs_true[!obs_true == '\r\n']

obs_pred = read_file("D:\\FER\\BIOINF2\\data\\three-even-islands\\sequences\\chr19_3.txt")
obs_pred = unlist(str_extract_all(obs_pred, boundary("character")))
obs_pred = toupper(obs_pred)
obs_pred = obs_pred[!obs_pred == '\r\n']

#if learning
bw = baumWelch(hmm, obs_true, 1)
print(bw$hmm)

viterb1 = viterbi(bw$hmm, obs_pred)
print(viterb1)

fileConn<-file("chr19_3_viterbi.txt")
writeLines(paste(viterb1,collapse=''), fileConn)
close(fileConn)