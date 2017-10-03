require(haven)
require(Rcpp)
require(circular)

sourceCpp(file = "VenterMode.cpp") 
sourceCpp(file = "VenterModeCircular.cpp") 

source(file = "CircularMixedEffectsFunctionsII.R")

Exp1 <- read_spss(file = "Cognitive Maps/WarrenExp1.sav")

Exp1$Error_rad <- (Exp1$Error/180)*pi
Exp1$Learning_c <- Exp1$Learning-mean(Exp1$Learning)
Exp1$MazeTrialInt <- Exp1$Maze*Exp1$Trial_type
Exp1$MazeLearning <- Exp1$Maze*Exp1$Learning_c
Exp1$TrialLearning <- Exp1$Trial_type*Exp1$Learning_c
Exp1$TrialMazeLearning <- Exp1$Trial_type * Exp1$Maze * Exp1$Learning_c

plot(circular(Exp1$Error_rad[Exp1$Maze == 0]))
plot(circular(Exp1$Error_rad[Exp1$Maze == 1]))
plot(circular(Exp1$Error_rad[Exp1$Trial_type == 0]))
plot(circular(Exp1$Error_rad[Exp1$Trial_type == 1]))
plot(circular(Exp1$Error_rad[Exp1$Maze == 0 & Exp1$Trial_type == 0]))
plot(circular(Exp1$Error_rad[Exp1$Maze == 1 & Exp1$Trial_type == 0]))
plot(circular(Exp1$Error_rad[Exp1$Maze == 0 & Exp1$Trial_type == 1]))
plot(circular(Exp1$Error_rad[Exp1$Maze == 1 & Exp1$Trial_type == 1]))

plot(circular(Exp1$Error_rad[Learning_c < 0]))
plot(circular(Exp1$Error_rad[Learning_c > 0]))

plot(circular(Exp1$Error_rad[Exp1$Maze == 0 & Exp1$Trial_type == 0 & Learning_c < 0]))
plot(circular(Exp1$Error_rad[Exp1$Maze == 1 & Exp1$Trial_type == 0 & Learning_c < 0]))
plot(circular(Exp1$Error_rad[Exp1$Maze == 0 & Exp1$Trial_type == 1 & Learning_c < 0]))
plot(circular(Exp1$Error_rad[Exp1$Maze == 1 & Exp1$Trial_type == 1 & Learning_c < 0]))

plot(circular(Exp1$Error_rad[Exp1$Maze == 0 & Exp1$Trial_type == 0 & Learning_c > 0]))
plot(circular(Exp1$Error_rad[Exp1$Maze == 1 & Exp1$Trial_type == 0 & Learning_c > 0]))
plot(circular(Exp1$Error_rad[Exp1$Maze == 0 & Exp1$Trial_type == 1 & Learning_c > 0]))
plot(circular(Exp1$Error_rad[Exp1$Maze == 1 & Exp1$Trial_type == 1 & Learning_c > 0]))

head(Exp1)

OutputIO <- me_model(Exp1, Error_rad, Subject, Trial_no,
                   ~1, ~1,
                   ~1, ~1,
                   its = 10000, burn = 1, n.lag = 1, priors = "default")

Output1 <- me_model(Exp1, Error_rad, Subject, Trial_no,
                   ~1 + Maze + Trial_type, ~1,
                   ~1 + Maze + Trial_type, ~1,
                   its = 10000, burn = 1, n.lag = 1, priors = "default")

Output1a <- me_model(Exp1, Error_rad, Subject, Trial_no,
                     ~1 + Learning_c, ~1,
                     ~1 + Learning_c, ~1,
                     its = 10000, burn = 1, n.lag = 1, priors = "default")

Output1b <- me_model(Exp1, Error_rad, Subject, Trial_no,
                     ~1 + Maze + Trial_type + Learning_c + MazeLearning, ~1,
                     ~1 + Maze + Trial_type + Learning_c + MazeLearning, ~1,
                     its = 10000, burn = 1, n.lag = 1, priors = "default")
  
Output1c <- me_model(Exp1, Error_rad, Subject, Trial_no,
                     ~1 + Maze + Trial_type + Learning_c + TrialLearning, ~1,
                     ~1 + Maze + Trial_type + Learning_c + TrialLearning, ~1,
                     its = 10000, burn = 1, n.lag = 1, priors = "default") 

Output2 <- me_model(Exp1, Error_rad, Subject, Trial_no,
                   ~1 + Maze + Trial_type + MazeTrialInt, ~1,
                   ~1 + Maze + Trial_type + MazeTrialInt, ~1,
                   its = 10000, burn = 1, n.lag = 1, priors = "default")

Output3 <- me_model(Exp1, Error_rad, Subject, Trial_no,
                    ~1 + Maze + Trial_type + MazeTrialInt + Learning_c, ~1,
                    ~1 + Maze + Trial_type + MazeTrialInt + Learning_c, ~1,
                    its = 10000, burn = 1, n.lag = 1, priors = "default")

Output4 <- me_model(Exp1, Error_rad, Subject, Trial_no,
                     ~1 + Maze + Trial_type + MazeTrialInt + TrialLearning + MazeLearning + TrialMazeLearning + Learning, ~1,
                     ~1 + Maze + Trial_type + MazeTrialInt + TrialLearning + MazeLearning + TrialMazeLearning+ Learning, ~1,
                     its = 10000, burn = 1, n.lag = 1, priors = "default")

OutputIO$summary
Output1$summary
Output1a$summary
Output1b$summary
Output1c$summary
Output2$summary
Output3$summary
Output4$summary

res <- Output2$output

Euclidean_standard <- atan2(res$Beta.II[101:10000, 1], res$Beta.I[101:10000, 1])

Euclidean_probe <- atan2(res$Beta.II[101:10000, 1] + res$Beta.II[101:10000, 3],
                         res$Beta.I[101:10000, 1] + res$Beta.I[101:10000, 3])

Wormhole_standard <- atan2(res$Beta.II[101:10000, 1] + res$Beta.II[101:10000, 2],
                           res$Beta.I[101:10000, 1] + res$Beta.I[101:10000, 2])

Wormhole_probe <- atan2(res$Beta.II[101:10000, 1] + res$Beta.II[101:10000, 3] + res$Beta.II[101:10000, 4],
                         res$Beta.I[101:10000, 1] + res$Beta.I[101:10000, 3] + res$Beta.I[101:10000, 4])

modeC(Euclidean_standard)
modeC(Euclidean_probe)
modeC(Wormhole_standard)
modeC(Wormhole_probe)

hpdC(Euclidean_standard)
hpdC(Euclidean_probe)
hpdC(Wormhole_standard)
hpdC(Wormhole_probe)

