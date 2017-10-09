Motor <- haven::read_spss(file = "data-raw/Motor Cognition/PhaseDifference.sav")
Motor <- Motor[,1:3]

Motor <- cbind(Motor, Motor$PhaseDiff*(pi/180))
colnames(Motor) <- c("Cond", colnames(Motor[,-c(1,4)]), "Phaserad")

Motor$Cond <- factor(Motor$Cond, labels = c("imp", "semi.imp", "exp"))
Motor$AvAmp <- Motor$AvAmp - mean(Motor$AvAmp)

devtools::use_data(Motor)


