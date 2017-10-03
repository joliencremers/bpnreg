Motor <- haven::read_spss(file = "data-raw/Motor Cognition/PhaseDifference.sav")
Motor <- Motor[,1:3]

Motor <- cbind(Motor, Motor$PhaseDiff*(pi/180))
colnames(Motor) <- c(colnames(Motor[,-4]), "Phaserad")

Motor$Condition <- factor(Motor$Condition, labels = c("implicit", "semi.implicit", "explicit"))
Motor$AvAmp <- Motor$AvAmp - mean(Motor$AvAmp)

devtools::use_data(Motor)


