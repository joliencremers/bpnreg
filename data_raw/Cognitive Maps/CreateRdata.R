Maps <- haven::read_spss(file = "data-raw/Cognitive Maps/WarrenExp1.sav")

Maps$Error.rad <- (Maps$Error/180)*pi
Maps$Learn.c <- Maps$Learn-mean(Maps$Learn)
Maps$Maze <- as.factor(Maps$Maze)
Maps$Trial.type <- as.factor(Maps$Trial.type)

devtools::use_data(Maps)

