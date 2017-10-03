Maps <- haven::read_spss(file = "data-raw/Cognitive Maps/WarrenExp1.sav")

Maps$Error.rad <- (Maps$Error/180)*pi
Maps$Learning.c <- Maps$Learning-mean(Maps$Learning)
Maps$Maze <- as.factor(Maps$Maze)
Maps$Trial.type <- as.factor(Maps$Trial.type)

devtools::use_data(Maps)

