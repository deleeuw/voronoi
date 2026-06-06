data(galo, package = "homals")

galo <- galo[, 1:4]

ncat_galo <- apply(galo, 2, function(x) length(unique(x)))

lab1 <- c("F", "M")
lab2 <- c("1", "2", "3", "4", "5", "6", "7", "8", "9")
lab3 <- c("Agr", "Ext", "Gen", "Grls", "Man", "None", "Uni")
lab4 <- c("LowWC", "MidWC", "Prof", "Shop", "Skil", "Unsk")

labs_galo <- c(lab1, lab2, lab3, lab4)