data(galo, package = "homals")

galo <- as.data.frame(galo[, 1:4])

lab1 <- c("F", "M")
lab2 <- c("1", "2", "3", "4", "5", "6", "7", "8", "9")
lab3 <- c("Agr", "Ext", "Gen", "Grls", "Man", "None", "Uni")
lab4 <- c("LwWC", "MdWC", "Prof", "Shop", "Skil", "Unsk")

labs_galo <- c(lab1, lab2, lab3, lab4)
ncat_galo <- c(2, 9, 7, 6)