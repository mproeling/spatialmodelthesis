
missing.persons = as.data.frame(seq(1, 1899))
colnames(missing.persons) = c("id")

missing.persons10 = as.data.frame(cbind(unique(missing.sample.10$sampleN), 1))
missing.persons25 = as.data.frame(cbind(unique(missing.sample.25$sampleN), 1))
missing.persons35 = as.data.frame(cbind(unique(missing.sample.35$sampleN), 1))
missing.persons50 = as.data.frame(cbind(unique(missing.sample.50$sampleN), 1))

colnames(missing.persons10) = c("id", "missing")
colnames(missing.persons25) = c("id", "missing")
colnames(missing.persons35) = c("id", "missing")
colnames(missing.persons50) = c("id", "missing")

missing.persons.all = merge(missing.persons, missing.persons10, by = "id", all = TRUE)
missing.persons.all = merge(missing.persons.all, missing.persons25, by = "id", all = TRUE)
missing.persons.all = merge(missing.persons.all, missing.persons35, by = "id", all = TRUE)
missing.persons.all = merge(missing.persons.all, missing.persons50, by = "id", all = TRUE)
colnames(missing.persons.all) = c("id", "missing10", "missing25", "missing35", "missing50")
missing.persons.all[is.na(missing.persons.all$missing10), 2] <- 0
missing.persons.all[is.na(missing.persons.all$missing25), 3] <- 0
missing.persons.all[is.na(missing.persons.all$missing35), 4] <- 0
missing.persons.all[is.na(missing.persons.all$missing50), 5] <- 0

write.table(missing.persons.all, "missing.persons.all.csv", quote = F, sep = ",", col.names = T, row.names = F)
