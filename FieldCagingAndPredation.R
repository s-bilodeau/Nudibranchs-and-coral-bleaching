#clear workspace
rm(list=ls())


########################## Fisher's Exact Tests for Caged/Uncaged Eggs and Nudibranchs ##########################

CageNudi <- matrix(c(6, 1, 98, 191), nrow = 2,
                   dimnames =
                     list(c("Caged", "No Cage"),
                          c("Nudi", "No Nudi"))) #includes corals with both nudibranchs and eggs

fisher.test(CageNudi) 

CageEggs <- matrix(c(7, 0, 97, 192), nrow = 2,
                   dimnames =
                     list(c("Caged", "No Cage"),
                          c("Eggs", "No Eggs"))) #includes corals with both nudibranchs and eggs

fisher.test(CageEggs) 



########################## Fisher's Exact Tests for Predation ##########################

NudiPred <- matrix(c(20, 0, 20, 0), nrow = 2,
                   dimnames =
                     list(c("Eaten", "Rejected"),
                          c("Nudi", "Control"))) #comparing palatability of nudibranchs and controls (tuna)

fisher.test(NudiPred) #p=1 means no difference between slugs and tuna

EggPred <- matrix(c(8, 0, 8, 0), nrow = 2,
                   dimnames =
                     list(c("Eaten", "Rejected"),
                          c("Eggs", "Control"))) #comparing palatability of eggs and controls (tuna)

fisher.test(EggPred) #p=1 means no difference between eggs and tuna

