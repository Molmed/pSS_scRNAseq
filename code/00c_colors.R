#00b_pheno.R

####----- create color objects for plotting
# create named color vector
colors_cellType <- c(naive= rev(brewer.pal(7, "Greens")),
                     dn = brewer.pal(5, "Blues"),
                     memory = brewer.pal(6, "YlOrRd"),
                     memory_stressed = c("gray80", "gray60", "gray40"),
                     memory_platelet = brewer.pal(3, "RdPu")[2],
                     plasma_cells = brewer.pal(3, "Purples")[3])

#show_col(colors_cellType)

# clusters

# 23 cell clusters:
# Naive 7 - green
# DN 5 - blue
# Memory_IgM - ?
# Memory_Classical - ?
# Memory_Stressed 3 - gray
# Memory_Platelet 1 - pink
# Plasma_cells 1 - purple

colors_clusters <- colors_cellType[c("naive1", "naive2", "naive3",
                                     "memory1", "naive4", "memory2",
                                     "naive4", "dn1", "memory_stressed1",
                                     "dn2", "naive5", "dn3",
                                     "memory3", "dn4", "memory4",
                                     "dn5", "memory5", "memory4",
                                     "memory_stressed2", "memory_platelet", "naive6",
                                     "memory_stressed3", "plasma_cells")]

# cellTypeFine (16)
colors_cellTypeFine <- colors_cellType[c("memory5",
                                         "dn1", "dn2", "dn3", "dn4", "dn5",
                                         "memory1", "memory3", "memory4",
                                         "memory_platelet",
                                         "memory_stressed1", "memory_stressed2",
                                         "naive2", "naive4", "naive6",
                                         "plasma_cells")]

# cellTypeMain (7)
colors_cellTypeMain <- colors_cellType[c("memory5", "dn3", "memory3", "memory_platelet", "memory_stressed1",
                                         "naive4", "plasma_cells", "naive6")]

# cellType scpred
colors_scpred1 <- colors_cellType[c("memory5", "dn3", "memory3",
                                    "naive2", "naive4", "memory_stressed1")] #6

colors_scpred2 <- colors_cellType[c("memory5", "memory4",
                                    "dn2", "dn3", "dn4", "dn5",
                                    "memory2", "memory1",
                                    "naive2", "naive4", "memory_stressed1")] #11
