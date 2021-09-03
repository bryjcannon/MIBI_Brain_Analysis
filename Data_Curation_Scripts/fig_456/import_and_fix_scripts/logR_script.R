# sample (sample_n) for later use by object type
n_sub_fraction <- 0.7

#tau_datsets
train_mglia_T <- sample_frac(filter(disease_cell_overlap, FlowSOM_ids == 'microglia' & tau_obj_POS == 0), n_sub_fraction)
train_mglia_T <- rbind(train_mglia_T, sample_frac(filter(disease_cell_overlap, FlowSOM_ids == 'microglia' & tau_obj_POS == 1), n_sub_fraction))
test_mglia_T <- setdiff(filter(disease_cell_overlap, FlowSOM_ids == 'microglia'), train_mglia_T)

train_neuron_T <- sample_frac(filter(disease_cell_overlap, FlowSOM_ids == 'neurons' & tau_obj_POS == 0), n_sub_fraction)
train_neuron_T <- rbind(train_neuron_T, sample_frac(filter(disease_cell_overlap, FlowSOM_ids == 'neurons' & tau_obj_POS == 1), n_sub_fraction))
test_neuron_T <- setdiff(filter(disease_cell_overlap, FlowSOM_ids == 'neurons'), train_neuron_T)

#amyloid_datsets
train_mglia_A <- sample_frac(filter(disease_cell_overlap, FlowSOM_ids == 'microglia' & amyloid_obj_POS == 0), n_sub_fraction)
train_mglia_A <- rbind(train_mglia_A, sample_frac(filter(disease_cell_overlap, FlowSOM_ids == 'microglia' & amyloid_obj_POS == 1), n_sub_fraction))
test_mglia_A <- setdiff(filter(disease_cell_overlap, FlowSOM_ids == 'microglia'), train_mglia_A)

train_neuron_A <- sample_frac(filter(disease_cell_overlap, FlowSOM_ids == 'neurons' & amyloid_obj_POS == 0), n_sub_fraction)
train_neuron_A <- rbind(train_neuron_A, sample_frac(filter(disease_cell_overlap, FlowSOM_ids == 'neurons' & amyloid_obj_POS == 1), n_sub_fraction))
test_neuron_A <- setdiff(filter(disease_cell_overlap, FlowSOM_ids == 'neurons'), train_neuron_A)

# models
mglia_tau_logR.fits <- glm(formula = tau_obj_POS ~ Iba1 + CD45 + ApoE4 + CD33Lyo + CD47 + EEA1 + MFN2, family = 'binomial', data = train_mglia_T)
mglia_amyloid_logR.fits <- glm(formula = amyloid_obj_POS ~ Iba1 + CD45 + ApoE4 + CD33Lyo + CD47 + EEA1 + MFN2, family = 'binomial', data = train_mglia_A)
neuro_tau_logR.fits <- glm(formula = tau_obj_POS ~ Parvalbumin + PanGAD6567 + VGAT + VGLUT1  + VGLUT2 + PSD95 + Synaptophysin + MFN2 + CD47 + Presenilin1NTF +  X8OHGuano + PolyubiK63 + PolyubiK48 + MAP, family = 'binomial', data = train_neuron_T)
neuro_amyloid_logR.fits <- glm(formula = amyloid_obj_POS ~ Parvalbumin + PanGAD6567 + VGAT + VGLUT1  + VGLUT2 + PSD95 + Synaptophysin + MFN2 + CD47 + Presenilin1NTF +  X8OHGuano + PolyubiK63 + PolyubiK48 + MAP, family = 'binomial', data = train_neuron_A)

summary(neuro_tau_logR.fits)

# tests - probs
mglia_tau_logR.probs <- predict(mglia_tau_logR.fits, test_mglia_T, type='response')
mglia_amyloid_logR.probs <- predict(mglia_amyloid_logR.fits, test_mglia_A, type='response')
neuro_tau_logR.probs <- predict(neuro_tau_logR.fits, test_neuron_T, type='response')
neuro_amyloid_logR.probs <- predict(neuro_amyloid_logR.fits, test_neuron_A, type='response')

# tests - predict
mglia_tau_logR.pred <- rep("0", dim(test_mglia_T)[1])
mglia_amyloid_logR.pred <- rep("0", dim(test_mglia_A)[1])
neuro_tau_logR.pred <- rep("0", dim(test_neuron_T)[1])
neuro_amyloid_logR.pred <- rep("0", dim(test_neuron_A)[1])

mglia_tau_logR.pred[mglia_tau_logR.probs >.4]="1"
mglia_amyloid_logR.pred[mglia_amyloid_logR.probs >.45]="1"
neuro_tau_logR.pred[neuro_tau_logR.probs >.4]="1"
neuro_amyloid_logR.pred[neuro_amyloid_logR.probs >.4]="1"

# tests - predict
tau_obj_check <- test_neuron_T$tau_obj_POS
table(neuro_tau_logR.pred, tau_obj_check)
mean(neuro_tau_logR.pred == tau_obj_check)

amyloid_obj_check <- test_neuron_A$amyloid_obj_POS
table(neuro_amyloid_logR.pred, amyloid_obj_check)
mean(neuro_amyloid_logR.pred == amyloid_obj_check)

tau_obj_check <- test_mglia_T$tau_obj_POS
table(mglia_tau_logR.pred, tau_obj_check)
mean(mglia_tau_logR.pred == tau_obj_check)

amyloid_obj_check <- test_mglia_A$amyloid_obj_POS
table(mglia_amyloid_logR.pred, amyloid_obj_check)
mean(mglia_amyloid_logR.pred == amyloid_obj_check)


#pcs
pca_res <- prcomp(filter(disease_cell_overlap, FlowSOM_ids == 'microglia')[glia_struct3], scale. = TRUE)
> autoplot(pca_res)

