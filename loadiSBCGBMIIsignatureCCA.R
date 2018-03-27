### This file loads SBC + CCA gene signature ####
loadiSBCGBMIIsignatureCCA = function(Y1, Y2, Y1.new, Y2.new){
################################################## CCA analysis ###################
Y1.first <- Y1
Y2.first <- Y2
res.cc <- cc(Y1.first,Y2.first)
plot(res.cc$cor,type="b")


Y1.first.test <- Y1.new
Y2.first.test <- Y2.new

############ Apparently 50 percent of the Variance is explained by the first two canonical axis
sum(res.cc$cor[1:10])/sum(res.cc$cor)

Y1 <- res.cc$scores$xscores[,1:10]
Y2 <- res.cc$scores$yscores[,1:10]


#################################################################################
res.cc.test <- comput(Y1.first.test, Y2.first.test, res.cc)
Y1.test <- res.cc.test$xscores[,1:10]
Y2.test <- res.cc.test$yscores[,1:10]

#################### SAVING THE DATA SETS ################################################
relev <- list('Y1' =Y1, 'Y2' = Y2, 'Y1.new' = Y1.test, 'Y2.new' = Y2.test)

}
