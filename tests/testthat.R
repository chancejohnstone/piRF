library(testthat)
library(piRF)

set.seed(2020)
airfoil <- airfoil[sample.int(n=nrow(airfoil), nrow(airfoil)), ]

output <- rfint(pressure~.,
                     train_data = airfoil[1:500, ],
                     test_data = airfoil[500:550, ],
                     method = c("Zhang","quantile", "HDI", "Roy", "Ghosal", "Romano", "Tung"),
                     alpha = 0.1,
                     num_threads = 1)

y <- airfoil$pressure[500:550]
mean(output$int$Zhang[,1] < y & output$int$Zhang[,2] > y)
mean(output$int$quantile[,1] < y & output$int$quantile[,2] > y)
mean(output$int$HDI[,1] < y & output$int$HDI[,2] > y)
mean(output$int$Roy[,1] < y & output$int$Roy[,2] > y)
mean(output$int$Ghosal[,1] < y & output$int$Ghosal[,2] > y)
mean(output$int$Romano[,1] < y & output$int$Romano[,2] > y)
mean(output$int$Tung[,1] < y & output$int$Tung[,2] > y)
