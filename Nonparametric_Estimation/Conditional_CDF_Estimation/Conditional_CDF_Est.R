rm(list=ls())

library(MASS);library(plot3D)

unsmooth_CDF = function(grid_X,grid_Y,sample_X,sample_Y,h_X){
  F_hat = matrix(0,length(grid_X),length(grid_Y))
  for (j in 1:length(grid_Y)){
    y_idx = (sample_Y <= grid_Y[j])+0
    for (i in 1:length(grid_X)){
      F_hat[i,j] = (mean(dnorm((sample_X-grid_X[i])/h_X)*y_idx)/mean(dnorm((sample_X-grid_X[i])/h_X)))
    }
  }
  return(F_hat)
}

unsmooth_CDF_pred = function(grid_X,grid_Y,sample_X,sample_Y,h_X){
  for (j in 1:length(grid_Y)){
    y_idx = (sample_Y <= grid_Y[j])+0
    for (i in 1:length(grid_X)){
      F_hat = (mean(dnorm((sample_X-grid_X[i])/h_X)*y_idx)/mean(dnorm((sample_X-grid_X[i])/h_X)))
    }
  }
  return(F_hat)
}

smooth_CDF = function(grid_X,grid_Y,sample_X,sample_Y,h_X,h_Y){
  F_hat = matrix(0,length(grid_X),length(grid_Y))
  for (j in 1:length(grid_Y)){
    p = pnorm(grid_Y[j],sample_Y,h_Y)
    for (i in 1:length(grid_X)){
      F_hat[i,j] = (mean(dnorm((sample_X-grid_X[i])/h_X)*p)/mean(dnorm((sample_X-grid_X[i])/h_X)))
    }
  }
  return(F_hat)
}

smooth_CDF_pred = function(grid_X,grid_Y,sample_X,sample_Y,h_X,h_Y){
  F_hat = matrix(0,length(grid_X),length(grid_Y))
  for (j in 1:length(grid_Y)){
    p = pnorm(grid_Y[j],sample_Y,h_Y)
    for (i in 1:length(grid_X)){
      F_hat = (mean(dnorm((sample_X-grid_X[i])/h_X)*p)/mean(dnorm((sample_X-grid_X[i])/h_X)))
    }
  }
  return(F_hat)
}


plot(a_X,a_Y, xlab="X",ylab="Y", main="Scatter plot of X,Y")
plot(b_X,b_Y, xlab="X",ylab="Y", main="Scatter plot of X,Y")

#------------------------------ Data Generate
datasize = 200
a_mu = c(0,0)
a_Sigma = as.matrix(rbind(c(1,0.8),c(0.8,1)))
a_XY = mvrnorm(datasize, mu=a_mu, Sigma=a_Sigma)
a_X = a_XY[,1]
a_Y = a_XY[,2]

b_X = rnorm(datasize,mean=0,sd=1)
b_eps = rnorm(datasize,mean=0,sd=1)
b_Y = sin(pi*b_X) + exp(-b_X^2)*b_eps


#------------------------------ problem a 
# Ready
grid_a_X = seq(min(a_X),max(a_X),0.1)
grid_a_Y = seq(min(a_Y),max(a_Y),0.1)

h_X = 0.36

F_hat = unsmooth_CDF(grid_a_X,grid_a_Y,a_X,a_Y,h_X)

h_X = 0.21
h_Y = 0.33
F_hat2 = smooth_CDF(grid_a_X,grid_a_Y,a_X,a_Y,h_X,h_Y)

persp(grid_a_X,grid_a_Y,F_hat, theta = 0, phi = 20, expand = 0.5, col = "lightblue", xlab="X",ylab="Y",zlab="F_hat",
      main="Unsmooth Estimation (h_x=0.36)")
persp(grid_a_X,grid_a_Y,F_hat2, theta = 0, phi = 20, expand = 0.5, col = "lightblue", xlab="X",ylab="Y",zlab="F_hat",
      main="Smooth Estimation (h_x=0.21, h_y=0.33)")


h_vec = seq(0.01,0.5,0.01)
val3.vec=NULL
for (h in 1:length(h_vec)){
  print(paste("h: ", h))
  val2.vec=NULL
  for (y in 1:length(grid_a_Y)){
    val.vec=NULL
    for (i in 1:datasize){
      val = ((a_Y[i] <= grid_a_Y[y])+0 - unsmooth_CDF_pred(a_X[i],grid_a_Y[y],a_X[-i],a_Y[-i],h_vec[h]))^2
      val.vec = c(val.vec,val)
    }
    val2 = mean(val.vec)
    val2.vec = c(val2.vec,val2)
  }
  val3 = mean(val2.vec)
  val3.vec = c(val3.vec,val3)
}

plot(h_vec,val3.vec, xlab = "h", ylab = "CV(h)", main="Cross-validation curve (optimal h = 0.36)")
opt = which.min(val3.vec[-1])
h_vec[opt+1] # opt_h_x = 0.36

h_vec = seq(0.01,0.5,0.01)
h_df = expand.grid(h_vec,h_vec)
val3.vec=NULL
for (h_x in 1:length(h_vec)){
  print(paste("h_x: ", h_x))
  for (h_y in 1:length(h_vec)){
    print(paste("h_y: ", h_y))
    val2.vec=NULL
    for (y in 1:length(grid_a_Y)){
      val.vec=NULL
      for (i in 1:datasize){
        val = ((a_Y[i] <= grid_a_Y[y])+0 - smooth_CDF_pred(a_X[i],grid_a_Y[y],a_X[-i],a_Y[-i],h_vec[h_x],h_vec[h_y]))^2
        val.vec = c(val.vec,val)
      }
      val2 = mean(val.vec)
      val2.vec = c(val2.vec,val2)
    }
    val3 = mean(val2.vec)
    val3.vec = c(val3.vec,val3)
  }
}

cbind(h_df,val3.vec)

opt = which.min(val3.vec)
h_df[opt,] # h_X=0.21, h_Y=0.33

scatter3D(h_df$Var1,h_df$Var2,val3.vec, phi=0, theta=225, xlab="h_x", ylab="h_y", zlab="CV(h_x,h_y)",
          main="Cross-validation plot (optimal h_x=0.21, h_y=0.33)")


#------------------------------ problem b
# Ready
grid_b_X = b_X
grid_b_Y = b_Y

grid_b_X = seq(min(b_X),max(b_X),0.1)
grid_b_Y = seq(min(b_Y),max(b_Y),0.1)

h_x = 0.1
h_y = 0.1

b_F_hat = b_F_hat2 = matrix(0, length(grid_b_X),length(grid_b_Y) )

# Estimation
for (j in 1:length(grid_b_Y)){
  b_Y_idx = (b_Y <= grid_b_Y[j])+0
  for (i in 1:length(grid_b_X)){
    b_F_hat[i,j] = (mean(dnorm((b_X-grid_b_X[i])/h_x)*b_Y_idx)/mean(dnorm((b_X-grid_b_X[i])/h_x)))
  }
}

for (j in 1:length(grid_b_Y)){
  p = pnorm(grid_b_Y[j],b_Y,h_y)
  for (i in 1:length(grid_b_X)){
    b_F_hat2[i,j] = (mean(dnorm((b_X-grid_b_X[i])/h_x)*p)/mean(dnorm((b_X-grid_b_X[i])/h_x)))
  }
}

persp(grid_b_X,grid_b_Y,b_F_hat, theta = 0, phi = 20, expand = 0.5, col = "lightblue")
persp(grid_b_X,grid_b_Y,b_F_hat2, theta = 0, phi = 20, expand = 0.5, col = "lightblue")




