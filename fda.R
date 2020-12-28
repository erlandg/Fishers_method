# 1-2: overlooked
# 3-(p-1): predictor values
# p: class
data = data.frame(read.delim('Data/t11-9.dat', header=FALSE, sep = "",
                             dec = ".", fill = TRUE, comment.char = "",
                             stringsAsFactors=FALSE))
# how many of the first variables to overlook
ovl = 2
p = ncol(data)-ovl-1
n = nrow(data)
# Three-class case
g = 3

# Not multinormal data (not neeeded for FDA)
library(MVN)
mvn(data.matrix(data[,(ovl+1):(p+ovl)]))
mu = colMeans(data.matrix(data[,(ovl+1):(p+ovl)]))

S.pooled = 0
# between-class sep.
B = 0
# Within-class sep.
W = 0
# iterate through classes
for (i in g) {
  temp.data = data[data[,(p+ovl+1)] == i,]
  temp.z = data.matrix(temp.data[,(ovl+1):(p+ovl)])
  S.pooled = S.pooled + ((nrow(temp.data)-1)*cov(temp.z)/
                           (nrow(data)-length(unique(data[,(p+ovl+1)]))))
  B = B + (colMeans(temp.z) - mu)%*%t(colMeans(temp.z) - mu)
  W = W + (nrow(temp.data)-1)*cov(temp.z)
}

eig = eigen(solve(W)%*%B)

# g = 3, p = 8
# s <= min(3-1, 8)
s = 2
# Discriminants
a = Re(eig$vectors[,1:s])

y = data.frame(cbind(data.matrix(data[,(ovl+1):(p+ovl)])%*%a, data[,(p+ovl+1)]))
colnames(y) = c('F1', 'F2', 'Group')

# Expectations E[Y|pi_i]. Equal for all Y
E1 = t(a)%*%colMeans(data.matrix(data[,(ovl+1):(p+ovl)])[data[,(p+ovl+1)] == 1,])
E2 = t(a)%*%colMeans(data.matrix(data[,(ovl+1):(p+ovl)])[data[,(p+ovl+1)] == 2,])
E3 = t(a)%*%colMeans(data.matrix(data[,(ovl+1):(p+ovl)])[data[,(p+ovl+1)] == 3,])
# may add arbitrary no. of expectations to extend to g > 3


# new sample vector
new.X = c(120,3,2,200,1,15,10,100)
Y = t(a)%*%new.X
print(Y)

Y_E1 = sqrt(sum((Y-E1)^2))
Y_E2 = sqrt(sum((Y-E2)^2))
Y_E3 = sqrt(sum((Y-E3)^2))
# Fisher's classification rule
classification.idx = print(which.min(c(Y_E1, Y_E2, Y_E3)))


# plot segregation in reduced dimensionality
res = 200
x.ax = c()
xex = c(min(c(y[,1], Y[1])), max(c(y[,1], Y[1])))
xex = c(xex[1]-.05*diff(xex), xex[2]+.05*diff(xex))
y.ax = c()
yex = c(min(c(y[,2], Y[2])), max(c(y[,2], Y[2])))
yex = c(yex[1]-.05*diff(yex), yex[2]+.05*diff(yex))
class.ax = c()
# Iterate through a 200 x 200 space to find segregation regions
for (x_ in seq(xex[1], xex[2], length.out=res)) {
  for (y_ in seq(yex[1], yex[2], length.out=res)) {
    x.ax = c(x.ax, x_)
    y.ax = c(y.ax, y_)
    Y.temp = c(x_, y_)
    Y_E1.temp = sqrt(sum((Y.temp-E1)^2))
    Y_E2.temp = sqrt(sum((Y.temp-E2)^2))
    Y_E3.temp = sqrt(sum((Y.temp-E3)^2))
    class.ax = c(class.ax, which.min(c(Y_E1.temp, Y_E2.temp, Y_E3.temp)))
  }
}
df = data.frame(x.ax, y.ax, class.ax)

# create figure with plotly
library(plotly)
fig = plot_ly(data=df, x=~x.ax, y=~y.ax, z=~class.ax,
              type='contour', contours=list(start=1,end=3,size=1)) %>%
  colorbar(title='Class') %>%
  add_trace(x=y[data[,(p+ovl+1)]==3,1], y=y[data[,(p+ovl+1)]==3,2], name='Class 3 (Q)',
            marker=list(color='green', size=10), type='scatter') %>%
  add_trace(x=y[data[,(p+ovl+1)]==2,1], y=y[data[,(p+ovl+1)]==2,2], name='Class 2 (K)',
            marker=list(color='blue', size=10), type='scatter') %>%
  add_trace(x=y[data[,(p+ovl+1)]==1,1], y=y[data[,(p+ovl+1)]==1,2], name='Class 1 (G)',
            marker=list(color='black', size=10), type='scatter')
fig


# Confusion matrix
conf = matrix(rep(0,9), nrow=3)
for (idx in seq(nrow(data))) {
  Y.temp = y[idx,1:2]
  Y_E1.temp = sqrt(sum((Y.temp-E1)^2))
  Y_E2.temp = sqrt(sum((Y.temp-E2)^2))
  Y_E3.temp = sqrt(sum((Y.temp-E3)^2))
  cl.temp = which.min(c(Y_E1.temp, Y_E2.temp, Y_E3.temp))
  conf[as.integer(data[idx,(p+ovl+1)]),as.integer(cl.temp)] = 1 +
    conf[as.integer(data[idx,(p+ovl+1)]),as.integer(cl.temp)]
}
# Confusion matrix
print(conf)
# Error rate
print(1 - sum(diag(conf))/sum(conf))
