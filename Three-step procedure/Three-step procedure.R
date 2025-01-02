library(glmnet)
library(SLOPE)
library(ggplot2)
library(dplyr)
graphics.off()
set.seed(2)
save.fig=TRUE
n=100;
p=200;
beta.val <- c(20,10)

plot_fig <- function(beta.est, beta, name){
  
  
  df <- data.frame(
    index = seq_along(beta.est),
    betadebvec = beta.est,
    beta = beta
  )
  RMSE = sqrt(mean((beta-beta.est)^2))
  # Identify continuous segments of each unique beta value
  df_segments <- df %>%
    arrange(index) %>%
    group_by(beta) %>%
    mutate(segment = cumsum(c(1, diff(index) != 1))) %>%
    ungroup()
  
  fig <- ggplot(df, aes(x = index)) +
    # Plot points first
    geom_point(aes(y = betadebvec), color = "black", size = 2, alpha = 0.7) +
    # Plot lines on top of the points
    geom_line(
      data = df_segments,
      aes(y = beta, group = interaction(beta, segment), color = factor(beta)),
      size = 1.2
    ) +
    labs(
      title = paste0(name," (RMSE = ", round(RMSE, 3), ")"),
      x = "Index",
      y = "Estimated Value"
    ) +
    theme_minimal(base_size = 14) +
    # Remove the legend
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      legend.position = "none"
    )
  return(fig)
}


####How the data were generated

kor=0.8
Covar1<-matrix(kor,10,10)
diag(Covar1)=1
C1<-diag(20);
Covar<-kronecker(C1,Covar1);
diag(Covar)=1

Ch=t(chol(Covar))


X=matrix(rnorm(n*p),n,p) %*% t(Ch)



beta<-rep(0,p)
beta[1:65]<- beta.val[1]
beta[66:117]<- beta.val[2]
beta_unique <- unique(beta)
n.beta <- length(unique(beta))

Y<-X%*%beta+0.8*rnorm(n);
C<-t(X)%*%X/n;


###
# estimating through SLOPE
##
q=0.5;
seq1<-seq(1:p);
lambda<-qnorm(1-q*seq1/2/p);

####alpha was selected manually as to minimize MSE

B<-SLOPE(X,Y, lambda=lambda, alpha=0.07/sqrt(n),intercept=FALSE, scale='none', solver='admm')
betahat=coefficients(B);

fig <- plot_fig(betahat, beta, 'SLOPE')

if(save.fig){
  ggsave("slope_estimate.pdf",fig)
}else{
  print(fig)
}



####truncated SLOPE
library(grpSLOPE);




lambdatest=(42/sqrt(n))*lambda;

betacor<-grpSLOPE::prox_sorted_L1(betahat, lambda=lambdatest)


fig <- plot_fig(betacor, beta, 'Truncated SLOPE')
if(save.fig){
  ggsave("truncated_slope_estimate.pdf",fig)
}else{
  print(fig)
}



####debiased SLOPE
val_ <- unique(betacor)
val_ <- val_[val_>0]
n.clust <- length(val_)
tildeX <- matrix(0,nrow=n, ncol = n.clust)
for(i in 1:length(val_)){
  ind <- which(betacor==val_[i])
  tildeX[,i] =  rowSums(X[,ind])
}

betadeb<-solve(t(tildeX)%*%tildeX)%*%t(tildeX)%*%Y;
betadebvec <- rep(0,p)
for(i in 1:length(val_)){
  ind <- which(betacor==val_[i])
  betadebvec[ind] <- betadeb[i]
}

fig <- plot_fig(betadebvec, beta, 'Debiased truncated SLOPE')

if(save.fig){
  ggsave("debiased_truncated_slope_estimate.pdf",fig)
}else{
  print(fig)
}
