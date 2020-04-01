
##### Determine the best E #####
TestE <- function(time_series,
                  d_name,
                  E_range,
                  nn = "e+1",
                  lib_type = "half",
                  lib_l = NULL,
                  pred_l = NULL,
                  visualize = TRUE,
                  out_pdf = TRUE,
                  E_only = FALSE){
  ts.l <- length(time_series)
  
  # Select library reconstruction strategies
  if(lib_type == "half"){
    # Half time series used to reconstruct library
    lib_l <- c(1, floor(ts.l/2))
    pred_l <- c(floor(ts.l/2) + 1, ts.l)
  }else if(lib_type == "full"){
    # Full time series used to reconstruct library
    lib_l <- c(1, ts.l)
    pred_l <- c(1, ts.l)
  }else if(lib_type == "manual"){
    lib_l <- lib_l
    pred_l <- pred_l
  }else{
    warning("Invalid lib_type option!")
  }
  
  # Do simplex projection to determine the best E
  simplex.out <- simplex(time_series, lib_l, pred_l, num_neighbors = nn, E = E_range, silent = T)
  
  if(visualize){
    # Visualize results
    op <- par(mfrow=c(1,3))
    plot(simplex.out$E, simplex.out$rho, type = "b")
    plot(simplex.out$E, simplex.out$mae, type = "b")
    plot(simplex.out$E, simplex.out$rmse, type = "b")
    par(op)
  }
  
  if(out_pdf){
  pdf(paste0("01_Embedding_Out/E_",d_name,".pdf"),width=1.6,height=1.5)
  par(mar = c(2.1, 2.6, 1.2, 0.5), mgp = c(2.5, 0.25, 0),ps=6)
  plot(simplex.out$E, simplex.out$rmse, type = "o", pch=16, cex=0.5, axes=F,ann=F,
       las=2,tcl=0.2)
  axis(side=1,tcl=-0.2,at=seq(4,24,4),labels=F)
  mtext(c("4","8","12","16","20","24"),at=seq(4,24,4),side=1,line=-0.15)
  axis(side=2,tcl=-0.2,las=2)
  mtext("Embedding dimension (E)",side=1,line=0.4)
  mtext("Forecast skill (RMSE)",side=2,line=1.2)
  mtext(d_name,side=3,line=0.1)
  box()
  dev.off()
  }
  
  if(E_only){
    return(simplex.out[which.min(simplex.out$rmse),"E"])
  }else{
    return(simplex.out)
  }
}



###### CCM with 100 seasonal surrogate function #####
CCMwSurrogate <- function(x, y, Ex, Ey, lag){
  x_xmap_y_mean <- rep(NA, length=length(lag))
  y_xmap_x_mean <- rep(NA, length=length(lag))
  
  for(i in 1:length(lag)){
    x_xmap_y <- ccm(d, E = Ex, lib_column = x, target_column = y, lib_sizes = 50, tp = lag[i], silent = T)
    y_xmap_x <- ccm(d, E = Ey, lib_column = y, target_column = x, lib_sizes = 50, tp = lag[i], silent = T)
    x_xmap_y_mean[i] <- ccm_means(x_xmap_y)$rho
    y_xmap_x_mean[i] <- ccm_means(y_xmap_x)$rho
  }
  
  
  # Seasonal surrogate 1 (x xmap y)
  num_surr <- 100
  surr_y <- make_surrogate_data(d[,y], method = "seasonal", T_period = 24, num_surr = num_surr)
  # Do CCM for all surrogate time series
  x_xmap_surr_y <- matrix(NA, nrow = length(lag), ncol = num_surr)
  for (i in 1:num_surr) {
    for(j in 1:length(lag)){
      x_xmap_surr_y[j,i] <- ccm_means(ccm(cbind(d[,x], surr_y[,i]), E = Ex, lib_sizes = 50, tp=lag[j], silent = T))$rho
    }
  }
  # Extract quantiles
  x_xmap_surr_y_qs <- t(apply(x_xmap_surr_y, 1, function(rhos) quantile(rhos, p = c(0.025, 0.5, 0.975))))
  # Combine resutls
  x_xmap_y_summary <- cbind(lag, x_xmap_y_mean, x_xmap_surr_y_qs)
  colnames(x_xmap_y_summary) <- c("lag_1", "rho_1", "lower2.5_1", "median_1", "upper97.5_1")
  
  # Seasonal surrogate 2 (y xmap x)
  surr_x <- make_surrogate_data(d[,x], method = "seasonal", T_period = 24, num_surr = num_surr)
  # Do CCM for all surrogate time series
  y_xmap_surr_x <- matrix(NA, nrow = length(lag), ncol = num_surr)
  for (i in 1:num_surr) {
    for(j in 1:length(lag)){
      y_xmap_surr_x[j,i] <- ccm_means(ccm(cbind(d[,y], surr_x[,i]), E = Ey, lib_sizes = 50, tp=lag[j], silent = T))$rho
    }
  }
  # Extract quantiles
  y_xmap_surr_x_qs <- t(apply(y_xmap_surr_x, 1, function(rhos) quantile(rhos, p = c(0.025, 0.5, 0.975))))
  # Combine resutls
  y_xmap_x_summary <- cbind(lag, y_xmap_x_mean, y_xmap_surr_x_qs)
  colnames(y_xmap_x_summary) <- c("lag_2", "rho_2", "lower2.5_2", "median_2", "upper97.5_2")
  
  # Return all result
  return(as.data.frame(cbind(x_xmap_y_summary, y_xmap_x_summary)))
}





###### CCM with 100 seasonal surrogate function for temperature #####
CCMwSurrogateTemp <- function(x, y, Ex, lag){
  x_xmap_y_mean <- rep(NA, length=length(lag))
  
  for(i in 1:length(lag)){
    x_xmap_y <- ccm(d, E = Ex, lib_column = x, target_column = y, lib_sizes = 50, tp = lag[i], silent = T)
    x_xmap_y_mean[i] <- ccm_means(x_xmap_y)$rho
  }
  
  
  # Seasonal surrogate 1 (x xmap y)
  num_surr <- 100
  surr_y <- make_surrogate_data(d[,y], method = "seasonal", T_period = 24, num_surr = num_surr)
  # Do CCM for all surrogate time series
  x_xmap_surr_y <- matrix(NA, nrow = length(lag), ncol = num_surr)
  for (i in 1:num_surr) {
    for(j in 1:length(lag)){
      x_xmap_surr_y[j,i] <- ccm_means(ccm(cbind(d[,x], surr_y[,i]), E = Ex, lib_sizes = 50, tp=lag[j], silent = T))$rho
    }
  }
  # Extract quantiles
  x_xmap_surr_y_qs <- t(apply(x_xmap_surr_y, 1, function(rhos) quantile(rhos, p = c(0.025, 0.5, 0.975))))
  # Combine resutls
  x_xmap_y_summary <- cbind(lag, x_xmap_y_mean, x_xmap_surr_y_qs)
  colnames(x_xmap_y_summary) <- c("lag_1", "rho_1", "lower2.5_1", "median_1", "upper97.5_1")
  
  # Return all result
  return(as.data.frame(x_xmap_y_summary))
}





##### CCM plot function #####
CCMplot <- function(ccm_res, x_name, y_name){
  
  g2<-ggplot(data = ccm_res, aes(x = lag_1)) + 
    geom_line(aes(y=rho_1, colour = "rho_1")) +
    geom_line(aes(y=rho_2, colour = "rho_2")) +
    # Add 95% confidence intervals
    geom_ribbon(ymin = ccm_res$lower2.5_1, ymax = ccm_res$upper97.5_1, alpha = 0.2, fill = rgb(1,0,1)) +
    geom_ribbon(ymin = ccm_res$lower2.5_2, ymax = ccm_res$upper97.5_2, alpha = 0.2, fill = rgb(0.55,0,0)) +
    geom_vline(xintercept = 0, linetype="dashed", color = "black", size=0.3) +
    scale_colour_manual("", breaks=c("rho_1","rho_2"), 
                        values=c(rgb(1,0,1),rgb(0.55,0,0)), 
                        labels=c(paste0(x_name," xmap ",y_name), 
                                 paste0(y_name," xmap ",x_name))) +
    scale_x_continuous(breaks=seq(ccm_res$lag_1[1],ccm_res$lag_1[nrow(ccm_res)],2)) +
    scale_y_continuous(breaks=seq(0,1,0.2)) +
    coord_cartesian(ylim=c(0,1)) +
    theme_bw() + 
    theme(legend.position=c(0.5,1.31),legend.key.height=unit(.3, "cm"),
          legend.key.width = unit(.4, "cm"),legend.direction="vertical",
          legend.background = element_rect(fill="transparent"),
          legend.spacing.x = unit(0.05, 'cm'),
          plot.margin=unit(c(1, 1, 1, 1),"lines"), 
          legend.text=element_text(size=6,hjust=0.5,margin=margin(b=2)), 
          axis.title=element_text(size=6), axis.text=element_text(size=6), 
          axis.title.x=element_text(margin = margin(t=-0.5)), 
          axis.title.y=element_text(margin=margin(r=-1)), 
          axis.text.x=element_text(margin=margin(t=0.5)), 
          axis.text.y=element_text(margin=margin(r=0.5))) + 
    labs(title="", x=expression(paste("Time to prediction (",italic(tp),")")), 
         y=expression(paste("Cross-map skill (", italic(rho), ")")))
  ggsave(paste0("02_CCM_Out/",x_name,"_",y_name,".pdf"), g2, 
         height = 1.7, width = 1.6)
  
}

CCMplotMagnif <- function(ccm_res, x_name, y_name){
  
  g2<-ggplot(data = ccm_res, aes(x = lag_1)) + 
    geom_line(aes(y=rho_1, colour = "rho_1")) +
    geom_line(aes(y=rho_2, colour = "rho_2")) +
    # Add 95% confidence intervals
    geom_ribbon(ymin = ccm_res$lower2.5_1, ymax = ccm_res$upper97.5_1, alpha = 0.2, fill = rgb(1,0,1)) +
    geom_ribbon(ymin = ccm_res$lower2.5_2, ymax = ccm_res$upper97.5_2, alpha = 0.2, fill = rgb(0.55,0,0)) +
    geom_vline(xintercept = 0, linetype="dashed", color = "black", size=0.3) +
    scale_colour_manual("", breaks=c("rho_1","rho_2"), 
                        values=c(rgb(1,0,1),rgb(0.55,0,0)), 
                        labels=c(paste0(x_name," xmap ",y_name), 
                                 paste0(y_name," xmap ",x_name))) +
    scale_x_continuous(breaks=seq(ccm_res$lag_1[1],ccm_res$lag_1[nrow(ccm_res)],2)) +
    scale_y_continuous(breaks=seq(0.8,1,0.05)) +
    coord_cartesian(ylim=c(0.8,1)) +
    theme_bw() + 
    theme(legend.position=c(0.5,1.31),legend.key.height=unit(.3, "cm"),
          legend.key.width = unit(.4, "cm"),legend.direction="vertical",
          legend.background = element_rect(fill="transparent"), 
          legend.spacing.x = unit(0.05, 'cm'),
          plot.margin=unit(c(1, 1, 1, 1),"lines"), 
          legend.text=element_text(size=6,hjust=0.5,margin=margin(b=2)), 
          axis.title=element_text(size=6), axis.text=element_text(size=6), 
          axis.title.x=element_text(margin = margin(t=-0.5)), 
          axis.title.y=element_text(margin=margin(r=-1)), 
          axis.text.x=element_text(margin=margin(t=0.5)), 
          axis.text.y=element_text(margin=margin(r=0.5))) + 
    labs(title="", x=expression(paste("Time to prediction (",italic(tp),")")),  
         y=expression(paste("Cross-map skill (", italic(rho), ")")))
  ggsave(paste0("02_CCM_Out/",x_name,"_",y_name,"magnif.pdf"), g2, 
         height = 1.7, width = 1.6)
  
}





##### CCM plot function for temperature #####
CCMplotTemp <- function(ccm_res, x_name, y_name){
  
  g2<-ggplot(data = ccm_res, aes(x = lag_1)) + 
    geom_line(aes(y=rho_1, colour = "rho_1")) +
    # Add 95% confidence intervals
    geom_ribbon(ymin = ccm_res$lower2.5_1, ymax = ccm_res$upper97.5_1, alpha = 0.2, fill = rgb(1,0,1)) +
    geom_vline(xintercept = 0, linetype="dashed", color = "black", size=0.3) +
    scale_colour_manual("", breaks=c("rho_1"), 
                        values=c(rgb(1,0,1)), 
                        labels=c(paste0(x_name," xmap ",y_name) 
                        )) +
    scale_x_continuous(breaks=seq(ccm_res$lag_1[1],ccm_res$lag_1[nrow(ccm_res)],2)) +
    scale_y_continuous(breaks=seq(0,1,0.2)) +
    coord_cartesian(ylim=c(0,1)) +
    theme_bw() + 
    theme(legend.position=c(0.5,1.22),legend.key.height=unit(.3, "cm"),
          legend.key.width = unit(.4, "cm"),legend.direction="vertical",
          legend.background = element_rect(fill="transparent"),
          legend.spacing.x = unit(0.05, 'cm'),
          plot.margin=unit(c(1, 1, 1, 1),"lines"), 
          legend.text=element_text(size=6,hjust=0.5,margin=margin(b=2)), 
          axis.title=element_text(size=6), axis.text=element_text(size=6), 
          axis.title.x=element_text(margin = margin(t=-0.5)), 
          axis.title.y=element_text(margin=margin(r=-1)), 
          axis.text.x=element_text(margin=margin(t=0.5)), 
          axis.text.y=element_text(margin=margin(r=0.5))) + 
    labs(title="", x=expression(paste("Time to prediction (",italic(tp),")")), 
         y=expression(paste("Cross-map skill (", italic(rho), ")")))
  ggsave(paste0("02_CCM_Out/",x_name,"_",y_name,".pdf"), g2, 
         height = 1.7, width = 1.6)
  
}

CCMplotMagnifTemp <- function(ccm_res, x_name, y_name){
  
  g2<-ggplot(data = ccm_res, aes(x = lag_1)) + 
    geom_line(aes(y=rho_1, colour = "rho_1")) +
    # Add 95% confidence intervals
    geom_ribbon(ymin = ccm_res$lower2.5_1, ymax = ccm_res$upper97.5_1, alpha = 0.2, fill = rgb(1,0,1)) +
    geom_vline(xintercept = 0, linetype="dashed", color = "black", size=0.3) +
    scale_colour_manual("", breaks=c("rho_1"), 
                        values=c(rgb(1,0,1)), 
                        labels=c(paste0(x_name," xmap ",y_name)
                        )) +
    scale_x_continuous(breaks=seq(ccm_res$lag_1[1],ccm_res$lag_1[nrow(ccm_res)],2)) +
    scale_y_continuous(breaks=seq(0.8,1,0.05)) +
    coord_cartesian(ylim=c(0.8,1)) +
    theme_bw() + 
    theme(legend.position=c(0.5,1.22),legend.key.height=unit(.3, "cm"),
          legend.key.width = unit(.4, "cm"),legend.direction="vertical",
          legend.background = element_rect(fill="transparent"), 
          legend.spacing.x = unit(0.05, 'cm'),
          plot.margin=unit(c(1, 1, 1, 1),"lines"), 
          legend.text=element_text(size=6,hjust=0.5,margin=margin(b=2)), 
          axis.title=element_text(size=6), axis.text=element_text(size=6), 
          axis.title.x=element_text(margin = margin(t=-0.5)), 
          axis.title.y=element_text(margin=margin(r=-1)), 
          axis.text.x=element_text(margin=margin(t=0.5)), 
          axis.text.y=element_text(margin=margin(r=0.5))) + 
    labs(title="", x=expression(paste("Time to prediction (",italic(tp),")")),  
         y=expression(paste("Cross-map skill (", italic(rho), ")")))
  ggsave(paste0("02_CCM_Out/",x_name,"_",y_name,"magnif.pdf"), g2, 
         height = 1.7, width = 1.6)
  
}





##### Check convergence #####
Convergence <- function(data, x, y, Ex, Ey, tpx, tpy){
  
  libx <- seq(Ex+1, 50, by = 1)
  liby <- seq(Ey+1, 50, by = 1)
  x_xmap_y <- ccm(data, E = Ex, lib_column = x, target_column = y, 
                  lib_sizes = libx, tp=tpx, silent=T)
  y_xmap_x <- ccm(data, E = Ey, lib_column = y, target_column = x, 
                  lib_sizes = liby, tp=tpy, silent=T)
  x_xmap_y_mean <- ccm_means(x_xmap_y)$rho
  y_xmap_x_mean <- ccm_means(y_xmap_x)$rho

  # Seasonal surrogate 1 (x xmap y)
  num_surr <- 100
  surr_y <- make_surrogate_data(data[,y], method = "seasonal", 
                                 T_period = 24, num_surr = num_surr)
  # Do CCM for all surrogate time series
  x_xmap_surr_y <- matrix(NA, nrow = length(libx), ncol = num_surr)
  for (i in 1:num_surr) {
    x_xmap_surr_y[,i] <- ccm_means(ccm(cbind(data[,x], surr_y[,i]), E=Ex,
                                      lib_column = 1, target_column = 2,
                                      lib_sizes = libx, tp=tpx, silent=T))$rho
  }

  # Extract quantiles
  x_xmap_surr_y_qs <- t(apply(x_xmap_surr_y, 1, function(rhos) quantile(rhos, p = c(0.025, 0.5, 0.975),na.rm=T)))
  # Combine resutls
  x_xmap_y_summary <- cbind(libx, x_xmap_y_mean, x_xmap_surr_y_qs)
  colnames(x_xmap_y_summary) <- c("lib_x", "rho_x", "lower2.5_x", "median_x", "upper97.5_x")
  
  
  # Seasonal surrogate 2 (y xmap x)
  num_surr <- 100
  surr_x <- make_surrogate_data(data[,x], method = "seasonal", 
                                T_period = 24, num_surr = num_surr)
  # Do CCM for all surrogate time series
  y_xmap_surr_x <- matrix(NA, nrow = length(liby), ncol = num_surr)
  for (i in 1:num_surr) {
    y_xmap_surr_x[,i] <- ccm_means(ccm(cbind(data[,y], surr_x[,i]), E=Ey,
                                       lib_column = 1, target_column = 2,
                                       lib_sizes = liby, tp=tpy, silent=T))$rho
  }
  
  # Extract quantiles
  y_xmap_surr_x_qs <- t(apply(y_xmap_surr_x, 1, function(rhos) quantile(rhos, p = c(0.025, 0.5, 0.975),na.rm=T)))
  # Combine resutls
  y_xmap_x_summary <- cbind(liby, y_xmap_x_mean, y_xmap_surr_x_qs)
  colnames(y_xmap_x_summary) <- c("lib_y", "rho_y", "lower2.5_y", "median_y", "upper97.5_y")

  # Return all result
  return(list(as.data.frame(x_xmap_y_summary), as.data.frame(y_xmap_x_summary)))
}





##### Check convergence for temperature #####
Convergence_Temp <- function(data, x, y, Ex, tpx){
  
  libx <- seq(Ex+1, 50, by = 1)
  x_xmap_y <- ccm(data, E = Ex, lib_column = x, target_column = y, 
                  lib_sizes = libx, tp=tpx, silent=T)
  x_xmap_y_mean <- ccm_means(x_xmap_y)$rho
  
  # Seasonal surrogate 1 (x xmap y)
  num_surr <- 100
  surr_y <- make_surrogate_data(data[,y], method = "seasonal", 
                                T_period = 24, num_surr = num_surr)
  # Do CCM for all surrogate time series
  x_xmap_surr_y <- matrix(NA, nrow = length(libx), ncol = num_surr)
  for (i in 1:num_surr) {
    x_xmap_surr_y[,i] <- ccm_means(ccm(cbind(data[,x], surr_y[,i]), E=Ex,
                                       lib_column = 1, target_column = 2,
                                       lib_sizes = libx, tp=tpx, silent=T))$rho
  }
  
  # Extract quantiles
  x_xmap_surr_y_qs <- t(apply(x_xmap_surr_y, 1, function(rhos) quantile(rhos, p = c(0.025, 0.5, 0.975),na.rm=T)))
  # Combine resutls
  x_xmap_y_summary <- cbind(libx, x_xmap_y_mean, x_xmap_surr_y_qs)
  colnames(x_xmap_y_summary) <- c("lib_x", "rho_x", "lower2.5_x", "median_x", "upper97.5_x")
  
  # Return all result
  return(list(as.data.frame(x_xmap_y_summary)))
}





##### Convergence plot function #####
Convergence_plot <- function(ccm_res, x, y, tpx, tpy){
  
  library(ggplot2)
  library(reshape2)
  g2 <- ggplot(data = ccm_res[[1]], aes(x = lib_x)) + 
    geom_ribbon(ymin = ccm_res[[1]]$lower2.5_x, ymax = ccm_res[[1]]$upper97.5_x, 
                alpha = 0.2, fill = rgb(1,0,1)) +
    geom_line(aes(y = rho_x, colour = "rho_x")) +
    scale_colour_manual("", breaks="rho_x", values=rgb(1,0,1)) +
    scale_x_continuous(breaks=seq(0,50,10)) +
    scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) +
    theme_bw() + 
    theme(legend.position="none", plot.margin=unit(c(1, 1, 1, 1),"lines"), 
         plot.title=element_text(size=6,hjust=0.5,margin=margin(b=2)), 
         axis.title=element_text(size=6), axis.text=element_text(size=6), 
         axis.title.x=element_text(margin = margin(t=-0.4)), 
         axis.title.y=element_text(margin=margin(r=-1)), 
         axis.text.x=element_text(margin=margin(t=0.5)), 
         axis.text.y=element_text(margin=margin(r=0.5))) + 
    labs(title=c(paste0(x," xmap ", y, "\n (tp = ",tpx,")")), x="Library size", 
         y=expression(paste("Cross-map skill (", italic(rho), ")")))
  ggsave(paste0("03_Convergence_Out/Convergence_",x,"xmap",y,"_tp=",tpx,".pdf"), g2, height = 1.6, width = 1.4)

  g2 <- ggplot(data = ccm_res[[2]], aes(x = lib_y)) + 
    geom_ribbon(ymin = ccm_res[[2]]$lower2.5_y, ymax = ccm_res[[2]]$upper97.5_y, 
                alpha = 0.2, fill = rgb(0.55,0,0)) +
    geom_line(aes(y = rho_y, colour = "rho_y")) +
    scale_colour_manual("", breaks="rho_y", values=rgb(0.55,0,0)) +
    scale_x_continuous(breaks=seq(0,50,10)) +
    scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) +
    theme_bw() + 
    theme(legend.position="none", plot.margin=unit(c(1, 1, 1, 1),"lines"), 
          plot.title=element_text(size=6,hjust=0.5,margin=margin(b=2)), 
          axis.title=element_text(size=6), axis.text=element_text(size=6), 
          axis.title.x=element_text(margin = margin(t=-0.4)), 
          axis.title.y=element_text(margin=margin(r=-1)), 
          axis.text.x=element_text(margin=margin(t=0.5)), 
          axis.text.y=element_text(margin=margin(r=0.5))) + 
    labs(title=c(paste0(y," xmap ", x, "\n (tp = ",tpy,")")), x="Library size", 
         y=expression(paste("Cross-map skill (", italic(rho), ")")))
  ggsave(paste0("03_Convergence_Out/Convergence_",y,"xmap",x,"_tp=",tpy,".pdf"), g2, height = 1.6, width = 1.4)
  
}





##### Convergence plot function for temperature #####
Convergence_plot_Temp <- function(ccm_res, x, y, tpx){
  
  library(ggplot2)
  library(reshape2)
  g2 <- ggplot(data = ccm_res[[1]], aes(x = lib_x)) + 
    geom_ribbon(ymin = ccm_res[[1]]$lower2.5_x, ymax = ccm_res[[1]]$upper97.5_x, 
                alpha = 0.2, fill = rgb(1,0,1)) +
    geom_line(aes(y = rho_x, colour = "rho_x")) +
    scale_colour_manual("", breaks="rho_x", values=rgb(1,0,1)) +
    scale_x_continuous(breaks=seq(0,50,10)) +
    scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) +
    theme_bw() + 
    theme(legend.position="none", plot.margin=unit(c(1, 1, 1, 1),"lines"), 
          plot.title=element_text(size=6,hjust=0.5,margin=margin(b=2)), 
          axis.title=element_text(size=6), axis.text=element_text(size=6), 
          axis.title.x=element_text(margin = margin(t=-0.4)), 
          axis.title.y=element_text(margin=margin(r=-1)), 
          axis.text.x=element_text(margin=margin(t=0.5)), 
          axis.text.y=element_text(margin=margin(r=0.5))) + 
    labs(title=c(paste0(x," xmap ", y, "\n (tp = ",tpx,")")), x="Library size", 
         y=expression(paste("Cross-map skill (", italic(rho), ")")))
  ggsave(paste0("03_Convergence_Out/Convergence_",x,"xmap",y,"_tp=",tpx,".pdf"), g2, height = 1.6, width = 1.4)
  
}


