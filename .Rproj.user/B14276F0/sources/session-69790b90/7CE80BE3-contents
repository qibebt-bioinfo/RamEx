
#' Detect outliers in a dataset using various quality control methods
#'
#' This function applies multiple quality control methods to detect outliers in a dataset.
#' It uses Hotelling's T2, SNR, MCD, and distance methods to identify and mark outliers.
#'
#' @param matrix A matrix of spectral data where each row represents a single spectrum.
#' @param var_tol Tolerance for variation in the data. Defaults to 0.5.
#' @return A list containing the outliers and the run time for each method.
#' @export detect_all
#' @importFrom disprofas get_hotellings
#' @importFrom stats quantile
#' @importFrom stats mahalanobis
detect_all <- function(matrix, var_tol=0.5){
  data <- new('Ramanome', datasets = list(data=matrix), wavenumber=as.numeric(colnames(matrix)))
  data %<>% pre.spike() %>% pre.smooth(.,m = 0, p = 5, w = 7, delta.wav = 2) %>% pre.baseline %>% pre.normalize(.,'CH')
  matrix <- data@datasets$normalized.data
  outliers_all_ <- outliers_all <- data.frame(matrix('Bad',ncol=5, nrow=nrow(matrix)))
  colnames(outliers_all_) <- colnames(outliers_all) <- c('T2','SNR','mcd','dis','jump') #
  cal_time <- list()
  outlier_na <- which(is.na(rowSums(matrix)))
  if(length(outlier_na)!=0) {matrix=matrix[-outlier_na,];outliers_all=outliers_all_[-outlier_na,]}
  print('Start Hotelling T2!')
  start_time <- Sys.time()
  outliers_all$T2 <- detect_T2(matrix)$out %>% gsub('TRUE','Good',.) %>% gsub('FALSE','Bad',.)
  end_time <- Sys.time()
  cal_time$T2 <- as.numeric(difftime(end_time, start_time, units = "secs"))
  print('Start SNR!')
  start_time <- Sys.time()
  outliers_all$SNR <- spc_QC_SNR(matrix, level = 'medium') %>% gsub('TRUE','Good',.) %>% gsub('FALSE','Bad',.)
  end_time <- Sys.time()
  cal_time$SNR <- as.numeric(difftime(end_time, start_time, units = "secs"))
  print('Start PC-MCD!')
  start_time <- Sys.time()
  data.red <- prcomp_irlba(matrix, n = 20, center = TRUE, scale. = TRUE)$x[,1:5]
  outliers_all$mcd <- rep('Good', nrow(outliers_all))
  outliers_all$mcd[outliers_mcd_my(data.red,index_good=outliers_all$mcd=='Good', alpha = 0.001)$outliers_pos] <- 'Bad'
  end_time <- Sys.time()
  cal_time$mcd <- as.numeric(difftime(end_time, start_time, units = "secs"))
  print('Start Distance!')
  start_time <- Sys.time()
  outliers_all$dis <- detect_dis(matrix, min.dis=ncol(matrix)*0.01)$out %>% gsub('TRUE','Good',.) %>% gsub('FALSE','Bad',.)
  end_time <- Sys.time()
  cal_time$dis <- as.numeric(difftime(end_time, start_time, units = "secs"))
  # print('Start Valley-MCD!')
  # start_time <- Sys.time()
  # outliers_all$valley_mcd <- qc_jump_peaks_mcd(matrix)$index_good %>% gsub('TRUE','Good',.) %>% gsub('FALSE','Bad',.)
  # end_time <- Sys.time()
  # cal_time$valley_mcd <- as.numeric(difftime(end_time, start_time, units = "secs"))
  print('Start Jump!')
  start_time <- Sys.time()
  outliers_all$jump <- qc_jump(matrix, var_tol=var_tol, max_iterations=100)$index_good %>% gsub('TRUE','Good',.) %>% gsub('FALSE','Bad',.)
  end_time <- Sys.time()
  cal_time$jump <- as.numeric(difftime(end_time, start_time, units = "secs"))

  if(length(outlier_na)!=0) outliers_all_[-outlier_na,]=outliers_all else outliers_all_=outliers_all
  return(list(outliers=outliers_all_, run_time=cal_time))
}


#' Detect outliers in a dataset using Hotelling's T2 method
#'
#' This function identifies outliers in a dataset by applying Hotelling's T2 statistical method.
#' It uses parallel processing to speed up the computation.
#'
#' @param x A matrix of spectral data where each row represents a single spectrum.
#' @return A data frame containing two columns: 'out' indicating whether each spectrum is an outlier,
#' and 'dis' containing the corresponding Hotelling's T2 distance values.
#' @export detect_T2
#' @importFrom parallel makeCluster parApply
#' @importFrom disprofas get_hotellings
detect_T2 <- function(x){
  pred_outliers <- rep(TRUE, nrow(x))
  n_out <- 0
  temp_outliers <- pred_outliers
  cl <- makeCluster(getOption("cl.cores", detectCores() ))
  i <- 1
  while(TRUE){
    mean_spec <- as.matrix(parApply(cl, x[temp_outliers,],2, median))
    temp_outliers <- pred_outliers
    hotel_p <- parApply(cl,as.matrix(x),1,function(x,spec=mean_spec){
      return(disprofas::get_hotellings(as.matrix(x), spec, 0.05)$Parameters['p.F'])
    })
    temp_outliers[hotel_p < 0.05] <- FALSE
    if(n_out==length(temp_outliers[!temp_outliers]))
      break
    if(i > 30)
      break
    i <- i +1
    n_out <- length(pred_outliers[!temp_outliers])
    print(n_out)
  }
  return(data.frame(out=temp_outliers, dis=hotel_p))
}

#' Detect outliers in a dataset based on distance from the mean spectrum
#'
#' This function identifies outliers in a dataset by calculating the distance of each spectrum
#' from the mean spectrum of non-outlier spectra. The distance is calculated as the square
#' root of the sum of the squared differences between each spectrum and the mean spectrum.
#' Spectra with a distance greater than a specified threshold are marked as outliers.
#'
#' @param x A matrix of spectral data where each row represents a single spectrum.
#' @param min.dis The minimum distance threshold for marking a spectrum as an outlier.
#'                 Defaults to 1.
#' @return A data frame containing two columns: 'out' indicating whether each spectrum is an
#'         outlier, and 'dis' containing the corresponding distance values.
#' @export detect_dis
#' @importFrom stats quantile
detect_dis <- function(x, min.dis=1){
  pred_outliers <- rep(TRUE, nrow(x))
  n_out <- length(pred_outliers[!pred_outliers])
  mean_spec <- colMeans(x[pred_outliers,])
  dis_matrix <- apply(x, 1, function(row) sqrt(sum((row - mean_spec)^2)))
  temp.dis <- quantile(dis_matrix, probs = 0.9)
  while(TRUE){
    pred_outliers[dis_matrix>max(temp.dis, min.dis)] <- FALSE
    mean_spec <- colMeans(x[pred_outliers,])
    dis_matrix <- apply(x, 1, function(row) sqrt(sum((row - mean_spec)^2)))
    if(n_out==length(pred_outliers[!pred_outliers]))  {temp.dis <- temp.dis*0.95;pred_outliers[dis_matrix>max(temp.dis, min.dis)] <- FALSE}
    else{temp.dis <- quantile(dis_matrix, probs = 0.9)}
    print(temp.dis)
    n_out <- length(pred_outliers[!pred_outliers])
    if(all(dis_matrix[pred_outliers]<min.dis))
      break
  }

  return(data.frame(out=pred_outliers, dis=dis_matrix))
}
#' Detect outliers using the Minimum Covariance Determinant (MCD) method
#'
#' This function extends the outliers_mcdEst function to detect outliers in a dataset
#' using the Minimum Covariance Determinant method. It modifies the output to include
#' additional information such as the distance and center of the outliers.
#'
#' @param x The input dataset.
#' @param index_good A vector indicating the indices of good observations.
#' @param h The fraction of the data to use for computing the MCD.
#' @param alpha The significance level for the MCD.
#' @param na.rm A logical value indicating whether to remove NA values.
#' @return An object of class "outliers_mcd_my" containing the following components:
#'   - distance: The maximum distance threshold for outliers.
#'   - center: The center of the outliers.
#'   - outliers_pos: The positions of the outliers.
#'   - nb: A named number giving the total number of outliers detected.
#' @export outliers_mcd_my
#' @importFrom stats qchisq
outliers_mcd_my <- function(x,index_good,h = .5,alpha = .01, na.rm = TRUE){
  out <- outliers_mcdEst(x,index_good,h,alpha,na.rm)
  out$distance <- out$MaxDist
  out$center <- out$center
  out$call <- match.call()
  out$nb <- c(total = length(out$outliers_pos))

  class(out) <- "outliers_mcd_my"
  out
}


#' Estimate outliers using the Minimum Covariance Determinant (MCD) method
#'
#' This function implements the Minimum Covariance Determinant method to estimate
#' outliers in a dataset. It calculates the MCD estimate of scatter and location,
#' and identifies observations that are likely to be outliers based on their distance
#' from the centroid.
#'
#' @param x The input dataset.
#' @param index_good A vector indicating the indices of good observations to use for the estimation.
#' @param h The fraction of the data to use for computing the MCD (between 0 and 1). Defaults to 0.75.
#' @param alpha The significance level for the MCD. Defaults to 0.01.
#' @param na.rm A logical value indicating whether to remove NA values. Defaults to TRUE.
#' @return A list containing the following components:
#'   - MaxDist: The maximum distance threshold for outliers.
#'   - center: The estimated centroid of the data.
#'   - outliers_pos: The positions of the outliers.
#'   - outliers_val: The values of the outliers.
#'   - dist_from_center: The distances of each observation from the centroid.
#' @export outliers_mcdEst
#' @importFrom stats qchisq
#' @importFrom stats mahalanobis
outliers_mcdEst <- function(x,index_good,
                            h = .75, # fraction of data we wanna keep
                            # to compute the MCD (between 0 and 1)
                            alpha = .01,
                            na.rm = TRUE){

  if (na.rm == TRUE ) {
    data <- na.omit(x)
  } else {data <- x}

  for (i in seq_len(ncol(data))){
    if(inherits(data[,i],c("numeric","integer")) == FALSE)
      stop("Data are neither numeric nor integer")
  }

  #Creating covariance matrix for Minimum Covariance Determinant
  # by default, use the "best" method = exhaustive method
  print(nrow(data))
  # output <- fastMCD(data, h=0)
  output <- covFastMCD(data[index_good,], alpha = 0.6, m=10, l=1, delta = 0.05)
  print('FastMCD finished!')
  # output <- cov.mcd(data,cor = FALSE,quantile.used = nrow(data)*h)

  cutoff <- (qchisq(p = 1-alpha, df = ncol(data)))
  # cor = FALSE to avoid useless output(correlation matrix)

  #Distances from centroid for each matrix
  cat('Calculating the mahalanobis distance \n')
  dist <- mahalanobis(data,output$center,output$cov, tol=-1) # distance

  #Detecting outliers
  names_outliers <- which(dist > cutoff)
  coordinates <- data.frame(
    matrix(NA, nrow = length(names_outliers), ncol = ncol(data))
  )
  for (k in seq_len(ncol(data))){
    coordinates[,k] <- data[,k][dist > cutoff]
  }


  # print results
  meth <- "Minimum Covariance Determinant estimator"

  # Return results in list()
  invisible(
    list(MaxDist = cutoff,
         center = output$center,
         outliers_pos = names_outliers,
         outliers_val=coordinates,
         dist_from_center = dist)
  )
}


#' Fast Minimum Covariance Determinant Estimation
#'
#' This function implements the Fast Minimum Covariance Determinant (MCD) method to estimate
#' the covariance matrix and identify potential outliers in the data. The MCD method is
#' robust to outliers and can be used for multivariate data.
#'
#' @param x The input dataset, a matrix or data frame.
#' @param alpha The significance level for the MCD estimation.
#' @param m The number of initial subsets to consider.
#' @param l The number of final subsets to consider.
#' @param delta The quantile for the reweighting step.
#' @return A list containing the raw center, raw covariance matrix, best subset indices,
#' weights for reweighting, and the reweighted center and covariance matrix.
#' @export covFastMCD
#' @importFrom stats qchisq
#' @importFrom stats pgamma
#' @importFrom stats mahalanobis
covFastMCD = function(x, alpha, m, l, delta){
  x = data.frame(as.matrix(x,nrow = nrow(x),ncol = ncol(x)))
  p = ncol(x) #number of varaibles in data
  n = nrow(x) #number of observations in data
  h = h.alpha.n(alpha, n, p) #computing the subsetsize h using above p and n and given alpha
  determ = as.vector(rep(0,m)) #initialization of determinant vector
  q.alpha = qchisq(alpha,df = p)
  q.delta = qchisq((1-delta),df = p)
  c.a = alpha/pgamma(q.alpha/2, shape = p/2+1, scale = 1)
  c.delta = (1-delta)/pgamma(q.delta/2, shape = p/2+1, scale = 1)
  cnp = .MCDcnp2(p, n, alpha)
  #Before we start we might wanna check whther some conditions for inputes are satisifiede
  #condition of n&p
  if (n <= p){
    stop("Error: n <= p  the sample size is too small and MCD can not be performed!")

  }
  if(n == p + 1){
    stop("Error: n == p+1  is too small sample size for MCD")

  }
  if(n<2*p){
    warning("n<2*p this might be a too small sample size")
  }
  #condition of nh&alpha
  if(h > n ){
    stop("Error: h value is too large: h > n choose lower alpha ")
  }
  # if(h < (n + p + 1)/2 ){
  #   print("Error: h value is too small: h < (n+p+1)/2 choose higher alpha")
  #   break
  # }
  if(h == n){ #in this case the MCD location estimator T is average of the whole dataset and scatter MCD is the cov matrix
    T_mcd = apply(x,2,mean)
    S_mcd = (h-1)/h*var(x)
    mcdh = list(centerh = T_mcd,scatterh = S_mcd)
    return(mcdh)

  }
  # Determining top 10 initial subsets
  #initialization of all required lists in which we will store the generated samples,matrices and vectors
  T_0 = list()
  S_0 = list()
  H_0 = list()
  T_1 = list()
  S_1 = list()
  H_1 = list()
  T_2 = list()
  S_2 = list()
  H_2 = list()
  T_3 = list()
  S_3 = list()
  H_3 = list()
  H_final = list()
  list_sets = list()
  S_final = list()
  sigma = list()
  center = list()

  for(i in 1:m){

    H_0[[i]] = x[sample(1:n,p+1),] #randomly selecting p+1 observations from data
    T_0[[i]] = apply(H_0[[i]],2,mean)
    S_0[[i]] = p/(p+1)*var(H_0[[i]]) #correcting for 1/n-1 to have 1/n (multiplying (n-1)/n)

    #if the det(S_0) is 0 we should keep adding observations to random subset untill it becomes nonzero
    while (det(S_0[[i]]) == 0){
      for (k in 1:(n-p-1)){
        subsetnew = x[sample(n,(p+1)+k),]
      }
      H_0[[i]] = subsetnew
      T_0[[i]] = apply(subsetnew,2,mean)
      S_0[[i]] = (p+k)/(p+1+k)*var(subsetnew)
    }

    mah_dist0 = mahalanobis(x,T_0[[i]],S_0[[i]], tol=-1)
    d_sqr0 = sort(mah_dist0,decreasing = FALSE,index.return = TRUE)
    H1_index = d_sqr0$ix[1:h]
    H_1[[i]] = x[H1_index,] #the first subset with h components that will be used in C1-step
    #C1-step looking at H1
    T_1[[i]] = apply(H_1[[i]],2,mean) #center based on H1
    S_1[[i]] = (h-1)/h*var(H_1[[i]]) #scatter based on H1
    mah_dist1 = mahalanobis(x,T_1[[i]],S_1[[i]], tol=-1)
    d_sqr1 = sort(mah_dist1,decreasing = FALSE,index.return = TRUE)
    H2_index = d_sqr1$ix[1:h] #indeces of h small d^2's observations
    H_2[[i]] = x[H2_index,] #H2 subset
    T_2[[i]] = apply(H_2[[i]],2,mean) #center based on H2
    S_2[[i]] = (h-1)/h*var(H_2[[i]]) #scatter based on H2

    #C2-step looking at H2
    mah_dist2 = mahalanobis(x,T_2[[i]],S_2[[i]], tol=-1)
    d_sqr2= sort(mah_dist2,decreasing = FALSE,index.return = TRUE)
    H3_index = d_sqr2$ix[1:h] #indeces of h small d^2's observations
    H_3[[i]] = x[H3_index,] #H3 subset
    T_3[[i]] = apply(H_3[[i]],2,mean) #center based on H3
    S_3[[i]] = (h-1)/h*var(H_3[[i]]) #scatter based on H3
    determ[i] = det(S_3[[i]]) #we take the det of the best subset
  }


  D = sort(determ, decreasing = FALSE, index.return = TRUE) #sorting the deteminant of all m subsets  and picking smallest l determinant
  det10 = D$ix[1:l] #getting indexes of these 10 low det sets

  #getting the best busets and the corresponding covariances
  for(r in 1:l){
    list_sets[[r]] = H_3[[det10[r]]] #in this case H3 has been used
    S_final[[r]] = S_3[[det10[r]]]
  }

  objectmcd = list(call = match.call(),H0sets = H_0,H1sets = H_1,H2sets = H_2, H3sets = H_3,final10_set = list_sets, finalS = S_final, index_10 = det10)
  # while loop for C-steps of convergence for top 10 subsets from initial step
  for(k in 1:l){
    H_new = data.frame(as.matrix(objectmcd$final10_set[[k]],nrow = h,ncol = p)) #H1 for k =1 is the first subset of 10 finial subsets
    term1 = 1 #determinant of new one
    term2 = 2 #determinant of old one

    while(term1 != term2 && term1 != 0 ){ #continue c-steps until det(S_new) = det(S_old) or det(S_new) = 0

      T_old = apply(H_new,2,mean)
      S_old = (h-1)/h*var(H_new)
      term2 = det(S_old) #determinant of the initial subset
      mah_distnew = mahalanobis(x,T_old,S_old, tol=-1)
      d_ksqr = sort(mah_distnew,decreasing = FALSE,index.return = TRUE)
      index_new = d_ksqr$ix[1:h]
      #updating infromation based on new indices
      H_old = H_new  #new new sumbset becomes old subset to compare it weith new one
      H_new = x[index_new,] #updating H_new using the indices determined by H_old
      T_new = apply(H_new,2,mean)
      S_new = (h-1)/h*var(H_new)
      term1 = det(S_new) #updating new determinant

    }

    center[[k]] = T_new # storing the final center estimate in list of this final 10 centers
    sigma[[k]] = S_new  #store in the final 10 scatters
    H_final[[k]] = H_new  #store in the final 10 subsets

  }
  #deciding which one has the minimum determiniant of these final 10 sets
  d = rep(0,l)
  for(i in 1:l){
    d[i] = det(sigma[[i]])
  }
  det10= sort(d,decreasing = FALSE,index.return = TRUE) #det10 is sorted vector of top 10 final final determinant
  #index of the one with smallest determinant from all 10, where first one is the one with smallest det
  mindet = det10$ix[1] #the index of the best set  corresponding the smallest det
  # the final T and S estimators with index mindet
  T_mcd = center[[mindet]] #center estimater of the minimal det set
  S_mcd = c.a*cnp*sigma[[mindet]] #scatter estimator of the minimal det set
  H_mcd = H_final[[mindet]] #final subset with smallest det
  indeces_H_mcd = as.integer(rownames(H_mcd)) #indeces of the final subset
  result1 = list(raw.center  = T_mcd, raw.cov = S_mcd , best = indeces_H_mcd)
  # Reweighting step
  #(1-delta quantile of Chi-square dist. we will use 0.025 to expect that 97.5% of data will be included)
  w = rep(0,n) #initialization of weight vector
  d_sqr_r = mahalanobis(x,result1$raw.center,result1$raw.cov, tol=-1) #mah. distanaces for all observations
  cutoff = qchisq((1-delta), df = p) #quantile of chi-square distribution with o degrees of freedom
  for(i in 1:n){
    if(d_sqr_r[i]<=cutoff){
      w[i] = 1
    }
  }
  weights = t(w) #weights matrix
  # Reweighted MCD location estimator
  x = as.matrix(x,nrow=n,ncol = p)
  T_rwgt = 1/sum(weights)*(weights%*%x)
  #Reweighted MCD scatter estimator
  S = diag(0,p) #zero diagonal matrix
  sigmas = list() #list to stor each ith sigma matrix
  T = as.matrix(T_rwgt,nrow = nrow(T),ncol =ncol(T))
  for(i in 1:n){
    sigmas[[i]] = t(x[i,]-T)%*%(x[i,]-T)*weights[i]
    S = S + sigmas[[i]] #adding each time the new sigma
  }
  S_rwgt = S/sum(weights) #the sum of all covariance matrices devided by sum of weights
  Srwgt = c.delta*cnp*S_rwgt
  #the output adding the newresults to result1 with raw estimates
  list = list.append(result1, weights = weights,center = T_rwgt,cov = Srwgt )
  return(list)
}

#' Calculate the Signal-to-Noise Ratio (SNR) for a dataset
#'
#' This function calculates the SNR for a dataset by comparing the signal and noise levels.
#'
#' @param data A matrix of spectral data where each row represents a single spectrum.
#' @param level The level of strictness for the SNR calculation. Defaults to "medium".
#' @return A vector indicating whether each spectrum passes the SNR test.
#' @export spc_QC_SNR
#' @importFrom stats sd
spc_QC_SNR <- function(data, level = "medium")
{
  # library('baseline')
  # library('signal')
  data <- as.matrix(data)
  # data <- data[,apply(data,2,min)>15]
  # ba_li <- baseline.modpolyfit(t(data), degree = 3)$corrected
  # # smooth_data <- t(data)
  # for( i in 1:nrow(ba_li))
  # {
  #   smooth_data[i,] <- as.matrix(sgolayfilt(ba_li[i,],m=0,p=3,n=9))
  # }
  # nor_data <- smooth_data/apply(smooth_data, 1, max)
  #
  waves <- as.numeric(colnames(data))
  band_noise <- c(which((waves > 1800) & (waves < 1900)))
  band_signal <- c(which(waves > 600 & waves < 3050))
  if (level == 'hard')
  {
    QC1 <- 0.05
    QC2 <- 0.03
    QC3 <- 10
  }
  if (level == 'medium')
  {
    QC1 <- 0.03
    QC2 <- 0.03
    QC3 <- 5
  }
  if (level == 'easy')
  {
    QC1 <- 0.3
    QC2 <- 0.3
    QC3 <- 3
  }
  IBG <- apply(data[,band_noise],1,mean)
  IRAM <- apply(data[,band_signal],1,max)-IBG
  sdBG <- apply(data[,band_signal],1,sd)
  SNR <- IRAM/sdBG


  index_good <- apply(abs(data[,band_noise]),1,function(x) {return((sd(x)<QC2)&(mean(x)<QC1))}) &
    SNR > QC3 &     !duplicated(data)
  # data_good <- data[index_good,]
  index_good[is.na(index_good)] <- 'FALSE'
  return(index_good)
}
