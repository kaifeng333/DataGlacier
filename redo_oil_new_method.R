##test git


# library(extraDistr)    ##library to generate Dirichlet random variables
library(MCMCpack)
library(Rcpp)
sourceCpp("oil_decision_VPIVOI.cpp")

##define a function to compute VOI
VOI_compute = function(f_good,pmf_prob,p_over_c,n_wells){
  x_wells_list = seq(0,n_wells)
  p_l_bin_mat = matrix(NaN,nrow = length(f_good),ncol = length(x_wells_list))
  E_p_over_c_vec = rep(NaN,length(x_wells_list))
  p_l_bin_check = rep(NaN,length(x_wells_list))
  for (j in 1:length(x_wells_list)){
    p_l_bin_mat[,j] = dbinom(x_wells_list[j],n_wells,f_good)
    p_l_bin_check[j] = sum(p_l_bin_mat[,j]*pmf_prob)
    E_p_over_c_vec[j] = pmax(sum(p_l_bin_mat[,j]*pmf_prob/
                                   p_l_bin_check[j]*p_over_c),0)
  }
  VOI = sum(E_p_over_c_vec*p_l_bin_check)-pmax(sum(pmf_prob*p_over_c),0)
  return(VOI)
}



###############
##test codes


##compute based on utility values

delta_f = 1e-1
##key variable: f,frequency of good wells in the play
f_all = seq(0,1,delta_f)

f_good = f_all[1:(length(f_all)-1)] + diff(f_all)/2
length(f_good)

##variable for sensitivity analysis: r/c, return on investment ratio
##how much return per cost or investment
# r_over_c = unique(c(seq(1,2,0.05),seq(2,10,0.1),seq(10,20,0.5),seq(10,100,10)))
# r_over_c = c(1.05,1.1,2,5)
r_over_c = 1.25
length(r_over_c)

#######################################################################
##computation
##compute profit if drill a well
##negative choose not develop, positive choose develop
##profit over cost, non-dimensional
##cost is positive value; profit has signs
p_over_c = matrix(NaN,nrow = length(f_good),ncol = length(r_over_c))

for (j in 1:length(r_over_c)){
  p_over_c[,j] = -1 + f_good*r_over_c[j]
}

p_over_c = data.frame(p_over_c)
colnames(p_over_c) = r_over_c

summary(p_over_c)
length(p_over_c)

##visualize when r/c = 5 (utility values)
# plot(f_good,p_over_c[,which(r_over_c==5)],type = 'l',lwd = 3,ylim = c(-1,4))



##############################################
##for all outcomes, Monte Carlo
tol_VOI = 1e-10
num_R = 1e+6

set.seed(1234)


####
##define a function to generate random probabilities
##num_R is the number of realizations of MC
##f_good is the frequency of good wells: # of outcomes for each alternative
prob_mat_fxn = function(num_R,f_good){
  prob_mat = matrix(NaN,nrow = num_R,ncol = length(f_good))

  # ##method 1
  # alpha_vec = rep(1,length(f_good))
  # for (j in 1:length(f_good)){
  #   prob_mat[,j] = rgamma(num_R,alpha_vec[j],1)
  #   prob_mat[,j] = rgamma(num_R,alpha_vec[j],1)
  # }
  # prob_mat = prob_mat/rowSums(prob_mat)

  ##method 2 (not work for lots of realizations, i.e., 10 million)
  prob_mat = rdirichlet(num_R,rep(1,length(f_good)))

  # prob_mat = data.frame(prob_mat)
  return(prob_mat)
}

t_begin = Sys.time()
prob_mat_A1 = prob_mat_fxn(num_R,f_good)
Sys.time()-t_begin

summary(prob_mat_A1)

# prob_mat_A1[1:10,] = diag(x=1,10)
###############################################
##VPI analysis
##VPI for all r_over_c values
acc_times_A1 = 1e+4

num_dis_list = seq(11,101,10)
# num_dis_list = 21

VPI_realization_A1_mat = matrix(NaN,nrow = num_R,ncol = length(r_over_c))
##prob for num of realizations of
p_VPI_vec_A1_mat = matrix(NaN,nrow = num_R,ncol = length(r_over_c))

p_VPI_check_A1 = rep(0,length(r_over_c))

##output matrix
mean_VPI_A1_vec = rep(NaN,length(r_over_c))
E_joint_A1_mat = matrix(NaN,nrow = length(f_good),ncol = length(r_over_c))

# pro_VPI = txtProgressBar(1,length(r_over_c),style = 3)
t_begin = Sys.time()
for (j in 1:length(r_over_c)){
  # setTxtProgressBar(pro_VPI,j)
  ##each r_over_c, use different random probabilities
  prob_mat_A1 = prob_mat_fxn(num_R,f_good)

  VPI_realization_A1_mat[,j] = rowSums(t(t(prob_mat_A1)*pmax(0,p_over_c[,j]))) -
    pmax(rowSums(t(t(prob_mat_A1)*p_over_c[,j])),0)

  min_VPI_A1 = floor(min(VPI_realization_A1_mat[,j])*acc_times_A1)/acc_times_A1
  max_VPI_A1 = ceiling(max(VPI_realization_A1_mat[,j])*acc_times_A1)/acc_times_A1

  if (abs(max_VPI_A1 - min_VPI_A1)<0.05){
    p_VPI_vec_A1 = rep(1/num_R,num_R)
    p_VPI_check_A1[j] = sum(p_VPI_vec_A1)
  }
  else{
    k = length(num_dis_list)
    ##p_VPI_vec is the corrected f(p11,p12)
    while(abs(p_VPI_check_A1[j] - 1) > 1e-6){
      num_dis = num_dis_list[k]
      VPI_pos_bin_A1 = seq(min_VPI_A1,max_VPI_A1,length.out = num_dis)
      p_mass_pos_VPI_A1 = p_mass_pos_fxn(VPI_pos_bin_A1, min(VPI_realization_A1_mat[,j]))

      p_VPI_vec_A1_mat[,j] = p_potential(VPI_realization_A1_mat[,j],VPI_pos_bin_A1,p_mass_pos_VPI_A1,num_dis)

      p_VPI_check_A1[j] = sum(p_VPI_vec_A1_mat[,j])
      k = k-1
    }
  }
  ##outputs
  mean_VPI_A1_vec[j] = sum(VPI_realization_A1_mat[,j]*p_VPI_vec_A1_mat[,j])
  E_joint_A1_mat[,j] = colSums(prob_mat_A1*p_VPI_vec_A1_mat[,j])
}
Sys.time()-t_begin

summary(p_VPI_check_A1)
summary(VPI_realization_A1_mat)

VPI_pos_bin_A1
p_mass_pos_VPI_A1

# ##check the expected joint (make sure all 1's)
# colSums(E_joint_A1_mat)
# dim(E_joint_A1_mat)
# 
# mean_VPI_A1_vec
# 
# summary(E_joint_A1_mat)
# 
# 
# 
# colSums(p_VPI_vec_A1_mat)
# summary(VPI_realization_A1_mat)

##################################################################
##VPI analysis
##compare with C++ codes
num_dis_list = seq(11,101,10)


##compute VOI first
##likelihood function
n_wells = 5
# n_wells = 10
x_wells_list = seq(0,n_wells)
p_l_bin_mat = matrix(NaN,nrow = length(f_good),ncol = length(x_wells_list))
p_l_bin_mat_check = matrix(NaN,nrow = num_R,ncol = n_wells + 1)
for (j in 1:length(x_wells_list)){
  p_l_bin_mat[,j] = dbinom(x_wells_list[j],n_wells,f_good)
}

##posterior and VOI
prob_mat_A1_post = matrix(NaN,nrow = num_R,ncol = length(f_good))
E_p_over_c_post = matrix(NaN,nrow = num_R,ncol = n_wells+1)
for (j in 1:length(x_wells_list)){
  p_l_bin_mat_check[,j] = colSums(p_l_bin_mat[,j]*t(prob_mat_A1))
  prob_mat_A1_post = t(p_l_bin_mat[,j]*t(prob_mat_A1))/p_l_bin_mat_check[,j]
  E_p_over_c_post[,j] = colSums(t(prob_mat_A1_post)*p_over_c[,1])
}
##check
summary(rowSums(p_l_bin_mat_check))

VOI_realization = rowSums(pmax(E_p_over_c_post,0)*p_l_bin_mat_check)-
  pmax(rowSums(t(t(prob_mat_A1)*p_over_c[,1])),0)
VOI_realization[abs(VOI_realization)<1e-10] = 0

summary(VPI_realization_A1_mat)
summary(VOI_realization)



#############################################
##adjust prob
VPI_realization = VPI_realization_A1_mat[,1]
VPI_pos_bin = VPI_pos_bin_A1
p_mass_pos_VPI = p_mass_pos_VPI_A1



num_dis_VOI_list = c(3,6,seq(11,101,10))
VOI_dis_mat = matrix(NaN,nrow = max(num_dis_VOI_list),ncol = length(p_mass_pos_VPI))
res_num_dis_VPI = rep(NaN,length(p_mass_pos_VPI))

num_min_VOI_not0 = rep(NaN,length(p_mass_pos_VPI))
p_check = rep(0,length(p_mass_pos_VPI))
VOI_ocheck = rep(NaN,length(p_mass_pos_VPI))
VOI_num_check = rep(NaN,length(p_mass_pos_VPI))
p_VOI_vec = rep(NaN,length(VPI_realization))
pro_voi = txtProgressBar(2,length(VPI_pos_bin),style = 3)
t_begin = Sys.time()
##first assign prob for VPI=VOIO=0 (if there is VPI=0)
if (min(VPI_realization)==0){
  p_VOI_vec[which(VPI_realization==0)] = 0.5/sum(VPI_realization==0)
}
for (i in 2:length(VPI_pos_bin)){
  setTxtProgressBar(pro_voi,i)
  
  if (i==2 & (min(VPI_realization)>0)){
    index = which(VPI_realization>=VPI_pos_bin[i-1] & VPI_realization<=VPI_pos_bin[i])
  }
  else{
    index = which(VPI_realization>VPI_pos_bin[i-1] & VPI_realization<=VPI_pos_bin[i])
  }
  weight_i = p_mass_pos_VPI[i-1]
  temp_i = VOI_realization[index]
  
  ##output the matrix of each VPI_pos_bin
  
  ##check decisions of each possible X
  E_p_over_c_post_temp_i = E_p_over_c_post[index,]
  decision_mat = pmax(E_p_over_c_post_temp_i,0)
  decision_mat = E_p_over_c_post_temp_i - decision_mat
  decision_mat[decision_mat<0] = 1
  decision_mat = decision_mat + 1
  
  if (length(index)>1){
    ##note this fxn asks index at least 2
    index_same = union(which(rowSums(decision_mat)==6),which(rowSums(decision_mat)==12))
    
    ##temporary vector to store assigned probabilities
    p_VOI_vec_temp = rep(NaN,length(index))
    
    ##if there are VOI=0 or min(VOI) that decision not change 
    ##no matter what k1 measurements
    
    if ((0 < length(index_same)) & (length(index_same) < length(index)) & (length(table(temp_i))>=5)){
      ##prob for same indices
      p_VOI_vec_temp = rep(NaN,length(index))
      p_VOI_vec_temp[index_same] = 0.5/length(index_same)*weight_i
      ##prob for different indices (VOI is not minimum and meaningful)
      temp_i_remain = temp_i[-index_same]
      min_VOI_remain = floor(min(temp_i_remain)*10000)/10000
      max_VOI_remain = ceiling(max(temp_i_remain)*10000)/10000
      c(min_VOI_remain,max_VOI_remain)
      k = length(num_dis_VOI_list)
      if ((max_VOI_remain-min_VOI_remain)/min_VOI_remain>5e-2){
        while (p_check[i-1]!=1){
          num_dis_VOI = num_dis_VOI_list[k]
          k = k-1
          VOI_pos_bin = seq(min_VOI_remain,max_VOI_remain,length.out = num_dis_VOI)
          p_mass_pos_VOI = p_mass_pos_fxn_VOI(VOI_pos_bin)
          
          p_VOI_vec_temp[-index_same] = p_potential_VOI(temp_i_remain,VOI_pos_bin,p_mass_pos_VOI,num_dis_VOI)*weight_i
          
          p_VOI_vec[index] = p_VOI_vec_temp
          p_check[i-1] = round(sum(p_VOI_vec[index]/weight_i),5)
          
          if (k<0){
            break
          }
        }
      }
      else{
        p_VOI_vec_temp[-index_same] = 0.5/length(temp_i_remain)*weight_i
        
        p_VOI_vec[index] = p_VOI_vec_temp
        p_check[i-1] = round(sum(p_VOI_vec[index]/weight_i),5)
      }
      
    }
    else if (length(index_same)==length(index) || (length(table(temp_i))<5)){
      p_VOI_vec[index] = 1.0/length(temp_i)*weight_i
      p_check[i-1] = round(sum(p_VOI_vec[index]/weight_i),5)
      num_dis_VOI = 0
      VOI_pos_bin = rep(NaN,5)
    }
    ##last scenario: valuable info, decisions different
    else{
      min_VOI = floor(min(temp_i)*10000)/10000
      max_VOI = ceiling(max(temp_i)*10000)/10000
      c(min_VOI,max_VOI)
      k = length(num_dis_VOI_list)
      while (p_check[i-1]!=1){
        num_dis_VOI = num_dis_VOI_list[k]
        k = k-1
        VOI_pos_bin = seq(min_VOI,max_VOI,length.out = num_dis_VOI)
        p_mass_pos_VOI = p_mass_pos_fxn(VOI_pos_bin,min(temp_i))
        
        p_VOI_vec[index] = p_potential(temp_i,VOI_pos_bin,p_mass_pos_VOI,num_dis_VOI)*weight_i
        
        p_check[i-1] = round(sum(p_VOI_vec[index]/weight_i),5)
        
        if (k<0){
          break
        }
      }
    }
  }
  else{
    p_VOI_vec[index] = 1.0/length(temp_i)*weight_i
    p_check[i-1] = round(sum(p_VOI_vec[index]/weight_i),5)
    num_dis_VOI = 0
    VOI_pos_bin = rep(NaN,5)
  }
  if (min(temp_i)>0){
    num_min_VOI_not0[i-1] = sum(abs(temp_i-min(temp_i))<1e-10)
    # print(table(decision_k1))
    # print(table(decision_k2))
  }
  res_num_dis_VPI[i-1] = num_dis_VOI
  VOI_ocheck[i-1] = length(which(temp_i==0))
  VOI_num_check[i-1] = length(temp_i)
  VOI_dis_mat[1:length(VOI_pos_bin),i-1] = VOI_pos_bin
}
t_end = Sys.time()
t_end-t_begin
sum(p_VOI_vec)

num_min_VOI_not0
p_check
VOI_ocheck
VOI_dis_mat


#####################################
##output
mean_VPI = sum(VPI_realization*p_VPI_vec_A1_mat[,1])
mean_VPI
mean_VOI = sum(VOI_realization*p_VOI_vec)
mean_VOI


##PMF
PMF_VOI = colSums(prob_mat_A1*p_VOI_vec)
PMF_VOI
sum(PMF_VOI)

##plug PMF_VOI into VOI fxn
VOI_compute(f_good,PMF_VOI,p_over_c[,1],n_wells)




##Prior PMF plots based on VPI analysis adn VOI analysis
barplot(rbind(E_joint_A1_mat[,1],PMF_VOI),
        names.arg = f_good,beside = TRUE,ylim = c(0,1),col = c('gray40','blue'),
        ylab = 'PMF',cex.lab = 1.2,cex.axis = 1.2,cex.names = 1.2)
title(xlab = 'f: frequency of good wells in the play',cex.lab = 1.2)
text(20,0.9,paste0('r/c = ',r_over_c),cex = 1.2)
legend('topleft',c('PMF based on VPI Analysis','PMF based on VOI Analysis'),
       col = c('gray40','blue'),cex = 1.2,
       bty = 'n',pch = c(15,15))

# ##pdf
# plot(f_good,E_joint_A1_mat[,1]/diff(f_good)[1],type = 'l',lwd = 3,xlab = '',
#      ylab = 'PDF',cex.axis = 1.2,cex.lab = 1.2,ylim = c(0,3))
# title(xlab = 'f: frequency of good wells in the play',cex.lab = 1.2)
# lines(f_good,PMF_VOI/diff(f_good)[1],lwd = 3,lty = 2,col = 'blue')
# text(0.05,2.8,paste0('r/c = ',r_over_c),cex = 1.2)
# legend('topright',c('PMF based on VPI Analysis','PMF based on VOI Analysis'),
#        col = c('gray40','blue'),cex = 1.2,
#        bty = 'n',pch = c(15,15))


#################################################
##plot PDF of VPI vs VOI
#################################################
VPI_pos_bin_A1
p_mass_pos_VPI_A1

##number of zeros and minimum non-zeros of each VPI bin
num_min_VOI_not0
p_check
VOI_ocheck
VOI_dis_mat



VPI_pos_bin_mid = VPI_pos_bin_A1[1:(length(VPI_pos_bin_A1)-1)]+
  diff(VPI_pos_bin_A1)/2
VPI_pos_bin_mid   ## this is the VPI-axis for Matlab 3-D plots
delta_VPI = diff(VPI_pos_bin_A1)[1]
delta_VPI

delta_VOI = rep(1,length(VPI_pos_bin_mid))

##1st: if we only output non-zero VOI realizations of each VPI bin
summary(p_VOI_vec)
sum(p_VOI_vec)

p_VOI_vec_pdf = rep(NaN,length(p_VOI_vec))

##check
p_VOI_pdf_check = rep(NaN,length(VPI_pos_bin_mid))
p_VOI_pdf_check_non_zero = rep(NaN,length(VPI_pos_bin_mid))

p_VOI_pmf_check_zero = rep(NaN,length(VPI_pos_bin_mid))

##matlab input
plt_matlab_VOI = rep(NaN,3)
plt_matlab_VOI_zero = rep(NaN,3)

pro_VOI_pdf = txtProgressBar(1,length(VPI_pos_bin_mid),style = 3)

t_begin = Sys.time()
for (i in 1:length(VPI_pos_bin_mid)){
  setTxtProgressBar(pro_VOI_pdf,i)
  
  if (i==1){
    index = which(VPI_realization>=VPI_pos_bin[i] & VPI_realization<=VPI_pos_bin[i+1])
  }
  else{
    index = which(VPI_realization>VPI_pos_bin[i] & VPI_realization<=VPI_pos_bin[i+1])
  }
  
  VOI_temp_i = VOI_realization[index]
  p_VOI_vec_temp_i = p_VOI_vec[index]
  p_VOI_vec_pdf_temp_i = p_VOI_vec_pdf[index]
  ##if VPI bin all zero-VOIs
  
  
  ##if the VPI bin not all zero-VOIs
  if (max(VOI_temp_i)>0){
    ##delta_VOI
    if (!is.na(VOI_dis_mat[1,i])){
      delta_VOI[i] = diff(VOI_dis_mat[,i])[1]    ##use VOI-bins while VOI-analysis
      VOI_pos_bin_temp_i = na.omit(VOI_dis_mat[,i])
      VOI_pos_bin_temp_i_mid = VOI_pos_bin_temp_i[1:(length(VOI_pos_bin_temp_i)-1)]+
        diff(VOI_pos_bin_temp_i)/2
      
      mat_temp = cbind(rep(VPI_pos_bin_mid[i],length(VOI_pos_bin_temp_i_mid)),
                       VOI_pos_bin_temp_i_mid,
                       rep(1,length(VOI_pos_bin_temp_i_mid)))
    }
    else{
      VOI_pos_bin_temp_i_mid = mean(VOI_temp_i)
      
      mat_temp = cbind(VPI_pos_bin_mid[i],
                       VOI_pos_bin_temp_i_mid,
                       1)
    }
    
    
    index_zero_i = which(VOI_temp_i==0)
    if (length(index_zero_i)>0){
      # p_VOI_vec_pdf_temp_i[-index_zero_i] = p_VOI_vec_temp_i[-index_zero_i]/delta_VOI[i]/delta_VPI
      p_VOI_pmf_check_zero[i] = sum(p_VOI_vec_temp_i[index_zero_i])
      p_VOI_pdf_check_non_zero[i] = sum(p_VOI_vec_temp_i[-index_zero_i])
    }
    else{
      # p_VOI_vec_pdf_temp_i = p_VOI_vec_temp_i/delta_VOI[i]/delta_VPI
      p_VOI_pdf_check_non_zero[i] = sum(p_VOI_vec_temp_i)
    }
    mat_temp[,3] = rep(p_VOI_pdf_check_non_zero[i]/length(VOI_pos_bin_temp_i_mid),
                       dim(mat_temp)[1])
    plt_matlab_VOI = rbind(plt_matlab_VOI,mat_temp)
  }
  # p_VOI_vec_pdf[index] = p_VOI_vec_pdf_temp_i
  
  p_VOI_pdf_check[i] = sum(p_VOI_vec_temp_i)
}
t_end = Sys.time()
t_end - t_begin

p_VOI_pdf_check
p_VOI_pdf_check_non_zero
p_VOI_pmf_check_zero

# summary(p_VOI_vec_pdf)
sum(VOI_ocheck)

plt_matlab_VOI = na.omit(plt_matlab_VOI)
summary(plt_matlab_VOI)
dim(plt_matlab_VOI)

##for zero VOis
plt_matlab_VOI_zero = cbind(VPI_pos_bin_mid,rep(0,length(VPI_pos_bin_mid)),
                            p_VOI_pmf_check_zero)

##output data for Matlab
# write.csv(plt_matlab_VOI,file = 'plt_matlab_VOI.csv')

# 
# ##prepare for 3D surafae plots
# ##change pdf to pmf
# plt_matlab_VOI_pmf = plt_matlab_VOI
# pmf_check = rep(NaN,length(VPI_pos_bin_mid))
# for (i in 1:length(VPI_pos_bin_mid)){
#   index_i = which(plt_matlab_VOI[,1]==VPI_pos_bin_mid[i])
#   plt_matlab_VOI_pmf[index_i,3] = plt_matlab_VOI[index_i,3]*delta_VPI*delta_VOI[i]
#   # pmf_check[i] = sum(plt_matlab_VOI_pmf[index_i,3])
# }
# 
# # pmf_check
# summary(plt_matlab_VOI_pmf)
# dim(plt_matlab_VOI_pmf)
# 
# summary(VOI_realization)

plt_matlab_VOI_pmf = plt_matlab_VOI
######
##VOI_pos_bin_plt only used for 3D surface plots in Matlab
num_dis_plt = 101

VPI_pos_bin_plt_mid = unique(plt_matlab_VOI_pmf[,1])
VPI_pos_bin_plt_mid


VOI_pos_bin_plt = seq(0,max(VOI_realization),length.out = num_dis_plt)
##make sure last element is max
VOI_pos_bin_plt[length(VOI_pos_bin_plt)] = max(VOI_realization)

VOI_pos_bin_plt_mid = VOI_pos_bin_plt[1:(length(VOI_pos_bin_plt)-1)]+diff(VOI_pos_bin_plt)/2
VOI_pos_bin_plt_mid



plt_matlab_VOI_surf = matrix(0,nrow = length(VPI_pos_bin_plt_mid),ncol = length(VOI_pos_bin_plt_mid))

pro_plt = txtProgressBar(1,dim(plt_matlab_VOI_surf)[1],style = 3)
for (i in 1:dim(plt_matlab_VOI_surf)[1]){
  setTxtProgressBar(pro_plt,i)
  # VOI_plt_temp = plt_matlab_VOI_pmf[which(plt_matlab_VOI_pmf[,1]==VPI_pos_bin_plt_mid[i]),]
  VOI_plt_temp = rep(NaN,3)
  VOI_plt_temp = na.omit(rbind(VOI_plt_temp,
                               plt_matlab_VOI_pmf[which(plt_matlab_VOI_pmf[,1]==VPI_pos_bin_plt_mid[i]),]))
  for (j in 1:dim(VOI_plt_temp)[1]){
    index = floor((VOI_plt_temp[j,2]-1e-16)/diff(VOI_pos_bin_plt)[1])+1
    plt_matlab_VOI_surf[i,index] = plt_matlab_VOI_surf[i,index]+VOI_plt_temp[j,3]
  }
}
summary(plt_matlab_VOI_surf) 
  
##pdf
plt_matlab_VOI_surf = plt_matlab_VOI_surf/diff(VPI_pos_bin_plt_mid)[1]/diff(VOI_pos_bin_plt_mid)[1]
max(plt_matlab_VOI_surf) 

##check dimensions
c(length(VPI_pos_bin_plt_mid),length(VOI_pos_bin_plt_mid))
dim(plt_matlab_VOI_surf)

write.csv(rbind(c(-1,VOI_pos_bin_plt_mid),cbind(VPI_pos_bin_plt_mid,plt_matlab_VOI_surf)),
          file = 'pdf VOIvs VPI.csv')


plt_matlab_VOI_zero[is.na(plt_matlab_VOI_zero)] = -1
write.csv(plt_matlab_VOI_zero,file = 'pmf VOIvs VPI.csv')


####################################
##posterior PMF plots for different X out of n values
##prior
PMF_VOI
##likelihood
p_l_bin_mat
dim(p_l_bin_mat_check)

##posterior for expected frequencies of f_good
PMF_VOI_post_mat = matrix(NaN,nrow = length(f_good),ncol = length(x_wells_list))
for (j in 1:length(x_wells_list)){
  ##compute posterior based on expected prior frequencies (wrong!!!)
  # PMF_VOI_post_mat[,j] = p_l_bin_mat[,j]*PMF_VOI/sum(p_l_bin_mat[,j]*PMF_VOI)
  ##posterior for each realization; then take expected value over all realizations with p_VOI_vec
  PMF_VOI_post_mat[,j] = colSums(t(p_l_bin_mat[,j]*t(prob_mat_A1))/p_l_bin_mat_check[,j]*p_VOI_vec)
}
colSums(PMF_VOI_post_mat)



VOI_realization_check = sum(rowSums(pmax(E_p_over_c_post,0)*p_l_bin_mat_check*p_VOI_vec)-
  pmax(rowSums(t(t(prob_mat_A1)*p_over_c[,1])*p_VOI_vec),0))
VOI_realization_check




##compute E(VOI) with PMF_VOI
##wrong way
# PMF_VOI*p_l_bin_mat[,1]/sum(PMF_VOI*p_l_bin_mat[,1])


##right way
# colSums(t(t(prob_mat_A1)*p_l_bin_mat[,1])/colSums(t(prob_mat_A1)*p_l_bin_mat[,1])*p_VOI_vec)






##plot the likelihood
##bar chart
for (j in 1:length(x_wells_list)){
  barplot(p_l_bin_mat[,j],names.arg = f_good,ylim = c(0,1),
          col = 'orange',ylab = 'Likelihoods',cex.lab = 1.2,cex.axis = 1.2,
          cex.names = 1.2)
  title(xlab = 'f: frequency of good wells in the play',cex.lab = 1.2)
  text(5,0.8,paste0(x_wells_list[j],' good wells out of 5 test wells'),cex = 1.2)
}



# ##continuous
# plot(f_good,p_l_bin_mat[,1],type = 'l',lwd = 3,ylim = c(0,1),xlab = '',
#      ylab = 'Likelihoods',cex.axis = 1.2,cex.lab = 1.2,xlim = c(0,1))
# title(xlab = 'f: frequency of good wells in the play',cex.lab = 1.2)
# for (j in 2:length(x_wells_list)){
#   lines(f_good,p_l_bin_mat[,j],lwd = 3,lty = j)
# }
# legend('top',c('0 good wells out of 5 test wells',
#                   '1 good wells out of 5 test wells',
#                   '2 good wells out of 5 test wells',
#                   '3 good wells out of 5 test wells',
#                   '4 good wells out of 5 test wells',
#                   '5 good wells out of 5 test wells'),cex = 1.2,bty = 'n',
#        lwd = 3,lty = seq(1,6))


##plot the posterior
for (j in 1:length(x_wells_list)){
  barplot(rbind(PMF_VOI,PMF_VOI_post_mat[,j]),
          names.arg = f_good,beside = TRUE,ylim = c(0,1),col = c('gray40','blue'),
          ylab = 'PMF',cex.lab = 1.2,cex.axis = 1.2,cex.names = 1.2)
  title(xlab = 'f: frequency of good wells in the play',cex.lab = 1.2)
  text(20,0.9,paste0('r/c = ',r_over_c),cex = 1.2)
  legend('topleft',c('Prior PMF','Posterior PMF'),col = c('gray40','blue'),cex = 1.2,
         bty = 'n',pch = c(15,15))
  text(20,0.6,paste0(x_wells_list[j],' good wells out of 5 test wells'),cex = 1.2)
}


# ##PDF
# plot(f_good,E_joint_A1_mat[,1]/diff(f_good)[1],type = 'l',lwd = 3,xlab = '',
#      ylab = 'PDF',cex.axis = 1.2,cex.lab = 1.2,ylim = c(0,15),col = 'red')
# title(xlab = 'f: frequency of good wells in the play',cex.lab = 1.2)
# for (j in 1:length(x_wells_list)){
#   lines(f_good,PMF_VOI_post_mat[,j]/diff(f_good)[1],lwd = 3,lty = j)
# }
# legend('top',c('0 good wells out of 5 test wells',
#                '1 good wells out of 5 test wells',
#                '2 good wells out of 5 test wells',
#                '3 good wells out of 5 test wells',
#                '4 good wells out of 5 test wells',
#                '5 good wells out of 5 test wells'),cex = 1.2,bty = 'n',
#        lwd = 3,lty = seq(1,6))
# text(0.15,4.8,paste0('r/c = ',r_over_c),cex = 1.2)







#######
##E(VPI) versus r/c
plot(1/r_over_c,mean_VPI_A1_vec/r_over_c,type = 'l',lwd = 3,xlab = '',ylab = '',
     cex.lab = 1.2,cex.axis = 1.2,ylim = c(0,0.14),)
title(xlab = expression(paste(f^{'*'},": break-even frequency of good wells in the play")),
      cex.lab = 1.2)
title(ylab = expression(paste('Value of Perfect Information / ',r[max])),cex.lab = 1.2,
      line = 2.5)

######
##PDF
plot(f_good,E_joint_A1_mat[,which(r_over_c==5)]/diff(f_good)[1],type = 'l',lwd = 3,xlab = '',
     ylab = 'PDF',cex.lab = 1.2,cex.axis = 1.2,ylim = c(0,4))
title(xlab = 'f: frequency of good wells in the play',cex.lab = 1.2)
lines(f_good,E_joint_A1_mat[,which(r_over_c==2)]/diff(f_good)[1],lwd = 3,lty = 2)
lines(f_good,E_joint_A1_mat[,which(r_over_c==1)]/diff(f_good)[1],lwd = 3,lty = 3)
legend('topright',c('r/c = 5','r/c = 2','r/c = 1'),lwd = 3,lty = seq(1,3),cex = 1.2,
       bty = 'n')

###############
##PMF
# plot(f_good,E_joint_A1_mat[,which(r_over_c==5)],type = 'h',lwd = 3,xlab = '',
#      ylab = 'PDF',cex.lab = 1.2,cex.axis = 1.2,ylim = c(0,4))
# title(xlab = 'f: frequency of good wells in the play',cex.lab = 1.2)
# lines(f_good,E_joint_A1_mat[,which(r_over_c==2)]/diff(f_good)[1],lwd = 3,lty = 2)
# lines(f_good,E_joint_A1_mat[,which(r_over_c==1)]/diff(f_good)[1],lwd = 3,lty = 3)
# legend('topright',c('r/c = 5','r/c = 2','r/c = 1'),lwd = 3,lty = seq(1,3),cex = 1.2,
#        bty = 'n')








##E(f) for different r_over_c
mean_f_good = colSums(f_good*E_joint_A1_mat)
mean_f_good

plot(r_over_c,mean_f_good,type = 'l')



##if 2 good wells out of 5
n_wells = 5
x_good_wells = 2
p_l_bin = choose(5,2)*f_good^(x_good_wells)*(1-f_good)*(n_wells-x_good_wells)
E_joint_post_A1_mat = sweep(p_l_bin*E_joint_A1_mat,MARGIN = 2,
                            colSums(p_l_bin*E_joint_A1_mat),FUN = '/')


##likelihood
plot(f_good,p_l_bin,lwd = 3,type = 'l',cex.lab = 1.2,cex.axis = 1.2,
     xlab = '',ylab = '',xlim = c(0,1))
title(xlab = 'f: frequency of good wells in the play',cex.lab = 1.2)
title(ylab = 'Likelihood Function: 2 good wells out of 5',cex.lab = 1.2)


##posterior
plot(f_good,E_joint_post_A1_mat[,which(r_over_c==5)]/diff(f_good)[1],
     type = 'l',xlab = '',ylab = 'Posterior PDF',cex.lab = 1.2,cex.axis = 1.2,
     lwd = 3,xlim = c(0,1),ylim = c(0,2))
title(xlab = 'f: frequency of good wells in the play',cex.lab = 1.2)
lines(f_good,E_joint_post_A1_mat[,which(r_over_c==2)]/diff(f_good)[1],lwd = 3,lty = 2)
lines(f_good,E_joint_post_A1_mat[,which(r_over_c==1)]/diff(f_good)[1],lwd = 3,lty = 3)
legend('topleft',c('r/c = 5','r/c = 2','r/c = 1'),lwd = 3,lty = seq(1,3),cex = 1.2,
       bty = 'n')



#######################################
##VOI analysis
##n test wells, x good wells
n_wells = 5

###figure 9
##VOI vs n_wells
r_over_c

n_wells_list = seq(1,n_wells)
VOI_list_n_wells = matrix(NaN,nrow = length(n_wells_list),ncol = length(r_over_c))
p_likelihood_check = matrix(NaN,nrow = length(n_wells_list),ncol = length(r_over_c))

E_joint_post_A1_mat = matrix(NaN,nrow = length(f_good),
                             ncol = length(r_over_c))

pro_VOI = txtProgressBar(1,length(r_over_c),style = 3)
t_begin = Sys.time()
for (k in 1:length(r_over_c)){
  setTxtProgressBar(pro_VOI,k)
  for (j in 1:length(n_wells_list)){
    n_wells = n_wells_list[j]
    x_good_wells_list = seq(0,n_wells)
    
    E_p_over_c_list = rep(NaN,length(x_good_wells_list))
    p_l_bin_list_A1 = rep(NaN,length(x_good_wells_list))
    
    for (i in 1:length(x_good_wells_list)){
      p_l_bin = choose(n_wells,x_good_wells_list[i])*
        f_good^(x_good_wells_list[i])*
        (1-f_good)^(n_wells-x_good_wells_list[i])
      p_l_bin_list_A1[i] = sum(p_l_bin*E_joint_A1_mat[,k])
      E_joint_post_A1_mat[,k] = p_l_bin*E_joint_A1_mat[,k]/p_l_bin_list_A1[i]
      
      E_p_over_c_list[i] = max(sum(E_joint_post_A1_mat[,k]*p_over_c[,k]),0)
    }
    VOI_list_n_wells[j,k] = sum(p_l_bin_list_A1*E_p_over_c_list) - 
      max(sum(E_joint_A1_mat[,k]*p_over_c[,k]),0)
    
    p_likelihood_check[j,k] = sum(p_l_bin_list_A1)
  }
}
Sys.time() - t_begin

##check
VOI_list_n_wells

summary(VOI_list_n_wells)

summary(p_likelihood_check)


##plot VOI vs n_test wells
plot(c(0,n_wells_list),c(0,VOI_list_n_wells[,which(r_over_c==5)]/5),type = 'l',lwd = 3,
     ylim = c(0,0.15),xlab = 'Number of test wells',ylab = '',cex.axis = 1.2,
     cex.lab = 1.2)
title(ylab = expression(paste('Value of Imperfect Information / ',r[max])),
      cex.lab = 1.2,line = 2.5)
lines(c(0,n_wells_list),c(0,VOI_list_n_wells[,which(r_over_c==1.25)]/1.25),lty = 2)
lines(c(0,n_wells_list),c(0,VOI_list_n_wells[,which(r_over_c==2)]/2),lty = 3,lwd = 3)
legend('right',c('r/c = 5','r/c = 1.25','r/c = 2'),lwd = 3,lty = seq(1,2,3),cex = 1.2,
       bty = 'n')



plot(c(0,n_wells_list),c(0,VOI_list_n_wells[,which(r_over_c==1.25)]/1.25),type = 'l',lwd = 3,
     ylim = c(0,0.15),xlab = 'Number of test wells',ylab = '',cex.axis = 1.2,lty = 2,
     cex.lab = 1.2)

#####
##VOI adjusted prob
E_joint_post_A1_mat

colSums(E_joint_post_A1_mat)

# index_bar_plt = c(1,1.1,2,5,seq(18,19,.5),50,100)
index_bar_plt  =c(1.25,5)
# for (i in 1:length(index_bar_plt)){
#   index_col = which(r_over_c==index_bar_plt[i])
#   # barplot(E_joint_A1_mat[,index_col],name= 'f: frequency of good wells in the play',cex.lab = 1.2)
#   plot(f_good,E_joint_A1_mat[,index_col]/diff(f_good)[1],type = 'l',ylim = c(0,2))
#   text(0.8,1.8,paste0('r/c = ',index_bar_plt[i]),cex = 1.2)
# }
for (i in 1:length(index_bar_plt)){
  index_col = which(r_over_c==index_bar_plt[i])
  barplot(rbind(E_joint_A1_mat[,index_col],E_joint_post_A1_mat[,index_col]),
          names.arg = f_good,beside = TRUE,ylim = c(0,1),col = c('gray40','blue'),
          ylab = 'PMF',cex.lab = 1.2,cex.axis = 1.2,cex.names = 1.2)
  title(xlab = 'f: frequency of good wells in the play',cex.lab = 1.2)
  text(20,0.9,paste0('r/c = ',index_bar_plt[i]),cex = 1.2)
  legend('topleft',c('Prior PMF','Posterior PMF'),col = c('gray40','blue'),cex = 1.2,
         bty = 'n',pch = c(15,15))
}








