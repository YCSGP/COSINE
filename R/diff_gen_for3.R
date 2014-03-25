diff_gen_for3 <-
function(data1,data2,data3){

   num_sample<-dim(data1)[1]
   num_gene<-dim(data1)[2]
   diff_expr<-rep(0,num_gene)
   diff_coex<-matrix(0,ncol=num_gene,nrow=num_gene)

   # calculate the statistics measuing the differential expression of each gene between the 2 groups

   for(i in 1:num_gene){
       data<-c(data1[,i],data2[,i],data3[,i])
       type<-c(rep(1,num_sample),rep(2,num_sample),rep(3,num_sample))
       diff_expr[i]<-f.test(data,type)
   }

   # calculate the statistics measuing the differential coexpression of each gene-pair between the 2 groups

   for(i in 1:(num_gene-1)){
        for(j in (i+1):num_gene){
           data.x<-c(data1[,i],data2[,i],data3[,i])
           data.y<-c(data1[,j],data2[,j],data3[,j])
           cond_fyx<-cond.fyx(data.y,data.x,type)
           cond_fxy<-cond.fyx(data.x,data.y,type)
           diff_coex[i,j]<-(cond_fyx+cond_fxy)/2
           diff_coex[j,i]<-diff_coex[i,j]
        }
   }

   return(list(diff_expr,diff_coex))
}

