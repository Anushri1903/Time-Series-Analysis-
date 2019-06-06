# applying Inbuild GA on Rossler dataset

library(tseriesChaos)
require(graphics)
ap<-sim.cont(rossler.syst, start=0, end= 650, dt=0.1,start.x=c(0,0,0), parms=c(0.15,0.20,10))   #sunspots dataset
ap<-ap[1:6500]
convo<-c()
t<-c()
t1<-c()
t2<-c()
nr<-c()
nr1<-c()
nr2<-c()
mse<-c()
fit<-c()
offsp_delay<-c()
offsp_dimen<-c()
intervtrain<-list()
intervtest<-list()
del<-c()
dimn<-c()

l<-1625
ntrain<-1219
ntest<-406
sig<-c()
std1<-c()
n_step=4
l1<-3/4  # part of data taken as training set
l2<-1/4  # part of data taken as test set

ranges<-list()
sumfit<-rep(0,5)

d=1
Rastrigin<-function(dl1,dm1,dl2,dm2,dl3,dm3,dl4,dm4)
{
  actions_r<-c(dl1,dm1,dl2,dm2,dl3,dm3,dl4,dm4)
  print(actions_r)
  zero<-1
  for( stp in 1:n_step)
  {  
    
    t[stp]<-zero # start with t=0
    nr[stp]<-1625# random number generated
    t1[stp]<-t[stp] # training set of l1 fraction
    nr1[stp]<-l1*nr[stp]
    t2[stp]<-t1[stp]+ceiling(nr1[stp])
    nr2[stp]<-l2*nr[stp]
    
    intervtrain[[stp]]<-c(t1[stp],ceiling(t1[stp]+nr1[stp]))
    training_set<-ap[intervtrain[[stp]][1]:intervtrain[[stp]][2]] # training set
    intervtest[[stp]]<-c(t2[stp],ceiling(t2[stp]+nr2[stp]))
    testing_set<-ap[intervtest[[stp]][1]:intervtest[[stp]][2]] # testing set
    t<-d+1
    mutut<-actions_r[d:t]
    Initial_val<-fitness(mutut,training_set,testing_set)
    convo<-cbind(convo,Initial_val)
    zero<-zero+ceiling(nr1[stp])  # moving on data 
    d=d+2
  }
  sum_of_fitness<-sum(convo)
  return(sum_of_fitness)
}

input_to_Rf<-function(out,train,delay,dim,m) # build training set for Random Forest
{
  input<-matrix(0,length(out),dim)
  
  for(k in 1:length(out))
  { 
    
    for(j in 1:dim)
    { 
      
      input[k,j]<-train[(m+k)-(j*delay)] 
    }
  }
  
  return(input)
}
fitness<-function(mutut,train,test) # fitness function
{
    actions<-round(mutut)
    m<-(actions[1]*actions[2])                          # m= delay*dimension
    trainout<-train[(m+1):ntrain]
    testout<-test[(m+1):ntest]
    inputtrain<-data.frame()
    inputtest<-data.frame()
    inputtrain<-input_to_Rf(trainout,train,actions[1],actions[2],m)  # input to NN is generated based on delay and dimension
    inputtest<-input_to_Rf(testout,test,actions[1],actions[2],m)
    inputrnn<-cbind(inputtrain,trainout)
    inputsnn<-cbind(inputtest,testout)
    l<-actions[2]
    inputrnn<-as.data.frame(inputrnn)
    xname<-paste0(names(inputrnn[1:l]), collapse = '+') # for NN
    formula = as.formula(paste0("trainout~",xname))
    model<-randomForest(formula = formula,data=inputrnn)
    pred<-predict(model,inputtest)
    mse<-mean((pred-testout) ^2)   # fitness function
   #fitval<-mse/mean(testout)
   #print(fitval)
  return(mse)  # returns best fittest parents
}
GA <- ga(type = "real-valued",
         fitness = function(x) -Rastrigin(x[1], x[2], x[3],x[4],x[5],x[6],x[7],x[8]),
         lower = c(1,1,1,1,1,1,1,1), upper = c(10,15,10,15,10,15,10,15),
         popSize =15, maxiter = 200,
         pmutation = ga_pmutation)
plot(GA)
sink("GA_statRossler.txt")
print(summary(GA))
sink()