library(ncdf4)

TAS_Models=list()
TAS_Models_hist_rcp=list()
Models_bdd=list.files(path='/bdd/CMIP5/output/')
m=1
for (i in 1:length(Models_bdd)){
  Model=Models_bdd[i]
  
  #Listing model versions
  Model_versions_path=paste0('/bdd/CMIP5/output/',Model,'/')
  Model_versions=list.files(path=Model_versions_path)
  
  for (j in 1:length(Model_versions)){
    Model_version=Model_versions[j]
    
    #Listing experiments
    Experiments_path=paste0('/bdd/CMIP5/output/',Model,'/',Model_version,'/')
    Experiments=list.files(path=Experiments_path)
    Historical_filepath_tas=paste0('/bdd/CMIP5/output/',Model,'/',Model_version,'/historical/mon/atmos/Amon/r1i1p1/latest/tas/')
    RCP85_filepath_tas=paste0('/bdd/CMIP5/output/',Model,'/',Model_version,'/rcp85/mon/atmos/Amon/r1i1p1/latest/tas/')
    
    if(dir.exists(Historical_filepath_tas)==TRUE & dir.exists(RCP85_filepath_tas)==TRUE)    {    
      TAS_Models[[m]]=Model
      TAS_Models_hist_rcp[[m]]=Model_version                
      m=m+1
    }
  }
}

names(TAS_Models_hist_rcp)=TAS_Models

TAS_MODELS=matrix(0,nrow=length(TAS_Models),ncol=2)

for (i in 1:length(TAS_Models)){
  TAS_MODELS[i,1]=TAS_Models[[i]];TAS_MODELS[i,2]=TAS_Models_hist_rcp[[i]]
}


saveRDS(TAS_MODELS, file = 'models_list.rds')


