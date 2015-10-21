# function to dicotamise a vector based on whether each value is below or above the sample median

DicotamiseOnMedian <- function(data){

  med <- median(data)
  val <- vector()

  for(i in 1:length(data)){
    
    if(data[i] < med){
      val[i] <- 0
       
    }else{
      
      val[i] <- 1
    
    }
    
  }

  return(val)
  
  } 