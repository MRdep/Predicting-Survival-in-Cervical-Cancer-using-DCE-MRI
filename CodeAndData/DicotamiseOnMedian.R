# function to dicotamise a vector based on whether each value is below or above the sample median

DicotamiseOnIndex <- function(data){

  splitpoint <- median(data)
  val <- vector()

  for(i in 1:length(data)){
    
    if(data[i] < splitpoint){
      val[i] <- 0
       
    }else{
      
      val[i] <- 1
    
    }
    
  }

  return(val)
  
  } 