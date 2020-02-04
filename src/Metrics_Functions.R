#author: Reinier de Vries
#Aim: function descriptions for LiDAR metrics 

# Function descriptions
# Calculations of single-grid cell (1m) measures

Heightmetrics = function(z)
{
  heightmetric = list(
    zmax = max(z), # single highest point
    zmean = mean(z),
    z090quantile = quantile(z, 0.90)
  )
  return(heightmetric)
}


Amplmetrics = function(i)
{
  amplmetric = list(
    istd = sd(i),
    imean = mean(i)
  )
  return(amplmetric)
}


#VertDistr_Metrics = function(z)
#{
#  library("e1071")
#  
#  p=proportion(z, by = 1)
#  p_whnull=p[p>0]
#  
#  vertdistr_metrics = list(
#    zstd = sd(z),
#    zkurto = kurtosis(z)
#  )
#  return(vertdistr_metrics)
#}

# Inputs for functions: LidR-LAS-class