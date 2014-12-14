//
// test the ExpressionInternal class' gaussian random number generator
// Check mean, standard deviation, skew and kurtosis
// skew should be close to zero.  kurtosis should be close to 3.0
// mean and standard deviation should be close to specified values.
// "Close to" is determined by rectal extraction.
//
#include <Xyce_config.h>
#include <N_UTL_Misc.h>
#include <N_UTL_fwd.h>
#include <N_UTL_RandomNumbers.h>

#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_BlockMap.h>

#include <iostream>
#include <vector>
#include <cmath>

int main(int argc, char* argv[])
{

  Xyce::Util::RandomNumbers theRandomNumbers;
  double desiredMean;
  double desiredSD;


  double sum, sumsq, temp,stddev, mean, meansq;
  double sumcub,sumquad, meancub, meanquad;
  double mu3, mu4;
  sum=0; sumsq=0;
  sumcub=0; sumquad=0;

  theRandomNumbers.seedRandom(0);

  desiredMean=3;
  desiredSD=0.1;

  for (int i=0; i<1000000; i++)
  {
    temp=theRandomNumbers.gaussianRandom(desiredMean,desiredSD);
    sum +=temp;
    sumsq += temp*temp;
    sumcub += temp*temp*temp;
    sumquad += temp*temp*temp*temp;
    mean=sum/(i+1.0);
    meansq=sumsq/(i+1.0);
    meancub=sumcub/(i+1.0);
    meanquad=sumquad/(i+1.0);
    stddev = sqrt(meansq-mean*mean);
    mu3=meancub-3*mean*meansq+2*mean*mean*mean;
    mu4=meanquad-4*mean*meancub+6*mean*mean*meansq-3*mean*mean*mean*mean;

    std::cout << "mean:" << mean  << "    std_dev:" << stddev  
              << "    skew:" << mu3/(stddev*stddev*stddev)
              << "    kurt:" << mu4/(stddev*stddev*stddev*stddev)
              << "\r";
  }
  std::cout << "\n";

  if (fabs(mean-desiredMean)<1e-3 && 
      fabs(stddev-desiredSD)<1e-3 && 
      fabs(mu3/(stddev*stddev*stddev)) <1e-2 && 
      fabs(mu4/(stddev*stddev*stddev*stddev)-3)<1e-2)
     exit(0);
  else
     exit(1);
}
