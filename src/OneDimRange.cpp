// @author Wushi Dong

// OneDimRange.cpp

#include "OneDimRange.h"

OneDimRange::~OneDimRange(){}

void OneDimRange::Print() const
{
  for(int i = 0; i < mN; ++i)
  {
    std::cout << mStart + (double)i * mSpacing << " ";
    std::cout << std::endl;
  }
}

//int main()
//{
//  OneDimRange <double> one_dim_range_double(0., 10., 2.5);
//  cout << one_dim_range_double.mStart << endl;
//  cout << one_dim_range_double.mEnd << endl;
//  cout << one_dim_range_double.mN << endl;
//  cout << one_dim_range_double.mSpacing << endl;
//
//  OneDimRange <double> one_dim_range_int(0, 10, 1);
//  cout << one_dim_range_int.mStart << endl;
//  cout << one_dim_range_int.mEnd << endl;
//  cout << one_dim_range_int.mN << endl;
//  cout << one_dim_range_int.mSpacing << endl;
//  
//  return 1;
//}

