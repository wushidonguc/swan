// @author Wushi Dong

// OneDimRange.h
// Purpose: Class of one-dimensional range containing min, max, spacing, and 
// number of points. It can also serve as base class for specified
// one-dimensional properties.

#ifndef ONEDIMRANGE_H
#define ONEDIMRANGE_H

#include <iostream>

class OneDimRange
{
  public:
    // Constructor
    OneDimRange(double start, double end, int n): mStart(start), mEnd(end),
    mN(n), mSpacing((end - start) / double(n)){}
    OneDimRange(double start, double end, double spacing): mStart(start),
    mEnd(end), mN(int((end - start) / spacing) + 1), mSpacing(spacing){}

    // Destructor
    ~OneDimRange();

    double mStart;
    double mEnd;
    int mN;
    double mSpacing;

    void Print() const;
};

#endif

