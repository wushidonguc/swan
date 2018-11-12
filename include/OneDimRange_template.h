//
#ifndef ONEDIMRANGE_H
#define ONEDIMRANGE_H

#include "header.h"

/* Potential profile class */
template <class numtype> class OneDimRange
{
  public:
//    OneDimRange(numtype start, numtype end, int n): mStart(start), mEnd(end), mN(n), mSpacing((end - start) / double(n)){}
    OneDimRange(numtype start, numtype end, numtype spacing): mStart(start), mEnd(end), mN(int((end - start) / spacing)), mSpacing(spacing){}
//    ~OneDimRange();

    numtype mStart;
    numtype mEnd;
    int mN;
    numtype mSpacing;

//    void Print() const;
};

#endif

