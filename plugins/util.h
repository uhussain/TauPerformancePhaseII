// -*- C++ -*-                
//                            
// Package:     -
// File:        util.h       
//                            
//                            

#include<cmath>

bool pairCompare(const std::pair<float, unsigned>& firstElem, const std::pair<float, unsigned>& secondElem) {
    return firstElem.first < secondElem.first;

}
