#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

bool reinjection(double xn1, double xn, double ubr, double lbr) {

    if (xn1 < ubr && xn1 > lbr && xn > ubr && xn < lbr){
        return true;
    }
    else {
        return false;
    }
}

bool ejection(double xn1, double xn, double ubr, double lbr) {

    if (xn1 > ubr && xn1 < lbr && xn < ubr && xn > lbr) {
        return true;
    }
    else {
        return false;
    }
}