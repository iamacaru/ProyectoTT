#include "../include/Sign_.h"

using namespace std;

double sign_(double a, double b) {
    if (b >= 0.0) {
        return abs(a);
    } else {
        return - abs(a);
    }
}