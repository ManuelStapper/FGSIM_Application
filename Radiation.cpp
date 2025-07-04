#include <Rcpp.h>
using namespace Rcpp;

// Helper function to sum matrix elements within a specific range using a loop
// Switch x and y to try...
double sumMatrixRange(NumericMatrix popMat, int col, int rowL, int rowU) {
  int rowMax = popMat.nrow();
  if (rowL < 1) rowL = 1;
  if (rowU > rowMax) rowU = rowMax;
  
  // Adjust 1-based indexing to 0-based for C++
  col = col - 1;
  rowL = rowL - 1;
  rowU = rowU - 1;
  
  double total = 0.0;
  
  // Loop through the y range and sum up the values in the matrix
  for (int row = rowL; row <= rowU; ++row) {
    total += popMat(row, col);
  }
  
  return total;
}

// [[Rcpp::export]]
double computePop(int rowO, int colO, int r, NumericMatrix popMat) {
  int colMax = popMat.ncol();
  
  double t1 = r / 16.0;
  double t2 = 0;
  int x  = r;
  int y  = 0;
  IntegerVector p1 = {r, 0};
  IntegerVector p2 = {0, r};
  
  NumericVector done(r + 1, false);
  double out = 0.0;
  
  // Sum for the initial column at xO (the center of the circle)
  out += sumMatrixRange(popMat, colO, rowO - r, rowO + r);
  
  while (x >= y) {
    IntegerVector p1New = {x, y};
    IntegerVector p2New = {y, x};
    
    if (p1New[0] != p1[0]) {
      if ((colO + p1[0]) <= colMax) {
        out += sumMatrixRange(popMat, colO + p1[0], rowO - p1[1], rowO + p1[1]);
      }
      if ((colO - p1[0]) > 0) {
        out += sumMatrixRange(popMat, colO - p1[0], rowO - p1[1], rowO + p1[1]);
      }
      done[p1[0]] = true;
    }
    
    if ((colO + p2New[0]) <= colMax) {
      out += sumMatrixRange(popMat, colO + p2New[0], rowO - p2New[1], rowO + p2New[1]);
    }
    
    if ((colO - p2New[0]) > 0) {
      out += sumMatrixRange(popMat, colO - p2New[0], rowO - p2New[1], rowO + p2New[1]);
    }
    done[p2New[0]] = true;
    
    p1 = p1New;
    p2 = p2New;
    y  = y + 1;
    t1 = t1 + y;
    t2 = t1 - x;
    if (t2 >= 0) {
      t1 = t2;
      x  = x - 1;
    }
  }
  
  // Final iteration to make sure no points are left unsummed
  IntegerVector p1New = {x, y};
  IntegerVector p2New = {y, x};
  
  if (p1New[0] != p1[0] && !done[p1[0]]) {
    if ((colO + p1[0]) <= colMax) {
      out += sumMatrixRange(popMat, colO + p1[0], rowO - p1[1], rowO + p1[1]);
    }
    
    if ((colO - p1[0]) > 0) {
      out += sumMatrixRange(popMat, colO - p1[0], rowO - p1[1], rowO + p1[1]);
    }
    done[p1[0]] = true;
  }
  
  if (!done[p2New[0]]) {
    if ((colO + p2New[0]) <= colMax) {
      out += sumMatrixRange(popMat, colO + p2New[0], rowO - p2New[1], rowO + p2New[1]);
    }
    
    if ((colO - p2New[0]) > 0) {
      out += sumMatrixRange(popMat, colO - p2New[0], rowO - p2New[1], rowO + p2New[1]);
    }
  }
  
  return out;
}